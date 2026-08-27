#!/usr/bin/env python3
"""Audit the immutable Qualification v2 panel without network access.

Scientific engines run only after every adopted record has exact source and
expectation metadata. Until then this runner emits a deterministic blocked
report rather than converting missing public cases into favorable evidence.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
from pathlib import Path
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
REQUIRED_RECORD_FIELDS = {
    "record_id", "record_kind", "stratum", "split", "expected_result",
    "source_release", "checksum", "license", "citation", "artifact",
    "exclusion_rationale", "pdb_entry_id", "homolog_group",
}


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical(value))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_artifact(relative: str) -> Path | None:
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        return None
    target = (HERE / path).resolve()
    try:
        target.relative_to(HERE)
    except ValueError:
        return None
    return target


def verify_source_lock(lock: Mapping[str, Any]) -> dict[str, Any]:
    prior = lock.get("prior_collection", {})
    prior_path = (HERE / str(prior.get("manifest", ""))).resolve()
    prior_ok = prior_path.is_file() and sha256(prior_path) == str(prior.get("sha256", ""))
    records = []
    for source in lock.get("sources", []):
        target = safe_artifact(str(source.get("artifact", "")))
        observed = sha256(target) if target is not None and target.is_file() else None
        records.append({
            "source_id": str(source.get("source_id", "")),
            "artifact": str(source.get("artifact", "")),
            "expected_sha256": str(source.get("sha256", "")),
            "observed_sha256": observed,
            "passed": observed == str(source.get("sha256", "")),
        })
    return {
        "prior_collection_checksum_valid": prior_ok,
        "adopted_source_count": len(records),
        "adopted_sources_valid": all(row["passed"] for row in records),
        "sources": records,
    }


def matches(record: Mapping[str, Any], requirement: Mapping[str, Any]) -> bool:
    return all(record.get(key) == value for key, value in requirement.items() if key != "count")


def validate_panel(panel: Mapping[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    records = panel.get("records", [])
    errors: list[str] = []
    if not isinstance(records, list):
        records = []
        errors.append("records must be a list")
    seen_ids: set[str] = set()
    seen_split_entries: dict[str, set[str]] = {"development": set(), "held_out": set()}
    seen_split_groups: dict[str, set[str]] = {"development": set(), "held_out": set()}
    rows = []
    for record in records:
        missing = sorted(REQUIRED_RECORD_FIELDS - set(record)) if isinstance(record, Mapping) else sorted(REQUIRED_RECORD_FIELDS)
        if missing:
            errors.append(f"record is missing fields: {', '.join(missing)}")
            continue
        record_id = str(record["record_id"])
        if record_id in seen_ids:
            errors.append(f"duplicate record_id: {record_id}")
        seen_ids.add(record_id)
        split = str(record["split"])
        if split in seen_split_entries:
            seen_split_entries[split].add(str(record["pdb_entry_id"]).upper())
            seen_split_groups[split].add(str(record["homolog_group"]))
        target = safe_artifact(str(record["artifact"]))
        observed = sha256(target) if target is not None and target.is_file() else None
        checksum_ok = observed == str(record["checksum"])
        if not checksum_ok:
            errors.append(f"{record_id}: source artifact is missing or checksum-mismatched")
        rows.append({
            "workflow": panel["workflow"], "panel_id": panel["panel_id"],
            "record_id": record_id, "record_kind": record["record_kind"],
            "stratum": record["stratum"], "split": split,
            "pdb_entry_id": record["pdb_entry_id"], "source_release": record["source_release"],
            "checksum_valid": checksum_ok, "execution_state": "not_executed_panel_audit_only",
        })
    leaked_entries = sorted(seen_split_entries["development"] & seen_split_entries["held_out"])
    leaked_groups = sorted(seen_split_groups["development"] & seen_split_groups["held_out"])
    if leaked_entries:
        errors.append("PDB entry leakage across development and held-out splits: " + ", ".join(leaked_entries))
    if leaked_groups:
        errors.append("homolog-group leakage across development and held-out splits: " + ", ".join(leaked_groups))
    requirements = []
    for requirement in panel.get("requirements", []):
        observed = sum(1 for record in records if isinstance(record, Mapping) and matches(record, requirement))
        required = int(requirement["count"])
        requirements.append({
            **requirement, "observed_count": observed,
            "missing_count": max(required - observed, 0), "passed": observed == required,
        })
        if observed != required:
            errors.append(
                f"{requirement['record_kind']}:{requirement['stratum']}:{requirement['split']} "
                f"requires {required}, observed {observed}"
            )
    return {
        "workflow": panel["workflow"], "panel_id": panel["panel_id"],
        "state": "ready_for_execution" if not errors else "blocked_panel_incomplete",
        "record_count": len(records), "requirements": requirements,
        "errors": sorted(errors), "gates": panel.get("gates", {}),
    }, rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def render_report(result: Mapping[str, Any]) -> str:
    sections = []
    for panel in result["panels"]:
        requirement_rows = "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(row[key]))}</td>" for key in ("record_kind", "stratum", "split", "count", "observed_count", "missing_count")) + "</tr>"
            for row in panel["requirements"]
        )
        sections.append(
            f"<section><h2>{html.escape(panel['workflow'])}</h2><p><strong>{html.escape(panel['state'])}</strong></p>"
            "<table><thead><tr><th>Record</th><th>Stratum</th><th>Split</th><th>Required</th><th>Observed</th><th>Missing</th></tr></thead>"
            f"<tbody>{requirement_rows}</tbody></table></section>"
        )
    return """<!doctype html><html lang="en"><head><meta charset="utf-8"><title>YAUVI Mark 1 Qualification v2</title>
<style>body{font:14px/1.5 system-ui;max-width:1050px;margin:32px auto;padding:0 20px;color:#17211e}h1,h2{color:#064f43}section{break-inside:avoid;margin:24px 0}table{width:100%;border-collapse:collapse}th,td{border:1px solid #b9c8c2;padding:6px;text-align:left}th{background:#eaf3ef}.notice{padding:12px;border-left:5px solid #b36b00;background:#fff5df}@media print{button{display:none}body{margin:0}}</style></head><body>
<button onclick="window.print()">Print / Save as PDF</button><h1>YAUVI Mark 1 Scientific Qualification v2</h1>
<p class="notice">This is a panel-composition audit, not a scientific qualification pass. Missing public cases remain blocked evidence.</p>
<p>Release state: <strong>""" + html.escape(str(result["overall_state"])) + "</strong></p>" + "".join(sections) + """
<h2>Non-claims</h2><ul><li>Orientation is not native exposure.</li><li>Conformational resemblance is not biochemical activity.</li><li>A mapped site is not observed catalysis.</li><li>An interface is not affinity.</li><li>Similarity is not exact functional transfer.</li></ul></body></html>"""


def main() -> int:
    manifest = read_json(HERE / "PANEL_MANIFEST.json")
    source_lock = verify_source_lock(read_json(HERE / "SOURCE_LOCK.json"))
    panels = []; case_rows = []
    for panel in manifest["panels"]:
        summary, rows = validate_panel(panel)
        panels.append(summary); case_rows.extend(rows)
    ready = source_lock["prior_collection_checksum_valid"] and source_lock["adopted_sources_valid"] and all(
        panel["state"] == "ready_for_execution" for panel in panels
    )
    result = {
        "schema_version": "2.0", "collection_id": manifest["collection_id"],
        "overall_state": "ready_for_scientific_execution" if ready else "blocked_panel_incomplete",
        "panel_composition_ready": ready,
        "scientific_execution_performed": False,
        "software_tests_are_separate": True,
        "source_lock": source_lock, "panels": panels,
        "limitations": [
            "This audit validates panel composition and source binding only.",
            "No scope becomes qualified until its scientific engines and predeclared gates execute successfully.",
        ],
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_json(RESULTS / "QUALIFICATION_V2_STATUS.json", result)
    write_tsv(RESULTS / "CASE_STATUS.tsv", sorted(case_rows, key=lambda row: (row["workflow"], row["record_id"])), [
        "workflow", "panel_id", "record_id", "record_kind", "stratum", "split",
        "pdb_entry_id", "source_release", "checksum_valid", "execution_state",
    ])
    stratum_rows = [
        {"workflow": panel["workflow"], **requirement}
        for panel in panels for requirement in panel["requirements"]
    ]
    write_tsv(RESULTS / "STRATUM_STATUS.tsv", stratum_rows, [
        "workflow", "record_kind", "stratum", "split", "count", "observed_count", "missing_count", "passed",
    ])
    (RESULTS / "QUALIFICATION_REPORT.html").write_text(render_report(result), encoding="utf-8")
    checksums = {
        path.name: sha256(path)
        for path in sorted(RESULTS.iterdir())
        if path.is_file() and path.name != "CHECKSUMS.json"
    }
    write_json(RESULTS / "CHECKSUMS.json", checksums)
    print(result["overall_state"])
    # Composition readiness is not scientific qualification. This runner only
    # audits panel composition and source binding; a later immutable executor
    # version must run the canonical engines and evaluate every gate. The exit
    # code therefore reports composition state, and can never report a scope as
    # qualified: 0 means the panel is composed and ready to be executed, 1 means
    # it is not. `scientific_execution_performed` stays false either way.
    return 0 if result.get("panel_composition_ready") else 1


if __name__ == "__main__":
    raise SystemExit(main())
