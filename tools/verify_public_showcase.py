#!/usr/bin/env python3
"""Verify deterministic evidence, redaction, and public-showcase boundaries."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / "yauvi-structural-workbench" / "public-showcase"
TECHNICAL = ROOT / "yauvi-structural-workbench" / "showcase" / "five-human-use-cases"
SFCSA = ROOT / "yauvi-structural-workbench" / "showcase" / "sfcsa-ceiling-case"
TEXT_SUFFIXES = {".html", ".css", ".js", ".json", ".tsv", ".md", ".txt"}


def fail(message: str) -> None:
    raise SystemExit(f"public showcase verification failed: {message}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_data() -> dict:
    text = (PUBLIC / "data.js").read_text(encoding="utf-8")
    prefix, suffix = "window.YAUVI_PUBLIC_SHOWCASE = ", ";\n"
    if not text.startswith(prefix) or not text.endswith(suffix):
        fail("data.js is not a deterministic embedded JSON assignment")
    return json.loads(text[len(prefix):-len(suffix)])


def main() -> int:
    required = {
        "index.html", "styles.css", "app.js", "data.js", "README.md",
        "PUBLIC_SHOWCASE_MANIFEST.json", "CHECKSUMS.json",
    }
    if not required.issubset({path.name for path in PUBLIC.iterdir() if path.is_file()}):
        fail("one or more public package files are missing")
    data = load_data()
    workflows, cases = data.get("workflows", []), data.get("cases", [])
    if len(workflows) != 6 or len(cases) != 6:
        fail("expected six explained workflows and six traceable cases")
    case_types = {case.get("analysis_type") for case in cases}
    expected_cases = {
        "structure_qc", "membrane_orientation", "conformational_state",
        "functional_site_state", "assembly_interface", "sf_csa",
    }
    if case_types != expected_cases:
        fail("the six case identities are incomplete")
    synthetic_cases = [case for case in cases if case.get("analysis_type") != "sf_csa"]
    if len(synthetic_cases) != 5 or any(case.get("test_state") != "passed_synthetic_case" for case in synthetic_cases):
        fail("the five synthetic analysis cases are incomplete")
    sf_case = next(case for case in cases if case.get("analysis_type") == "sf_csa")
    if sf_case.get("test_state") != "passed_stubbed_pipeline_case":
        fail("SF-CSA is not labeled as a stubbed pipeline case")
    disclosure = (sf_case.get("runtime_disclosure") or "").lower()
    if "test doubles" not in disclosure or "computed no alignments" not in disclosure:
        fail("SF-CSA runtime boundary is not explicit")
    sf_csa = next((item for item in workflows if item.get("analysis_type") == "sf_csa"), None)
    if not sf_csa or sf_csa.get("showcase_state") != "passed_stubbed_pipeline_case":
        fail("SF-CSA showcase state does not match the stubbed case")
    if "computed no alignments" not in sf_csa.get("showcase_note", ""):
        fail("SF-CSA workflow card hides its no-alignment boundary")
    baseline = json.loads((ROOT / "yauvi-structural-workbench" / "BASELINE.json").read_text(encoding="utf-8"))
    release = json.loads((ROOT / "yauvi-structural-workbench" / "RELEASE_STATUS.json").read_text(encoding="utf-8"))
    identity = json.loads((ROOT / "yauvi-structural-workbench" / "PLATFORM_IDENTITY.json").read_text(encoding="utf-8"))
    if data.get("product") != identity.get("display_name") or data.get("platform_identity") != identity:
        fail("Mark 1 platform identity is missing or drifted")
    if identity.get("edition") != "Mark 1" or identity.get("publication_authorized"):
        fail("Mark 1 identity has an invalid edition or publication state")
    if data["baseline"]["total_passed"] != baseline["total_passed"]:
        fail("software test total was hard-coded or drifted from BASELINE.json")
    if data["release"]["release_state"] != release["release_state"]:
        fail("release status drifted")
    if release["publication_authorized"] or data["release"]["publication_authorized"]:
        fail("public showcase incorrectly authorizes publication")
    qualification = data.get("qualification", {})
    if qualification.get("workflow_counts") != {"blocked": 0, "failed": 0, "partial": 2, "passed": 4}:
        fail("public qualification counts drifted")
    qualification_cases = qualification.get("cases", [])
    if len(qualification_cases) != 6:
        fail("six public qualification narratives are required")
    qualification_v2 = data.get("qualification_v2", {})
    # These were pinned to the frozen values "blocked_panel_incomplete" and
    # scientific_execution_performed is False. That guarded honestly against
    # over-claiming while nothing had executed, but it also fixed the showcase
    # to a state the project has since moved past, and it would have failed on
    # the day the last panel is adopted and the audit correctly reports a
    # composed panel. The check is now consistency with the evidence: the
    # showcase must say exactly what the executed results and the composition
    # audit say, and must never claim a qualified scope.
    v2_status = json.loads((ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
                            / "results" / "QUALIFICATION_V2_STATUS.json").read_text(encoding="utf-8"))
    summary = json.loads((ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
                          / "results" / "EXECUTION_SUMMARY.json").read_text(encoding="utf-8"))
    if qualification_v2.get("overall_state") != v2_status["overall_state"]:
        fail("Qualification v2 composition state drifted from the audit")
    execution = qualification_v2.get("scientific_execution", {})
    for field in ("panels_executed", "panels_total", "cases_passed", "cases_required",
                  "every_executed_panel_passed", "workflows_executed", "workflows_not_executed",
                  "all_release_blocking_scopes_qualified", "second_machine_reproduction"):
        if execution.get(field) != summary[field]:
            fail(f"showcase execution field '{field}' drifted from EXECUTION_SUMMARY.json")
    if qualification_v2.get("scientific_execution_performed") != summary["scientific_execution_performed"]:
        fail("showcase execution flag drifted from the executed evidence")
    # The gate that must never close from generated data.
    if execution.get("all_release_blocking_scopes_qualified") is not False:
        fail("public showcase claims Qualification v2 scopes are qualified")
    if execution.get("second_machine_reproduction") != "not_recorded":
        fail("public showcase claims an independent second-machine reproduction")
    if not str(execution.get("scope_qualification_note", "")).strip():
        fail("execution counts are published without the note that bounds them")
    if int(qualification_v2.get("missing_records", 0)) <= 0 and summary["cases_passed"] >= summary["cases_required"]:
        fail("Qualification v2 missing-record count is not visible")

    # RELEASE_STATUS.json records a sha256 beside every evidence document it
    # cites, and until now nothing compared them to the files. The membrane
    # result was rewritten by the drift change and its recorded digest went
    # stale in the public repository without any check noticing, which is the
    # whole point of recording a digest. Verify the chain it claims.
    #
    # QUALIFICATION_V2_STATUS.json is deliberately excluded. run_qualification.py
    # regenerates it, and it embeds an observed_sha256 per locked source that is
    # null wherever the artifacts have not been acquired -- so its bytes differ
    # between a machine with sources/ populated and the reviewer gate, which
    # audits composition offline. Its recorded digest cannot be a contract, and
    # asserting it here passed locally and failed CI for exactly that reason.
    workbench = ROOT / "yauvi-structural-workbench"
    v2_evidence = release["qualification_evidence"]["current_v2"]
    digests = [
        (v2_evidence["panel_manifest"], v2_evidence["panel_manifest_sha256"]),
        (v2_evidence["source_lock"], v2_evidence["source_lock_sha256"]),
        (v2_evidence["execution_summary"], v2_evidence["execution_summary_sha256"]),
    ] + [
        (path, v2_evidence["execution_results_sha256"][workflow])
        for workflow, path in sorted(v2_evidence["execution_results"].items())
    ]
    for relative, recorded in digests:
        if digest(workbench / relative) != recorded:
            fail(f"RELEASE_STATUS.json records a stale digest for {relative}")
    case_statuses = {item.get("analysis_type"): item.get("status") for item in qualification_cases}
    expected_statuses = {
        "structure_qc": "passed", "membrane_orientation": "partial",
        "conformational_state": "partial", "functional_site_state": "passed",
        "assembly_interface": "passed", "sf_csa": "passed",
    }
    if case_statuses != expected_statuses:
        fail("qualification states do not match the canonical result")
    workflow_benchmarks = {item["analysis_type"]: item.get("external_benchmark") for item in workflows}
    expected_benchmarks = {
        key: "public_case_passed" if value == "passed" else "partial_public_case"
        for key, value in expected_statuses.items()
    }
    if workflow_benchmarks != expected_benchmarks:
        fail("workflow benchmark labels do not match public qualification")
    for workflow in workflows:
        inputs = workflow.get("inputs", [])
        if not inputs or not any(item.get("required") for item in inputs):
            fail(f'{workflow.get("analysis_type")} lacks required file guidance')
        for item in inputs:
            if not item.get("label") or not item.get("extensions") or not item.get("absence_effect"):
                fail(f'{workflow.get("analysis_type")} contains incomplete file guidance')
    for case in qualification_cases:
        if not case.get("finding") or not case.get("biological_context") or not case.get("remaining_limit"):
            fail(f'{case.get("analysis_type")} qualification narrative is incomplete')
        for link in case.get("source_links", []):
            if not re.fullmatch(r"https://[^\s]+", link.get("url", "")):
                fail(f'{case.get("analysis_type")} contains an unsafe source link')

    roadmap = data.get("publication_roadmap", {})
    if roadmap.get("current_phase") != "local_hardening" or len(roadmap.get("phases", [])) != 4:
        fail("publication roadmap is missing or drifted")
    if not any(item.get("gate_id") == "publication_approval" and item.get("state") == "blocked"
               for item in roadmap.get("gates", [])):
        fail("publication approval boundary is not visible")
    if len(data.get("reviewer_files", [])) != 4:
        fail("reviewer pack is incomplete")

    for case in cases:
        if not case.get("input_sha256") or not case.get("evidence_files"):
            fail(f'{case.get("case_id")} lacks checksum-bound evidence')
        for value in case["input_sha256"].values():
            if not re.fullmatch(r"[0-9a-f]{64}", value):
                fail(f'{case["case_id"]} contains an invalid input digest')
        for record in case["evidence_files"]:
            relative = record.get("path", "")
            path = (PUBLIC / relative).resolve()
            try:
                path.relative_to(PUBLIC.resolve())
            except ValueError:
                fail(f"unsafe evidence path: {relative}")
            if path.suffix.lower() not in {".json", ".tsv"} or not path.is_file():
                fail(f"missing or unsupported evidence file: {relative}")
            if digest(path) != record.get("sha256"):
                fail(f"public evidence checksum mismatch: {relative}")

    manifest = json.loads((PUBLIC / "PUBLIC_SHOWCASE_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("workflow_count") != 6 or manifest.get("executed_case_count") != 6:
        fail("manifest counts are incorrect")
    if manifest.get("pending_public_cases") != []:
        fail("manifest incorrectly marks a generated public case pending")
    if manifest.get("network_dependencies") != [] or manifest.get("external_uploads") != []:
        fail("public package declares a network or upload dependency")
    expected_sources = {
        "SHOWCASE.json": digest(TECHNICAL / "SHOWCASE.json"),
        "SF_CSA_CASE.json": digest(SFCSA / "CASE.json"),
        "BASELINE.json": digest(ROOT / "yauvi-structural-workbench" / "BASELINE.json"),
        "RELEASE_STATUS.json": digest(ROOT / "yauvi-structural-workbench" / "RELEASE_STATUS.json"),
        "PLATFORM_IDENTITY.json": digest(ROOT / "yauvi-structural-workbench" / "PLATFORM_IDENTITY.json"),
        "START_HERE.md": digest(ROOT / "yauvi-structural-workbench" / "START_HERE.md"),
        "JOSS_PUBLICATION_ROADMAP.json": digest(ROOT / "yauvi-structural-workbench" / "JOSS_PUBLICATION_ROADMAP.json"),
        "QUALIFICATION_RESULTS.json": digest(ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v1" / "results" / "QUALIFICATION_RESULTS.json"),
        "SOURCE_VERIFICATION.json": digest(ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v1" / "results" / "SOURCE_VERIFICATION.json"),
        "SOURCE_LOCK.json": digest(ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v1" / "SOURCE_LOCK.json"),
        "QUALIFICATION_V2_STATUS.json": digest(ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2" / "results" / "QUALIFICATION_V2_STATUS.json"),
        "QUALIFICATION_V2_PANEL_MANIFEST.json": digest(ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2" / "PANEL_MANIFEST.json"),
    }
    for name, value in expected_sources.items():
        if manifest.get("source_sha256", {}).get(name) != value:
            fail(f"source digest drifted: {name}")
    for relative, expected in manifest.get("generated_sha256", {}).items():
        path = PUBLIC / relative
        if not path.is_file() or digest(path) != expected:
            fail(f"generated digest drifted: {relative}")

    checksums = json.loads((PUBLIC / "CHECKSUMS.json").read_text(encoding="utf-8"))
    expected_files = {
        path.relative_to(PUBLIC).as_posix() for path in PUBLIC.rglob("*")
        if path.is_file() and path.name != "CHECKSUMS.json"
    }
    if set(checksums) != expected_files:
        fail("CHECKSUMS.json does not cover the complete public package")
    for relative, expected in checksums.items():
        if digest(PUBLIC / relative) != expected:
            fail(f"package checksum mismatch: {relative}")

    index = (PUBLIC / "index.html").read_text(encoding="utf-8")
    for phrase in (
        "From protein coordinates to", "Explore six traceable cases", "Start Mark 1 locally",
        "Inspect raw evidence", "Print / Save as PDF", "Six scientific questions",
        "Six traceable evidence cases", "Independent public qualification",
        "What agreed—and what did not", "JOSS publication path",
        "YAUVI Structural Biology Platform — Mark 1", "Share and start",
        "share/PLATFORM_IDENTITY.json", "share/START_HERE.md",
    ):
        if phrase not in index:
            fail(f"public narrative is missing: {phrase}")
    if re.search(r"(?:src|href)=[\"']https?://", index, flags=re.IGNORECASE):
        fail("index.html depends on a remote asset or link")
    app = (PUBLIC / "app.js").read_text(encoding="utf-8")
    if 'fetch("/api/structural-tools"' not in app or "location.protocol" not in app or "location.hostname" not in app:
        fail("same-origin loopback detection is incomplete")
    if re.search(r"fetch\(\s*[`\"']https?://", app):
        fail("app.js attempts an external fetch")
    for phrase in (
        "Synthetic case passed", "Pipeline passed · engines stubbed",
        "Research benefits and scientific limits", "Open raw-file inventory",
        "Public case passed", "Public case partial", "Files you provide",
    ):
        if phrase not in app:
            fail(f"progressive-disclosure UX is missing: {phrase}")

    forbidden = re.compile(
        r"(?:/Users/|/private/var/|file://|20\d\d-\d\d-\d\dT\d\d:|"
        r"YAUVI-PeriodontalPathogens|OralomeVax|K9-PeriodontalVax|unpublished[_ -]sequence)",
        flags=re.IGNORECASE,
    )
    for path in PUBLIC.rglob("*"):
        if path.is_file() and path.suffix.lower() in TEXT_SUFFIXES:
            if forbidden.search(path.read_text(encoding="utf-8", errors="replace")):
                fail(f"private path, timestamp, or campaign identifier leaked into {path.relative_to(PUBLIC)}")
    print(f"verified public evidence showcase: 6 workflows, 5 synthetic analysis cases, "
          f"1 stubbed SF-CSA pipeline case, 4 passed and 2 partial public qualification cases, "
          f"{baseline['total_passed']} baseline tests reported")
    return 0


if __name__ == "__main__":
    sys.exit(main())
