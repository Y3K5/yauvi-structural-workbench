#!/usr/bin/env python3
"""Execute an adopted Qualification v2 stratum and evaluate its gates.

The companion `run_qualification.py` audits panel *composition* only. This
runner executes the canonical engine for each adopted case and compares the
result against the expectation the case declares.

Scope discipline, deliberately narrow:

* Executing a stratum qualifies nothing. A workflow scope becomes qualified
  only when every stratum in its panel passes and reproduces on a second
  machine. This runner reports one stratum and says so.
* Gate semantics come from the panel's declared `gate_semantics`, not from this
  file. Choosing a mode is an adoption decision recorded in data, so the runner
  cannot quietly reinterpret a frozen gate.
* Network access is forbidden during execution. Artifacts must already be
  acquired and checksum-verified.

Usage:  python run_execution.py [--panel ADOPTION_DRAFT_XRAY.json]
Exit:   0 every case passed, 1 a case failed or the panel could not be executed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def close(a: Any, b: Any, tol: float) -> bool:
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return a == b


def invoke(inv: Mapping[str, Any], exe: str, out: Path):
    """Run the engine for one case. The single place a case is executed."""
    if out.exists():
        shutil.rmtree(out)
    cmd = [exe, "run", "--structure", str(HERE / inv["structure"]),
           "--reference-fasta", str(HERE / inv["reference_fasta"]),
           "--validation-report", str(HERE / inv["validation_report"]),
           "--out", str(out)]
    if inv.get("provenance"):
        cmd += ["--provenance", str(HERE / inv["provenance"])]
    if inv.get("chain"):
        cmd += ["--chain", inv["chain"]]
    return subprocess.run(cmd, capture_output=True, text=True), out


def measure(record: Mapping[str, Any], exe: str, out_root: Path) -> dict[str, Any]:
    """Emit an expectation from an actual run.

    Expectations are recorded by executing, never transcribed by hand. An
    earlier hand-recorded expectation was wrong because the command used to read
    it was piped, so the shell reported the pipe's exit status instead of the
    engine's.
    """
    inv = record["expected_result"]["invocation"]
    proc, out = invoke(inv, exe, out_root / record["record_id"])
    ev = json.loads((out / "STRUCTURE_EVIDENCE.json").read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
    c, cs = ev["completeness"], ev["chain_summaries"][0]
    missing = manifest.get("missing_evidence") or []
    return {
        "cli_exit_code": proc.returncode,
        "exit_code_reason": ("complete: every declared evidence leg was supplied"
                             if not missing else
                             f"scientifically incomplete, missing: {', '.join(missing)}"),
        "invocation": inv,
        "residue_identity": {k: c[k] for k in
            ("identity_fraction", "coverage_fraction", "mapped_residues",
             "reference_length", "coordinate_residues", "state")},
        "chain_summary": {k: cs[k] for k in
            ("chain_breaks", "missing_backbone_residues", "residues")},
        "official_metric_import": {"state": ev["external_validation"]["state"],
                                   "values": ev["external_validation"]["metrics"]},
        "gemmi_coordinate_validation": ev["coordinate"]["parser"]["gemmi_validation"],
        "missing_evidence": missing,
    }


def run_case(record: Mapping[str, Any], exe: str, out_root: Path) -> dict[str, Any]:
    expected = record["expected_result"]
    inv = expected["invocation"]
    checks: list[dict[str, Any]] = []

    artifact = HERE / record["artifact"]
    observed_digest = sha256(artifact) if artifact.is_file() else None
    checks.append({"check": "source_artifact_checksum", "required": True,
                   "expected": record["checksum"], "observed": observed_digest,
                   "passed": observed_digest == record["checksum"]})
    if observed_digest != record["checksum"]:
        return {"record_id": record["record_id"], "passed": False, "checks": checks,
                "note": "artifact missing or checksum-mismatched; engine not invoked"}

    proc, out = invoke(inv, exe, out_root / record["record_id"])

    checks.append({"check": "cli_exit_code", "required": True,
                   "expected": expected["cli_exit_code"], "observed": proc.returncode,
                   "passed": proc.returncode == expected["cli_exit_code"]})

    evidence_path = out / "STRUCTURE_EVIDENCE.json"
    if not evidence_path.is_file():
        checks.append({"check": "evidence_written", "required": True, "expected": "STRUCTURE_EVIDENCE.json",
                       "observed": None, "passed": False, "stderr": proc.stderr[-400:]})
        return {"record_id": record["record_id"], "passed": False, "checks": checks}
    ev = json.loads(evidence_path.read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))

    checks.append({"check": "gemmi_coordinate_validation", "required": True,
                   "expected": expected["gemmi_coordinate_validation"],
                   "observed": ev["coordinate"]["parser"]["gemmi_validation"],
                   "passed": ev["coordinate"]["parser"]["gemmi_validation"] == expected["gemmi_coordinate_validation"]})

    # --- gate: residue identity -------------------------------------------
    gates = record["_gate_semantics"]
    ri_mode = gates["residue_identity"]["mode"]
    tol = float(gates["residue_identity"].get("tolerance", 1e-6))
    obs_c, exp_c = ev["completeness"], expected["residue_identity"]
    if ri_mode == "equals_one":
        ok = close(obs_c["identity_fraction"], 1.0, tol)
        checks.append({"check": f"residue_identity[{ri_mode}]", "required": True, "expected": 1.0,
                       "observed": obs_c["identity_fraction"], "passed": ok})
    else:
        ok = close(obs_c["identity_fraction"], exp_c["identity_fraction"], tol)
        checks.append({"check": f"residue_identity[{ri_mode}]", "required": True,
                       "expected": exp_c["identity_fraction"], "observed": obs_c["identity_fraction"], "passed": ok})
    for field in ("coverage_fraction", "mapped_residues", "reference_length"):
        checks.append({"check": f"completeness.{field}", "required": True,
                       "expected": exp_c[field], "observed": obs_c[field],
                       "passed": close(obs_c[field], exp_c[field], tol)})

    # --- gate: official metric import -------------------------------------
    exp_m = expected["official_metric_import"]
    obs_ev = ev["external_validation"]
    checks.append({"check": "official_validation_state", "required": True,
                   "expected": exp_m["state"], "observed": obs_ev.get("state"),
                   "passed": obs_ev.get("state") == exp_m["state"]})
    obs_metrics = obs_ev.get("metrics") or {}
    for name in gates["official_metric_import"]["required_metrics"]:
        exp_v = (exp_m.get("values") or {}).get(name)
        checks.append({"check": f"metric.{name}", "required": True, "expected": exp_v,
                       "observed": obs_metrics.get(name),
                       "passed": name in obs_metrics and close(obs_metrics[name], exp_v, tol)})

    # --- gate: fail-closed behaviour --------------------------------------
    missing = manifest.get("missing_evidence") or []
    fail_closed_ok = (proc.returncode == 0) if not missing else (proc.returncode == 1)
    checks.append({"check": "missing_evidence_fail_closed", "required": True,
                   "expected": "exit 1 when evidence is missing, 0 otherwise",
                   "observed": {"exit_code": proc.returncode, "missing_evidence": missing},
                   "passed": fail_closed_ok})

    return {"record_id": record["record_id"], "pdb_entry_id": record["pdb_entry_id"],
            "stratum": record["stratum"], "split": record["split"],
            "passed": all(c["passed"] for c in checks if c["required"]), "checks": checks}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, default=HERE / "ADOPTION_DRAFT_XRAY.json")
    ap.add_argument("--out", type=Path, default=RESULTS / "execution")
    ap.add_argument("--record", action="store_true",
                    help="Re-measure every case and write the expectations back into the panel.")
    args = ap.parse_args(argv)

    exe = shutil.which("structqc")
    if exe is None:
        print("structqc is not on PATH; install the distribution first", file=sys.stderr)
        return 1
    panel = json.loads(args.panel.read_text(encoding="utf-8"))
    semantics = panel.get("gate_semantics")
    if not semantics:
        print(f"{args.panel.name} declares no gate_semantics; refusing to assume any", file=sys.stderr)
        return 1

    args.out.mkdir(parents=True, exist_ok=True)

    if args.record:
        for record in panel["records"] + panel.get("controls", []):
            record["expected_result"] = measure(record, exe, args.out)
            print(f"  recorded {record['record_id']}: "
                  f"exit={record['expected_result']['cli_exit_code']} "
                  f"identity={record['expected_result']['residue_identity']['identity_fraction']}")
        args.panel.write_text(json.dumps(panel, indent=2) + "\n", encoding="utf-8")
        print(f"wrote measured expectations into {args.panel.name}")
        return 0

    cases = []
    for record in panel["records"]:
        cases.append(run_case({**record, "_gate_semantics": semantics}, exe, args.out))
    # Controls are not counted toward the panel's required case counts. They
    # exist so a gate that would otherwise be trivially satisfied is exercised
    # against a deliberately incomplete run.
    controls = [run_case({**r, "_gate_semantics": semantics}, exe, args.out)
                for r in panel.get("controls", [])]

    passed = sum(1 for c in cases if c["passed"])
    controls_passed = sum(1 for c in controls if c["passed"])
    result = {
        "schema_version": "1.0",
        "collection_id": panel.get("panel_id"),
        "stratum": panel.get("stratum"),
        "gate_semantics": semantics,
        "scientific_execution_performed": True,
        "case_counts": {"total": len(cases), "passed": passed, "failed": len(cases) - passed},
        "stratum_state": ("passed" if passed == len(cases) and controls_passed == len(controls)
                          else "failed"),
        "control_counts": {"total": len(controls), "passed": controls_passed},
        "controls": controls,
        "scope_qualified": False,
        "scope_qualification_note": (
            "Executing one stratum qualifies no workflow scope. The StructQC panel also requires the "
            "cryo_em, nmr, and alphafold strata, and every Mark 1 scope additionally requires "
            "independent reproduction on a second machine."),
        "runtime": {"python": platform.python_version(), "machine": platform.machine(),
                    "platform": platform.system()},
        "cases": cases,
    }
    for name, value in (("EXECUTION_STATUS.json", result),):
        path = RESULTS / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(canonical(value))
    print(f"{panel.get('stratum')} stratum: {passed}/{len(cases)} cases passed, "
          f"{controls_passed}/{len(controls)} controls passed "
          f"[residue_identity={semantics['residue_identity']['mode']}]")
    for c in cases + controls:
        if not c["passed"]:
            bad = [k["check"] for k in c["checks"] if k["required"] and not k["passed"]]
            print(f"  FAILED {c['record_id']}: {', '.join(bad)}")
    print("scope_qualified: false (one stratum of four; second-machine reproduction outstanding)")
    return 0 if passed == len(cases) and controls_passed == len(controls) else 1


if __name__ == "__main__":
    raise SystemExit(main())
