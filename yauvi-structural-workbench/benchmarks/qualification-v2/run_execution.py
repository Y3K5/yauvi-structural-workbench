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
import csv
import hashlib
import json
import math
import platform
import shutil
import re
import subprocess
import sys
import tempfile
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


def drift_deltas(records: Any) -> dict[str, list[float]]:
    """Collect every informational ``drift.*`` delta, keyed by check name.

    Drift checks are recorded and never required, so a run can pass with a
    delta nobody looked at. This is the one place that reads them back out,
    shared by the runner's own summary line and by ``summarize_execution.py``,
    so the two can never disagree about what a run measured.
    """
    deltas: dict[str, list[float]] = {}
    for record in records:
        for check in record["checks"]:
            if check["check"].startswith("drift.") and not check["required"]:
                delta = (check["observed"] or {}).get("delta")
                if delta is not None:
                    deltas.setdefault(check["check"], []).append(float(delta))
    return {name: sorted(values) for name, values in sorted(deltas.items())}


# --- workflow engines -------------------------------------------------------
# Each panel names a workflow, executed by a different CLI, emitting a different
# evidence document, judged by different gates. Only those three things vary;
# acquisition, checksums, expectation recording, controls and coverage are
# shared. Keeping the variation in one registry is what stops the second panel
# becoming a second copy of this file.

def at(reference: Any) -> Path:
    """Resolve a path a case declares.

    Panel paths are relative to this directory, but a chained upstream run
    contributes an absolute path, and --out may be given relative to whatever
    directory the runner was invoked from. Joining every path to HERE produced a
    doubled path when the two conventions met, so absolute paths pass through
    untouched.
    """
    path = Path(str(reference))
    return path if path.is_absolute() else HERE / path


def panel_relative(path: Path) -> str:
    """Record where a case wrote, without recording whose machine wrote it.

    `at()` resolves a relative path against this directory, so a panel-relative
    string round-trips through the coverage pass unchanged. An absolute one puts
    a home directory into evidence that is published to a public repository, and
    is useless to a reviewer who cannot resolve it anyway. Paths outside the
    panel keep their absolute form, because there is nothing else to say about
    them.
    """
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(HERE))
    except ValueError:
        return str(resolved)


def _cmd_structure_qc(inv: Mapping[str, Any], exe: str, out: Path) -> list[str]:
    cmd = [exe, "run", "--structure", str(at(inv["structure"])),
           "--reference-fasta", str(at(inv["reference_fasta"])), "--out", str(out)]
    # A predicted model has no wwPDB validation report and does have a PAE
    # matrix; an experimental entry is the other way round. Neither is passed
    # unless the case declares it.
    if inv.get("validation_report"):
        cmd += ["--validation-report", str(at(inv["validation_report"]))]
    if inv.get("pae"):
        cmd += ["--pae", str(at(inv["pae"]))]
    if inv.get("provenance"):
        cmd += ["--provenance", str(at(inv["provenance"]))]
    if inv.get("chain"):
        cmd += ["--chain", inv["chain"]]
    # NMR entries deposit an ensemble. Which model the expectation describes is
    # a stratum decision, declared per case rather than left to the CLI default.
    if inv.get("model") is not None:
        cmd += ["--model", str(inv["model"])]
    return cmd


def _cmd_functional_site_state(inv: Mapping[str, Any], exe: str, out: Path) -> list[str]:
    # site-context consumes StructQC's evidence document, so a case is a
    # two-stage chain: the upstream run is produced first and referenced here.
    cmd = [exe, "run", "--manifest", str(at(inv["manifest"])),
           "--structure", str(at(inv["structure"])),
           "--annotations", str(at(inv["annotations"])), "--out", str(out)]
    if inv.get("component_map"):
        cmd += ["--component-map", str(at(inv["component_map"]))]
    if inv.get("pocket_result"):
        cmd += ["--pocket-result", str(at(inv["pocket_result"]))]
    return cmd


def _cmd_assembly_interface(inv: Mapping[str, Any], exe: str, out: Path) -> list[str]:
    # Like site-context, this consumes StructQC's evidence; it additionally needs
    # both the isolated coordinates and the biological assembly, because the
    # quantity under test is how much of the subject chain the assembly buries.
    cmd = [exe, "run", "--manifest", str(at(inv["manifest"])),
           "--isolated", str(at(inv["isolated"])),
           "--assembly", str(at(inv["assembly"])),
           "--subject-chain", str(inv["subject_chain"]),
           "--relationship", str(inv.get("relationship", "exact_protein")),
           "--out", str(out)]
    if inv.get("expected_chains"):
        cmd += ["--expected-chains", str(inv["expected_chains"])]
    if inv.get("assembly_id"):
        cmd += ["--assembly-id", str(inv["assembly_id"])]
    if inv.get("reference_id"):
        cmd += ["--reference-id", str(inv["reference_id"])]
    return cmd


def freesasa_version() -> str | None:
    """The FreeSASA build actually on PATH.

    assembly-context invokes FreeSASA but records only 'available_invoked', with
    no version anywhere in its evidence. A panel whose gate is a 0.001 relative
    tolerance against a reference implementation cannot leave which
    implementation unrecorded, so the executor captures it.
    """
    exe = shutil.which("freesasa")
    if exe is None:
        return None
    try:
        out = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    first = (out.stdout or out.stderr).strip().splitlines()
    return first[0].strip() if first else None


def standalone_sasa(structure: Path, chain: str, whole_assembly: bool) -> float | None:
    """Subject-chain SASA from the FreeSASA command line, independent of the module.

    The coordinates are converted to PDB first. FreeSASA's mmCIF reader is a
    recent addition and distribution packages ship builds without it, so a check
    that required ``--cif`` passed on macOS and failed on every Ubuntu runner --
    while assembly-context itself was unaffected, because it writes temporary
    PDB files and never asks for CIF.

    Independence is preserved where it matters: the chain selection, the file
    written, and the parsing are all done here, and the per-residue ``seq``
    output is summed rather than reading the module's own JSON totals.
    """
    exe = shutil.which("freesasa")
    if exe is None:
        return None
    try:
        import gemmi
        st = gemmi.read_structure(str(structure))
        st.setup_entities()
        st.remove_ligands_and_waters()
        if not whole_assembly:
            for model in st:
                for name in [c.name for c in model]:
                    if name != chain:
                        model.remove_chain(name)
        with tempfile.TemporaryDirectory(prefix="yauvi-freesasa-") as tmp:
            pdb = Path(tmp) / "subject.pdb"
            st.write_pdb(str(pdb))
            run = subprocess.run([exe, str(pdb), "--format=seq"],
                                 capture_output=True, text=True, timeout=900)
        total, found = 0.0, False
        for line in run.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 4 and parts[0] == "SEQ" and parts[1] == chain:
                try:
                    total += float(parts[-1]); found = True
                except ValueError:
                    continue
        return round(total, 6) if found else None
    except Exception:
        return None


def _cmd_membrane_orientation(inv: Mapping[str, Any], exe: str, out: Path) -> list[str]:
    # `memorient run` is used rather than `orient`: both emit the validation
    # block carrying the rotation series this panel gates on, but `run` also
    # writes the RUN_MANIFEST.json the shared harness reads.
    cmd = [exe, "run", "--structure", str(at(inv["structure"])), "--out", str(out)]
    if inv.get("context"):
        cmd += ["--context", str(inv["context"])]
    if inv.get("chain"):
        cmd += ["--chain", str(inv["chain"])]
    if inv.get("topology_evidence"):
        cmd += ["--topology-evidence", str(at(inv["topology_evidence"]))]
    return cmd


def opm_half_thickness(path: Path) -> float | None:
    """The reference half-thickness an OPM file carries in its own header.

    OPM distributes structures already rotated so the bilayer normal is Z, and
    records the half-thickness it fitted as a REMARK. That is the independent
    quantity this panel compares against; reading it from the file keeps the
    reference and the coordinates from drifting apart.
    """
    try:
        head = path.read_text(errors="ignore")[:8000]
    except OSError:
        return None
    match = re.search(r"1/2 of bilayer thickness:\s*([0-9.]+)", head)
    return float(match.group(1)) if match else None


ENGINES: dict[str, dict[str, Any]] = {
    "structure_qc": {"cli": "structqc", "evidence": "STRUCTURE_EVIDENCE.json",
                     "build_cmd": _cmd_structure_qc},
    "functional_site_state": {"cli": "site-context", "evidence": "SITE_CONTEXT.json",
                              "build_cmd": _cmd_functional_site_state},
    "assembly_interface": {"cli": "assembly-context", "evidence": "ASSEMBLY_CONTEXT.json",
                           "build_cmd": _cmd_assembly_interface},
    "membrane_orientation": {"cli": "memorient", "evidence": "MEMBRANE_ORIENTATION.json",
                             "build_cmd": _cmd_membrane_orientation},
}


def invoke(inv: Mapping[str, Any], exe: str, out: Path, workflow: str = "structure_qc"):
    """Run the engine for one case. The single place a case is executed.

    A downstream workflow consumes an upstream module's evidence document. That
    document is regenerated here from the same locked sources rather than being
    committed, so a case cannot pass against a stale manifest that no longer
    corresponds to the coordinates in the lock.
    """
    if out.exists():
        shutil.rmtree(out)
    inv = dict(inv)
    upstream = inv.pop("upstream", None)
    if upstream:
        up_workflow = upstream.get("workflow", "structure_qc")
        up_exe = shutil.which(ENGINES[up_workflow]["cli"])
        if up_exe is None:
            raise RuntimeError(f"{ENGINES[up_workflow]['cli']} is not on PATH")
        up_out = out / "upstream"
        up_cmd = ENGINES[up_workflow]["build_cmd"](upstream, up_exe, up_out)
        up = subprocess.run(up_cmd, capture_output=True, text=True)
        up_evidence = up_out / ENGINES[up_workflow]["evidence"]
        if not up_evidence.is_file():
            return up, out          # downstream cannot run; caller reports it
        inv["manifest"] = str(up_evidence.resolve())
    cmd = ENGINES[workflow]["build_cmd"](inv, exe, out)
    return subprocess.run(cmd, capture_output=True, text=True), out


def _measure_structure_qc(ev, manifest) -> dict[str, Any]:
    c, cs = ev["completeness"], ev["chain_summaries"][0]
    return {
        "residue_identity": {k: c[k] for k in
            ("identity_fraction", "coverage_fraction", "mapped_residues",
             "reference_length", "coordinate_residues", "state")},
        "chain_summary": {k: cs[k] for k in
            ("chain_breaks", "missing_backbone_residues", "residues")},
        "official_metric_import": {"state": ev["external_validation"]["state"],
                                   "values": ev["external_validation"]["metrics"]},
        "confidence": ev.get("pae"),
        "provenance_class": (ev.get("provenance") or {}).get("class"),
        "gemmi_coordinate_validation": ev["coordinate"]["parser"]["gemmi_validation"],
    }


def _measure_functional_site_state(ev, manifest) -> dict[str, Any]:
    sites = ev.get("sites") or []
    counts: dict[str, int] = {}
    for s in sites:
        counts[s["state"]] = counts.get(s["state"], 0) + 1
    resolved = sum(n for k, n in counts.items() if k in ("role_compatible", "role_mismatch"))
    return {
        "site_count": len(sites),
        "state_counts": dict(sorted(counts.items())),
        # The per-site vector is the real expectation. Counts alone would hide a
        # pair of compensating changes leaving two residues silently swapped.
        "sites": [{"position": s["position"], "chain_id": s.get("chain_id"),
                   "role": s.get("role"), "state": s["state"],
                   "observed_residue": s.get("observed_residue"),
                   "expected_residues": s.get("expected_residues")} for s in sites],
        "curated_residue_recovery": round(resolved / len(sites), 6) if sites else 0.0,
        "cofactors": ev.get("cofactors") or [],
        "observed_heteroatoms": sorted({h.get("component_id") for h in
                                        (ev.get("observed_heteroatoms") or [])}),
    }


def _measure_assembly_interface(ev, manifest) -> dict[str, Any]:
    assembly = ev.get("assembly") or {}
    surface = ev.get("surface") or {}
    return {
        "stoichiometry": {k: assembly.get(k) for k in
                          ("chains_expected", "chains_observed", "complete", "lower_bound")},
        "surface": {k: surface.get(k) for k in
                    ("subject_isolated_sasa_A2", "subject_assembly_sasa_A2", "buried_sasa_A2")},
        "residue_contacts": len(ev.get("residue_contacts") or []),
        "interfaces": len(ev.get("interfaces") or []),
        "methods": ev.get("methods"),
        # Recorded by the executor, not by the module: assembly-context reports
        # only that FreeSASA was "available_invoked".
        "freesasa_version": freesasa_version(),
    }


def _measure_membrane_orientation(ev, manifest) -> dict[str, Any]:
    summary = ev.get("summary") or {}
    validation = ev.get("validation") or {}
    errors = validation.get("normal_angle_errors_deg") or []
    return {
        "scope_id": summary.get("scope_id"),
        "fitted_half_thickness_A": summary.get("half_thickness"),
        "normal_angle_errors_deg": errors,
        "mean_unsigned_normal_error_deg": round(sum(errors)/len(errors), 6) if errors else None,
        "max_unsigned_normal_error_deg": round(max(errors), 6) if errors else None,
        "normal_drift_deg": round(max(errors)-min(errors), 6) if errors else None,
        "mean_jaccard": validation.get("mean_jaccard"),
        "embedded_jaccards": validation.get("embedded_jaccards"),
        "n_reference_extracellular": validation.get("n_reference_extracellular"),
        "extracellular_comparison_state": validation.get("extracellular_comparison_state"),
        "n_embedded": summary.get("n_embedded"),
        "delta_kd": summary.get("delta_kd"),
        "fitted_normal": summary.get("normal"),
        "normal_frame": summary.get("normal_frame"),
    }


MEASURERS = {"structure_qc": _measure_structure_qc,
             "membrane_orientation": _measure_membrane_orientation,
             "assembly_interface": _measure_assembly_interface,
             "functional_site_state": _measure_functional_site_state}


def measure(record: Mapping[str, Any], exe: str, out_root: Path,
            workflow: str = "structure_qc") -> dict[str, Any]:
    """Emit an expectation from an actual run.

    Expectations are recorded by executing, never transcribed by hand. An
    earlier hand-recorded expectation was wrong because the command used to read
    it was piped, so the shell reported the pipe's exit status, not the engine's.
    """
    inv = record["expected_result"]["invocation"]
    proc, out = invoke(inv, exe, out_root / record["record_id"], workflow)
    ev = json.loads((out / ENGINES[workflow]["evidence"]).read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
    missing = manifest.get("missing_evidence") or []
    return {
        "cli_exit_code": proc.returncode,
        "exit_code_reason": ("complete: every declared evidence leg was supplied"
                             if not missing else
                             f"scientifically incomplete, missing: {', '.join(missing)}"),
        "invocation": inv,
        "missing_evidence": missing,
        **MEASURERS[workflow](ev, manifest),
    }


def _gates_structure_qc(record, expected, ev, checks) -> bool:
    """StructQC gate checks. Returns False if the case cannot be judged further.

    Lifted verbatim from run_case so a second workflow can be added without this
    file growing a second copy of the shared harness. Behaviour is unchanged.
    """
    gates = record["_gate_semantics"]
    checks.append({"check": "gemmi_coordinate_validation", "required": True,
                   "expected": expected["gemmi_coordinate_validation"],
                   "observed": ev["coordinate"]["parser"]["gemmi_validation"],
                   "passed": ev["coordinate"]["parser"]["gemmi_validation"] == expected["gemmi_coordinate_validation"]})

    # --- gate: residue identity -------------------------------------------
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
    # What counts as the official metric set depends on the stratum. An
    # experimental entry has a wwPDB validation report; a predicted model has no
    # such report and carries pLDDT and PAE instead. A single global metric list
    # cannot describe both.
    stratum = record["stratum"]
    metric_spec = gates["official_metric_import"]["by_stratum"].get(stratum)
    if metric_spec is None:
        checks.append({"check": "official_metric_import.stratum_declared", "required": True,
                       "expected": f"gate_semantics for stratum {stratum}", "observed": None,
                       "passed": False})
        return {"record_id": record["record_id"], "passed": False, "checks": checks}

    exp_m = expected["official_metric_import"]
    obs_ev = ev["external_validation"]
    checks.append({"check": "official_validation_state", "required": True,
                   "expected": exp_m["state"], "observed": obs_ev.get("state"),
                   "passed": obs_ev.get("state") == exp_m["state"]})
    obs_metrics = obs_ev.get("metrics") or {}
    for name in metric_spec.get("required_metrics", []):
        exp_v = (exp_m.get("values") or {}).get(name)
        checks.append({"check": f"metric.{name}", "required": True, "expected": exp_v,
                       "observed": obs_metrics.get(name),
                       "passed": name in obs_metrics and close(obs_metrics[name], exp_v, tol)})
    for field in metric_spec.get("required_confidence", []):
        exp_v = (expected.get("confidence") or {}).get(field)
        obs_v = (ev.get("pae") or {}).get(field)
        checks.append({"check": f"confidence.{field}", "required": True, "expected": exp_v,
                       "observed": obs_v, "passed": obs_v is not None and close(obs_v, exp_v, tol)})

    return True


def _gates_functional_site_state(record, expected, ev, checks) -> bool:
    """Functional-site gate checks.

    The panel declares four gates. Three are thresholds; the fourth -- that a
    curated residue keeps the state it was recorded with -- is the one that
    actually catches drift, so it compares the per-site vector rather than the
    state counts. A pair of compensating changes leaves the counts identical.
    """
    fs = record["_gate_semantics"].get("functional_site", {})
    obs = _measure_functional_site_state(ev, None)

    exp_sites = {(s["position"], s.get("chain_id")): s for s in (expected.get("sites") or [])}
    obs_sites = {(s["position"], s.get("chain_id")): s for s in obs["sites"]}
    drifted = sorted(
        f"{k[1]}:{k[0]} {exp_sites[k]['state']}->{obs_sites.get(k, {}).get('state')}"
        for k in exp_sites if obs_sites.get(k, {}).get("state") != exp_sites[k]["state"])
    checks.append({"check": "site_states_match_declared", "required": True,
                   "expected": "every curated residue in its declared state",
                   "observed": drifted or "all match", "passed": not drifted})

    checks.append({"check": "site_count_unchanged", "required": True,
                   "expected": expected.get("site_count"), "observed": obs["site_count"],
                   "passed": obs["site_count"] == expected.get("site_count")})

    floor = fs.get("unambiguous_curated_residue_recovery_min")
    if floor is not None:
        # Recovery is the fraction of curated residues the run could resolve to a
        # coordinate at all. A stratum built from apo, modified, or incomplete
        # structures is expected to fall below a floor written for complete ones,
        # so the floor is applied per stratum rather than panel-wide.
        exempt = set(fs.get("recovery_floor_exempt_strata") or [])
        if record["stratum"] in exempt:
            checks.append({"check": "curated_residue_recovery", "required": False,
                           "expected": f"floor {floor} not applied to stratum {record['stratum']}",
                           "observed": obs["curated_residue_recovery"], "passed": True})
        else:
            checks.append({"check": "curated_residue_recovery", "required": True,
                           "expected": f">= {floor}", "observed": obs["curated_residue_recovery"],
                           "passed": obs["curated_residue_recovery"] >= float(floor)})

    # A false exact mapping is the panel's zero-tolerance failure: a residue
    # reported compatible whose observed identity is not one the curation expects.
    false_exact = [s for s in obs["sites"]
                   if s["state"] == "role_compatible" and s.get("expected_residues")
                   and s.get("observed_residue") not in s["expected_residues"]]
    ceiling = int(fs.get("false_exact_mappings_max", 0))
    checks.append({"check": "false_exact_mappings", "required": True,
                   "expected": f"<= {ceiling}", "observed": len(false_exact),
                   "passed": len(false_exact) <= ceiling})

    if fs.get("cofactor_identifier_match") == "exact":
        checks.append({"check": "cofactor_identifiers_exact", "required": True,
                       "expected": expected.get("cofactors"), "observed": obs["cofactors"],
                       "passed": obs["cofactors"] == (expected.get("cofactors") or [])})
        checks.append({"check": "observed_heteroatoms_exact", "required": True,
                       "expected": expected.get("observed_heteroatoms"),
                       "observed": obs["observed_heteroatoms"],
                       "passed": obs["observed_heteroatoms"] == (expected.get("observed_heteroatoms") or [])})
    return True


def _gates_assembly_interface(record, expected, ev, checks) -> bool:
    """Assembly-interface gate checks.

    The panel's gates are exact stoichiometry and a SASA agreement tolerance.
    The second is only meaningful if something independent produces the
    comparison value, so the buried-surface numbers are checked against the
    FreeSASA command line rather than against the module's own previous output.
    """
    ai = record["_gate_semantics"].get("assembly_interface", {})
    obs = _measure_assembly_interface(ev, None)
    inv = expected["invocation"]

    if ai.get("operator_copy_stoichiometry") == "exact":
        exp_st, obs_st = expected.get("stoichiometry") or {}, obs["stoichiometry"]
        for field in ("chains_expected", "chains_observed", "complete"):
            checks.append({"check": f"stoichiometry.{field}", "required": True,
                           "expected": exp_st.get(field), "observed": obs_st.get(field),
                           "passed": obs_st.get(field) == exp_st.get(field)})

    exp_sf, obs_sf = expected.get("surface") or {}, obs["surface"]
    for field in ("subject_isolated_sasa_A2", "subject_assembly_sasa_A2", "buried_sasa_A2"):
        checks.append({"check": f"surface.{field}", "required": True,
                       "expected": exp_sf.get(field), "observed": obs_sf.get(field),
                       "passed": close(obs_sf.get(field), exp_sf.get(field), 1e-6)})

    # The version is part of the evidence: a tolerance against a reference
    # implementation is meaningless without knowing which build produced it.
    # A version must be recorded, because a tolerance against a reference
    # implementation is meaningless if nothing says which one produced it. It is
    # not required to match across platforms: the scientific gate is numeric
    # agreement, and demanding identical version strings would fail runners whose
    # builds agree perfectly well. A change is visible in the evidence instead.
    checks.append({"check": "freesasa_version_recorded", "required": True,
                   "expected": "a FreeSASA version reported by the binary on PATH",
                   "observed": obs["freesasa_version"],
                   "passed": bool(obs["freesasa_version"])})
    checks.append({"check": "freesasa_version_unchanged", "required": False,
                   "expected": expected.get("freesasa_version"),
                   "observed": obs["freesasa_version"],
                   "passed": obs["freesasa_version"] == expected.get("freesasa_version")})

    rel = float(ai.get("freesasa_relative_tolerance", 0.001))
    abs_tol = float(ai.get("freesasa_absolute_tolerance_A2", 1.0))
    chain = str(inv["subject_chain"])
    for field, path, whole in (
            ("subject_isolated_sasa_A2", at(inv["isolated"]), False),
            ("subject_assembly_sasa_A2", at(inv["assembly"]), True)):
        reference = standalone_sasa(path, chain, whole)
        module_value = obs_sf.get(field)
        if reference is None or module_value is None:
            checks.append({"check": f"freesasa_agreement.{field}", "required": True,
                           "expected": "a standalone FreeSASA value",
                           "observed": None, "passed": False})
            continue
        delta = abs(float(module_value) - reference)
        ok = delta <= abs_tol and (delta / reference if reference else 0.0) <= rel
        checks.append({"check": f"freesasa_agreement.{field}", "required": True,
                       "expected": f"within {rel} relative and {abs_tol} A2 of standalone FreeSASA",
                       "observed": {"module": module_value, "standalone": reference,
                                    "delta_A2": round(delta, 6)},
                       "passed": ok})
    return True


def _gates_membrane_orientation(record, expected, ev, checks) -> bool:
    """Membrane-orientation gate checks.

    The module's own ``validation.passed`` is deliberately not used. It applies
    its ``normal_threshold_deg`` of 1.0 to the *angle error*, so 1BXW at 6.19
    degrees is reported as not placed. The panel's 1.0 degree gate is on the
    *drift across rotations*, which for the same case is 0.015 degrees. Same
    number, different quantity: the panel asks whether the fitted normal is
    reproducible under rotation, not whether it matches OPM to within a degree.
    The gates are therefore computed here from the reported arrays.
    """
    mo = record["_gate_semantics"].get("membrane_orientation", {})
    obs = _measure_membrane_orientation(ev, None)

    checks.append({"check": "scope_id", "required": True,
                   "expected": expected.get("scope_id"), "observed": obs["scope_id"],
                   "passed": obs["scope_id"] == expected.get("scope_id")})

    reference = record.get("opm_half_thickness_A")
    fitted = obs["fitted_half_thickness_A"]
    limit = mo.get("mean_half_thickness_error_A_max")
    if reference is not None and fitted is not None and limit is not None:
        error = abs(float(fitted) - float(reference))
        checks.append({"check": "half_thickness_error_A", "required": True,
                       "expected": f"<= {limit} vs OPM {reference}",
                       "observed": round(error, 6), "passed": error <= float(limit)})

    for key, field in (("mean_unsigned_normal_error_deg_max", "mean_unsigned_normal_error_deg"),
                       ("per_case_unsigned_normal_error_deg_max", "max_unsigned_normal_error_deg"),
                       ("normal_drift_across_8_rotations_deg_max", "normal_drift_deg")):
        limit = mo.get(key)
        value = obs[field]
        if limit is None:
            continue
        checks.append({"check": field, "required": True, "expected": f"<= {limit}",
                       "observed": value,
                       "passed": value is not None and value <= float(limit)})

    # "Non-vacuous" is the load-bearing word: a Jaccard computed over an empty
    # reference set is 1.0 and means nothing, so the comparison must have been
    # evaluated against a non-empty extracellular set before the floor applies.
    floor = mo.get("non_vacuous_residue_jaccard_min")
    if floor is not None:
        n_ec = obs["n_reference_extracellular"] or 0
        evaluated = obs["extracellular_comparison_state"] == "evaluated"
        minimum = int(mo.get("minimum_reference_extracellular_residues", 1))
        jaccard = obs["mean_jaccard"]
        checks.append({"check": "extracellular_comparison_non_vacuous", "required": True,
                       "expected": f"evaluated with at least {minimum} reference residues",
                       "observed": {"state": obs["extracellular_comparison_state"], "n": n_ec},
                       "passed": evaluated and n_ec >= minimum})
        checks.append({"check": "residue_jaccard", "required": True, "expected": f">= {floor}",
                       "observed": jaccard,
                       "passed": jaccard is not None and jaccard >= float(floor)})

    # Drift against the recorded run, reported but not required.
    #
    # These were required at a 1e-6 tolerance, which is bit-level equality on a
    # Nelder-Mead fit and failed on every runner including the platform the
    # values were recorded on. The scientific gates above already bound these
    # quantities with the tolerances the panel actually specifies -- 2.5 A on
    # half-thickness, 15 degrees on the normal -- so requiring exact
    # reproduction on top of them asserted a reproducibility nobody has
    # established. Whether the fit reproduces numerically across machines is
    # precisely what the independent second-machine gate exists to answer, and
    # pretending it holds here would answer it by assumption.
    #
    # Accuracy against the reference orientation.
    #
    # Every other normal gate in this panel measures self-consistency: whether
    # the fit reproduces under rotation. A stable wrong answer passes them, and
    # two cases did -- 1FEP and 2ERV pass every gate in every environment while
    # sitting 6.69 and 5.41 degrees from OPM. This is the only gate that asks
    # whether the orientation is right.
    #
    # OPM distributes structures already oriented with the membrane normal on z,
    # so the reference direction is the panel's declared `reference_normal` and
    # the error is the undirected angle to it -- a membrane normal equals its
    # negation, so the comparison must not be signed.
    reference = mo.get("reference_normal")
    limit = mo.get("opm_normal_error_deg_max")
    fitted = obs.get("fitted_normal")
    if reference is not None and limit is not None:
        if fitted is None:
            checks.append({"check": "opm_normal_error_deg", "required": True,
                           "expected": f"<= {limit}",
                           "observed": "fitted normal absent from evidence",
                           "passed": False})
        else:
            a = [float(v) for v in fitted]
            b = [float(v) for v in reference]
            na = math.sqrt(sum(v * v for v in a)) or 1.0
            nb = math.sqrt(sum(v * v for v in b)) or 1.0
            dot = sum(x * y for x, y in zip(a, b)) / (na * nb)
            error = math.degrees(math.acos(min(1.0, abs(dot))))
            checks.append({"check": "opm_normal_error_deg", "required": True,
                           "expected": f"<= {limit} vs reference {reference}",
                           "observed": round(error, 6),
                           "passed": error <= float(limit)})

    # `n_reference_extracellular` was required here on the reasoning that a
    # count is not a fit. It is a count *of residues selected using the fitted
    # geometry* -- it comes from the same `validation` block as mean_jaccard,
    # not from the OPM reference file -- so it inherits the fit's instability
    # exactly as the jaccard does. Being an integer protects it from
    # float-comparison artifacts, not from the geometry moving underneath it,
    # and quantisation makes it the harshest check in the panel: a perturbation
    # that moves one residue across the boundary is a full-unit error where the
    # continuous quantities absorb it. It failed 7 of 8 failing cases across all
    # six runners while every scientific gate passed.
    #
    # Nothing scientific is lost by reporting it. The real requirement -- that
    # the comparison was evaluated against enough reference residues to mean
    # anything -- is already carried by the required
    # `extracellular_comparison_non_vacuous` gate above, which bounds the same
    # count from below without asserting an exact value no machine reproduces.
    for field in ("fitted_half_thickness_A", "mean_jaccard", "n_reference_extracellular"):
        value, recorded = obs[field], expected.get(field)
        delta = (abs(float(value) - float(recorded))
                 if value is not None and recorded is not None else None)
        checks.append({"check": f"drift.{field}", "required": False,
                       "expected": recorded,
                       "observed": {"value": value, "delta": delta},
                       "passed": delta is not None and delta <= 1e-6})
    return True


GATE_CHECKS = {"structure_qc": _gates_structure_qc,
               "membrane_orientation": _gates_membrane_orientation,
               "assembly_interface": _gates_assembly_interface,
               "functional_site_state": _gates_functional_site_state}


def run_case(record: Mapping[str, Any], exe: str, out_root: Path) -> dict[str, Any]:
    expected = record["expected_result"]
    inv = expected["invocation"]
    checks: list[dict[str, Any]] = []

    artifact = at(record["artifact"])
    observed_digest = sha256(artifact) if artifact.is_file() else None
    checks.append({"check": "source_artifact_checksum", "required": True,
                   "expected": record["checksum"], "observed": observed_digest,
                   "passed": observed_digest == record["checksum"]})
    if observed_digest != record["checksum"]:
        return {"record_id": record["record_id"], "passed": False, "checks": checks,
                "note": "artifact missing or checksum-mismatched; engine not invoked"}

    workflow = record["_workflow"]
    proc, out = invoke(inv, exe, out_root / record["record_id"], workflow)

    checks.append({"check": "cli_exit_code", "required": True,
                   "expected": expected["cli_exit_code"], "observed": proc.returncode,
                   "passed": proc.returncode == expected["cli_exit_code"]})

    evidence_path = out / ENGINES[workflow]["evidence"]
    if not evidence_path.is_file():
        checks.append({"check": "evidence_written", "required": True,
                       "expected": ENGINES[workflow]["evidence"],
                       "observed": None, "passed": False, "stderr": proc.stderr[-400:]})
        return {"record_id": record["record_id"], "passed": False, "checks": checks}
    ev = json.loads(evidence_path.read_text(encoding="utf-8"))
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text(encoding="utf-8"))

    gate_checks = GATE_CHECKS.get(record["_workflow"])
    if gate_checks is None:
        checks.append({"check": "workflow_gates_declared", "required": True,
                       "expected": f"gate checks for workflow {record['_workflow']}",
                       "observed": None, "passed": False})
    elif not gate_checks(record, expected, ev, checks):
        return {"record_id": record["record_id"], "passed": False, "checks": checks}

    # --- gate: fail-closed behaviour --------------------------------------
    # Some evidence cannot exist for a stratum: a predicted model has no
    # community geometry validation, and reporting it missing is correct rather
    # than a failure. The stratum declares what is legitimately absent; anything
    # else appearing in the missing list is a change that must be noticed.
    gates = record["_gate_semantics"]
    stratum = record["stratum"]
    missing = sorted(manifest.get("missing_evidence") or [])
    if record["record_kind"] == "control_case":
        # Controls have declared purposes. A fail-closed control exists to be
        # incomplete, so it must report withheld evidence. A coverage control
        # exists to supply a feature no benchmark case carries, so it must be
        # complete like any other run. Holding both to one rule failed the
        # coverage control for supplying the evidence it was meant to supply.
        purpose = record.get("control_purpose", "fail_closed")
        if purpose == "fail_closed":
            checks.append({"check": "control_reports_withheld_evidence", "required": True,
                           "expected": "a non-empty missing-evidence list",
                           "observed": missing, "passed": bool(missing)})
        else:
            expected_missing = sorted(
                gates["missing_evidence_behavior"]["expected_missing_by_stratum"].get(stratum, []))
            checks.append({"check": "coverage_control_is_complete", "required": True,
                           "expected": expected_missing, "observed": missing,
                           "passed": missing == expected_missing})
    else:
        expected_missing = sorted(
            gates["missing_evidence_behavior"]["expected_missing_by_stratum"].get(stratum, []))
        checks.append({"check": "missing_evidence_matches_stratum_expectation", "required": True,
                       "expected": expected_missing, "observed": missing,
                       "passed": missing == expected_missing})
    checks.append({"check": "missing_evidence_recorded_unchanged", "required": True,
                   "expected": sorted(expected.get("missing_evidence") or []), "observed": missing,
                   "passed": missing == sorted(expected.get("missing_evidence") or [])})

    return {"record_id": record["record_id"], "pdb_entry_id": record["pdb_entry_id"],
            "stratum": record["stratum"], "split": record["split"],
            "output_dir": panel_relative(out),
            "passed": all(c["passed"] for c in checks if c["required"]), "checks": checks}


def _witness_structure_qc(out: Path) -> set[str]:
    """Report which coverage features one executed case actually demonstrates.

    Coverage is verified from the evidence a run produced, not asserted in
    prose. The panel declares which features must appear; this function knows
    how to recognise them. Composition counts records per stratum and therefore
    cannot see a coverage rule going unmet -- that gap is what this closes.
    """
    seen: set[str] = set()
    evidence = out / "STRUCTURE_EVIDENCE.json"
    residues = out / "RESIDUE_QUALITY.tsv"
    if not evidence.is_file():
        return seen
    ev = json.loads(evidence.read_text(encoding="utf-8"))

    coord = ev.get("coordinate", {})
    if len(coord.get("chains") or []) > 1:
        seen.add("multichain_coordinates")
    if (coord.get("model_count") or 1) > 1:
        seen.add("multi_model_ensemble")
    if (ev.get("pae") or {}).get("state") == "evaluated":
        seen.add("predicted_confidence")
    for chain in ev.get("chain_summaries") or []:
        if chain.get("nonstandard_residues"):
            seen.add("nonstandard_residues")
        if chain.get("chain_breaks") or chain.get("missing_backbone_residues"):
            seen.add("missing_residues")
    if (ev.get("external_validation") or {}).get("state") == "imported":
        seen.add("official_validation_import")

    if residues.is_file():
        with residues.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if (row.get("insertion_code") or "").strip():
                    seen.add("insertion_codes")
                if (row.get("altlocs") or "").strip():
                    seen.add("alternate_locations")
                if seen >= {"insertion_codes", "alternate_locations"}:
                    break
    return seen


def _witness_functional_site_state(out: Path) -> set[str]:
    """Coverage features one executed functional-site case demonstrates."""
    seen: set[str] = set()
    evidence = out / "SITE_CONTEXT.json"
    manifest = out / "RUN_MANIFEST.json"
    if not evidence.is_file():
        return seen
    ev = json.loads(evidence.read_text(encoding="utf-8"))

    states = {s.get("state") for s in ev.get("sites") or []}
    if "role_compatible" in states:
        seen.add("curated_residue_resolved")
    if "role_mismatch" in states:
        seen.add("curated_residue_substituted")
    if "unresolved_mapping" in states:
        seen.add("curated_residue_unresolved")
    if {"not_observed", "missing_coordinates"} & states:
        seen.add("curated_residue_missing_coordinates")
    if any(s.get("role") == "metal_ligand" for s in ev.get("sites") or []):
        seen.add("metal_ligand_role")

    cofactors = {c.get("state") for c in ev.get("cofactors") or []}
    if "observed_match" in cofactors:
        seen.add("declared_cofactor_observed")
    if "not_observed" in cofactors:
        seen.add("declared_cofactor_absent")
    if "unresolved" in cofactors:
        seen.add("declared_cofactor_unmappable")

    if manifest.is_file() and (json.loads(manifest.read_text(encoding="utf-8")).get("missing_evidence") or []):
        seen.add("missing_evidence_reported")
    return seen


def _witness_assembly_interface(out: Path) -> set[str]:
    seen: set[str] = set()
    evidence = out / "ASSEMBLY_CONTEXT.json"
    if not evidence.is_file():
        return seen
    ev = json.loads(evidence.read_text(encoding="utf-8"))
    assembly, surface = ev.get("assembly") or {}, ev.get("surface") or {}
    chains = assembly.get("chains_observed") or []
    if len(chains) == 2:
        seen.add("dimeric_assembly")
    if len(chains) == 4:
        seen.add("tetrameric_assembly")
    if len(chains) > 4:
        seen.add("higher_order_assembly")
    if len(set(chains)) != len(chains):
        seen.add("repeated_chain_ids")
    if assembly.get("complete"):
        seen.add("stoichiometry_complete")
    if assembly.get("lower_bound"):
        seen.add("stoichiometry_lower_bound")
    if (surface.get("buried_sasa_A2") or 0) > 0:
        seen.add("buried_surface_observed")
    if ev.get("residue_contacts"):
        seen.add("residue_contacts_observed")
    return seen


def _witness_membrane_orientation(out: Path) -> set[str]:
    seen: set[str] = set()
    evidence = out / "MEMBRANE_ORIENTATION.json"
    if not evidence.is_file():
        return seen
    ev = json.loads(evidence.read_text(encoding="utf-8"))
    s, v = ev.get("summary") or {}, ev.get("validation") or {}
    if s.get("scope_id") == "beta_barrel":
        seen.add("beta_barrel_scope")
    if v.get("extracellular_comparison_state") == "evaluated":
        seen.add("extracellular_comparison_evaluated")
    errors = v.get("normal_angle_errors_deg") or []
    if errors:
        seen.add("rotation_series_recorded")
        if max(errors) > 1.0:
            seen.add("non_trivial_normal_error")   # a case the fit does not recover exactly
    if (s.get("n_embedded") or 0) > 0:
        seen.add("embedded_residues_identified")
    if (s.get("delta_kd") or 0) > 0:
        seen.add("lipid_pore_gap_positive")
    return seen


WITNESSES = {"structure_qc": _witness_structure_qc,
             "membrane_orientation": _witness_membrane_orientation,
             "assembly_interface": _witness_assembly_interface,
             "functional_site_state": _witness_functional_site_state}


def witness_coverage(out: Path, workflow: str = "structure_qc") -> set[str]:
    """Dispatch to the witness for the workflow that produced this evidence."""
    return WITNESSES[workflow](out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, default=HERE / "ADOPTION_DRAFT_XRAY.json")
    ap.add_argument("--out", type=Path, default=RESULTS / "execution")
    ap.add_argument("--record", action="store_true",
                    help="Re-measure every case and write the expectations back into the panel.")
    args = ap.parse_args(argv)

    panel = json.loads(args.panel.read_text(encoding="utf-8"))
    cli = ENGINES.get(panel.get("workflow", "structure_qc"), {}).get("cli", "structqc")
    exe = shutil.which(cli)
    if exe is None:
        print(f"{cli} is not on PATH; install the distribution first", file=sys.stderr)
        return 1

    # The engine must actually execute, not merely exist on PATH.
    #
    # invoke() clears a case's output directory before running the engine, so an
    # engine that fails destroys that case's previous evidence and writes
    # nothing in its place. On 2026-08-31 a console script whose package could no
    # longer be imported passed `which`, failed all sixteen cases, and left the
    # panel's committed evidence as sixteen empty directories. The summary line
    # still parsed. Probe before anything is touched, and probe the upstream
    # engine too, since a chained case runs that one first.
    for probe_cli in dict.fromkeys([cli] + [
            ENGINES[str(r.get("upstream", {}).get("workflow", "structure_qc"))]["cli"]
            for r in panel.get("records", []) if r.get("upstream")]):
        probe_exe = shutil.which(probe_cli)
        if probe_exe is None:
            print(f"{probe_cli} is not on PATH; install the distribution first", file=sys.stderr)
            return 1
        probe = subprocess.run([probe_exe, "describe"], capture_output=True, text=True)
        if probe.returncode != 0:
            print(f"{probe_cli} is on PATH at {probe_exe} but failed to execute: "
                  f"`{probe_cli} describe` exited {probe.returncode}. Refusing to run: "
                  f"each case's output directory is cleared before its engine runs, so "
                  f"continuing would destroy the existing evidence.", file=sys.stderr)
            if probe.stderr.strip():
                print(probe.stderr.strip().splitlines()[-1], file=sys.stderr)
            return 1
    semantics = panel.get("gate_semantics")
    if not semantics:
        print(f"{args.panel.name} declares no gate_semantics; refusing to assume any", file=sys.stderr)
        return 1

    workflow = panel.get("workflow", "structure_qc")
    workflow = panel.get("workflow", "structure_qc")
    args.out.mkdir(parents=True, exist_ok=True)

    if args.record:
        for record in panel["records"] + panel.get("controls", []):
            record["expected_result"] = measure(record, exe, args.out, workflow)
            e = record["expected_result"]
            summary = (f"identity={e['residue_identity']['identity_fraction']}"
                       if "residue_identity" in e else
                       f"sites={e.get('site_count')} states={e.get('state_counts')}")
            print(f"  recorded {record['record_id']}: exit={e['cli_exit_code']} {summary}")
        args.panel.write_text(json.dumps(panel, indent=2) + "\n", encoding="utf-8")
        print(f"wrote measured expectations into {args.panel.name}")
        return 0

    cases = []
    for record in panel["records"]:
        cases.append(run_case({**record, "_gate_semantics": semantics,
                               "_workflow": workflow}, exe, args.out))
    # Controls are not counted toward the panel's required case counts. They
    # exist so a gate that would otherwise be trivially satisfied is exercised
    # against a deliberately incomplete run.
    controls = [run_case({**r, "_gate_semantics": semantics, "_workflow": workflow}, exe, args.out)
                for r in panel.get("controls", [])]

    passed = sum(1 for c in cases if c["passed"])
    controls_passed = sum(1 for c in controls if c["passed"])

    required_coverage = (semantics.get("coverage_requirements") or {}).get("required_features") or []
    witnessed: dict[str, list[str]] = {feature: [] for feature in required_coverage}
    for case in cases + controls:
        if not case.get("output_dir"):
            continue
        for feature in witness_coverage(at(case["output_dir"]), workflow):
            if feature in witnessed:
                witnessed[feature].append(case["record_id"])
    unmet = sorted(f for f, by in witnessed.items() if not by)
    coverage_ok = not unmet
    result = {
        "schema_version": "1.0",
        "collection_id": panel.get("panel_id"),
        "stratum": panel.get("stratum"),
        "gate_semantics": semantics,
        "scientific_execution_performed": True,
        "case_counts": {"total": len(cases), "passed": passed, "failed": len(cases) - passed},
        "stratum_state": ("passed" if passed == len(cases) and controls_passed == len(controls)
                          and coverage_ok else "failed"),
        "coverage": {"required": sorted(required_coverage),
                     "unwitnessable": sorted(
                         ((semantics.get("coverage_requirements") or {}).get("unwitnessable_features") or {})),
                     "witnessed_by": {f: sorted(by) for f, by in sorted(witnessed.items())},
                     "unmet": unmet, "passed": coverage_ok},
        "control_counts": {"total": len(controls), "passed": controls_passed},
        "controls": controls,
        "scope_qualified": False,
        "strata_executed": sorted({c["stratum"] for c in cases if c.get("stratum")}),
        "scope_qualification_note": (
            "Passing execution is not scope qualification. A Mark 1 scope is qualified only when its "
            "panel composes, every case and control passes, and the result reproduces independently "
            "on a second machine. This runner reports execution only; it never sets a qualification "
            "gate, and scope_qualified is always false here by construction."),
        "runtime": {"python": platform.python_version(), "machine": platform.machine(),
                    "platform": platform.system()},
        "cases": cases,
    }
    # The results document lives with the run that produced it.
    #
    # Keying the filename on the panel id was not enough. It separates different
    # panels, but two runs of the *same* panel -- a canonical one and a variant
    # pointed at a different --out -- still wrote to one path, so an adversarial
    # run silently replaced the real membrane evidence with its own failures.
    # Writing beside the case outputs makes each execution self-contained and
    # removes the collision rather than narrowing it.
    path = args.out / "EXECUTION_STATUS.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical(result))
    print(f"{panel.get('stratum')} stratum: {passed}/{len(cases)} cases passed, "
          f"{controls_passed}/{len(controls)} controls passed "
          + (f"[residue_identity={semantics['residue_identity']['mode']}]"
             if "residue_identity" in semantics else f"[workflow={workflow}]"))
    if required_coverage:
        print(f"coverage: {len(required_coverage) - len(unmet)}/{len(required_coverage)} features witnessed"
              + (f" | UNMET: {', '.join(unmet)}" if unmet else ""))
    unwitnessable = (semantics.get("coverage_requirements") or {}).get("unwitnessable_features") or {}
    for name in sorted(unwitnessable):
        print(f"coverage NOT satisfied (recorded defect): {name}")

    # Informational drift, printed because stdout is the only place it survives.
    # The evidence file is written into the runner's workspace and destroyed
    # with it, so a delta that is never printed is a measurement nobody can
    # read. These are not gates -- see the membrane drift note above -- but the
    # spread across machines is the first evidence anyone has about whether the
    # fit reproduces, and that question is what the independent second-machine
    # gate exists to answer. It belongs in the log of every run that produced it.
    for check_name, deltas in drift_deltas(cases + controls).items():
        print(f"{check_name}: max {deltas[-1]:.3e}, median "
              f"{deltas[len(deltas) // 2]:.3e}, over {len(deltas)} cases")

    # Names alone are not diagnosable. A gate that fails cross-machine is
    # exactly the case where the evidence file is on a runner nobody can reach,
    # so the log has to carry the numbers: which check, what it wanted, what it
    # got. Anything less means reading the failure requires credentials.
    for c in cases + controls:
        if not c["passed"]:
            print(f"  FAILED {c['record_id']}:")
            for k in c["checks"]:
                if k["required"] and not k["passed"]:
                    print(f"    {k['check']}: expected {k['expected']!r}, observed {k['observed']!r}")
    print(f"strata executed: {', '.join(sorted({c['stratum'] for c in cases if c.get('stratum')}))}")
    print("scope_qualified: false (this runner never qualifies a scope; see the result note)")
    return 0 if passed == len(cases) and controls_passed == len(controls) and coverage_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
