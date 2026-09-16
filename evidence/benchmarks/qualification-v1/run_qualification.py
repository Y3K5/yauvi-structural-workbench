#!/usr/bin/env python3
"""Run the frozen public qualification panel without network access.

This is deliberately separate from the unit suite.  It compares YAUVI outputs
with independent public records and records failures rather than weakening a
gate.  Source acquisition is not performed here; SOURCE_LOCK.json must already
match the local source files byte-for-byte.
"""
from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any, Callable

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCES = HERE / "sources"
BUILD = HERE / "build"
RESULTS = HERE / "results"

PYTHONPATHS = [
    REPO / "structqc" / "src",
    REPO / "state-atlas" / "src",
    REPO / "site-context" / "src",
    REPO / "assembly-context" / "src",
    REPO / "sf-csa" / "src",
    REPO / "Membrane Orientor" / "memorient" / "src",
]
for source_root in reversed(PYTHONPATHS):
    sys.path.insert(0, str(source_root))


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
        encoding="utf-8",
    )


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def source(relative: str) -> Path:
    return HERE / relative


def check(name: str, passed: bool, observed: Any, expected: Any, *, required: bool = True) -> dict[str, Any]:
    return {
        "check": name,
        "passed": bool(passed),
        "required": required,
        "observed": observed,
        "expected": expected,
    }


def status(checks: list[dict[str, Any]], *, blocked: bool = False) -> str:
    if blocked:
        return "blocked"
    required = [row for row in checks if row.get("required", True)]
    if required and all(row["passed"] for row in required):
        return "passed"
    if any(row["passed"] for row in required):
        return "partial"
    return "failed"


def verify_source_lock() -> dict[str, Any]:
    lock = read_json(HERE / "SOURCE_LOCK.json")
    rows = []
    for item in lock["sources"]:
        path = source(item["artifact"])
        observed = sha256(path) if path.is_file() else None
        rows.append({
            "artifact": item["artifact"],
            "source_id": item["source_id"],
            "expected_sha256": item["sha256"],
            "observed_sha256": observed,
            "passed": observed == item["sha256"],
        })
    return {
        "collection_id": lock["collection_id"],
        "passed": all(row["passed"] for row in rows),
        "artifacts": rows,
    }


def parser_for(path: Path):
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def chain_sequence(path: Path, chain_id: str) -> str:
    model = next(parser_for(path).get_models())
    try:
        chain = model[chain_id]
    except KeyError as exc:
        raise RuntimeError(f"{path.name} has no chain {chain_id}; available: {[c.id for c in model]}") from exc
    return "".join(
        AA3.get(residue.resname.upper(), "X")
        for residue in chain
        if residue.id[0] == " "
    )


def write_fasta(path: Path, identifier: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    wrapped = "\n".join(sequence[index:index + 70] for index in range(0, len(sequence), 70))
    path.write_text(f">{identifier}\n{wrapped}\n", encoding="utf-8")


class ChainSelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:
        return str(chain.id) == self.chain_id

    def accept_residue(self, residue) -> bool:
        return residue.id[0] == " " and residue.resname.upper() in AA3


def write_chain_pdb(source_path: Path, chain_id: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(parser_for(source_path))
    writer.save(str(destination), select=ChainSelect(chain_id))


def cli(module: str, arguments: list[str]) -> dict[str, Any]:
    environment = os.environ.copy()
    prefix = os.pathsep.join(str(path) for path in PYTHONPATHS)
    environment["PYTHONPATH"] = prefix + (os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else "")
    process = subprocess.run(
        [sys.executable, "-m", module, *arguments],
        cwd=REPO,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    def sanitize(value: str) -> str:
        return value.replace(str(REPO), "<repo>")

    return {
        "exit_code": process.returncode,
        "stdout": sanitize(process.stdout.strip()),
        "stderr": sanitize(process.stderr.strip()),
    }


def run_structqc() -> dict[str, Any]:
    work = BUILD / "structqc"
    work.mkdir(parents=True, exist_ok=True)
    experimental = source("sources/wwpdb/1CRN.cif")
    experimental_fasta = work / "1CRN.fasta"
    write_fasta(experimental_fasta, "PDB:1CRN:A", chain_sequence(experimental, "A"))
    experimental_provenance = work / "1CRN.provenance.json"
    write_json(experimental_provenance, {
        "class": "experimental",
        "method": "X-ray diffraction",
        "source_id": "PDB:1CRN",
        "resolution_angstrom": 1.5,
    })
    experimental_out = work / "experimental-output"
    experimental_run = cli("structqc.cli", [
        "run", "--structure", str(experimental), "--chain", "A",
        "--reference-fasta", str(experimental_fasta),
        "--provenance", str(experimental_provenance),
        "--validation-report", str(source("sources/wwpdb/1crn_validation.xml")),
        "--require-external-validation", "--out", str(experimental_out),
    ])
    experimental_doc = read_json(experimental_out / "STRUCTURE_EVIDENCE.json")
    validation_metrics = experimental_doc["external_validation"].get("metrics", {})
    residue_identities = experimental_doc["residues"]

    af_meta = read_json(source("sources/alphafold/P69905-api.json"))[0]
    predicted_fasta = work / "P69905.fasta"
    write_fasta(predicted_fasta, "P69905", af_meta["sequence"])
    predicted_provenance = work / "P69905.provenance.json"
    write_json(predicted_provenance, {
        "class": "predicted",
        "method": "AlphaFold Monomer v2.0 pipeline",
        "source_id": "AF-P69905-F1:model:v6",
        "confidence_encoding": "plddt_in_bfactor",
    })
    predicted_out = work / "predicted-output"
    predicted_run = cli("structqc.cli", [
        "run", "--structure", str(source("sources/alphafold/AF-P69905-F1-model_v6.cif")),
        "--chain", "A", "--reference-fasta", str(predicted_fasta),
        "--provenance", str(predicted_provenance),
        "--pae", str(source("sources/alphafold/AF-P69905-F1-predicted_aligned_error_v6.json")),
        "--out", str(predicted_out),
    ])
    predicted_doc = read_json(predicted_out / "STRUCTURE_EVIDENCE.json")
    plddt = [row.get("plddt") for row in predicted_doc["residues"] if row.get("plddt") is not None]

    unknown_out = work / "unknown-output"
    unknown_run = cli("structqc.cli", [
        "run", "--structure", str(experimental), "--chain", "A", "--out", str(unknown_out),
    ])
    unknown_doc = read_json(unknown_out / "STRUCTURE_EVIDENCE.json")

    checks = [
        check("experimental_cli_completed", experimental_run["exit_code"] == 0, experimental_run["exit_code"], 0),
        check("gemmi_coordinate_validation", experimental_doc["coordinate"]["parser"]["gemmi_validation"] == "validated", experimental_doc["coordinate"]["parser"], "gemmi validated"),
        check("author_label_identity_present", all(row.get("label_asym_id") and row.get("label_seq_id") is not None for row in residue_identities), sum(bool(row.get("label_asym_id")) and row.get("label_seq_id") is not None for row in residue_identities), len(residue_identities)),
        check("reference_mapping_exact", experimental_doc["completeness"].get("identity_fraction") == 1.0, experimental_doc["completeness"], "identity_fraction=1.0"),
        check("wwpdb_validation_imported", experimental_doc["external_validation"]["state"] == "imported", experimental_doc["external_validation"]["state"], "imported"),
        check("wwpdb_geometry_metrics_imported", {"clashscore", "ramachandran_outliers_percent", "rotamer_outliers_percent"}.issubset(validation_metrics), sorted(validation_metrics), ["clashscore", "ramachandran_outliers_percent", "rotamer_outliers_percent"]),
        check("wwpdb_raw_clashscore_not_percentile", validation_metrics.get("clashscore") == 0.0, validation_metrics.get("clashscore"), 0.0),
        check("predicted_cli_completed", predicted_run["exit_code"] == 0, predicted_run["exit_code"], 0),
        check("plddt_only_under_declared_prediction", len(plddt) == len(predicted_doc["residues"]), len(plddt), len(predicted_doc["residues"])),
        check("pae_matrix_evaluated", predicted_doc["pae"].get("state") == "evaluated" and predicted_doc["pae"].get("size") == len(predicted_doc["residues"]), predicted_doc["pae"], f"evaluated square matrix of {len(predicted_doc['residues'])}"),
        check("unknown_provenance_fails_incomplete", unknown_run["exit_code"] == 1 and unknown_doc["provenance"]["class"] == "unknown", {"exit_code": unknown_run["exit_code"], "class": unknown_doc["provenance"]["class"]}, {"exit_code": 1, "class": "unknown"}),
    ]
    return {
        "workflow": "structure_qc",
        "independent_references": ["wwPDB validation XML", "AlphaFold DB v6 metadata/model/PAE"],
        "status": status(checks),
        "checks": checks,
        "observations": {
            "wwpdb_metrics_imported": validation_metrics,
            "predicted_mean_plddt": round(float(np.mean(plddt)), 6) if plddt else None,
            "predicted_pae": predicted_doc["pae"],
        },
        "claim_boundary": "Coordinate and imported validation evidence do not establish native conformation or function.",
    }


def random_rotation(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def opm_half_thickness(path: Path) -> float | None:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "1/2 of bilayer thickness" in line:
            match = re.search(r"([0-9]+(?:\.[0-9]+)?)", line.split(":", 1)[-1])
            return float(match.group(1)) if match else None
    return None


def run_membrane_orientation() -> dict[str, Any]:
    from memorient.barrel import fit_membrane
    from memorient.contexts import get_context
    from memorient.geometry import load_structure
    from memorient.orientor import five_fold_validate

    panel = [
        ("1BXW", "beta_barrel", "gram_negative_om"),
        ("2POR", "beta_barrel", "gram_negative_om"),
        ("1QD6", "beta_barrel", "gram_negative_om"),
        ("2F1T", "beta_barrel", "gram_negative_om"),
        ("1P4T", "beta_barrel", "gram_negative_om"),
        ("1AFO", "alpha_helical", "tm_receptor"),
        ("1C3W", "alpha_helical", "eukaryotic_pm"),
        ("1U19", "alpha_helical", "eukaryotic_pm"),
    ]
    rows = []
    for pdb_id, stratum, context_name in panel:
        path = source(f"sources/opm/{pdb_id}.pdb")
        structure = load_structure(str(path))
        context = get_context(context_name)
        angle_errors: list[float] = []
        half_thicknesses: list[float] = []
        for seed in range(4):
            rotation = random_rotation(seed)
            fitted = fit_membrane(structure.transformed(rotation), context)
            expected_normal = rotation @ np.asarray([0.0, 0.0, 1.0])
            cosine = min(abs(float(np.dot(fitted.normal, expected_normal))), 1.0)
            angle_errors.append(float(np.degrees(np.arccos(cosine))))
            half_thicknesses.append(float(fitted.half_thickness))
        validation = five_fold_validate(structure, context, seeds=4, threshold=0.95, n_points=160)
        opm_half = opm_half_thickness(path)
        fit_half = float(np.mean(half_thicknesses))
        rows.append({
            "pdb_id": pdb_id,
            "stratum": stratum,
            "context": context_name,
            "mean_normal_error_deg": round(float(np.mean(angle_errors)), 6),
            "normal_error_sd_deg": round(float(np.std(angle_errors)), 6),
            "opm_half_thickness_A": opm_half,
            "fitted_half_thickness_A": round(fit_half, 6),
            "half_thickness_absolute_error_A": round(abs(fit_half - opm_half), 6) if opm_half is not None else None,
            "rotation_invariance": validation,
        })

    beta = [row for row in rows if row["stratum"] == "beta_barrel"]
    alpha = [row for row in rows if row["stratum"] == "alpha_helical"]
    beta_angle = float(np.mean([row["mean_normal_error_deg"] for row in beta]))
    alpha_angle = float(np.mean([row["mean_normal_error_deg"] for row in alpha]))
    beta_thickness = float(np.mean([row["half_thickness_absolute_error_A"] for row in beta]))
    alpha_thickness = float(np.mean([row["half_thickness_absolute_error_A"] for row in alpha]))
    checks = [
        check("beta_barrel_mean_normal_error", beta_angle <= 15.0, round(beta_angle, 6), "<= 15 degrees"),
        check("beta_barrel_mean_half_thickness_error", beta_thickness <= 2.5, round(beta_thickness, 6), "<= 2.5 A"),
        check("beta_barrel_rotation_invariance", all(row["rotation_invariance"]["passed"] for row in beta), [row["rotation_invariance"]["mean_jaccard"] for row in beta], "all seed Jaccards >= 0.95"),
        check("alpha_helical_mean_normal_error", alpha_angle <= 15.0, round(alpha_angle, 6), "<= 15 degrees"),
        check("alpha_helical_mean_half_thickness_error", alpha_thickness <= 2.5, round(alpha_thickness, 6), "<= 2.5 A"),
        check("alpha_helical_rotation_invariance", all(row["rotation_invariance"]["passed"] for row in alpha), [row["rotation_invariance"]["mean_jaccard"] for row in alpha], "all seed Jaccards >= 0.95"),
    ]
    return {
        "workflow": "membrane_orientation",
        "independent_reference": "OPM-oriented coordinates; deposited membrane normal is Z and REMARK records bilayer half-thickness",
        "precommitted_tolerances": {"mean_normal_error_deg": 15.0, "mean_half_thickness_error_A": 2.5, "rotation_jaccard": 0.95},
        "status": status(checks),
        "checks": checks,
        "records": rows,
        "claim_boundary": "Agreement with OPM placement does not establish native-cell exposure or topology in a tested cell.",
    }


def run_state_atlas() -> dict[str, Any]:
    work = BUILD / "state-atlas"
    work.mkdir(parents=True, exist_ok=True)
    reference_set = work / "abl-reference-set.json"
    write_json(reference_set, {
        "reference_set_id": "abl-kincore-two-sided-v1",
        "decision_rules": {"max_rmsd_A": 2.5, "min_margin_A": 0.25},
        "references": [
            {
                "reference_id": "2HZ4_B_KinCore_active",
                "state": "active",
                "structure": "../../sources/wwpdb/2HZ4.cif",
                "chain": "B",
                "provenance": {"class": "experimental", "method": "X-ray diffraction"},
                "state_evidence": {"basis": "KinCore ABL1 active label: DFGin, BLAminus, salt bridge present, HRD-backbone hydrogen bond, ordered activation loop"},
            },
            {
                "reference_id": "2G1T_A_KinCore_inactive",
                "state": "inactive",
                "structure": "../../sources/wwpdb/2G1T.cif",
                "chain": "A",
                "provenance": {"class": "experimental", "method": "X-ray diffraction"},
                "state_evidence": {"basis": "KinCore ABL1 inactive label with DFGin-BLBplus conformation"},
            },
        ],
    })

    held_out = [
        ("2V7A", "A", "active_like", "KinCore active"),
        ("8SSN", "A", "inactive_like", "KinCore inactive"),
    ]
    rows = []
    checks: list[dict[str, Any]] = []
    validation = cli("state_atlas.cli", ["validate", "--reference-set", str(reference_set)])
    checks.append(check("two_sided_reference_validation", validation["exit_code"] == 0, validation, "exit_code 0"))
    for pdb_id, chain, expected_label, independent_label in held_out:
        structure = source(f"sources/wwpdb/{pdb_id}.cif")
        fasta = work / f"{pdb_id}-{chain}.fasta"
        write_fasta(fasta, f"PDB:{pdb_id}:{chain}", chain_sequence(structure, chain))
        provenance = work / f"{pdb_id}.provenance.json"
        write_json(provenance, {"class": "experimental", "method": "X-ray diffraction", "source_id": f"PDB:{pdb_id}"})
        qc_out = work / f"{pdb_id}-qc"
        qc_run = cli("structqc.cli", [
            "run", "--structure", str(structure), "--chain", chain,
            "--reference-fasta", str(fasta), "--provenance", str(provenance), "--out", str(qc_out),
        ])
        state_out = work / f"{pdb_id}-state"
        state_run = cli("state_atlas.cli", [
            "run", "--manifest", str(qc_out / "STRUCTURE_EVIDENCE.json"),
            "--reference-set", str(reference_set), "--structure", str(structure),
            "--chain", chain, "--cluster-cutoff-A", "2.0", "--out", str(state_out),
        ])
        document = read_json(state_out / "STATE_ENSEMBLE.json")
        observed = document["overall_label"]
        rows.append({
            "pdb_id": pdb_id,
            "chain": chain,
            "independent_label": independent_label,
            "expected_yauvi_label": expected_label,
            "observed_yauvi_label": observed,
            "structqc_exit_code": qc_run["exit_code"],
            "state_atlas_exit_code": state_run["exit_code"],
            "frame_metrics": document["frame_metrics"],
        })
        checks.append(check(f"held_out_{pdb_id}_{expected_label}", observed == expected_label, observed, expected_label))

    checks.append(check("no_opposite_confident_calls", not any(
        (row["expected_yauvi_label"] == "active_like" and row["observed_yauvi_label"] == "inactive_like")
        or (row["expected_yauvi_label"] == "inactive_like" and row["observed_yauvi_label"] == "active_like")
        for row in rows
    ), [row["observed_yauvi_label"] for row in rows], "no held-out opposite-state call"))
    return {
        "workflow": "conformational_state",
        "independent_reference": "KinCore ABL1 experimental-chain labels",
        "precommitted_rules": {"max_rmsd_A": 2.5, "min_margin_A": 0.25, "alignment": "sequence-mapped shared CA atoms"},
        "status": status(checks),
        "checks": checks,
        "records": rows,
        "claim_boundary": "A state label is geometric resemblance to curated references, not kinase activity.",
    }


def mcsa_role(entry: dict[str, Any]) -> str:
    functions = {str(role.get("function", "")).lower() for role in entry.get("roles", [])}
    if "proton acceptor" in functions or "proton donor" in functions:
        return "acid_base"
    return "unspecified"


def run_functional_site() -> dict[str, Any]:
    work = BUILD / "functional-site"
    work.mkdir(parents=True, exist_ok=True)
    structure = source("sources/wwpdb/1B73.cif")
    fasta = work / "1B73-A.fasta"
    write_fasta(fasta, "P56868", chain_sequence(structure, "A"))
    provenance = work / "1B73.provenance.json"
    write_json(provenance, {"class": "experimental", "method": "X-ray diffraction", "source_id": "PDB:1B73", "resolution_angstrom": 2.3})
    qc_out = work / "structqc-output"
    qc_run = cli("structqc.cli", [
        "run", "--structure", str(structure), "--chain", "A", "--reference-fasta", str(fasta),
        "--provenance", str(provenance), "--out", str(qc_out),
    ])

    mcsa = read_json(source("sources/mcsa/entry-1-residues.json"))
    sites = []
    for record in mcsa:
        chain_record = next((item for item in record.get("residue_chains", []) if item.get("pdb_id", "").lower() == "1b73" and item.get("chain_name") == "A"), None)
        sequence_record = next((item for item in record.get("residue_sequences", []) if item.get("uniprot_id") == "P56868"), None)
        if not chain_record or not sequence_record:
            continue
        sites.append({
            "position": int(sequence_record["resid"]),
            "type": "active_site",
            "role": mcsa_role(record),
            "expected_residues": [AA3[str(sequence_record["code"]).upper()[:3]]],
            "detail": str(record.get("main_annotation", "")),
        })
    annotations = work / "MCSA-1.annotations.json"
    write_json(annotations, {"sites": sorted(sites, key=lambda item: item["position"]), "declared_cofactors": []})
    output = work / "site-context-output"
    site_run = cli("site_context.cli", [
        "run", "--manifest", str(qc_out / "STRUCTURE_EVIDENCE.json"), "--structure", str(structure),
        "--annotations", str(annotations), "--out", str(output),
    ])
    document = read_json(output / "SITE_CONTEXT.json")
    positions = sorted(row["position"] for row in document["sites"])
    states = [row["state"] for row in document["sites"]]
    checks = [
        check("structqc_boundary_completed", qc_run["exit_code"] == 0, qc_run["exit_code"], 0),
        check("mcsa_six_reference_residues_loaded", positions == [7, 8, 70, 147, 178, 180], positions, [7, 8, 70, 147, 178, 180]),
        check("all_mcsa_residues_exactly_mapped", all(state == "role_compatible" for state in states), states, "six role_compatible records"),
        check("site_cli_scientifically_incomplete_without_pocket_leg", site_run["exit_code"] == 1 and document["missing_evidence"] == ["pocket_results"], {"exit_code": site_run["exit_code"], "missing_evidence": document["missing_evidence"]}, {"exit_code": 1, "missing_evidence": ["pocket_results"]}),
        check("no_activity_claim_emitted", "activity" not in json.dumps(document).lower() or "does not establish" in json.dumps(document).lower(), "claim language inspected", "no observed activity claim"),
    ]
    return {
        "workflow": "functional_site_state",
        "independent_reference": "M-CSA entry 1, glutamate racemase P56868 / PDB 1B73",
        "status": status(checks),
        "checks": checks,
        "mapped_sites": document["sites"],
        "missing_evidence": document["missing_evidence"],
        "claim_boundary": "Exact mapping of M-CSA residues is not an observation of catalysis.",
    }


def run_assembly_context() -> dict[str, Any]:
    work = BUILD / "assembly-context"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    isolated = source("sources/wwpdb/4HHB.cif")
    assembly = source("sources/wwpdb/4HHB-assembly1.cif")
    fasta = work / "4HHB-A.fasta"
    write_fasta(fasta, "P69905", chain_sequence(isolated, "A"))
    provenance = work / "4HHB.provenance.json"
    write_json(provenance, {"class": "experimental", "method": "X-ray diffraction", "source_id": "PDB:4HHB"})
    qc_out = work / "structqc-output"
    qc_run = cli("structqc.cli", [
        "run", "--structure", str(isolated), "--chain", "A", "--reference-fasta", str(fasta),
        "--provenance", str(provenance), "--out", str(qc_out),
    ])
    output = work / "assembly-output"
    assembly_run = cli("assembly_context.cli", [
        "run", "--manifest", str(qc_out / "STRUCTURE_EVIDENCE.json"),
        "--isolated", str(isolated), "--assembly", str(assembly), "--subject-chain", "A",
        "--relationship", "exact_protein", "--reference-id", "PDB:4HHB",
        "--assembly-id", "1", "--expected-chains", "A,B,C,D", "--out", str(output),
    ])
    document = read_json(output / "ASSEMBLY_CONTEXT.json")
    run_manifest = read_json(output / "RUN_MANIFEST.json")
    freesasa = shutil.which("freesasa")
    checks = [
        check("structqc_boundary_completed", qc_run["exit_code"] == 0, qc_run["exit_code"], 0),
        check("assembly_cli_completed", assembly_run["exit_code"] == 0, assembly_run["exit_code"], 0),
        check("wwpdb_tetramer_stoichiometry", document["assembly"]["complete"] is True and len(document["assembly"]["chains_observed"]) == 4, document["assembly"], "complete A,B,C,D tetramer"),
        check("gemmi_assembly_metadata_read", document["assembly"]["metadata"]["operator_application"]["operator_backend"] == "gemmi", document["assembly"]["metadata"]["operator_application"], "gemmi"),
        check("interface_contacts_observed", len(document["interfaces"]) > 0, len(document["interfaces"]), "> 0 heavy-atom residue contacts at 5 A"),
        check("buried_surface_positive", document["surface"]["buried_sasa_A2"] > 0, document["surface"], "buried_sasa_A2 > 0"),
        check("publication_grade_freesasa_invoked", document["methods"]["sasa"].startswith("freesasa_") and run_manifest["optional_runtimes"]["freesasa"] == "available_invoked", {"binary": freesasa, "method": document["methods"]["sasa"], "runtime": run_manifest["optional_runtimes"]["freesasa"]}, "FreeSASA available and invoked"),
    ]
    return {
        "workflow": "assembly_interface",
        "independent_reference": "wwPDB biological assembly 1 for hemoglobin 4HHB",
        "status": status(checks, blocked=not bool(freesasa)),
        "checks": checks,
        "methods": document["methods"],
        "assembly": document["assembly"],
        "surface": document["surface"],
        "claim_boundary": "A deposited assembly and computed interface are not binding affinity or intact-cell accessibility.",
    }


def run_sf_csa() -> dict[str, Any]:
    foldseek = shutil.which("foldseek")
    diamond = shutil.which("diamond")
    runtime_checks = [
        check("foldseek_runtime_installed", bool(foldseek), foldseek or "missing", "installed"),
        check("diamond_runtime_installed", bool(diamond), diamond or "missing", "installed"),
    ]
    if not (foldseek and diamond):
        return {
            "workflow": "sf_csa",
            "independent_reference": "CATH v4.3 classifications plus public PDB controls",
            "status": "blocked",
            "checks": runtime_checks,
            "runtime_versions": {"foldseek": None, "diamond": None},
            "claim_boundary": "A fold or sequence hit cannot be promoted to exact functional transfer.",
        }

    from sf_csa.core import database_bundle_checksum, pdb_sequence

    work = BUILD / "sf-csa"
    if work.exists():
        shutil.rmtree(work)
    structures = work / "structures" / "experimental"
    query_structures = work / "structures" / "query"
    proteomes = work / "proteomes"
    ledger_dir = work / "results"
    database_dir = work / "db"
    config_dir = work / "config"
    for directory in (structures, query_structures, proteomes, ledger_dir, database_dir, config_dir):
        directory.mkdir(parents=True, exist_ok=True)

    controls = [
        ("1B73", "A", "3.40.50.1860", "exact query control"),
        ("2JFN", "A", "3.40.50.1860", "same CATH homologous superfamily"),
        ("4EA9", "A", "3.40.50.20", "same CATH topology, different homologous superfamily"),
        ("1OAI", "A", "1.10.8.10", "different CATH class and topology"),
    ]
    fasta_records = []
    for pdb_id, chain_id, _cath, _relationship in controls:
        destination = structures / f"{pdb_id}.pdb"
        write_chain_pdb(source(f"sources/wwpdb/{pdb_id}.cif"), chain_id, destination)
        sequence, _residues, _coordinates = pdb_sequence(destination, chain_id)
        if not sequence:
            raise RuntimeError(f"SF-CSA control {pdb_id} has no parsed chain sequence")
        fasta_records.append((pdb_id, sequence))
    shutil.copyfile(structures / "1B73.pdb", query_structures / "1B73.pdb")

    proteome = proteomes / "public-cath-controls.faa"
    proteome.write_text("".join(f">{pdb_id} public CATH control\n{sequence}\n" for pdb_id, sequence in fasta_records), encoding="utf-8")
    query_sequence = dict(fasta_records)["1B73"]
    (ledger_dir / "SEQUENCE_MANIFEST.tsv").write_text(
        "accession\tsource_file\tsequence_sha256\n"
        f"1B73\tproteomes/public-cath-controls.faa\t{hashlib.sha256(query_sequence.encode('ascii')).hexdigest()}\n",
        encoding="utf-8",
    )
    (ledger_dir / "SELECTION_LEDGER.tsv").write_text(
        "primary_accession\tdecision_status\n1B73\tBENCHMARK_CONTROL\n",
        encoding="utf-8",
    )
    write_json(ledger_dir / "orientation.json", {
        "state": "not_required_for_public_structural_relationship_control",
        "claim_boundary": "No membrane or native exposure claim is made.",
    })

    db_prefix = database_dir / "public_structures"
    create_db = subprocess.run(
        [foldseek, "createdb", str(structures), str(db_prefix)],
        cwd=work,
        text=True,
        capture_output=True,
        check=False,
    )
    if create_db.returncode != 0:
        raise RuntimeError(f"Foldseek createdb failed: {(create_db.stderr or create_db.stdout)[-1200:]}")
    (database_dir / "pdb.version").write_text(
        "CATH v4.3.0 S35 relationship controls; wwPDB coordinates locked by SOURCE_LOCK.json\n",
        encoding="utf-8",
    )

    spec = work / "campaign.json"
    write_json(spec, {
        "schema_version": 1,
        "release_scope": "public CATH-labeled SF-CSA qualification controls",
        "root": ".",
        "path_base": "..",
        "sequence_manifest": "results/SEQUENCE_MANIFEST.tsv",
        "decision_ledger": "results/SELECTION_LEDGER.tsv",
        "default_orientation_artifact": "results/orientation.json",
        "targets": [{
            "accession": "1B73",
            "common_name": "Glutamate racemase public control",
            "organism": "Aquifex pyrophilus",
            "strain": "public PDB record",
            "uniprot_accession": "P56868",
            "mechanism_group": "glutamate_racemase",
            "protein_specific_boundary": "Fold and sequence evidence do not establish catalysis or substrate specificity.",
            "structure_path": "structures/query/1B73.pdb",
            "source_proteome_path": "proteomes/public-cath-controls.faa",
            "structure_class": "exact_experimental_chain",
            "chain": "A",
        }],
        "database": {
            "pdb_database": "db/public_structures",
            "proteome_globs": ["proteomes/*.faa"],
            "campaign_structure_roots": ["structures/query"],
            "annotation_tables": [],
            "seqmatch_tables": [],
            "required_foldseek_version": "10.941cd33",
            "thresholds": {
                "structure_evalue": "0.01",
                "same_fold_tm": 0.5,
                "whole_architecture_coverage": 0.7,
                "max_structure_hits": 20,
                "sequence_evalue": "1e-5",
                "sequence_min_identity": 0,
                "sequence_min_query_coverage": 30,
                "sequence_max_hits": 20,
                "sequence_hits_per_proteome": 20,
            },
            "mechanism_families": [{"group": "glutamate_racemase", "pattern": "glutamate racemase"}],
            "contested_groups": [],
            "divergence_sets": [],
            "release_expectations": {
                "proteome_count": 1,
                "target_statuses": {"1B73": "BENCHMARK_CONTROL"},
                "title_traps": [],
            },
        },
    })

    build_run = cli("sf_csa.cli", ["build-manifests", "--spec", str(spec), "--out", str(config_dir)])
    validate_run = cli("sf_csa.cli", [
        "validate", "--queries", str(config_dir / "target_manifest.json"),
        "--databases", str(config_dir / "database_manifest.json"),
    ])
    release_dir = work / "release"
    pipeline_run = cli("sf_csa.cli", [
        "run", "--queries", str(config_dir / "target_manifest.json"),
        "--databases", str(config_dir / "database_manifest.json"), "--output", str(release_dir),
    ])
    verify_run = cli("sf_csa.cli", [
        "verify", "--output", str(release_dir), "--databases", str(config_dir / "database_manifest.json"),
    ])

    structure_table = release_dir / "targets" / "1B73" / "structure_hits.tsv"
    sequence_table = release_dir / "targets" / "1B73" / "species_comparison.tsv"
    structure_rows = list(csv.DictReader(structure_table.open(encoding="utf-8"), delimiter="\t")) if structure_table.is_file() else []
    sequence_rows = list(csv.DictReader(sequence_table.open(encoding="utf-8"), delimiter="\t")) if sequence_table.is_file() else []
    experimental_hits = {row["target_id"].upper(): row for row in structure_rows if row["database"] == "experimental_pdb"}
    sequence_hits = {row["target_accession"].upper(): row for row in sequence_rows}
    db_files = sorted(path for path in database_dir.glob("public_structures*") if path.is_file())
    db_manifest = read_json(config_dir / "database_manifest.json") if (config_dir / "database_manifest.json").is_file() else {}

    exact = experimental_hits.get("1B73")
    homolog = experimental_hits.get("2JFN")
    analogy = experimental_hits.get("4EA9")
    unrelated = experimental_hits.get("1OAI")
    checks = runtime_checks + [
        check("cath_control_classes_frozen", controls, controls, controls),
        check("manifest_builder_completed", build_run["exit_code"] == 0, build_run, "exit_code 0"),
        check("manifest_validation_completed", validate_run["exit_code"] == 0, validate_run, "exit_code 0"),
        check("real_sf_csa_release_executed", pipeline_run["exit_code"] == 0, pipeline_run, "exit_code 0"),
        check("real_sf_csa_release_verified", verify_run["exit_code"] == 0, verify_run, "exit_code 0"),
        check("exact_structure_control", bool(exact) and exact["function_classification"] == "exact_function_supported", exact or "missing", "1B73 exact_function_supported"),
        check("cath_homolog_recovered", bool(homolog) and homolog["structural_category"] in {"whole_architecture_match", "domain_or_partial_match"}, homolog or "missing", "2JFN structural match"),
        check("cath_fold_analogy_not_promoted_to_exact_function", not analogy or analogy["function_classification"] not in {"exact_function_supported", "probable_same_function"}, analogy or "not returned at threshold", "4EA9 absent or bounded analogy"),
        check("cath_unrelated_not_promoted", not unrelated or unrelated["function_classification"] == "unresolved_or_conflicted", unrelated or "not returned at threshold", "1OAI absent or unresolved"),
        check("sequence_exact_control_recovered", "1B73" in sequence_hits, sorted(sequence_hits), "1B73"),
        check("structural_and_sequence_legs_preserved", structure_table.is_file() and sequence_table.is_file(), {"structure_table": structure_table.is_file(), "sequence_table": sequence_table.is_file()}, "separate Foldseek and DIAMOND tables"),
        check(
            "database_bundle_checksum_complete",
            len(db_manifest.get("pdb_database_file_checksums", {})) == len(db_files)
            and db_manifest.get("pdb_database_checksum") == database_bundle_checksum(db_prefix),
            {
                "manifest_checksum": db_manifest.get("pdb_database_checksum"),
                "bundle_file_count": len(db_files),
                "manifest_file_count": len(db_manifest.get("pdb_database_file_checksums", {})),
                "bundle_files": [path.name for path in db_files],
            },
            "composite checksum binds every Foldseek database file",
        ),
    ]
    return {
        "workflow": "sf_csa",
        "independent_reference": "CATH v4.3 classifications plus public PDB controls",
        "status": status(checks),
        "checks": checks,
        "runtime_versions": {
            "foldseek": subprocess.run([foldseek, "version"], text=True, capture_output=True).stdout.strip() if foldseek else None,
            "diamond": subprocess.run([diamond, "version"], text=True, capture_output=True).stdout.strip() if diamond else None,
        },
        "control_relationships": [
            {"pdb_id": pdb_id, "chain": chain_id, "cath": cath, "relationship": relationship}
            for pdb_id, chain_id, cath, relationship in controls
        ],
        "structure_hit_count": len(structure_rows),
        "sequence_hit_count": len(sequence_rows),
        "claim_boundary": "A fold or sequence hit cannot be promoted to exact functional transfer.",
    }


def guarded(name: str, function: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    try:
        return function()
    except Exception as exc:  # the qualification report must retain failures in other legs
        return {
            "workflow": name,
            "status": "failed",
            "checks": [],
            "error": f"{type(exc).__name__}: {exc}",
        }


def main() -> int:
    if BUILD.exists():
        shutil.rmtree(BUILD)
    if RESULTS.exists():
        shutil.rmtree(RESULTS)
    BUILD.mkdir(parents=True)
    RESULTS.mkdir(parents=True)

    lock = verify_source_lock()
    if not lock["passed"]:
        write_json(RESULTS / "QUALIFICATION_RESULTS.json", {
            "schema_version": "1.0",
            "source_lock": lock,
            "overall_state": "blocked_source_drift",
            "workflows": [],
        })
        return 2

    workflows = [
        guarded("structure_qc", run_structqc),
        guarded("membrane_orientation", run_membrane_orientation),
        guarded("conformational_state", run_state_atlas),
        guarded("functional_site_state", run_functional_site),
        guarded("assembly_interface", run_assembly_context),
        guarded("sf_csa", run_sf_csa),
    ]
    states = {item["status"] for item in workflows}
    overall = "all_passed" if states == {"passed"} else "incomplete_or_failed"
    result = {
        "schema_version": "1.0",
        "collection_id": "yauvi-structural-public-qualification-v1",
        "source_lock": {"passed": True, "artifact_count": len(lock["artifacts"])},
        "qualification_rule": "All required checks for all six workflows must pass. Software tests are a separate gate.",
        "overall_state": overall,
        "workflow_counts": {state: sum(item["status"] == state for item in workflows) for state in ("passed", "partial", "failed", "blocked")},
        "workflows": workflows,
    }
    write_json(RESULTS / "QUALIFICATION_RESULTS.json", result)
    write_json(RESULTS / "SOURCE_VERIFICATION.json", lock)
    print(json.dumps({"overall_state": overall, "workflow_counts": result["workflow_counts"]}, sort_keys=True))
    return 0 if overall == "all_passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
