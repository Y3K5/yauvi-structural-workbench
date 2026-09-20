#!/usr/bin/env python3
"""Run the hash-locked public CA II StructQC comparison into a new output folder."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parent
LOCK_PATH = ROOT / "public_research_sources.json"
PROTOCOL_PATH = ROOT / "PUBLIC_RESEARCH_PROTOCOL.md"
ACQUISITION_PATH = ROOT / "acquire_public_inputs.py"


def find_source_root(script_root: Path) -> Path:
    """Find the nearest full workbench source tree, or use this kit when copied alone."""
    script_root = script_root.resolve()
    for candidate in (script_root, *script_root.parents):
        if (candidate / "pyproject.toml").is_file() and (candidate / "software").is_dir():
            return candidate.resolve()
    return script_root


SOURCE_TREE_ROOT = find_source_root(ROOT)


def assert_outside_source(path: Path, source_root: Path, label: str) -> None:
    path = path.expanduser().resolve()
    source_root = source_root.resolve()
    if path == source_root or source_root in path.parents:
        raise ValueError(f"{label} must be outside the source tree: {path}")


def verify_input(path: Path, source: dict) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"required input missing: {source['path']}")
    actual_hash = sha256(path)
    actual_size = path.stat().st_size
    if actual_hash != source["sha256"] or actual_size != source["byte_count"]:
        raise RuntimeError(
            f"source-lock mismatch for {source['path']}: "
            f"sha256={actual_hash}, bytes={actual_size}"
        )
    return {"path": source["path"], "sha256": actual_hash, "byte_count": actual_size}


def validate_model_chain(entry_id: str, models: list, model_index: int, chain_id: str) -> None:
    if model_index < 0 or model_index >= len(models):
        raise RuntimeError(f"model {model_index} not found for {entry_id}")
    available = sorted(str(chain.id) for chain in models[model_index])
    if chain_id not in available:
        raise RuntimeError(
            f"selected chain {chain_id!r} not found for {entry_id}; "
            f"available chains: {available}"
        )


def require_new_output(path: Path) -> None:
    if path.exists():
        raise ValueError(f"output path already exists; choose a new directory: {path}")


def distribution_identity(package_name: str, module_file: Path) -> dict[str, object]:
    """Report distributions that provide a package, distinguishing source overrides."""
    candidates = importlib.metadata.packages_distributions().get(package_name, [])
    distributions = []
    matching_distributions = []
    for name in candidates:
        try:
            dist = importlib.metadata.distribution(name)
        except importlib.metadata.PackageNotFoundError:
            continue
        distributions.append({"name": dist.metadata.get("Name", name), "version": dist.version})
        for item in dist.files or ():
            if Path(dist.locate_file(item)).resolve() == module_file.resolve():
                matching_distributions.append({"name": dist.metadata.get("Name", name), "version": dist.version})
                break
    return {
        "distributions_providing_package": distributions,
        "module_distributions": matching_distributions,
        "module_origin_state": "installed_distribution" if matching_distributions else "source_tree_or_editable",
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path, help="directory containing the locked public files")
    parser.add_argument("--out", required=True, type=Path, help="new output directory outside this source tree")
    args = parser.parse_args()
    inputs = args.inputs.expanduser().resolve()
    out = args.out.expanduser().resolve()
    assert_outside_source(inputs, SOURCE_TREE_ROOT, "input directory")
    assert_outside_source(out, SOURCE_TREE_ROOT, "output directory")
    if not inputs.is_dir():
        parser.error(f"input directory does not exist: {inputs}")
    require_new_output(out)

    try:
        from structqc.core import analyze, read_fasta, write_outputs
        import Bio
        import numpy
        from Bio.PDB import MMCIFParser
    except ImportError as exc:
        parser.error(f"StructQC and its runtime dependencies must be installed or on PYTHONPATH: {exc}")

    lock = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    reference = lock["reference"]
    structures = lock["structures"]
    checked = []
    for source in [reference, *structures]:
        checked.append(verify_input(inputs / source["path"], source))

    fasta_path = inputs / reference["path"]
    ref_id, ref_sequence = read_fasta(fasta_path)
    if ref_id != reference["id"]:
        raise RuntimeError(f"reference identifier mismatch: expected {reference['id']}, observed {ref_id}")
    # Check every fixed chain before creating the result directory. A missing
    # chain must stop the case cleanly, before a partial comparison is written.
    for source in structures:
        coordinate_path = inputs / source["path"]
        parsed = MMCIFParser(QUIET=True).get_structure(source["entry_id"], str(coordinate_path))
        models = list(parsed.get_models())
        validate_model_chain(source["entry_id"], models, int(source["model"]), source["chain"])
    out.mkdir(parents=True)
    write_json(out / "INPUTS_CHECKED.json", {"case_id": lock["case_id"], "files": checked})
    results = []
    for source in structures:
        coord_path = inputs / source["path"]
        provenance = {
            "class": "experimental",
            "method": source["method"],
            "source_id": f"RCSB PDB:{source['entry_id']}",
            "resolution_angstrom": source["resolution_angstrom"],
        }
        module_result = analyze(
            coord_path,
            subject_id=source["uniprot_accession"],
            provenance=provenance,
            reference_sequence=ref_sequence,
            reference_id=ref_id,
            model_index=source["model"],
            chain=source["chain"],
        )
        entry_out = out / "structqc" / source["entry_id"]
        write_outputs(entry_out, module_result)
        comp = module_result["completeness"]
        seen = {
            int(row["sequence_index"])
            for row in module_result["residues"]
            if row.get("sequence_index") is not None
        }
        missing = sorted(set(range(1, int(comp["reference_length"]) + 1)) - seen)
        b_values = [float(row["mean_b_factor"]) for row in module_result["residues"] if row.get("mean_b_factor") is not None]
        b_summary = {
            "n_residues": len(b_values),
            "minimum": round(min(b_values), 6) if b_values else None,
            "median": round(statistics.median(b_values), 6) if b_values else None,
            "maximum": round(max(b_values), 6) if b_values else None,
        }
        results.append({
            "entry_id": source["entry_id"],
            "entry_url": source["entry_url"],
            "title": source["title"],
            "ligand_context": source["ligand_context"],
            "method": source["method"],
            "resolution_angstrom": source["resolution_angstrom"],
            "selected_model": source["model"],
            "selected_chain": source["chain"],
            "structure_sha256": module_result["coordinate"]["sha256"],
            "reference_id": ref_id,
            "reference_length": comp["reference_length"],
            "coordinate_residues": comp["coordinate_residues"],
            "mapped_residues": comp["mapped_residues"],
            "identity_fraction": comp["identity_fraction"],
            "coverage_fraction": comp["coverage_fraction"],
            "missing_reference_positions": missing,
            "per_residue_mean_b_factor": b_summary,
            "external_validation_state": module_result["external_validation"]["state"],
            "pae_state": module_result["pae"]["state"],
            "warnings": module_result["warnings"],
        })

    module_file = Path(sys.modules["structqc.core"].__file__).resolve()
    distribution = distribution_identity("structqc", module_file)
    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        **distribution,
        "structqc_schema_version": "1.1",
        "structqc_core_sha256": sha256(module_file),
        "case_driver_sha256": sha256(Path(__file__).resolve()),
        "acquisition_script_sha256": sha256(ACQUISITION_PATH),
        "biopython": Bio.__version__,
        "numpy": numpy.__version__,
    }
    comparison = {
        "schema_version": "1.0",
        "case_id": lock["case_id"],
        "protocol_sha256": sha256(PROTOCOL_PATH),
        "source_lock_sha256": sha256(LOCK_PATH),
        "environment": environment,
        "results": results,
        "claim_limits": lock["claim_limits"],
    }
    write_json(out / "CASE_RESULTS.json", comparison)
    write_json(out / "RUN_ENVIRONMENT.json", environment)
    report = [
        "# Public CA II structure coverage comparison",
        "",
        "This report describes two hash-locked public coordinate files using StructQC. It does not test ligand effects, activity, affinity, native conformation, or which model is better.",
        "",
        "## Results",
        "",
        "| PDB entry | Context | Resolution (Å) | Resolved residues | Identity | Reference coverage | Missing P00918 positions | Median per-residue mean B factor |",
        "|---|---|---:|---:|---:|---:|---|---:|",
    ]
    for row in results:
        missing_text = ", ".join(map(str, row["missing_reference_positions"])) or "none"
        report.append(
            f"| [{row['entry_id']}]({row['entry_url']}) | {row['ligand_context']} | "
            f"{row['resolution_angstrom']:.2f} | {row['coordinate_residues']} | "
            f"{row['identity_fraction']:.3f} | {row['coverage_fraction']:.3f} | "
            f"{missing_text} | {row['per_residue_mean_b_factor']['median']:.3f} |"
        )
    report += [
        "",
        "## Interpretation limits",
        "",
        "Coverage and identity are computed against UniProt P00918 by the StructQC sequence mapping. B factors are shown only as deposited coordinate descriptors; different resolution, crystal environment, refinement, and ligand context prevent treating their difference as a controlled quality or ligand effect. No external geometry validation report or PAE was supplied, so these remain missing or unevaluated in the module evidence.",
        "",
        "The exact input and software identities are in `INPUTS_CHECKED.json`, `RUN_ENVIRONMENT.json`, and each `structqc/<PDB>/RUN_MANIFEST.json`. The complete machine-readable comparison is `CASE_RESULTS.json`.",
        "",
        "No participant, private use case, or independent researcher was involved. This is an example workflow, not evidence of independent use.",
        "",
    ]
    (out / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"public_research_case: {exc}", file=sys.stderr)
        raise SystemExit(2)
