from __future__ import annotations

import argparse
from importlib.resources import files
import json
from pathlib import Path
import sys

from .core import InputError, analyze, read_fasta, read_json, read_validation_report, write_outputs


def describe() -> dict:
    return {
        "module_id": "structure_quality", "component_id": "structqc", "package": "structqc",
        "version": "0.1.0", "display_name": "Structure Provenance and Quality",
        "one_line": "Are these coordinates identity-bound, provenance-declared, and structurally inspectable?",
        "commands": ["describe", "validate", "fetch", "run"],
        "inputs": [{"name": "structure", "format": "PDB or mmCIF", "required": True},
                   {"name": "reference_fasta", "format": "FASTA", "required": False},
                   {"name": "provenance", "format": "JSON", "required": False},
                   {"name": "pae", "format": "JSON", "required": False},
                   {"name": "validation_report", "format": "wwPDB XML or validator JSON", "required": False}],
        "outputs": [
            {"name": "STRUCTURE_EVIDENCE.json", "contract": "structure_evidence_manifest", "format": "json"},
            {"name": "RESIDUE_QUALITY.tsv", "contract": "residue_quality_table", "format": "tsv"},
            {"name": "STRUCTURE_LAYER.json", "contract": "structure_layer_bundle", "format": "json"},
            {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest", "format": "json"},
        ],
        "runtimes": ["yauvi-python"], "optional_runtimes": ["gemmi", "MolProbity", "Phenix"],
        "deterministic_output": True,
        "limitations": [
            "Coordinate quality does not establish native conformation.",
            "Unknown provenance remains unknown.",
            "Completeness is unevaluated without a reference sequence.",
        ],
    }


def _analysis(args):
    reference_id, sequence = read_fasta(args.reference_fasta)
    return analyze(
        args.structure, subject_id=args.subject_id, provenance=read_json(args.provenance),
        reference_sequence=sequence, reference_id=reference_id, pae=read_json(args.pae),
        validation_report=read_validation_report(args.validation_report),
        model_index=args.model, chain=args.chain,
    )


def _add_inputs(parser):
    parser.add_argument("--structure", required=True)
    parser.add_argument("--subject-id")
    parser.add_argument("--reference-fasta")
    parser.add_argument("--provenance")
    parser.add_argument("--pae")
    parser.add_argument("--validation-report")
    parser.add_argument("--require-external-validation", action="store_true")
    parser.add_argument("--model", type=int, default=0)
    parser.add_argument("--chain")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="structqc")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("describe")
    validate = sub.add_parser("validate")
    _add_inputs(validate)
    fetch = sub.add_parser("fetch")
    fetch.add_argument("--plan", action="store_true", required=True)
    run = sub.add_parser("run")
    _add_inputs(run)
    run.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "describe":
            print(json.dumps(describe(), indent=2, sort_keys=True))
            return 0
        if args.command == "fetch":
            print((files("structqc") / "sources.yaml").read_text(encoding="utf-8"), end="")
            return 0
        document = _analysis(args)
        if args.command == "validate":
            print(json.dumps({"valid": True, "subject": document["subject"], "warnings": document["warnings"]}, sort_keys=True))
        else:
            write_outputs(args.out, document)
        incomplete = document["provenance"]["class"] == "unknown"
        if args.require_external_validation and document["external_validation"]["state"] != "imported":
            incomplete = True
        return 1 if incomplete else 0
    except (InputError, OSError, ValueError) as exc:
        print(f"structqc: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
