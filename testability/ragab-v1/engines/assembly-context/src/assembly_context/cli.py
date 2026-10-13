from __future__ import annotations

import argparse
from importlib.resources import files
import json
import sys

from .core import InputError, analyze, read_json, write_outputs


def describe() -> dict:
    return {
        "module_id": "assembly_context", "component_id": "assembly-context", "package": "assembly_context",
        "version": "0.1.0", "display_name": "Biological Assembly and Interface Context",
        "one_line": "Which residues are contacted or buried when an isolated chain is placed in a declared assembly?",
        "commands": ["describe", "validate", "fetch", "run"],
        "inputs": [{"name": "structure_evidence_manifest", "contract": "structure_evidence_manifest", "required": True},
                   {"name": "isolated", "format": "PDB or mmCIF", "required": True},
                   {"name": "assembly", "format": "expanded PDB or mmCIF assembly", "required": True}],
        "outputs": [{"name": "ASSEMBLY_CONTEXT.json", "contract": "assembly_context_summary", "format": "json"},
                    {"name": "INTERFACES.tsv", "contract": "assembly_interface_table", "format": "tsv"},
                    {"name": "ASSEMBLY_LAYER.json", "contract": "structure_layer_bundle", "format": "json"},
                    {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest", "format": "json"}],
        "runtimes": ["yauvi-python"], "optional_runtimes": ["freesasa", "gemmi"],
        "deterministic_output": True,
        "limitations": ["One coordinate state is not native exposure.", "Homolog assemblies transfer architecture only.",
                        "Partial assemblies yield lower bounds."],
    }


def _inputs(parser):
    parser.add_argument("--manifest", required=True); parser.add_argument("--isolated", required=True)
    parser.add_argument("--assembly", required=True); parser.add_argument("--subject-chain", required=True)
    parser.add_argument("--relationship", required=True,
                        choices=["exact_protein", "homolog_assembly", "architecture_analogy", "unresolved"])
    parser.add_argument("--reference-id", default=""); parser.add_argument("--assembly-id")
    parser.add_argument("--expected-chains", default="", help="comma-separated chain ids")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="assembly-context"); sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("describe"); fetch = sub.add_parser("fetch"); fetch.add_argument("--plan", action="store_true", required=True)
    validate = sub.add_parser("validate"); _inputs(validate)
    run = sub.add_parser("run"); _inputs(run); run.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.command == "describe": print(json.dumps(describe(), indent=2, sort_keys=True)); return 0
        if args.command == "fetch": print((files("assembly_context") / "sources.yaml").read_text(), end=""); return 0
        document = analyze(read_json(args.manifest), args.isolated, args.assembly,
                           subject_chain=args.subject_chain, relationship=args.relationship,
                           reference_id=args.reference_id, assembly_id=args.assembly_id,
                           expected_chains=[x for x in args.expected_chains.split(",") if x])
        if args.command == "validate":
            print(json.dumps({"valid": True, "chains": document["assembly"]["chains_observed"]}, sort_keys=True))
        else: write_outputs(args.out, document)
        return 1 if document["assembly"]["lower_bound"] or args.relationship == "unresolved" else 0
    except (InputError, OSError, ValueError) as exc:
        print(f"assembly-context: {exc}", file=sys.stderr); return 2


if __name__ == "__main__": raise SystemExit(main())
