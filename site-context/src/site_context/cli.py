from __future__ import annotations

import argparse
from importlib.resources import files
import json
import sys

from .core import InputError, analyze, read_annotations, read_json, write_outputs


def describe() -> dict:
    return {
        "module_id": "site_context", "component_id": "site-context", "package": "site_context", "version": "0.1.0",
        "display_name": "Functional Sites, Cofactors, and Pockets",
        "one_line": "Where do declared functional residues, observed cofactors, and predicted pockets lie in this structure?",
        "commands": ["describe", "validate", "fetch", "run"],
        "inputs": [{"name": "structure_evidence_manifest", "contract": "structure_evidence_manifest", "required": True},
                   {"name": "annotations", "format": "JSON or TSV", "required": True},
                   {"name": "pocket_result", "format": "method-declared JSON", "required": False}],
        "outputs": [{"name": "SITE_CONTEXT.json", "contract": "site_context_summary", "format": "json"},
                    {"name": "SITE_RESIDUES.tsv", "contract": "site_residue_table", "format": "tsv"},
                    {"name": "POCKETS.tsv", "contract": "pocket_evidence_table", "format": "tsv"},
                    {"name": "SITE_LAYER.json", "contract": "structure_layer_bundle", "format": "json"},
                    {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest", "format": "json"}],
        "runtimes": ["yauvi-python"], "optional_runtimes": ["fpocket", "p2rank"], "deterministic_output": True,
        "limitations": ["Evidence legs remain separate.", "Pocket scores are method-specific.",
                        "Proximity does not establish function or affinity."],
    }


def _inputs(p):
    p.add_argument("--manifest", required=True); p.add_argument("--structure", required=True)
    p.add_argument("--annotations", required=True); p.add_argument("--component-map")
    p.add_argument("--pocket-result", action="append", default=[])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="site-context"); sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("describe"); fetch = sub.add_parser("fetch"); fetch.add_argument("--plan", action="store_true", required=True)
    validate = sub.add_parser("validate"); _inputs(validate)
    run = sub.add_parser("run"); _inputs(run); run.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.command == "describe": print(json.dumps(describe(), indent=2, sort_keys=True)); return 0
        if args.command == "fetch": print((files("site_context") / "sources.yaml").read_text(), end=""); return 0
        doc = analyze(read_json(args.manifest), args.structure, read_annotations(args.annotations),
                      component_map=read_json(args.component_map), pocket_results=[read_json(p) for p in args.pocket_result])
        if args.command == "validate": print(json.dumps({"valid": True, "sites": len(doc["sites"])}, sort_keys=True))
        else: write_outputs(args.out, doc)
        unresolved = any(s["state"].startswith("unresolved") or s["state"] == "missing_coordinates" for s in doc["sites"])
        return 1 if unresolved or doc["missing_evidence"] else 0
    except (InputError, OSError, ValueError) as exc:
        print(f"site-context: {exc}", file=sys.stderr); return 2


if __name__ == "__main__": raise SystemExit(main())
