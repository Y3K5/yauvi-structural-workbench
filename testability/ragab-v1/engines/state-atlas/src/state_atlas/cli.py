from __future__ import annotations

import argparse
from importlib.resources import files
import json
from pathlib import Path
import sys

from .core import InputError, analyze, read_json, sha256, validate_reference_set, write_outputs


def describe() -> dict:
    return {
        "module_id": "conformational_state", "component_id": "state-atlas", "package": "state_atlas", "version": "0.2.0",
        "display_name": "Conformational State and MD Ensemble Atlas",
        "one_line": "Which experimentally bounded conformational states do these structures or trajectory frames resemble?",
        "commands": ["describe", "validate", "fetch", "run"],
        "inputs": [{"name": "structure_evidence_manifest", "contract": "structure_evidence_manifest", "required": True},
                   {"name": "reference_set", "format": "JSON", "required": True},
                   {"name": "alignment_map", "format": "JSON", "required": False,
                    "note": "required by Reference Set v2 and all Mark 1 ABL-family calls"},
                   {"name": "structure_or_trajectory", "format": "PDB/mmCIF or XTC/DCD/TRR plus topology", "required": True}],
        "outputs": [{"name": "STATE_ENSEMBLE.json", "contract": "state_ensemble_summary", "format": "json"},
                    {"name": "FRAME_METRICS.tsv", "contract": "state_frame_table", "format": "tsv"},
                    {"name": "CLUSTERS.tsv", "contract": "state_cluster_table", "format": "tsv"},
                    {"name": "STATE_LAYER.json", "contract": "structure_layer_bundle", "format": "json"},
                    {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest", "format": "json"}],
        "runtimes": ["yauvi-python"], "optional_runtimes": ["mdanalysis"], "deterministic_output": True,
        "limitations": ["Labels are structural resemblance, not activity.", "Populations depend on trajectory preparation.",
                        "Unresolved frames remain visible.",
                        "Mark 1 qualification is limited to checksum-bound ABL-family Reference Set v2 analyses."],
        "scientific_scopes": [
            {"scope_id": "abl_family", "scientific_state": "prototype",
             "benchmark_collection": "yauvi-structural-public-qualification-v2", "release_blocking": True,
             "supported_subject_class": "ABL-family conformational resemblance",
             "known_limitations": ["Becomes conditionally qualified only after the v2 held-out gate passes."]},
            {"scope_id": "other_proteins", "scientific_state": "prototype",
             "benchmark_collection": "not_qualified", "release_blocking": False,
             "supported_subject_class": "user-declared two-sided reference families"},
        ],
    }


def _inputs(p):
    p.add_argument("--manifest", required=True); p.add_argument("--reference-set", required=True)
    p.add_argument("--alignment-map")
    p.add_argument("--structure"); p.add_argument("--topology"); p.add_argument("--trajectory")
    p.add_argument("--chain"); p.add_argument("--selection", default="protein and name CA")
    p.add_argument("--stride", type=int, default=1); p.add_argument("--pbc", choices=["none", "unwrap"])
    p.add_argument("--cluster-cutoff-A", type=float, default=2.0); p.add_argument("--collective-variables")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="state-atlas"); sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("describe"); fetch = sub.add_parser("fetch"); fetch.add_argument("--plan", action="store_true", required=True)
    validate = sub.add_parser("validate"); validate.add_argument("--reference-set", required=True); validate.add_argument("--alignment-map")
    run = sub.add_parser("run"); _inputs(run); run.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.command == "describe": print(json.dumps(describe(), indent=2, sort_keys=True)); return 0
        if args.command == "fetch": print((files("state_atlas") / "sources.yaml").read_text(), end=""); return 0
        reference_path = Path(args.reference_set); refs = read_json(reference_path)
        alignment_path = Path(args.alignment_map) if getattr(args, "alignment_map", None) else None
        alignment = read_json(alignment_path) if alignment_path else None
        alignment_digest = sha256(alignment_path) if alignment_path else None
        if args.command == "validate":
            errors = validate_reference_set(
                refs, base=reference_path.parent, alignment_map=alignment,
                alignment_map_sha256=alignment_digest,
            )
            print(json.dumps({"valid": not errors, "errors": errors}, indent=2, sort_keys=True)); return 0 if not errors else 1
        doc = analyze(read_json(args.manifest), refs, reference_base=reference_path.parent,
                      structure_path=args.structure, topology_path=args.topology, trajectory_path=args.trajectory,
                      chain=args.chain, selection=args.selection, stride=args.stride, pbc=args.pbc,
                      cluster_cutoff_A=args.cluster_cutoff_A,
                      collective_variables=(read_json(args.collective_variables) or [] if args.collective_variables else []),
                      alignment_map=alignment, alignment_map_sha256=alignment_digest)
        write_outputs(args.out, doc)
        return 1 if doc["overall_label"] == "unresolved" else 0
    except (InputError, OSError, ValueError) as exc:
        print(f"state-atlas: {exc}", file=sys.stderr); return 2


if __name__ == "__main__": raise SystemExit(main())
