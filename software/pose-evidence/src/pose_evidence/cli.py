from __future__ import annotations

import argparse
from importlib.resources import files
import json
from pathlib import Path
import sys

from .core import InputError, analyze, read_json, sha256, write_outputs


def describe() -> dict:
    return {"module_id": "pose_evidence", "component_id": "pose-evidence", "package": "pose_evidence", "version": "0.1.0",
            "display_name": "Docking Pose Evidence", "one_line": "What geometry, clustering, consistency, and unresolved evidence is present across imported poses?",
            "commands": ["describe", "validate", "fetch", "run"],
            "inputs": [{"name": "pose_config", "format": "json", "required": True}, {"name": "docking_job_bundle", "contract": "docking_job_bundle", "required": False}],
            "outputs": [{"name": "POSE_EVIDENCE.json", "contract": "pose_evidence_summary", "format": "json"},
                        {"name": "POSE_CLUSTERS.tsv", "contract": "pose_cluster_table", "format": "tsv"},
                        {"name": "POSE_CONTACTS.tsv", "contract": "pose_contact_table", "format": "tsv"},
                        {"name": "POSE_LAYER.json", "contract": "structure_layer_bundle", "format": "json"},
                        {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest", "format": "json"}],
            "runtimes": ["yauvi-python"], "optional_runtimes": ["freesasa"], "deterministic_output": True,
            "limitations": ["Scores and units remain engine-specific.", "Pose evidence is not affinity, engagement, activity, or mechanism.",
                            "Hydrogen bonds, salt bridges, metal contacts, and clashes are geometry proxies."]}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="pose-evidence"); sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("describe"); fetch = sub.add_parser("fetch"); fetch.add_argument("--plan", action="store_true", required=True)
    for name in ("validate", "run"):
        command = sub.add_parser(name); command.add_argument("--config", required=True)
        if name == "run": command.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "describe": print(json.dumps(describe(), indent=2, sort_keys=True)); return 0
        if args.command == "fetch": print((files("pose_evidence") / "sources.yaml").read_text(), end=""); return 0
        config_path = Path(args.config).resolve(); document, digests = analyze(read_json(config_path), base_dir=config_path.parent)
        if args.command == "validate": print(json.dumps({"valid": True, "poses": len(document["poses"]), "clusters": len({r['cluster_id'] for r in document['clusters']})}, sort_keys=True))
        else: write_outputs(args.out, document, input_sha256=digests, config_sha256=sha256(config_path))
        return 1 if document["buried_sasa"]["state"] == "unevaluated" else 0
    except (InputError, OSError, ValueError) as exc:
        print(f"pose-evidence: {exc}", file=sys.stderr); return 2


if __name__ == "__main__": raise SystemExit(main())
