from __future__ import annotations

import argparse
from importlib.resources import files
import json
from pathlib import Path
import subprocess
import sys

from .core import InputError, ScientificBlock, execute_job, prepare_job, read_json, sha256, write_job


def describe() -> dict:
    return {
        "module_id": "docking_preparation", "component_id": "dockprep", "package": "dockprep", "version": "0.1.0",
        "display_name": "Docking Preparation and Execution", "one_line": "Can a checksum-approved receptor and ligand be prepared and run locally without exceeding the evidence?",
        "commands": ["describe", "validate", "fetch", "run", "execute"],
        "inputs": [{"name": "job_spec", "format": "json", "required": True},
                   {"name": "protein_case_manifest", "contract": "protein_case_manifest", "required": False},
                   {"name": "residue_identity_map", "contract": "residue_identity_map", "required": False}],
        "outputs": [{"name": "DOCKING_JOB.json", "contract": "docking_job_bundle", "format": "json"},
                    {"name": "DOCKING_EXECUTION.json", "contract": "docking_execution_manifest", "format": "json"},
                    {"name": "RUN_MANIFEST.json", "contract": "structure_analysis_run_manifest", "format": "json"}],
        "runtimes": ["yauvi-python"], "optional_runtimes": ["autodock-vina", "meeko", "haddock3"], "deterministic_output": True,
        "limitations": ["Runtime absence blocks only that docking leg.", "Docking output is pose evidence, not affinity or mechanism.",
                        "Chemistry, source terms, exact mappings, and human approval are hard gates."],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="dockprep")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("describe")
    fetch = sub.add_parser("fetch"); fetch.add_argument("--plan", action="store_true", required=True)
    for name in ("validate", "run"):
        command = sub.add_parser(name); command.add_argument("--config", required=True)
        if name == "run": command.add_argument("--out", required=True)
    execute = sub.add_parser("execute"); execute.add_argument("--job", required=True); execute.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "describe": print(json.dumps(describe(), indent=2, sort_keys=True)); return 0
        if args.command == "fetch": print((files("dockprep") / "sources.yaml").read_text(), end=""); return 0
        if args.command == "execute":
            _manifest, code = execute_job(args.job, args.out)
            return code
        spec_path = Path(args.config).resolve()
        if args.command == "validate":
            # Validation uses a private scratch directory beside no scientific output.
            import tempfile
            with tempfile.TemporaryDirectory() as temporary:
                document = prepare_job(read_json(spec_path), base_dir=spec_path.parent, out_dir=temporary)
            print(json.dumps({"valid": True, "job_id": document["job_id"], "missing_evidence": document["missing_evidence"]}, sort_keys=True))
            return 1 if document["missing_evidence"] else 0
        document = prepare_job(read_json(spec_path), base_dir=spec_path.parent, out_dir=args.out)
        write_job(args.out, document, config_sha256=sha256(spec_path))
        return 1 if document["missing_evidence"] else 0
    except ScientificBlock as exc:
        print(f"dockprep: scientifically blocked: {exc}", file=sys.stderr); return 1
    except (InputError, OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"dockprep: {exc}", file=sys.stderr); return 2


if __name__ == "__main__": raise SystemExit(main())
