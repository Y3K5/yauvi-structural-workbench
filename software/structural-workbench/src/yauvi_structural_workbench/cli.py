"""Structural-only YAUVI command line used by the JOSS reviewer distribution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

from yauvi_platform.structural_workbench import (
    AnalysisError, StructuralAnalysisStore, analysis_definitions,
)
from . import __version__


def _workspace(value: str | None) -> Path:
    if value:
        return Path(value).resolve()
    return Path.cwd().resolve()


def _print(value: Any) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="yauvi",
        description="Local, evidence-bounded structural protein analysis.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--workspace", help="Local analysis workspace (default: current directory).")
    groups = parser.add_subparsers(dest="group", required=True)

    analysis = groups.add_parser(
        "analysis", help="Create, validate, run, and export structural analyses.",
        description=(
            "Create and run local structural analyses. Use 'analysis types' to list "
            "analysis types and 'analysis inputs --type TYPE' to see allowed input "
            "roles, file formats, and requirements."
        ),
    )
    actions = analysis.add_subparsers(dest="action", required=True)
    definitions = analysis_definitions()
    type_ids = [item["analysis_type"] for item in definitions]
    types = actions.add_parser("types", help="List supported analysis types and their claim limits.")
    inputs = actions.add_parser("inputs", help="Show allowed input roles and formats for an analysis type.")
    inputs.add_argument("--type", required=True, choices=type_ids, help="Analysis type from 'analysis types'.")
    create = actions.add_parser("create", help="Create an analysis case.")
    create.add_argument("--analysis", required=True, help="Local case identifier.")
    create.add_argument("--type", required=True, choices=type_ids, help="Analysis type from 'analysis types'.")
    create.add_argument("--question", required=True)
    create.add_argument("--subject-id", default="")
    add = actions.add_parser("add", help="Attach one local input file to an analysis.")
    add.add_argument("--analysis", required=True)
    add.add_argument("--role", required=True)
    add.add_argument("--file", required=True)
    add.description = "Use 'analysis inputs --type TYPE' to see allowed roles and file formats."
    for action in ("validate", "run"):
        command = actions.add_parser(action)
        command.add_argument("--analysis", required=True)
    export = actions.add_parser("export")
    export.add_argument("--analysis", required=True)
    export.add_argument("--out", required=True)

    params = actions.add_parser("parameters", help="Create a revision from a JSON parameter file.")
    params.add_argument("--analysis", required=True)
    params.add_argument("--file", required=True)
    example = groups.add_parser("example", help="Create a bundled synthetic StructQC analysis offline.")
    example.add_argument("--analysis", default="qc-example")
    example.add_argument("--without-validation", action="store_true")

    workbench = groups.add_parser("workbench", help="Serve the loopback-only browser workbench.")
    actions = workbench.add_subparsers(dest="action", required=True)
    for action in ("serve", "open"):
        command = actions.add_parser(action)
        command.add_argument("--host", default="127.0.0.1")
        command.add_argument("--port", default=8947, type=int)
        command.add_argument("--allow-reference-fetch", action="store_true")
        command.add_argument("--label", default="Local workspace", help="Name shown for this local workspace.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    workspace = _workspace(args.workspace)
    try:
        if args.group == "example":
            from .example import create_example
            _print(create_example(workspace, args.analysis, without_validation=args.without_validation))
            return 0
        if args.group == "analysis":
            definitions = analysis_definitions()
            by_type = {item["analysis_type"]: item for item in definitions}
            if args.action == "types":
                _print([
                    {
                        "analysis_type": item["analysis_type"],
                        "title": item["title"],
                        "question": item["question"],
                        "readiness": item["readiness"],
                        "claim_ceiling": item["claim_ceiling"],
                    }
                    for item in definitions
                ])
                return 0
            if args.action == "inputs":
                definition = by_type[args.type]
                _print({
                    "analysis_type": definition["analysis_type"],
                    "inputs": [
                        {
                            key: role[key]
                            for key in (
                                "role", "label", "required", "multiple", "minimum_files",
                                "accepted_artifact_types", "accepted_extensions", "format_guide",
                                "source_ids", "absence_effect", "sensitivity",
                            )
                            if key in role
                        }
                        for role in definition["inputs"]
                    ],
                })
                return 0
            store = StructuralAnalysisStore(workspace)
            if args.action == "create":
                _print(store.create(
                    args.analysis, analysis_type=args.type,
                    question=args.question, subject_id=args.subject_id,
                ))
                return 0
            if args.action == "add":
                _print(store.add_file(args.analysis, role=args.role, path=args.file))
                return 0
            if args.action == "parameters":
                _print(store.update_parameters(args.analysis, json.loads(Path(args.file).read_text())))
                return 0
            if args.action == "validate":
                result = store.preflight(args.analysis)
                _print(result)
                return 0 if result["valid"] else 1
            if args.action == "run":
                result = store.run(args.analysis)
                _print(result)
                if result["status"] == "completed":
                    return 0
                return 1 if result["status"] in {"scientifically_incomplete", "blocked"} else 2
            _print(store.export(args.analysis, args.out))
            return 0

        if args.host not in {"127.0.0.1", "localhost", "::1"}:
            raise AnalysisError("YAUVI is local-only; --host must be a loopback address")
        if args.action == "open":
            from .launcher import open_workbench
            return open_workbench(workspace, args.host, args.port, args.allow_reference_fetch, label=args.label)
        from .server import serve
        return serve(workspace, args.host, args.port, args.allow_reference_fetch, label=args.label)

    except (AnalysisError, OSError, ValueError) as exc:
        print(f"YAUVI blocked: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
