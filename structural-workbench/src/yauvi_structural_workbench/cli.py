"""Structural-only YAUVI command line used by the JOSS reviewer distribution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

from yauvi_platform.structural_workbench import AnalysisError, StructuralAnalysisStore


def _workspace(value: str | None) -> Path:
    if value:
        return Path(value).resolve()
    current = Path.cwd().resolve()
    for candidate in (current, *current.parents):
        if (candidate / "apps" / "yauvi").is_dir() and (candidate / "catalogs").is_dir():
            return candidate
    return current


def _print(value: Any) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="yauvi",
        description="Local, evidence-bounded structural protein analysis.",
    )
    parser.add_argument("--workspace", help="Structural Workbench repository root.")
    groups = parser.add_subparsers(dest="group", required=True)

    analysis = groups.add_parser("analysis", help="Create, validate, run, and export structural analyses.")
    actions = analysis.add_subparsers(dest="action", required=True)
    create = actions.add_parser("create")
    create.add_argument("--analysis", required=True)
    create.add_argument("--type", required=True)
    create.add_argument("--question", required=True)
    create.add_argument("--subject-id", default="")
    add = actions.add_parser("add")
    add.add_argument("--analysis", required=True)
    add.add_argument("--role", required=True)
    add.add_argument("--file", required=True)
    for action in ("validate", "run"):
        command = actions.add_parser(action)
        command.add_argument("--analysis", required=True)
    export = actions.add_parser("export")
    export.add_argument("--analysis", required=True)
    export.add_argument("--out", required=True)

    workbench = groups.add_parser("workbench", help="Serve the loopback-only browser workbench.")
    serve = workbench.add_subparsers(dest="action", required=True).add_parser("serve")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", default=8931, type=int)
    serve.add_argument("--allow-reference-fetch", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    workspace = _workspace(args.workspace)
    try:
        if args.group == "analysis":
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
        server = workspace / "apps" / "yauvi" / "controller" / "server.py"
        if not server.is_file():
            raise AnalysisError(
                "the browser controller is not present; install from the reviewer repository checkout"
            )
        command = [sys.executable, str(server), "--host", args.host, "--port", str(args.port)]
        if args.allow_reference_fetch:
            command.append("--allow-reference-fetch")
        return subprocess.run(command, cwd=server.parent, check=False).returncode
    except (AnalysisError, OSError, ValueError) as exc:
        print(f"YAUVI blocked: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
