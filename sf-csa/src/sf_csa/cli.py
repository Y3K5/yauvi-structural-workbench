#!/usr/bin/env python3
"""`sf-csa` — structure/function comparison across species.

    sf-csa build-manifests --spec campaign.json --out config/
    sf-csa run       --queries q.json --databases db.json --output results/
    sf-csa verify    --output results/ --databases db.json
    sf-csa describe                     machine-readable IO contract
    sf-csa validate  --queries q.json --databases db.json
    sf-csa fetch --plan                 what raw files are needed, and where from

`run` and `verify` keep their original flags. The rest are the uniform module
contract every package in the platform answers.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

from .core import SFCSError, run_pipeline, verify_release

EXIT_OK = 0
EXIT_FAILED = 1
EXIT_BLOCKED = 2


def cmd_run(args) -> int:
    manifest = run_pipeline(Path(args.queries), Path(args.databases), Path(args.output))
    print(
        f"SF-CSA release {manifest['release_id']}: "
        f"{manifest['query_count']} queries, {manifest['proteome_count']} proteomes"
    )
    return EXIT_OK


def cmd_verify(args) -> int:
    errors = verify_release(Path(args.output), Path(args.databases))
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return EXIT_FAILED
    print("SF-CSA release verification: PASS")
    return EXIT_OK


def cmd_build_manifests(args) -> int:
    from .manifests import SpecError, build

    try:
        target_path, database_path = build(args.spec, args.out)
    except SpecError as exc:
        print(f"campaign spec error: {exc}", file=sys.stderr)
        return EXIT_BLOCKED
    print(f"wrote {target_path}")
    print(f"wrote {database_path}")
    return EXIT_OK


def cmd_describe(args) -> int:
    from .module_contract import describe

    print(json.dumps(describe(), indent=2, sort_keys=True))
    return EXIT_OK


def cmd_validate(args) -> int:
    from .module_contract import validate_manifests

    checks = validate_manifests(args.queries, args.databases)
    failed = 0
    for name, ok, detail in checks:
        if ok:
            print(f"  [ok]      {name}" + (f"  {detail}" if detail else ""))
        else:
            failed += 1
            print(f"  [FAILED]  {name}  {detail}")
    return EXIT_OK if not failed else EXIT_FAILED


def cmd_fetch(args) -> int:
    from .module_contract import SOURCES_MANIFEST

    if not SOURCES_MANIFEST.is_file():
        print("no source manifest shipped with this install", file=sys.stderr)
        return EXIT_BLOCKED
    if args.plan_only:
        print(SOURCES_MANIFEST.read_text(encoding="utf-8"))
        return EXIT_OK
    try:
        from yauvi_sources.cli import main as fetch_main
    except ImportError:
        print(
            "acquisition is provided by yauvi-sources, which is not installed.\n"
            "  pip install 'yauvi-sources[fetch]'\n"
            f"  then: yauvi-fetch plan --for sf_csa --manifest {SOURCES_MANIFEST}\n"
            "Or run `sf-csa fetch --plan` to print the declared sources.",
            file=sys.stderr,
        )
        return EXIT_FAILED
    return fetch_main(["plan", "--for", "sf_csa", "--manifest", str(SOURCES_MANIFEST)])


def build_parser() -> argparse.ArgumentParser:
    from . import __version__

    parser = argparse.ArgumentParser(
        prog="sf-csa", description="Structure-Function Comparative Species Analysis"
    )
    parser.add_argument("--version", action="version", version=f"sf-csa {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the comparison and write a release")
    run.add_argument("--queries", required=True)
    run.add_argument("--databases", required=True)
    run.add_argument("--output", required=True)
    run.set_defaults(func=cmd_run)

    verify = sub.add_parser("verify", help="audit an existing release")
    verify.add_argument("--output", required=True)
    # The release is audited against the manifest it was configured with, not
    # against the shape it recorded for itself.
    verify.add_argument("--databases", required=True)
    verify.set_defaults(func=cmd_verify)

    build_cmd = sub.add_parser(
        "build-manifests", help="build checksum-pinned manifests from a campaign spec"
    )
    build_cmd.add_argument("--spec", required=True, help="campaign spec JSON")
    build_cmd.add_argument("--out", required=True, help="directory to write the manifests into")
    build_cmd.set_defaults(func=cmd_build_manifests)

    describe = sub.add_parser("describe", help="print the machine-readable IO contract")
    describe.set_defaults(func=cmd_describe)

    validate = sub.add_parser("validate", help="check the manifests without running")
    validate.add_argument("--queries", required=True)
    validate.add_argument("--databases", required=True)
    validate.set_defaults(func=cmd_validate)

    fetch = sub.add_parser("fetch", help="what raw files this module needs, and where from")
    fetch.add_argument(
        "--plan", dest="plan_only", action="store_true", help="print the declared sources only"
    )
    fetch.set_defaults(func=cmd_fetch)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except SFCSError as exc:
        # Fail-closed is the module's contract: a precondition that is not met
        # blocks, and says so, rather than producing a partial release.
        print(f"SF-CSA fail-closed: {exc}", file=sys.stderr)
        return EXIT_BLOCKED
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
