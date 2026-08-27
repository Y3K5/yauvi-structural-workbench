#!/usr/bin/env python3
"""Run the public/reviewer Structural Workbench tests in isolated processes.

The repository contains legacy test modules with duplicate basenames, so one
large pytest collection is not reliable.  This runner mirrors the reviewer CI,
keeps network and optional-adapter tests off by default, and reports aggregate
counts without including unrelated private YAUVI platform tests.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SUITES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("structqc", ("structqc/tests",)),
    ("memorient", ("Membrane Orientor/memorient/tests",)),
    ("state-atlas", ("state-atlas/tests",)),
    ("site-context", ("site-context/tests",)),
    ("activity-state", ("activity-state/tests",)),
    ("assembly-context", ("assembly-context/tests",)),
    ("sf-csa", ("sf-csa/tests",)),
    ("sf-csa fixture", ("tools/fixtures/sfcsa",)),
    ("structural workbench", ("platform/tests/test_structural_workbench.py",)),
    ("qualification v2", ("platform/tests/test_qualification_v2.py",)),
    ("source registry", ("sources/tests",)),
)
COUNT = re.compile(r"(?P<count>\d+) (?P<kind>passed|failed|skipped|deselected|error|errors)\b")


def parse_counts(output: str) -> dict[str, int]:
    counts = {"passed": 0, "failed": 0, "skipped": 0, "deselected": 0, "errors": 0}
    for match in COUNT.finditer(output):
        kind = match.group("kind")
        if kind == "error":
            kind = "errors"
        counts[kind] = max(counts[kind], int(match.group("count")))
    return counts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path, help="Optional path for the machine-readable summary.")
    parser.add_argument("--verbose", action="store_true", help="Show each complete pytest output.")
    parser.add_argument(
        "--core-only", action="store_true",
        help="Accepted for compatibility. Every suite in this distribution is already core.",
    )
    args = parser.parse_args(argv)

    records: list[dict[str, Any]] = []
    totals = {"passed": 0, "failed": 0, "skipped": 0, "deselected": 0, "errors": 0}
    failed = False
    for label, paths in SUITES:
        command = [
            sys.executable, "-m", "pytest", *paths, "-q",
            "-m", "not network and not adapter", "--tb=short",
        ]
        run = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
        output = run.stdout + run.stderr
        counts = parse_counts(output)
        for key in totals:
            totals[key] += counts[key]
        record = {
            "suite": label,
            "paths": list(paths),
            "exit_code": run.returncode,
            **counts,
        }
        records.append(record)
        failed = failed or run.returncode != 0
        state = "PASS" if run.returncode == 0 else "FAIL"
        detail = ", ".join(f"{value} {key}" for key, value in counts.items() if value)
        print(f"[{state}] {label}: {detail or 'no test count parsed'}")
        if args.verbose or run.returncode != 0:
            print(output.rstrip())

    summary = {
        "schema_version": "1.0",
        "python": sys.version.split()[0],
        "selection": "not network and not adapter",
        "scope": (
            "JOSS reviewer distribution"
        ),
        "totals": totals,
        "suites": records,
        "passed": not failed,
        "scientific_boundary": (
            "Software test success does not establish external scientific qualification, "
            "biochemical activity, native exposure, binding, safety, or efficacy."
        ),
    }
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"passed": summary["passed"], "totals": totals}, sort_keys=True))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
