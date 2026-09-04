#!/usr/bin/env python3
"""Run the public/reviewer Structural Workbench tests in isolated processes.

The repository contains legacy test modules with duplicate basenames, so one
large pytest collection is not reliable.  This runner mirrors the reviewer CI,
keeps network and optional-adapter tests off by default, and reports aggregate
counts without including unrelated private YAUVI platform tests.
"""
from __future__ import annotations

import argparse
import importlib.metadata as importlib_metadata
import json
from pathlib import Path
import re
import subprocess
import sys
import tomllib
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

# Every pyproject in the distribution, so a floor is read from the package that
# declares it rather than from a list here that can drift away from them.
PYPROJECTS: tuple[str, ...] = (
    "pyproject.toml", "structqc/pyproject.toml", "Membrane Orientor/memorient/pyproject.toml",
    "state-atlas/pyproject.toml", "site-context/pyproject.toml", "activity-state/pyproject.toml",
    "assembly-context/pyproject.toml", "sf-csa/pyproject.toml",
)
REQUIREMENT = re.compile(r"^\s*([A-Za-z0-9._-]+)\s*>=\s*([0-9][0-9A-Za-z.\-]*)")


def _version_tuple(text: str) -> tuple[int, ...]:
    parts: list[int] = []
    for chunk in text.split("."):
        digits = re.match(r"\d+", chunk)
        if digits is None:
            break
        parts.append(int(digits.group()))
    return tuple(parts)


def declared_floors() -> dict[str, tuple[str, str]]:
    """distribution name -> (floor, the pyproject that asks for it)."""
    floors: dict[str, tuple[str, str]] = {}
    for relative in PYPROJECTS:
        path = ROOT / relative
        if not path.is_file():
            continue
        data = tomllib.loads(path.read_text(encoding="utf-8"))
        for spec in data.get("project", {}).get("dependencies", []):
            match = REQUIREMENT.match(spec)
            if match is None:
                continue
            name, floor = match.group(1).lower(), match.group(2)
            current = floors.get(name)
            if current is None or _version_tuple(floor) > _version_tuple(current[0]):
                floors[name] = (floor, relative)
    return floors


def preflight() -> list[str]:
    """Check the interpreter satisfies what the packages declare they need.

    Without this an out-of-spec environment reports as a scientific failure. On
    2026-09-03 Biopython 1.78 sat under a declared floor of 1.81, and the two
    visible symptoms were `No module named 'Bio.PDB.SASA'` and a site-context
    cofactor asserting `observed_match` and getting `not_observed` -- the second
    of which reads as a result about proteins and is not one. It is a parser
    behaviour: 1.78 leaves `residue.resname` unstripped, so a PDB-padded ` ZN`
    never matches a declared `ZN`.

    Compiled extensions are checked by importing them, because a wheel built for
    the wrong architecture installs without complaint and fails only on import.
    """
    problems: list[str] = []
    for name, (floor, source) in sorted(declared_floors().items()):
        try:
            installed = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            problems.append(f"{name}: not installed, {source} requires >={floor}")
            continue
        if _version_tuple(installed) < _version_tuple(floor):
            problems.append(f"{name}: {installed} installed, {source} requires >={floor}")
    for module, package in (("Bio.PDB.PDBParser", "biopython"), ("Bio.Align", "biopython"),
                            ("numpy", "numpy"), ("scipy", "scipy"), ("gemmi", "gemmi")):
        try:
            __import__(module)
        except ImportError as exc:
            problems.append(f"{module} ({package}) will not import: {exc.__class__.__name__}: {exc}")
    return problems



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
    parser.add_argument(
        "--skip-preflight", action="store_true",
        help="Run the suites even if the environment does not satisfy the declared floors.",
    )
    args = parser.parse_args(argv)

    problems = preflight()
    if problems:
        for problem in problems:
            print(f"[ENV] {problem}")
        if not args.skip_preflight:
            print("\nThe environment does not satisfy what the packages declare. Test failures "
                  "here would describe this environment, not the software. Fix the environment, "
                  "or pass --skip-preflight to run anyway.")
            if args.json_out:
                args.json_out.parent.mkdir(parents=True, exist_ok=True)
                args.json_out.write_text(json.dumps(
                    {"passed": False, "environment_problems": problems}, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
            return 1
        print("[ENV] --skip-preflight given; results below describe this environment.\n")

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
