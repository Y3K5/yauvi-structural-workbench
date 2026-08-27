#!/usr/bin/env python3
"""Verify the five-case structural showcase and its human-readable boundaries."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
SHOWCASE = ROOT / "yauvi-structural-workbench" / "showcase" / "five-human-use-cases"


def fail(message: str) -> None:
    raise SystemExit(f"showcase verification failed: {message}")


def main() -> int:
    summary = json.loads((SHOWCASE / "SHOWCASE.json").read_text(encoding="utf-8"))
    cases = summary.get("cases", [])
    if summary.get("test_cases") != 5 or summary.get("passed") != 5 or len(cases) != 5:
        fail("expected five passed cases")
    if [case.get("case_id") for case in cases] != [f"HUC-{index:02d}" for index in range(1, 6)]:
        fail("case identifiers are incomplete or out of order")
    for case in cases:
        if case.get("test_state") != "passed" or case.get("exit_code") != 0:
            fail(f'{case.get("case_id")} did not complete with exit code 0')
        for field in ("human_label", "human_question", "observed_result", "human_benefits", "non_claim"):
            if not case.get(field):
                fail(f'{case["case_id"]} is missing {field}')
        for relative in case.get("evidence_files", []):
            path = (SHOWCASE / relative).resolve()
            try:
                path.relative_to(SHOWCASE.resolve())
            except ValueError:
                fail(f"unsafe evidence path: {relative}")
            if not path.is_file():
                fail(f"missing evidence file: {relative}")

    records = {case["case_id"]: case for case in cases}
    if records["HUC-01"]["measurements"][0]["value"] != "100.0%":
        fail("coordinate mapping coverage changed")
    if records["HUC-02"]["measurements"][0]["value"] != "tm_helix_experimental":
        fail("alpha-helical prototype label changed")
    if records["HUC-02"]["measurements"][1]["value"] != "tm_helix_axis_v2":
        fail("checksum-bound alpha-helical axis method was not used")
    if "Mark 1 alpha-helical qualification" not in records["HUC-02"]["non_claim"]:
        fail("alpha-helical qualification boundary is missing")
    if records["HUC-03"]["measurements"][0]["value"] != "active_like":
        fail("state resemblance fixture changed")
    if records["HUC-04"]["measurements"][1]["value"] != "3":
        fail("functional-site role mapping changed")
    if float(records["HUC-05"]["measurements"][2]["value"].split()[0]) <= 0:
        fail("assembly burial is no longer positive")

    checksums = json.loads((SHOWCASE / "CHECKSUMS.json").read_text(encoding="utf-8"))
    for relative, expected in checksums.items():
        path = SHOWCASE / relative
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            fail(f"checksum mismatch: {relative}")
    forbidden = re.compile(r"(?:/Users/|/private/var/|file://|20\d\d-\d\d-\d\dT\d\d:)")
    for path in SHOWCASE.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".json", ".tsv", ".txt", ".html", ".md"}:
            if forbidden.search(path.read_text(encoding="utf-8", errors="replace")):
                fail(f"absolute path or timestamp leaked into {path.relative_to(SHOWCASE)}")
    print("verified five human use cases: 5/5 passed; checksums and claim boundaries intact")
    return 0


if __name__ == "__main__":
    sys.exit(main())
