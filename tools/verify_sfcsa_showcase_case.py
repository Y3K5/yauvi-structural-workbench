#!/usr/bin/env python3
"""Verify the SF-CSA stubbed-pipeline showcase case and its claim boundary."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
SHOWCASE = ROOT / "yauvi-structural-workbench" / "showcase" / "sfcsa-ceiling-case"
ADOPTION = ROOT / "yauvi-structural-workbench" / "docs" / "SFCSA_ARTIFACT_ADOPTION.json"
ARTIFACT = ROOT / "artifacts" / "protein-platform-modularization-reviewed.zip"


def fail(message: str) -> None:
    raise SystemExit(f"SF-CSA showcase verification failed: {message}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    case = json.loads((SHOWCASE / "CASE.json").read_text(encoding="utf-8"))
    adoption = json.loads(ADOPTION.read_text(encoding="utf-8"))
    # The reviewed archive is not redistributed with this distribution, so in a
    # published checkout this check has no artifact to hash. Say that, rather
    # than raising FileNotFoundError from the middle of a checksum comparison.
    # It still fails: an absent archive is unverified provenance, and the one
    # thing this repository may never do is let missing evidence pass as
    # favorable evidence.
    if not ARTIFACT.is_file():
        fail(
            f"the reviewed source archive is not present at {ARTIFACT.relative_to(ROOT)}. "
            "It is not redistributed with this distribution, so this check can only run "
            "in the development tree that holds it. Its recorded digest is "
            f"{adoption.get('source_archive', {}).get('sha256', 'unrecorded')}."
        )
    if sha256(ARTIFACT) != adoption.get("source_archive", {}).get("sha256"):
        fail("reviewed artifact checksum no longer matches the adoption record")
    for relative, expected in adoption.get("canonical_file_sha256", {}).items():
        if not (ROOT / relative).is_file() or sha256(ROOT / relative) != expected:
            fail(f"canonical adopted component drifted: {relative}")
    if not adoption.get("not_adopted") or not adoption.get("open_scientific_gates"):
        fail("artifact exclusions or open scientific gates are missing")
    if case.get("case_id") != "HUC-06" or case.get("analysis_type") != "sf_csa":
        fail("case identity is incorrect")
    if case.get("test_state") != "passed_stubbed_pipeline_case":
        fail("case is not explicitly labeled as a stubbed pipeline execution")
    if case.get("exit_codes") != {"run": 0, "verify": 0}:
        fail("canonical run or release audit did not pass")
    disclosure = (case.get("runtime_disclosure") or "").lower()
    if "computed no alignments" not in disclosure or "test doubles" not in disclosure:
        fail("stub-runtime disclosure is incomplete")
    if case.get("scientific_qualification_state") != "external_benchmark_pending":
        fail("stubbed execution was promoted to external scientific qualification")
    if case.get("measurements", [])[-1].get("value") != "stubbed":
        fail("alignment engines are not visibly labeled stubbed")
    for relative, expected in {**case.get("input_sha256", {}), **case.get("evidence_sha256", {})}.items():
        path = (SHOWCASE / relative).resolve()
        try:
            path.relative_to(SHOWCASE.resolve())
        except ValueError:
            fail(f"unsafe case path: {relative}")
        if not path.is_file() or sha256(path) != expected:
            fail(f"evidence checksum mismatch: {relative}")
    checksums = json.loads((SHOWCASE / "CHECKSUMS.json").read_text(encoding="utf-8"))
    expected_files = {
        path.relative_to(SHOWCASE).as_posix() for path in SHOWCASE.rglob("*")
        if path.is_file() and path.name != "CHECKSUMS.json"
    }
    if set(checksums) != expected_files:
        fail("CHECKSUMS.json does not cover the complete case")
    for relative, expected in checksums.items():
        if sha256(SHOWCASE / relative) != expected:
            fail(f"package checksum mismatch: {relative}")
    forbidden = re.compile(r"(?:/Users/|/private/var/|file://|20\d\d-\d\d-\d\dT\d\d:)")
    for path in SHOWCASE.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".json", ".tsv", ".txt", ".md", ".html"}:
            if forbidden.search(path.read_text(encoding="utf-8", errors="replace")):
                fail(f"absolute path or timestamp leaked into {path.relative_to(SHOWCASE)}")
    print("verified SF-CSA showcase: canonical run and audit passed; engines disclosed as test doubles")
    return 0


if __name__ == "__main__":
    sys.exit(main())
