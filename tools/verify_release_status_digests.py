#!/usr/bin/env python3
"""Compare every digest recorded in RELEASE_STATUS.json against its file.

`RELEASE_STATUS.json` records a sha256 beside each Qualification v2 evidence
document so a reader can tell that the numbers quoted on the reviewer surfaces
came from the files in the tree. Nothing compared them, and the adoption
protocol names the consequence directly (rule 2): a stale membrane digest sat
in the public repository unnoticed, because a digest nothing checks is a claim,
not a check.

On 2026-09-02 three of eight were stale. The manifest had been revised twice
(collection 2.4 and the 2026-09-01 ABL requirement change, 114 required cases
to 110) and neither the composition audit nor the execution summary had been
regenerated, so the recorded totals described a panel that no longer existed.

Verify mode is the default and is what the test suite runs. `--update`
recomputes and rewrites the recorded digests, and is the correct way to record
evidence you have deliberately regenerated -- never a way to silence a
mismatch you have not explained.

Usage:  python tools/verify_release_status_digests.py [--update]
Exit:   0 every recorded digest matches, 1 one does not or its file is missing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / "yauvi-structural-workbench" / "RELEASE_STATUS.json"
BASE = STATUS.parent


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def recorded_pairs(status: dict) -> list[tuple[str, str, str]]:
    """Every (label, relative path, recorded digest) the status file declares.

    Two shapes: a `<key>` / `<key>_sha256` pair, and the `execution_results` /
    `execution_results_sha256` mapping keyed by workflow. Both are discovered
    rather than listed, so a new evidence document is covered the moment it is
    recorded in the same shape.
    """
    v2 = status["qualification_evidence"]["current_v2"]
    pairs: list[tuple[str, str, str]] = []
    for key, digest in v2.items():
        if not key.endswith("_sha256"):
            continue
        source = key[: -len("_sha256")]
        value = v2.get(source)
        if isinstance(value, str):
            pairs.append((source, value, digest))
        elif isinstance(value, dict) and isinstance(digest, dict):
            for workflow, path in value.items():
                if workflow in digest:
                    pairs.append((f"{source}:{workflow}", path, digest[workflow]))
    return pairs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update", action="store_true",
        help="Rewrite the recorded digests from the files. For evidence you meant to regenerate.",
    )
    args = parser.parse_args(argv)

    status = json.loads(STATUS.read_text(encoding="utf-8"))
    v2 = status["qualification_evidence"]["current_v2"]
    stale: list[str] = []
    missing: list[str] = []

    for label, relative, digest in recorded_pairs(status):
        actual = sha256(BASE / relative)
        if actual is None:
            missing.append(f"{label}: {relative} does not exist")
            continue
        if actual == digest:
            continue
        stale.append(f"{label}: {relative}\n    recorded {digest}\n    actual   {actual}")
        if args.update:
            key, _, workflow = label.partition(":")
            if workflow:
                v2[f"{key}_sha256"][workflow] = actual
            else:
                v2[f"{key}_sha256"] = actual

    if missing:
        for line in missing:
            print(f"MISSING {line}")
    if stale:
        for line in stale:
            print(f"{'UPDATED' if args.update else 'STALE'} {line}")

    if args.update and stale and not missing:
        STATUS.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        print(f"\nrewrote {len(stale)} digest(s) in {STATUS.relative_to(ROOT)}")
        return 0
    if not stale and not missing:
        print(f"{len(recorded_pairs(status))} recorded digests match their files")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
