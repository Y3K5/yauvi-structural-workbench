#!/usr/bin/env python3
"""Acquire the artifacts named in SOURCE_LOCK.json and verify their digests.

Acquisition is the only step permitted to touch the network. Execution reads
the files this writes and never contacts a provider, which is what
`execution_policy.network_access: forbidden` requires.

Artifacts are never committed: `ships_public_records` is false, so the lock
records provider, URL, and SHA-256 and every consumer re-acquires the identical
bytes.

Usage:
  python acquire_sources.py                # download what is missing, verify all
  python acquire_sources.py --verify-only  # verify only, download nothing
Exit: 0 every locked artifact is present and matches, 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import gzip
import http.client
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOCK = HERE / "SOURCE_LOCK.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, expected: str = "", attempts: int = 6) -> None:
    """Download one artifact, retrying until its digest matches.

    Retrying only on transport errors is not enough: a provider under load can
    truncate a response in ways that surface as a short read, and the larger
    coordinate files here (several megabytes) hit that often enough to break a
    whole acquisition. The digest is the real success criterion, so the loop
    retries until the bytes on disk match the lock, streaming to a partial file
    so a failed attempt never leaves a plausible-looking artifact behind.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "yauvi-qualification/2.0"})
            with urllib.request.urlopen(req, timeout=300) as r:
                declared = r.headers.get("Content-Length")
                with part.open("wb") as fh:
                    shutil.copyfileobj(r, fh, 1 << 16)
            size = part.stat().st_size
            if declared is not None and size != int(declared):
                raise http.client.IncompleteRead(b"", int(declared) - size)

            data = part.read_bytes()
            if url.endswith(".gz") and not dest.name.endswith(".gz"):
                data = gzip.decompress(data)
                part.write_bytes(data)

            if expected:
                observed = hashlib.sha256(part.read_bytes()).hexdigest()
                if observed != expected:
                    raise ValueError(f"digest mismatch (got {observed[:12]}...)")
            part.replace(dest)
            return
        except Exception as exc:  # transport, truncation, or digest mismatch
            last = exc
            part.unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(min(2 ** attempt, 20))
    raise RuntimeError(f"{attempts} attempts failed: {last}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify-only", action="store_true",
                    help="Do not download; only check what is already present.")
    args = ap.parse_args(argv)

    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    sources = lock.get("sources", [])
    if not sources:
        print("SOURCE_LOCK.json adopts no sources; nothing to acquire.", file=sys.stderr)
        return 1

    problems: list[str] = []
    for entry in sources:
        artifact = HERE / entry["artifact"]
        expected = entry.get("sha256", "")
        url = entry.get("url")
        if not artifact.is_file():
            if entry.get("acquisition") == "committed_in_repository":
                problems.append(f"missing from the repository: {entry['artifact']}")
                continue
            if args.verify_only:
                problems.append(f"missing: {entry['artifact']}")
                continue
            if not url:
                problems.append(f"no url to acquire: {entry['artifact']}")
                continue
            try:
                fetch(url, artifact, expected)
            except Exception as exc:
                problems.append(f"download failed: {entry['artifact']}: {exc}")
                continue
        observed = sha256(artifact)
        if observed != expected:
            problems.append(f"checksum mismatch: {entry['artifact']}\n"
                            f"    expected {expected}\n    observed {observed}")

    if problems:
        print(f"{len(problems)} of {len(sources)} locked artifacts failed:", file=sys.stderr)
        print("\n".join("  " + p for p in problems), file=sys.stderr)
        return 1
    print(f"all {len(sources)} locked artifacts present and checksum-verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
