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
import http.client
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


def fetch(url: str, dest: Path, attempts: int = 4) -> None:
    """Download one artifact, retrying transient transport failures.

    Providers truncate connections under load. Without a retry the whole
    acquisition fails on a hiccup that has nothing to do with the science, and
    a scheduled qualification run turns into noise. The digest check downstream
    is what actually establishes correctness, so retrying is safe.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "yauvi-qualification/2.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                declared = r.headers.get("Content-Length")
                data = r.read()
            if declared is not None and len(data) != int(declared):
                raise http.client.IncompleteRead(data, int(declared) - len(data))
            if url.endswith(".gz") and not dest.name.endswith(".gz"):
                import gzip
                data = gzip.decompress(data)
            dest.write_bytes(data)
            return
        except (urllib.error.URLError, http.client.HTTPException, OSError, EOFError) as exc:
            # http.client.IncompleteRead derives from HTTPException, not from
            # URLError or OSError, so it escapes the obvious except clause and a
            # truncated download reads as a hard failure. Providers truncate
            # under load often enough that this is the common retry case.
            last = exc
            if attempt < attempts:
                time.sleep(2 ** attempt)
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
                fetch(url, artifact)
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
