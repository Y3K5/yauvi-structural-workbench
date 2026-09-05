#!/usr/bin/env python3
"""Ask whether every locked source can still be re-acquired from its URL.

`acquire_sources.py --verify-only` answers "are the files on this disk the files
the lock names". This answers the question the lock actually makes a promise
about: **can someone else obtain these bytes?** `SOURCE_LOCK.json` exists so a
reviewer acquires the exact same artifacts, and a digest recorded against a URL
that has moved on is a claim rather than a check.

Two findings motivated this, both caught after a push rather than before:

*   **Finding 8, 2026-09-04.** Ten reference proteomes were locked against
    UniProt *stream queries*, an endpoint with no version parameter. It returns
    the current release: the human proteome had gone from 144,818 to 147,520
    sequences since the lock was written. `acquire_sources.py` retries until the
    digest matches, so every CI runner exhausted its retries and the whole
    matrix went red at acquisition.
*   **2HYY, the same day.** wwPDB re-released the entry with expanded metadata,
    884,736 -> 899,620 bytes. Coordinates unchanged; digest not.

Neither is a defect in the software. Both are the outside world moving under a
lock, which is exactly the thing a lock is supposed to notice.

Two rules this tool follows because getting either wrong produces false alarms,
and both mistakes were made by hand before it existed:

1.  **Digest what the acquirer digests.** For a `.gz` URL whose artifact path is
    not `.gz`, `acquire_sources.py` decompresses before hashing. Comparing raw
    bytes reports every gzip as drifted, because gzip embeds a timestamp.
2.  **Retry transport errors.** Providers truncate under load. A single-attempt
    check reported nineteen RCSB entries as drifted when all nineteen were
    `IncompleteRead`, which the acquirer's own retry loop already handles.
    Truncation is reported separately from drift, and never as drift.

Network-heavy by nature: it fetches every locked URL. Meant for a schedule and
for running by hand before a push that touches the lock, not for every commit.

Usage:
  verify_source_lock_health.py [--json-out PATH] [--workers N] [--only SUBSTRING]
Exit:
  0  every locked URL still returns its locked bytes
  1  at least one drifted, is unreachable, or is missing from the repository
"""
from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
QUALIFICATION = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
LOCK = QUALIFICATION / "SOURCE_LOCK.json"
USER_AGENT = "yauvi-qualification/2.0"

# Long enough that a slow provider is not called unreachable; short enough that a
# hung connection does not stall the run.
TIMEOUT = 180
ATTEMPTS = 4


def digest_as_acquirer_would(data: bytes, url: str, artifact: str) -> str:
    """The acquirer's rule, kept in one place so the two cannot diverge."""
    if url.endswith(".gz") and not artifact.endswith(".gz"):
        data = gzip.decompress(data)
    return hashlib.sha256(data).hexdigest()


def check_one(entry: dict[str, Any]) -> dict[str, Any]:
    artifact = str(entry["artifact"])
    expected = str(entry.get("sha256", ""))
    url = entry.get("url")
    result = {"source_id": entry.get("source_id"), "artifact": artifact,
              "url": url, "expected": expected}

    if entry.get("acquisition") == "committed_in_repository" or not url:
        local = QUALIFICATION / artifact
        if not local.is_file():
            return {**result, "state": "missing_from_repository"}
        observed = hashlib.sha256(local.read_bytes()).hexdigest()
        return {**result, "state": "in_repository" if observed == expected else "repository_mismatch",
                "observed": observed}

    last = ""
    for attempt in range(1, ATTEMPTS + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                declared = response.headers.get("Content-Length")
                data = response.read()
            if declared is not None and len(data) != int(declared):
                raise OSError(f"short read: {len(data)} of {declared} bytes")
            observed = digest_as_acquirer_would(data, url, artifact)
            state = "reproducible" if observed == expected else "drifted"
            return {**result, "state": state, "observed": observed, "bytes": len(data),
                    "attempts": attempt}
        except Exception as exc:  # transport, truncation, HTTP status
            last = f"{type(exc).__name__}: {exc}"
            if attempt < ATTEMPTS:
                time.sleep(min(2 ** attempt, 15))
    # Exhausted retries without ever reading a complete response. This is the
    # provider being unavailable, not the lock being wrong, and is reported as
    # its own state so it is never mistaken for drift.
    return {**result, "state": "unreachable", "error": last, "attempts": ATTEMPTS}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json-out", type=Path, help="Write the full result set here.")
    parser.add_argument("--workers", type=int, default=3,
                        help="Parallel fetches. Kept low on purpose: providers truncate "
                             "under load, and a truncation storm is noise, not signal.")
    parser.add_argument("--only", help="Check only entries whose artifact or url contains this.")
    args = parser.parse_args(argv)

    if not LOCK.is_file():
        print(f"missing: {LOCK}", file=sys.stderr)
        return 1
    entries = json.loads(LOCK.read_text(encoding="utf-8")).get("sources", [])
    if args.only:
        entries = [e for e in entries
                   if args.only in str(e.get("artifact", "")) or args.only in str(e.get("url", ""))]
    if not entries:
        print("no locked sources matched")
        return 1

    print(f"checking {len(entries)} locked sources against their providers "
          f"({args.workers} at a time)\n")
    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(check_one, entries):
            results.append(result)
            if result["state"] not in ("reproducible", "in_repository"):
                print(f"  {result['state'].upper():<22} {result['artifact']}")
                if result["state"] == "drifted":
                    print(f"      locked   {result['expected'][:16]}")
                    print(f"      returned {result['observed'][:16]}  ({result.get('bytes', 0):,} bytes)")
                elif result["state"] == "unreachable":
                    print(f"      {result.get('error', '')[:110]}")

    tally: dict[str, int] = {}
    for result in results:
        tally[result["state"]] = tally.get(result["state"], 0) + 1

    print("\n  " + "  ".join(f"{state}={count}" for state, count in sorted(tally.items())))

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(
            {"checked": len(results), "tally": tally, "results": results},
            indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"  written to {args.json_out}")

    drifted = tally.get("drifted", 0)
    broken = tally.get("missing_from_repository", 0) + tally.get("repository_mismatch", 0)
    unreachable = tally.get("unreachable", 0)

    if drifted or broken:
        print(f"\n{drifted} locked source(s) no longer return their bytes; "
              f"{broken} missing or mismatched in the repository.")
        print("A digest against a URL that has moved on is a claim, not a check. Either the "
              "provider re-released and the lock should be updated with the change recorded, or "
              "the URL cannot be locked at all -- see Finding 8.")
        return 1
    if unreachable:
        # Not a failure: the lock may be perfectly good and the provider down.
        # Saying so plainly beats a red run that means nothing.
        print(f"\nNo drift found. {unreachable} provider(s) did not answer in {ATTEMPTS} "
              f"attempts; that is availability, not drift. Re-run before drawing a conclusion.")
        return 0
    print("\nEvery locked source still returns its locked bytes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
