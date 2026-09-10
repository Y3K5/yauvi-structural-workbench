#!/usr/bin/env python3
"""Add a verified mirror URL to matching SOURCE_LOCK.json entries.

A mirror is only worth having if it serves the exact bytes the lock already
records, so this refuses to write one it has not proven. For each matching
entry it downloads the candidate mirror, applies the same decompression rule the
acquirer applies, digests the result, and requires it to equal the recorded
sha256. One mismatch and nothing is written at all -- a lock half-pointed at a
mirror that disagrees is worse than a lock with no mirror.

`url` is left untouched. It stays the archive of record: the citable location,
the one the manifest names, and the one that has to resolve years from now. The
mirror only changes which host carries the traffic.

Typical use, after publishing the ten proteomes as a Hugging Face dataset and
pinning the commit that contains them:

    python tools/add_source_mirror.py \\
        --match sources/proteomes/ \\
        --template 'https://huggingface.co/datasets/OWNER/NAME/resolve/COMMIT_SHA/{filename}'

`{filename}` expands to the basename of the artifact plus the archive URL's
compression suffix, so `sources/proteomes/UP000000579.fasta` locked against a
`.fasta.gz` archive URL resolves to `UP000000579.fasta.gz`. Use `{artifact}`
for the full locked path, or `{basename}` for the name with no suffix added.

Pin a commit SHA rather than a branch. A branch name is a moving target, which
is the whole reason these files could not be locked against a live UniProt query
in the first place; a commit SHA is immutable and is a stronger version
parameter than a record id.

Usage:
    python tools/add_source_mirror.py --match PREFIX --template URL [--dry-run]
    python tools/add_source_mirror.py --match PREFIX --remove-mirrors
Exit: 0 written (or nothing to do), 1 a candidate failed verification, 2 bad usage.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2" / "SOURCE_LOCK.json"
USER_AGENT = "yauvi-qualification/2.0"


def expand(template: str, entry: dict) -> str:
    artifact = entry["artifact"]
    basename = Path(artifact).name
    url = entry.get("url") or ""
    suffix = ".gz" if url.endswith(".gz") and not basename.endswith(".gz") else ""
    return template.format(artifact=artifact, basename=basename, filename=basename + suffix)


def digest_as_acquirer_would(data: bytes, url: str, artifact: str) -> str:
    """Hash what the acquirer would store, not what the wire carried.

    `acquire_sources.py` decompresses a `.gz` URL whose artifact path is not
    `.gz` before hashing. gzip embeds an mtime, so a mirror re-compressed at a
    different moment carries different bytes while holding identical content --
    comparing the archives would reject a perfectly good mirror.
    """
    if url.endswith(".gz") and not artifact.endswith(".gz"):
        data = gzip.decompress(data)
    return hashlib.sha256(data).hexdigest()


def verify(url: str, entry: dict) -> str | None:
    """Return None if the mirror serves the locked bytes, else why it does not."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            data = response.read()
    except urllib.error.HTTPError as exc:
        return f"HTTP {exc.code} {exc.reason}"
    except Exception as exc:  # transport, DNS, TLS
        return f"{type(exc).__name__}: {exc}"

    try:
        observed = digest_as_acquirer_would(data, url, entry["artifact"])
    except (OSError, EOFError) as exc:
        return f"not the archive it claims to be ({len(data)} bytes): {exc}"

    expected = entry.get("sha256", "")
    if observed != expected:
        return f"digest mismatch: expected {expected[:12]}..., served {observed[:12]}..."
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--match", required=True,
                        help="Add the mirror to entries whose artifact path starts with this.")
    parser.add_argument("--template",
                        help="Mirror URL, with {filename}, {basename} or {artifact}.")
    parser.add_argument("--remove-mirrors", action="store_true",
                        help="Drop every mirror from matching entries instead of adding one.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Verify and report; write nothing.")
    args = parser.parse_args(argv)

    if bool(args.template) == bool(args.remove_mirrors):
        parser.error("give either --template or --remove-mirrors")

    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    matched = [e for e in lock.get("sources", []) if e.get("artifact", "").startswith(args.match)]
    if not matched:
        print(f"No locked artifact starts with {args.match!r}.", file=sys.stderr)
        return 2

    if args.remove_mirrors:
        removed = sum(1 for e in matched if e.pop("mirrors", None) is not None)
        if not args.dry_run and removed:
            LOCK.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
        print(f"{'would remove' if args.dry_run else 'removed'} mirrors from {removed} entries")
        return 0

    print(f"verifying {len(matched)} candidate mirrors against their locked digests")
    planned: list[tuple[dict, str]] = []
    failures: list[str] = []
    for entry in matched:
        url = expand(args.template, entry)
        problem = verify(url, entry)
        if problem:
            failures.append(f"{entry['artifact']}\n    {url}\n    {problem}")
            print(f"  FAILED  {entry['artifact']}")
        else:
            planned.append((entry, url))
            print(f"  ok      {entry['artifact']}")

    if failures:
        print(f"\n{len(failures)} of {len(matched)} candidates did not serve the locked bytes. "
              f"Nothing was written.", file=sys.stderr)
        print("\n".join("  " + f for f in failures), file=sys.stderr)
        return 1

    if args.dry_run:
        print(f"\nDry run. {len(planned)} verified mirrors would be recorded.")
        return 0

    for entry, url in planned:
        mirrors = entry.get("mirrors") or []
        if isinstance(mirrors, str):
            mirrors = [mirrors]
        if url not in mirrors:
            mirrors.append(url)
        entry["mirrors"] = mirrors
    LOCK.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")

    print(f"\nRecorded {len(planned)} verified mirrors in {LOCK.name}.")
    print("The archive of record in `url` is unchanged, and no recorded checksum moved: "
          "a mirror adds a retrieval path, not an adopted source, so the manifest's "
          "immutability policy is untouched.")
    print("Next, in order:")
    print(f"  python tools/verify_source_lock_health.py --only "
          f"{args.match.rstrip('/').split('/')[-1]}")
    print("  python tools/verify_release_status_digests.py --update   "
          "# the lock's own digest moved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
