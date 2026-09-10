#!/usr/bin/env python3
"""Add a verified mirror URL to matching SOURCE_LOCK.json entries.

A mirror is only worth having if it serves the exact bytes the lock already
records, so this refuses to write one it has not proven. For each matching
entry it downloads the candidate mirror *through the acquirer's own `fetch`* --
same retries, same decompression rule -- digests the result, and requires it to
equal the recorded sha256. Sharing that download is the point: a verifier with a
weaker fetch than the acquirer refuses mirrors that would have worked. One mismatch and nothing is written at all -- a lock half-pointed at a
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
import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COLLECTION = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
LOCK = COLLECTION / "SOURCE_LOCK.json"
ACQUIRER = COLLECTION / "acquire_sources.py"


def expand(template: str, entry: dict) -> str:
    artifact = entry["artifact"]
    basename = Path(artifact).name
    url = entry.get("url") or ""
    suffix = ".gz" if url.endswith(".gz") and not basename.endswith(".gz") else ""
    return template.format(artifact=artifact, basename=basename, filename=basename + suffix)


_acquire_sources = None


def acquirer():
    """The acquirer's own downloader, loaded from the collection it belongs to.

    A mirror is worth recording only if acquisition can actually use it, so it
    is proven by downloading it exactly the way acquisition will: the same retry
    policy, the same decompression rule, the same treatment of a definitive
    refusal.

    Reimplementing that here is how the two came to disagree. This verifier used
    to read the whole body in a single call with no retry, while
    `acquire_sources.fetch` streams, checks the declared length, and retries --
    because a provider under load truncates large responses. Verifying the ten
    proteome mirrors, three passes over the identical commit failed six, then
    one, then two of them, a different subset each time and every one an
    `IncompleteRead` a few kilobytes from the end. Not one was a digest
    mismatch. The mirror was good; the verifier was less robust than the thing
    it was verifying for, and so refused mirrors acquisition would have used.
    """
    global _acquire_sources
    if _acquire_sources is None:
        spec = importlib.util.spec_from_file_location("acquire_sources", ACQUIRER)
        if spec is None or spec.loader is None:  # pragma: no cover - packaging error
            raise RuntimeError(f"cannot load the acquirer from {ACQUIRER}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _acquire_sources = module
    return _acquire_sources


def verify(url: str, entry: dict) -> str | None:
    """Return None if the mirror serves the locked bytes, else why it does not.

    The digest is compared here rather than passed to `fetch` as its `expected`
    argument, deliberately. `fetch` treats a mismatch as retryable, which is
    right for acquisition -- a short read can produce one -- but wrong here: a
    mirror pointed at the wrong file would spend the whole retry budget on every
    artifact before saying so, and this tool exists to say so immediately. The
    retries stay where they belong, on the transport.
    """
    acquire = acquirer()
    with tempfile.TemporaryDirectory() as tmp:
        # The name decides decompression: `fetch` gunzips a `.gz` URL whose
        # destination is not `.gz`, which is what makes the digest comparable to
        # a lock that records decompressed content.
        dest = Path(tmp) / Path(entry["artifact"]).name
        try:
            acquire.fetch(url, dest)
        except acquire.Unacquirable as exc:
            return str(exc)  # a definitive answer: retrying cannot change it
        except Exception as exc:  # transport, truncation, TLS, DNS
            return f"{type(exc).__name__}: {exc}"
        observed = hashlib.sha256(dest.read_bytes()).hexdigest()

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
