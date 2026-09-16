#!/usr/bin/env python3
"""Acquire the artifacts named in SOURCE_LOCK.json and verify their digests.

Acquisition is the only step permitted to touch the network. Execution reads
the files this writes and never contacts a provider, which is what
`execution_policy.network_access: forbidden` requires.

Artifacts are never committed: `ships_public_records` is false, so the lock
records provider, URL, and SHA-256 and every consumer re-acquires the identical
bytes.

An entry may also carry `mirrors`: further URLs serving the identical bytes,
tried before `url`. `url` stays the archive of record -- the citable location a
reader is sent to -- while a mirror carries the traffic. Every candidate is held
to the same recorded digest, so a mirror can never weaken the lock; it can only
make acquisition survive one provider having a bad day.

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


#: HTTP statuses that will answer the same way however many times they are asked.
#: Retrying them wastes the backoff budget and, worse, buries the status: a run
#: that failed on a 404 and a run that failed on a truncated read both ended up
#: reporting "6 attempts failed", which is what made a broken acquisition step
#: unreadable in CI.
TERMINAL_STATUS = frozenset({400, 401, 403, 404, 405, 410, 451})

#: Statuses that mean "ask again later", and mean it on the provider's schedule
#: rather than ours. A repository serving a 150 MB set to a serial client will
#: throttle it, and the exponential backoff here tops out at twenty seconds --
#: shorter than the window a throttle usually asks for. Honouring `Retry-After`
#: is the difference between waiting the time the provider named and spending the
#: whole retry budget discovering that twenty seconds was not enough.
THROTTLED_STATUS = frozenset({429, 503})

#: A provider may name a very long window. Wait a bounded amount of it, so a
#: single artifact cannot hold the job to its 45-minute limit on its own.
MAX_RETRY_AFTER = 120


class Unacquirable(Exception):
    """The provider answered definitively, so retrying cannot help."""


def _retry_after(response_headers) -> float | None:
    """Seconds the provider asked us to wait, if it named a number."""
    value = response_headers.get("Retry-After") if response_headers else None
    if not value:
        return None
    try:
        return min(float(value.strip()), MAX_RETRY_AFTER)
    except ValueError:
        return None  # HTTP-date form; fall back to the normal backoff


def _describe(response_body: bytes, content_type: str | None) -> str:
    """A short, safe characterisation of what arrived instead of the artifact."""
    head = response_body[:80]
    if head.startswith(b"\x1f\x8b"):
        return "gzip"
    if head.lstrip()[:1] in (b"<",):
        return f"markup ({content_type or 'no content-type'})"
    if head.lstrip()[:1] in (b"{", b"["):
        return f"json ({content_type or 'no content-type'})"
    return f"{content_type or 'no content-type'}, first bytes {head[:24]!r}"


def fetch(url: str, dest: Path, expected: str = "", attempts: int = 6) -> None:
    """Download one artifact, retrying until its digest matches.

    Retrying only on transport errors is not enough: a provider under load can
    truncate a response in ways that surface as a short read, and the larger
    coordinate files here (several megabytes) hit that often enough to break a
    whole acquisition. The digest is the real success criterion, so the loop
    retries until the bytes on disk match the lock, streaming to a partial file
    so a failed attempt never leaves a plausible-looking artifact behind.

    A definitive answer is not retried. A 404 is not a busy provider; it is a URL
    that no longer names the bytes the lock recorded, and treating it as transient
    both hides the status and spends fifty seconds per artifact discovering the
    same thing six times. The same applies to a response that arrives intact but
    is not the artifact -- an HTML preview page where a gzip stream was expected
    means the URL addresses a landing page, not a file, and the message now says
    so instead of reporting a generic decompression error.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    last: Exception | None = None
    asked_to_wait: float | None = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "yauvi-qualification/2.0"})
            try:
                response = urllib.request.urlopen(req, timeout=300)
            except urllib.error.HTTPError as exc:
                if exc.code in TERMINAL_STATUS:
                    raise Unacquirable(f"HTTP {exc.code} {exc.reason} for {url}") from exc
                if exc.code in THROTTLED_STATUS:
                    asked_to_wait = _retry_after(exc.headers)
                raise
            with response as r:
                declared = r.headers.get("Content-Length")
                content_type = r.headers.get("Content-Type")
                with part.open("wb") as fh:
                    shutil.copyfileobj(r, fh, 1 << 16)
            size = part.stat().st_size
            if declared is not None and size != int(declared):
                raise http.client.IncompleteRead(b"", int(declared) - size)

            data = part.read_bytes()
            if url.endswith(".gz") and not dest.name.endswith(".gz"):
                try:
                    data = gzip.decompress(data)
                except (OSError, EOFError) as exc:
                    raise Unacquirable(
                        f"{url} ends .gz but served {_describe(data, content_type)} "
                        f"({size} bytes): {exc}"
                    ) from exc
                part.write_bytes(data)

            if expected:
                observed = hashlib.sha256(part.read_bytes()).hexdigest()
                if observed != expected:
                    raise ValueError(f"digest mismatch (got {observed[:12]}...)")
            part.replace(dest)
            return
        except Unacquirable:
            part.unlink(missing_ok=True)
            raise
        except Exception as exc:  # transport, truncation, or digest mismatch
            last = exc
            part.unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(asked_to_wait if asked_to_wait else min(2 ** attempt, 20))
                asked_to_wait = None
    raise RuntimeError(f"{attempts} attempts failed: {last}")


def candidates(entry: dict) -> list[str]:
    """Every URL that may serve this artifact, in the order to try them.

    `url` is the archive of record: the location the manifest cites, the one a
    reader is pointed at, and the one that has to still resolve in five years.
    `mirrors` are operational copies of the identical bytes.

    Mirrors are tried *first*, and the archive last. That looks backwards until
    you notice what broke: seven runners re-downloading a 146 MB set on every
    push is roughly a gigabyte a day aimed at a preservation archive, and the
    archive started returning gateway timeouts. The thing you cite and the thing
    you hammer should not be the same thing.

    Ordering cannot affect correctness, only speed and politeness: every
    candidate is held to the same recorded sha256, so a mirror that serves
    anything else fails exactly as a bad archive copy would. That is what makes
    a mirror free to add -- the lock already treats every host as untrusted.
    """
    mirrors = entry.get("mirrors") or []
    if isinstance(mirrors, str):
        mirrors = [mirrors]
    url = entry.get("url")
    ordered = [*mirrors, url] if url else list(mirrors)
    # Preserve order while dropping a mirror that merely repeats the archive URL.
    return list(dict.fromkeys(u for u in ordered if u))


def acquire(entry: dict, artifact: Path) -> tuple[str, list[str]]:
    """Fetch one artifact from the first candidate that serves its locked bytes.

    Returns the URL that worked and the failures met on the way, so a run that
    succeeded on a mirror still says the archive was unreachable rather than
    reporting a silent success. A degraded provider that nothing reports is how
    this became a three-day outage.
    """
    expected = entry.get("sha256", "")
    attempted: list[str] = []
    for url in candidates(entry):
        try:
            fetch(url, artifact, expected)
            return url, attempted
        except Exception as exc:
            attempted.append(f"{url}: {exc}")
    raise RuntimeError("; ".join(attempted) if attempted else "no url to acquire")


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
    degraded: list[str] = []
    restored = 0
    for entry in sources:
        artifact = HERE / entry["artifact"]
        expected = entry.get("sha256", "")
        if not artifact.is_file():
            if entry.get("acquisition") == "committed_in_repository":
                problems.append(f"missing from the repository: {entry['artifact']}")
                continue
            if args.verify_only:
                problems.append(f"missing: {entry['artifact']}")
                continue
            if not candidates(entry):
                problems.append(f"no url to acquire: {entry['artifact']}")
                continue
            try:
                served_by, attempted = acquire(entry, artifact)
            except Exception as exc:
                problems.append(f"download failed: {entry['artifact']}: {exc}")
                continue
            if attempted:
                degraded.append(f"{entry['artifact']}\n"
                                f"    served by {served_by}\n"
                                + "\n".join(f"    after {a}" for a in attempted))
        else:
            restored += 1
        observed = sha256(artifact)
        if observed != expected:
            problems.append(f"checksum mismatch: {entry['artifact']}\n"
                            f"    expected {expected}\n    observed {observed}")

    # A run that fell back is a working run and a warning at the same time. Say
    # so on success, or the provider that is failing stays invisible until the
    # day the fallback fails too.
    if degraded:
        print(f"{len(degraded)} artifact(s) came from a fallback candidate:")
        print("\n".join("  " + d for d in degraded))

    if problems:
        print(f"{len(problems)} of {len(sources)} locked artifacts failed:", file=sys.stderr)
        print("\n".join("  " + p for p in problems), file=sys.stderr)
        return 1
    print(f"all {len(sources)} locked artifacts present and checksum-verified "
          f"({restored} already present, {len(sources) - restored} acquired)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
