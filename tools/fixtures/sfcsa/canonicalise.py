#!/usr/bin/env python3
"""Rewrite an SF-CSA release into a machine-independent form for golden diffing.

A release as produced by `sf_csa.cli run` is *not* byte-comparable across
machines, and pretending otherwise is how golden fixtures rot into
"regenerate it until the diff goes away". What varies:

1.  `SF_CSA_RELEASE_MANIFEST.json` records `release_id`, which the pipeline
    takes from the basename of `--output`. A run into `/tmp/xyz` and a run
    into `./out` produce different manifests from identical inputs.
2.  `CHECKSUMS.json` digests that file, so it inherits the same variation.

**Closed 2026-09-04.** `proteome_denominator.json` used to be a third source of
variation: `build_proteome_universe` stored `str(path)` after `Path.resolve()`,
so it embedded the checkout location and had to be rewritten here. It now
records paths relative to the database root, so no substitution is needed and
the file is portable as produced. That was a privacy fix rather than a fixture
one -- the absolute form reached release evidence staged for a public
repository -- and the property is held by
`sf-csa/tests/test_release_paths_are_portable.py`. This note stays because a
canonicaliser that silently keeps substituting a field the pipeline no longer
emits is how the next reader concludes the pipeline still emits it.

This script copies a release into a canonical tree with those substitutions
applied and the checksums recomputed over the canonical bytes. Recomputing is
deliberate: the digests stay load-bearing in the golden tree, so a genuine
content change in any release file still shows up twice -- once in the file and
once in `CHECKSUMS.json` -- rather than being papered over by a placeholder.

`work/` is excluded. It holds pipeline scratch (`rbh_source_<hash>.dmnd`, where
the hash is `sha256` of the *absolute* source-proteome path) and is not part of
the release contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

FIXTURE_ROOT_TOKEN = "<FIXTURE_ROOT>"
RELEASE_ID_TOKEN = "<RELEASE_ID>"

# Relative to the release root. Everything else is copied byte-for-byte.
EXCLUDED_DIRS = ("work",)
CHECKSUM_FILE = "CHECKSUMS.json"


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _rewrite(text: str, fixture_root: Path, release_id: str) -> str:
    """Apply the two machine-dependent substitutions.

    The fixture root is replaced in both its resolved and its literal form,
    because a release may have been produced through a symlinked path (on macOS
    `/tmp` resolves to `/private/tmp`) and only one of the two will appear.

    The two forms are applied **longest first**, and that is load-bearing rather
    than tidy. When one form is a suffix of the other -- `/tmp/x` and
    `/private/tmp/x` -- replacing the short one first rewrites the tail of the
    long one and leaves the orphaned prefix behind, so the output reads
    `/private<FIXTURE_ROOT>` and the golden comparison fails on a machine
    detail. Iterating a set here made that failure depend on hash order: it did
    not reproduce in the checkout the fixture was authored in and appeared on
    the first run from a symlinked directory.
    """
    forms = sorted({str(fixture_root), str(fixture_root.resolve())}, key=len, reverse=True)
    for form in forms:
        text = text.replace(form, FIXTURE_ROOT_TOKEN)
    if release_id:
        text = text.replace(f'"release_id": "{release_id}"', f'"release_id": "{RELEASE_ID_TOKEN}"')
    return text


def canonicalise(release: Path, dest: Path, fixture_root: Path) -> list[str]:
    """Write a canonical copy of `release` into `dest`; return the file list."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)

    release_id = release.resolve().name
    written: list[str] = []

    for src in sorted(release.rglob("*")):
        if not src.is_file():
            continue
        rel = src.relative_to(release)
        if rel.parts and rel.parts[0] in EXCLUDED_DIRS:
            continue
        if rel.as_posix() == CHECKSUM_FILE:
            continue  # rebuilt at the end, over canonical bytes
        out = dest / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        raw = src.read_bytes()
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            out.write_bytes(raw)  # binary payload: copied untouched
        else:
            out.write_text(_rewrite(text, fixture_root, release_id), encoding="utf-8")
        written.append(rel.as_posix())

    # Rebuild CHECKSUMS.json over the canonical tree, keeping the original's
    # key set so a file that stopped being checksummed is still a visible diff.
    original = release / CHECKSUM_FILE
    if original.exists():
        keys = sorted(json.loads(original.read_text(encoding="utf-8")))
        rebuilt = {}
        for key in keys:
            target = dest / key
            rebuilt[key] = _sha256_bytes(target.read_bytes()) if target.exists() else "MISSING"
        (dest / CHECKSUM_FILE).write_text(
            json.dumps(rebuilt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        written.append(CHECKSUM_FILE)

    return sorted(written)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--release", required=True, help="release directory produced by sf_csa run")
    parser.add_argument("--dest", required=True, help="directory to write the canonical copy into")
    parser.add_argument(
        "--fixture-root",
        default=str(Path(__file__).resolve().parent),
        help="path whose occurrences are replaced with %s (default: this script's directory)"
        % FIXTURE_ROOT_TOKEN,
    )
    args = parser.parse_args(argv)

    files = canonicalise(Path(args.release), Path(args.dest), Path(args.fixture_root))
    print(f"canonicalised {len(files)} file(s) into {args.dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
