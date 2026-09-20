#!/usr/bin/env python3
"""Acquire the public files named in public_research_sources.json."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import urllib.request

ROOT = Path(__file__).resolve().parent
LOCK = ROOT / "public_research_sources.json"


def find_source_root(script_root: Path) -> Path:
    script_root = script_root.resolve()
    for candidate in (script_root, *script_root.parents):
        if (candidate / "pyproject.toml").is_file() and (candidate / "software").is_dir():
            return candidate.resolve()
    return script_root


SOURCE_TREE_ROOT = find_source_root(ROOT)


def assert_outside_source(path: Path, source_root: Path = SOURCE_TREE_ROOT) -> None:
    path = path.expanduser().resolve()
    source_root = source_root.resolve()
    if path == source_root or source_root in path.parents:
        raise ValueError(f"input output directory must be outside the source tree: {path}")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="new input directory outside this source tree")
    args = parser.parse_args()
    out = args.out.expanduser().resolve()
    try:
        assert_outside_source(out)
    except ValueError as exc:
        parser.error(str(exc))
    if out.exists() and any(out.iterdir()):
        parser.error(f"output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    sources = [lock["reference"], *lock["structures"]]
    downloaded = []
    for source in sources:
        path = out / source["path"]
        request = urllib.request.Request(source.get("coordinate_url", source.get("url")), headers={"User-Agent": "YAUVI-public-research-example/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            raw = response.read()
            final_url = response.geturl()
        if not raw:
            raise RuntimeError(f"empty response for {source['path']}")
        path.write_bytes(raw)
        actual = digest(path)
        expected = source.get("sha256") or ""
        if expected and actual != expected:
            path.unlink()
            raise RuntimeError(f"SHA-256 mismatch for {source['path']}: expected {expected}, received {actual}")
        expected_size = source.get("byte_count")
        if expected_size is not None and len(raw) != expected_size:
            path.unlink()
            raise RuntimeError(f"byte-count mismatch for {source['path']}: expected {expected_size}, received {len(raw)}")
        downloaded.append({"path": source["path"], "sha256": actual, "byte_count": len(raw), "retrieved_url": final_url})
    print(json.dumps({"input_directory": str(out), "files": downloaded}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"acquire_public_inputs: {exc}", file=sys.stderr)
        raise SystemExit(2)
