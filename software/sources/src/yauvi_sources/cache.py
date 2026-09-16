"""The local store of acquired input files.

Layout, rooted at the cache directory (default `~/.cache/yauvi/sources`, override
with `YAUVI_SOURCE_CACHE` or `--cache`):

    <cache>/<source_id>/<sha256>/<filename>
    <cache>/<source_id>/MANIFEST.json

Content addressing by digest is what makes a run reproducible: two files that
differ in a single byte occupy different directories, so a silently-updated
remote release can never overwrite the copy an earlier run was built on. The
per-source `MANIFEST.json` is an append-only list of every entry ever acquired,
newest last, each recording where it came from and when.

`retrieved_at` is the one field here that is not derived from content. It is
recorded because a source whose upstream carries no version — the registry says
of UniProt, "the release is whatever UniProt serves that day" — can only be
pinned by the date it was taken plus the digest of what arrived.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Iterator, Sequence

DEFAULT_CACHE_ENV = "YAUVI_SOURCE_CACHE"
MANIFEST_NAME = "MANIFEST.json"


class CacheError(RuntimeError):
    pass


@dataclass(frozen=True)
class CacheEntry:
    """One acquired file, identified by its content."""

    source_id: str
    filename: str
    sha256: str
    bytes: int
    retrieved_at: str
    origin: str          # URL fetched, or 'staged:<original path>' for manual files
    version: str = ""    # upstream release string when the endpoint reports one
    note: str = ""

    @property
    def relative_path(self) -> str:
        return f"{self.source_id}/{self.sha256}/{self.filename}"


def default_cache_dir() -> Path:
    override = os.environ.get(DEFAULT_CACHE_ENV)
    if override:
        return Path(override).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    root = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return root / "yauvi" / "sources"


def sha256_file(path: str | Path, *, chunk: int = 1 << 20) -> tuple[str, int]:
    """Digest a file without reading it all into memory. Returns (hexdigest, bytes)."""
    digest = hashlib.sha256()
    total = 0
    with Path(path).open("rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
            total += len(block)
    return digest.hexdigest(), total


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class SourceCache:
    """Content-addressed storage for acquired source files."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root).expanduser() if root else default_cache_dir()

    # -- paths ------------------------------------------------------------

    def source_dir(self, source_id: str) -> Path:
        return self.root / source_id

    def manifest_path(self, source_id: str) -> Path:
        return self.source_dir(source_id) / MANIFEST_NAME

    def path_for(self, entry: CacheEntry) -> Path:
        return self.root / entry.relative_path

    # -- manifest ---------------------------------------------------------

    def entries(self, source_id: str) -> list[CacheEntry]:
        path = self.manifest_path(source_id)
        if not path.is_file():
            return []
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise CacheError(f"cache manifest for {source_id} is unreadable: {exc}") from exc
        known = set(CacheEntry.__dataclass_fields__)
        return [CacheEntry(**{k: v for k, v in item.items() if k in known}) for item in raw]

    def latest(self, source_id: str) -> CacheEntry | None:
        """The most recently acquired entry for a source, if any is present on disk."""
        for entry in reversed(self.entries(source_id)):
            if self.path_for(entry).is_file():
                return entry
        return None

    def has(self, source_id: str) -> bool:
        return self.latest(source_id) is not None

    def _write_manifest(self, source_id: str, entries: Sequence[CacheEntry]) -> None:
        path = self.manifest_path(source_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps([asdict(e) for e in entries], indent=2, sort_keys=True) + "\n"
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(payload, encoding="utf-8")
        tmp.replace(path)

    # -- writing ----------------------------------------------------------

    def store(
        self,
        source_id: str,
        payload: bytes,
        *,
        filename: str,
        origin: str,
        version: str = "",
        note: str = "",
    ) -> CacheEntry:
        """Write bytes into the cache under their own digest and record them."""
        if not payload:
            raise CacheError(f"refusing to cache an empty payload for {source_id}")
        digest = hashlib.sha256(payload).hexdigest()
        target = self.root / source_id / digest / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            tmp = target.with_name(target.name + ".part")
            tmp.write_bytes(payload)
            tmp.replace(target)

        entry = CacheEntry(
            source_id=source_id,
            filename=filename,
            sha256=digest,
            bytes=len(payload),
            retrieved_at=utc_now(),
            origin=origin,
            version=version,
            note=note,
        )
        existing = self.entries(source_id)
        # Re-acquiring identical content is a no-op in the manifest: the digest
        # already pins it, and a duplicate row would imply a change that did not
        # happen.
        if not any(e.sha256 == digest and e.filename == filename for e in existing):
            existing.append(entry)
            self._write_manifest(source_id, existing)
            return entry
        return next(e for e in existing if e.sha256 == digest and e.filename == filename)

    def stage(self, source_id: str, path: str | Path, *, note: str = "") -> CacheEntry:
        """Adopt a file the user acquired by hand, hashing it on the way in."""
        source_path = Path(path).expanduser()
        if not source_path.is_file():
            raise CacheError(f"cannot stage {source_id}: no such file: {source_path}")
        payload = source_path.read_bytes()
        return self.store(
            source_id,
            payload,
            filename=source_path.name,
            origin=f"staged:{source_path}",
            note=note or "manually acquired",
        )

    # -- verification -----------------------------------------------------

    def verify(self, source_id: str | None = None) -> Iterator[tuple[CacheEntry, bool, str]]:
        """Re-hash cached files against their manifests. Yields (entry, ok, detail)."""
        source_ids = [source_id] if source_id else self.source_ids()
        for sid in source_ids:
            for entry in self.entries(sid):
                path = self.path_for(entry)
                if not path.is_file():
                    yield entry, False, "missing from disk"
                    continue
                digest, size = sha256_file(path)
                if digest != entry.sha256:
                    yield entry, False, f"digest mismatch: on disk {digest}"
                elif size != entry.bytes:
                    yield entry, False, f"size mismatch: on disk {size} bytes"
                else:
                    yield entry, True, "ok"

    def source_ids(self) -> list[str]:
        if not self.root.is_dir():
            return []
        return sorted(p.name for p in self.root.iterdir() if (p / MANIFEST_NAME).is_file())
