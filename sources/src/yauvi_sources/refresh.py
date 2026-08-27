"""Opt-in refresh of public oral-reference caches.

This layer never reads Protein Cases, atlas observations, local omics, or user
targets. It uses only a fixed bulk-source list, updates the content-addressed
source cache, and leaves adoption into an atlas revision to an explicit build.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
from typing import Any, Callable, Mapping

from .cache import SourceCache
from .fetchers.http import _public_get


REFRESH_SPECS: Mapping[str, Mapping[str, Any]] = {
    "homd": {
        "url": "https://v31a.homd.org/download", "filename": "ehomd-download-index.html",
        "max_bytes": 4 * 1024 * 1024, "signature": lambda body: b"homd" in body.lower() and b"download" in body.lower(),
    },
    "human_protein_atlas": {
        "url": "https://www.proteinatlas.org/download/proteinatlas.tsv.zip", "filename": "proteinatlas.tsv.zip",
        "max_bytes": 128 * 1024 * 1024, "signature": lambda body: body.startswith(b"PK"),
    },
    "rhea": {
        "url": "https://ftp.expasy.org/databases/rhea/tsv/rhea-tsv.tar.gz", "filename": "rhea-tsv.tar.gz",
        "max_bytes": 128 * 1024 * 1024, "signature": lambda body: body.startswith(b"\x1f\x8b"),
    },
    "proteomexchange": {
        "url": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD006367&outputMode=JSON",
        "filename": "PXD006367.json", "max_bytes": 8 * 1024 * 1024,
        "signature": lambda body: b"PXD006367" in body,
    },
}


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _stamp(value: datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


class ReferenceRefreshManager:
    """Bounded status and refresh coordinator for fixed public providers."""

    def __init__(self, cache: SourceCache | str | Path | None = None, *,
                 now: Callable[[], datetime] = _utc_now) -> None:
        self.cache = cache if isinstance(cache, SourceCache) else SourceCache(cache)
        self.now = now
        self.state_path = self.cache.root / "ORAL_REFERENCE_REFRESH.json"

    def _state(self) -> dict[str, Any]:
        if not self.state_path.is_file():
            return {"schema_version": "1.0", "sources": {}}
        try:
            value = json.loads(self.state_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {"schema_version": "1.0", "sources": {}}
        return value if isinstance(value, dict) else {"schema_version": "1.0", "sources": {}}

    def _write_state(self, value: Mapping[str, Any]) -> None:
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.state_path.with_name(f".{self.state_path.name}.{os.getpid()}.tmp")
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, self.state_path)

    def status(self, *, online_enabled: bool) -> dict[str, Any]:
        state = self._state()
        result = {}
        for source_id in REFRESH_SPECS:
            latest = self.cache.latest(source_id)
            saved = state.get("sources", {}).get(source_id, {})
            status = saved.get("status", "offline" if not online_enabled else "source_failed")
            if not online_enabled:
                status = "offline"
            result[source_id] = {
                "status": status, "checked_at": saved.get("checked_at", ""),
                "sha256": latest.sha256 if latest else "", "version": latest.version if latest else "",
                "detail": saved.get("detail", "Reference refresh is disabled." if not online_enabled else "No successful check recorded."),
            }
        return {"online_enabled": online_enabled, "sources": result,
                "limitations": ["Refresh updates only the source cache; atlas revisions adopt new checksums explicitly."]}

    def refresh(self, *, force: bool = False) -> dict[str, Any]:
        state = self._state()
        sources = state.setdefault("sources", {})
        now = self.now()
        for source_id, spec in REFRESH_SPECS.items():
            prior = sources.get(source_id, {})
            if not force and prior.get("checked_at"):
                try:
                    checked = datetime.strptime(prior["checked_at"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
                except ValueError:
                    checked = None
                if checked is not None and (now - checked).total_seconds() < 24 * 3600:
                    continue
            latest = self.cache.latest(source_id)
            headers = {}
            if prior.get("etag"):
                headers["If-None-Match"] = prior["etag"]
            if prior.get("last_modified"):
                headers["If-Modified-Since"] = prior["last_modified"]
            response, reason = _public_get(spec["url"], headers=headers, max_bytes=int(spec["max_bytes"]))
            checked_at = _stamp(now)
            if reason == "not_modified":
                sources[source_id] = {**prior, "status": "current", "checked_at": checked_at, "detail": "Provider reported no change."}
                continue
            if response is None:
                sources[source_id] = {**prior, "status": "stale_cache" if latest else ("offline" if reason.startswith("network_") else "source_failed"),
                                      "checked_at": checked_at, "detail": reason}
                continue
            body = response.content
            if not spec["signature"](body):
                sources[source_id] = {**prior, "status": "stale_cache" if latest else "source_failed",
                                      "checked_at": checked_at, "detail": "content_signature_mismatch"}
                continue
            entry = self.cache.store(source_id, body, filename=spec["filename"], origin=response.url,
                                     version=str(response.headers.get("Last-Modified", "")), note="opt-in oral reference refresh")
            changed = latest is not None and latest.sha256 != entry.sha256
            sources[source_id] = {
                "status": "update_available" if changed else "current", "checked_at": checked_at,
                "detail": "New content is cached and requires explicit atlas revision adoption." if changed else "Cached public reference is current.",
                "sha256": entry.sha256, "etag": str(response.headers.get("ETag", "")),
                "last_modified": str(response.headers.get("Last-Modified", "")),
            }
        self._write_state(state)
        return self.status(online_enabled=True)
