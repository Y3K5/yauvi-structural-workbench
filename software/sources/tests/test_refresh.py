from __future__ import annotations

from datetime import datetime, timezone

from yauvi_sources import ReferenceRefreshManager, SourceCache


class Response:
    def __init__(self, content: bytes, url: str, headers=None):
        self.content = content
        self.url = url
        self.headers = headers or {}


def test_refresh_is_opt_in_and_never_contacts_network_for_status(tmp_path, monkeypatch):
    def blocked(*args, **kwargs):
        raise AssertionError("network must not be contacted")

    monkeypatch.setattr("yauvi_sources.refresh._public_get", blocked)
    manager = ReferenceRefreshManager(tmp_path)
    status = manager.status(online_enabled=False)
    assert all(row["status"] == "offline" for row in status["sources"].values())


def test_refresh_caches_content_and_skips_recent_checks(tmp_path, monkeypatch):
    calls = []
    payloads = {
        "homd.org": b"HOMD download release",
        "proteinatlas.org": b"PK\x03\x04hpa",
        "expasy.org": b"\x1f\x8brhea",
        "proteomexchange.org": b'{"dataset":"PXD006367"}',
    }

    def fake(url, **kwargs):
        calls.append({"url": url, **kwargs})
        key = next(key for key in payloads if key in url)
        return Response(payloads[key], url, {"ETag": f'"{key}"'}), ""

    monkeypatch.setattr("yauvi_sources.refresh._public_get", fake)
    now = datetime(2026, 8, 22, 12, tzinfo=timezone.utc)
    manager = ReferenceRefreshManager(tmp_path, now=lambda: now)
    first = manager.refresh()
    assert len(calls) == 4
    assert all(row["status"] == "current" for row in first["sources"].values())
    manager.refresh()
    assert len(calls) == 4


def test_changed_cache_is_update_available_and_failure_uses_stale_cache(tmp_path, monkeypatch):
    cache = SourceCache(tmp_path)
    cache.store("rhea", b"\x1f\x8bold", filename="rhea-tsv.tar.gz", origin="old")
    def fake(url, **kwargs):
        if "expasy" in url:
            return Response(b"\x1f\x8bnew", url), ""
        return None, "network_unreachable:connection_error"

    monkeypatch.setattr("yauvi_sources.refresh._public_get", fake)
    manager = ReferenceRefreshManager(cache, now=lambda: datetime(2026, 8, 22, tzinfo=timezone.utc))
    status = manager.refresh(force=True)
    assert status["sources"]["rhea"]["status"] == "update_available"
    assert status["sources"]["rhea"]["sha256"] != ""
    assert status["sources"]["homd"]["status"] == "offline"


def test_refresh_sends_conditional_headers(tmp_path, monkeypatch):
    calls = []
    def fake(url, **kwargs):
        calls.append(kwargs.get("headers", {}))
        return Response(b"HOMD download" if "homd" in url else b"PK" if "proteinatlas" in url else b"\x1f\x8b" if "expasy" in url else b"PXD006367",
                        url, {"ETag": "v1", "Last-Modified": "date"}), ""
    monkeypatch.setattr("yauvi_sources.refresh._public_get", fake)
    manager = ReferenceRefreshManager(tmp_path, now=lambda: datetime(2026, 8, 22, tzinfo=timezone.utc))
    manager.refresh(force=True)
    manager.now = lambda: datetime(2026, 8, 24, tzinfo=timezone.utc)
    manager.refresh(force=True)
    assert all(headers.get("If-None-Match") == "v1" for headers in calls[4:])
