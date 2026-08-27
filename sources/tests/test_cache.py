"""The cache is content-addressed, and it notices when a file changes."""
from __future__ import annotations

import pytest

from yauvi_sources import CacheError, SourceCache
from yauvi_sources.cache import default_cache_dir

PAYLOAD = b">p1\nMKVLAAG\n>p2\nMQQPTRA\n"


def test_store_puts_content_under_its_digest(cache):
    entry = cache.store("open_db", PAYLOAD, filename="db.fasta", origin="https://example.invalid/x")
    assert entry.sha256 == __import__("hashlib").sha256(PAYLOAD).hexdigest()
    assert entry.bytes == len(PAYLOAD)
    assert cache.path_for(entry).read_bytes() == PAYLOAD
    assert entry.sha256 in str(cache.path_for(entry))


def test_differing_content_never_overwrites(cache):
    first = cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o1")
    second = cache.store("open_db", PAYLOAD + b">p3\nMMM\n", filename="db.fasta", origin="o2")
    assert first.sha256 != second.sha256
    # Both remain retrievable: an updated upstream cannot destroy the copy an
    # earlier run was built on.
    assert cache.path_for(first).is_file()
    assert cache.path_for(second).is_file()
    assert len(cache.entries("open_db")) == 2


def test_reacquiring_identical_content_does_not_duplicate_the_record(cache):
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o1")
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o1")
    assert len(cache.entries("open_db")) == 1


def test_latest_returns_the_most_recent_present_entry(cache):
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o1")
    second = cache.store("open_db", PAYLOAD + b"x", filename="db.fasta", origin="o2")
    assert cache.latest("open_db").sha256 == second.sha256


def test_latest_skips_entries_whose_file_was_deleted(cache):
    first = cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o1")
    second = cache.store("open_db", PAYLOAD + b"x", filename="db.fasta", origin="o2")
    cache.path_for(second).unlink()
    assert cache.latest("open_db").sha256 == first.sha256


def test_empty_payload_is_refused(cache):
    with pytest.raises(CacheError, match="empty payload"):
        cache.store("open_db", b"", filename="db.fasta", origin="o")


def test_has_is_false_before_anything_is_stored(cache):
    assert not cache.has("open_db")
    assert cache.latest("open_db") is None


def test_stage_adopts_a_manual_file_and_records_its_origin(cache, tmp_path):
    path = tmp_path / "deg_bacteria.faa"
    path.write_bytes(PAYLOAD)
    entry = cache.stage("gated_db", path, note="downloaded after registering")
    assert entry.filename == "deg_bacteria.faa"
    assert entry.origin == f"staged:{path}"
    assert entry.note == "downloaded after registering"
    assert cache.path_for(entry).read_bytes() == PAYLOAD


def test_staging_a_missing_file_is_an_error(cache, tmp_path):
    with pytest.raises(CacheError, match="no such file"):
        cache.stage("gated_db", tmp_path / "absent.faa")


def test_verify_passes_on_an_untouched_cache(cache):
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o")
    results = list(cache.verify())
    assert results and all(ok for _, ok, _ in results)


def test_verify_detects_a_modified_file(cache):
    entry = cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o")
    path = cache.path_for(entry)
    path.write_bytes(PAYLOAD + b"tampered")
    results = list(cache.verify())
    assert [detail for _, ok, detail in results if not ok][0].startswith("digest mismatch")


def test_verify_detects_a_deleted_file(cache):
    entry = cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o")
    cache.path_for(entry).unlink()
    assert [detail for _, ok, detail in cache.verify() if not ok] == ["missing from disk"]


def test_manifest_records_the_upstream_version_when_one_is_reported(cache):
    entry = cache.store(
        "open_db", PAYLOAD, filename="db.fasta", origin="o", version="2026_03"
    )
    assert cache.entries("open_db")[0].version == "2026_03"
    assert entry.retrieved_at.endswith("Z")


def test_source_ids_lists_only_sources_with_a_manifest(cache):
    cache.store("open_db", PAYLOAD, filename="a", origin="o")
    cache.store("gated_db", PAYLOAD, filename="b", origin="o")
    assert cache.source_ids() == ["gated_db", "open_db"]


def test_cache_root_honours_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("YAUVI_SOURCE_CACHE", str(tmp_path / "elsewhere"))
    assert default_cache_dir() == tmp_path / "elsewhere"
    assert SourceCache().root == tmp_path / "elsewhere"
