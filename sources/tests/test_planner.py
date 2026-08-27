"""Joining requirements, policy, and cache into an answer."""
from __future__ import annotations

import pytest

from yauvi_sources import FetchClass, RegistryError, build_plan, render_plan
from yauvi_sources.manifest import parse_manifest

PAYLOAD = b">p\nMKVLAAG\n"


def test_plan_covers_every_declared_requirement(manifest, registry, cache):
    plan = build_plan(manifest, registry, cache)
    assert [i.source_id for i in plan.items] == manifest.source_ids()


def test_a_required_absent_source_blocks(manifest, registry, cache):
    plan = build_plan(manifest, registry, cache)
    assert not plan.satisfied
    assert [i.source_id for i in plan.blocking()] == ["open_db"]


def test_runtimes_and_heuristics_never_block(manifest, registry, cache):
    """A binary is resolved by runtime preflight, not by this layer."""
    plan = build_plan(manifest, registry, cache)
    blocking = {i.source_id for i in plan.blocking()}
    assert "a_binary" not in blocking
    assert "a_heuristic" not in blocking


def test_caching_the_required_source_satisfies_the_plan(manifest, registry, cache):
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o")
    plan = build_plan(manifest, registry, cache)
    assert plan.satisfied
    assert plan.blocking() == []


def test_optional_absent_sources_do_not_block(registry, cache):
    manifest = parse_manifest(
        {"module_id": "m", "requires": [{"source_id": "gated_db", "required": False}]}
    )
    assert build_plan(manifest, registry, cache).satisfied


def test_fetchable_lists_only_uncached_open_sources(manifest, registry, cache):
    assert [i.source_id for i in build_plan(manifest, registry, cache).fetchable()] == ["open_db"]
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o")
    assert build_plan(manifest, registry, cache).fetchable() == []


def test_manual_lists_what_a_human_must_supply(manifest, registry, cache):
    manual = {i.source_id for i in build_plan(manifest, registry, cache).manual()}
    assert manual == {"gated_db", "a_table"}


def test_a_staged_source_stops_being_manual(manifest, registry, cache):
    cache.store("gated_db", PAYLOAD, filename="deg.faa", origin="staged:/x")
    manual = {i.source_id for i in build_plan(manifest, registry, cache).manual()}
    assert "gated_db" not in manual


def test_item_status_words(manifest, registry, cache):
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o")
    statuses = {i.source_id: i.status for i in build_plan(manifest, registry, cache).items}
    assert statuses == {
        "open_db": "present",
        "gated_db": "manual",
        "a_table": "manual",
        "a_binary": "runtime",
        "a_heuristic": "internal",
    }


def test_plan_rejects_a_module_naming_an_undeclared_source(registry, cache):
    manifest = parse_manifest({"module_id": "m", "requires": ["not_in_registry"]})
    with pytest.raises(RegistryError, match="not declared"):
        build_plan(manifest, registry, cache)


def test_render_names_the_blocking_sources(manifest, registry, cache):
    text = render_plan(build_plan(manifest, registry, cache))
    assert "NOT SATISFIED" in text
    assert "open_db" in text
    # The manual instructions must be printed, not merely counted.
    assert "yauvi-fetch stage gated_db" in text


def test_render_reports_satisfaction_and_the_recorded_digest(manifest, registry, cache):
    cache.store("open_db", PAYLOAD, filename="db.fasta", origin="o", version="2026_03")
    text = render_plan(build_plan(manifest, registry, cache), verbose=True)
    assert "SATISFIED" in text
    assert "version: 2026_03" in text
    assert "CC BY 4.0." in text


def test_fetch_class_may_download_is_true_only_for_open(registry, cache, manifest):
    for item in build_plan(manifest, registry, cache).items:
        assert item.can_fetch == (item.fetch_class is FetchClass.OPEN_FETCHABLE)
