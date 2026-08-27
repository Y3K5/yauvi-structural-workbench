"""Fetch policy: what code may and may not go and get."""
from __future__ import annotations

import pytest

from yauvi_sources import FetchClass, PolicyError, classify, instructions_for
from yauvi_sources.policy import ACCESS_TO_CLASS
from yauvi_sources.registry import Source


def test_open_source_is_fetchable(registry):
    assert classify(registry.get("open_db")) is FetchClass.OPEN_FETCHABLE
    assert classify(registry.get("open_db")).may_download


def test_licence_gated_source_is_never_downloadable(registry):
    fetch_class = classify(registry.get("gated_db"))
    assert fetch_class is FetchClass.LICENSE_GATED
    assert not fetch_class.may_download


def test_binary_is_a_runtime_not_a_file(registry):
    assert classify(registry.get("a_binary")) is FetchClass.RUNTIME


def test_internal_heuristic_has_nothing_to_fetch(registry):
    assert classify(registry.get("a_heuristic")) is FetchClass.INTERNAL


def test_table_only_status_overrides_a_fetchable_transport(registry):
    """A source reachable over an API is still not automated when its status says so.

    This is the IEDB case: `access: network_api_and_download` would make it
    downloadable on transport alone, but the registry records its status as
    `table_only` and states that no code path fetches it.
    """
    source = registry.get("reachable_over_api")
    assert source.access == "network_api_and_download"
    assert ACCESS_TO_CLASS[source.access] is FetchClass.OPEN_FETCHABLE
    assert classify(source) is FetchClass.TABLE_ONLY
    assert not classify(source).may_download


def test_unmapped_access_mode_raises_instead_of_defaulting():
    """An access mode nobody has decided about must not become downloadable."""
    invented = Source(
        source_id="mystery",
        display_name="Mystery",
        kind="sequence_db",
        channel="safety",
        status="proposed",
        access="carrier_pigeon",
    )
    with pytest.raises(PolicyError, match="no fetch policy"):
        classify(invented)


def test_licence_gated_instructions_name_the_staging_command(registry):
    text = instructions_for(registry.get("gated_db"))
    assert "yauvi-fetch stage gated_db" in text
    assert "Registration required." in text


def test_table_only_instructions_describe_the_expected_export(registry):
    text = instructions_for(registry.get("a_table"))
    assert "never produces" in text
    assert "yauvi-fetch stage a_table" in text


def test_runtime_instructions_point_at_the_runtime_registry(registry):
    assert "runtime-registry.yaml" in instructions_for(registry.get("a_binary"))


# --- the real file -------------------------------------------------------


def test_every_real_source_has_a_policy(real_registry):
    """No source in the shipped registry may be un-classifiable."""
    for source in real_registry:
        classify(source)  # raises PolicyError if unmapped


def test_real_licence_gated_sources_are_not_downloadable(real_registry):
    for source_id in ("deg", "drugbank"):
        assert not classify(real_registry.get(source_id)).may_download


def test_iedb_is_not_automated(real_registry):
    assert classify(real_registry.get("iedb")) is FetchClass.TABLE_ONLY
