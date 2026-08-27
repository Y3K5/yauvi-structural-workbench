"""The registry loader refuses to guess."""
from __future__ import annotations

import textwrap

import pytest

from yauvi_sources import RegistryError, SourceRegistry


def test_loads_every_synthetic_source(registry):
    assert len(registry) == 6
    assert "open_db" in registry
    assert registry.catalog_id == "test-registry"


def test_get_returns_the_declared_fields(registry):
    source = registry.get("open_db")
    assert source.display_name == "An Open Database"
    assert source.channel == "conservation"
    assert source.license_note == "CC BY 4.0."
    assert source.url == "https://example.invalid/open_db.fasta"


def test_unknown_source_raises_rather_than_returning_none(registry):
    # A typo must not become a silently-unfetched input.
    with pytest.raises(RegistryError, match="not declared"):
        registry.get("open_bd")


def test_resolve_many_reports_every_miss_at_once(registry):
    with pytest.raises(RegistryError) as excinfo:
        registry.resolve_many(["open_db", "nope_one", "nope_two"])
    message = str(excinfo.value)
    assert "nope_one" in message and "nope_two" in message


def _write(tmp_path, body: str):
    path = tmp_path / "sources.yaml"
    path.write_text(textwrap.dedent(body).strip(), encoding="utf-8")
    return path


def test_missing_required_field_is_an_error(tmp_path):
    path = _write(
        tmp_path,
        """
        sources:
          - source_id: half_declared
            display_name: Missing Its Access Mode
            kind: sequence_db
            channel: safety
            status: wired
        """,
    )
    with pytest.raises(RegistryError, match="missing required field"):
        SourceRegistry.load(path)


def test_unknown_status_is_an_error(tmp_path):
    # Drift between the file's vocabulary and this loader must be loud.
    path = _write(
        tmp_path,
        """
        sources:
          - source_id: odd
            display_name: Odd Status
            kind: sequence_db
            channel: safety
            status: probably_fine
            access: internal
        """,
    )
    with pytest.raises(RegistryError, match="unknown status"):
        SourceRegistry.load(path)


def test_duplicate_source_id_is_an_error(tmp_path):
    path = _write(
        tmp_path,
        """
        sources:
          - source_id: twice
            display_name: First
            kind: sequence_db
            channel: safety
            status: wired
            access: internal
          - source_id: twice
            display_name: Second
            kind: sequence_db
            channel: safety
            status: wired
            access: internal
        """,
    )
    with pytest.raises(RegistryError, match="duplicate"):
        SourceRegistry.load(path)


def test_empty_registry_is_an_error(tmp_path):
    path = _write(tmp_path, "sources: []")
    with pytest.raises(RegistryError, match="declares no sources"):
        SourceRegistry.load(path)


def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(RegistryError, match="not found"):
        SourceRegistry.load(tmp_path / "absent.yaml")


# --- the real file -------------------------------------------------------


def test_real_registry_loads(real_registry):
    """The whole point of this package: catalogs/sources.yaml must stay readable."""
    assert len(real_registry) >= 30
    assert real_registry.catalog_id == "yauvi-evidence-sources"


def test_real_registry_sources_all_carry_provenance(real_registry):
    for source in real_registry:
        assert source.display_name, f"{source.source_id} has no display name"
        assert source.channel, f"{source.source_id} has no channel"
        assert source.access, f"{source.source_id} has no access mode"


def test_real_registry_key_sources_present(real_registry):
    for source_id in ("uniprot_proteomes", "deg", "diamond", "iedb", "alphafold_db", "pdb"):
        assert source_id in real_registry, f"{source_id} disappeared from the registry"
