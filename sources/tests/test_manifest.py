"""A module's declaration of what it needs."""
from __future__ import annotations

import pytest

from yauvi_sources import ManifestError, load_manifest, resolve_manifest
from yauvi_sources.manifest import parse_manifest


def test_parses_full_requirements(manifest):
    assert manifest.module_id == "testmod"
    assert manifest.source_ids()[0] == "open_db"
    assert len(manifest.required_only()) == 2


def test_string_shorthand_is_accepted():
    parsed = parse_manifest({"module_id": "m", "requires": ["a", "b"]})
    assert parsed.source_ids() == ["a", "b"]
    assert all(r.required for r in parsed.requires)


def test_missing_module_id_is_an_error():
    with pytest.raises(ManifestError, match="missing module_id"):
        parse_manifest({"requires": []})


def test_duplicate_requirement_is_an_error():
    with pytest.raises(ManifestError, match="twice"):
        parse_manifest({"module_id": "m", "requires": ["a", "a"]})


def test_requirement_without_source_id_is_an_error():
    with pytest.raises(ManifestError, match="without source_id"):
        parse_manifest({"module_id": "m", "requires": [{"role": "orphan"}]})


def test_requires_must_be_a_list():
    with pytest.raises(ManifestError, match="must be a list"):
        parse_manifest({"module_id": "m", "requires": {"a": 1}})


def test_a_module_may_declare_nothing():
    assert parse_manifest({"module_id": "m"}).source_ids() == []


def test_load_from_disk(tmp_path):
    path = tmp_path / "sources.yaml"
    path.write_text("module_id: m\nrequires:\n  - source_id: a\n", encoding="utf-8")
    assert load_manifest(path).module_id == "m"


def test_missing_manifest_file_is_an_error(tmp_path):
    with pytest.raises(ManifestError, match="not found"):
        load_manifest(tmp_path / "absent.yaml")


def test_resolve_reports_where_to_look_when_nothing_is_found(tmp_path):
    with pytest.raises(ManifestError, match="no source manifest found"):
        resolve_manifest("not_a_module", workspace=tmp_path)


def test_resolve_prefers_an_explicit_path(tmp_path):
    path = tmp_path / "custom.yaml"
    path.write_text("module_id: custom\nrequires: []\n", encoding="utf-8")
    assert resolve_manifest("subproteo", explicit_path=path).module_id == "custom"


# --- the shipped manifests ----------------------------------------------


@pytest.mark.parametrize("module_id", ["subproteo", "memorient"])
def test_shipped_manifests_resolve_against_the_real_registry(
    module_id, workspace, real_registry
):
    # Skip rather than fail when the module is simply absent: this suite runs in
    # environments where only yauvi-sources is installed, and "that package is
    # not here" is not a finding about this package.
    try:
        manifest = resolve_manifest(module_id, workspace=workspace)
    except ManifestError:
        pytest.skip(f"{module_id} is neither installed nor present in this tree")
    # Raises if the module names a source the registry does not declare.
    real_registry.resolve_many(manifest.source_ids())
    assert manifest.module_id == module_id
    assert manifest.source_ids(), f"{module_id} declares no sources"
