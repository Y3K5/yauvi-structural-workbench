"""The campaign spec is data, and the builder refuses to guess."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from sf_csa.manifests import SpecError, build, build_database_manifest, build_target_manifest, load_spec
from conftest import write_spec


def test_a_valid_spec_loads(campaign):
    spec, root = load_spec(campaign)
    assert spec["release_scope"] == "one synthetic target"
    assert root.is_dir()


def test_build_produces_both_manifests(campaign, tmp_path):
    target_path, database_path = build(campaign, tmp_path / "out")
    assert target_path.is_file() and database_path.is_file()
    targets = json.loads(target_path.read_text(encoding="utf-8"))
    assert len(targets["queries"]) == 1
    assert targets["queries"][0]["accession"] == "EX0001"


def test_every_query_is_checksummed(campaign, tmp_path):
    target_path, _ = build(campaign, tmp_path / "out")
    query = json.loads(target_path.read_text(encoding="utf-8"))["queries"][0]
    assert len(query["structure_sha256"]) == 64
    assert len(query["sequence_sha256"]) == 64


def test_paths_are_recorded_relative_to_the_campaign_root(campaign, tmp_path):
    target_path, _ = build(campaign, tmp_path / "out")
    query = json.loads(target_path.read_text(encoding="utf-8"))["queries"][0]
    for key in ("fasta_path", "structure_path"):
        assert not query[key].startswith("/"), f"{key} is absolute"


def test_the_decision_status_comes_from_the_ledger(campaign, tmp_path):
    target_path, _ = build(campaign, tmp_path / "out")
    assert json.loads(target_path.read_text(encoding="utf-8"))["queries"][0]["decision_status"] == "SELECTED"


def test_the_default_orientation_artifact_is_applied(campaign, tmp_path):
    target_path, _ = build(campaign, tmp_path / "out")
    query = json.loads(target_path.read_text(encoding="utf-8"))["queries"][0]
    assert query["orientation_artifact"] == "results/orient_manifest.json"


def test_a_per_target_orientation_override_wins(campaign, tmp_path, spec_document):
    spec_document["targets"][0]["orientation_artifact"] = "results/special.json"
    write_spec(campaign, spec_document)
    target_path, _ = build(campaign, tmp_path / "out")
    query = json.loads(target_path.read_text(encoding="utf-8"))["queries"][0]
    assert query["orientation_artifact"] == "results/special.json"


# --- refusals ------------------------------------------------------------


def test_a_missing_spec_is_an_error(tmp_path):
    with pytest.raises(SpecError, match="not found"):
        load_spec(tmp_path / "absent.json")


def test_a_spec_with_no_targets_is_an_error(campaign, spec_document):
    spec_document["targets"] = []
    write_spec(campaign, spec_document)
    with pytest.raises(SpecError, match="declares no targets"):
        load_spec(campaign)


def test_a_target_missing_a_required_field_is_an_error(campaign, spec_document):
    del spec_document["targets"][0]["mechanism_group"]
    write_spec(campaign, spec_document)
    with pytest.raises(SpecError, match="mechanism_group"):
        load_spec(campaign)


def test_a_duplicated_target_is_an_error(campaign, spec_document):
    spec_document["targets"].append(dict(spec_document["targets"][0]))
    write_spec(campaign, spec_document)
    with pytest.raises(SpecError, match="declared twice"):
        load_spec(campaign)


def test_malformed_json_is_an_error(campaign):
    campaign.write_text("{not json", encoding="utf-8")
    with pytest.raises(SpecError, match="not valid JSON"):
        load_spec(campaign)


def test_sequence_ledger_drift_blocks_the_build(campaign, tmp_path):
    """The ledger is the authority; a changed FASTA must stop the build."""
    spec, root = load_spec(campaign)
    (root / "proteomes" / "example.faa").write_text(
        ">sp|EX0001|EXAMPLE_ORG\nMKVLAAGIVGLTTHAADQPRSTWYCHANGED\n", encoding="utf-8"
    )
    with pytest.raises(SpecError, match="sequence ledger drift"):
        build_target_manifest(spec, root)


def test_a_missing_structure_blocks_the_build(campaign):
    spec, root = load_spec(campaign)
    (root / "results" / "EX0001.pdb").unlink()
    with pytest.raises(SpecError, match="structure not found"):
        build_target_manifest(spec, root)


def test_an_accession_absent_from_the_ledger_blocks_the_build(campaign):
    spec, root = load_spec(campaign)
    (root / "results" / "SEQUENCE_MANIFEST.tsv").write_text(
        "accession\tsource_file\tsequence_sha256\n", encoding="utf-8"
    )
    with pytest.raises(SpecError, match="not in the sequence manifest"):
        build_target_manifest(spec, root)


def test_a_missing_structure_database_blocks_the_build(campaign):
    spec, root = load_spec(campaign)
    (root / "db" / "structdb").unlink()
    with pytest.raises(SpecError, match="structure database not found"):
        build_database_manifest(spec, root, query_count=1)


def test_an_unversioned_database_blocks_the_build(campaign):
    """An unversioned database cannot be cited, so it may not be used."""
    spec, root = load_spec(campaign)
    (root / "db" / "pdb.version").unlink()
    with pytest.raises(SpecError, match="no version recorded"):
        build_database_manifest(spec, root, query_count=1)


def test_a_spec_without_a_database_block_is_an_error(campaign, spec_document):
    del spec_document["database"]
    write_spec(campaign, spec_document)
    spec, root = load_spec(campaign)
    with pytest.raises(SpecError, match="no `database` block"):
        build_database_manifest(spec, root, query_count=1)


# --- what the manifest carries -------------------------------------------


def test_the_campaign_tables_are_written_into_the_manifest(campaign, tmp_path):
    """The rules that used to be literals in core.py must be recorded per release."""
    _, database_path = build(campaign, tmp_path / "out")
    document = json.loads(database_path.read_text(encoding="utf-8"))
    for key in ("mechanism_families", "contested_groups", "divergence_sets"):
        assert document[key], f"{key} is not recorded in the database manifest"


def test_a_spec_may_override_the_default_tables(campaign, tmp_path, spec_document):
    spec_document["database"]["contested_groups"] = [
        {"group": "my_group", "reason": "locally contested"}
    ]
    write_spec(campaign, spec_document)
    _, database_path = build(campaign, tmp_path / "out")
    document = json.loads(database_path.read_text(encoding="utf-8"))
    assert document["contested_groups"][0]["group"] == "my_group"


def test_the_database_checksum_and_version_are_recorded(campaign, tmp_path):
    _, database_path = build(campaign, tmp_path / "out")
    document = json.loads(database_path.read_text(encoding="utf-8"))
    assert len(document["pdb_database_checksum"]) == 64
    assert document["pdb_database_version"] == "2026-01-01; frozen snapshot"


def test_foldseek_database_checksum_binds_every_prefix_sidecar(campaign, tmp_path):
    spec, root = load_spec(campaign)
    sidecar = root / "db" / "structdb.index"
    sidecar.write_text("first index\n", encoding="utf-8")
    first = build_database_manifest(spec, root, query_count=1)
    assert sorted(first["pdb_database_file_checksums"]) == ["structdb", "structdb.index"]
    sidecar.write_text("drifted index\n", encoding="utf-8")
    second = build_database_manifest(spec, root, query_count=1)
    assert first["pdb_database_checksum"] != second["pdb_database_checksum"]


def test_the_query_count_expectation_matches_the_targets(campaign, tmp_path):
    _, database_path = build(campaign, tmp_path / "out")
    document = json.loads(database_path.read_text(encoding="utf-8"))
    assert document["release_expectations"]["query_count"] == 1


def test_the_closed_vocabulary_is_recorded_by_default(campaign, tmp_path):
    _, database_path = build(campaign, tmp_path / "out")
    document = json.loads(database_path.read_text(encoding="utf-8"))
    assert len(document["classification_vocabulary"]) == 6
