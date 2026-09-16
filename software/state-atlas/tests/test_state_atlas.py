from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from state_atlas.core import InputError, analyze, validate_reference_set, write_outputs


BASE = [(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2)]
INACTIVE = [(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 5)]


def pdb(points, models=None):
    collections = models or [points]; lines = []
    for mi, coords in enumerate(collections, start=1):
        if len(collections) > 1: lines.append(f"MODEL     {mi:4d}")
        for i, (x, y, z) in enumerate(coords, start=1):
            lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
        if len(collections) > 1: lines.append("ENDMDL")
    lines.append("END"); return "\n".join(lines) + "\n"


def setup(tmp_path: Path):
    active = tmp_path / "active.pdb"; inactive = tmp_path / "inactive.pdb"; query = tmp_path / "query.pdb"
    active.write_text(pdb(BASE)); inactive.write_text(pdb(INACTIVE)); query.write_text(pdb(BASE, [BASE, INACTIVE]))
    refs = {"reference_set_id": "synthetic-two-state", "decision_rules": {"max_rmsd_A": 0.5, "min_margin_A": 0.2},
            "references": [
                {"reference_id": "ACTIVE", "state": "active", "structure": "active.pdb",
                 "provenance": {"class": "experimental", "method": "synthetic coordinates"},
                 "state_evidence": {"basis": "synthetic compact fourth residue"}},
                {"reference_id": "INACTIVE", "state": "inactive", "structure": "inactive.pdb",
                 "provenance": {"class": "experimental", "method": "synthetic coordinates"},
                 "state_evidence": {"basis": "synthetic displaced fourth residue"}},
            ]}
    digest = hashlib.sha256(query.read_bytes()).hexdigest()
    manifest = {"subject": {"id": "Q"}, "coordinate": {"sha256": digest}}
    return query, refs, manifest


def test_mixed_static_ensemble(tmp_path):
    query, refs, manifest = setup(tmp_path)
    doc = analyze(manifest, refs, reference_base=tmp_path, structure_path=query, cluster_cutoff_A=0.5)
    assert doc["overall_label"] == "mixed"
    assert doc["frames_total"] == 2
    assert doc["populations"]["unresolved"]["count"] == 0
    out = tmp_path / "out"; write_outputs(out, doc)
    assert (out / "STATE_LAYER.json").is_file()


def test_unresolved_frames_remain_in_total_population(tmp_path):
    query, refs, manifest = setup(tmp_path)
    midpoint = [tuple((np.asarray(a) + np.asarray(b)) / 2) for a, b in zip(BASE, INACTIVE)]
    query.write_text(pdb(BASE, [BASE, INACTIVE, midpoint]))
    manifest["coordinate"]["sha256"] = hashlib.sha256(query.read_bytes()).hexdigest()
    document = analyze(manifest, refs, reference_base=tmp_path, structure_path=query, cluster_cutoff_A=0.5)
    assert [row["call"] for row in document["frame_metrics"]] == ["active_like", "inactive_like", "unresolved"]
    assert document["frames_total"] == 3
    assert document["frames_interpretable"] == 2
    assert document["populations"]["unresolved"]["fraction_total"] == pytest.approx(1 / 3, abs=1e-6)


def test_one_sided_reference_set_is_refused(tmp_path):
    _query, refs, _manifest = setup(tmp_path)
    refs["references"] = refs["references"][:1]
    assert any("active and inactive" in e for e in validate_reference_set(refs, base=tmp_path))


def test_trajectory_requires_explicit_topology_and_pbc(tmp_path):
    query, refs, manifest = setup(tmp_path)
    fake = tmp_path / "run.xtc"; fake.write_bytes(b"not a trajectory")
    with pytest.raises(InputError, match="topology and explicit"):
        analyze(manifest, refs, reference_base=tmp_path, trajectory_path=fake)


def test_rigid_transform_and_model_order_preserve_population_calls(tmp_path):
    query, refs, manifest = setup(tmp_path)
    baseline = analyze(manifest, refs, reference_base=tmp_path, structure_path=query, cluster_cutoff_A=0.5)
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    shift = np.asarray([7.0, -4.0, 3.0])
    active_moved = [tuple(np.asarray(point) @ rotation + shift) for point in BASE]
    inactive_moved = [tuple(np.asarray(point) @ rotation + shift) for point in INACTIVE]
    moved = tmp_path / "moved.pdb"; moved.write_text(pdb(active_moved, [active_moved, inactive_moved]))
    moved_manifest = {"subject": {"id": "Q"}, "coordinate": {"sha256": hashlib.sha256(moved.read_bytes()).hexdigest()}}
    transformed = analyze(moved_manifest, refs, reference_base=tmp_path, structure_path=moved, cluster_cutoff_A=0.5)

    reversed_models = tmp_path / "reversed.pdb"; reversed_models.write_text(pdb(INACTIVE, [INACTIVE, BASE]))
    reversed_manifest = {"subject": {"id": "Q"}, "coordinate": {"sha256": hashlib.sha256(reversed_models.read_bytes()).hexdigest()}}
    reordered = analyze(reversed_manifest, refs, reference_base=tmp_path, structure_path=reversed_models, cluster_cutoff_A=0.5)
    assert transformed["populations"] == baseline["populations"]
    assert transformed["overall_label"] == baseline["overall_label"]
    assert reordered["populations"] == baseline["populations"]
    assert reordered["overall_label"] == baseline["overall_label"]


def _v2_fixture(tmp_path: Path):
    active_points = [
        (1.2 * (index % 17), 1.2 * (index // 17), 0.3 * ((index * 7) % 11))
        for index in range(254)
    ]
    inactive_points = [
        (x + (2.5 if index >= 190 else 0.0), y, z + (4.0 if index >= 190 else 0.0))
        for index, (x, y, z) in enumerate(active_points)
    ]
    active = tmp_path / "active-v2.pdb"
    active_two = tmp_path / "active-v2-two.pdb"
    inactive = tmp_path / "inactive-v2.pdb"
    inactive_two = tmp_path / "inactive-v2-two.pdb"
    query = tmp_path / "long-query.pdb"
    active.write_text(pdb(active_points), encoding="utf-8")
    active_two.write_text(pdb(active_points), encoding="utf-8")
    inactive.write_text(pdb(inactive_points), encoding="utf-8")
    inactive_two.write_text(pdb(inactive_points), encoding="utf-8")
    query_points = [(80, 80, 80), *active_points, (-80, -80, -80)]
    query.write_text(pdb(query_points), encoding="utf-8")
    rows_query = [
        {"uniprot_position": 242 + index, "chain_id": "A", "auth_seq_id": 2 + index, "mapping_state": "exact"}
        for index in range(254)
    ]
    rows_reference = [
        {"uniprot_position": 242 + index, "chain_id": "A", "auth_seq_id": 1 + index, "mapping_state": "exact"}
        for index in range(254)
    ]
    alignment = {
        "schema_version": "1.0", "coordinate_system": "uniprot",
        "source": {"id": "SIFTS:synthetic", "citation": "synthetic exact mapping", "sha256": "1" * 64},
        "domain": {"uniprot_start": 242, "uniprot_end": 495},
        "query": rows_query,
        "references": {
            "ACTIVE": rows_reference, "ACTIVE_TWO": rows_reference,
            "INACTIVE": rows_reference, "INACTIVE_TWO": rows_reference,
        },
    }
    alignment_path = tmp_path / "alignment-map.json"
    alignment_path.write_text(json.dumps(alignment, sort_keys=True) + "\n", encoding="utf-8")
    alignment_digest = hashlib.sha256(alignment_path.read_bytes()).hexdigest()
    references = {
        "schema_version": "2.0", "reference_set_id": "synthetic-abl-v2",
        "qualification_scope": "abl_family", "subject_family": "ABL1",
        "alignment_map_sha256": alignment_digest,
        "alignment_mask": {"coordinate_system": "uniprot", "uniprot_start": 242, "uniprot_end": 495, "minimum_coverage": 0.9},
        "decision_rules": {"max_rmsd_A": 2.5, "min_margin_A": 0.25},
        "references": [
            {"reference_id": "ACTIVE", "state": "active", "structure": active.name,
             "structure_sha256": hashlib.sha256(active.read_bytes()).hexdigest(), "chain": "A", "pdb_entry_id": "SYN1",
             "provenance": {"class": "experimental", "method": "synthetic coordinates"},
             "state_evidence": {"basis": "compact fourth residue", "citation": "synthetic active evidence"}},
            {"reference_id": "ACTIVE_TWO", "state": "active", "structure": active_two.name,
             "structure_sha256": hashlib.sha256(active_two.read_bytes()).hexdigest(), "chain": "A", "pdb_entry_id": "SYN2",
             "provenance": {"class": "experimental", "method": "synthetic coordinates"},
             "state_evidence": {"basis": "compact domain control", "citation": "synthetic active evidence"}},
            {"reference_id": "INACTIVE", "state": "inactive", "structure": inactive.name,
             "structure_sha256": hashlib.sha256(inactive.read_bytes()).hexdigest(), "chain": "A", "pdb_entry_id": "SYN3",
             "provenance": {"class": "experimental", "method": "synthetic coordinates"},
             "state_evidence": {"basis": "extended fourth residue", "citation": "synthetic inactive evidence"}},
            {"reference_id": "INACTIVE_TWO", "state": "inactive", "structure": inactive_two.name,
             "structure_sha256": hashlib.sha256(inactive_two.read_bytes()).hexdigest(), "chain": "A", "pdb_entry_id": "SYN4",
             "provenance": {"class": "experimental", "method": "synthetic coordinates"},
             "state_evidence": {"basis": "extended domain control", "citation": "synthetic inactive evidence"}},
        ],
    }
    manifest = {"subject": {"id": "ABL1_SYNTHETIC"}, "coordinate": {"sha256": hashlib.sha256(query.read_bytes()).hexdigest()}}
    return query, references, alignment, alignment_digest, manifest


def test_reference_set_v2_uses_exact_domain_map_and_excludes_construct_extensions(tmp_path):
    query, references, alignment, digest, manifest = _v2_fixture(tmp_path)
    document = analyze(
        manifest, references, reference_base=tmp_path, structure_path=query, chain="A",
        alignment_map=alignment, alignment_map_sha256=digest, cluster_cutoff_A=0.5,
    )
    assert document["overall_label"] == "active_like"
    assert document["scientific_scope"]["scope_id"] == "abl_family"
    assert document["scientific_scope"]["scientific_state"] == "prototype"
    assert document["scientific_scope"]["qualification_gate"] == "qualification_v2_held_out_panel_pending"
    for reference in document["reference_set"]["references"]:
        assert reference["mapping"]["positions_mapped"] == 254
        assert reference["mapping"]["coverage"] == 1.0
    assert document["config"]["ensemble_alignment_residue_count"] == 254


def test_reference_set_v2_requires_checksum_bound_map(tmp_path):
    _query, references, alignment, digest, _manifest = _v2_fixture(tmp_path)
    assert any("checksum" in error for error in validate_reference_set(
        references, base=tmp_path, alignment_map=alignment, alignment_map_sha256="0" * 64,
    ))
    assert any("requires an alignment map" in error for error in validate_reference_set(references, base=tmp_path))


def test_reference_set_v2_freezes_abl_scope_and_requires_multiple_references(tmp_path):
    _query, references, alignment, digest, _manifest = _v2_fixture(tmp_path)
    narrowed = copy.deepcopy(references)
    narrowed["references"] = [
        next(row for row in narrowed["references"] if row["reference_id"] == "ACTIVE"),
        next(row for row in narrowed["references"] if row["reference_id"] == "INACTIVE"),
    ]
    errors = validate_reference_set(
        narrowed, base=tmp_path, alignment_map=alignment, alignment_map_sha256=digest,
    )
    assert sum("requires multiple" in error for error in errors) == 2

    changed = copy.deepcopy(references)
    changed["alignment_mask"]["uniprot_end"] = 494
    changed["decision_rules"]["max_rmsd_A"] = 3.0
    errors = validate_reference_set(
        changed, base=tmp_path, alignment_map=alignment, alignment_map_sha256=digest,
    )
    assert any("242-495" in error for error in errors)
    assert any("frozen at 2.5" in error for error in errors)


def test_reference_set_v2_reports_ambiguous_positions_without_using_them(tmp_path):
    query, references, alignment, _digest, manifest = _v2_fixture(tmp_path)
    alignment["query"][0]["mapping_state"] = "ambiguous"
    alignment_path = tmp_path / "ambiguous-alignment-map.json"
    alignment_path.write_text(json.dumps(alignment, sort_keys=True) + "\n", encoding="utf-8")
    digest = hashlib.sha256(alignment_path.read_bytes()).hexdigest()
    references["alignment_map_sha256"] = digest
    document = analyze(
        manifest, references, reference_base=tmp_path, structure_path=query, chain="A",
        alignment_map=alignment, alignment_map_sha256=digest, cluster_cutoff_A=0.5,
    )
    for reference in document["reference_set"]["references"]:
        mapping = reference["mapping"]
        assert mapping["positions_mapped"] == 253
        assert mapping["ambiguous_uniprot_positions"] == [242]
        assert 242 not in mapping["included_uniprot_positions"]


def test_static_selection_is_rejected_instead_of_silently_ignored(tmp_path):
    query, refs, manifest = setup(tmp_path)
    with pytest.raises(InputError, match="applies only to trajectories"):
        analyze(
            manifest, refs, reference_base=tmp_path, structure_path=query,
            selection="protein and resid 2:4", cluster_cutoff_A=0.5,
        )
