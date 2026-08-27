"""Unified orientor: routing, alpha-helical orienter (P2), localization seam, invariance."""

from __future__ import annotations

import os

import numpy as np
import pytest

from memorient.contexts import get_context
from memorient.orientor import (
    LocalizationCall,
    five_fold_validate,
    orient_structure,
)

from synthetic import make_barrel, make_tm_helix, make_soluble_blob

GN = get_context("gram_negative_om")
TM = get_context("tm_receptor")
SOL = get_context("soluble_secreted")


def tm_topology():
    return {
        "spans": [{"chain_id": "A", "start_auth_seq_id": 11, "end_auth_seq_id": 31}],
        "sidedness": {"extracellular_residue": {"chain_id": "A", "auth_seq_id": 41}},
        "source": {"id": "synthetic-topology", "citation": "test fixture"},
    }


# -- routing -------------------------------------------------------------------------------


def test_barrel_routes_to_barrel_normal_and_classifies():
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=10, peri_loop_len=2, seed=0)
    r = orient_structure(s, GN, n_points=120, validate=False)
    assert r.method == "barrel_normal"
    assert r.label == "barrel"
    assert r.fit is not None and r.projection is not None
    assert r.fit.delta_kd > 1.0
    # oriented frame: membrane centred near origin, extracellular residues on +Z
    ec = [l for l in r.labels.labels if l.extracellular]
    assert len(ec) > 0
    assert np.mean([r.structure.ca[i, 2] for i, l in enumerate(r.labels.labels) if l.extracellular]) > 0


def test_tm_helix_routes_to_belt_and_orients_by_positive_inside():
    h = make_tm_helix(seed=1)
    r = orient_structure(h, TM, n_points=120, validate=False, topology_evidence=tm_topology())
    assert r.method == "tm_helix_axis_v2"
    assert r.label == "tm_helix_experimental"
    assert r.scope_id == "alpha_helical"
    assert r.scientific_readiness == "prototype"
    assert r.scientific_state == "placement_evaluated"
    assert r.side.votes["declared_topology"] != 0


def test_tm_helix_without_spans_is_incomplete_and_has_no_sided_claim():
    h = make_tm_helix(seed=1)
    result = orient_structure(h, TM, n_points=120, validate=False)
    assert result.method == "tm_helix_belt_legacy_experimental"
    assert result.scientific_state == "insufficient_topology"
    assert result.extracellular_resids() == []
    assert result.labels.surface_set == []
    assert result.host_antibody_accessible is False


def test_topology_mapping_distinguishes_author_insertion_codes():
    helix = make_tm_helix(seed=1)
    helix.resids[11] = helix.resids[10]
    helix.icodes[11] = "A"
    residues = [
        {
            "chain_id": str(helix.chains[index]),
            "auth_seq_id": int(helix.resids[index]),
            "insertion_code": str(helix.icodes[index]),
        }
        for index in range(10, 31)
    ]
    topology = {
        "spans": [{"residues": residues}],
        "sidedness": {
            "extracellular_residue": {
                "chain_id": "A", "auth_seq_id": 41, "insertion_code": "",
            }
        },
        "source": {"id": "synthetic-topology", "citation": "test fixture"},
    }
    result = orient_structure(
        helix, TM, n_points=80, validate=False, topology_evidence=topology,
    )
    duplicated = int(helix.resids[10])
    mapped = [row for row in result.residue_table() if row["resid"] == duplicated]
    assert {row["insertion_code"] for row in mapped} == {"", "A"}
    exact_keys = result.to_dict()["surface_residue_keys"]
    assert all("insertion_code" in key for key in exact_keys)


def test_soluble_routes_to_sasa_only_and_has_no_membrane():
    b = make_soluble_blob(n_res=120, seed=2)
    r = orient_structure(b, SOL, n_points=120, validate=False)
    assert r.method == "sasa_only"
    assert r.label == "soluble"
    assert r.projection is None          # NO membrane zones for a soluble protein
    assert r.fit is None
    assert len(r.labels.surface_set) > 0  # exposed residues are the surface set


# -- P2: alpha-helical orienter has no lipid/pore metric ----------------------------------


def test_tm_context_does_not_compute_lipid_pore_gap():
    h = make_tm_helix(seed=1)
    r = orient_structure(h, TM, n_points=120, validate=False, topology_evidence=tm_topology())
    # tm_receptor context does not declare lipid_pore_gap -> it must not appear in metrics
    assert "lipid_pore_gap" not in r.metrics
    assert "hydrophobic_belt" in r.metrics


# -- P1 seam: localization can veto geometry ----------------------------------------------


def test_localization_veto_flips_host_accessibility():
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=16, peri_loop_len=2, seed=0)
    # geometry alone: surface-exposed
    r_open = orient_structure(s, GN, n_points=120, validate=False)
    # biology says periplasmic/shielded -> not host-accessible despite geometry
    veto = LocalizationCall(localization="periplasmic", surface_exposed=False, source="test")
    r_veto = orient_structure(s, GN, localization=veto, n_points=120, validate=False)
    assert r_veto.host_antibody_accessible is False
    if len(r_open.labels.surface_set) > 0:
        assert r_open.host_antibody_accessible is True


# -- rotation invariance ------------------------------------------------------------------


def test_synthetic_barrel_is_broadly_rotation_invariant():
    """A perfectly cylindrical synthetic barrel has degenerate in-plane axes (the geometry
    ruling): the normal can tilt slightly under rotation with nothing to pin it, so a handful
    of boundary residues flicker. We require the extracellular set to be broadly stable here
    and enforce the strict >=0.95 bar on a real OMP (below), which has loop asymmetry."""
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=10, peri_loop_len=2, seed=0)
    v = five_fold_validate(s, GN, n_points=160)
    assert v["mean_jaccard"] >= 0.90
    assert min(v["jaccards"]) >= 0.80


def test_soluble_is_rotation_invariant():
    b = make_soluble_blob(n_res=120, seed=2)
    v = five_fold_validate(b, SOL, n_points=160)
    assert v["passed"]                # soluble surface set is frame-independent
    assert v["mean_jaccard"] >= 0.95


def test_empty_extracellular_sets_are_not_a_passing_jaccard():
    helix = make_tm_helix(seed=1)
    validation = five_fold_validate(helix, TM, seeds=2, n_points=80)
    assert validation["extracellular_comparison_state"] == "not_applicable_empty_or_unresolved"
    assert validation["jaccards"] == [None, None]
    assert validation["sidedness_passed"] is False
    assert validation["passed"] is False


@pytest.mark.network
def test_real_omp_is_strictly_rotation_invariant():
    """Real OMPs (with genuine loop asymmetry) must hit the strict >=0.95 Jaccard bar."""
    import urllib.request
    from memorient.geometry import load_structure

    path = "/tmp/1bxw_test.pdb"
    if not os.path.exists(path):
        urllib.request.urlretrieve("https://files.rcsb.org/download/1BXW.pdb", path)
    s = load_structure(path, chain="A")
    v = five_fold_validate(s, GN, n_points=240)
    assert v["passed"]
    assert v["mean_jaccard"] >= 0.95
