"""Membrane-zone + context-gating tests."""

from __future__ import annotations

import numpy as np
import pytest

from memorient.barrel import fit_membrane
from memorient.contexts import Metric, get_context
from memorient.membrane import (
    ACC_ANTIBODY,
    ACC_LPS_SHIELDED,
    ZONE_CORE,
    ZONE_EXTRACELLULAR,
    ZONE_PERIPLASMIC,
    context_metrics,
    project_membrane,
)
from memorient.sasa import compute_sasa

from synthetic import make_barrel

GN = get_context("gram_negative_om")
EUK = get_context("eukaryotic_pm")
SOL = get_context("soluble_secreted")


def _fit_and_project(ctx, ec_sign=1):
    s = make_barrel(n_strands=12, strand_len=10, seed=0)
    fit = fit_membrane(s, ctx)
    rsa = compute_sasa(s, n_points=120)["rsa"]
    proj = project_membrane(s, fit, ctx, ec_sign=ec_sign, rsa=rsa)
    return s, fit, proj


def test_zones_span_the_membrane():
    s, fit, proj = _fit_and_project(GN)
    zones = set(proj.zone.tolist())
    assert ZONE_CORE in zones
    # both extracellular and periplasmic zones present in a full barrel
    assert ZONE_EXTRACELLULAR in zones
    assert ZONE_PERIPLASMIC in zones


def test_core_residues_have_facing():
    s, fit, proj = _fit_and_project(GN)
    core = proj.zone == ZONE_CORE
    facings = set(proj.facing[core].tolist())
    assert "lipid" in facings and "pore" in facings


def test_soluble_context_refuses_membrane_projection():
    s = make_barrel(seed=1)
    fit = fit_membrane(s, GN)  # fit exists, but soluble context must refuse zones
    with pytest.raises(ValueError):
        project_membrane(s, fit, SOL, ec_sign=1)


def test_lps_shielding_only_in_gram_negative():
    # gram-negative OM: some extracellular-interface residues are LPS-shielded
    _, _, proj_gn = _fit_and_project(GN)
    assert np.any(proj_gn.accessibility == ACC_LPS_SHIELDED)

    # eukaryotic PM: no LPS band, so no residue is LPS-shielded
    _, _, proj_euk = _fit_and_project(EUK)
    assert not np.any(proj_euk.accessibility == ACC_LPS_SHIELDED)


def test_context_metrics_are_gated():
    s, fit, proj = _fit_and_project(GN)
    m_gn = context_metrics(s, fit, GN, proj)
    # gram-negative: has lipid/pore gap + LPS shielding, NOT positive-inside
    assert Metric.LIPID_PORE_GAP in m_gn
    assert Metric.LPS_SHIELDING in m_gn
    assert Metric.POSITIVE_INSIDE not in m_gn

    s2, fit2, proj2 = _fit_and_project(EUK)
    m_euk = context_metrics(s2, fit2, EUK, proj2)
    # eukaryotic PM: has positive-inside, NOT LPS shielding
    assert Metric.POSITIVE_INSIDE in m_euk
    assert Metric.LPS_SHIELDING not in m_euk


def test_barrel_signature_metrics_have_right_sign():
    s, fit, proj = _fit_and_project(GN)
    m = context_metrics(s, fit, GN, proj)
    # lipid face more hydrophobic than pore (positive gap)
    assert m[Metric.LIPID_PORE_GAP] > 0
    # embedded belt more hydrophobic than flanks
    assert m[Metric.HYDROPHOBIC_BELT] > 0
    # aromatic girdle enriched at interface (>= 0, planted in synthetic)
    assert m[Metric.AROMATIC_GIRDLE] >= 0
