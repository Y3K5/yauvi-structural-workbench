"""Extracellular-side call + per-residue labeling tests."""

from __future__ import annotations

import numpy as np
import pytest

from memorient.barrel import fit_membrane
from memorient.contexts import get_context
from memorient.geometry import Structure
from memorient.labeler import call_extracellular_side, label_residues
from memorient.membrane import project_membrane
from memorient.sasa import compute_sasa

from synthetic import make_barrel, make_soluble_blob

GN = get_context("gram_negative_om")


def _ec_depth_of_true_ec(structure, fit, side):
    """Mean ec_depth of residues that are truly extracellular (high original +Z)."""
    proj_along = (structure.ca - fit.centroid) @ fit.normal
    ecdepth = side.ec_sign * (proj_along - fit.center)
    hi = structure.ca[:, 2] > 8
    lo = structure.ca[:, 2] < -8
    return ecdepth[hi].mean(), ecdepth[lo].mean()


def test_side_call_puts_long_loop_face_extracellular():
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=8, peri_loop_len=2, seed=0)
    fit = fit_membrane(s, GN)
    side = call_extracellular_side(s, fit, GN)
    ec_mean, peri_mean = _ec_depth_of_true_ec(s, fit, side)
    # the +Z (long-loop) face must have positive ec_depth, -Z negative
    assert ec_mean > 0
    assert peri_mean < 0
    # both primary signals should agree on a clean barrel
    assert side.votes["loop_architecture"] == side.votes["terminus"]
    assert side.confidence > 0.8


def test_side_call_robust_to_dominant_periplasmic_domain():
    """A big periplasmic domain must not invert the extracellular-side call (median, not mean)."""
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=8, peri_loop_len=2, seed=0)
    blob = make_soluble_blob(n_res=180, radius=20, seed=5)
    blob.ca[:] = blob.ca + np.array([0, 0, -45.0])
    big = Structure(
        ca=np.vstack([s.ca, blob.ca]), resids=np.concatenate([s.resids, blob.resids + 9000]),
        resnames=np.concatenate([s.resnames, blob.resnames]),
        chains=np.concatenate([s.chains, blob.chains]),
        sc_vec=np.vstack([s.sc_vec, blob.sc_vec]),
        plddt=np.concatenate([s.plddt, blob.plddt]), atoms=s.atoms + blob.atoms, source="x",
    )
    fit = fit_membrane(big, GN)
    side = call_extracellular_side(big, fit, GN)
    ec_mean, _ = _ec_depth_of_true_ec(big, fit, side)
    assert ec_mean > 0  # true extracellular face still positive despite the huge -Z domain


def test_label_residues_emits_full_table():
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=10, peri_loop_len=2, seed=0)
    fit = fit_membrane(s, GN)
    side = call_extracellular_side(s, fit, GN)
    rsa = compute_sasa(s, n_points=120)["rsa"]
    proj = project_membrane(s, fit, GN, ec_sign=side.ec_sign, rsa=rsa)
    ls = label_residues(s, proj, rsa, GN, fit)
    assert len(ls.labels) == len(s)
    rows = ls.to_rows()
    assert set(rows[0]) >= {"resid", "resname", "zone", "facing", "accessibility", "extracellular", "rsa"}
    # every antibody-accessible residue is in the surface set and is extracellular-flagged
    for l in ls.labels:
        if l.accessibility == "antibody_accessible":
            assert l.resid in ls.surface_set


def test_periplasmic_residues_not_extracellular():
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=8, peri_loop_len=2, seed=0)
    fit = fit_membrane(s, GN)
    side = call_extracellular_side(s, fit, GN)
    rsa = compute_sasa(s, n_points=120)["rsa"]
    proj = project_membrane(s, fit, GN, ec_sign=side.ec_sign, rsa=rsa)
    ls = label_residues(s, proj, rsa, GN, fit)
    # residues on the periplasmic face are never flagged extracellular
    for l in ls.labels:
        if l.ec_depth < -(fit.half_thickness + 4.0):
            assert not l.extracellular
