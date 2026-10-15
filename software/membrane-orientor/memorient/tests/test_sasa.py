"""SASA tests: geometric invariants + agreement with BioPython's Shrake-Rupley (if online)."""

from __future__ import annotations

import numpy as np
import pytest

from memorient.geometry import structure_from_string
from memorient.sasa import _fibonacci_sphere, atom_sasa, compute_sasa

from synthetic import make_soluble_blob


def test_fibonacci_sphere_is_unit_and_uniform():
    pts = _fibonacci_sphere(500)
    assert np.allclose(np.linalg.norm(pts, axis=1), 1.0, atol=1e-9)
    # roughly zero-centred
    assert np.allclose(pts.mean(axis=0), 0.0, atol=0.05)


def test_isolated_atom_sasa_equals_full_sphere():
    # one atom, radius 1.7, probe 1.4 -> full inflated sphere area 4*pi*(3.1)^2
    r = 1.70
    probe = 1.40
    sa = atom_sasa(np.zeros((1, 3)), np.array([r]), n_points=1000, probe=probe)[0]
    expected = 4.0 * np.pi * (r + probe) ** 2
    assert sa == pytest.approx(expected, rel=1e-6)


def test_two_far_atoms_are_both_fully_exposed():
    coords = np.array([[0.0, 0, 0], [50.0, 0, 0]])
    radii = np.array([1.7, 1.7])
    sa = atom_sasa(coords, radii, n_points=500)
    full = 4.0 * np.pi * (1.7 + 1.4) ** 2
    assert np.allclose(sa, full, rtol=1e-6)


def test_buried_atom_has_less_sasa_than_isolated():
    # central atom surrounded by 6 neighbours on axes
    d = 2.0
    coords = np.array([[0, 0, 0], [d, 0, 0], [-d, 0, 0], [0, d, 0], [0, -d, 0], [0, 0, d], [0, 0, -d]], dtype=float)
    radii = np.full(7, 1.7)
    sa = atom_sasa(coords, radii, n_points=1000)
    full = 4.0 * np.pi * (1.7 + 1.4) ** 2
    assert sa[0] < 0.5 * full  # heavily occluded
    assert all(sa[1:] > 0)


def test_rsa_in_unit_range_on_synthetic():
    s = make_soluble_blob(seed=1)
    res = compute_sasa(s, n_points=120)
    assert res["rsa"].min() >= 0.0
    assert res["rsa"].max() <= 1.5  # RSA can slightly exceed 1 for extended residues
    assert len(res["rsa"]) == len(s)


def test_sasa_invariant_under_rotation():
    """Total SASA is frame-invariant; per-residue converges as sampling density rises.

    The Fibonacci test points are fixed in the lab frame, so at finite n_points a rotated
    molecule is sampled slightly differently. Total SASA is robust to this; per-residue has
    a few A^2 of sampling noise that shrinks with more points.
    """
    from synthetic import random_rotation
    s = make_soluble_blob(seed=2)
    a = compute_sasa(s, n_points=960)
    s2 = s.transformed(random_rotation(3), t=np.array([10.0, -4, 7]))
    b = compute_sasa(s2, n_points=960)
    # total SASA nearly identical
    assert a["sasa"].sum() == pytest.approx(b["sasa"].sum(), rel=0.01)
    # per-residue: mean absolute difference is small relative to typical residue SASA
    mad = np.mean(np.abs(np.sort(a["sasa"]) - np.sort(b["sasa"])))
    assert mad < 3.0


@pytest.mark.network
def test_agrees_with_biopython_on_crambin():
    import io
    import urllib.request

    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley

    try:
        txt = urllib.request.urlopen("https://files.rcsb.org/download/1CRN.pdb", timeout=30).read().decode()
    except Exception as e:  # pragma: no cover - network guard
        pytest.skip(f"network unavailable: {e}")

    s = structure_from_string(txt, fmt="pdb")
    mine = compute_sasa(s, n_points=240)["sasa"].sum()

    model = list(PDBParser(QUIET=True).get_structure("s", io.StringIO(txt)).get_models())[0]
    ShrakeRupley(probe_radius=1.40, n_points=240).compute(model, level="R")
    bp = sum(
        r.sasa for ch in model for r in ch
        if r.id[0].strip() == "" and any(a.get_name() == "CA" for a in r)
    )
    assert mine == pytest.approx(bp, rel=0.03)
