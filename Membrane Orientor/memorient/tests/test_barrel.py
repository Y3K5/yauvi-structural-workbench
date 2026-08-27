"""Membrane-fit + classifier tests on synthetic structures (network-free) and OMP checks."""

from __future__ import annotations

import numpy as np
import pytest

from memorient.barrel import classify_membrane_protein, fit_membrane
from memorient.contexts import get_context
from memorient.geometry import Structure

from synthetic import make_barrel, make_soluble_blob, random_rotation

GN = get_context("gram_negative_om")


def _axis_angle(n, ref):
    n = np.asarray(n, float) / np.linalg.norm(n)
    ref = np.asarray(ref, float) / np.linalg.norm(ref)
    return np.degrees(np.arccos(np.clip(abs(np.dot(n, ref)), 0, 1)))


def test_fit_recovers_synthetic_barrel_normal():
    s = make_barrel(n_strands=12, strand_len=10, seed=0)
    fit = fit_membrane(s, GN)
    assert _axis_angle(fit.normal, [0, 0, 1]) < 12.0
    assert fit.delta_kd > 1.5          # strong lipid/pore gap
    assert 8.0 < fit.half_thickness < 18.0
    assert fit.n_embedded >= 40
    assert fit.inner_frac < 0.25       # hollow


def test_synthetic_barrel_classifies_as_barrel():
    s = make_barrel(seed=1)
    cls = classify_membrane_protein(s, GN)
    assert cls.label == "barrel"
    assert cls.confidence > 0.7


def test_soluble_blob_classifies_as_surface():
    s = make_soluble_blob(n_res=150, radius=16, seed=2)
    cls = classify_membrane_protein(s, GN)
    assert cls.label in ("surface", "uncertain")
    assert not (cls.label == "barrel")


def test_fit_robust_to_dominant_extramembrane_domain():
    """A big soluble domain hung off one end must not capture the normal (the whole point)."""
    s = make_barrel(n_strands=12, strand_len=10, seed=0)
    blob = make_soluble_blob(n_res=220, radius=22, seed=9)
    blob.ca[:] = blob.ca + np.array([0, 0, 60.0])
    big = Structure(
        ca=np.vstack([s.ca, blob.ca]),
        resids=np.concatenate([s.resids, blob.resids + 10000]),
        resnames=np.concatenate([s.resnames, blob.resnames]),
        chains=np.concatenate([s.chains, blob.chains]),
        sc_vec=np.vstack([s.sc_vec, blob.sc_vec]),
        plddt=np.concatenate([s.plddt, blob.plddt]),
        atoms=s.atoms + blob.atoms, source="barrel+blob",
    )
    fit = fit_membrane(big, GN)
    # normal still along Z, not dragged toward the blob offset axis
    assert _axis_angle(fit.normal, [0, 0, 1]) < 15.0
    assert fit.delta_kd > 1.0


@pytest.mark.parametrize("seed", range(5))
def test_fit_normal_tracks_input_rotation(seed):
    s = make_barrel(n_strands=12, strand_len=10, seed=0)
    R = random_rotation(seed)
    fit = fit_membrane(s.transformed(R), GN)
    ref = R @ np.array([0, 0, 1.0])
    assert _axis_angle(fit.normal, ref) < 15.0


@pytest.mark.network
def test_real_omp_classifies_barrel():
    import urllib.request

    from memorient.geometry import structure_from_string

    try:
        txt = urllib.request.urlopen("https://files.rcsb.org/download/1BXW.pdb", timeout=30).read().decode()
    except Exception as e:  # pragma: no cover
        pytest.skip(f"network unavailable: {e}")
    s = structure_from_string(txt, fmt="pdb")
    cls = classify_membrane_protein(s, GN)
    assert cls.label == "barrel"
    assert cls.fit.half_thickness < 18.0
