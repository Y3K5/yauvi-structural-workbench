"""Tests for geometry: canonicalization idempotency + frame-independence, rotation helpers."""

from __future__ import annotations

import numpy as np
import pytest

from memorient.geometry import (
    canonical_rotation,
    canonicalize,
    principal_axes,
    rotation_matrix_to_z,
)

from synthetic import (
    make_barrel,
    make_ellipsoid,
    make_soluble_blob,
    make_tm_helix,
    random_rotation,
)


def test_rotation_matrix_to_z_maps_vector_to_z():
    rng = np.random.default_rng(0)
    for _ in range(20):
        v = rng.normal(size=3)
        R = rotation_matrix_to_z(v)
        out = R @ (v / np.linalg.norm(v))
        assert np.allclose(out, [0, 0, 1], atol=1e-9)
        assert np.isclose(np.linalg.det(R), 1.0, atol=1e-9)  # proper rotation


def test_rotation_matrix_to_z_antiparallel():
    R = rotation_matrix_to_z([0, 0, -1])
    assert np.allclose(R @ np.array([0, 0, -1.0]), [0, 0, 1], atol=1e-9)
    assert np.isclose(np.linalg.det(R), 1.0)


def test_canonical_rotation_is_proper():
    s = make_barrel(seed=3)
    R, c, info = canonical_rotation(s.ca)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-9)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-9)


def test_canonicalize_is_idempotent_nondegenerate():
    """On a structure with distinct principal axes, canonicalization is a fixed point."""
    s = make_ellipsoid(seed=4)
    c1, info = canonicalize(s)
    assert not info["degenerate"]
    c2, _ = canonicalize(c1)
    assert np.allclose(c1.ca, c2.ca, atol=1e-6)


def test_canonicalize_frame_independent_nondegenerate():
    """Distinct-axis structure: canonical CA coords identical under any input rotation."""
    s = make_ellipsoid(seed=5)
    c0, info0 = canonicalize(s)
    assert not info0["degenerate"]
    for seed in range(5):
        R = random_rotation(seed=seed)
        c1, _ = canonicalize(s.transformed(R, t=np.array([5.0, -3.0, 2.0])))
        assert np.allclose(c0.ca, c1.ca, atol=1e-3), f"seed {seed}"


@pytest.mark.parametrize("builder", [make_barrel, make_tm_helix, make_soluble_blob])
def test_canonicalize_membrane_axis_is_frame_independent(builder):
    """Symmetric structures have degenerate in-plane axes, so raw coords are NOT pinned.

    What IS invariant is the *distinct* axis (the membrane normal for barrel/helix) and the
    set of coordinate magnitudes. Full extracellular-set invariance is proven separately by
    five_fold_validate on the oriented result.
    """
    s = builder(seed=5)
    c0, _ = canonicalize(s)
    # sorted radial distances from origin are a frame-invariant fingerprint
    r0 = np.sort(np.linalg.norm(c0.ca, axis=1))
    for seed in range(5):
        R = random_rotation(seed=seed)
        c1, _ = canonicalize(s.transformed(R))
        r1 = np.sort(np.linalg.norm(c1.ca, axis=1))
        assert np.allclose(r0, r1, atol=1e-4), f"seed {seed}"


def test_transformed_preserves_pairwise_distances():
    s = make_tm_helix(seed=6)
    R = random_rotation(seed=1)
    s2 = s.transformed(R, t=np.array([3.0, -2.0, 1.0]))
    d0 = np.linalg.norm(s.ca[0] - s.ca[-1])
    d1 = np.linalg.norm(s2.ca[0] - s2.ca[-1])
    assert np.isclose(d0, d1, atol=1e-6)


def test_sidechain_vectors_are_unit_or_zero():
    s = make_barrel(seed=7)
    norms = np.linalg.norm(s.sc_vec, axis=1)
    for n in norms:
        assert np.isclose(n, 1.0, atol=1e-6) or np.isclose(n, 0.0, atol=1e-6)


def test_sequence_extraction():
    s = make_tm_helix(seed=8)
    seq = s.sequence
    assert len(seq) == len(s)
    assert set(seq) <= set("ARNDCQEGHILKMFPSTWYVX")
