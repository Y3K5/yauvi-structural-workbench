"""Rigid-transform contracts for the membrane-fit coordinate search.

These tests establish numerical repeatability only.  They do not compare a fitted plane with
an experimental membrane orientation and therefore are not evidence of biological accuracy.
"""

from __future__ import annotations

import numpy as np
import pytest

from memorient.barrel import fit_membrane, fit_membrane_on_normal
from memorient.contexts import get_context
from memorient.geometry import (
    canonical_axis_sign,
    canonical_rotation,
    ordered_intrinsic_rotation,
)

from synthetic import make_barrel, make_tm_helix, random_rotation


GN = get_context("gram_negative_om")
TM = get_context("tm_receptor")


def _axis_angle_deg(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float) / np.linalg.norm(left)
    right = np.asarray(right, dtype=float) / np.linalg.norm(right)
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(left, right)), 0.0, 1.0))))


def test_ordered_intrinsic_frame_is_equivariant_for_pca_degenerate_barrel():
    """The search gauge remains pinned when a barrel's in-plane PCA axes do not."""
    structure = make_barrel(n_strands=12, strand_len=10, seed=17)
    _, _, pca_info = canonical_rotation(structure.ca)
    assert pca_info["degenerate"]

    frame, centroid, _ = ordered_intrinsic_rotation(structure.ca)
    reference = (structure.ca - centroid) @ frame.T
    translation = np.array([19.25, -7.5, 3.125])
    for seed in range(5):
        rotated = structure.transformed(random_rotation(seed), t=translation)
        moved_frame, moved_centroid, _ = ordered_intrinsic_rotation(rotated.ca)
        observed = (rotated.ca - moved_centroid) @ moved_frame.T
        assert np.allclose(observed, reference, rtol=0.0, atol=1e-10), seed


def test_fit_membrane_is_rotation_and_translation_equivariant():
    """A rigid copy has the same fitted physical slab and residue membership."""
    structure = make_barrel(n_strands=12, strand_len=10, seed=0)
    reference = fit_membrane(structure, GN)
    translation = np.array([13.0, -9.0, 4.0])

    for seed in range(5):
        rotation = random_rotation(seed)
        moved = fit_membrane(structure.transformed(rotation, t=translation), GN)
        assert _axis_angle_deg(rotation.T @ moved.normal, reference.normal) < 1e-5, seed
        assert moved.center == pytest.approx(reference.center, abs=1e-10)
        assert moved.half_thickness == pytest.approx(reference.half_thickness, abs=1e-10)
        assert moved.score == pytest.approx(reference.score, abs=1e-12)
        assert moved.components == pytest.approx(reference.components, abs=1e-12)
        assert np.array_equal(moved.embedded_mask, reference.embedded_mask)


def test_fixed_helix_normal_fit_is_rotation_and_translation_equivariant():
    """The fixed-normal path used by mapped alpha helices obeys the same slab contract."""
    structure = make_tm_helix(seed=23)
    normal = np.array([0.0, 0.0, 1.0])
    reference = fit_membrane_on_normal(structure, TM, normal)
    translation = np.array([-4.0, 8.5, 12.0])

    for seed in range(5):
        rotation = random_rotation(seed)
        moved = fit_membrane_on_normal(
            structure.transformed(rotation, t=translation), TM, rotation @ normal,
        )
        assert _axis_angle_deg(rotation.T @ moved.normal, reference.normal) < 1e-5, seed
        assert moved.center == pytest.approx(reference.center, abs=1e-10)
        assert moved.half_thickness == pytest.approx(reference.half_thickness, abs=1e-10)
        assert moved.score == pytest.approx(reference.score, abs=1e-12)
        assert np.array_equal(moved.embedded_mask, reference.embedded_mask)


def test_axis_sign_fallback_uses_intrinsic_order_for_symmetric_skew():
    """Zero-skew sign selection rotates with the points instead of the lab axes."""
    centred = np.array([
        [2.0, 0.0, 0.0],
        [-2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ])
    axis = np.array([1.0, 0.0, 0.0])
    signed = canonical_axis_sign(-axis, centred)
    assert np.allclose(signed, axis)

    for seed in range(5):
        rotation = random_rotation(seed)
        moved_signed = canonical_axis_sign(-(rotation @ axis), centred @ rotation.T)
        assert np.allclose(moved_signed, rotation @ signed, rtol=0.0, atol=1e-12), seed


@pytest.mark.parametrize(
    "coords, message",
    [
        (np.empty((0, 3)), "at least two"),
        (np.array([[0.0, 0.0, 0.0]]), "at least two"),
        (np.zeros((6, 3)), "coincident"),
        (np.column_stack([np.zeros(6), np.zeros(6), np.arange(6.0)]), "collinear"),
        (np.array([[0.0, 0.0, 0.0], [np.nan, 1.0, 2.0]]), "finite"),
    ],
)
def test_intrinsic_frame_rejects_geometrically_undefined_inputs(coords, message):
    """Undefined gauges fail explicitly instead of borrowing a lab-frame direction."""
    with pytest.raises(ValueError, match=message):
        ordered_intrinsic_rotation(coords)


def test_fit_rejects_collinear_structure_instead_of_using_lab_frame():
    structure = make_tm_helix(seed=31)
    structure.ca = np.column_stack([
        np.zeros(len(structure)), np.zeros(len(structure)), np.arange(len(structure), dtype=float),
    ])
    with pytest.raises(ValueError, match="collinear"):
        fit_membrane(structure, GN)


def test_fit_rejects_invalid_scan_and_fixed_normal():
    barrel = make_barrel(seed=37)
    with pytest.raises(ValueError, match="positive integer"):
        fit_membrane(barrel, GN, n_scan=0)

    helix = make_tm_helix(seed=41)
    for invalid in (np.zeros(3), np.array([np.nan, 0.0, 1.0])):
        with pytest.raises(ValueError, match="finite and non-zero"):
            fit_membrane_on_normal(helix, TM, invalid)
