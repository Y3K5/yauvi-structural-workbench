"""Tests for the membrane-context registry (stdlib-only, no compute deps needed)."""

from __future__ import annotations

import dataclasses

import pytest

from memorient.contexts import (
    MembraneContext,
    MembraneModel,
    Metric,
    OrientationMethod,
    ThicknessPrior,
    REGISTRY,
    default_context,
    get_context,
    list_contexts,
)

EXPECTED = {
    "gram_negative_om",
    "eukaryotic_pm",
    "tm_receptor",
    "gram_positive_surface",
    "soluble_secreted",
}


def test_registry_has_expected_contexts():
    assert set(REGISTRY) == EXPECTED
    assert {c.name for c in list_contexts()} == EXPECTED


def test_every_context_validates_and_is_serializable():
    import json

    for c in list_contexts():
        d = c.to_dict()
        # round-trips through JSON
        assert json.loads(json.dumps(d))["name"] == c.name
        # controlled vocabularies respected
        assert c.membrane_model in MembraneModel.ALL
        assert c.orientation_method in OrientationMethod.ALL
        assert all(m in Metric.ALL for m in c.metrics)


def test_rotation_invariance_always_active():
    # The one metric that is meaningful in every context.
    for c in list_contexts():
        assert c.is_metric_active(Metric.ROTATION_INVARIANCE)


def test_metrics_are_not_a_blanket():
    """The whole point: metrics differ by physics."""
    gn = get_context("gram_negative_om")
    euk = get_context("eukaryotic_pm")
    sol = get_context("soluble_secreted")
    tm = get_context("tm_receptor")

    # LPS shielding is gram-negative-OM-only.
    assert gn.is_metric_active(Metric.LPS_SHIELDING)
    assert not euk.is_metric_active(Metric.LPS_SHIELDING)
    assert not sol.is_metric_active(Metric.LPS_SHIELDING)

    # positive-inside orients alpha-helical membranes, not barrels or soluble.
    assert euk.is_metric_active(Metric.POSITIVE_INSIDE)
    assert not gn.is_metric_active(Metric.POSITIVE_INSIDE)

    # a single-pass TM receptor has no lumen -> no lipid/pore gap.
    assert not tm.is_metric_active(Metric.LIPID_PORE_GAP)
    assert gn.is_metric_active(Metric.LIPID_PORE_GAP)

    # soluble runs only rotation-invariance.
    assert sol.metrics == (Metric.ROTATION_INVARIANCE,)


def test_membrane_model_and_sides_consistency():
    gn = get_context("gram_negative_om")
    assert gn.has_bilayer and gn.is_asymmetric and gn.lps_shielding

    euk = get_context("eukaryotic_pm")
    assert euk.has_bilayer and not euk.is_asymmetric and not euk.lps_shielding

    for name in ("gram_positive_surface", "soluble_secreted"):
        c = get_context(name)
        assert not c.has_bilayer
        assert c.thickness_prior is None
        assert not c.has_membrane_sides


def test_thickness_prior_penalty_is_zero_at_mean():
    tp = ThicknessPrior(mean=13.0, sd=2.0)
    assert tp.penalty(13.0) == pytest.approx(0.0)
    assert tp.penalty(15.0) == pytest.approx(1.0)  # one sd out
    assert tp.penalty(11.0) == pytest.approx(1.0)


def test_get_context_unknown_raises_helpful():
    with pytest.raises(KeyError) as ei:
        get_context("mitochondrial_cristae")
    assert "unknown membrane context" in str(ei.value)


def test_default_context_is_gram_negative():
    assert default_context().name == "gram_negative_om"


def test_context_construction_rejects_bad_values():
    # rotation-invariance omitted
    with pytest.raises(ValueError):
        MembraneContext(
            name="bad", description="", membrane_model=MembraneModel.NONE,
            orientation_method=OrientationMethod.SASA_ONLY, thickness_prior=None,
            metrics=(Metric.AROMATIC_GIRDLE,),
        )
    # thickness prior on a non-membrane context
    with pytest.raises(ValueError):
        MembraneContext(
            name="bad", description="", membrane_model=MembraneModel.NONE,
            orientation_method=OrientationMethod.SASA_ONLY,
            thickness_prior=ThicknessPrior(10.0, 1.0),
            metrics=(Metric.ROTATION_INVARIANCE,),
        )
    # LPS shielding without asymmetric model
    with pytest.raises(ValueError):
        MembraneContext(
            name="bad", description="", membrane_model=MembraneModel.SYMMETRIC_PHOSPHOLIPID,
            orientation_method=OrientationMethod.TM_HELIX_BELT,
            thickness_prior=ThicknessPrior(15.0, 2.0),
            metrics=(Metric.ROTATION_INVARIANCE,), lps_shielding=True,
        )


def test_contexts_are_frozen():
    c = default_context()
    with pytest.raises(dataclasses.FrozenInstanceError):
        c.name = "nope"  # type: ignore[misc]
