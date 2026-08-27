"""Interpretation is bounded, and the campaign rules are configuration.

The two rules that used to be literals inside `classify_hit` — the contested
group, and the set of groups that share a framework but not a substrate — are now
manifest entries. These tests hold both halves: the defaults still behave exactly
as the hardcoded version did, and a campaign against other organisms can replace
them without editing source.
"""
from __future__ import annotations

import pytest

from sf_csa.core import (
    CLASSIFICATION_VOCABULARY,
    DEFAULT_CONTESTED_GROUPS,
    DEFAULT_DIVERGENCE_SETS,
    DEFAULT_MECHANISM_FAMILIES,
    classify_hit,
    classify_title,
    structural_category,
)


def query(accession="Q1", group="omp85_bama"):
    return {"accession": accession, "mechanism_group": group}


def hit(target="T1", header="BamA outer membrane protein assembly factor"):
    return {"target": target, "theader": header}


# --- the closed vocabulary ------------------------------------------------


def test_the_vocabulary_has_exactly_six_labels():
    assert len(CLASSIFICATION_VOCABULARY) == 6
    assert len(set(CLASSIFICATION_VOCABULARY)) == 6


@pytest.mark.parametrize(
    ("qgroup", "tgroup", "category"),
    [
        ("omp85_bama", "omp85_bama", "whole_architecture_match"),
        ("omp85_bama", "t9ss_porg", "domain_or_partial_match"),
        ("msp_contested", "msp_contested", "whole_architecture_match"),
        ("susc_raga_importer", "generic_om_barrel", "whole_architecture_match"),
        ("unknown", "unknown", "below_structural_similarity_threshold"),
    ],
)
def test_every_outcome_is_inside_the_vocabulary(qgroup, tgroup, category):
    label, *_ = classify_hit(
        query(group=qgroup), hit(), category, {"mechanism_group": tgroup}
    )
    assert label in CLASSIFICATION_VOCABULARY


# --- behaviour preserved from the hardcoded version -----------------------


def test_a_self_match_is_a_control_not_evidence():
    label, basis, _, boundary = classify_hit(
        query(accession="Q1"), hit(target="Q1"), "whole_architecture_match"
    )
    assert label == "exact_function_supported"
    assert "control" in boundary


def test_below_threshold_never_gets_a_fold_interpretation():
    label, *_ = classify_hit(
        query(), hit(), "below_structural_similarity_threshold", {"mechanism_group": "omp85_bama"}
    )
    assert label == "unresolved_or_conflicted"


def test_a_contested_group_is_never_promoted():
    """However good the structural match, a contested group stays unresolved."""
    label, basis, _, _ = classify_hit(
        query(group="msp_contested"),
        hit(target="T9"),
        "whole_architecture_match",
        {"mechanism_group": "msp_contested", "rbh": True},
    )
    assert label == "unresolved_or_conflicted"
    assert "contested" in basis


def test_rbh_plus_whole_architecture_reaches_probable_same_function():
    label, *_ = classify_hit(
        query(), hit(), "whole_architecture_match", {"mechanism_group": "omp85_bama", "rbh": True}
    )
    assert label == "probable_same_function"


def test_same_group_without_rbh_stops_at_mechanism_class():
    label, _, _, boundary = classify_hit(
        query(), hit(), "whole_architecture_match", {"mechanism_group": "omp85_bama"}
    )
    assert label == "same_mechanism_class"
    assert "not transferred" in boundary


def test_the_default_divergence_set_still_applies():
    label, *_ = classify_hit(
        query(group="susc_raga_importer"),
        hit(),
        "whole_architecture_match",
        {"mechanism_group": "tonb_heme_receptor"},
    )
    assert label == "candidate_functional_divergence"


def test_an_unrelated_cross_group_match_is_analogy_only():
    label, *_ = classify_hit(
        query(group="omp85_bama"), hit(), "whole_architecture_match",
        {"mechanism_group": "t9ss_porg"},
    )
    assert label == "structural_analogy_only"


# --- the rules are now configuration --------------------------------------


def test_a_campaign_can_declare_its_own_contested_group():
    label, basis, _, _ = classify_hit(
        query(group="my_family"),
        hit(),
        "whole_architecture_match",
        {"mechanism_group": "my_family"},
        None,
        [{"group": "my_family", "reason": "locally contested"}],
    )
    assert label == "unresolved_or_conflicted"
    assert basis == "locally contested"


def test_supplying_an_empty_contested_list_disables_the_default():
    """A campaign with no contested groups must not inherit periodontal biology."""
    label, *_ = classify_hit(
        query(group="msp_contested"),
        hit(),
        "whole_architecture_match",
        {"mechanism_group": "msp_contested"},
        None,
        [],
    )
    assert label == "same_mechanism_class"


def test_a_campaign_can_declare_its_own_divergence_set():
    label, basis, _, _ = classify_hit(
        query(group="alpha"),
        hit(),
        "whole_architecture_match",
        {"mechanism_group": "beta"},
        None,
        None,
        [{"groups": ["alpha", "beta"], "reason": "shared scaffold, different substrate"}],
    )
    assert label == "candidate_functional_divergence"
    assert basis == "shared scaffold, different substrate"


def test_supplying_an_empty_divergence_list_disables_the_default():
    label, *_ = classify_hit(
        query(group="susc_raga_importer"),
        hit(),
        "whole_architecture_match",
        {"mechanism_group": "tonb_heme_receptor"},
        None,
        None,
        [],
    )
    assert label == "structural_analogy_only"


def test_the_defaults_are_documented_as_overridable():
    for table in (DEFAULT_CONTESTED_GROUPS, DEFAULT_DIVERGENCE_SETS):
        assert table, "a default table is empty"
        for entry in table:
            assert entry.get("reason"), "every default entry must state its reason"


# --- title classification -------------------------------------------------


def test_a_title_is_matched_to_a_mechanism_family():
    assert classify_title("BamA outer membrane protein assembly factor") == "omp85_bama"


def test_refinement_narrows_a_matched_family():
    """A TonB-dependent transporter becomes a heme receptor when the title says so."""
    assert classify_title("TonB-dependent hemoglobin receptor") == "tonb_heme_receptor"
    assert classify_title("TonB-dependent transporter") == "susc_raga_importer"


def test_an_unmatched_title_is_unknown():
    assert classify_title("hypothetical protein of unknown function") == "unknown"


def test_a_campaign_can_supply_its_own_families():
    families = [{"group": "my_family", "pattern": r"widget"}]
    assert classify_title("widget-binding protein", families) == "my_family"
    assert classify_title("BamA assembly factor", families) == "unknown"


def test_default_families_are_regex_patterns_that_compile():
    import re

    for family in DEFAULT_MECHANISM_FAMILIES:
        re.compile(family["pattern"])
        for narrower in family.get("refine", ()):
            re.compile(narrower["pattern"])


# --- structural category --------------------------------------------------


def test_whole_architecture_needs_score_and_coverage_on_both_sides():
    assert structural_category(
        {"alntmscore": "0.8", "qcov": "0.9", "tcov": "0.9"}, 0.7, 0.5
    ) == "whole_architecture_match"


def test_good_score_with_poor_coverage_is_only_partial():
    assert structural_category(
        {"alntmscore": "0.8", "qcov": "0.3", "tcov": "0.9"}, 0.7, 0.5
    ) == "domain_or_partial_match"


def test_a_low_score_is_below_threshold():
    assert structural_category(
        {"alntmscore": "0.2", "qcov": "0.9", "tcov": "0.9"}, 0.7, 0.5
    ) == "below_structural_similarity_threshold"


def test_missing_scores_are_treated_as_zero_not_as_a_match():
    assert structural_category({}, 0.7, 0.5) == "below_structural_similarity_threshold"
