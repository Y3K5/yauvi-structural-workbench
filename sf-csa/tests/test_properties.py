#!/usr/bin/env python3
"""Property-based tests for the sf-csa interpretation rules.

`test_classification.py` checks that the default rules behave as the hardcoded
version did, example by example. These tests check the claims that must hold for
every hit, which is what "closed vocabulary" and "bounded interpretation"
actually assert:

  P1  `classify_hit` returns a 4-tuple whose interpretation is always in
      `CLASSIFICATION_VOCABULARY`, for every hit and every configuration
  P2  a hit below the structural threshold is `unresolved_or_conflicted`,
      whatever else the hit says
  P3  a contested group is never promoted, at any TM score or coverage
  P4  a missing, blank or unparseable similarity score reads as below threshold
      — never as a pass
  P5  `classify_title` is total and deterministic, and its family ordering is a
      contract rather than an accident
  P6  the structural and sequence legs never merge into one score

The generated space deliberately includes hostile input — empty strings, `None`,
non-numeric scores, out-of-range coverages — because the real inputs are parsed
out of a Foldseek TSV, and a field that arrives blank must not become a pass.
"""
from __future__ import annotations

import inspect

import pytest

from hypothesis import assume, given, settings, strategies as st

from sf_csa.core import (
    CLASSIFICATION_VOCABULARY,
    DEFAULT_CONTESTED_GROUPS,
    DEFAULT_DIVERGENCE_SETS,
    DEFAULT_MECHANISM_FAMILIES,
    DEFAULT_TITLE_TRAPS,
    classify_hit,
    classify_title,
    norm_id,
    structural_category,
)

settings.register_profile("sfcsa", settings(max_examples=200, deadline=None))
settings.load_profile("sfcsa")

# The claims each vocabulary member makes about shared function, strongest
# first. Used only to express "never promoted above", not as a total order.
CLAIM_STRENGTH = {
    "exact_function_supported": 5,
    "probable_same_function": 4,
    "same_mechanism_class": 3,
    "candidate_functional_divergence": 2,
    "structural_analogy_only": 1,
    "unresolved_or_conflicted": 0,
}

# The two labels that assert shared function. The contested-group rule and the
# title traps both exist to keep specific inputs out of this set.
FUNCTION_CLAIMS = frozenset({"exact_function_supported", "probable_same_function"})

KNOWN_GROUPS = sorted(
    {family["group"] for family in DEFAULT_MECHANISM_FAMILIES}
    | {
        narrower["group"]
        for family in DEFAULT_MECHANISM_FAMILIES
        for narrower in family.get("refine", ())
    }
    | {"unknown"}
)

CONTESTED_GROUP_NAMES = frozenset(entry["group"] for entry in DEFAULT_CONTESTED_GROUPS)

# Thresholds used by the shipped campaign config, as the defaults to vary around.
WHOLE_COVERAGE = 0.8
SAME_FOLD_TM = 0.5


# -- strategies -----------------------------------------------------------

scores = st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False)

# What a Foldseek TSV field can actually contain by the time it reaches here:
# a number, a numeric string, or nothing at all.
hostile_scores = st.one_of(
    scores,
    st.just(""),
    st.just(None),
    st.just("0.9"),
    st.just("nan"),
    st.integers(min_value=0, max_value=1),
)

titles = st.one_of(
    st.text(max_size=60),
    st.sampled_from(
        [
            "BamA outer membrane protein assembly factor",
            "RagA TonB-dependent receptor",
            "TonB-dependent heme receptor",
            "PorG type IX secretion protein",
            "major surface protein MSP",
            "outer-membrane porin barrel",
            "toluene monooxygenase",
            "hypothetical protein",
            "",
        ]
    ),
)


@st.composite
def hits(draw):
    return {
        "target": draw(st.text(alphabet="ABCDEFGH0123456789|._", max_size=14)),
        "theader": draw(titles),
        "alntmscore": draw(hostile_scores),
        "qcov": draw(hostile_scores),
        "tcov": draw(hostile_scores),
        "evalue": draw(st.one_of(st.just(""), st.floats(0.0, 10.0))),
    }


@st.composite
def queries(draw):
    return {
        "accession": draw(st.text(alphabet="ABCDEFGH0123456789", min_size=1, max_size=8)),
        "mechanism_group": draw(st.sampled_from(KNOWN_GROUPS)),
    }


categories = st.sampled_from(
    ["whole_architecture_match", "domain_or_partial_match", "below_structural_similarity_threshold"]
)


@st.composite
def target_metas(draw):
    return draw(
        st.one_of(
            st.none(),
            st.builds(
                lambda group, rbh: {"mechanism_group": group, "rbh": rbh},
                st.sampled_from(KNOWN_GROUPS),
                st.booleans(),
            ),
            st.just({}),
        )
    )


# -- P1: the vocabulary is closed ----------------------------------------


@given(queries(), hits(), categories, target_metas())
def test_every_interpretation_is_in_the_vocabulary(query, hit, category, target_meta):
    interpretation, reason, evidence, boundary = classify_hit(
        query, hit, category, target_meta
    )
    assert interpretation in CLASSIFICATION_VOCABULARY


@given(queries(), hits(), categories, target_metas())
def test_classify_hit_is_total(query, hit, category, target_meta):
    """No hit in the generated space raises. A crash is not a classification."""
    result = classify_hit(query, hit, category, target_meta)
    assert isinstance(result, tuple) and len(result) == 4


@given(queries(), hits(), categories, target_metas())
def test_every_interpretation_carries_a_reason_and_a_boundary(
    query, hit, category, target_meta
):
    """The boundary is the scientific content: what the label does *not* license."""
    interpretation, reason, evidence, boundary = classify_hit(
        query, hit, category, target_meta
    )
    assert reason.strip()
    assert evidence.strip()
    assert boundary.strip()


@given(queries(), hits(), categories, target_metas())
def test_no_free_text_from_the_hit_reaches_the_interpretation(
    query, hit, category, target_meta
):
    """A PDB title is data. It is never promoted into a label."""
    interpretation, _, _, _ = classify_hit(query, hit, category, target_meta)
    assert interpretation in CLASSIFICATION_VOCABULARY
    assert interpretation != hit["theader"]


# -- P2: below threshold is decisive -------------------------------------


@given(queries(), hits(), target_metas())
def test_below_threshold_is_always_unresolved(query, hit, target_meta):
    """No amount of matching annotation can rescue a failed structural test."""
    interpretation, _, _, _ = classify_hit(
        query, hit, "below_structural_similarity_threshold", target_meta
    )
    # The one exception is an exact self-match, which is a control rather than
    # an independent observation, and is checked separately below.
    if norm_id(hit.get("target", "")) != query["accession"]:
        assert interpretation == "unresolved_or_conflicted"


@given(queries(), hits())
def test_below_threshold_never_claims_shared_function_for_a_different_accession(query, hit):
    assume(norm_id(hit.get("target", "")) != query["accession"])
    interpretation, _, _, _ = classify_hit(
        query, hit, "below_structural_similarity_threshold", {"mechanism_group": query["mechanism_group"], "rbh": True}
    )
    assert interpretation not in FUNCTION_CLAIMS


@given(scores, scores, scores)
def test_structural_category_respects_its_own_thresholds(tm, qcov, tcov):
    """The category boundaries are exactly what the thresholds say."""
    category = structural_category(
        {"alntmscore": tm, "qcov": qcov, "tcov": tcov}, WHOLE_COVERAGE, SAME_FOLD_TM
    )
    if tm < SAME_FOLD_TM:
        assert category == "below_structural_similarity_threshold"
    elif qcov >= WHOLE_COVERAGE and tcov >= WHOLE_COVERAGE:
        assert category == "whole_architecture_match"
    else:
        assert category == "domain_or_partial_match"


@given(scores, scores, scores)
def test_raising_the_threshold_never_strengthens_the_category(tm, qcov, tcov):
    """Monotonicity in the threshold: a stricter run cannot claim more."""
    rank = {
        "below_structural_similarity_threshold": 0,
        "domain_or_partial_match": 1,
        "whole_architecture_match": 2,
    }
    hit = {"alntmscore": tm, "qcov": qcov, "tcov": tcov}
    lenient = structural_category(hit, WHOLE_COVERAGE, SAME_FOLD_TM)
    strict = structural_category(hit, min(1.0, WHOLE_COVERAGE + 0.15), min(1.0, SAME_FOLD_TM + 0.2))
    assert rank[strict] <= rank[lenient]


# -- P3: a contested group is never promoted -----------------------------


@given(hits(), categories, st.booleans(), st.sampled_from(sorted(CONTESTED_GROUP_NAMES)))
def test_a_contested_group_is_never_promoted(hit, category, rbh, group):
    """Whatever the structural evidence, a contested group stays unresolved.

    This is the rule most at risk from a well-meaning refactor: the match looks
    excellent by every numeric measure, and the reason it may not be promoted is
    biological rather than computational.
    """
    query = {"accession": "Q_CONTESTED", "mechanism_group": group}
    target_meta = {"mechanism_group": group, "rbh": rbh}
    assume(norm_id(hit.get("target", "")) != query["accession"])
    interpretation, _, _, _ = classify_hit(query, hit, category, target_meta)
    assert interpretation not in FUNCTION_CLAIMS


@given(scores, scores, scores, st.booleans())
def test_a_contested_group_is_unresolved_at_every_score(tm, qcov, tcov, rbh):
    group = sorted(CONTESTED_GROUP_NAMES)[0]
    query = {"accession": "Q_CONTESTED", "mechanism_group": group}
    hit = {"target": "T_OTHER", "theader": "", "alntmscore": tm, "qcov": qcov, "tcov": tcov}
    category = structural_category(hit, WHOLE_COVERAGE, SAME_FOLD_TM)
    interpretation, _, _, _ = classify_hit(
        query, hit, category, {"mechanism_group": group, "rbh": rbh}
    )
    assert interpretation == "unresolved_or_conflicted"


@given(hits(), categories)
def test_an_empty_contested_list_is_not_treated_as_the_default(hit, category):
    """Passing `[]` means "nothing is contested", not "use the defaults".

    The distinction matters: `contested = DEFAULT if contested is None else
    contested` is correct, and `contested or DEFAULT` would silently restore the
    campaign's biology into a run that explicitly disabled it.
    """
    group = sorted(CONTESTED_GROUP_NAMES)[0]
    query = {"accession": "Q", "mechanism_group": group}
    assume(norm_id(hit.get("target", "")) != "Q")
    assume(category != "below_structural_similarity_threshold")
    with_default = classify_hit(query, hit, category, {"mechanism_group": group, "rbh": True})
    with_none = classify_hit(
        query, hit, category, {"mechanism_group": group, "rbh": True}, contested=[]
    )
    assert with_default[0] == "unresolved_or_conflicted"
    assert with_none[0] != "unresolved_or_conflicted"


# -- P4: a missing score is not a pass -----------------------------------


@given(st.sampled_from(["", None, "  ", "\t", "NA", "-"]), hostile_scores, hostile_scores)
def test_a_blank_similarity_score_never_becomes_a_pass(tm, qcov, tcov):
    """Fail closed on a field that never arrived — by either mechanism.

    The invariant that matters is that a missing score cannot yield a structural
    claim. This module reaches that two different ways depending on the exact
    bytes in the field (see the next test), so the property is stated as "below
    threshold or a hard error", which is what fail-closed actually requires.
    """
    try:
        category = structural_category(
            {"alntmscore": tm, "qcov": qcov, "tcov": tcov}, WHOLE_COVERAGE, SAME_FOLD_TM
        )
    except ValueError:
        return  # a loud refusal is a fail-closed outcome
    assert category == "below_structural_similarity_threshold"


def test_a_blank_score_and_a_whitespace_score_fail_closed_differently():
    """Documenting an inconsistency, rather than papering over it.

    `float(hit.get("alntmscore") or 0)` makes `""` and `None` falsy, so they
    become `0.0` — a silent below-threshold hit that is counted as a real
    comparison in the release manifest. A whitespace-only field is truthy, so it
    reaches `float()` and raises, aborting the run.

    Both are fail-closed in the sense that neither becomes a pass, but the
    difference is decided by whether an upstream tool emitted `""` or `" "` for
    the same missing value. One outcome is a row in the output, the other is no
    output at all. See REVIEW.md §7.
    """
    empty = structural_category({"alntmscore": ""}, WHOLE_COVERAGE, SAME_FOLD_TM)
    assert empty == "below_structural_similarity_threshold"

    with pytest.raises(ValueError):
        structural_category({"alntmscore": " "}, WHOLE_COVERAGE, SAME_FOLD_TM)


@given(hostile_scores, hostile_scores)
def test_a_missing_score_key_is_below_threshold(qcov, tcov):
    """A hit dict with no `alntmscore` key at all must not pass."""
    category = structural_category({"qcov": qcov, "tcov": tcov}, WHOLE_COVERAGE, SAME_FOLD_TM)
    assert category == "below_structural_similarity_threshold"


@given(scores, st.floats(min_value=WHOLE_COVERAGE, max_value=1.0))
def test_missing_coverage_cannot_reach_whole_architecture(tm, present):
    """Blank coverage may not be read as full coverage — on either side.

    Each coverage field is blanked independently while the other is left at a
    passing value, so a default that reads one missing field as full coverage is
    caught even though the other field still fails. Mutation testing (mutant S5)
    showed an earlier version blanked both at once and so missed exactly that.
    """
    for qcov, tcov in (("", present), (present, ""), ("", ""), (None, present), (present, None)):
        category = structural_category(
            {"alntmscore": tm, "qcov": qcov, "tcov": tcov}, WHOLE_COVERAGE, SAME_FOLD_TM
        )
        assert category != "whole_architecture_match"


@given(queries(), titles)
def test_a_hit_with_no_scores_at_all_never_claims_shared_function(query, title):
    hit = {"target": "T_UNRELATED", "theader": title}
    category = structural_category(hit, WHOLE_COVERAGE, SAME_FOLD_TM)
    assume(norm_id(hit["target"]) != query["accession"])
    interpretation, _, _, _ = classify_hit(query, hit, category)
    assert interpretation == "unresolved_or_conflicted"


def test_a_non_numeric_score_is_a_hard_failure_not_a_silent_pass():
    """`float('abc')` raises rather than defaulting — the fail-closed choice.

    Pinned because a `try/except: return 0` "fix" here would be a regression: a
    corrupted TSV field must stop the run, not quietly become a below-threshold
    hit that gets counted as a real comparison.
    """
    with pytest.raises(ValueError):
        structural_category({"alntmscore": "not-a-number"}, WHOLE_COVERAGE, SAME_FOLD_TM)


# -- P5: classify_title is total, deterministic, and order-dependent ------


@given(titles)
def test_classify_title_is_total(title):
    assert classify_title(title) in KNOWN_GROUPS


@given(titles)
def test_classify_title_is_deterministic(title):
    assert classify_title(title) == classify_title(title)


@given(st.text(max_size=200))
def test_arbitrary_text_never_invents_a_group(title):
    assert classify_title(title) in KNOWN_GROUPS


@given(titles)
def test_classify_title_is_case_insensitive(title):
    assert classify_title(title.upper()) == classify_title(title.lower())


@given(titles)
def test_an_empty_family_list_yields_unknown(title):
    """`families=[]` means no families, not "fall back to the defaults"."""
    assert classify_title(title, families=[]) == "unknown"


def test_family_order_is_a_contract_not_an_accident():
    """A title matching two families resolves to whichever is declared first.

    `DEFAULT_MECHANISM_FAMILIES` is scanned in order and the first match wins, so
    the list order is load-bearing: "TonB-dependent receptor porin" matches both
    `susc_raga_importer` and the catch-all `generic_om_barrel`, and only the
    declaration order decides. Reordering the list silently reclassifies hits, so
    the ordering is asserted here rather than left implicit.
    """
    title = "TonB-dependent receptor porin barrel"
    assert classify_title(title) == "susc_raga_importer"

    reordered = [
        family
        for name in ("generic_om_barrel", "susc_raga_importer")
        for family in DEFAULT_MECHANISM_FAMILIES
        if family["group"] == name
    ]
    assert classify_title(title, families=reordered) == "generic_om_barrel"


def test_a_refinement_wins_over_its_parent_family():
    """The heme refinement must beat the generic importer, or the narrower
    mechanism claim is lost."""
    assert classify_title("TonB-dependent heme receptor") == "tonb_heme_receptor"
    assert classify_title("RagA TonB-dependent receptor") == "susc_raga_importer"


# -- P6: the two legs never merge ----------------------------------------


def test_no_signature_in_the_interpretation_path_combines_the_two_legs():
    """Structural and sequence evidence are reported side by side, never summed.

    Asserted structurally: `classify_hit` and `structural_category` receive
    structural evidence and manifest configuration, and nothing carrying a
    sequence-identity or bitscore field. A parameter named for sequence
    similarity appearing here would mean the legs had been merged.
    """
    forbidden = ("bitscore", "pident", "fident", "sequence_identity", "combined", "score_sum")
    for function in (classify_hit, structural_category):
        parameters = set(inspect.signature(function).parameters)
        assert not parameters & set(forbidden)


@given(queries(), hits(), categories, target_metas())
def test_the_interpretation_never_returns_a_numeric_score(query, hit, category, target_meta):
    """The output is a claim with a boundary, not a number to be ranked."""
    for field in classify_hit(query, hit, category, target_meta):
        assert isinstance(field, str)
        assert not field.replace(".", "", 1).isdigit()


@given(queries(), hits(), categories)
def test_sequence_evidence_only_enters_through_the_rbh_flag(query, hit, category):
    """The only sequence input to interpretation is a boolean, by design.

    Reciprocal-best-hit status is a yes/no fact about the sequence leg. It can
    promote a match to `probable_same_function`, but no continuous sequence
    similarity is blended into the structural claim — so flipping the flag can
    change the label by at most one step in the vocabulary.
    """
    assume(norm_id(hit.get("target", "")) != query["accession"])
    group = query["mechanism_group"]
    without = classify_hit(query, hit, category, {"mechanism_group": group, "rbh": False})[0]
    with_rbh = classify_hit(query, hit, category, {"mechanism_group": group, "rbh": True})[0]
    assert CLAIM_STRENGTH[with_rbh] - CLAIM_STRENGTH[without] <= 1


@given(queries(), hits(), categories)
def test_rbh_alone_never_promotes_a_partial_match_to_shared_function(query, hit, category):
    """A reciprocal best hit on a partial structural match is not shared function."""
    assume(category == "domain_or_partial_match")
    assume(norm_id(hit.get("target", "")) != query["accession"])
    group = query["mechanism_group"]
    assume(group not in CONTESTED_GROUP_NAMES)
    interpretation, _, _, _ = classify_hit(
        query, hit, category, {"mechanism_group": group, "rbh": True}
    )
    assert interpretation != "probable_same_function"


# -- the divergence rule -------------------------------------------------


@given(hits(), categories, st.booleans())
def test_a_cross_group_match_inside_a_divergence_set_is_divergence(hit, category, rbh):
    """Shared framework, different substrate: divergence, not analogy."""
    entry = DEFAULT_DIVERGENCE_SETS[0]
    qgroup, tgroup = entry["groups"][0], entry["groups"][1]
    query = {"accession": "Q_DIVERGE", "mechanism_group": qgroup}
    assume(norm_id(hit.get("target", "")) != query["accession"])
    assume(category != "below_structural_similarity_threshold")
    interpretation, _, _, _ = classify_hit(
        query, hit, category, {"mechanism_group": tgroup, "rbh": rbh}
    )
    assert interpretation == "candidate_functional_divergence"
    assert interpretation not in FUNCTION_CLAIMS


@given(hits(), categories)
def test_divergence_is_never_a_shared_function_claim(hit, category):
    entry = DEFAULT_DIVERGENCE_SETS[0]
    for qgroup in entry["groups"]:
        for tgroup in entry["groups"]:
            if qgroup == tgroup:
                continue
            query = {"accession": "Q_D", "mechanism_group": qgroup}
            assume(norm_id(hit.get("target", "")) != "Q_D")
            interpretation, _, _, _ = classify_hit(
                query, hit, category, {"mechanism_group": tgroup, "rbh": True}
            )
            assert interpretation not in FUNCTION_CLAIMS


# -- the self-match control ----------------------------------------------


@given(hits(), categories)
def test_a_self_match_is_labelled_a_control(hit, category):
    """An accession matching itself is `exact_function_supported`, and its
    boundary must say it is a control rather than evidence."""
    accession = "P00001"
    query = {"accession": accession, "mechanism_group": "omp85_bama"}
    hit = dict(hit, target=accession)
    interpretation, reason, evidence, boundary = classify_hit(query, hit, category)
    assert interpretation == "exact_function_supported"
    assert "control" in boundary.lower()


@given(st.text(alphabet="ABCDEFGH0123456789", min_size=1, max_size=8))
def test_a_self_match_survives_id_decoration(accession):
    """`norm_id` strips `.pdb`, `_oriented` and database prefixes, so a
    decorated self-match is still recognised as a control rather than being
    reported as an independent confirmation."""
    for decorated in (f"{accession}.pdb", f"{accession}_oriented", f"sp|{accession}|NAME"):
        assert norm_id(decorated) == accession
        query = {"accession": accession, "mechanism_group": "omp85_bama"}
        interpretation, _, _, boundary = classify_hit(
            query, {"target": decorated, "theader": ""}, "whole_architecture_match"
        )
        assert interpretation == "exact_function_supported"
        assert "control" in boundary.lower()


# -- the title traps -----------------------------------------------------


@given(st.sampled_from([t["substring"] for t in DEFAULT_TITLE_TRAPS]))
def test_a_trapped_substring_has_a_declared_prohibition(substring):
    """Each trap names the labels it forbids, and they are real labels."""
    trap = next(t for t in DEFAULT_TITLE_TRAPS if t["substring"] == substring)
    assert trap["must_not_promote_to"]
    for label in trap["must_not_promote_to"]:
        assert label in CLASSIFICATION_VOCABULARY


def test_the_title_trap_is_a_release_audit_not_a_classifier_guard():
    """Documenting where the trap does and does not act.

    The traps are enforced by `verify_release`, which audits a finished release.
    `classify_hit` itself does not consult them: a hit whose title contains a
    trapped substring can be classified `probable_same_function` in-process, and
    the violation is caught only when the release is verified.

    That is a defensible design — the audit checks the artefact rather than
    trusting the code that wrote it — but it means the trap is not a guard on the
    classifier, and a caller using `classify_hit` directly gets no protection
    from it. Pinned so the distinction stays visible. See REVIEW.md §7.
    """
    trapped = DEFAULT_TITLE_TRAPS[0]["substring"]
    query = {"accession": "Q_TRAP", "mechanism_group": "omp85_bama"}
    hit = {
        "target": "T_TRAP",
        "theader": f"BamA outer membrane protein assembly factor with {trapped}",
        "alntmscore": 0.9,
        "qcov": 0.95,
        "tcov": 0.95,
    }
    interpretation, _, _, _ = classify_hit(
        query, hit, "whole_architecture_match", {"mechanism_group": "omp85_bama", "rbh": True}
    )
    # The classifier promotes it; only the release audit would object.
    assert interpretation == "probable_same_function"
    assert interpretation in set(DEFAULT_TITLE_TRAPS[0]["must_not_promote_to"])


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
