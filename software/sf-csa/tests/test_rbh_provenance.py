"""`probable_same_function` must rest on computed evidence, not on curator input.

These tests were written red, before the fix, to hold three separate defects that
`SF_CSA_PREADOPTION_FINDINGS.md` and the code review of 2026-09-01 established:

1. **Fail-open.** A curator writing `"rbh": true` into a query record reaches
   `classify_hit` through `target_meta` and emits the strongest non-exact
   functional label in the closed vocabulary, carrying the basis string "RBH plus
   compatible whole architecture and mechanism annotation" — for a pair on which
   no reciprocal best hit was ever computed.
2. **Dead computed path.** The real reciprocal-best-hit runs after classification
   and writes to `h["orthology_status"]` on sequence rows, so the label is never
   reached by evidence.
3. **Not pairwise.** `target_meta` is keyed by target accession alone, so an
   `rbh` flag living on a target row is query-independent by construction: one
   query establishing RBH would flag that target for every other query.

The panel gate `analogy_or_unrelated_promoted_to_function_max = 0` is only a real
gate once all three hold. Until then it can be satisfied, or violated, without
any computed evidence existing.
"""
from __future__ import annotations

import pytest

from sf_csa.core import classify_hit


WHOLE = "whole_architecture_match"


def query(accession="Q1", group="omp85_bama"):
    return {"accession": accession, "mechanism_group": group}


def hit(target="T1", header="BamA outer membrane protein assembly factor"):
    return {"target": target, "theader": header}


def target_row(group="omp85_bama", **extra):
    """A validated query row as it appears in `target_meta`: curator keys survive."""
    return {"accession": "T1", "mechanism_group": group, **extra}


# --- defect 1: asserted evidence must not promote --------------------------


def test_curator_asserted_rbh_does_not_promote_to_a_functional_claim():
    """`"rbh": true` in a manifest is an assertion, not a measurement."""
    label, basis, _, _ = classify_hit(
        query(), hit(), WHOLE, target_row(rbh=True),
    )
    assert label != "probable_same_function", (
        "a curator-supplied rbh field promoted a pair to the strongest "
        f"non-exact functional label, with basis {basis!r}"
    )
    assert label == "same_mechanism_class"


def test_asserted_rbh_cannot_be_laundered_through_an_unknown_key():
    """Unknown keys survive `{**q, ...}` into target_meta; none may reach the label."""
    for key in ("rbh", "RBH", "reciprocal_best_hit"):
        label, _, _, _ = classify_hit(
            query(), hit(), WHOLE, target_row(**{key: True}),
        )
        assert label != "probable_same_function", f"{key!r} promoted the pair"


# --- defect 2: the computed path must reach the label ----------------------


def test_computed_rbh_promotes_when_passed_explicitly():
    """The label is reachable by evidence, through an explicit argument."""
    label, basis, evidence, limit = classify_hit(
        query(), hit(), WHOLE, target_row(), rbh=True,
    )
    assert label == "probable_same_function"
    assert evidence == "E3_curated_plus_E4_computational"
    assert "still require direct validation" in limit


def test_computed_rbh_does_not_promote_below_whole_architecture():
    """RBH alone is not enough; the structural category still bounds the claim."""
    label, _, _, _ = classify_hit(
        query(), hit(), "domain_or_partial_match", target_row(), rbh=True,
    )
    assert label != "probable_same_function"


def test_computed_rbh_does_not_promote_across_mechanism_groups():
    label, _, _, _ = classify_hit(
        query(), hit(), WHOLE, target_row(group="t9ss_porg"), rbh=True,
    )
    assert label != "probable_same_function"


# --- defect 3: the relation is pairwise ------------------------------------


def test_rbh_is_pairwise_not_a_property_of_the_target():
    """The same target row must classify differently for different queries.

    This is the defect that survives a naive fix. Moving the RBH computation
    earlier and writing it into `target_meta` produces a working path that is
    still wrong, because `target_meta` is keyed by target accession alone.
    """
    shared_target = target_row()
    promoted, _, _, _ = classify_hit(
        query("Q1"), hit(), WHOLE, shared_target, rbh=True,
    )
    not_promoted, _, _, _ = classify_hit(
        query("Q2"), hit(), WHOLE, shared_target, rbh=False,
    )
    assert promoted == "probable_same_function"
    assert not_promoted == "same_mechanism_class"


def test_default_is_unpromoted():
    """Omitting the flag must never promote: absence of evidence is not evidence."""
    label, _, _, _ = classify_hit(query(), hit(), WHOLE, target_row())
    assert label == "same_mechanism_class"


# --- the manifest gate -----------------------------------------------------


def test_query_records_reject_reserved_computed_fields():
    """Curator manifests may not carry fields the pipeline computes."""
    from sf_csa.manifests import RESERVED_COMPUTED_FIELDS, reject_reserved_fields

    assert "rbh" in RESERVED_COMPUTED_FIELDS
    reject_reserved_fields({"accession": "Q1", "mechanism_group": "x"}, "query Q1")
    with pytest.raises(Exception) as exc:
        reject_reserved_fields({"accession": "Q1", "rbh": True}, "query Q1")
    assert "rbh" in str(exc.value)
