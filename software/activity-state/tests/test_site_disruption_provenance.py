"""`active_site_disrupted` may only be reached from position-specific evidence.

The screen behind that label asked one question: is the residue at an annotated
ACT_SITE position outside `CATALYTICALLY_COMPETENT`, a 13-letter set covering
nucleophiles, acid/base pairs, metal ligands, and the two that act through
backbone geometry. Membership in a set that broad is not a position-specific
chemistry test, and the label it fed is the strongest negative claim in the
closed vocabulary. Two consequences, in opposite directions:

* It over-claims. Seven residues fall outside the set, so any of them at an
  annotated position produced `active_site_disrupted` -- without knowing which
  residue the annotation expects there, without the role, and without comparing
  against an experimentally validated ortholog. A sequence that is the wrong
  isoform, or numbered against a different entry, is indistinguishable from a
  pseudoenzyme by this test.
* It is silent on the more common degradation. A catalytic Cys to Ser, or His
  to Asn, stays inside the set, so a genuinely dead site reads as `supported`.

So the strong label now requires an expected residue for that position, and
without one a non-competent residue caps the label at `indeterminate` -- the
observation is kept, in the signal and in the rationale, and the claim is not
made. Recorded in CHANGES.md.
"""
from __future__ import annotations

import pytest

from actstate.core import ProteinRecord, assess, completeness_signal

from conftest import TRIAD_ACT_SITE, sequence_with


def record_with(**kwargs) -> ProteinRecord:
    return ProteinRecord(
        accession=kwargs.pop("accession", "P"),
        sequence=kwargs.pop("sequence", sequence_with()),
        act_site_raw=kwargs.pop("act_site_raw", TRIAD_ACT_SITE),
        **kwargs,
    )


# -- the over-claim ------------------------------------------------------


def test_a_non_competent_residue_alone_does_not_establish_disruption():
    # Position 5 holds Ala, which is outside the competence set. Before this
    # change that alone produced `active_site_disrupted`.
    result = assess(record_with(sequence=sequence_with(p5="A")))
    assert result.label == "indeterminate"
    assert result.label != "active_site_disrupted"


def test_the_observation_survives_the_cap():
    """Capping the label must not discard the finding. Rule 3."""
    result = assess(record_with(sequence=sequence_with(p5="A")))
    signal = result.signal("completeness")
    assert signal.state == "contradicted"
    assert "5A" in signal.detail
    assert signal.values["basis"] == "generic_competence_set"
    # The rationale has to say what would raise this to the strong label,
    # otherwise the reader cannot tell a capped claim from an absent one.
    assert "expected residue" in result.rationale
    assert "5" in result.rationale


def test_a_capped_claim_is_never_lifted_into_a_positive_one(write_pdb, clustered_triad):
    """The half of P3 that must survive: no favourable evidence rescues it."""
    from actstate.structure import read_structure

    structure = read_structure(
        write_pdb("e.pdb", clustered_triad, header="HEADER    HYDROLASE   01-JAN-26   0AAA\n")
    )
    result = assess(
        record_with(sequence=sequence_with(p5="A")),
        structure=structure,
        reference_comparison={"reference": "1ABC", "state": "active"},
        fold_state={"state": "active_assembly"},
    )
    assert result.label == "indeterminate"
    assert result.label not in ("active_state_supported", "probable_active", "apo_but_competent")


# -- the position-specific path ------------------------------------------


def test_an_expected_residue_mismatch_does_establish_disruption():
    result = assess(
        record_with(sequence=sequence_with(p5="A")),
        expected_residues={5: "H", 9: "D", 14: "S"},
    )
    assert result.label == "active_site_disrupted"
    signal = result.signal("completeness")
    assert signal.values["basis"] == "position_specific_expected_residue"
    assert "expected H" in signal.detail and "observed A" in signal.detail


def test_a_substitution_inside_the_competence_set_is_caught_when_expected_is_known():
    """The false negative the generic screen cannot see: catalytic Ser to Thr."""
    generic = assess(record_with(sequence=sequence_with(p14="T")))
    assert generic.signal("completeness").state == "supported"

    specific = assess(
        record_with(sequence=sequence_with(p14="T")),
        expected_residues={14: "S"},
    )
    assert specific.label == "active_site_disrupted"
    assert "14" in specific.signal("completeness").detail


def test_matching_the_expectation_is_recorded_as_position_specific():
    result = assess(record_with(), expected_residues={5: "H", 9: "D", 14: "S"})
    signal = result.signal("completeness")
    assert signal.state == "supported"
    assert signal.values["position_specific_confirmed"] == 3
    assert signal.values["basis"] == "position_specific_expected_residue"


def test_an_expectation_outside_the_competence_set_still_governs():
    """The ortholog is the authority on its own site, not the residue set."""
    result = assess(
        record_with(sequence=sequence_with(p5="A")),
        expected_residues={5: "A"},
    )
    assert result.label != "active_site_disrupted"
    assert result.signal("completeness").state == "supported"


# -- what a curator may not do -------------------------------------------


@pytest.mark.parametrize(
    "expected",
    [
        {99: "H"},          # not an annotated catalytic position
        {5: "B"},           # not a standard amino acid
        {5: "HIS"},         # three-letter code
        {5: ""},            # empty
        {5: None},          # not a string
        {"five": "H"},      # not a position
    ],
)
def test_unusable_expectations_are_rejected_and_recorded(expected):
    result = assess(record_with(sequence=sequence_with(p5="A")), expected_residues=expected)
    signal = result.signal("completeness")
    # Rejected outright: it cannot license the strong label, and it is named
    # rather than dropped, so a curator sees that their entry did nothing.
    assert result.label != "active_site_disrupted"
    assert signal.values["rejected_expectations"]


def test_assess_stays_total_on_a_malformed_expectation():
    """P1: assess raises on no input. Validation belongs at the IO boundary."""
    assess(record_with(), expected_residues={"x": object()})
    assess(record_with(), expected_residues="not a mapping")
    assess(record_with(), expected_residues=None)


def test_expectations_are_keyed_per_position_not_per_record():
    """One matching position cannot vouch for another that mismatches."""
    result = assess(
        record_with(sequence=sequence_with(p5="H", p14="A")),
        expected_residues={5: "H", 14: "S"},
    )
    assert result.label == "active_site_disrupted"
    assert "14" in result.signal("completeness").detail


# -- the signal keeps its other guarantees --------------------------------


def test_no_annotation_is_still_indeterminate_for_its_own_reason():
    result = assess(ProteinRecord(accession="P", sequence="AAAA"), expected_residues={1: "H"})
    assert result.label == "indeterminate"
    assert result.signal("completeness").state == "unevaluated"
    assert "no ACT_SITE annotation" in result.signal("completeness").detail


def test_completeness_signal_takes_the_expectation_directly():
    record = record_with(sequence=sequence_with(p9="A"))
    signal = completeness_signal(record, record.features(), expected_residues={9: "D"})
    assert signal.state == "contradicted"
    assert signal.values["basis"] == "position_specific_expected_residue"
