"""Absent density may not improve a geometry call.

`geometry_signal` bailed out only when fewer than two declared catalytic positions
were resolved. With three of four present it measured the cluster over those three
and returned `supported`, so dropping the one residue that sits far from the others
turned a failure into a pass. On a protein whose signature is a mobile catalytic
loop that is usually unmodelled, that scores a site as better the less of it was
observed.

The case that exposed it: human PRPS1 entries 8DBG, 8DBI and 8DBL declare catalytic
positions 130, 196, 220 and 221. R196 has no density in any chain of those three, and
in every entry that does resolve it, it sits ~20 A from D221. Dropping it turned a
20.64 A failure into a 10.48 A pass, and the label read "the site is intact".
"""
from __future__ import annotations

from actstate.core import SIGNAL_STATES, ProteinRecord, assess
from actstate.structure import read_structure
from conftest import sequence_with


# --- the regression ------------------------------------------------------------

def test_dropping_a_distant_catalytic_residue_does_not_manufacture_a_cluster(
    intact_record, write_pdb
):
    """The heart of it: the resolved pair clusters, the declared set does not."""
    # 5 and 14 sit 4 A apart. 9 is declared catalytic and simply absent — exactly
    # the shape of 8DBG, where 130/220/221 cluster and 196 has no density.
    partial = read_structure(write_pdb("partial.pdb", [
        (5, "HIS", (0.0, 0.0, 0.0)),
        (14, "SER", (2.0, 3.5, 0.0)),
    ]))
    signal = assess(intact_record, structure=partial).signal("geometry")
    assert signal.state == "unevaluable", (
        "a site measured over part of its declared set must not return a verdict")
    assert signal.values["missing_from_structure"] == [9]
    assert signal.values["declared_set_coverage"] == {"declared": 3, "resolved": 2}
    assert "max_separation_angstrom" not in signal.values


def test_the_label_does_not_claim_an_intact_site_on_partial_coverage(
    intact_record, write_pdb
):
    """`probable_active` reads "the site is intact". It is a claim about coordinates."""
    partial = read_structure(write_pdb("partial.pdb", [
        (5, "HIS", (0.0, 0.0, 0.0)),
        (14, "SER", (2.0, 3.5, 0.0)),
    ]))
    assessment = assess(intact_record, structure=partial)
    assert assessment.label == "indeterminate"
    assert "the site is intact" not in assessment.rationale
    assert "declared catalytic position" in assessment.rationale


def test_partial_coverage_cannot_be_lifted_by_other_signals(intact_record, write_pdb):
    """The branch returns rather than capping, so nothing downstream can raise it.

    This is what makes the ordering constraint real: supplying a conformational
    reference set removes one of the reasons that previously capped these entries at
    `probable_active`, and without this branch they would climb to
    `active_state_supported` — the strongest label in the vocabulary — on coordinates
    with no density at the catalytic residue.
    """
    partial = read_structure(write_pdb("partial.pdb", [
        (5, "HIS", (0.0, 0.0, 0.0)),
        (14, "SER", (2.0, 3.5, 0.0)),
    ]))
    assessment = assess(
        intact_record,
        structure=partial,
        reference_comparison={"reference": "0REF", "state": "active", "score": 0.99},
        fold_state={"state": "active_assembly", "source": "test"},
    )
    assert assessment.label == "indeterminate"


# --- the behaviour that must not change ----------------------------------------

def test_a_complete_declared_set_still_returns_a_verdict(intact_record, write_pdb,
                                                         clustered_triad):
    structure = read_structure(write_pdb("all.pdb", clustered_triad))
    signal = assess(intact_record, structure=structure).signal("geometry")
    assert signal.state == "supported"
    assert signal.values["missing_from_structure"] == []
    assert signal.values["declared_set_coverage"] == {"declared": 3, "resolved": 3}
    assert signal.values["max_separation_angstrom"] > 0


def test_a_complete_but_dispersed_set_is_still_contradicted(intact_record, write_pdb,
                                                            dispersed_triad):
    structure = read_structure(write_pdb("far.pdb", dispersed_triad))
    assessment = assess(intact_record, structure=structure)
    assert assessment.signal("geometry").state == "contradicted"
    assert assessment.label == "inactive_conformation"


def test_no_structure_is_still_unavailable_not_unevaluable(intact_record):
    """Scope guard: the two states mean different things and must stay distinct."""
    signal = assess(intact_record, structure=None).signal("geometry")
    assert signal.state == "unavailable"
    assert "no structure supplied" in signal.detail


def test_no_annotated_site_is_still_unevaluated(write_pdb, clustered_triad):
    record = ProteinRecord(accession="P_NOSITE", sequence=sequence_with(), act_site_raw="")
    structure = read_structure(write_pdb("all.pdb", clustered_triad))
    assert assess(record, structure=structure).signal("geometry").state == "unevaluated"


def test_unevaluable_is_in_the_closed_vocabulary():
    assert "unevaluable" in SIGNAL_STATES
    assert "unavailable" in SIGNAL_STATES
