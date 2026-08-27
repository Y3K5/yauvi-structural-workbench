#!/usr/bin/env python3
"""Property-based tests for the activity-state classifier.

The example-based suites in this package check that named fixtures produce named
labels. These tests check the things that must hold for *every* input, which is
what the module's scientific claims actually are: a closed vocabulary is a claim
about all inputs, not about six of them.

Six invariants from the module docstring and the platform's fail-closed policy
are encoded here as properties:

  P1  the label is always drawn from the closed vocabulary, and `assess` is total
      (no input raises)
  P2  every signal reports a state from `SIGNAL_STATES` and explains itself
  P3  a degraded catalytic position is decisive: no other evidence, however
      favourable, can lift the label off `active_site_disrupted`
  P4  no annotated catalytic site always yields `indeterminate` — never any
      `inactive` label (absence of annotation is not absence of function)
  P5  a predicted model can never reach `active_state_supported`, and neither can
      any input with a non-informative signal
  P6  result documents are byte-identical across runs on the same input

On monotonicity
---------------
The brief asks for "removing a signal never improves a label". That is *false*
of this classifier as written, and deliberately so — see
`test_removing_contradicting_evidence_does_raise_the_label`, which pins the
behaviour rather than asserting it away. Dropping a structure that showed a
dispersed site moves `inactive_conformation` to `probable_active`, because the
reason for the negative claim went with it.

So monotonicity is asserted where it is meaningful: on inputs with no
contradicted signal, removing evidence never *increases* positive claim
strength (`test_removing_evidence_never_strengthens_a_positive_claim`), and
removing evidence can never *reach* the strongest label
(`test_removing_evidence_never_reaches_the_strongest_label`). The distinction is
not pedantry: the first is a property of the ordering, the second is invariant 5
of the brief, and only the second is a fail-closed guarantee.
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from hypothesis import assume, given, settings, strategies as st

from actstate.core import (
    CATALYTICALLY_COMPETENT,
    LABELS,
    SIGNAL_STATES,
    SITE_CLUSTER_MAX_ANGSTROM,
    ProteinRecord,
    assess,
    geometry_signal,
)
from actstate.features import parse_features
from actstate.io import build_document, write_json
from actstate.structure import THREE_TO_ONE, Heteroatom, Residue, Structure

# Deterministic runs: a property suite that finds a different counterexample on
# every CI run is not a regression test.
settings.register_profile("actstate", settings(max_examples=150, deadline=None))
settings.load_profile("actstate")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items() if len(k) == 3 and v != "M"}
ONE_TO_THREE["M"] = "MET"

# How strongly each label asserts that the protein is in a functionally
# competent state. This is *not* a total order on the vocabulary and is not
# claimed to be one: `inactive_conformation` and `active_site_disrupted` are
# both "no positive claim", and they differ in what they rule out, not in
# strength. The ranking exists only to make "never stronger than" testable.
POSITIVE_STRENGTH = {
    "active_state_supported": 4,
    "probable_active": 3,
    "apo_but_competent": 2,
    "indeterminate": 1,
    "inactive_conformation": 0,
    "active_site_disrupted": 0,
}

# Labels that assert the protein is *not* in a working state. P4 forbids
# reaching either of these from missing annotation.
NEGATIVE_LABELS = frozenset({"inactive_conformation", "active_site_disrupted"})

COMPETENT = "".join(sorted(CATALYTICALLY_COMPETENT))
INCOMPETENT = "".join(sorted(set(AMINO_ACIDS) - set(CATALYTICALLY_COMPETENT)))


# -- strategies -----------------------------------------------------------


def act_site_string(positions, *, experimental=True) -> str:
    """Build a UniProt-style ACT_SITE column for the given 1-based positions."""
    code = "ECO:0000269|PubMed:1" if experimental else "ECO:0000255"
    return "; ".join(
        f'ACT_SITE {p}; /note="Generated"; /evidence="{code}"' for p in positions
    )


@st.composite
def sequences(draw, min_size=4, max_size=40):
    return draw(st.text(alphabet=AMINO_ACIDS, min_size=min_size, max_size=max_size))


@st.composite
def catalytic_positions(draw, sequence, *, min_count=1, max_count=4):
    """1-based positions inside the sequence, sorted and distinct."""
    count = draw(st.integers(min_value=min_count, max_value=max_count))
    return sorted(
        draw(
            st.lists(
                st.integers(min_value=1, max_value=len(sequence)),
                min_size=count,
                max_size=count,
                unique=True,
            )
        )
    )


@st.composite
def records(draw, *, annotated=None, competent=None):
    """A ProteinRecord, optionally forced to be annotated and/or site-competent.

    `annotated=None` lets Hypothesis choose, so the unconstrained properties see
    both regimes.
    """
    sequence = list(draw(sequences()))
    if annotated is None:
        annotated = draw(st.booleans())

    positions: list[int] = []
    if annotated:
        positions = draw(catalytic_positions(sequence))
        for position in positions:
            if competent is True:
                sequence[position - 1] = draw(st.sampled_from(COMPETENT))
            elif competent is False:
                sequence[position - 1] = draw(st.sampled_from(INCOMPETENT))

    cofactor = draw(st.sampled_from(["", "COFACTOR: Name=Zn(2+); Xref=CHEBI:29105;"]))
    return ProteinRecord(
        accession=draw(st.text(alphabet="ABCDEFGHIJ0123456789_", min_size=1, max_size=8)),
        sequence="".join(sequence),
        act_site_raw=act_site_string(positions) if positions else "",
        cofactor_raw=cofactor,
    ), positions


@st.composite
def structures(draw, record, positions, *, clustered=None, predicted=None, heteroatoms=None):
    """Coordinates for a record, with the catalytic residues placed on purpose.

    Single chain throughout: multi-chain behaviour is a separate question (see
    REVIEW.md) and mixing it in here would make every property about chains.
    """
    if clustered is None:
        clustered = draw(st.booleans())
    if predicted is None:
        predicted = draw(st.booleans())
    if heteroatoms is None:
        heteroatoms = draw(st.booleans())

    # Clustered: inside a box well under the bound. Dispersed: far enough apart
    # that the widest pair must exceed it regardless of which pair is widest.
    spread = 3.0 if clustered else 4.0 * SITE_CLUSTER_MAX_ANGSTROM

    residues = []
    for index, position in enumerate(positions or [1]):
        letter = record.sequence[position - 1] if position <= len(record.sequence) else "A"
        residues.append(
            Residue(
                seq_id=position,
                name=ONE_TO_THREE.get(letter, "ALA"),
                chain="A",
                x=float(index * spread),
                y=0.0,
                z=0.0,
                atom_used="CB",
                b_factor=draw(st.floats(min_value=0.0, max_value=100.0)),
            )
        )

    hets = ()
    if heteroatoms:
        hets = (Heteroatom(name="ZN", chain="A", seq_id=900, atom_count=1),)

    return Structure(
        identifier=record.accession,
        residues=tuple(residues),
        heteroatoms=hets,
        is_predicted=predicted,
        source_note="generated" if predicted else "",
    )


@st.composite
def full_case(draw):
    """A record, its coordinates, and both sidecars — every knob free."""
    record, positions = draw(records())
    has_structure = draw(st.booleans())
    structure = draw(structures(record, positions)) if has_structure else None
    comparison = draw(
        st.one_of(
            st.none(),
            st.just({"reference": "0REF", "state": "active", "score": 0.9}),
            st.just({"reference": "0REF", "state": "inactive", "score": 0.9}),
            st.just({"reference": "", "state": "sideways"}),
        )
    )
    fold_state = draw(
        st.one_of(
            st.none(),
            st.just({"state": "active_assembly"}),
            st.just({"state": "isolated_fold"}),
            st.just({"state": "something_else"}),
        )
    )
    return record, structure, comparison, fold_state


# -- P1: the vocabulary is closed and `assess` is total --------------------


@given(full_case())
def test_label_is_always_in_the_closed_vocabulary(case):
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    assert result.label in LABELS


@given(full_case())
def test_assess_is_total(case):
    """No input in the generated space raises. Fail-closed means a label, not a crash."""
    record, structure, comparison, fold_state = case
    assess(record, structure=structure, reference_comparison=comparison, fold_state=fold_state)


@given(full_case())
def test_every_signal_is_reported_exactly_once(case):
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    names = [s.name for s in result.signals]
    assert sorted(names) == ["assembly", "completeness", "conformation", "geometry", "occupancy"]


@given(st.text(max_size=200))
def test_arbitrary_feature_text_never_escapes_the_vocabulary(raw):
    """Free text in the annotation column is never promoted into a label."""
    result = assess(ProteinRecord(accession="P", sequence="ACDEFGHIK", act_site_raw=raw))
    assert result.label in LABELS


@given(st.text(alphabet=AMINO_ACIDS, max_size=30), st.text(max_size=100))
def test_arbitrary_cofactor_text_never_escapes_the_vocabulary(sequence, cofactor):
    result = assess(
        ProteinRecord(accession="P", sequence=sequence, act_site_raw="", cofactor_raw=cofactor)
    )
    assert result.label in LABELS


# -- P2: every signal has a state and a reason ----------------------------


@given(full_case())
def test_signal_states_are_from_the_closed_set(case):
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    for signal in result.signals:
        assert signal.state in SIGNAL_STATES


@given(full_case())
def test_every_signal_explains_itself(case):
    """An unavailable signal that does not say why is indistinguishable from a bug."""
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    for signal in result.signals:
        assert signal.detail.strip()


@given(full_case())
def test_a_label_always_carries_a_rationale(case):
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    assert result.rationale.strip()


@given(full_case())
def test_no_signal_is_ever_a_neutral_placeholder(case):
    """Fail-closed: a signal is informative or it is explicitly not evaluated.

    There is no fifth state meaning "assume fine", and an unavailable signal
    never carries values that would let a reader treat it as a pass.
    """
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    for signal in result.signals:
        if signal.state in ("unavailable", "unevaluated"):
            assert not signal.informative


# -- P3: a degraded catalytic position is decisive ------------------------


@given(records(annotated=True, competent=False))
def test_a_degraded_site_is_disrupted_whatever_else_is_supplied(case):
    """The strongest possible supporting evidence cannot rescue a dead site."""
    record, positions = case
    assume(positions)
    structure = Structure(
        identifier=record.accession,
        residues=tuple(
            Residue(
                seq_id=p,
                name=ONE_TO_THREE.get(record.sequence[p - 1], "ALA"),
                chain="A",
                x=float(i * 2.0),
                y=0.0,
                z=0.0,
                atom_used="CB",
            )
            for i, p in enumerate(positions)
        ),
        heteroatoms=(Heteroatom(name="ZN", chain="A", seq_id=900, atom_count=1),),
        is_predicted=False,
    )
    result = assess(
        record,
        structure=structure,
        reference_comparison={"reference": "0REF", "state": "active", "score": 1.0},
        fold_state={"state": "active_assembly"},
    )
    assert result.label == "active_site_disrupted"


@given(records(annotated=True, competent=False))
def test_a_degraded_site_names_the_offending_position(case):
    record, positions = case
    assume(positions)
    signal = assess(record).signal("completeness")
    assert signal.state == "contradicted"
    assert any(str(p) in signal.detail for p in positions)


# -- P4: absence of annotation is not absence of function -----------------


@given(records(annotated=False))
def test_no_annotation_is_always_indeterminate(case):
    record, _ = case
    assume(not record.features().catalytic_positions())
    assert assess(record).label == "indeterminate"


@given(records(annotated=False))
def test_no_annotation_never_yields_a_negative_label(case):
    """The rule that separates this module from a naive pseudoenzyme caller."""
    record, _ = case
    assume(not record.features().catalytic_positions())
    structure = Structure(
        identifier="S",
        residues=(Residue(seq_id=1, name="ALA", chain="A", x=0.0, y=0.0, z=0.0, atom_used="CB"),),
        heteroatoms=(),
        is_predicted=False,
    )
    for comparison in (None, {"reference": "0REF", "state": "inactive"}):
        result = assess(record, structure=structure, reference_comparison=comparison)
        assert result.label not in NEGATIVE_LABELS
        assert result.label == "indeterminate"


@given(sequences(), st.integers(min_value=1, max_value=6))
def test_annotation_beyond_the_sequence_is_indeterminate_not_disrupted(sequence, overshoot):
    """A position off the end means annotation and sequence disagree.

    That is a mismatch to report, not evidence of a broken site.
    """
    position = len(sequence) + overshoot
    record = ProteinRecord(
        accession="P", sequence=sequence, act_site_raw=act_site_string([position])
    )
    result = assess(record)
    assert result.label == "indeterminate"
    assert result.signal("completeness").state == "unavailable"


# -- P5: predictions and missing signals cannot reach the top label -------


@given(records(annotated=True, competent=True))
def test_a_predicted_model_never_reaches_the_strongest_label(case):
    record, positions = case
    assume(positions)
    structure = Structure(
        identifier=record.accession,
        residues=tuple(
            Residue(seq_id=p, name=ONE_TO_THREE.get(record.sequence[p - 1], "ALA"),
                    chain="A", x=float(i * 2.0), y=0.0, z=0.0, atom_used="CB")
            for i, p in enumerate(positions)
        ),
        heteroatoms=(Heteroatom(name="ZN", chain="A", seq_id=900, atom_count=1),),
        is_predicted=True,
        source_note="generated",
    )
    result = assess(
        record,
        structure=structure,
        reference_comparison={"reference": "0REF", "state": "active", "score": 1.0},
        fold_state={"state": "active_assembly"},
    )
    assert result.label != "active_state_supported"


@given(full_case())
def test_the_strongest_label_requires_nothing_unavailable_or_contradicted(case):
    """`active_state_supported` is only reachable when nothing is missing or against.

    Note the exact claim: no signal may be `unavailable` (could not be checked)
    and none `contradicted`. `unevaluated` (not applicable to this entry) does
    *not* block — which is a weaker guarantee than the module docstring's "every
    signal available and consistent" suggests. See
    `test_the_strongest_label_can_be_reached_without_checking_geometry`.
    """
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    if result.label == "active_state_supported":
        assert not any(s.state == "unavailable" for s in result.signals)
        assert not any(s.state == "contradicted" for s in result.signals)
        assert structure is not None and not structure.is_predicted


def test_the_strongest_label_can_be_reached_without_checking_geometry():
    """Documenting a gap Hypothesis found, rather than asserting it away.

    `assign_label` blocks the strongest label on any `unavailable` signal but not
    on an `unevaluated` one. For occupancy that is defensible: an entry with no
    declared cofactor and no heteroatom present poses no occupancy question.

    For geometry it is harder to defend. A single annotated catalytic position
    yields `unevaluated` ("no geometry to check against itself"), so an entry
    with one ACT_SITE reaches `active_state_supported` with the geometry signal
    never evaluated — the very signal that would detect a dispersed site. Two of
    five signals can be skipped this way while the label still claims "every
    evaluated signal is consistent".

    The rationale text is literally true (every *evaluated* signal was
    consistent) and the docstring's "every signal available and consistent" is
    not. This test pins the behaviour so that closing the gap is a deliberate
    decision with a visible diff. See REVIEW.md §4.
    """
    record = ProteinRecord(
        accession="ONE_SITE", sequence="AAASAAA", act_site_raw=act_site_string([4])
    )
    structure = Structure(
        identifier="ONE_SITE",
        residues=(
            Residue(seq_id=4, name="SER", chain="A", x=0.0, y=0.0, z=0.0, atom_used="CB"),
        ),
        heteroatoms=(),
        is_predicted=False,
    )
    result = assess(
        record,
        structure=structure,
        reference_comparison={"reference": "0REF", "state": "active"},
        fold_state={"state": "active_assembly"},
    )
    assert result.label == "active_state_supported"
    assert result.signal("geometry").state == "unevaluated"
    assert result.signal("occupancy").state == "unevaluated"
    # Two of five signals carried no information, yet the strongest label stands.
    assert sum(1 for s in result.signals if not s.informative) == 2


def test_an_isolated_fold_cannot_reach_the_strongest_label():
    """Regression: the assembly signal must actually reach the label.

    Found by `test_the_strongest_label_requires_nothing_unavailable_or_contradicted`.
    `assembly_signal` was computed and reported, but `assign_label` never read it,
    so an entry whose fold_state said `isolated_fold` — the evidence describes the
    monomer, not the working assembly — was labelled `active_state_supported`
    with the rationale "every evaluated signal is consistent with a working
    state". The assembly signal said the opposite, in the same document.

    A contradicted assembly signal caps the claim at `probable_active`. It does
    not invert it to `inactive_conformation`: fold_state reports what the
    evidence describes, not that the monomer fold is a non-functional pose.
    """
    record = ProteinRecord(
        accession="MONOMER", sequence="AAAAHAAADAAAASAAAAA",
        act_site_raw=act_site_string([5, 9, 14]),
    )
    structure = Structure(
        identifier="MONOMER",
        residues=(
            Residue(seq_id=5, name="HIS", chain="A", x=0.0, y=0.0, z=0.0, atom_used="CB"),
            Residue(seq_id=9, name="ASP", chain="A", x=4.0, y=0.0, z=0.0, atom_used="CB"),
            Residue(seq_id=14, name="SER", chain="A", x=2.0, y=3.5, z=0.0, atom_used="CB"),
        ),
        heteroatoms=(),
        is_predicted=False,
    )
    comparison = {"reference": "0REF", "state": "active"}

    isolated = assess(
        record, structure=structure, reference_comparison=comparison,
        fold_state={"state": "isolated_fold"},
    )
    assert isolated.signal("assembly").state == "contradicted"
    assert isolated.label == "probable_active"
    # The rationale must name the reason, not contradict it.
    assert "isolated monomer" in isolated.rationale

    # The same entry in its working assembly does reach the strongest label,
    # so the gate is specific rather than a blanket downgrade.
    assembled = assess(
        record, structure=structure, reference_comparison=comparison,
        fold_state={"state": "active_assembly"},
    )
    assert assembled.label == "active_state_supported"


@given(full_case())
def test_a_contradicted_signal_is_always_visible_in_the_label(case):
    """No signal may argue against the protein while the label claims support."""
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    if any(s.state == "contradicted" for s in result.signals):
        assert result.label != "active_state_supported"


@given(full_case())
def test_the_strongest_label_requires_experimental_coordinates(case):
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    if result.label == "active_state_supported":
        assert structure is not None
        assert not structure.is_predicted


@given(full_case())
def test_an_unavailable_signal_is_named_in_the_rationale(case):
    """A signal that could not be evaluated never silently drops out."""
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    if result.label == "probable_active":
        missing = [s.name for s in result.signals if s.state == "unavailable"]
        if missing and "predicted model" not in result.rationale:
            assert any(name in result.rationale for name in missing)


# -- monotonicity, stated precisely --------------------------------------


@given(full_case())
def test_removing_evidence_never_reaches_the_strongest_label(case):
    """Invariant 5, as a property: less evidence cannot produce the top label.

    Every one-piece-removed variant of every generated case is checked, so this
    covers dropping the structure, either sidecar, or the sequence.
    """
    record, structure, comparison, fold_state = case
    variants = [
        assess(record, structure=None, reference_comparison=comparison, fold_state=fold_state),
        assess(record, structure=structure, reference_comparison=None, fold_state=fold_state),
        assess(record, structure=structure, reference_comparison=comparison, fold_state=None),
        assess(
            ProteinRecord(
                accession=record.accession,
                sequence="",
                act_site_raw=record.act_site_raw,
                cofactor_raw=record.cofactor_raw,
            ),
            structure=structure,
            reference_comparison=comparison,
            fold_state=fold_state,
        ),
    ]
    for result in variants:
        assert result.label != "active_state_supported"


@given(full_case())
def test_removing_evidence_never_strengthens_a_positive_claim(case):
    """Monotonicity, on the sublattice where it is meaningful.

    Restricted to cases with no contradicted signal — where nothing is arguing
    against the protein, removing evidence can only weaken the claim. The
    unrestricted version is false by design; see the next test.
    """
    record, structure, comparison, fold_state = case
    full = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    assume(all(s.state != "contradicted" for s in full.signals))

    for reduced in (
        assess(record, structure=None, reference_comparison=comparison, fold_state=fold_state),
        assess(record, structure=structure, reference_comparison=None, fold_state=fold_state),
        assess(record, structure=structure, reference_comparison=comparison, fold_state=None),
        assess(record, structure=None, reference_comparison=None, fold_state=None),
    ):
        assume(all(s.state != "contradicted" for s in reduced.signals))
        assert POSITIVE_STRENGTH[reduced.label] <= POSITIVE_STRENGTH[full.label]


def test_removing_contradicting_evidence_does_raise_the_label():
    """Pinning the intended non-monotonicity, so a future change is a decision.

    A dispersed active site yields `inactive_conformation`. Drop the coordinates
    and the label rises to `probable_active`, because the evidence for the
    negative claim went with them. This is correct — the module reports what it
    can see — but it means "more evidence is always a weaker claim" is not a
    property of this classifier, and a reader comparing two runs with different
    inputs cannot treat the labels as ordered.
    """
    sequence = "AAAAHAAADAAAASAAAAA"
    record = ProteinRecord(
        accession="P", sequence=sequence, act_site_raw=act_site_string([5, 9, 14])
    )
    dispersed = Structure(
        identifier="P",
        residues=(
            Residue(seq_id=5, name="HIS", chain="A", x=0.0, y=0.0, z=0.0, atom_used="CB"),
            Residue(seq_id=9, name="ASP", chain="A", x=60.0, y=0.0, z=0.0, atom_used="CB"),
            Residue(seq_id=14, name="SER", chain="A", x=0.0, y=60.0, z=0.0, atom_used="CB"),
        ),
        heteroatoms=(),
        is_predicted=False,
    )
    with_structure = assess(record, structure=dispersed)
    without = assess(record, structure=None)

    assert with_structure.label == "inactive_conformation"
    assert without.label == "probable_active"
    assert POSITIVE_STRENGTH[without.label] > POSITIVE_STRENGTH[with_structure.label]


# -- P6: determinism ------------------------------------------------------


@given(st.lists(full_case(), min_size=1, max_size=4))
def test_result_documents_are_byte_identical_across_runs(cases):
    """Determinism is a documented guarantee with tests on it; this generalises them."""
    def run():
        assessments = [
            assess(r, structure=s, reference_comparison=c, fold_state=f)
            for r, s, c, f in cases
        ]
        return json.dumps(
            build_document(assessments, config={"chain": "all"}),
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        )

    assert run() == run()


@given(st.lists(full_case(), min_size=1, max_size=3))
def test_written_documents_are_byte_identical_across_runs(cases):
    # A function-scoped tmp_path fixture cannot be combined with @given (the
    # directory would be shared across examples), so the scratch directory is
    # created per example instead.
    directory = Path(tempfile.mkdtemp())
    payloads = []
    for run in ("a", "b"):
        assessments = [
            assess(r, structure=s, reference_comparison=c, fold_state=f)
            for r, s, c, f in cases
        ]
        path = write_json(
            directory / f"{run}.json", build_document(assessments, config={"chain": "all"})
        )
        payloads.append(path.read_bytes())
    shutil.rmtree(directory, ignore_errors=True)
    assert payloads[0] == payloads[1]


@given(full_case())
def test_assessment_order_does_not_affect_the_document(case):
    """Sorting by accession means input order cannot change the output."""
    record, structure, comparison, fold_state = case
    other = ProteinRecord(accession="ZZZ_LAST", sequence="ACDEFGHIK")
    one = assess(record, structure=structure, reference_comparison=comparison, fold_state=fold_state)
    two = assess(other)
    forward = build_document([one, two], config={})
    backward = build_document([two, one], config={})
    assert forward == backward


# -- the geometry bound behaves like a bound -----------------------------


@given(records(annotated=True, competent=True), st.floats(min_value=0.5, max_value=8.0))
def test_a_tight_cluster_is_never_reported_as_dispersed(case, spacing):
    record, positions = case
    assume(len(positions) >= 2)
    # Place residues on a line with spacing chosen so the widest pair is inside
    # the bound; that is the definition of clustered, so it must not contradict.
    span = spacing * (len(positions) - 1)
    assume(span <= SITE_CLUSTER_MAX_ANGSTROM)
    structure = Structure(
        identifier="S",
        residues=tuple(
            Residue(seq_id=p, name=ONE_TO_THREE.get(record.sequence[p - 1], "ALA"),
                    chain="A", x=float(i) * spacing, y=0.0, z=0.0, atom_used="CB")
            for i, p in enumerate(positions)
        ),
        heteroatoms=(),
        is_predicted=False,
    )
    assert assess(record, structure=structure).signal("geometry").state != "contradicted"


@given(records(annotated=True, competent=True), st.floats(min_value=1.01, max_value=20.0))
def test_a_site_wider_than_the_bound_is_always_contradicted(case, factor):
    record, positions = case
    assume(len(positions) >= 2)
    separation = SITE_CLUSTER_MAX_ANGSTROM * factor
    structure = Structure(
        identifier="S",
        residues=tuple(
            Residue(seq_id=p, name=ONE_TO_THREE.get(record.sequence[p - 1], "ALA"),
                    chain="A", x=float(i) * separation, y=0.0, z=0.0, atom_used="CB")
            for i, p in enumerate(positions)
        ),
        heteroatoms=(),
        is_predicted=False,
    )
    result = assess(record, structure=structure)
    assert result.signal("geometry").state == "contradicted"
    assert result.label == "inactive_conformation"


@given(records(annotated=True, competent=True))
def test_a_single_catalytic_position_has_no_geometry(case):
    """One residue cannot be dispersed, so the signal must not claim it is."""
    record, positions = case
    assume(len(positions) == 1)
    structure = Structure(
        identifier="S",
        residues=(
            Residue(seq_id=positions[0], name="SER", chain="A",
                    x=0.0, y=0.0, z=0.0, atom_used="CB"),
        ),
        heteroatoms=(),
        is_predicted=False,
    )
    assert assess(record, structure=structure).signal("geometry").state == "unevaluated"


# -- occupancy and the conformation signal -------------------------------


@given(records(annotated=True, competent=True))
def test_an_uninterpretable_reference_comparison_is_unavailable_not_a_pass(case):
    """Malformed input to the conformation signal fails closed."""
    record, positions = case
    assume(positions)
    # The token list is exactly ("active", "inactive"). Every other token,
    # including the plausible-looking ones a curator might write, must land in
    # `unavailable` rather than being interpreted. Mutation testing (mutant A9)
    # showed an earlier version of this test missed a widened token set,
    # because it never tried a token that a widened set would accept.
    for comparison in (
        {},
        {"reference": "0REF"},
        {"state": "active"},
        {"reference": "", "state": "active"},
        {"reference": "0REF", "state": "maybe"},
        {"reference": "0REF", "state": "ACTIVE"},
        {"reference": "0REF", "state": "unknown"},
        {"reference": "0REF", "state": "indeterminate"},
        {"reference": "0REF", "state": "none"},
        {"reference": "0REF", "state": ""},
        {"reference": "0REF", "state": "active_assembly"},
    ):
        signal = assess(record, reference_comparison=comparison).signal("conformation")
        assert signal.state == "unavailable"


@given(records(annotated=True, competent=True))
def test_an_unrecognised_fold_state_is_unavailable_not_a_pass(case):
    record, positions = case
    assume(positions)
    for fold_state in ({}, {"state": ""}, {"state": "monomer"}, {"state": "ACTIVE_ASSEMBLY"}):
        signal = assess(record, fold_state=fold_state).signal("assembly")
        assert signal.state == "unavailable"


@given(records(annotated=True, competent=True))
def test_solvent_alone_is_never_cofactor_occupancy(case):
    """Waters and buffer components are not evidence of an occupied site."""
    record, positions = case
    assume(positions)
    record = ProteinRecord(
        accession=record.accession,
        sequence=record.sequence,
        act_site_raw=record.act_site_raw,
        cofactor_raw="COFACTOR: Name=Zn(2+); Xref=CHEBI:29105;",
    )
    structure = Structure(
        identifier="S",
        residues=tuple(
            Residue(seq_id=p, name="SER", chain="A", x=float(i) * 2.0, y=0.0, z=0.0,
                    atom_used="CB")
            for i, p in enumerate(positions)
        ),
        heteroatoms=(
            Heteroatom(name="HOH", chain="A", seq_id=901, atom_count=1),
            Heteroatom(name="GOL", chain="A", seq_id=902, atom_count=6),
            Heteroatom(name="SO4", chain="A", seq_id=903, atom_count=5),
        ),
        is_predicted=False,
    )
    signal = assess(record, structure=structure).signal("occupancy")
    assert signal.state == "contradicted"


@given(full_case())
def test_a_declared_cofactor_is_always_reported(case):
    record, structure, comparison, fold_state = case
    result = assess(
        record, structure=structure, reference_comparison=comparison, fold_state=fold_state
    )
    assert tuple(result.declared_cofactors) == record.declared_cofactors()


# -- what the parser could not read is always surfaced -------------------


@given(st.lists(st.sampled_from(["?", "<5", "12..?", ">400", "?..9"]), min_size=1, max_size=4))
def test_unreadable_positions_are_counted_never_guessed(tokens):
    raw = "; ".join(f"ACT_SITE {t}" for t in tokens)
    result = assess(ProteinRecord(accession="P", sequence="ACDEFGHIKLMNP", act_site_raw=raw))
    assert len(result.unparsed_features) == len(tokens)
    # Nothing was parsed, so there is no site to judge.
    assert result.label == "indeterminate"


@given(sequences(), st.integers(min_value=1, max_value=10))
def test_a_readable_and_an_unreadable_position_both_surface(sequence, extra):
    position = min(extra, len(sequence))
    raw = f'{act_site_string([position])}; ACT_SITE ?'
    result = assess(ProteinRecord(accession="P", sequence=sequence, act_site_raw=raw))
    assert position in result.catalytic_positions
    assert result.unparsed_features


# -- a structure with several chains is ambiguous, not a free choice ------
#
# REVIEW.md §10. `Structure.by_seq_id` keys on seq_id alone, so before the fix a
# homodimer's chain B silently overwrote chain A and geometry was asserted from
# whichever chain the parser happened to read last — a confident claim measured
# on half the file. Ambiguity now fails closed.


TRIAD_POSITIONS = (5, 9, 14)
# H at 5, D at 9, S at 14 — all in CATALYTICALLY_COMPETENT, so completeness is
# never the reason a label is capped in the tests below.
SEQUENCE_WITH_DEFAULTS = "AAAAHAAADAAAASAAAAA"


def _multichain(chains, positions=TRIAD_POSITIONS):
    """The same triad repeated in each chain, translated well apart."""
    names = {5: "HIS", 9: "ASP", 14: "SER"}
    offsets = {chain: 60.0 * index for index, chain in enumerate(chains)}
    return Structure(
        identifier="MULTI",
        residues=tuple(
            Residue(
                seq_id=position,
                name=names.get(position, "SER"),
                chain=chain,
                x=offsets[chain] + (2.0 if index else 0.0),
                y=3.5 if index == 2 else 0.0,
                z=0.0,
                atom_used="CB",
            )
            for chain in chains
            for index, position in enumerate(positions)
        ),
        heteroatoms=(),
        is_predicted=False,
    )


@given(st.integers(min_value=2, max_value=6))
@settings(max_examples=25, deadline=None)
def test_catalytic_positions_in_several_chains_never_yield_a_geometry_verdict(count):
    """With no chain selected, multi-chain geometry is unavailable, not decided."""
    chains = [chr(ord("A") + index) for index in range(count)]
    signal = geometry_signal(parse_features(act_site_string(TRIAD_POSITIONS)),
                             _multichain(chains))
    assert signal.state == "unavailable"
    # The reason must name the chains, so a user can act on it.
    for chain in chains:
        assert chain in signal.detail


@given(st.integers(min_value=2, max_value=4))
@settings(max_examples=20, deadline=None)
def test_selecting_a_chain_restores_a_verdict(count):
    """The ambiguity is about the *choice*, not the coordinates."""
    chains = [chr(ord("A") + index) for index in range(count)]
    features = parse_features(act_site_string(TRIAD_POSITIONS))
    structure = _multichain(chains)
    for chain in chains:
        signal = geometry_signal(features, structure, chain=chain)
        assert signal.state == "supported", f"chain {chain}: {signal.detail}"


def test_a_homodimer_cannot_reach_the_strongest_label_without_a_chain():
    """The end-to-end consequence: an ambiguous structure caps the claim."""
    record = ProteinRecord(
        accession="DIMER",
        sequence=SEQUENCE_WITH_DEFAULTS,
        act_site_raw=act_site_string(TRIAD_POSITIONS),
    )
    result = assess(
        record,
        structure=_multichain(["A", "B"]),
        reference_comparison={"reference": "0REF", "state": "active"},
        fold_state={"state": "active_assembly"},
    )
    assert result.label == "probable_active"
    assert "geometry" in result.rationale
    # And with the ambiguity resolved, the strongest label is reachable again.
    resolved = assess(
        record,
        structure=_multichain(["A", "B"]),
        chain="A",
        reference_comparison={"reference": "0REF", "state": "active"},
        fold_state={"state": "active_assembly"},
    )
    assert resolved.label == "active_state_supported"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
