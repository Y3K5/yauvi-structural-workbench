"""The five signals, and the rules that turn them into a label."""
from __future__ import annotations

import pytest

from actstate.core import LABELS, ProteinRecord, assess
from actstate.structure import read_structure
from conftest import TRIAD_ACT_SITE, sequence_with


# --- completeness --------------------------------------------------------


def test_intact_site_supports_completeness(intact_record):
    signal = assess(intact_record).signal("completeness")
    assert signal.state == "supported"
    assert "5H" in signal.detail and "9D" in signal.detail


def test_a_degraded_position_contradicts_completeness():
    record = ProteinRecord(
        accession="P", sequence=sequence_with(p5="A", p9="A"), act_site_raw=TRIAD_ACT_SITE
    )
    signal = assess(record).signal("completeness")
    assert signal.state == "contradicted"
    assert signal.values["degraded"] == {"5": "A", "9": "A"}


def test_no_annotation_leaves_completeness_unevaluated():
    signal = assess(ProteinRecord(accession="P", sequence="A" * 20)).signal("completeness")
    assert signal.state == "unevaluated"


def test_no_sequence_makes_completeness_unavailable():
    record = ProteinRecord(accession="P", act_site_raw=TRIAD_ACT_SITE)
    assert assess(record).signal("completeness").state == "unavailable"


def test_annotation_beyond_the_sequence_is_a_mismatch_not_a_disruption():
    """A short sequence and a long annotation are not the same entry."""
    record = ProteinRecord(accession="P", sequence="AAA", act_site_raw=TRIAD_ACT_SITE)
    signal = assess(record).signal("completeness")
    assert signal.state == "unavailable"
    assert signal.values["out_of_range"] == [5, 9, 14]


def test_experimental_evidence_is_counted(intact_record):
    assert assess(intact_record).signal("completeness").values["experimentally_evidenced"] == 2


# --- geometry ------------------------------------------------------------


def test_clustered_residues_support_geometry(intact_record, write_pdb, clustered_triad):
    structure = read_structure(write_pdb("x.pdb", clustered_triad))
    signal = assess(intact_record, structure=structure).signal("geometry")
    assert signal.state == "supported"
    assert signal.values["max_separation_angstrom"] < 16.0


def test_dispersed_residues_contradict_geometry(intact_record, write_pdb, dispersed_triad):
    structure = read_structure(write_pdb("x.pdb", dispersed_triad))
    signal = assess(intact_record, structure=structure).signal("geometry")
    assert signal.state == "contradicted"
    assert signal.values["max_separation_angstrom"] > 16.0


def test_no_structure_makes_geometry_unavailable(intact_record):
    assert assess(intact_record).signal("geometry").state == "unavailable"


def test_a_single_catalytic_position_has_no_geometry(write_pdb):
    record = ProteinRecord(accession="P", sequence=sequence_with(), act_site_raw="ACT_SITE 5")
    structure = read_structure(write_pdb("x.pdb", [(5, "HIS", (0.0, 0.0, 0.0))]))
    assert assess(record, structure=structure).signal("geometry").state == "unevaluated"


def test_unresolved_positions_block_geometry_rather_than_biasing_it(
    intact_record, write_pdb
):
    """Two of three residues missing from the coordinates is not a verdict."""
    structure = read_structure(write_pdb("x.pdb", [(5, "HIS", (0.0, 0.0, 0.0))]))
    signal = assess(intact_record, structure=structure).signal("geometry")
    assert signal.state == "unavailable"
    assert signal.values["missing_from_structure"] == [9, 14]


def test_the_cluster_bound_is_configurable(intact_record, write_pdb, dispersed_triad):
    structure = read_structure(write_pdb("x.pdb", dispersed_triad))
    signal = assess(intact_record, structure=structure, max_separation=100.0).signal("geometry")
    assert signal.state == "supported"


# --- occupancy -----------------------------------------------------------


def test_declared_cofactor_absent_from_coordinates_is_apo(write_pdb, clustered_triad):
    record = ProteinRecord(
        accession="P",
        sequence=sequence_with(),
        act_site_raw=TRIAD_ACT_SITE,
        cofactor_raw="COFACTOR: Name=Zn(2+); Xref=ChEBI:CHEBI:29105;",
    )
    structure = read_structure(write_pdb("x.pdb", clustered_triad))
    assert assess(record, structure=structure).signal("occupancy").state == "contradicted"


def test_buffer_only_structure_is_still_apo(write_pdb, clustered_triad):
    record = ProteinRecord(
        accession="P",
        sequence=sequence_with(),
        act_site_raw=TRIAD_ACT_SITE,
        cofactor_raw="COFACTOR: Name=Zn(2+);",
    )
    structure = read_structure(
        write_pdb(
            "x.pdb", clustered_triad, heteroatoms=[("GOL", (9.0, 9.0, 9.0), "C")]
        )
    )
    assert assess(record, structure=structure).signal("occupancy").state == "contradicted"


def test_declared_and_present_cofactor_supports_occupancy(write_pdb, clustered_triad):
    record = ProteinRecord(
        accession="P",
        sequence=sequence_with(),
        act_site_raw=TRIAD_ACT_SITE,
        cofactor_raw="COFACTOR: Name=Zn(2+);",
    )
    structure = read_structure(
        write_pdb("x.pdb", clustered_triad, heteroatoms=[("ZN", (2.0, 1.5, 1.0), "ZN")])
    )
    assert assess(record, structure=structure).signal("occupancy").state == "supported"


def test_no_cofactor_declared_or_present_is_not_a_question(
    intact_record, write_pdb, clustered_triad
):
    structure = read_structure(write_pdb("x.pdb", clustered_triad))
    assert assess(intact_record, structure=structure).signal("occupancy").state == "unevaluated"


# --- conformation and assembly ------------------------------------------


def test_conformation_is_unavailable_without_a_comparison(intact_record):
    """The fail-closed case: no aligner run means no conclusion, not a neutral one."""
    signal = assess(intact_record).signal("conformation")
    assert signal.state == "unavailable"
    assert "aligner" in signal.detail


def test_an_active_reference_supports_conformation(intact_record):
    signal = assess(
        intact_record, reference_comparison={"reference": "1ABC", "state": "active", "score": 0.9}
    ).signal("conformation")
    assert signal.state == "supported"
    assert signal.values["reference"] == "1ABC"


def test_an_inactive_reference_contradicts_conformation(intact_record):
    signal = assess(
        intact_record, reference_comparison={"reference": "2XYZ", "state": "inactive"}
    ).signal("conformation")
    assert signal.state == "contradicted"


def test_an_uninterpretable_comparison_is_unavailable(intact_record):
    signal = assess(intact_record, reference_comparison={"score": 0.9}).signal("conformation")
    assert signal.state == "unavailable"


def test_assembly_reads_fold_state(intact_record):
    assert (
        assess(intact_record, fold_state={"state": "active_assembly"}).signal("assembly").state
        == "supported"
    )
    assert (
        assess(intact_record, fold_state={"state": "isolated_fold"}).signal("assembly").state
        == "contradicted"
    )
    assert assess(intact_record).signal("assembly").state == "unavailable"


# --- labelling rules -----------------------------------------------------


def test_no_annotation_is_indeterminate_never_inactive():
    """Rule 1. Absence of annotation is not evidence of absence of function."""
    result = assess(ProteinRecord(accession="P", sequence="A" * 20))
    assert result.label == "indeterminate"


def test_a_disrupted_site_is_decisive_even_with_good_geometry(write_pdb, clustered_triad):
    record = ProteinRecord(
        accession="P", sequence=sequence_with(p5="A", p9="A"), act_site_raw=TRIAD_ACT_SITE
    )
    structure = read_structure(write_pdb("x.pdb", clustered_triad))
    result = assess(record, structure=structure)
    assert result.label == "active_site_disrupted"
    # Geometry still reports separately; it is not overwritten by the label.
    assert result.signal("geometry").state == "supported"


def test_a_predicted_model_can_never_be_active_state_supported(
    intact_record, write_pdb, clustered_triad
):
    """Rule 2."""
    structure = read_structure(
        write_pdb(
            "AF-P-F1.pdb", clustered_triad, header="HEADER    PREDICTED\nREMARK ALPHAFOLD\n"
        )
    )
    result = assess(
        intact_record,
        structure=structure,
        reference_comparison={"reference": "1ABC", "state": "active"},
        fold_state={"state": "active_assembly"},
    )
    assert result.label == "probable_active"
    assert "predicted model" in result.rationale


def test_full_evidence_on_experimental_coordinates_is_supported(write_pdb, clustered_triad):
    record = ProteinRecord(
        accession="P",
        sequence=sequence_with(),
        act_site_raw=TRIAD_ACT_SITE,
        cofactor_raw="COFACTOR: Name=Zn(2+);",
    )
    structure = read_structure(
        write_pdb(
            "e.pdb",
            clustered_triad,
            header="HEADER    HYDROLASE   01-JAN-26   0AAA\n",
            heteroatoms=[("ZN", (2.0, 1.5, 1.0), "ZN")],
        )
    )
    result = assess(
        record,
        structure=structure,
        reference_comparison={"reference": "1ABC", "state": "active"},
        fold_state={"state": "active_assembly"},
    )
    assert result.label == "active_state_supported"


def test_an_unavailable_signal_caps_the_label_at_probable(
    intact_record, write_pdb, clustered_triad
):
    """Rule 3: an unevaluated signal never quietly stops mattering."""
    structure = read_structure(
        write_pdb("e.pdb", clustered_triad, header="HEADER    HYDROLASE   01-JAN-26   0AAA\n")
    )
    result = assess(intact_record, structure=structure)
    assert result.label == "probable_active"
    assert "conformation" in result.rationale


def test_dispersed_geometry_is_an_inactive_conformation(
    intact_record, write_pdb, dispersed_triad
):
    structure = read_structure(write_pdb("d.pdb", dispersed_triad))
    assert assess(intact_record, structure=structure).label == "inactive_conformation"


def test_every_label_is_in_the_closed_vocabulary(intact_record):
    assert assess(intact_record).label in LABELS


def test_an_invented_label_is_rejected():
    from actstate.core import ActivityAssessment

    with pytest.raises(ValueError, match="closed vocabulary"):
        ActivityAssessment(accession="P", label="probably_fine", signals=())


def test_a_signal_must_explain_itself():
    from actstate.core import Signal

    with pytest.raises(ValueError, match="must explain itself"):
        Signal("completeness", "supported", "")


def test_an_unknown_signal_state_is_rejected():
    from actstate.core import Signal

    with pytest.raises(ValueError, match="unknown signal state"):
        Signal("completeness", "probably", "detail")


def test_every_assessment_reports_all_five_signals(intact_record):
    names = [s.name for s in assess(intact_record).signals]
    assert names == ["completeness", "geometry", "occupancy", "conformation", "assembly"]


def test_unparsed_features_are_surfaced_on_the_assessment():
    record = ProteinRecord(
        accession="P", sequence=sequence_with(), act_site_raw="ACT_SITE ?; ACT_SITE 5"
    )
    assert assess(record).unparsed_features == ("ACT_SITE ?",)
