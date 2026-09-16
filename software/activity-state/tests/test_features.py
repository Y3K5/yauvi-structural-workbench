"""Parsing UniProt feature columns."""
from __future__ import annotations

from actstate.features import parse_cofactors, parse_features
from conftest import TRIAD_ACT_SITE


def test_parses_point_features():
    parsed = parse_features(TRIAD_ACT_SITE)
    assert parsed.catalytic_positions() == (5, 9, 14)
    assert len(parsed.catalytic) == 3
    assert parsed.has_catalytic_annotation


def test_keeps_the_note_and_the_evidence_code():
    first = parse_features(TRIAD_ACT_SITE).catalytic[0]
    assert first.note == "Charge relay system"
    assert "ECO:0000269" in first.evidence


def test_distinguishes_experimental_from_inferred_evidence():
    catalytic = parse_features(TRIAD_ACT_SITE).catalytic
    # ECO:0000269 is experimental; ECO:0000255 is inferred by similarity.
    assert [f.experimentally_evidenced for f in catalytic] == [True, False, True]


def test_parses_a_span_into_every_position():
    parsed = parse_features('BINDING 57..59; /ligand="ATP"; /ligand_id="ChEBI:CHEBI:30616"')
    assert parsed.ligand_positions() == (57, 58, 59)
    assert parsed.ligands() == ("ATP",)


def test_an_unreadable_position_is_reported_not_guessed():
    """UniProt uses '?', '<' and '>' for uncertain positions."""
    parsed = parse_features('ACT_SITE ?; /note="unknown"; ACT_SITE 42; /note="known"')
    assert parsed.catalytic_positions() == (42,)
    assert parsed.unparsed == ("ACT_SITE ?",)


def test_open_ended_positions_are_unparsed():
    parsed = parse_features("ACT_SITE <1; ACT_SITE 5..>9")
    assert parsed.catalytic_positions() == ()
    assert len(parsed.unparsed) == 2


def test_a_reversed_span_is_unparsed():
    assert parse_features("ACT_SITE 90..10").unparsed == ("ACT_SITE 90..10",)


def test_empty_and_missing_columns_are_safe():
    for raw in (None, "", "   "):
        parsed = parse_features(raw)
        assert parsed.features == () and parsed.unparsed == ()
        assert not parsed.has_catalytic_annotation


def test_unknown_feature_keys_are_ignored():
    parsed = parse_features('SIGNAL 1..22; ACT_SITE 30; /note="x"')
    assert parsed.catalytic_positions() == (30,)


def test_site_and_binding_are_kept_apart():
    parsed = parse_features('ACT_SITE 10; BINDING 20; /ligand="Mg"; SITE 30; /note="cleavage"')
    assert parsed.catalytic_positions() == (10,)
    assert parsed.ligand_positions() == (20,)
    assert len(parsed.of_kind("SITE")) == 1


def test_duplicate_positions_are_collapsed():
    assert parse_features("ACT_SITE 5; ACT_SITE 5").catalytic_positions() == (5,)


def test_cofactor_names_are_extracted():
    raw = (
        "COFACTOR: Name=Zn(2+); Xref=ChEBI:CHEBI:29105; Evidence=ECO:0000250; "
        "COFACTOR: Name=Mg(2+); Xref=ChEBI:CHEBI:18420;"
    )
    assert parse_cofactors(raw) == ("Mg(2+)", "Zn(2+)")


def test_no_cofactor_column_is_empty():
    assert parse_cofactors(None) == () and parse_cofactors("") == ()
