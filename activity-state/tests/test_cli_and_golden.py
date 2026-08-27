"""Contract (T2) and golden (T3) tiers.

T2 asserts that what the module *says* it does matches what it does. T3 asserts
that running it on the bundled fixtures produces the same bytes every time.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from actstate import __version__
from actstate.cli import RESULT_NAME, TABLE_NAME, EXIT_BLOCKED, EXIT_OK, EXIT_USAGE, main
from actstate.core import LABELS, SIGNAL_STATES
from actstate.io import read_annotation_table, read_fasta

MODULE_SOURCES = Path(__file__).resolve().parent.parent / "src" / "actstate" / "sources.yaml"


# --- T2: the declared contract ------------------------------------------


@pytest.fixture
def described(capsys):
    assert main(["describe"]) == EXIT_OK
    return json.loads(capsys.readouterr().out)


def test_describe_is_valid_json_and_names_the_module(described):
    assert described["module_id"] == "activity_state"
    assert described["version"] == __version__


def test_describe_lists_the_closed_label_vocabulary(described):
    assert described["labels"] == list(LABELS)
    assert described["signal_states"] == list(SIGNAL_STATES)


def test_describe_names_every_signal_the_code_emits(described, intact_record):
    from actstate.core import assess

    declared = [s["name"] for s in described["signals"]]
    emitted = [s.name for s in assess(intact_record).signals]
    assert declared == emitted


def test_describe_declares_the_outputs_run_actually_writes(described):
    names = {o["name"] for o in described["outputs"]}
    assert names == {RESULT_NAME, TABLE_NAME}


def test_describe_states_its_limitations_honestly(described):
    text = " ".join(described["limitations"])
    assert "never be merged into a single" in text
    assert "indeterminate, never inactive" in text
    assert "never on its own yield active_state_supported" in text
    # The conformation signal has no reference set; that must be declared.
    assert "No such set ships with this module" in text


def test_the_source_manifest_ships_with_the_package():
    from importlib.resources import files

    assert (files("actstate") / "sources.yaml").is_file()


def test_fetch_plan_prints_the_declared_sources(capsys):
    assert main(["fetch", "--plan"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "module_id: actstate" in out
    assert "uniprot_proteomes" in out


# --- input handling ------------------------------------------------------


def test_column_aliases_are_accepted(tmp_path):
    """A UniProt export uses field names; the web UI uses display headers."""
    path = tmp_path / "a.tsv"
    path.write_text("Entry\tActive site\n P1 \tACT_SITE 5\n", encoding="utf-8")
    records = read_annotation_table(path)
    assert records[0].accession == "P1"
    assert records[0].act_site_raw == "ACT_SITE 5"


def test_a_table_without_an_accession_column_is_rejected(tmp_path):
    path = tmp_path / "a.tsv"
    path.write_text("name\tvalue\nx\ty\n", encoding="utf-8")
    with pytest.raises(Exception, match="accession column"):
        read_annotation_table(path)


def test_fasta_headers_in_uniprot_pipe_form_are_keyed_by_accession(tmp_path):
    path = tmp_path / "s.fasta"
    path.write_text(">sp|P12345|NAME_ORG\nMKV\nLAA\n>P99999\nMQQ\n", encoding="utf-8")
    sequences = read_fasta(path)
    assert sequences["P12345"] == "MKVLAA"
    assert sequences["P99999"] == "MQQ"


def test_missing_annotation_is_a_usage_error(tmp_path, capsys):
    assert main(["run", "--in", str(tmp_path), "--out", str(tmp_path / "o")]) == EXIT_USAGE
    assert "no annotation table found" in capsys.readouterr().err


# --- T3: golden -----------------------------------------------------------


@pytest.fixture
def run_examples(examples, tmp_path):
    def run(out_name="out"):
        out = tmp_path / out_name
        code = main(
            [
                "run",
                "--in", str(examples),
                "--out", str(out),
                "--fold-state", str(examples / "fold_state.json"),
                "--reference-comparison", str(examples / "reference_comparison.json"),
            ]
        )
        assert code == EXIT_OK
        return out

    return run


def test_the_fixtures_exercise_every_label(run_examples):
    document = json.loads((run_examples() / RESULT_NAME).read_text(encoding="utf-8"))
    produced = set(document["summary"]["labels"])
    assert produced == set(LABELS), f"fixtures do not cover: {set(LABELS) - produced}"


@pytest.mark.parametrize(
    ("accession", "expected"),
    [
        ("P_ACTIVE", "active_state_supported"),
        ("P_PREDICTED", "probable_active"),
        ("P_APO", "apo_but_competent"),
        ("P_DISPERSED", "inactive_conformation"),
        ("P_DISRUPTED", "active_site_disrupted"),
        ("P_NOSITE", "indeterminate"),
    ],
)
def test_each_fixture_gets_its_intended_label(run_examples, accession, expected):
    document = json.loads((run_examples() / RESULT_NAME).read_text(encoding="utf-8"))
    labels = {a["accession"]: a["label"] for a in document["assessments"]}
    assert labels[accession] == expected


def test_output_is_byte_identical_across_two_runs(run_examples):
    """Determinism is what makes recording a digest meaningful."""
    first, second = run_examples("a"), run_examples("b")
    for name in (RESULT_NAME, TABLE_NAME):
        digest_a = hashlib.sha256((first / name).read_bytes()).hexdigest()
        digest_b = hashlib.sha256((second / name).read_bytes()).hexdigest()
        assert digest_a == digest_b, f"{name} is not deterministic"


def test_the_document_carries_no_absolute_paths(run_examples):
    text = (run_examples() / RESULT_NAME).read_text(encoding="utf-8")
    assert "/Users/" not in text and "/tmp/" not in text


def test_the_document_carries_no_timestamp(run_examples):
    """A timestamp would make every run differ and defeat the golden test."""
    document = json.loads((run_examples() / RESULT_NAME).read_text(encoding="utf-8"))
    serialised = json.dumps(document)
    assert "retrieved_at" not in serialised and "generated_at" not in serialised


def test_assessments_are_sorted_by_accession(run_examples):
    document = json.loads((run_examples() / RESULT_NAME).read_text(encoding="utf-8"))
    accessions = [a["accession"] for a in document["assessments"]]
    assert accessions == sorted(accessions)


def test_the_table_has_one_row_per_protein(run_examples):
    lines = (run_examples() / TABLE_NAME).read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 7  # header + six proteins


def test_signals_are_reported_separately_in_the_output(run_examples):
    """The central design rule, asserted on the artifact itself."""
    document = json.loads((run_examples() / RESULT_NAME).read_text(encoding="utf-8"))
    for assessment in document["assessments"]:
        assert len(assessment["signals"]) == 5
        assert all(s["detail"] for s in assessment["signals"])
    # No aggregate score anywhere in the document.
    serialised = json.dumps(document)
    assert "activity_score" not in serialised and "total_score" not in serialised


def test_unevaluated_signals_are_recorded_not_dropped(run_examples):
    document = json.loads((run_examples() / RESULT_NAME).read_text(encoding="utf-8"))
    nosite = next(a for a in document["assessments"] if a["accession"] == "P_NOSITE")
    states = {s["name"]: s["state"] for s in nosite["signals"]}
    assert states["conformation"] == "unavailable"
    assert states["assembly"] == "unavailable"


# --- validate ------------------------------------------------------------


def test_validate_reports_the_inputs(examples, capsys):
    assert main(["validate", "--in", str(examples)]) == EXIT_OK
    out = capsys.readouterr().out
    assert "6 protein(s)" in out
    assert "entries with an annotated catalytic site: 5" in out


def test_validate_blocks_when_nothing_is_annotated(tmp_path, capsys):
    path = tmp_path / "annotations.tsv"
    path.write_text("Entry\tft_act_site\nP1\t\n", encoding="utf-8")
    assert main(["validate", "--in", str(tmp_path)]) == EXIT_BLOCKED
    assert "No entry carries an ACT_SITE annotation" in capsys.readouterr().out
