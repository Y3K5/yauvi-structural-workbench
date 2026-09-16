from __future__ import annotations

import json
from pathlib import Path

from structqc.cli import main
from structqc.core import analyze, read_validation_report, write_outputs


PDB = """\
ATOM      1  N   ALA A   1      -1.200   0.000   0.000  1.00 91.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 91.00           C
ATOM      3  C   ALA A   1       1.300   0.000   0.000  1.00 91.00           C
ATOM      4  O   ALA A   1       2.100   0.300   0.000  1.00 91.00           O
ATOM      5  N   GLY A   2       1.500   1.200   0.000  1.00 82.00           N
ATOM      6  CA  GLY A   2       2.700   1.500   0.000  1.00 82.00           C
ATOM      7  C   GLY A   2       3.600   0.400   0.000  1.00 82.00           C
ATOM      8  O   GLY A   2       4.800   0.600   0.000  1.00 82.00           O
TER
END
"""


def fixture(tmp_path: Path) -> Path:
    path = tmp_path / "model.pdb"
    path.write_text(PDB)
    return path


def test_unknown_provenance_stays_unknown(tmp_path):
    doc = analyze(fixture(tmp_path))
    assert doc["provenance"]["class"] == "unknown"
    assert doc["completeness"]["state"] == "unevaluated"


def test_predicted_confidence_is_declared(tmp_path):
    doc = analyze(fixture(tmp_path), provenance={"class": "predicted", "method": "AlphaFold"})
    assert [r["plddt"] for r in doc["residues"]] == [91.0, 82.0]


def test_outputs_are_deterministic_and_relative(tmp_path):
    source = fixture(tmp_path)
    doc = analyze(source, provenance={"class": "predicted", "method": "AlphaFold"})
    first, second = tmp_path / "a", tmp_path / "b"
    write_outputs(first, doc)
    write_outputs(second, doc)
    assert {p.name: p.read_bytes() for p in first.iterdir()} == {p.name: p.read_bytes() for p in second.iterdir()}
    assert str(tmp_path) not in (first / "STRUCTURE_EVIDENCE.json").read_text()


def test_cli_contract(tmp_path, capsys):
    assert main(["describe"]) == 0
    assert json.loads(capsys.readouterr().out)["module_id"] == "structure_quality"
    source = fixture(tmp_path)
    provenance = tmp_path / "provenance.json"
    provenance.write_text('{"class":"predicted","method":"AlphaFold"}')
    assert main(["run", "--structure", str(source), "--provenance", str(provenance), "--out", str(tmp_path / "out")]) == 0


def test_external_validation_is_checksum_bound_and_can_be_required(tmp_path):
    source = fixture(tmp_path)
    provenance = tmp_path / "provenance.json"
    provenance.write_text('{"class":"predicted","method":"AlphaFold"}')
    report = tmp_path / "validation.xml"
    report.write_text('<validation clashscore="2.5" rama_outlier_percent="0.2"/>')
    imported = read_validation_report(report)
    assert imported["metrics"]["clashscore"] == 2.5
    assert main(["run", "--structure", str(source), "--provenance", str(provenance),
                 "--validation-report", str(report), "--require-external-validation",
                 "--out", str(tmp_path / "validated")]) == 0
    assert main(["run", "--structure", str(source), "--provenance", str(provenance),
                 "--require-external-validation", "--out", str(tmp_path / "incomplete")]) == 1


def test_wwpdb_validation_attributes_do_not_confuse_percentiles_with_metrics(tmp_path):
    report = tmp_path / "wwpdb-validation.xml"
    report.write_text(
        '<wwPDB-validation-information><Entry PDB-resolution="1.50" clashscore="0.00" '
        'percent-rama-outliers="0.20" percent-rota-outliers="1.30" '
        'absolute-percentile-clashscore="100.0"/></wwPDB-validation-information>'
    )
    metrics = read_validation_report(report)["metrics"]
    assert metrics == {
        "clashscore": 0.0,
        "ramachandran_outliers_percent": 0.2,
        "resolution_angstrom": 1.5,
        "rotamer_outliers_percent": 1.3,
    }


def test_atom_order_does_not_change_residue_quality(tmp_path):
    original = fixture(tmp_path)
    lines = original.read_text().splitlines()
    reordered = tmp_path / "reordered.pdb"
    reordered.write_text("\n".join([*reversed(lines[:4]), *reversed(lines[4:8]), *lines[8:], ""]))
    provenance = {"class": "predicted", "method": "AlphaFold"}
    left = analyze(original, provenance=provenance)
    right = analyze(reordered, provenance=provenance)
    assert left["chain_summaries"] == right["chain_summaries"]
    assert left["residues"] == right["residues"]
