from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from assembly_context import core
from assembly_context.core import analyze, write_outputs


def transform(text: str) -> str:
    output = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            x, y, z = -y + 11.0, x - 7.0, z + 3.0
            line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
        output.append(line)
    return "\n".join(output) + "\n"


def reverse_chain_blocks(text: str) -> str:
    lines = text.splitlines()
    split = lines.index("TER")
    second_end = len(lines) - 2
    return "\n".join([*lines[split + 1:second_end + 1], *lines[:split + 1], "END", ""])


def pdb(chain_b_x=3.0):
    return f"""\
ATOM      1  N   ALA A   1      -1.200   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       1.300   0.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       2.100   0.300   0.000  1.00 20.00           O
TER
ATOM      5  N   GLY B   1       {chain_b_x:5.1f}   0.000   0.000  1.00 20.00           N
ATOM      6  CA  GLY B   1       {chain_b_x+1.2:5.1f}   0.000   0.000  1.00 20.00           C
ATOM      7  C   GLY B   1       {chain_b_x+2.5:5.1f}   0.000   0.000  1.00 20.00           C
ATOM      8  O   GLY B   1       {chain_b_x+3.3:5.1f}   0.300   0.000  1.00 20.00           O
TER
END
"""


def test_contacts_and_burial(tmp_path: Path):
    isolated = tmp_path / "isolated.pdb"; assembly = tmp_path / "assembly.pdb"
    isolated.write_text(pdb(30.0)); assembly.write_text(pdb(3.0))
    digest = hashlib.sha256(isolated.read_bytes()).hexdigest()
    manifest = {"subject": {"id": "Q"}, "coordinate": {"sha256": digest}}
    doc = analyze(manifest, isolated, assembly, subject_chain="A", relationship="exact_protein", expected_chains=["A", "B"])
    assert doc["interfaces"]
    assert doc["surface"]["buried_sasa_A2"] > 0
    assert doc["assembly"]["lower_bound"] is False


def test_partial_assembly_is_a_lower_bound(tmp_path: Path):
    isolated = tmp_path / "isolated.pdb"; assembly = tmp_path / "assembly.pdb"
    isolated.write_text(pdb(30.0)); assembly.write_text(pdb(3.0))
    digest = hashlib.sha256(isolated.read_bytes()).hexdigest()
    doc = analyze({"subject": {"id": "Q"}, "coordinate": {"sha256": digest}}, isolated, assembly,
                  subject_chain="A", relationship="homolog_assembly", expected_chains=["A", "B", "C"])
    assert doc["assembly"]["complete"] is False
    out = tmp_path / "out"; write_outputs(out, doc)
    assert (out / "ASSEMBLY_LAYER.json").is_file()


def test_rigid_transform_and_chain_order_do_not_change_interface_geometry(tmp_path):
    isolated = tmp_path / "isolated.pdb"; assembly = tmp_path / "assembly.pdb"
    isolated.write_text(pdb(30.0)); assembly.write_text(pdb(3.0))
    manifest = {"subject": {"id": "Q"}, "coordinate": {"sha256": hashlib.sha256(isolated.read_bytes()).hexdigest()}}
    original = analyze(manifest, isolated, assembly, subject_chain="A", relationship="exact_protein", expected_chains=["A", "B"])

    moved_isolated = tmp_path / "moved_isolated.pdb"; moved_assembly = tmp_path / "moved_assembly.pdb"
    moved_isolated.write_text(transform(isolated.read_text())); moved_assembly.write_text(transform(assembly.read_text()))
    moved_manifest = {"subject": {"id": "Q"}, "coordinate": {"sha256": hashlib.sha256(moved_isolated.read_bytes()).hexdigest()}}
    moved = analyze(moved_manifest, moved_isolated, moved_assembly, subject_chain="A", relationship="exact_protein", expected_chains=["A", "B"])

    reversed_assembly = tmp_path / "reversed_assembly.pdb"; reversed_assembly.write_text(reverse_chain_blocks(assembly.read_text()))
    reversed_doc = analyze(manifest, isolated, reversed_assembly, subject_chain="A", relationship="exact_protein", expected_chains=["A", "B"])
    assert original["interfaces"] == moved["interfaces"] == reversed_doc["interfaces"]
    assert original["surface"] == moved["surface"] == reversed_doc["surface"]


def _operator_mmcif(path: Path) -> None:
    gemmi = pytest.importorskip("gemmi")
    structure = gemmi.Structure()
    structure.name = "synthetic_operator_fixture"
    structure.cell = gemmi.UnitCell(100, 100, 100, 90, 90, 90)
    structure.spacegroup_hm = "P 1"
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.seqid = gemmi.SeqId(1, " ")
    for serial, (name, element, xyz) in enumerate((
        ("N", "N", (0.0, 0.0, 0.0)),
        ("CA", "C", (1.2, 0.0, 0.0)),
        ("C", "C", (2.4, 0.0, 0.0)),
        ("O", "O", (3.1, 0.4, 0.0)),
    ), 1):
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element(element)
        atom.pos = gemmi.Position(*xyz)
        atom.serial = serial
        residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    structure.assign_subchains()
    structure.assign_label_seq_id()
    subchain = structure[0][0][0].subchain
    assembly = gemmi.Assembly("1")
    generator = gemmi.Assembly.Gen()
    generator.subchains = [subchain]
    for name, offset in (("1", 0.0), ("2", 4.0)):
        operator = gemmi.Assembly.Operator()
        operator.name = name
        transform_value = gemmi.Transform()
        transform_value.vec.x = offset
        operator.transform = transform_value
        generator.operators.append(operator)
    assembly.generators.append(generator)
    structure.assemblies.append(assembly)
    structure.make_mmcif_document().write_file(str(path))


def test_gemmi_applies_declared_mmcif_assembly_operators(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(core.shutil, "which", lambda _name: None)
    isolated = tmp_path / "isolated.pdb"
    isolated.write_text(pdb(30.0))
    assembly = tmp_path / "assembly.cif"
    _operator_mmcif(assembly)
    manifest = {
        "subject": {"id": "Q"},
        "coordinate": {"sha256": hashlib.sha256(isolated.read_bytes()).hexdigest()},
    }
    document = analyze(
        manifest, isolated, assembly, subject_chain="A", relationship="exact_protein",
        assembly_id="1", expected_chains=["A"],
    )
    application = document["assembly"]["metadata"]["operator_application"]
    assert application["operator_backend"] == "gemmi"
    assert application["operator_state"] == "applied"
    assert document["config"]["resolved_subject_chain"] == "A1"
    assert document["assembly"]["chains_observed"] == ["A1", "A2"]
    assert document["interfaces"]


def test_available_freesasa_runtime_is_invoked_and_named(tmp_path: Path, monkeypatch):
    executable = tmp_path / "freesasa"
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import json, pathlib, sys\n"
        "assert '--depth=chain' in sys.argv\n"
        "assert '--output-depth=chain' not in sys.argv\n"
        "isolated = pathlib.Path(sys.argv[-1]).name == 'isolated.pdb'\n"
        "total = 100.0 if isolated else 150.0\n"
        "chains = [{'label':'A','area':{'total':100.0 if isolated else 60.0}}]\n"
        "if not isolated: chains.append({'label':'B','area':{'total':90.0}})\n"
        "print(json.dumps({'results':[{'structure':[{'area':{'total':total},'chains':chains}]}]}))\n"
    )
    executable.chmod(0o755)
    monkeypatch.setattr(core.shutil, "which", lambda name: str(executable) if name == "freesasa" else None)
    isolated = tmp_path / "isolated.pdb"
    assembly = tmp_path / "assembly.pdb"
    isolated.write_text(pdb(30.0))
    assembly.write_text(pdb(3.0))
    manifest = {
        "subject": {"id": "Q"},
        "coordinate": {"sha256": hashlib.sha256(isolated.read_bytes()).hexdigest()},
    }
    document = analyze(
        manifest, isolated, assembly, subject_chain="A", relationship="exact_protein",
        expected_chains=["A", "B"],
    )
    assert document["methods"]["sasa"] == "freesasa_lee_richards_default_single_thread"
    assert document["surface"] == {
        "subject_isolated_sasa_A2": 100.0,
        "subject_assembly_sasa_A2": 60.0,
        "buried_sasa_A2": 40.0,
    }
