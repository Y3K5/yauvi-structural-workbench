"""3D export (viz) + command-line tests."""

from __future__ import annotations

import json
import os

import pytest

from memorient.cli import main
from memorient.contexts import get_context
from memorient.orientor import orient_structure
from memorient.viz import display_oriented, write_3dmol_html, write_pymol_script

from synthetic import make_barrel, make_soluble_blob

GN = get_context("gram_negative_om")
SOL = get_context("soluble_secreted")


def test_display_oriented_has_slab_for_barrel():
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=16, peri_loop_len=2, seed=0)
    r = orient_structure(s, GN, n_points=120, validate=False)
    disp = display_oriented(r)
    # JSON-serializable
    json.dumps(disp)
    assert disp["membrane_slab"] is not None
    assert disp["membrane_slab"]["asymmetric"] is True
    assert "lps_upper_z" in disp["membrane_slab"]        # LPS buffer band present
    assert len(disp["residue_colors"]) == len(s)
    assert "pdb" in disp and disp["pdb"].startswith(("ATOM", "HEADER", "MODEL", "CRYST", "REMARK", "TITLE"))


def test_display_oriented_has_no_slab_for_soluble():
    b = make_soluble_blob(n_res=100, seed=2)
    r = orient_structure(b, SOL, n_points=120, validate=False)
    disp = display_oriented(r)
    assert disp["membrane_slab"] is None                 # no bilayer -> no slab in the viewer
    json.dumps(disp)


def test_write_pymol_script(tmp_path):
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=16, peri_loop_len=2, seed=0)
    r = orient_structure(s, GN, n_points=120, validate=False)
    pml = tmp_path / "view.pml"
    write_pymol_script(r, str(pml))
    text = pml.read_text()
    assert "load oriented.pdb" in text
    assert "color" in text
    assert "pseudoatom mem_ec" in text                   # slab planes drawn for a barrel


def test_write_3dmol_html(tmp_path):
    s = make_barrel(n_strands=12, strand_len=10, ec_loop_len=16, peri_loop_len=2, seed=0)
    r = orient_structure(s, GN, n_points=120, validate=False)
    html = tmp_path / "view.html"
    write_3dmol_html(r, str(html))
    text = html.read_text()
    assert "3Dmol" in text and "addModel" in text
    assert "addBox" in text                              # membrane slab drawn for a barrel
    assert "<!DOCTYPE html>" in text


def test_cli_contexts_runs(capsys):
    rc = main(["contexts"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "gram_negative_om" in out and "soluble_secreted" in out


def test_cli_contexts_json(capsys):
    rc = main(["contexts", "--json"])
    out = capsys.readouterr().out
    assert rc == 0
    data = json.loads(out)
    assert {c["name"] for c in data} >= {"gram_negative_om", "eukaryotic_pm", "soluble_secreted"}


def test_cli_describe_runs(capsys):
    rc = main(["describe", "eukaryotic_pm"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "positive_inside" in out


def test_cli_describe_unknown_context_errors(capsys):
    rc = main(["describe", "not_a_context"])
    assert rc == 2


def test_common_run_contract_writes_deterministic_module_outputs(tmp_path):
    structure = make_soluble_blob(n_res=40, seed=11)
    oriented = orient_structure(structure, SOL, n_points=120, validate=False)
    pdb = tmp_path / "synthetic.pdb"
    pdb.write_text(display_oriented(oriented)["pdb"])
    out = tmp_path / "out"
    assert main(["run", "--structure", str(pdb), "--context", "soluble_secreted", "--out", str(out)]) == 0
    expected = {
        "MEMBRANE_ORIENTATION.json", "ORIENTED_STRUCTURE.pdb", "MEMBRANE_LAYER.json",
        "MEMBRANE_ORIENTATION.pml", "RESIDUE_ORIENTATION.tsv", "RUN_MANIFEST.json",
    }
    assert expected == {path.name for path in out.iterdir()}
    manifest = json.loads((out / "RUN_MANIFEST.json").read_text())
    assert manifest["module_id"] == "membrane_orientation"
    assert "external_opm_ppm_benchmark_adoption" in manifest["missing_evidence"]


@pytest.mark.network
def test_cli_orient_writes_all_outputs(tmp_path):
    import urllib.request

    pdb = tmp_path / "1bxw.pdb"
    urllib.request.urlretrieve("https://files.rcsb.org/download/1BXW.pdb", str(pdb))
    out_json = tmp_path / "r.json"
    out_pdb = tmp_path / "r_oriented.pdb"
    out_viz = tmp_path / "r_viz.json"
    out_pml = tmp_path / "r.pml"
    rc = main([
        "orient", str(pdb), "--context", "gram_negative_om", "--chain", "A",
        "--out-json", str(out_json), "--out-pdb", str(out_pdb),
        "--out-viz", str(out_viz), "--out-pymol", str(out_pml), "--max-rows", "5",
    ])
    assert rc == 0
    for f in (out_json, out_pdb, out_viz, out_pml):
        assert f.exists() and f.stat().st_size > 0
    result = json.loads(out_json.read_text())
    assert result["summary"]["label"] == "barrel"
    assert result["summary"]["host_antibody_accessible"] is True
    assert len(result["surface_set"]) > 0
