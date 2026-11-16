"""Tests for ``memorient domain-orient`` (memorient.domain_report). Invented coordinates only; no network."""
import json
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.dirname(__file__))
from test_domain_orient import invented_barrel  # noqa: E402
from memorient.cli import main  # noqa: E402
from memorient.domain_report import membrane_frame, zone_for  # noqa: E402

ONE_TO_THREE = dict(L='LEU', S='SER', Y='TYR', G='GLY', D='ASP', N='ASN', T='THR', E='GLU', K='LYS', Q='GLN', A='ALA')


def write_pdb(path, ca, seq, chain='A'):
    lines = [f"ATOM  {i + 1:5d}  CA  {ONE_TO_THREE[a]} {chain}{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
             for i, ((x, y, z), a) in enumerate(zip(ca, seq))]
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def soluble_bundle():
    ca, seq = [], ''
    for h in range(4):
        for i in range(20):
            a = math.radians(100 * i); ca.append([2.3 * math.cos(a) + 10 * (h % 2), 2.3 * math.sin(a) + 10 * (h // 2), 1.5 * i]); seq += 'EKAQ'[i % 4]
    return np.array(ca), seq


def read_z(pdb_path):
    return [float(l[46:54]) for l in pdb_path.read_text().splitlines() if l.startswith('ATOM')]


def test_barrel_writes_four_files_with_membrane_along_z(tmp_path, capsys):
    ca, seq = invented_barrel()
    pdb = write_pdb(tmp_path / 'barrel.pdb', ca, seq)
    assert main(['domain-orient', str(pdb), '--out', str(tmp_path / 'out')]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed['method'].startswith('beta_barrel_core')
    out = tmp_path / 'out'
    assert sorted(p.name for p in out.iterdir()) == ['DOMAIN_ORIENTATION.json', 'MEMBRANE_LAYER.json', 'ORIENTED_STRUCTURE.pdb', 'RESIDUE_ZONES.tsv']
    doc = json.loads((out / 'DOMAIN_ORIENTATION.json').read_text())
    assert doc['scientific_qualification'] == 'not_established' and doc['claim_ceiling']
    assert doc['input']['file_name'] == 'barrel.pdb' and len(doc['input']['sha256']) == 64
    half = doc['result']['half']
    z = np.array(read_z(out / 'ORIENTED_STRUCTURE.pdb'))
    assert len(z) == len(ca)
    # Invented layout: 60-residue soluble domain first, then strands with long loops on the outside.
    loops_and_strands = z[60:]
    assert abs(np.median(loops_and_strands)) < half              # barrel sits across the slab, centred near 0
    assert np.mean(z[:60]) < -half                                # the soluble N-terminal domain is below the membrane (-Z)
    layer = json.loads((out / 'MEMBRANE_LAYER.json').read_text())
    assert layer['membrane_spanning'] is True and layer['label'] == 'barrel'
    assert layer['membrane_slab']['half_thickness'] == half and layer['orientation']['axis'] == '+z'
    assert len(layer['residue_colors']) == len(ca)
    rows = (out / 'RESIDUE_ZONES.tsv').read_text().splitlines()
    assert rows[0].split('\t')[:6] == ['chain', 'resnum', 'insertion_code', 'resname', 'z_A', 'zone'] and len(rows) == len(ca) + 1
    assert sum(doc['zone_counts'].values()) == len(ca)


def test_outputs_are_byte_identical_on_rerun(tmp_path):
    ca, seq = invented_barrel()
    pdb = write_pdb(tmp_path / 'barrel.pdb', ca, seq)
    assert main(['domain-orient', str(pdb), '--out', str(tmp_path / 'a')]) == 0
    assert main(['domain-orient', str(pdb), '--out', str(tmp_path / 'b')]) == 0
    for name in ('DOMAIN_ORIENTATION.json', 'MEMBRANE_LAYER.json', 'ORIENTED_STRUCTURE.pdb', 'RESIDUE_ZONES.tsv'):
        assert (tmp_path / 'a' / name).read_bytes() == (tmp_path / 'b' / name).read_bytes(), name


def test_soluble_structure_reports_no_core_and_writes_no_pose(tmp_path, capsys):
    ca, seq = soluble_bundle()
    pdb = write_pdb(tmp_path / 'soluble.pdb', ca, seq)
    assert main(['domain-orient', str(pdb), '--out', str(tmp_path / 'out')]) == 1
    out = tmp_path / 'out'
    assert sorted(p.name for p in out.iterdir()) == ['DOMAIN_ORIENTATION.json']
    doc = json.loads((out / 'DOMAIN_ORIENTATION.json').read_text())
    assert doc['result']['method'] == 'none' and any('No membrane core' in w for w in doc['warnings'])


def test_input_errors_exit_2_without_writing(tmp_path, capsys):
    ca, seq = invented_barrel()
    pdb = write_pdb(tmp_path / 'barrel.pdb', ca, seq)
    assert main(['domain-orient', str(pdb), '--chain', 'Z', '--out', str(tmp_path / 'none')]) == 2
    assert 'No usable C-alpha atoms' in capsys.readouterr().err and not (tmp_path / 'none').exists()
    (tmp_path / 'exists').mkdir()
    assert main(['domain-orient', str(pdb), '--out', str(tmp_path / 'exists')]) == 2


@pytest.mark.parametrize('normal', [[0, 0, 1], [1, 0, 0], [0.3, -0.5, 0.81], [-0.95, 0.1, 0.2]])
def test_membrane_frame_is_a_rotation_putting_the_normal_on_z(normal):
    R = membrane_frame(normal)
    n = np.asarray(normal, float) / np.linalg.norm(normal)
    assert np.allclose(R @ R.T, np.eye(3)) and math.isclose(np.linalg.det(R), 1.0)
    assert np.allclose(R @ n, [0, 0, 1])


def test_zone_boundaries():
    assert zone_for(0, 12) == zone_for(12, 12) == zone_for(-12, 12) == 'membrane_core'
    assert zone_for(16.9, 12) == zone_for(-17, 12) == 'interface'
    assert zone_for(17.1, 12) == 'positive_side' and zone_for(-17.1, 12) == 'negative_side'
