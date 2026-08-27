from __future__ import annotations

import hashlib
from pathlib import Path

from site_context.core import analyze, write_outputs


PDB = """\
ATOM      1  N   SER A   1      -1.200   0.000   0.000  1.00 20.00           N
ATOM      2  CA  SER A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   SER A   1       1.300   0.000   0.000  1.00 20.00           C
ATOM      4  O   SER A   1       2.100   0.300   0.000  1.00 20.00           O
ATOM      5  OG  SER A   1       0.000   1.400   0.000  1.00 20.00           O
ATOM      6  N   HIS A   2       1.500   1.200   0.000  1.00 20.00           N
ATOM      7  CA  HIS A   2       2.700   1.500   0.000  1.00 20.00           C
ATOM      8  C   HIS A   2       3.600   0.400   0.000  1.00 20.00           C
ATOM      9  O   HIS A   2       4.800   0.600   0.000  1.00 20.00           O
HETATM   10 ZN    ZN A 101       0.500   2.000   0.000  1.00 20.00          ZN
TER
END
"""


def inputs(tmp_path: Path):
    path = tmp_path / "site.pdb"; path.write_text(PDB)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {"subject": {"id": "ENZ"}, "coordinate": {"sha256": digest}, "residues": [
        {"chain_id": "A", "auth_seq_id": 1, "insertion_code": "", "one_letter": "S", "sequence_index": 1},
        {"chain_id": "A", "auth_seq_id": 2, "insertion_code": "", "one_letter": "H", "sequence_index": 2},
    ]}
    annotations = {"sites": [
        {"position": 1, "type": "active_site", "role": "nucleophile", "expected_residues": ["S", "C"]},
        {"position": 2, "type": "metal_ligand", "role": "metal_ligand", "expected_residues": ["H", "C", "D", "E"]},
    ], "declared_cofactors": [{"name": "zinc", "component_id": "ZN"}]}
    return path, manifest, annotations


def test_roles_cofactor_and_pockets_remain_separate(tmp_path):
    path, manifest, annotations = inputs(tmp_path)
    pocket = {"method": "fpocket", "pockets": [{"pocket_id": "P1", "score": 0.7,
              "residues": [{"chain_id": "A", "auth_seq_id": 1}]}]}
    doc = analyze(manifest, path, annotations, pocket_results=[pocket])
    assert {s["state"] for s in doc["sites"]} == {"role_compatible"}
    assert doc["cofactors"][0]["state"] == "observed_match"
    assert doc["pockets"][0]["method"] == "fpocket"
    out = tmp_path / "out"; write_outputs(out, doc)
    assert (out / "SITE_LAYER.json").is_file()


def test_unknown_cofactor_synonym_is_unresolved(tmp_path):
    path, manifest, annotations = inputs(tmp_path)
    annotations["declared_cofactors"] = [{"name": "zinc ion"}]
    doc = analyze(manifest, path, annotations)
    assert doc["cofactors"][0]["state"] == "unresolved"
