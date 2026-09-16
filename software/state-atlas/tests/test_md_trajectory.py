from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest

from state_atlas.core import analyze


mda = pytest.importorskip("MDAnalysis", reason="install yauvi-state-atlas[md] for trajectory coverage")

ACTIVE = np.asarray([(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2)], dtype=np.float32)
INACTIVE = np.asarray([(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 5)], dtype=np.float32)
UNRESOLVED = (ACTIVE + INACTIVE) / 2


def _pdb(points: np.ndarray) -> str:
    lines = [
        f"ATOM  {index:5d}  CA  ALA A{index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        for index, (x, y, z) in enumerate(points, start=1)
    ]
    return "\n".join([*lines, "END", ""])


def _trajectory(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "topology.pdb"
    topology.write_text(_pdb(ACTIVE), encoding="utf-8")
    universe = mda.Universe(str(topology))
    trajectory = tmp_path / "synthetic.dcd"
    with mda.Writer(str(trajectory), n_atoms=len(universe.atoms)) as writer:
        for coordinates in (ACTIVE, INACTIVE, UNRESOLVED, ACTIVE):
            universe.atoms.positions = coordinates
            writer.write(universe.atoms)
    return topology, trajectory


@pytest.mark.adapter
def test_real_dcd_preserves_active_inactive_and_unresolved_frames(tmp_path):
    topology, trajectory = _trajectory(tmp_path)
    active = tmp_path / "active.pdb"
    inactive = tmp_path / "inactive.pdb"
    active.write_text(_pdb(ACTIVE), encoding="utf-8")
    inactive.write_text(_pdb(INACTIVE), encoding="utf-8")
    references = {
        "reference_set_id": "synthetic-dcd-two-state",
        "decision_rules": {"max_rmsd_A": 0.45, "min_margin_A": 0.15},
        "references": [
            {
                "reference_id": "ACTIVE",
                "state": "active",
                "structure": active.name,
                "provenance": {"class": "experimental", "method": "synthetic coordinate fixture"},
                "state_evidence": {"basis": "compact fourth residue in synthetic fixture"},
            },
            {
                "reference_id": "INACTIVE",
                "state": "inactive",
                "structure": inactive.name,
                "provenance": {"class": "experimental", "method": "synthetic coordinate fixture"},
                "state_evidence": {"basis": "extended fourth residue in synthetic fixture"},
            },
        ],
    }
    manifest = {
        "subject": {"id": "SYNTHETIC_DCD"},
        "coordinate": {"sha256": hashlib.sha256(topology.read_bytes()).hexdigest()},
    }
    document = analyze(
        manifest,
        references,
        reference_base=tmp_path,
        topology_path=topology,
        trajectory_path=trajectory,
        selection="protein and name CA",
        stride=1,
        pbc="none",
        cluster_cutoff_A=0.45,
    )
    calls = [row["call"] for row in document["frame_metrics"]]
    assert calls == ["active_like", "inactive_like", "unresolved", "active_like"]
    assert document["overall_label"] == "mixed"
    assert document["frames_total"] == 4
    assert document["frames_interpretable"] == 3
    assert document["populations"]["unresolved"]["fraction_total"] == 0.25
    assert document["input_kind"] == "trajectory"
    assert document["config"]["pbc"] == "none"
    assert document["config"]["stride"] == 1
