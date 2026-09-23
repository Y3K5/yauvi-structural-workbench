"""Geometry controls: rigid motion, reflection, missing atoms and source drift."""
import copy
import io
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

_spec = importlib.util.spec_from_file_location("showcase_state_motion", Path(__file__).with_name("showcase_state_motion.py"))
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
angles, build_motion, fit = _module.angles, _module.build_motion, _module.fit


class GeometryControls(unittest.TestCase):
    def test_rigid_transform_is_not_conformational_change(self):
        fixed = np.array([[0., 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4], [2, 3, 5]])
        rotated = fixed @ np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]]) + [20, -3, 15]
        r, t = fit(rotated, fixed)
        np.testing.assert_allclose(rotated @ r + t, fixed, atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(r), 1.)

    def test_mirror_image_cannot_be_fit_as_a_rotation(self):
        fixed = np.array([[0., 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4], [2, 3, 5]])
        mirrored = fixed * [-1, 1, 1]
        r, t = fit(mirrored, fixed)
        self.assertAlmostEqual(np.linalg.det(r), 1.)
        self.assertGreater(np.linalg.norm(mirrored @ r + t - fixed), .5)


class AcquiredPairControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel = Path(__file__).resolve().parents[1] / "evidence/benchmarks/qualification-v2"
        if not (cls.panel / "sources/wwpdb/2GQG.cif").exists():
            raise unittest.SkipTest("Requires acquired, locked ABL sources")
        cls.motion = build_motion(cls.panel)

    def test_missing_and_modified_residues_are_not_filled(self):
        self.assertEqual(self.motion["excluded_positions"], [275, 393])
        self.assertEqual(self.motion["mapped_count"], 252)
        data = {r["position"]: r for r in self.motion["residues"]}
        self.assertIsNone(data[392]["active_angles"]["psi"])
        self.assertIsNone(data[394]["active_angles"]["phi"])

    def test_exported_coordinates_reproduce_displacement(self):
        parser = PDBParser(QUIET=True)
        a = parser.get_structure("a", io.StringIO(self.motion["active_pdb"]))[0]["A"]
        b = parser.get_structure("b", io.StringIO(self.motion["inactive_pdb"]))[0]["A"]
        for row in self.motion["residues"]:
            pos = row["position"]
            self.assertAlmostEqual(float(np.linalg.norm(a[pos]["CA"].coord-b[pos]["CA"].coord)), row["displacement"], delta=.002)
        # Torsions must not change when coordinates undergo a rigid transform.
        original = {r.id[1]: r for r in a}
        moved = copy.deepcopy(original)
        rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        for residue in moved.values():
            for atom in residue:
                atom.coord = atom.coord @ rotation + [11, -3, 8]
        self.assertEqual(angles(original, 382), angles(moved, 382))

    def test_changed_source_fails_before_geometry(self):
        # Minimal copy sufficient to reach the first digest check; no frozen file changed.
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            run = Path("results/execution-abl/v2-stateatlas-active-held_out-2GQG/RUN_MANIFEST.json")
            (root / run).parent.mkdir(parents=True)
            (root / run).write_bytes((self.panel / run).read_bytes())
            source = root / "sources/wwpdb/2GQG.cif"
            source.parent.mkdir(parents=True)
            source.write_text("changed coordinates")
            with self.assertRaisesRegex(ValueError, "input drift"):
                build_motion(root)


if __name__ == "__main__":
    unittest.main()
