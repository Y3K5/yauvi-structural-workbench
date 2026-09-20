"""Focused fail-closed checks for the public structural case runner."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

SCRIPT = Path(__file__).with_name("public_research_case.py")
SPEC = importlib.util.spec_from_file_location("public_research_case", SCRIPT)
assert SPEC and SPEC.loader
case = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(case)
ACQUIRE_SCRIPT = Path(__file__).with_name("acquire_public_inputs.py")
ACQUIRE_SPEC = importlib.util.spec_from_file_location("acquire_public_inputs", ACQUIRE_SCRIPT)
assert ACQUIRE_SPEC and ACQUIRE_SPEC.loader
acquire = importlib.util.module_from_spec(ACQUIRE_SPEC)
ACQUIRE_SPEC.loader.exec_module(acquire)


class PublicResearchCaseSafetyTests(unittest.TestCase):
    def test_finds_enclosing_workbench_source_tree_and_refuses_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            repo = Path(temporary) / "checkout"
            script_root = repo / "examples" / "structural-portfolio"
            (repo / "software").mkdir(parents=True)
            (repo / "pyproject.toml").write_text("[project]\nname='test'\n")
            script_root.mkdir(parents=True)
            self.assertEqual(case.find_source_root(script_root), repo.resolve())
            with self.assertRaisesRegex(ValueError, "outside the source tree"):
                case.assert_outside_source(repo / "software" / "analysis", repo, "output directory")
            with self.assertRaisesRegex(ValueError, "outside the source tree"):
                acquire.assert_outside_source(repo / "software" / "inputs", repo)
            case.assert_outside_source(Path(temporary) / "run", repo, "output directory")

    def test_standalone_copy_uses_case_directory_as_fallback_boundary(self):
        with tempfile.TemporaryDirectory() as temporary:
            standalone = Path(temporary) / "structural-portfolio"
            standalone.mkdir()
            self.assertEqual(case.find_source_root(standalone), standalone.resolve())
            with self.assertRaisesRegex(ValueError, "outside the source tree"):
                case.assert_outside_source(standalone / "results", standalone, "output directory")
            with self.assertRaisesRegex(ValueError, "outside the source tree"):
                acquire.assert_outside_source(standalone / "inputs", standalone)

    def test_source_lock_rejects_tampered_and_missing_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "coordinates.cif"
            path.write_bytes(b"locked bytes")
            source = {"path": path.name, "sha256": case.sha256(path), "byte_count": path.stat().st_size}
            self.assertEqual(case.verify_input(path, source)["sha256"], source["sha256"])
            path.write_bytes(b"tampered bytes")
            with self.assertRaisesRegex(RuntimeError, "source-lock mismatch"):
                case.verify_input(path, source)
            path.unlink()
            with self.assertRaisesRegex(RuntimeError, "required input missing"):
                case.verify_input(path, source)

    def test_missing_selected_chain_and_model_fail_closed(self):
        models = [[SimpleNamespace(id="A"), SimpleNamespace(id="B")]]
        case.validate_model_chain("TEST", models, 0, "B")
        with self.assertRaisesRegex(RuntimeError, "selected chain 'Z' not found"):
            case.validate_model_chain("TEST", models, 0, "Z")
        with self.assertRaisesRegex(RuntimeError, "model 1 not found"):
            case.validate_model_chain("TEST", models, 1, "A")

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            existing = Path(temporary) / "previous-run"
            existing.mkdir()
            marker = existing / "keep.txt"
            marker.write_text("preserve")
            with self.assertRaisesRegex(ValueError, "already exists"):
                case.require_new_output(existing)
            self.assertEqual(marker.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main()
