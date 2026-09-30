"""The bundled engines must not require nonexistent standalone metadata."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

spec = importlib.util.spec_from_file_location("distribution_runner", Path(__file__).with_name("run_structural_workbench_tests.py"))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


@pytest.mark.parametrize("name,version", [("memorient", "0.3.0"), ("yauvi-structqc", "0.1.0"), ("yauvi-assembly-context", "0.1.0")])
def test_engines_from_source_tree_satisfy_metadata_floor(monkeypatch, name, version):
    monkeypatch.setattr(runner.importlib.util, "find_spec", lambda module: SimpleNamespace(
        origin=str(next((runner.ROOT / "software").rglob(f"src/{module}/__init__.py")))))
    assert runner.bundled_version(name) == version


def test_external_dependency_cannot_be_satisfied_by_bundled_fallback():
    assert runner.bundled_version("numpy") is None


def test_unowned_engine_with_same_module_name_is_rejected(monkeypatch, tmp_path):
    monkeypatch.setattr(runner.importlib.util, "find_spec", lambda _: SimpleNamespace(origin=str(tmp_path / "__init__.py")))
    monkeypatch.setattr(runner.importlib_metadata, "distributions", lambda **_: [SimpleNamespace(files=[])])
    assert runner.bundled_version("memorient") is None


def test_installed_bundle_can_own_engine(monkeypatch, tmp_path):
    origin = tmp_path / "__init__.py"
    monkeypatch.setattr(runner.importlib.util, "find_spec", lambda _: SimpleNamespace(origin=str(origin)))
    monkeypatch.setattr(runner.importlib_metadata, "distributions", lambda **_: [
        SimpleNamespace(files=[]), SimpleNamespace(
            files=["memorient/__init__.py"], locate_file=lambda _: origin)])
    assert runner.bundled_version("memorient") == "0.3.0"
