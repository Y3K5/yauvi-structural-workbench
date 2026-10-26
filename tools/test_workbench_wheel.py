"""Exercise actual wheel admission and rejection with independent package fixtures."""
import importlib.util
from pathlib import Path
import zipfile

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 uses the already declared dependency.
    import tomli as tomllib


TOOLS = Path(__file__).resolve().parent
# Build test inputs from the declared distribution, never from the guard's list.
# The guard keeps its independent, hand-reviewed list; disagreement must fail.
CONFIG = tomllib.loads((TOOLS.parent / "pyproject.toml").read_text())
PACKAGE_ROOTS = {
    pattern.rstrip("*").split(".")[0]
    for pattern in CONFIG["tool"]["setuptools"]["packages"]["find"]["include"]
}


@pytest.fixture
def guard():
    spec = importlib.util.spec_from_file_location(
        "workbench_wheel_guard", TOOLS / "verify_structural_workbench_wheel.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_wheel(tmp_path, *, omit=(), extra=(), licenses=None):
    path = tmp_path / "example.whl"
    if licenses is None:
        licenses = ["example.dist-info/LICENSE"]
    with zipfile.ZipFile(path, "w") as archive:
        for root in sorted(PACKAGE_ROOTS - set(omit)):
            package = root + ("/structural_workbench" if root == "yauvi_platform" else "")
            archive.writestr(package + "/__init__.py", "")
        for name in extra:
            archive.writestr(name, "invented control")
        for name in licenses:
            archive.writestr(name, "invented license")
    return path


def test_current_distribution_is_accepted(guard, tmp_path):
    assert guard.main(["guard", str(make_wheel(tmp_path))]) == 0


@pytest.mark.parametrize("root", sorted(PACKAGE_ROOTS))
def test_every_required_package_must_be_present(guard, tmp_path, capsys, root):
    wheel = make_wheel(tmp_path, omit=[root])
    assert guard.main(["guard", str(wheel)]) == 1
    error = capsys.readouterr().err
    assert f"'missing': ['{root}']" in error


@pytest.mark.parametrize("root", ["accidental_module", "md_health_extra"])
def test_unlisted_package_still_fails(guard, tmp_path, capsys, root):
    wheel = make_wheel(tmp_path, extra=[root + "/__init__.py"])
    assert guard.main(["guard", str(wheel)]) == 1
    assert f"'unexpected': ['{root}']" in capsys.readouterr().err


@pytest.mark.parametrize("name", [
    "md_health/_archive/old.py",
    "md_health/oral_atlas/private.py",
    "md_health/subproteo/private.py",
])
def test_admitted_modules_cannot_hide_excluded_material(guard, tmp_path, capsys, name):
    wheel = make_wheel(tmp_path, extra=[name])
    assert guard.main(["guard", str(wheel)]) == 1
    assert "'forbidden': []" not in capsys.readouterr().err


@pytest.mark.parametrize("name", [
    "yauvi_platform/__init__.py",
    "yauvi_platform/unrelated/__init__.py",
])
def test_platform_namespace_remains_scoped(guard, tmp_path, capsys, name):
    wheel = make_wheel(tmp_path, extra=[name])
    assert guard.main(["guard", str(wheel)]) == 1
    assert "outside structural_workbench" in capsys.readouterr().err


def test_new_setuptools_license_layout_is_accepted(guard, tmp_path):
    wheel = make_wheel(tmp_path, licenses=[
        "example.dist-info/licenses/yauvi-structural-workbench/LICENSE"
    ])
    assert guard.main(["guard", str(wheel)]) == 0


@pytest.mark.parametrize("licenses", [[], [
    "example.dist-info/LICENSE",
    "example.dist-info/licenses/yauvi-structural-workbench/LICENSE",
]])
def test_exactly_one_license_is_required(guard, tmp_path, capsys, licenses):
    wheel = make_wheel(tmp_path, licenses=licenses)
    assert guard.main(["guard", str(wheel)]) == 1
    assert "exactly one complete project LICENSE" in capsys.readouterr().err
