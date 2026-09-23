"""The shared content screen catches what it must and never repeats a private term."""
from __future__ import annotations

import io
import json
from pathlib import Path
import sys
import zipfile

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import content_screen as cs  # noqa: E402

TERM = "Example Project 7"
# Built at run time so this file never contains a literal non-placeholder home path,
# which the repository-wide path verifier would rightly reject.
HOME = "/" + "Users/" + "alice"


@pytest.fixture(autouse=True)
def isolated(request, monkeypatch, tmp_path):
    """No real denylist or real home directory leaks into these tests."""
    if request.node.name == "test_public_tools_carry_no_term_from_a_local_denylist":
        return
    monkeypatch.delenv(cs.ENV_DENYLIST, raising=False)
    monkeypatch.setattr(cs, "DEFAULT_DENYLIST", tmp_path / "absent.json")
    monkeypatch.setattr(cs, "runtime_terms", lambda: [])


def categories(files, denylist=()):
    return {(f["file"], f["category"]) for f in cs.screen(files, list(denylist))["findings"]}


def test_private_term_matches_across_case_and_separators_without_being_reported():
    for spelling in ("example project 7", "EXAMPLE_PROJECT_7", "Example-Project-7", "example%20project%207"):
        report = cs.screen({"notes.txt": f"see {spelling} here".encode()}, [TERM])
        assert [f["category"] for f in report["findings"]] == ["private_term_1"]
        assert "project" not in json.dumps(report).lower()


def test_home_paths_are_caught_but_placeholders_and_bare_prefixes_are_not():
    assert categories({"a.json": ('{"p": "' + HOME + '/work/x.pdb"}').encode()}) == {("a.json", "absolute_local_path")}
    assert categories({"w.txt": ("C:" + "\\Users\\" + "alice\\data").encode()}) == {("w.txt", "absolute_local_path")}
    assert categories({"t.txt": ("/private" + "/tmp/run-1/out").encode()}) == {("t.txt", "absolute_local_path")}
    assert categories({"ok.txt": b"/Users/runner/work and /home/researcher/x"}) == set()
    assert categories({"test.py": b"assert '/Users/' not in text"}) == set()


def test_credentials_nested_archives_and_file_names_are_screened():
    inner = io.BytesIO()
    with zipfile.ZipFile(inner, "w") as z:
        z.writestr("pkg/secret.txt", "token ghp_" + "a" * 36)
    outer = io.BytesIO()
    with zipfile.ZipFile(outer, "w") as z:
        z.writestr("app.whl", inner.getvalue())
    found = categories({"kit.zip": outer.getvalue(), "example_project_7.txt": b"clean"}, [TERM])
    assert ("kit.zip!app.whl!pkg/secret.txt", "credential_shape") in found
    assert ("example_project_7.txt", "private_term_1") in found


def test_unsafe_members_and_opaque_archives_are_held_for_review():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        z.writestr("../escape.txt", "x")
    found = categories({"bad.zip": buffer.getvalue(), "data.tar.gz": b"\x1f\x8b"})
    assert ("bad.zip!../escape.txt", "unsafe_archive_member") in found
    assert ("data.tar.gz", "unsupported_archive_requires_review") in found


def test_missing_denylist_fails_closed_unless_opted_out():
    with pytest.raises(cs.DenylistError, match="no private denylist"):
        cs.load_denylist()
    assert cs.load_denylist(allow_missing=True) == ([], "not used (explicit opt-out)")


def test_denylist_sources_are_validated(tmp_path, monkeypatch):
    good = tmp_path / "deny.json"
    good.write_text(json.dumps([TERM]))
    assert cs.load_denylist(good)[0] == [TERM]
    monkeypatch.setenv(cs.ENV_DENYLIST, str(good))
    assert cs.load_denylist()[0] == [TERM]
    monkeypatch.setenv(cs.ENV_DENYLIST, str(tmp_path / "missing.json"))
    with pytest.raises(cs.DenylistError, match="does not exist"):
        cs.load_denylist()
    monkeypatch.delenv(cs.ENV_DENYLIST)
    for bad in ([], [""], ["  "], "not a list", [3]):
        good.write_text(json.dumps(bad))
        with pytest.raises(cs.DenylistError):
            cs.load_denylist(good)


def test_denylist_inside_the_repository_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(cs, "ROOT", tmp_path)
    inside = tmp_path / "tools" / "deny.json"
    inside.parent.mkdir()
    inside.write_text(json.dumps([TERM]))
    with pytest.raises(cs.DenylistError, match="outside the repository"):
        cs.load_denylist(inside)


def test_runtime_terms_skip_placeholder_accounts(monkeypatch):
    monkeypatch.undo()
    monkeypatch.setattr(cs.Path, "home", staticmethod(lambda: Path("/Users/runner")))
    monkeypatch.setattr(cs.getpass, "getuser", lambda: "runner")
    assert cs.runtime_terms() == []
    monkeypatch.setattr(cs.Path, "home", staticmethod(lambda: Path(HOME + "-example")))
    monkeypatch.setattr(cs.getpass, "getuser", lambda: "alice-example")
    assert cs.runtime_terms() == [HOME + "-example", "alice-example"]


def test_public_tools_carry_no_term_from_a_local_denylist():
    """Runs only where a maintainer's denylist exists; CI has none and skips."""
    try:
        terms, _ = cs.load_denylist()
    except cs.DenylistError:
        pytest.skip("no local private denylist on this machine")
    tools = Path(__file__).resolve().parent
    files = {p.name: p.read_bytes() for p in list(tools.glob("*.py")) + list(tools.glob("*kit/*")) if p.is_file()}
    hits = [f for f in cs.screen(files, terms)["findings"] if f["category"].startswith("private_term_")]
    assert hits == []
