"""Exit codes and output, because these commands are run from scripts."""
from __future__ import annotations

import json

import pytest

from yauvi_sources.cli import EXIT_OK, EXIT_UNSATISFIED, EXIT_USAGE, find_registry, main

PAYLOAD = b">p\nMKVLAAG\n"


@pytest.fixture
def cli(tmp_path, registry_path, monkeypatch):
    """Run the CLI against the synthetic registry and an isolated cache."""
    manifest_path = tmp_path / "mod_sources.yaml"
    manifest_path.write_text(
        "module_id: testmod\n"
        "requires:\n"
        "  - source_id: open_db\n"
        "    required: true\n"
        "  - source_id: gated_db\n"
        "    required: false\n",
        encoding="utf-8",
    )
    cache_dir = tmp_path / "cache"

    def run(*args: str) -> int:
        return main(
            ["--registry", str(registry_path), "--cache", str(cache_dir), *args]
        )

    run.manifest = str(manifest_path)
    run.cache_dir = cache_dir
    return run


def test_plan_exits_nonzero_when_a_required_source_is_absent(cli, capsys):
    code = cli("plan", "--for", "testmod", "--manifest", cli.manifest)
    assert code == EXIT_UNSATISFIED
    assert "NOT SATISFIED" in capsys.readouterr().out


def test_plan_exits_zero_once_the_requirement_is_staged(cli, capsys, tmp_path):
    staged = tmp_path / "db.fasta"
    staged.write_bytes(PAYLOAD)
    assert cli("stage", "open_db", str(staged)) == EXIT_OK
    capsys.readouterr()
    assert cli("plan", "--for", "testmod", "--manifest", cli.manifest) == EXIT_OK
    assert "SATISFIED" in capsys.readouterr().out


def test_plan_json_is_machine_readable(cli, capsys):
    cli("plan", "--for", "testmod", "--manifest", cli.manifest, "--json")
    payload = json.loads(capsys.readouterr().out)
    assert payload["module_id"] == "testmod"
    assert payload["satisfied"] is False
    ids = {item["source_id"]: item for item in payload["items"]}
    assert ids["open_db"]["blocks"] is True
    assert ids["gated_db"]["fetch_class"] == "license_gated"


def test_staging_an_undeclared_source_is_refused(cli, tmp_path, capsys):
    path = tmp_path / "x.faa"
    path.write_bytes(PAYLOAD)
    assert cli("stage", "not_declared", str(path)) == EXIT_USAGE
    assert "not declared" in capsys.readouterr().err


def test_where_reports_a_missing_source_on_stderr(cli, capsys):
    assert cli("where", "open_db") == EXIT_UNSATISFIED
    assert "not cached" in capsys.readouterr().err


def test_where_prints_the_path_once_staged(cli, tmp_path, capsys):
    staged = tmp_path / "db.fasta"
    staged.write_bytes(PAYLOAD)
    cli("stage", "open_db", str(staged))
    capsys.readouterr()
    assert cli("where", "open_db") == EXIT_OK
    assert capsys.readouterr().out.strip().endswith("db.fasta")


def test_verify_fails_after_the_cache_is_tampered_with(cli, tmp_path, capsys):
    staged = tmp_path / "db.fasta"
    staged.write_bytes(PAYLOAD)
    cli("stage", "open_db", str(staged))
    capsys.readouterr()
    cached = next(cli.cache_dir.rglob("db.fasta"))
    cached.write_bytes(PAYLOAD + b"tampered")
    assert cli("verify") == EXIT_UNSATISFIED
    assert "FAIL" in capsys.readouterr().out


def test_verify_on_an_empty_cache_is_not_a_failure(cli, capsys):
    assert cli("verify") == EXIT_OK
    assert "cache is empty" in capsys.readouterr().out


def test_sources_lists_the_registry(cli, capsys):
    assert cli("sources") == EXIT_OK
    out = capsys.readouterr().out
    assert "open_db" in out and "6 source(s)" in out


def test_sources_json_carries_the_fetch_class(cli, capsys):
    cli("sources", "--json")
    rows = {r["source_id"]: r for r in json.loads(capsys.readouterr().out)}
    assert rows["a_binary"]["fetch_class"] == "runtime"
    assert rows["reachable_over_api"]["fetch_class"] == "table_only"


def test_unknown_module_is_a_usage_error(cli, capsys):
    assert cli("plan", "--for", "no_such_module") == EXIT_USAGE
    assert "no source manifest found" in capsys.readouterr().err


def test_missing_registry_is_a_usage_error(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("YAUVI_SOURCES_REGISTRY", raising=False)
    monkeypatch.chdir(tmp_path)
    assert main(["sources"]) == EXIT_USAGE
    assert "could not locate" in capsys.readouterr().err


def test_find_registry_walks_upward(tmp_path, monkeypatch):
    (tmp_path / "catalogs").mkdir()
    registry = tmp_path / "catalogs" / "sources.yaml"
    registry.write_text("sources: []", encoding="utf-8")
    nested = tmp_path / "a" / "b"
    nested.mkdir(parents=True)
    monkeypatch.delenv("YAUVI_SOURCES_REGISTRY", raising=False)
    assert find_registry(None, start=nested) == registry
