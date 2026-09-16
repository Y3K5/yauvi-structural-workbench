"""The declared interface, and the checks that run without an external binary."""
from __future__ import annotations

import json

import pytest

from sf_csa import __version__
from sf_csa.cli import EXIT_BLOCKED, EXIT_FAILED, EXIT_OK, main
from sf_csa.core import CLASSIFICATION_VOCABULARY


@pytest.fixture
def described(capsys):
    assert main(["describe"]) == EXIT_OK
    return json.loads(capsys.readouterr().out)


def test_describe_names_the_module_and_the_package(described):
    assert described["module_id"] == "sf_csa"
    assert described["package"] == "sf_csa"
    assert described["version"] == __version__


def test_describe_publishes_the_closed_vocabulary(described):
    assert described["classification_vocabulary"] == list(CLASSIFICATION_VOCABULARY)


def test_describe_states_the_separation_rule(described):
    """The module's central design rule must appear in its declared limitations."""
    text = " ".join(described["limitations"])
    assert "must not be merged into a single similarity claim" in text
    assert "never as structural negatives" in text


def test_describe_admits_the_defaults_are_periodontal_biology(described):
    """The old descriptor called this out; the module must keep admitting it."""
    text = " ".join(described["limitations"])
    assert "periodontal-pathogen biology" in text


def test_describe_declares_the_external_runtimes(described):
    assert set(described["runtimes"]) >= {"foldseek", "diamond"}


def test_the_source_manifest_ships_with_the_package():
    from importlib.resources import files

    assert (files("sf_csa") / "sources.yaml").is_file()


def test_fetch_plan_prints_the_declared_sources(capsys):
    assert main(["fetch", "--plan"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "module_id: sf_csa" in out
    assert "pdb" in out


def test_the_package_ships_no_campaign_manifests():
    """Target dictionaries are private campaign material and must not be installed."""
    from importlib.resources import files

    package = files("sf_csa")
    assert not (package / "config").is_dir(), (
        "the package ships a config/ directory; campaign manifests belong to the project"
    )


def test_no_campaign_material_in_the_source():
    """The organisms, strains and antigens live in a project spec, not in code."""
    from importlib.resources import files

    for name in ("core.py", "manifests.py", "cli.py", "module_contract.py"):
        text = (files("sf_csa") / name).read_text(encoding="utf-8").lower()
        for token in ("gingivalis", "forsythia", "denticola", "oralome", "wp_0058", "wp_0142"):
            assert token not in text, f"{name} still names campaign material: {token}"


# --- validate -------------------------------------------------------------


def test_validate_passes_on_built_manifests(campaign, tmp_path, capsys):
    from sf_csa.manifests import build

    target_path, database_path = build(campaign, tmp_path / "out")
    code = main(["validate", "--queries", str(target_path), "--databases", str(database_path)])
    assert code == EXIT_OK
    assert "FAILED" not in capsys.readouterr().out


def test_validate_reports_a_missing_manifest(tmp_path, capsys):
    code = main(
        ["validate", "--queries", str(tmp_path / "a.json"), "--databases", str(tmp_path / "b.json")]
    )
    assert code == EXIT_FAILED
    assert "no such file" in capsys.readouterr().out


def test_validate_reports_malformed_json(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    assert main(["validate", "--queries", str(bad), "--databases", str(bad)]) == EXIT_FAILED
    assert "FAILED" in capsys.readouterr().out


def test_validate_flags_a_query_without_a_mechanism_group(tmp_path, capsys):
    queries = tmp_path / "q.json"
    queries.write_text(
        json.dumps({"queries": [{"accession": "X1"}]}), encoding="utf-8"
    )
    databases = tmp_path / "d.json"
    databases.write_text(json.dumps({"pdb_database": "x"}), encoding="utf-8")
    assert main(["validate", "--queries", str(queries), "--databases", str(databases)]) == EXIT_FAILED
    assert "mechanism_group" in capsys.readouterr().out


def test_validate_says_when_the_default_tables_will_be_used(campaign, tmp_path, capsys):
    """Silently inheriting the wrong biology is the failure this warns about."""
    queries = tmp_path / "q.json"
    queries.write_text(
        json.dumps({"queries": [{"accession": "X1", "mechanism_group": "g"}]}), encoding="utf-8"
    )
    databases = tmp_path / "d.json"
    databases.write_text(
        json.dumps(
            {
                "pdb_database": "x", "pdb_database_checksum": "y",
                "thresholds": {}, "classification_vocabulary": [],
            }
        ),
        encoding="utf-8",
    )
    main(["validate", "--queries", str(queries), "--databases", str(databases)])
    assert "the built-in periodontal-pathogen default will be used" in capsys.readouterr().out


# --- build-manifests ------------------------------------------------------


def test_build_manifests_through_the_cli(campaign, tmp_path, capsys):
    assert main(["build-manifests", "--spec", str(campaign), "--out", str(tmp_path / "cfg")]) == EXIT_OK
    assert (tmp_path / "cfg" / "target_manifest.json").is_file()
    assert "wrote" in capsys.readouterr().out


def test_build_manifests_reports_a_bad_spec(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text("{}", encoding="utf-8")
    assert main(["build-manifests", "--spec", str(bad), "--out", str(tmp_path / "cfg")]) == EXIT_BLOCKED
    assert "campaign spec error" in capsys.readouterr().err


def test_verify_requires_the_database_manifest_explicitly():
    """The release is audited against its manifest, not against its own recorded shape."""
    import argparse

    from sf_csa.cli import build_parser

    with pytest.raises(SystemExit):
        build_parser().parse_args(["verify", "--output", "somewhere"])
