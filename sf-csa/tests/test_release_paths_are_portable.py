"""Release evidence must not record the generating machine's filesystem layout.

`proteome_denominator.json` is part of a release and is published. It recorded
`str(path)` after `Path.resolve()` for every proteome FASTA, so it carried the
absolute checkout location — on 2026-09-03, ten occurrences of a home directory
in a file staged for a public repository. Three tracked `EXECUTION_STATUS.json`
files in the same repository already carried the same class through a different
writer, which is the point: fixing one writer does not close the class, so this
holds the property rather than the incident.

`stage_campaign` already recorded sources relative to the database root; this
brings `build_proteome_universe` to the same convention. The sha256 beside each
path is what identifies the file, so nothing identifying is lost.

Written to fail against the previous implementation: it asserted no absolute
path, which the first fix satisfied while silently recording only a basename,
because the local name `root` is rebound by the glob loop above. The
directory-component assertion is what distinguishes the two.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from sf_csa import core


@pytest.fixture()
def universe(tmp_path: Path):
    """A database root with two proteomes, laid out as a release expects."""
    root = tmp_path / "panel"
    proteomes = root / "sources" / "proteomes"
    proteomes.mkdir(parents=True)
    (proteomes / "UP000000001.fasta").write_text(">a|one desc\nMKV\n", encoding="utf-8")
    (proteomes / "UP000000002.fasta").write_text(">b|two desc\nMKW\n", encoding="utf-8")

    manifest_dir = root / "config" / "sf_csa"
    manifest_dir.mkdir(parents=True)
    manifest_path = manifest_dir / "database_manifest.json"
    db = {"path_base": "../..", "proteome_globs": ["sources/proteomes/*.fasta"]}
    manifest_path.write_text(json.dumps(db), encoding="utf-8")
    return db, manifest_path, root, tmp_path


def test_proteome_paths_are_relative_to_the_database_root(universe, tmp_path):
    db, manifest_path, root, _ = universe
    records, _index = core.build_proteome_universe(db, manifest_path, tmp_path / "u.fasta")

    assert len(records) == 2
    for record in records:
        path = record["path"]
        assert not Path(path).is_absolute(), f"absolute path in release evidence: {path}"
        assert str(root) not in path, f"database root leaked into release evidence: {path}"
        # The whole point of "relative to the root" rather than "basename": the
        # location inside the panel is still readable by a reviewer.
        assert path.startswith("sources/proteomes/"), (
            f"expected a root-relative path, got {path!r} — a bare basename means "
            "the relative branch did not fire"
        )
        assert record["sha256"], "the digest identifies the file and must survive"


def test_no_absolute_path_survives_anywhere_in_the_denominator(universe, tmp_path):
    """The published document as a whole, not just the field we remembered."""
    db, manifest_path, root, _ = universe
    records, _index = core.build_proteome_universe(db, manifest_path, tmp_path / "u.fasta")
    document = {
        "proteome_count": len(records),
        "protein_count": sum(r["protein_count"] for r in records),
        "proteomes": records,
    }
    serialized = json.dumps(document)
    assert str(root) not in serialized
    assert str(Path.home()) not in serialized


def test_a_proteome_outside_the_root_degrades_to_its_name(tmp_path):
    """Outside the root there is no relative form; a basename still leaks nothing."""
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    (outside / "UP000000003.fasta").write_text(">c|three d\nMKY\n", encoding="utf-8")

    root = tmp_path / "panel"
    manifest_dir = root / "config" / "sf_csa"
    manifest_dir.mkdir(parents=True)
    manifest_path = manifest_dir / "database_manifest.json"
    db = {"path_base": "../..", "proteome_globs": [f"{outside}/*.fasta"]}
    manifest_path.write_text(json.dumps(db), encoding="utf-8")

    records, _index = core.build_proteome_universe(db, manifest_path, tmp_path / "u.fasta")
    assert [r["path"] for r in records] == ["UP000000003.fasta"]
    assert not Path(records[0]["path"]).is_absolute()
