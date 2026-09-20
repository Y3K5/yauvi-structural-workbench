from __future__ import annotations

import hashlib
import importlib.util
import io
from pathlib import Path
import zipfile

import pytest

spec = importlib.util.spec_from_file_location("review_builder", Path(__file__).with_name("build_review_candidate.py"))
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def manifest(name, data):
    return {"files": [{"path": name, "sha256": hashlib.sha256(data).hexdigest()}]}


def test_manifest_binds_exact_bytes_and_ignores_unlisted_files(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "public.txt").write_bytes(b"reviewed")
    (source / "unlisted.txt").write_bytes(b"private")
    report = builder.build(source, manifest("public.txt", b"reviewed"), tmp_path / "out", [])
    assert report["screen_passed"]
    with zipfile.ZipFile(tmp_path / "out" / report["archive"]) as archive:
        assert set(archive.namelist()) == {"public.txt", "REVIEW_MANIFEST.json"}
    (source / "public.txt").write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed"):
        builder.build(source, manifest("public.txt", b"reviewed"), tmp_path / "again", [])


def test_private_content_in_nested_zip_is_held_without_revealing_term(tmp_path):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("report.txt", "confidential-example")
    source = tmp_path / "source"
    source.mkdir()
    data = buffer.getvalue()
    (source / "report.zip").write_bytes(data)
    report = builder.build(source, manifest("report.zip", data), tmp_path / "out", ["confidential-example"])
    assert not report["screen_passed"]
    assert any(x["file"] == "report.zip!report.txt" for x in report["findings"])
    assert "confidential-example" not in str(report)
    assert not (tmp_path / "out" / "yauvi-review-candidate.zip").exists()


@pytest.mark.parametrize("name", ["../escape", "/absolute", "personal/data.json", "a/../../escape", "evidence/benchmarks/qualification-v2/sources/raw.pdb"])
def test_disallowed_paths(name):
    with pytest.raises(ValueError):
        builder.safe_relative(name)


def test_symlink_to_even_an_inside_file_is_rejected(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "real").write_bytes(b"data")
    (source / "alias").symlink_to(source / "real")
    with pytest.raises(ValueError, match="symlinked"):
        builder.build(source, manifest("alias", b"data"), tmp_path / "out", [])


def test_reviewed_fixture_exception_is_visible_and_checksum_bound(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    data = b"/" + b"home/example-user/cache"
    (source / "fixture.txt").write_bytes(data)
    record = manifest("fixture.txt", data)
    record["reviewed_exceptions"] = [{"path": "fixture.txt",
        "sha256": hashlib.sha256(data).hexdigest(), "rule": "absolute_home_path",
        "count": 1, "reason": "Invented regression fixture", "approved_by": "test reviewer"}]
    report = builder.build(source, record, tmp_path / "out", [])
    assert report["screen_passed"]
    assert report["reviewed_findings"][0]["count"] == 1
    changed = data + b" changed"
    (source / "fixture.txt").write_bytes(changed)
    record["files"] = manifest("fixture.txt", changed)["files"]
    with pytest.raises(ValueError, match="stale"):
        builder.build(source, record, tmp_path / "changed", [])


def test_reviewed_fixture_does_not_suppress_private_terms(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    data = b"/" + b"home/example-user/cache confidential-example"
    (source / "fixture.txt").write_bytes(data)
    record = manifest("fixture.txt", data)
    record["reviewed_exceptions"] = [{"path": "fixture.txt",
        "sha256": hashlib.sha256(data).hexdigest(), "rule": "absolute_home_path",
        "count": 1, "reason": "Invented regression fixture", "approved_by": "test reviewer"}]
    report = builder.build(source, record, tmp_path / "out", ["confidential-example"])
    assert not report["screen_passed"]
    assert report["findings"][0]["rule"] == "private_term_1"
    assert not (tmp_path / "out/yauvi-review-candidate.zip").exists()
