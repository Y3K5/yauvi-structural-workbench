"""Tests for the SF-CSA offline fixture.

Run with:

    PY=$PWD/.venv/bin/python ../../.venv/bin/python -m pytest tools/fixtures/sfcsa/ -q

Three things are under test here, and they are not the same thing:

1.  The stubs themselves. `foldseek` and `diamond` are hand-written scripts
    that emit canned tables. If a stub silently stops emitting a hit, the
    pipeline still runs and still produces a clean release -- just a smaller
    one. So the stubs are tested directly, not only through the pipeline.
2.  The input generator. `build_inputs.py --check` proves the generated inputs
    are reproducible; without that, a fixture failure is ambiguous between
    "the code changed" and "the inputs drifted".
3.  The two scenarios, end to end, against the golden trees -- including that
    the trap scenario's audit *fails*. A fixture whose audit always passes
    does not demonstrate that the audit works, so "trap must fail" is an
    assertion, not an observation.

Everything is offline and stdlib + pytest only. `sf_csa` must be importable;
if it is not, the pipeline tests skip rather than error, because the stub and
generator tests are still meaningful on their own.
"""
from __future__ import annotations

import filecmp
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
STUB_BIN = HERE / "stub_bin"
INPUTS = HERE / "inputs"
GOLDEN = HERE / "golden"
HITS = json.loads((STUB_BIN / "hits.json").read_text(encoding="utf-8"))

sf_csa = pytest.importorskip("sf_csa", reason="sf-csa is not importable in this environment")


def stub_env(scenario: str) -> dict:
    """Environment with the stubs ahead of any real tools on PATH."""
    env = dict(os.environ)
    env["PATH"] = f"{STUB_BIN}{os.pathsep}{env.get('PATH', '')}"
    env["SFCSA_FIXTURE_SCENARIO"] = scenario
    return env


def run(args: list[str], scenario: str = "main", cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        args, env=stub_env(scenario), cwd=cwd, text=True, capture_output=True
    )


# --------------------------------------------------------------------------
# The stubs
# --------------------------------------------------------------------------

def test_the_stubs_are_executable():
    """A stub that is not +x is invisible to shutil.which, and the pipeline
    would fall through to a real foldseek if one happened to be installed."""
    for name in ("foldseek", "diamond"):
        path = STUB_BIN / name
        assert path.exists(), f"{name} stub missing"
        assert os.access(path, os.X_OK), f"{name} stub is not executable"


def test_the_stubs_shadow_any_real_tool():
    """The fixture is only offline if the stubs win the PATH lookup."""
    env = stub_env("main")
    for name in ("foldseek", "diamond"):
        found = shutil.which(name, path=env["PATH"])
        assert found is not None
        assert Path(found).resolve() == (STUB_BIN / name).resolve()


def test_the_foldseek_stub_reports_the_version_the_manifest_requires():
    """`run_pipeline` refuses to run when the reported version does not match
    `required_foldseek_version`, so this pairing is part of the contract."""
    required = json.loads((INPUTS / "database_manifest.json").read_text(encoding="utf-8"))[
        "required_foldseek_version"
    ]
    proc = run([str(STUB_BIN / "foldseek"), "version"])
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == required


def test_the_diamond_stub_reports_a_version():
    proc = run([str(STUB_BIN / "diamond"), "version"])
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == HITS["diamond_version"]


def test_the_foldseek_stub_emits_the_fields_it_was_asked_for(tmp_path):
    """The pipeline passes --format-output and then reads the columns
    positionally. A stub that emits its own column order would put TM scores in
    the coverage column and every classification would be wrong-but-plausible."""
    out = tmp_path / "hits.tsv"
    fields = HITS["fields"]
    proc = run(
        [
            str(STUB_BIN / "foldseek"), "easy-search",
            "QRY_A.pdb", "some/experimental_pdb", str(out), str(tmp_path / "tmp"),
            "--format-output", fields, "-e", "1e-3", "--max-seqs", "50",
        ]
    )
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in out.read_text(encoding="utf-8").splitlines() if line]
    assert rows, "stub emitted no hits for QRY_A"
    assert all(len(r) == len(fields.split(",")) for r in rows), "column count mismatch"


def test_the_foldseek_stub_routes_on_the_target_database(tmp_path):
    """Campaign-model hits and experimental-PDB hits are different tables. If
    the stub returned the same rows for both, the fixture would silently stop
    testing the campaign-overlay path."""
    fields = HITS["fields"]

    def targets_for(db_path: str) -> set[str]:
        out = tmp_path / f"{db_path.replace('/', '_')}.tsv"
        proc = run(
            [
                str(STUB_BIN / "foldseek"), "easy-search", "QRY_A.pdb", db_path,
                str(out), str(tmp_path / "tmp"), "--format-output", fields,
                "-e", "1e-3", "--max-seqs", "50",
            ]
        )
        assert proc.returncode == 0, proc.stderr
        return {line.split("\t")[1] for line in out.read_text().splitlines() if line}

    experimental = targets_for("databases/pdb_database")
    campaign = targets_for("work/campaign_structures")
    assert experimental and campaign
    assert experimental != campaign
    assert "QRY_A" in campaign, "self-match control missing from campaign_models"


def test_every_canned_hit_names_what_it_demonstrates():
    """An uncommented number in a fixture is a number nobody can safely change."""
    for hit in HITS["hits"]:
        assert hit.get("_comment"), f"hit {hit.get('target')} has no _comment"
    for hit in HITS["sequence_hits"]:
        assert hit.get("_comment"), f"sequence hit {hit.get('sseqid')} has no _comment"


def test_the_scenarios_partition_the_campaign_hits():
    """Both scenarios share one stub. A hit with no `_scenarios` key belongs to
    both; the trapped QRY_B hit must belong to `trap` alone, or the main
    scenario's audit would fail too and the fixture would prove nothing."""
    scoped = [h for h in HITS["hits"] if "_scenarios" in h]
    assert scoped, "no scenario-scoped hits: the two scenarios are identical"
    trap_only = [h for h in scoped if h["_scenarios"] == ["trap"]]
    assert trap_only, "no trap-only hit"
    assert all("toluene" in h["theader"].lower() for h in trap_only)


# --------------------------------------------------------------------------
# The input generator
# --------------------------------------------------------------------------

def test_the_generated_inputs_are_reproducible():
    """The query manifest carries SHA-256 checksums of the FASTA sequence and
    the PDB file, and the pipeline refuses to run on a mismatch. If the
    generator is not deterministic, the committed inputs and a fresh run
    disagree and the failure looks like a pipeline bug."""
    proc = run([sys.executable, str(HERE / "build_inputs.py"), "--check"])
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_both_scenarios_have_manifests():
    for name in (
        "query_manifest.json", "database_manifest.json",
        "query_manifest_trap.json", "database_manifest_trap.json",
    ):
        assert (INPUTS / name).exists(), f"{name} missing -- run build_inputs.py"


def test_no_manifest_declares_a_reserved_computed_field():
    """Inverted 2026-09-01, when the defect it recorded was fixed.

    This used to assert that the trap manifest declared `rbh: true` on QRY_B by
    hand, because nothing in the pipeline ever wrote the key and that was the
    only way to reach probable_same_function. Both halves changed: the RBH
    computation now runs before classification and reaches the label as a
    pairwise fact, and `reject_reserved_fields` refuses a curator-supplied `rbh`
    at manifest read time. A manifest declaring it no longer promotes a hit -- it
    fails the run. So the assertion is now the opposite one, and it guards the
    fix rather than recording the defect."""
    for name in ("query_manifest.json", "query_manifest_trap.json"):
        queries = json.loads((INPUTS / name).read_text(encoding="utf-8"))["queries"]
        flagged = [q["accession"] for q in queries if q.get("rbh")]
        assert not flagged, f"{name} declares rbh on {flagged}; the reader will reject the run"


def test_the_trap_manifest_declares_the_trap_it_tests():
    expectations = json.loads(
        (INPUTS / "database_manifest_trap.json").read_text(encoding="utf-8")
    )["release_expectations"]
    traps = expectations["title_traps"]
    assert traps, "trap scenario declares no title traps"
    assert any(t["substring"] == "toluene" for t in traps)


# --------------------------------------------------------------------------
# The two scenarios, end to end
# --------------------------------------------------------------------------

def pipeline(scenario: str, tmp_path: Path) -> tuple[Path, int]:
    """Run and verify one scenario; return the release path and the audit's exit code."""
    suffix = "_trap" if scenario == "trap" else ""
    out = tmp_path / scenario
    proc = run(
        [
            sys.executable, "-m", "sf_csa.cli", "run",
            "--queries", str(INPUTS / f"query_manifest{suffix}.json"),
            "--databases", str(INPUTS / f"database_manifest{suffix}.json"),
            "--output", str(out),
        ],
        scenario=scenario,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    audit = run(
        [
            sys.executable, "-m", "sf_csa.cli", "verify",
            "--output", str(out),
            "--databases", str(INPUTS / f"database_manifest{suffix}.json"),
        ],
        scenario=scenario,
    )
    return out, audit.returncode


def classifications(release: Path) -> dict[tuple[str, str], str]:
    import csv

    found = {}
    for path in sorted((release / "targets").glob("*/structure_hits.tsv")):
        for row in csv.DictReader(path.open(), delimiter="\t"):
            found[(row["query_id"], row["target_id"])] = row["function_classification"]
    return found


def test_the_main_scenario_audits_clean(tmp_path):
    release, rc = pipeline("main", tmp_path)
    assert rc == 0, "main scenario audit failed; the fixture's clean baseline is broken"
    assert (release / "SF_CSA_RELEASE_MANIFEST.json").exists()


def test_the_main_scenario_exercises_every_decision_path(tmp_path):
    """Five distinct outcomes on hand-placed numbers. If a threshold changes,
    one of these moves and the fixture says which."""
    release, _ = pipeline("main", tmp_path)
    labels = classifications(release)
    assert labels[("QRY_A", "QRY_A")] == "exact_function_supported"
    assert labels[("QRY_A", "SYN_WHOLE")] == "same_mechanism_class"
    assert labels[("QRY_A", "SYN_PART")] == "structural_analogy_only"
    assert labels[("QRY_B", "SYN_BELOW")] == "unresolved_or_conflicted"
    # Changed 2026-09-01 from same_mechanism_class. QRY_B now sits in proteome 1
    # as the top-scoring hit and points back at QRY_A in the reverse search, so
    # the pair is a computed reciprocal best hit and the structural leg is a
    # whole-architecture match in the same mechanism group. Both legs on one
    # target is what the label requires, and it is now reached by measurement
    # rather than by a manifest key.
    assert labels[("QRY_A", "QRY_B")] == "probable_same_function"


def test_the_main_scenario_reaches_probable_same_function_by_computation(tmp_path):
    """The signal this file said to watch for, arrived.

    Its predecessor asserted that QRY_A vs QRY_B stopped one rung short of
    probable_same_function because no code path set `rbh`, and said in as many
    words: "when that defect is fixed this test fails, which is the intended
    signal". The defect was fixed on 2026-09-01, so the assertion is inverted
    rather than deleted -- the fixture keeps its memory of what changed.

    What matters is *how* the label is reached. It must come from a reciprocal
    best hit the pipeline computed, not from a key a curator wrote, so the
    orthology status behind it is checked here too."""
    import csv as _csv

    release, _ = pipeline("main", tmp_path)
    assert classifications(release)[("QRY_A", "QRY_B")] == "probable_same_function"

    species = release / "targets" / "QRY_A" / "species_comparison.tsv"
    rows = list(_csv.DictReader(species.open(), delimiter="\t"))
    reciprocal = [r for r in rows
                  if r["target_accession"] == "QRY_B" and r["orthology_status"] == "reciprocal_best_hit"]
    assert reciprocal, "the label was reached without a computed reciprocal best hit behind it"


def test_the_trap_scenario_fails_its_audit(tmp_path):
    """The whole reason for a second scenario."""
    release, rc = pipeline("trap", tmp_path)
    assert rc != 0, "trap scenario audited clean: the title trap is not biting"
    assert classifications(release)[("QRY_A", "QRY_B")] == "probable_same_function"


def test_the_trap_scenario_names_the_substring_it_caught(tmp_path):
    _, rc = pipeline("trap", tmp_path)
    assert rc != 0
    audit = run(
        [
            sys.executable, "-m", "sf_csa.cli", "verify",
            "--output", str(tmp_path / "trap"),
            "--databases", str(INPUTS / "database_manifest_trap.json"),
        ],
        scenario="trap",
    )
    assert "toluene" in audit.stdout + audit.stderr


# --------------------------------------------------------------------------
# Golden comparison
# --------------------------------------------------------------------------

def canonical(release: Path, dest: Path) -> Path:
    proc = run(
        [
            sys.executable, str(HERE / "canonicalise.py"),
            "--release", str(release), "--dest", str(dest), "--fixture-root", str(HERE),
        ]
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return dest


def compare_trees(left: Path, right: Path) -> list[str]:
    """Recursive byte comparison; returns human-readable differences."""
    problems: list[str] = []
    left_files = {p.relative_to(left).as_posix() for p in left.rglob("*") if p.is_file()}
    right_files = {p.relative_to(right).as_posix() for p in right.rglob("*") if p.is_file()}
    for missing in sorted(left_files - right_files):
        problems.append(f"missing from run: {missing}")
    for extra in sorted(right_files - left_files):
        problems.append(f"not in golden: {extra}")
    for shared in sorted(left_files & right_files):
        if not filecmp.cmp(left / shared, right / shared, shallow=False):
            problems.append(f"content differs: {shared}")
    return problems


@pytest.mark.parametrize("scenario", ["main", "trap"])
def test_the_release_matches_its_golden_tree(scenario, tmp_path):
    """Compared on the canonical form: a raw release embeds its own output
    directory name in `release_id` and the absolute proteome paths in
    `proteome_denominator.json`, neither of which is a property of the code."""
    golden = GOLDEN / scenario
    if not golden.exists():
        pytest.skip(f"no golden tree for {scenario}; run run_fixture.sh --update-golden")
    release, _ = pipeline(scenario, tmp_path)
    produced = canonical(release, tmp_path / f"canon-{scenario}")
    problems = compare_trees(golden, produced)
    assert not problems, "golden mismatch:\n  " + "\n  ".join(problems)


def test_the_golden_trees_carry_live_checksums():
    """The canonicaliser recomputes CHECKSUMS.json over canonical bytes rather
    than blanking it, so a content change in a golden file shows up twice. If
    the digests were placeholders that second signal would be gone."""
    import hashlib

    for scenario in ("main", "trap"):
        root = GOLDEN / scenario
        if not root.exists():
            pytest.skip("no golden trees yet")
        digests = json.loads((root / "CHECKSUMS.json").read_text(encoding="utf-8"))
        assert digests, f"{scenario}: empty CHECKSUMS.json"
        for rel, digest in digests.items():
            target = root / rel
            assert target.exists(), f"{scenario}: {rel} checksummed but absent"
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            assert actual == digest, f"{scenario}: {rel} digest stale"


def test_the_two_golden_releases_actually_differ():
    """If main and trap produced identical releases, one of them would be
    decoration."""
    if not (GOLDEN / "main").exists() or not (GOLDEN / "trap").exists():
        pytest.skip("no golden trees yet")
    assert compare_trees(GOLDEN / "main", GOLDEN / "trap")


def test_the_canonicaliser_substitutes_the_longer_path_form_first():
    """A path substitution must not leave the prefix of a longer match behind.

    On macOS `/tmp` is a symlink to `/private/tmp`, so a fixture root at
    `/tmp/x` resolves to `/private/tmp/x` -- and the literal form is a *substring*
    of the resolved one. Replace the short form first and it rewrites the tail of
    the long one, leaving `/private<FIXTURE_ROOT>`: a machine-dependent string in
    a file whose whole purpose is to be machine-independent.

    This is a regression test for a real failure. The substitution iterated a
    `set`, so which form went first depended on hash order. The fixture passed in
    the directory it was authored in and failed on the first run from a symlinked
    path, which is the worst way for this to be discovered.

    The two forms are supplied through a stand-in rather than a real symlink,
    because reproducing a *prefix*-adding symlink (as opposed to a sibling one)
    requires writing at the filesystem root. The stand-in is what `_rewrite`
    consumes: something with `str()` and `.resolve()`.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_sfcsa_canonicalise_under_test", HERE / "canonicalise.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    class _RootWithSymlinkedPrefix:
        """Mimics /tmp/x, whose resolved form /private/tmp/x contains it."""

        def __init__(self, literal, resolved):
            self._literal = literal
            self._resolved = resolved

        def __str__(self):
            return self._literal

        def resolve(self):
            return self._resolved

    literal, resolved = "/tmp/fixture-root", "/private/tmp/fixture-root"
    assert literal in resolved, "test setup no longer reproduces the overlap"

    root = _RootWithSymlinkedPrefix(literal, resolved)

    for written_as in (literal, resolved):
        text = json.dumps({"path": f"{written_as}/inputs/proteomes/P1.faa"})
        out = module._rewrite(text, root, release_id="")
        assert "/private<FIXTURE_ROOT>" not in out, (
            f"written_as={written_as}: the short form was substituted first and "
            f"orphaned its prefix: {out}"
        )
        assert out == json.dumps(
            {"path": f"{module.FIXTURE_ROOT_TOKEN}/inputs/proteomes/P1.faa"}
        ), f"written_as={written_as}: {out}"
