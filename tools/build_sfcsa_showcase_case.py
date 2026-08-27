#!/usr/bin/env python3
"""Build one deterministic, public-safe SF-CSA pipeline showcase case.

The canonical SF-CSA package is executed through the reviewed offline fixture's
Foldseek and DIAMOND process-boundary test doubles. The stubs compute no
alignments. This case demonstrates orchestration, parsing, evidence separation,
classification boundaries, release verification, and deterministic artifacts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tools" / "fixtures" / "sfcsa"
DEFAULT_OUT = ROOT / "yauvi-structural-workbench" / "showcase" / "sfcsa-ceiling-case"
ARTIFACT = ROOT / "artifacts" / "protein-platform-modularization-reviewed.zip"


def canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sanitized(text: str, replacements: list[tuple[str, str]]) -> str:
    for source, target in sorted(replacements, key=lambda item: len(item[0]), reverse=True):
        text = text.replace(source, target)
    return text


def run(command: list[str], *, env: dict[str, str], replacements: list[tuple[str, str]]) -> tuple[int, str, str]:
    completed = subprocess.run(
        command, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False,
    )
    return (
        completed.returncode,
        sanitized(completed.stdout, replacements),
        sanitized(completed.stderr, replacements),
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build(output: Path, *, replace: bool = False) -> dict[str, Any]:
    output = output.resolve()
    if output.exists():
        if not replace:
            raise ValueError(f"output already exists: {output}; pass --replace")
        existing_case = output / "CASE.json"
        is_generated_case = False
        if existing_case.is_file():
            try:
                is_generated_case = json.loads(existing_case.read_text(encoding="utf-8")).get("case_id") == "HUC-06"
            except (OSError, json.JSONDecodeError):
                is_generated_case = False
        if output == ROOT or (ROOT not in output.parents and not is_generated_case):
            raise ValueError("refusing to replace a directory that is not a generated SF-CSA showcase")
        shutil.rmtree(output)
    output.mkdir(parents=True)

    inputs = output / "inputs"
    shutil.copytree(FIXTURE / "inputs", inputs)
    runtime = output / "runtime-fixture"
    runtime.mkdir()
    shutil.copyfile(FIXTURE / "stub_bin" / "hits.json", runtime / "hits.json")

    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "sf-csa" / "src") + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    env["PATH"] = str(FIXTURE / "stub_bin") + os.pathsep + env.get("PATH", "")
    env["SFCSA_FIXTURE_SCENARIO"] = "main"

    with tempfile.TemporaryDirectory(prefix="yauvi-sfcsa-showcase-") as temporary:
        scratch = Path(temporary)
        raw = scratch / "raw-release"
        canonical_release = scratch / "canonical-release"
        replacements = [(str(ROOT), "<repository>"), (str(output), "<showcase>"),
                        (str(scratch), "<scratch>")]
        run_code, run_stdout, run_stderr = run([
            sys.executable, "-m", "sf_csa.cli", "run",
            "--queries", str(inputs / "query_manifest.json"),
            "--databases", str(inputs / "database_manifest.json"),
            "--output", str(raw),
        ], env=env, replacements=replacements)
        verify_code, verify_stdout, verify_stderr = run([
            sys.executable, "-m", "sf_csa.cli", "verify",
            "--output", str(raw),
            "--databases", str(inputs / "database_manifest.json"),
        ], env=env, replacements=replacements)
        write_text(output / "SF_CSA_RUN_STDOUT.txt", run_stdout)
        write_text(output / "SF_CSA_RUN_STDERR.txt", run_stderr)
        write_text(output / "SF_CSA_VERIFY_STDOUT.txt", verify_stdout)
        write_text(output / "SF_CSA_VERIFY_STDERR.txt", verify_stderr)
        if run_code != 0 or verify_code != 0:
            raise RuntimeError(
                f"SF-CSA showcase fixture failed closed: run={run_code}, verify={verify_code}"
            )
        canonicalise = subprocess.run([
            sys.executable, str(FIXTURE / "canonicalise.py"),
            "--release", str(raw), "--dest", str(canonical_release),
            "--fixture-root", str(inputs),
        ], cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE,
           stderr=subprocess.PIPE, check=False)
        if canonicalise.returncode != 0:
            raise RuntimeError(f"SF-CSA canonicalisation failed: {canonicalise.stderr}")
        shutil.copytree(canonical_release, output / "release")

    release = output / "release"
    manifest = json.loads((release / "SF_CSA_RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    structural_rows: list[dict[str, str]] = []
    sequence_rows: list[dict[str, str]] = []
    for target in sorted((release / "targets").iterdir()):
        if target.is_dir():
            structural_rows.extend(read_tsv(target / "structure_hits.tsv"))
            sequence_rows.extend(read_tsv(target / "species_comparison.tsv"))
    labels: dict[str, int] = {}
    for row in structural_rows:
        label = row["function_classification"]
        labels[label] = labels.get(label, 0) + 1

    evidence_relatives = [
        "release/SF_CSA_RELEASE_MANIFEST.json",
        "release/RELEASE_COMPARISON_MATRIX.tsv",
        "release/targets/QRY_A/structure_hits.tsv",
        "release/targets/QRY_A/species_comparison.tsv",
        "release/targets/QRY_B/structure_hits.tsv",
        "release/targets/QRY_B/species_comparison.tsv",
        "release/CHECKSUMS.json",
    ]
    input_relatives = [
        "inputs/query_manifest.json", "inputs/database_manifest.json",
        "inputs/queries/QRY_A.faa", "inputs/queries/QRY_A.pdb",
        "inputs/queries/QRY_B.faa", "inputs/queries/QRY_B.pdb",
        "runtime-fixture/hits.json",
    ]
    case = {
        "schema_version": "1.0",
        "case_id": "HUC-06",
        "tool": "SF-CSA",
        "analysis_type": "sf_csa",
        "human_label": "How does structure-based evidence differ from sequence-based evidence?",
        "human_question": "Can the same invented proteins be compared through separate structural and sequence legs without turning similarity into exact functional transfer?",
        "test_state": "passed_stubbed_pipeline_case",
        "exit_codes": {"run": run_code, "verify": verify_code},
        "runtime_class": "deterministic_process_boundary_test_doubles",
        "observed_result": (
            f"{manifest['query_count']} invented queries completed a checksum-bound SF-CSA release; "
            "the release audit passed and structural classifications remained separate from sequence orthology candidates."
        ),
        "measurements": [
            {"label": "Release audit", "value": "passed", "help": "Canonical sf-csa verify returned exit code 0"},
            {"label": "Queries", "value": str(manifest["query_count"]), "help": "Invented checksum-bound query structures"},
            {"label": "Structural rows", "value": str(len(structural_rows)), "help": "Canned Foldseek-shaped rows interpreted by the real pipeline"},
            {"label": "Sequence rows", "value": str(len(sequence_rows)), "help": "Canned DIAMOND-shaped rows retained in a separate table"},
            {"label": "Evidence legs", "value": "2 separate", "help": "Structural similarity and sequence homology are not merged"},
            {"label": "Alignment engines", "value": "stubbed", "help": "No Foldseek or DIAMOND alignment was computed in this case"},
        ],
        "human_benefits": [
            "Shows where fold similarity and sequence homology agree or disagree.",
            "Prevents a strong structural match from silently becoming exact functional transfer.",
            "Preserves missing structures and unresolved relationships as visible evidence gaps.",
        ],
        "non_claim": (
            "Real alignment performance, biological function, orthology, substrate transfer, activity, "
            "pathogenic importance, or external scientific qualification."
        ),
        "runtime_disclosure": (
            "The canonical pipeline, subprocess construction, TSV parsing, classification, checksums, and "
            "release audit ran. Foldseek and DIAMOND were deterministic test doubles and computed no alignments."
        ),
        "known_findings": [
            "The title-trap protection currently acts during release verification, not direct classify_hit calls.",
            "The main pipeline does not currently feed computed sequence reciprocal-best-hit evidence back into structural classification.",
            "External mini-database runs with the installed Foldseek and DIAMOND binaries remain pending.",
        ],
        "classification_counts": dict(sorted(labels.items())),
        "input_sha256": {relative: sha256(output / relative) for relative in input_relatives},
        "evidence_files": evidence_relatives,
        "evidence_sha256": {relative: sha256(output / relative) for relative in evidence_relatives},
        "scientific_qualification_state": "external_benchmark_pending",
        "reviewed_artifact_sha256": sha256(ARTIFACT) if ARTIFACT.is_file() else "not_available",
    }
    write_text(output / "CASE.json", canonical(case))
    write_text(output / "README.md", """# SF-CSA public showcase case

This generated case runs the canonical SF-CSA pipeline through the reviewed
offline fixture's deterministic Foldseek and DIAMOND test doubles. The stubs
compute no alignments. The case demonstrates pipeline wiring, parsing, evidence
separation, classification, fail-closed release verification, and reproducible
artifacts—not biological performance or external scientific qualification.

Rebuild and verify from the repository root:

```bash
python tools/build_sfcsa_showcase_case.py --replace
python tools/verify_sfcsa_showcase_case.py
```
""")
    checksums = {
        path.relative_to(output).as_posix(): sha256(path)
        for path in sorted(output.rglob("*")) if path.is_file() and path.name != "CHECKSUMS.json"
    }
    write_text(output / "CHECKSUMS.json", canonical(checksums))
    return case


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    try:
        case = build(args.out, replace=args.replace)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"SF-CSA showcase build failed: {exc}", file=sys.stderr)
        return 2
    print(canonical({"showcase": str(args.out.resolve()), "case_id": case["case_id"],
                     "state": case["test_state"], "exit_codes": case["exit_codes"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
