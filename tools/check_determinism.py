#!/usr/bin/env python3
"""Run the shipped public case twice through the structural-workbench store.

The case is ``public-ca2-structure-coverage-v1``: chain A of RCSB PDB 2POW and
4RN4 against UniProt P00918. Inputs are the hash-locked files from
``examples/structural-portfolio/public_research_sources.json``. Each entry is
executed by ``StructuralAnalysisStore`` into its own workspace. Every emitted
file is byte-compared. Nothing is dropped for being a timestamp or a path.

A second mode, ``--compare``, reads one workspace from each CI runner and
reports cross-machine differences. A difference is accepted only when every
changed JSON value sits under a field listed in ``CROSS_ENVIRONMENT_FIELDS``.
Any other difference is a failure. Checksums and ``RAW_EVIDENCE.zip`` are
accepted only when the files they bind are themselves identical or explained
by that list.

Exit 0 when the two local runs match and, in ``--compare``, when every
cross-runner difference is explained. Exit 1 when bytes differ for a reason
that is not in the list. Exit 2 when the case cannot be run or the record is
unreadable.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "examples" / "structural-portfolio"
LOCK_PATH = CASE_DIR / "public_research_sources.json"
ACQUIRE = CASE_DIR / "acquire_public_inputs.py"
QUESTION = (
    "For chain A, what sequence identity and canonical-sequence coverage does "
    "StructQC report against UniProt P00918, and which canonical positions are "
    "unrepresented in the coordinate chain?"
)

# Within one machine the public case is byte-identical, including files that
# could have carried a clock or a path. This list is empty on purpose. Adding
# a name here drops that file from the byte compare, so the reason has to say
# why the variation is legitimate. An empty list means every emitted file is
# compared.
WITHIN_RUN_EXCLUSIONS: dict[str, str] = {}

# A cross-runner JSON difference is explained only when the path contains one
# of these markers. The compare prints the observed values; it does not delete
# them from the stored artifact.
CROSS_ENVIRONMENT_FIELDS: dict[str, str] = {
    "runtime_versions": (
        "Records the interpreter, operating system, CPU, and the numpy, scipy, "
        "gemmi, biopython, freesasa, foldseek, and diamond versions this process "
        "resolved. The run id is a hash of an identity that includes this object, "
        "so the id changes with it. These are not measurements of the structure."
    ),
    "run_id": (
        "First 16 hex digits of the SHA-256 of the case identity. The identity "
        "includes runtime_versions, so the directory name runs/run-<id> changes "
        "across the matrix. It is not a wall-clock timestamp."
    ),
    "gemmi_version": (
        "Parser build string stored beside the coordinate parse. Residue numbers "
        "in the evidence are rounded to 6 decimal places; this string is not one "
        "of those numbers. A difference in a coordinate or B-factor is not covered "
        "by this field and still fails the compare."
    ),
    "optional_runtimes.mkdssp": (
        "Whether the mkdssp executable is on PATH. This case does not invoke it. "
        "Runner images do not all ship it."
    ),
    "optional_runtimes.molprobity_or_phenix": (
        "Whether phenix.molprobity or clashscore is on PATH. This case does not "
        "invoke either."
    ),
}

RUN_ID_RE = re.compile(r"run-[0-9a-f]{16}")
GENERATOR_RE = re.compile(
    r"np\.random|numpy\.random|default_rng\(|os\.urandom\(|"
    r"secrets\.token_|secrets\.randbits\("
)

# Every generator hit under software/**/src and software/**/examples must match
# one of these substrings. A new hit fails the job instead of being folded into
# a seed field.
KNOWN_GENERATORS: dict[str, list[tuple[str, str]]] = {
    "software/membrane-orientor/memorient/src/memorient/orientor.py": [
        (
            "np.random.default_rng(0)",
            "MembraneOrient's rotation-invariance check draws rotations from "
            "numpy.random.default_rng(0). The 0 is a constant in that function, "
            "not an input. The public case is structure QC and does not call it. "
            "No seed field is added.",
        ),
    ],
    "software/membrane-orientor/memorient/examples/p4_benchmark.py": [
        (
            "np.random.default_rng(seed)",
            "Example benchmark script, not a store step and not the public case. "
            "Its seed argument is local to that script. No seed field is added "
            "to the case record.",
        ),
    ],
    "software/platform/src/yauvi_platform/structural_workbench/store.py": [
        (
            "os.urandom",
            "Ingest upload ids only. finalize_ingest deletes the ingest directory "
            "and does not write the id into the case. Repeated public-case runs "
            "are byte-identical, so this draw is not part of the emitted record.",
        ),
    ],
    "software/structural-workbench/src/yauvi_structural_workbench/server.py": [
        (
            "secrets.token_",
            "Loopback browser token and synthetic example ids. The public-case "
            "store run does not read them.",
        ),
    ],
    "software/structural-workbench/src/yauvi_structural_workbench/jobs.py": [
        (
            "secrets.token_",
            "Browser job ids. The public-case store run does not read them.",
        ),
    ],
    "software/platform/src/yauvi_platform/structural_workbench/protein_import.py": [
        (
            "secrets.token_",
            "Import record ids. Not part of the public-case store run.",
        ),
    ],
    "software/platform/src/yauvi_platform/structural_workbench/sources.py": [
        (
            "secrets.token_",
            "Source-acquisition ids. Not part of the public-case store run.",
        ),
    ],
}


@dataclass
class Sample:
    label: str
    files: dict[str, bytes]
    environment: dict[str, str] = field(default_factory=dict)


@dataclass
class Verdict:
    kind: str  # identical, explained, bug
    summary: str
    rows: list[str] = field(default_factory=list)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def logical_path(relative: str) -> str:
    """Align run directories whose names are environment-bound run ids."""
    return RUN_ID_RE.sub("run-ID", relative)


def snapshot(root: Path) -> dict[str, bytes]:
    files: dict[str, bytes] = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        files[path.relative_to(root).as_posix()] = path.read_bytes()
    return files


def logical_snapshot(files: dict[str, bytes]) -> dict[str, bytes]:
    out: dict[str, bytes] = {}
    for relative, data in files.items():
        key = logical_path(relative)
        if key in out:
            raise RuntimeError(f"two files collapsed to the same logical path: {key}")
        out[key] = data
    return out


def audit_generators() -> list[str]:
    """Search scientific source for generators. Unknown hits are failures."""
    unexpected: list[str] = []
    seen: set[str] = set()
    roots = [ROOT / "software"]
    for base in roots:
        for path in sorted(base.rglob("*.py")):
            relative = path.relative_to(ROOT).as_posix()
            if "/src/" not in f"/{relative}" and "/examples/" not in f"/{relative}":
                continue
            if "/tests/" in f"/{relative}":
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            for number, line in enumerate(text.splitlines(), 1):
                if GENERATOR_RE.search(line) is None:
                    continue
                seen.add(relative)
                allowed = KNOWN_GENERATORS.get(relative, [])
                if any(snippet in line for snippet, _reason in allowed):
                    continue
                unexpected.append(f"{relative}:{number}: {line.strip()}")
    missing = sorted(set(KNOWN_GENERATORS) - seen)
    for relative in missing:
        unexpected.append(f"{relative}: listed generator is no longer in the tree; update the audit")
    return unexpected


def print_generator_audit(unexpected: list[str]) -> None:
    print("Random-number audit (software/**/src and software/**/examples):")
    for relative, rows in KNOWN_GENERATORS.items():
        for _snippet, reason in rows:
            print(f"  {relative}")
            print(f"    {reason}")
    if unexpected:
        print("Unexpected generator usage:")
        for row in unexpected:
            print(f"  {row}")
    else:
        print("  No unlisted generator. The public case does not draw random numbers, so no seed is recorded.")


def _pointer(path: tuple[Any, ...]) -> str:
    text = ""
    for part in path:
        if isinstance(part, int):
            text += f"[{part}]"
        else:
            text += ("." if text else "") + str(part)
    return text or "$"


def _covered(path: tuple[Any, ...]) -> str | None:
    if "runtime_versions" in path:
        return "runtime_versions"
    if path and path[-1] == "run_id":
        return "run_id"
    if path and path[-1] == "gemmi_version":
        return "gemmi_version"
    if len(path) >= 2 and path[-2] == "optional_runtimes" and path[-1] in {"mkdssp", "molprobity_or_phenix"}:
        return f"optional_runtimes.{path[-1]}"
    return None


def _walk_diff(values: list[Any], path: tuple[Any, ...]) -> Iterable[tuple[tuple[Any, ...], list[Any]]]:
    if any(isinstance(value, dict) != isinstance(values[0], dict) for value in values) or any(
        isinstance(value, list) != isinstance(values[0], list) for value in values
    ):
        yield path, values
        return
    if isinstance(values[0], dict):
        keys: set[Any] = set()
        for value in values:
            keys.update(value)
        for key in sorted(keys, key=str):
            if any(key not in value for value in values):
                yield path + (key,), [value.get(key, "<missing>") for value in values]
                continue
            yield from _walk_diff([value[key] for value in values], path + (key,))
        return
    if isinstance(values[0], list):
        if any(len(value) != len(values[0]) for value in values):
            yield path, [f"<list length {len(value)}>" for value in values]
            return
        for index in range(len(values[0])):
            yield from _walk_diff([value[index] for value in values], path + (index,))
        return
    if any(value != values[0] for value in values):
        yield path, values


def _short(value: Any) -> str:
    text = repr(value)
    if len(text) > 180:
        return text[:177] + "..."
    return text


def _grouped(labels: list[str], values: list[Any]) -> str:
    groups: dict[str, list[str]] = {}
    for label, value in zip(labels, values):
        groups.setdefault(_short(value), []).append(label)
    return "; ".join(f"{text} [{', '.join(names)}]" for text, names in groups.items())


def classify_json(labels: list[str], blobs: list[bytes]) -> Verdict:
    if all(blob == blobs[0] for blob in blobs):
        return Verdict("identical", "byte-identical")
    try:
        values = [json.loads(blob) for blob in blobs]
    except json.JSONDecodeError as exc:
        return Verdict("bug", f"JSON does not parse: {exc}")
    diffs = list(_walk_diff(values, ()))
    if not diffs:
        return Verdict(
            "bug",
            "JSON bytes differ but the parsed values are equal. That is a formatting difference, not an environment field.",
        )
    rows: list[str] = []
    uncovered: list[str] = []
    for path, observed in diffs:
        marker = _covered(path)
        line = f"{_pointer(path)}: {_grouped(labels, observed)}"
        rows.append(line)
        if marker is None:
            uncovered.append(line)
    if uncovered:
        return Verdict("bug", "JSON values differ outside the stated environment fields.", uncovered)
    reasons = []
    used: set[str] = set()
    for path, _observed in diffs:
        marker = _covered(path)
        if marker and marker not in used:
            used.add(marker)
            reasons.append(f"{marker}: {CROSS_ENVIRONMENT_FIELDS[marker]}")
    return Verdict("explained", " ".join(reasons), rows)


def classify_html(labels: list[str], blobs: list[bytes]) -> Verdict:
    if all(blob == blobs[0] for blob in blobs):
        return Verdict("identical", "byte-identical")
    stripped = [RUN_ID_RE.sub("run-ID", blob.decode("utf-8", errors="replace")).encode("utf-8") for blob in blobs]
    if all(blob == stripped[0] for blob in stripped):
        return Verdict(
            "explained",
            "run_id: " + CROSS_ENVIRONMENT_FIELDS["run_id"],
            ["The only byte differences are run-[0-9a-f]{16} identifiers."],
        )
    return Verdict("bug", "HTML differs by more than the run id.", [_html_hint(blobs)])


def _html_hint(blobs: list[bytes]) -> str:
    texts = [blob.decode("utf-8", errors="replace") for blob in blobs]
    limit = min(len(text) for text in texts)
    for index in range(limit):
        chars = {text[index] for text in texts}
        if len(chars) > 1:
            start = max(0, index - 40)
            stop = min(limit, index + 40)
            return "first difference at " + str(index) + ": " + " || ".join(text[start:stop] for text in texts)
    return "files differ in length: " + ", ".join(str(len(blob)) for blob in blobs)


def classify_raw(labels: list[str], blobs: list[bytes]) -> Verdict:
    if all(blob == blobs[0] for blob in blobs):
        return Verdict("identical", "byte-identical")
    lengths = ", ".join(f"{label}={len(blob)}" for label, blob in zip(labels, blobs))
    return Verdict("bug", f"bytes differ ({lengths}). This file has no environment-field exception.", [_raw_hint(blobs)])


def _raw_hint(blobs: list[bytes]) -> str:
    limit = min(len(blob) for blob in blobs)
    for index in range(limit):
        if any(blob[index] != blobs[0][index] for blob in blobs):
            start = max(0, index - 24)
            stop = min(limit, index + 24)
            snippets = []
            for blob in blobs:
                snippets.append(blob[start:stop].decode("utf-8", errors="replace").replace("\n", "\\n"))
            return f"first difference at byte {index}: " + " || ".join(snippets)
    return "common prefix matches; lengths differ"


def classify_primary(labels: list[str], relative: str, blobs: list[bytes]) -> Verdict:
    name = Path(relative).name
    if name.endswith(".json"):
        return classify_json(labels, blobs)
    if name.endswith(".html"):
        return classify_html(labels, blobs)
    return classify_raw(labels, blobs)


def _zip_members(blob: bytes) -> dict[str, tuple[zipfile.ZipInfo, bytes]]:
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        return {info.filename: (info, archive.read(info)) for info in archive.infolist()}


def classify_checksum(labels: list[str], blobs: list[bytes], bound: dict[str, Verdict]) -> Verdict:
    """Accept a checksum change only when every changed entry names an explained file."""
    if all(blob == blobs[0] for blob in blobs):
        return Verdict("identical", "byte-identical")
    try:
        documents = [json.loads(blob) for blob in blobs]
    except json.JSONDecodeError as exc:
        return Verdict("bug", f"CHECKSUMS.json does not parse: {exc}")
    key_sets = [set(doc.get("files", {})) for doc in documents]
    if any(keys != key_sets[0] for keys in key_sets):
        return Verdict("bug", "CHECKSUMS.json does not name the same files on every runner.", [
            _grouped(labels, [sorted(keys) for keys in key_sets])
        ])
    without_files = []
    for doc in documents:
        copy = dict(doc)
        copy.pop("files", None)
        without_files.append(copy)
    meta = list(_walk_diff(without_files, ()))
    if meta:
        return Verdict("bug", "CHECKSUMS.json differs outside its file hash map.", [
            f"{_pointer(path)}: {_grouped(labels, observed)}" for path, observed in meta
        ])
    explained: list[str] = []
    problems: list[str] = []
    for name in sorted(key_sets[0]):
        hashes = [doc["files"][name] for doc in documents]
        target = bound.get(name)
        if len(set(hashes)) == 1:
            if target is not None and target.kind == "bug":
                problems.append(f"{name}: checksum is unchanged while the file differs unexpectedly")
            continue
        if target is None:
            problems.append(f"{name}: checksum changed and the compare has no classification for that file")
            continue
        if target.kind == "identical":
            problems.append(f"{name}: checksum changed even though the bound file is byte-identical")
            continue
        if target.kind == "bug":
            problems.append(f"{name}: checksum follows a file that differs outside the stated fields")
            continue
        explained.append(name)
    if problems:
        return Verdict("bug", "CHECKSUMS.json records a difference that is not explained by the bound files.", problems)
    if not explained:
        return Verdict("bug", "CHECKSUMS.json bytes differ without a hash change. That is formatting, not an environment field.")
    return Verdict(
        "explained",
        "Checksums bind the other artifacts. The entries that differ name files whose only differences are stated environment fields: "
        + ", ".join(explained),
        explained,
    )


def classify_zip(labels: list[str], blobs: list[bytes]) -> tuple[Verdict, dict[str, Verdict]]:
    if all(blob == blobs[0] for blob in blobs):
        members = _zip_members(blobs[0])
        return Verdict("identical", "byte-identical"), {name: Verdict("identical", "byte-identical") for name in members}
    try:
        opened = [_zip_members(blob) for blob in blobs]
    except zipfile.BadZipFile as exc:
        return Verdict("bug", f"RAW_EVIDENCE.zip is not a zip: {exc}"), {}
    names = [list(members) for members in opened]
    if any(name != names[0] for name in names):
        return Verdict("bug", "RAW_EVIDENCE.zip does not contain the same member names.", [
            _grouped(labels, names)
        ]), {}
    member_verdicts: dict[str, Verdict] = {}
    meta_problems: list[str] = []
    for name in names[0]:
        infos = [members[name][0] for members in opened]
        contents = [members[name][1] for members in opened]
        for attribute in ("date_time", "compress_type", "external_attr", "create_system"):
            observed = [getattr(info, attribute) for info in infos]
            if any(value != observed[0] for value in observed):
                meta_problems.append(f"{name} {attribute}: {_grouped(labels, observed)}")
        if name.endswith("CHECKSUMS.json"):
            continue
        member_verdicts[name] = classify_primary(labels, name, contents)
    checksum_names = [name for name in names[0] if name.endswith("CHECKSUMS.json")]
    for name in checksum_names:
        contents = [members[name][1] for members in opened]
        member_verdicts[name] = classify_checksum(labels, contents, member_verdicts)
    content_bugs = [f"{name}: {verdict.summary}" for name, verdict in member_verdicts.items() if verdict.kind == "bug"]
    if meta_problems or content_bugs:
        return Verdict(
            "bug",
            "RAW_EVIDENCE.zip differs in member metadata or in bytes that are not an environment field.",
            meta_problems + content_bugs,
        ), member_verdicts
    explained = [name for name, verdict in member_verdicts.items() if verdict.kind == "explained"]
    if not explained and any(blob != blobs[0] for blob in blobs):
        return Verdict(
            "bug",
            "RAW_EVIDENCE.zip bytes differ even though every member name, timestamp, attribute, and byte string matches. "
            "That is a container difference. It is reported rather than ignored; inspect the zip writers.",
            ["member metadata and member bytes match; container bytes do not"],
        ), member_verdicts
    detail = (
        "Archive members that differ do so only in stated environment fields: " + ", ".join(explained)
        if explained else "byte-identical"
    )
    return Verdict("explained" if explained else "identical", detail, explained), member_verdicts


def _run_prefix(files: dict[str, bytes], case_id: str) -> str | None:
    prefix = f"structural_analyses/cases/{case_id}/runs/run-ID/"
    if any(path.startswith(prefix) for path in files):
        return prefix
    return None


def compare_samples(samples: list[Sample]) -> tuple[list[tuple[str, Verdict]], int]:
    if len(samples) < 2:
        raise RuntimeError("at least two samples are required")
    labels = [sample.label for sample in samples]
    logical = [logical_snapshot(sample.files) for sample in samples]
    names = [set(files) for files in logical]
    findings: list[tuple[str, Verdict]] = []
    if any(name != names[0] for name in names):
        for sample, present in zip(samples, names):
            missing = sorted(names[0] - present) if sample is not samples[0] else []
            extra = sorted(present - names[0])
            if missing or extra:
                findings.append((
                    f"file set ({sample.label})",
                    Verdict("bug", "A runner emitted a different set of files.", [
                        f"missing: {missing}", f"extra: {extra}"
                    ]),
                ))
        common = sorted(set.intersection(*names)) if names else []
    else:
        common = sorted(names[0])
    primary = [
        path for path in common
        if not path.endswith("CHECKSUMS.json") and not path.endswith("RAW_EVIDENCE.zip")
    ]
    primary_verdicts: dict[str, Verdict] = {}
    for path in primary:
        verdict = classify_primary(labels, path, [files[path] for files in logical])
        primary_verdicts[path] = verdict
        findings.append((path, verdict))
    zip_paths = [path for path in common if path.endswith("RAW_EVIDENCE.zip")]
    zip_member_verdicts: dict[str, dict[str, Verdict]] = {}
    for path in zip_paths:
        verdict, members = classify_zip(labels, [files[path] for files in logical])
        zip_member_verdicts[path] = members
        findings.append((path, verdict))
    for path in [item for item in common if item.endswith("CHECKSUMS.json")]:
        case = _case_of(path)
        prefix = _run_prefix(logical[0], case) if case else None
        bound: dict[str, Verdict] = {}
        if prefix:
            for relative, verdict in primary_verdicts.items():
                if relative.startswith(prefix):
                    bound[relative[len(prefix):]] = verdict
            for zip_path, members in zip_member_verdicts.items():
                if prefix and zip_path.startswith(prefix):
                    bound["RAW_EVIDENCE.zip"] = next(item for item_path, item in findings if item_path == zip_path)
                    break
            export_zip = f"export/{case}/RAW_EVIDENCE.zip"
            if path.startswith("export/") and export_zip in zip_member_verdicts:
                bound["RAW_EVIDENCE.zip"] = next(item for item_path, item in findings if item_path == export_zip)
        verdict = classify_checksum(labels, [files[path] for files in logical], bound)
        findings.append((path, verdict))
    bugs = [item for item in findings if item[1].kind == "bug"]
    return findings, (1 if bugs else 0)


def _case_of(path: str) -> str | None:
    parts = path.split("/")
    if len(parts) >= 3 and parts[0] == "structural_analyses" and parts[1] == "cases":
        return parts[2]
    if len(parts) >= 2 and parts[0] == "export":
        return parts[1]
    return None


def print_findings(findings: list[tuple[str, Verdict]]) -> None:
    identical = [path for path, verdict in findings if verdict.kind == "identical"]
    explained = [(path, verdict) for path, verdict in findings if verdict.kind == "explained"]
    bugs = [(path, verdict) for path, verdict in findings if verdict.kind == "bug"]
    print(f"Identical files: {len(identical)}")
    for path in identical:
        print(f"  {path}")
    print(f"Files that differ for a stated reason: {len(explained)}")
    for path, verdict in explained:
        print(f"  {path}")
        print(f"    {verdict.summary}")
        for row in verdict.rows:
            print(f"    {row}")
    print(f"Unexpected differences: {len(bugs)}")
    for path, verdict in bugs:
        print(f"  {path}")
        print(f"    {verdict.summary}")
        for row in verdict.rows[:40]:
            print(f"    {row}")
        if len(verdict.rows) > 40:
            print(f"    ... {len(verdict.rows) - 40} more")


def within_run(first: dict[str, bytes], second: dict[str, bytes]) -> tuple[list[str], list[str]]:
    """Return (compared paths, unexpected differences). Exclusions are printed, not hidden."""
    print("Within-run exclusions:")
    if not WITHIN_RUN_EXCLUSIONS:
        print("  none. Every emitted file is byte-compared. No timestamp or absolute path was set aside.")
    else:
        for name, reason in sorted(WITHIN_RUN_EXCLUSIONS.items()):
            print(f"  {name}: {reason}")
    excluded = set(WITHIN_RUN_EXCLUSIONS)
    only_first = sorted(set(first) - set(second) - excluded)
    only_second = sorted(set(second) - set(first) - excluded)
    compared = sorted((set(first) & set(second)) - excluded)
    differed = [path for path in compared if first[path] != second[path]]
    problems = []
    if only_first:
        problems.append("only in the first run: " + ", ".join(only_first))
    if only_second:
        problems.append("only in the second run: " + ", ".join(only_second))
    problems.extend(differed)
    print(f"Files compared: {len(compared)}")
    for path in compared:
        print(f"  {path}")
    if problems:
        print("Within-run differences (these are bugs, not exclusions):")
        for item in problems:
            print(f"  {item}")
    else:
        print("Within-run result: byte-identical")
    return compared, problems


def _provenance(source: dict[str, Any]) -> dict[str, Any]:
    return {
        "class": "experimental",
        "method": source["method"],
        "resolution_angstrom": source["resolution_angstrom"],
        "source_id": f"RCSB PDB:{source['entry_id']}",
    }


def acquire_inputs(destination: Path) -> None:
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    completed = subprocess.run(
        [sys.executable, str(ACQUIRE), "--out", str(destination)],
        check=False,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"acquiring the locked public inputs failed with exit {completed.returncode}")


def run_workspace(workspace: Path, inputs: Path, lock: dict[str, Any], provenance_dir: Path) -> dict[str, Any]:
    from yauvi_platform.structural_workbench import StructuralAnalysisStore

    workspace.mkdir(parents=True)
    store = StructuralAnalysisStore(workspace)
    fasta = inputs / lock["reference"]["path"]
    cases = []
    for source in lock["structures"]:
        analysis_id = str(source["entry_id"]).lower()
        store.create(
            analysis_id,
            analysis_type="structure_qc",
            question=QUESTION,
            subject_id=str(source["uniprot_accession"]),
        )
        provenance = provenance_dir / f"{source['entry_id']}-provenance.json"
        provenance.write_text(
            json.dumps(_provenance(source), sort_keys=True, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
        store.add_file(analysis_id, role="structure", path=inputs / source["path"])
        store.add_file(analysis_id, role="reference_fasta", path=fasta)
        store.add_file(analysis_id, role="provenance", path=provenance)
        store.update_parameters(analysis_id, {"model": int(source["model"]), "chain": str(source["chain"])})
        preflight = store.preflight(analysis_id)
        if not preflight["valid"]:
            raise RuntimeError(f"preflight blocked {analysis_id}: {preflight}")
        record = store.run(analysis_id)
        if record["status"] not in {"completed", "scientifically_incomplete"}:
            raise RuntimeError(f"{analysis_id} did not produce a report bundle: {record['status']} {record.get('error', '')}")
        export_dir = workspace / "export" / analysis_id
        store.export(analysis_id, export_dir)
        cases.append({
            "analysis_id": analysis_id,
            "entry_id": source["entry_id"],
            "status": record["status"],
            "exit_code": record["exit_code"],
            "run_id": record["run_id"],
        })
        print(
            f"  {source['entry_id']}: status={record['status']} exit={record['exit_code']} "
            f"run_id={record['run_id']}"
        )
    return {"cases": cases}


def environment() -> dict[str, str]:
    return {
        "system": platform.system(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
    }


def write_record(out: Path, workspace: Path, compared: list[str], cases: list[dict[str, Any]]) -> None:
    if out.exists():
        shutil.rmtree(out)
    target = out / "workspace"
    shutil.copytree(workspace, target)
    manifest = {
        "schema_version": "1",
        "case_id": "public-ca2-structure-coverage-v1",
        "environment": environment(),
        "cases": cases,
        "files_compared": compared,
        "within_run_exclusions": WITHIN_RUN_EXCLUSIONS,
        "cross_environment_fields": CROSS_ENVIRONMENT_FIELDS,
        "file_sha256": {
            logical_path(path): sha256(data)
            for path, data in sorted(snapshot(workspace).items())
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def execute(out: Path, inputs: Path | None) -> int:
    unexpected = audit_generators()
    print_generator_audit(unexpected)
    if unexpected:
        return 1
    lock = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    if lock.get("case_id") != "public-ca2-structure-coverage-v1":
        raise RuntimeError(f"unexpected case lock: {lock.get('case_id')}")
    temporary = Path(tempfile.mkdtemp(prefix="yauvi-determinism-"))
    try:
        input_dir = inputs
        if input_dir is None:
            input_dir = temporary / "inputs"
            print(f"Acquiring locked public inputs into {input_dir}")
            acquire_inputs(input_dir)
        else:
            input_dir = input_dir.resolve()
        provenance_dir = temporary / "provenance"
        provenance_dir.mkdir()
        first_root = temporary / "run-1"
        second_root = temporary / "run-2"
        print(f"Run 1: {first_root}")
        first = run_workspace(first_root, input_dir, lock, provenance_dir)
        print(f"Run 2: {second_root}")
        second = run_workspace(second_root, input_dir, lock, provenance_dir)
        if first["cases"] != second["cases"]:
            print("The two runs did not record the same status, exit code, or run id.")
            print(json.dumps(first, indent=2, sort_keys=True))
            print(json.dumps(second, indent=2, sort_keys=True))
            return 1
        print(
            "Store status is scientifically_incomplete when no external validation "
            "report is attached: the registered structure_qc adapter passes "
            "--require-external-validation, and this public case does not supply one. "
            "Both runs must agree. Incomplete is a recorded result, not a hidden retry."
        )
        compared, problems = within_run(snapshot(first_root), snapshot(second_root))
        if problems:
            return 1
        write_record(out, first_root, compared, first["cases"])
        print(f"Wrote {out}")
        return 0
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def load_samples(directory: Path) -> list[Sample]:
    samples = []
    for child in sorted(path for path in directory.iterdir() if path.is_dir()):
        manifest_path = child / "manifest.json"
        workspace = child / "workspace"
        if not manifest_path.is_file() or not workspace.is_dir():
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        env = manifest.get("environment", {})
        label = f"{env.get('system', '?')}/{env.get('machine', '?')} py{env.get('python', '?')}"
        files = snapshot(workspace)
        recorded = manifest.get("file_sha256", {})
        logical = logical_snapshot(files)
        if recorded and recorded != {path: sha256(data) for path, data in logical.items()}:
            raise RuntimeError(f"{child.name}: manifest hashes do not match the uploaded workspace")
        samples.append(Sample(label=label, files=files, environment={str(k): str(v) for k, v in env.items()}))
    if len(samples) < 2:
        raise RuntimeError(f"need at least two determinism records in {directory}")
    return samples


def cross_compare(directory: Path) -> int:
    samples = load_samples(directory)
    print("Cross-environment samples:")
    for sample in samples:
        print(f"  {sample.label}")
    print("Fields that may differ, and why:")
    for name, reason in CROSS_ENVIRONMENT_FIELDS.items():
        print(f"  {name}: {reason}")
    findings, code = compare_samples(samples)
    print_findings(findings)
    return code


def _canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def _fixture_zip(members: list[tuple[str, bytes]]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, content in sorted(members):
            info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, content)
    return buffer.getvalue()


def _self_test_tree(run_id: str, python: str, mkdssp: str, residue: str) -> dict[str, bytes]:
    prefix = f"structural_analyses/cases/2pow/runs/{run_id}/"
    analysis = {
        "run_id": run_id,
        "runtime_versions": {"machine": "x86_64", "platform": "Linux", "python": python},
        "status": "scientifically_incomplete",
    }
    evidence = {"coordinate": {"parser": {"gemmi_version": "0.7.5"}}, "residues": [{"mean_b_factor": 1.0}]}
    module_manifest = {
        "optional_runtimes": {"gemmi": "available_invoked", "mkdssp": mkdssp},
        "runtime_versions": {"numpy": "2.5.4", "python": python},
    }
    report = {"documents": [{"document": evidence}, {"document": module_manifest}], "run_id": run_id}
    html = f"<p>Run: <code>{run_id}</code></p>\n".encode("utf-8")
    table = f"mean_b_factor\n{residue}\n".encode("utf-8")
    files = {
        prefix + "ANALYSIS_RUN.json": _canonical(analysis),
        prefix + "REPORT_DATA.json": _canonical(report),
        prefix + "REPORT.html": html,
        prefix + "outputs/structqc/STRUCTURE_EVIDENCE.json": _canonical(evidence),
        prefix + "outputs/structqc/RUN_MANIFEST.json": _canonical(module_manifest),
        prefix + "outputs/structqc/RESIDUE_QUALITY.tsv": table,
        "structural_analyses/cases/2pow/PREFLIGHT.json": _canonical({
            "readiness": {"optional_runtimes": {"gemmi": "available", "mkdssp": "missing", "molprobity_or_phenix": "missing"}},
            "valid": True,
        }),
    }
    members = [(path[len(prefix):], data) for path, data in files.items() if path.startswith(prefix)]
    digest_map = {name: sha256(data) for name, data in members}
    checksum = _canonical({"files": digest_map, "schema_version": "1.0"})
    files[prefix + "CHECKSUMS.json"] = checksum
    files[prefix + "RAW_EVIDENCE.zip"] = _fixture_zip(members + [("CHECKSUMS.json", checksum)])
    # Outside checksum also names the zip, matching the store.
    # The store records the zip in the checksum written after the archive.
    # The copy inside the archive does not name the archive or itself.
    outside = dict(digest_map)
    outside["RAW_EVIDENCE.zip"] = sha256(files[prefix + "RAW_EVIDENCE.zip"])
    files[prefix + "CHECKSUMS.json"] = _canonical({"files": outside, "schema_version": "1.0"})
    return files


def self_test() -> int:
    base = _self_test_tree("run-aaaaaaaaaaaaaaaa", "3.12.3", "not_available", "1.25")
    same = _self_test_tree("run-aaaaaaaaaaaaaaaa", "3.12.3", "not_available", "1.25")
    _findings, code = compare_samples([
        Sample("a", base), Sample("b", same),
    ])
    if code != 0:
        print("self-test failed: identical fixtures were reported as different")
        return 1
    shifted = _self_test_tree("run-bbbbbbbbbbbbbbbb", "3.10.16", "available_not_invoked", "1.25")
    # Rebuild is independent; mkdssp and python and run id change. Residue does not.
    # PREFLIGHT mkdssp stays "missing" in the fixture, so patch it.
    shifted["structural_analyses/cases/2pow/PREFLIGHT.json"] = _canonical({
        "readiness": {"optional_runtimes": {"gemmi": "available", "mkdssp": "available", "molprobity_or_phenix": "missing"}},
        "valid": True,
    })
    findings, code = compare_samples([Sample("linux", base), Sample("macos", shifted)])
    if code != 0:
        print("self-test failed: an environment-only change was reported as a bug")
        print_findings(findings)
        return 1
    if not any(verdict.kind == "explained" for _path, verdict in findings):
        print("self-test failed: environment change produced no explained difference")
        return 1
    mutated = _self_test_tree("run-aaaaaaaaaaaaaaaa", "3.12.3", "not_available", "9.99")
    findings, code = compare_samples([Sample("a", base), Sample("b", mutated)])
    if code == 0:
        print("self-test failed: a residue-table change was accepted")
        print_findings(findings)
        return 1
    whitespace = dict(base)
    evidence_path = "structural_analyses/cases/2pow/runs/run-aaaaaaaaaaaaaaaa/outputs/structqc/STRUCTURE_EVIDENCE.json"
    whitespace[evidence_path] = json.dumps(json.loads(base[evidence_path]), indent=2).encode("utf-8")
    _findings, code = compare_samples([Sample("a", base), Sample("b", whitespace)])
    if code == 0:
        print("self-test failed: JSON formatting difference was accepted")
        return 1
    gemmi_state = json.loads(base["structural_analyses/cases/2pow/PREFLIGHT.json"])
    gemmi_state["readiness"]["optional_runtimes"]["gemmi"] = "missing"
    changed = dict(base)
    changed["structural_analyses/cases/2pow/PREFLIGHT.json"] = _canonical(gemmi_state)
    _findings, code = compare_samples([Sample("a", base), Sample("b", changed)])
    if code == 0:
        print("self-test failed: a gemmi-availability change was treated as ignorable")
        return 1
    unexpected = audit_generators()
    if unexpected:
        print("self-test failed: generator audit is not clean against this tree")
        for row in unexpected:
            print(row)
        return 1
    print("self-test passed")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="directory for the hash record and the first-run workspace")
    parser.add_argument("--inputs", type=Path, help="directory that already holds the locked public files")
    parser.add_argument("--compare", type=Path, help="directory of downloaded determinism records")
    parser.add_argument("--self-test", action="store_true", help="check the comparer and the generator audit, then exit")
    args = parser.parse_args(argv)
    try:
        if args.self_test:
            return self_test()
        if self_test() != 0:
            return 1
        if args.compare:
            return cross_compare(args.compare)
        if args.out is None:
            parser.error("--out is required unless --compare or --self-test is set")
        return execute(args.out, args.inputs)
    except (OSError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"check_determinism: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
