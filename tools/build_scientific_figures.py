#!/usr/bin/env python3
"""Rebuild the JOSS paper figures from regenerated store runs.

``--check`` first reruns the shipped synthetic StructQC example and the public
CA2 case (PDB 2POW and 4RN4) through StructuralAnalysisStore. It then requires
each run id and artifact SHA-256 to match the pins in FIGURE_INPUTS.json.
The run directories are not committed. Schematic figures and the qualification
history receipt have no store run id.

    python tools/build_scientific_figures.py
    python tools/build_scientific_figures.py --check

The pins were recorded on Linux x86_64, CPython 3.12.3, with
``yauvi-structural-workbench/docs/requirements-linux-x86_64-py312.txt``.
A run id hashes that runtime, so another interpreter or library set fails
the pin check. CA2 acquisition needs network access to RCSB and UniProt.
The downloaded coordinates are not written into the repository.

``--check`` then rebuilds the figures and byte-compares them with the
committed copies, the same way ``tools/build_cli_reference.py --check``
compares the CLI reference. PNG bytes can still differ across FreeType or
libpng builds. When they do, the check compares the SVG and the extracted
plotted data and says why. A change in the plotted data fails.

Matplotlib is not part of the scientific install. The figure command needs
``matplotlib==3.11.2``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import statistics
import sys
import tempfile
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO / "yauvi-structural-workbench" / "paper"
FIGURES_DIR = PAPER_DIR / "figures"
DEFAULT_INPUTS = FIGURES_DIR / "FIGURE_INPUTS.json"
PROVENANCE_NAME = "FIGURE_PROVENANCE.json"
PAPER_MD = PAPER_DIR / "paper.md"
GENERATOR = Path(__file__)

SOURCE_DATE_EPOCH = "1700000000"
SVG_HASHSALT = "yauvi-paper-figures"
MATPLOTLIB_SPEC = "matplotlib==3.11.2"
RUN_ID_RE = re.compile(r"^run-[0-9a-f]{16}$")
FIGURE_ORDER = (
    "architecture",
    "worked-example",
    "qualification-history",
    "six-evidence-questions",
    "ca2-coverage",
)
NAVY = "#193641"
TEAL = "#26786b"
GOLD = "#a16929"

PNG_DIFFERENCE_REASON = (
    "PNG bytes rasterize text through FreeType and compress with libpng. "
    "Those libraries can emit different bytes on another platform even when "
    "the drawing commands are the same. Drawing uses the Agg backend. SVG "
    "text stays text (svg.fonttype none), clip-path ids use a fixed "
    "svg.hashsalt, and SOURCE_DATE_EPOCH fixes the embedded date. When PNG "
    "bytes differ, this check compares the SVG and the extracted plotted "
    "data. It fails when that plotted data changes."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def dump_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


_MATPLOTLIB = None


def configure_matplotlib():
    global _MATPLOTLIB
    if _MATPLOTLIB is not None:
        return _MATPLOTLIB
    os.environ["SOURCE_DATE_EPOCH"] = SOURCE_DATE_EPOCH
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "yauvi-figure-cache"))
    try:
        import matplotlib
    except ImportError:
        print(
            f"matplotlib is required to build the paper figures. Install {MATPLOTLIB_SPEC}.",
            file=sys.stderr,
        )
        raise SystemExit(2)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "svg.fonttype": "none",
        "svg.hashsalt": SVG_HASHSALT,
        "path.simplify": False,
        "agg.path.chunksize": 0,
    })
    _MATPLOTLIB = (plt, FancyBboxPatch)
    return _MATPLOTLIB


def save_figure(plt, fig, dest: Path, name: str) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest / f"{name}.png", dpi=200, bbox_inches="tight", facecolor="white")
    svg_path = dest / f"{name}.svg"
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")


def cited_in_paper(name: str) -> bool:
    return f"figures/{name}.png" in PAPER_MD.read_text(encoding="utf-8")


SYNTHETIC_ANALYSIS = "synthetic-structqc"
SYNTHETIC_QUESTION = "What does this two-residue fixture measure?"
SYNTHETIC_SUBJECT = "SYNTHETIC"
SYNTHETIC_EXAMPLES = REPO / "software" / "structqc" / "examples"
PINNED_LIBRARIES = ("numpy", "scipy", "gemmi", "biopython")
ABSENT_FROM_PIN = ("freesasa", "foldseek", "diamond")


def constraint_versions(path: Path) -> dict[str, str]:
    versions: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if "==" not in line:
            continue
        name, version = line.split("==", 1)
        versions[name.strip().lower()] = version.strip()
    return versions


def assert_pin_environment(manifest: dict[str, Any]) -> None:
    pin = manifest.get("pin_environment")
    if not isinstance(pin, dict):
        raise SystemExit("FIGURE_INPUTS.json has no pin_environment")
    problems: list[str] = []
    for key, actual in (
        ("system", platform.system()),
        ("machine", platform.machine()),
        ("python", platform.python_version()),
    ):
        if actual != pin.get(key):
            problems.append(f"{key}: pin {pin.get(key)!r}, this process {actual!r}")
    constraint_rel = pin.get("constraint")
    constraint = REPO / str(constraint_rel) if constraint_rel else None
    if constraint is None or not constraint.is_file():
        problems.append(f"missing constraint file {constraint_rel}")
    else:
        pinned = constraint_versions(constraint)
        import importlib.metadata as metadata

        for name in PINNED_LIBRARIES:
            expected = pinned.get(name)
            try:
                actual = metadata.version(name)
            except metadata.PackageNotFoundError:
                actual = "unavailable"
            if expected is None:
                problems.append(f"{name} is not pinned in {constraint_rel}")
            elif actual != expected:
                problems.append(f"{name}: constraint {expected}, installed {actual}")
    for binary in ABSENT_FROM_PIN:
        found = shutil.which(binary)
        if found:
            problems.append(f"{binary} is on PATH ({found}); the pins recorded it as unavailable")
    if problems:
        raise SystemExit(
            "Figure run-id pins were recorded on "
            f"{pin.get('system')} {pin.get('machine')}, CPython {pin.get('python')}, "
            f"with {pin.get('constraint')}.\n" + "\n".join(problems)
        )


def _run_dir(workspace: Path, analysis_id: str, run_id: str) -> Path:
    return workspace / "structural_analyses" / "cases" / analysis_id / "runs" / run_id


def _load_determinism():
    tools = str(REPO / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import check_determinism

    return check_determinism


def execute_public_runs(scratch: Path) -> dict[str, Path]:
    """Rerun the synthetic fixture and the public CA2 case into one temporary store."""
    from yauvi_platform.structural_workbench import StructuralAnalysisStore

    determinism = _load_determinism()
    workspace = scratch / "workspace"
    workspace.mkdir()
    store = StructuralAnalysisStore(workspace)
    print("Regenerating the synthetic StructQC run from software/structqc/examples", flush=True)
    store.create(
        SYNTHETIC_ANALYSIS,
        analysis_type="structure_qc",
        question=SYNTHETIC_QUESTION,
        subject_id=SYNTHETIC_SUBJECT,
    )
    for role, name in (
        ("structure", "model.pdb"),
        ("reference_fasta", "reference.fasta"),
        ("provenance", "provenance.json"),
        ("validation_report", "validation.json"),
    ):
        store.add_file(SYNTHETIC_ANALYSIS, role=role, path=SYNTHETIC_EXAMPLES / name)
    store.update_parameters(SYNTHETIC_ANALYSIS, {"model": 0, "chain": "A"})
    record = store.run(SYNTHETIC_ANALYSIS)
    print(f"  {SYNTHETIC_ANALYSIS}: status={record['status']} run_id={record['run_id']}", flush=True)
    if record["status"] != "completed":
        raise SystemExit(f"{SYNTHETIC_ANALYSIS} status {record['status']}: {record.get('error', '')}")
    runs = {SYNTHETIC_ANALYSIS: _run_dir(workspace, SYNTHETIC_ANALYSIS, record["run_id"])}

    lock = load_json(determinism.LOCK_PATH)
    if lock.get("case_id") != "public-ca2-structure-coverage-v1":
        raise SystemExit(f"unexpected public case lock: {lock.get('case_id')}")
    inputs = scratch / "inputs"
    provenance_dir = scratch / "provenance"
    provenance_dir.mkdir()
    print(f"Acquiring locked public CA2 inputs into {inputs}", flush=True)
    try:
        determinism.acquire_inputs(inputs)
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    try:
        fasta = inputs / lock["reference"]["path"]
        for source in lock["structures"]:
            analysis_id = str(source["entry_id"]).lower()
            store.create(
                analysis_id,
                analysis_type="structure_qc",
                question=determinism.QUESTION,
                subject_id=str(source["uniprot_accession"]),
            )
            provenance = provenance_dir / f"{source['entry_id']}-provenance.json"
            provenance.write_text(
                json.dumps(determinism._provenance(source), sort_keys=True, separators=(",", ":")) + "\n",
                encoding="utf-8",
            )
            store.add_file(analysis_id, role="structure", path=inputs / source["path"])
            store.add_file(analysis_id, role="reference_fasta", path=fasta)
            store.add_file(analysis_id, role="provenance", path=provenance)
            store.update_parameters(analysis_id, {"model": int(source["model"]), "chain": str(source["chain"])})
            preflight = store.preflight(analysis_id)
            if not preflight["valid"]:
                raise SystemExit(f"preflight blocked {analysis_id}: {preflight}")
            record = store.run(analysis_id)
            print(
                f"  {source['entry_id']}: status={record['status']} "
                f"exit={record['exit_code']} run_id={record['run_id']}",
                flush=True,
            )
            if record["status"] not in {"completed", "scientifically_incomplete"}:
                raise SystemExit(
                    f"{analysis_id} did not produce a report bundle: "
                    f"{record['status']} {record.get('error', '')}"
                )
            runs[analysis_id] = _run_dir(workspace, analysis_id, record["run_id"])
    finally:
        shutil.rmtree(inputs, ignore_errors=True)
        shutil.rmtree(provenance_dir, ignore_errors=True)
    return runs


def verify_pins(manifest: dict[str, Any], runs: dict[str, Path]) -> None:
    problems: list[str] = []
    seen: set[str] = set()
    for name in FIGURE_ORDER:
        spec = manifest["figures"][name]
        if spec["input_class"] != "store_run":
            continue
        for run_spec in spec["runs"]:
            analysis_id = run_spec["analysis_id"]
            seen.add(analysis_id)
            pinned_id = run_spec["run_id"]
            if not RUN_ID_RE.fullmatch(pinned_id):
                problems.append(f"{analysis_id} pin is not a run id: {pinned_id}")
                continue
            root = runs.get(analysis_id)
            if root is None:
                problems.append(f"no regenerated run for {analysis_id}")
                continue
            record_path = root / "ANALYSIS_RUN.json"
            if not record_path.is_file():
                problems.append(f"no ANALYSIS_RUN.json for {analysis_id} at {root}")
                continue
            record = load_json(record_path)
            actual_id = record.get("run_id")
            if actual_id != pinned_id or root.name != pinned_id:
                problems.append(
                    f"{analysis_id} regenerated {actual_id} does not match pin {pinned_id}; "
                    f"runtime_versions={json.dumps(record.get('runtime_versions'), sort_keys=True)}"
                )
            for artifact in run_spec["artifacts"]:
                path = root / artifact["path"]
                if not path.is_file():
                    problems.append(f"missing {analysis_id} {artifact['role']}: {artifact['path']}")
                    continue
                digest = sha256(path)
                if digest != artifact.get("sha256"):
                    problems.append(
                        f"{analysis_id} {artifact['role']} sha256 {digest} does not match pin {artifact.get('sha256')}"
                    )
    extra = set(runs) - seen
    if extra:
        problems.append("regenerated analyses are not pinned: " + ", ".join(sorted(extra)))
    if problems:
        raise SystemExit("Regenerated runs do not match the pins in FIGURE_INPUTS.json:\n" + "\n".join(problems))
    print("Regenerated run ids and artifact SHA-256s match FIGURE_INPUTS.json.", flush=True)


def prepare_runs(manifest: dict[str, Any]) -> tuple[Path, dict[str, Path]]:
    assert_pin_environment(manifest)
    scratch = Path(tempfile.mkdtemp(prefix="yauvi-figure-runs-"))
    try:
        runs = execute_public_runs(scratch)
        verify_pins(manifest, runs)
    except BaseException:
        shutil.rmtree(scratch, ignore_errors=True)
        raise
    return scratch, runs


def require_hash(path: Path, expected: str | None, role: str) -> str:
    if not path.is_file():
        raise SystemExit(f"missing {role}: {path}")
    actual = sha256(path)
    if expected and actual != expected:
        raise SystemExit(f"{role} sha256 mismatch for {path}: manifest {expected}, file {actual}")
    return actual


def resolve_repo_artifact(spec: dict[str, str]) -> dict[str, str]:
    path = REPO / spec["path"]
    digest = require_hash(path, spec.get("sha256"), spec["role"])
    return {"role": spec["role"], "path": spec["path"], "sha256": digest}


def resolve_run(spec: dict[str, Any], runs: dict[str, Path]) -> dict[str, Any]:
    analysis_id = spec["analysis_id"]
    root = runs[analysis_id]
    record = load_json(root / "ANALYSIS_RUN.json")
    run_id = record["run_id"]
    if run_id != spec["run_id"]:
        raise SystemExit(f"{analysis_id} run id {run_id} does not match pin {spec['run_id']}")
    artifacts = []
    for artifact in spec["artifacts"]:
        path = root / artifact["path"]
        digest = require_hash(path, artifact.get("sha256"), f"{analysis_id} {artifact['role']}")
        artifacts.append({
            "role": artifact["role"],
            "run_id": run_id,
            "analysis_id": analysis_id,
            "path": artifact["path"],
            "sha256": digest,
            "resolved_path": path,
        })
    return {
        "analysis_id": analysis_id,
        "entry_id": spec.get("entry_id", analysis_id),
        "run_id": run_id,
        "artifacts": artifacts,
        "run_dir": root,
    }


def artifact_path(run: dict[str, Any], role: str) -> Path:
    for artifact in run["artifacts"]:
        if artifact["role"] == role:
            return artifact["resolved_path"]
    raise SystemExit(f"{run['analysis_id']} has no {role} artifact")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def position_ranges(values: list[int]) -> str:
    if not values:
        return "none"
    ordered = sorted(values)
    chunks: list[str] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value == previous + 1:
            previous = value
            continue
        chunks.append(str(start) if start == previous else f"{start}–{previous}")
        start = previous = value
    chunks.append(str(start) if start == previous else f"{start}–{previous}")
    return ", ".join(chunks)


def draw_architecture(plt, FancyBboxPatch, dest: Path) -> dict[str, Any]:
    fig, ax = plt.subplots(num=1, figsize=(11, 3.6))
    ax.axis("off")
    ax.set_xlim(-0.16, 11)
    ax.set_ylim(-0.16, 3.6)
    ax.text(0, 3.3, "A shared path from input to interpretation", fontsize=18, fontweight="bold", color=NAVY)
    labels = [
        ("Input", "Structure + source\ncontent hashes"),
        ("Mapping", "Model, assembly,\nchain, residue"),
        ("Analysis", "Shared engine\nCLI or browser"),
        ("Evidence", "Measurements, gaps,\nimmutable run"),
        ("Interpretation", "Bounded claims\nportable report"),
    ]
    for index, (title, body) in enumerate(labels):
        x = index * 2.2
        ax.add_patch(FancyBboxPatch((x, 0.9), 1.9, 1.6, boxstyle="round,pad=.08", facecolor="#edf4f2", edgecolor=TEAL))
        ax.text(x + 0.95, 2.05, title, ha="center", fontweight="bold", color=NAVY)
        ax.text(x + 0.95, 1.45, body, ha="center", va="center", fontsize=10)
        if index < 4:
            ax.annotate("", xy=(x + 2.08, 1.7), xytext=(x + 1.95, 1.7), arrowprops={"arrowstyle": "->", "color": NAVY})
    status = "Execution outcome • Scientific qualification • Independent reproduction • Publication approval"
    limit = "Four separate states; success in one does not establish the others."
    ax.text(0, 0.3, status, fontsize=11, color=NAVY)
    ax.text(0, -0.1, limit, fontsize=10, color=GOLD)
    save_figure(plt, fig, dest, "architecture")
    return {"stages": [{"title": title, "body": body} for title, body in labels], "status_line": status, "limit": limit}


def draw_worked_example(plt, dest: Path, run: dict[str, Any]) -> dict[str, Any]:
    evidence = load_json(artifact_path(run, "structure_evidence"))
    rows = read_tsv(artifact_path(run, "residue_quality"))
    completeness = evidence["completeness"]
    clash = evidence["external_validation"]["metrics"]["clashscore"]
    if evidence["external_validation"]["state"] != "imported":
        raise SystemExit("worked-example expects an imported validation report on the synthetic run")
    if evidence["pae"]["state"] != "unevaluated":
        raise SystemExit("worked-example expects PAE to be unevaluated")
    points = []
    for row in rows:
        x_value, y_value, z_value = (float(part) for part in row["ca_xyz"].split(";"))
        points.append({
            "chain_id": row["chain_id"],
            "auth_seq_id": row["auth_seq_id"],
            "chemical_component_id": row["chemical_component_id"],
            "ca_xyz": row["ca_xyz"],
            "x": x_value,
            "y": y_value,
            "z": z_value,
        })
    lines = [
        f"Mapped residues: {completeness['mapped_residues']} of {completeness['reference_length']}",
        f"Sequence identity: {completeness['identity_fraction']:.0%} within this toy mapping",
        f"Imported clashscore: {clash} (synthetic fixture)",
        "PAE: not supplied; not inferred",
        "Input and validation records identified by SHA-256",
        "Limit: these two residues demonstrate software behavior.",
        "They do not establish structure validity or function.",
    ]
    fig, (ax, side) = plt.subplots(1, 2, num=2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [1, 1.5]})
    ax.plot([point["x"] for point in points], [point["y"] for point in points], color=TEAL, lw=2)
    for row, point in zip(rows, points):
        ax.scatter(point["x"], point["y"], s=120, color=TEAL)
        ax.annotate(
            f"{row['chain_id']}:{row['auth_seq_id']} {row['chemical_component_id']}",
            (point["x"], point["y"]),
            xytext=(8, -18),
            textcoords="offset points",
        )
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_title("Cα coordinates · synthetic example")
    ax.set_xlim(-0.6, 4)
    ax.set_ylim(-0.8, 2.4)
    ax.spines[["top", "right"]].set_visible(False)
    side.axis("off")
    side.text(0, 0.94, "What was actually measured?", fontsize=16, fontweight="bold", color=NAVY)
    for index, text in enumerate(lines):
        side.text(0, 0.78 - index * 0.105, text, color=NAVY if index < 5 else GOLD, fontsize=10.5)
    save_figure(plt, fig, dest, "worked-example")
    return {
        "run_id": run["run_id"],
        "analysis_id": run["analysis_id"],
        "points": points,
        "lines": lines,
    }


def qualification_actions_run(observation: dict[str, Any]) -> tuple[str, str]:
    workflows = [item for item in observation["workflows"] if item.get("name") == "qualification"]
    if len(workflows) != 1:
        raise SystemExit("CI observation must name exactly one qualification workflow")
    workflow = workflows[0]
    match = re.search(r"/runs/(\d+)", str(workflow.get("url", "")))
    if not match:
        raise SystemExit("qualification workflow URL has no actions run id")
    run_id = match.group(1)
    commit = str(observation["public_commit"])
    jobs = workflow.get("jobs") or []
    if not jobs or any(job.get("conclusion") != "success" for job in jobs):
        raise SystemExit("qualification workflow jobs did not all succeed")
    if any(job.get("head_sha") != commit for job in jobs):
        raise SystemExit("qualification jobs do not all record public_commit")
    panel_steps = [
        step.get("name")
        for job in jobs
        for step in job.get("steps") or []
        if str(step.get("name", "")).startswith("Execute the ")
    ]
    if not any("sf-csa" in name for name in panel_steps):
        raise SystemExit("qualification observation does not record an SF-CSA execution step")
    return run_id, commit


def qualification_table(summary: dict[str, Any], manifest: dict[str, Any], recon: dict[str, Any], observation: dict[str, Any]) -> dict[str, Any]:
    if manifest.get("collection_version") != "2.11" or recon["manifest"]["collection_version"] != "2.11":
        raise SystemExit("qualification figure expects collection 2.11")
    if recon["manifest"]["sha256"] != sha256(REPO / "evidence/benchmarks/qualification-v2/PANEL_MANIFEST.json"):
        raise SystemExit("reconciliation panel-manifest checksum does not match the file")
    adopted = {panel["workflow"]: panel["adopted_records"] for panel in recon["panels"]}
    by_workflow = {panel["workflow"]: panel for panel in summary["panels"]}
    actions_run, commit = qualification_actions_run(observation)
    order = (
        ("structure_qc", "StructQC"),
        ("conformational_state", "ABL state"),
        ("functional_site_state", "SiteContext"),
        ("assembly_interface", "AssemblyContext"),
        ("sf_csa", "SF-CSA"),
        ("membrane_orientation", "MembraneOrient"),
    )
    rows = []
    executed_passed = 0
    for workflow, label in order:
        if workflow == "sf_csa":
            if workflow not in summary["workflows_not_executed"]:
                raise SystemExit("retained summary unexpectedly executed SF-CSA")
            retained = "Not executed"
            adopted_cell = f"{adopted[workflow]}/{adopted[workflow]}"
            ci_cell = "Included in passed blocking aggregate"
        elif workflow == "membrane_orientation":
            panel = by_workflow[workflow]
            if panel["cases_adopted"] != adopted[workflow] or panel["release_blocking"] is not False:
                raise SystemExit("membrane panel does not match the retained nonblocking receipt")
            retained = f"{panel['cases']['passed']}/{panel['cases']['total']}; {panel['cases']['failed']} failed"
            adopted_cell = f"{adopted[workflow]}/{panel['cases_required']}; beta only"
            ci_cell = "Job completed; research-only, nonblocking"
            executed_passed += panel["cases"]["passed"]
        else:
            panel = by_workflow[workflow]
            controls = panel["controls"]
            retained = f"{panel['cases']['passed']}/{panel['cases']['total']}"
            if controls["total"]:
                noun = "control" if controls["total"] == 1 else "controls"
                retained = f"{retained}; {controls['passed']} {noun}"
            adopted_cell = f"{adopted[workflow]}/{adopted[workflow]}"
            ci_cell = "Included in passed blocking aggregate"
            executed_passed += panel["cases"]["passed"]
        rows.append([label, adopted_cell, retained, ci_cell])
    if executed_passed != summary["cases_passed"] or summary["cases_required"] != recon["required_records"]:
        raise SystemExit("retained pass count does not add up to the execution summary")
    if sum(adopted.values()) != recon["adopted_records"]:
        raise SystemExit("adopted panel counts do not add up to the reconciliation total")
    short_commit = commit[:12]
    footnotes = [
        (
            f"Collection v{manifest['collection_version']} adopts {recon['adopted_records']} of "
            f"{recon['required_records']} records. The retained {summary['cases_passed']}/"
            f"{summary['cases_required']} summary predates this reconciliation and omits SF-CSA execution."
        ),
        (
            f"Run {actions_run} passed its configured blocking gate on commit {short_commit}; "
            "its job metadata is not per-case scientific evidence."
        ),
        (
            f"Membrane coverage is beta-barrel only ({adopted['membrane_orientation']}/"
            f"{by_workflow['membrane_orientation']['cases_required']} adopted); execution remains "
            "research-only and nonblocking. No unified accuracy score or changed-candidate qualification is claimed."
        ),
    ]
    header = [
        "Workflow",
        "Adopted records\n(v2.11)",
        "Retained summary\n(single machine)",
        f"Run {actions_run}\n(configured gate)",
    ]
    return {
        "actions_run_id": actions_run,
        "actions_run_is_store_run": False,
        "public_commit": commit,
        "header": header,
        "rows": rows,
        "footnotes": footnotes,
        "label_note": (
            "The previously committed figure called this column Public CI #78. "
            "CI_OBSERVATION.json does not contain that number. The figure now uses "
            f"the qualification actions run id recorded there, {actions_run}."
        ),
    }


def draw_qualification(plt, dest: Path, table: dict[str, Any]) -> dict[str, Any]:
    fig, ax = plt.subplots(num=3, figsize=(12.8, 5.5))
    ax.axis("off")
    ax.set_xlim(-0.12, 12.8)
    ax.set_ylim(-1.03, 5.5)
    ax.text(0, 5.25, "Three records, three different questions", fontsize=17, fontweight="bold", color=NAVY)
    ax.text(0, 4.96, "Adopted panel coverage · retained single-machine execution · later public CI on a named commit", fontsize=10.5, color=NAVY)
    drawn = ax.table(
        cellText=table["rows"],
        colLabels=table["header"],
        loc="upper left",
        bbox=[0, 0.16, 1, 0.72],
        cellLoc="left",
        colLoc="left",
        colWidths=[0.20, 0.18, 0.23, 0.39],
    )
    drawn.auto_set_font_size(False)
    drawn.set_fontsize(9.5)
    for (row_index, _column), cell in drawn.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_facecolor(NAVY if row_index == 0 else "#edf4f2" if row_index % 2 else "#f7f9f8")
        cell.set_text_props(color="white" if row_index == 0 else NAVY)
    for index, footnote in enumerate(table["footnotes"]):
        ax.text(0, -0.00 - index * 0.13, footnote, transform=ax.transAxes, fontsize=9, color=GOLD if index == 2 else NAVY)
    save_figure(plt, fig, dest, "qualification-history")
    return {
        "actions_run_id": table["actions_run_id"],
        "actions_run_is_store_run": False,
        "public_commit": table["public_commit"],
        "header": table["header"],
        "rows": table["rows"],
        "footnotes": table["footnotes"],
        "label_note": table["label_note"],
    }


def draw_questions(plt, FancyBboxPatch, dest: Path) -> dict[str, Any]:
    fig, ax = plt.subplots(num=4, figsize=(11, 4.6))
    ax.axis("off")
    ax.set_xlim(-0.14, 11.12)
    ax.set_ylim(-0.95, 4.6)
    ax.text(0, 4.28, "Six structural evidence questions", fontsize=18, fontweight="bold", color=NAVY)
    ax.text(0, 3.98, "One gated manifest, five independent analyses", fontsize=10.5, color=NAVY)
    ax.add_patch(FancyBboxPatch((0, 0.58), 2.4, 3.27, boxstyle="round,pad=.06", facecolor="#edf4f2", edgecolor=TEAL))
    trust_lines = [
        (3.42, "Coordinate trust", {"fontweight": "bold", "fontsize": 12}),
        (3.14, "can I trust this file?", {"fontsize": 10}),
        (2.86, "structqc", {"fontsize": 10, "color": TEAL, "family": "DejaVu Sans Mono"}),
        (2.36, "identity · provenance", {"fontsize": 9.5}),
        (2.14, "SHA-256 · residue map", {"fontsize": 9.5}),
        (1.72, "provenance is recorded,", {"fontsize": 9, "style": "italic", "color": GOLD}),
        (1.52, "never inferred", {"fontsize": 9, "style": "italic", "color": GOLD}),
        (1.08, "no downstream module may", {"fontsize": 8.5}),
        (0.90, "interpret unbound coordinates", {"fontsize": 8.5}),
    ]
    for y_value, text, style in trust_lines:
        color = style.pop("color", NAVY)
        ax.text(1.2, y_value, text, ha="center", color=color, **style)
    ax.plot([0.35, 2.05], [2.62, 2.62], color="#b9c8c4", lw=0.8)
    scopes = [
        ("Membrane orientation", "memorient", "which way is out?", "geometry, not intact-cell exposure"),
        ("Conformational state", "state-atlas", "active-like or inactive-like?", "resemblance, not activity"),
        ("Functional site", "site-context · actstate", "is the site intact, and competent?", "annotation, not observed catalysis"),
        ("Assembly interface", "assembly-context", "what does it touch, and how much is buried?", "one coordinate state, not affinity"),
        ("Comparative analysis", "sf-csa", "what else looks like this — by fold, and by sequence?", "a shared fold is not a function"),
    ]
    ax.plot([2.4, 3.05], [2.215, 2.215], color=NAVY, lw=0.9)
    ax.plot([3.05, 3.05], [0.855, 3.575], color=NAVY, lw=0.9)
    renderer = fig.canvas.get_renderer()
    for index, (title, cli, question, limit) in enumerate(scopes):
        bottom = 3.30 - index * 0.68
        center = bottom + 0.275
        ax.add_patch(FancyBboxPatch((3.38, bottom), 7.62, 0.55, boxstyle="round,pad=.04", facecolor="#f7f9f8", edgecolor="#d3d9de"))
        ax.annotate("", xy=(3.36, center), xytext=(3.05, center), arrowprops={"arrowstyle": "->", "color": NAVY, "lw": 0.9})
        title_artist = ax.text(3.58, center + 0.10, title, fontsize=11.5, fontweight="bold", color=NAVY, va="center")
        ax.text(
            title_artist.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).x1 + 0.18,
            center + 0.10,
            cli,
            fontsize=8.5,
            color="#5f7570",
            family="DejaVu Sans Mono",
            va="center",
        )
        ax.text(3.58, center - 0.15, question, fontsize=9.5, color=NAVY, va="center")
        ax.text(10.85, center, limit, fontsize=8.5, style="italic", color=GOLD, ha="right", va="center")
    ax.add_patch(FancyBboxPatch((0, -0.52), 11, 0.44, boxstyle="round,pad=.02", facecolor="#edf4f2", edgecolor=TEAL))
    banner = "Evidence is reported separately — no combined score"
    ax.text(5.5, -0.30, banner, fontsize=12, fontweight="bold", color=NAVY, ha="center", va="center")
    footer = (
        "The five analyses are independently installable and consume the same checksum-bound manifest.\n"
        "Missing references, mappings or runtimes stay visible as missing, and no analysis feeds another a verdict."
    )
    ax.text(0, -0.78, footer, fontsize=9, color=GOLD, va="top")
    save_figure(plt, fig, dest, "six-evidence-questions")
    return {
        "coordinate_trust": [text for _y, text, _style in trust_lines],
        "scopes": [{"title": title, "cli": cli, "question": question, "limit": limit} for title, cli, question, limit in scopes],
        "banner": banner,
        "footer": footer,
    }


def coverage_from_evidence(evidence: dict[str, Any]) -> dict[str, Any]:
    completeness = evidence["completeness"]
    reference_length = int(completeness["reference_length"])
    mapped = int(completeness["mapped_residues"])
    present = {
        int(residue["sequence_index"])
        for residue in evidence["residues"]
        if residue.get("sequence_index") not in (None, "")
    }
    missing = [index for index in range(1, reference_length + 1) if index not in present]
    coverage = mapped / reference_length
    if abs(coverage - float(completeness["coverage_fraction"])) > 1e-5:
        raise SystemExit("coverage_fraction does not match mapped residues over reference length")
    b_factors = [float(residue["mean_b_factor"]) for residue in evidence["residues"] if residue.get("mean_b_factor") not in (None, "")]
    if not b_factors:
        raise SystemExit("residue evidence has no mean B factors")
    return {
        "chain": evidence["coordinate"].get("selected_chain"),
        "model": evidence["coordinate"].get("selected_model"),
        "structure_sha256": evidence["coordinate"].get("sha256"),
        "reference_id": evidence["reference"]["id"],
        "mapped_residues": mapped,
        "reference_length": reference_length,
        "coverage_fraction": f"{float(completeness['coverage_fraction']):.6f}",
        "identity_fraction": f"{float(completeness['identity_fraction']):.3f}",
        "missing_reference_positions": missing,
        "missing_label": position_ranges(missing),
        "median_mean_b_factor": f"{statistics.median(b_factors):.3f}",
        "external_validation_state": evidence["external_validation"]["state"],
        "pae_state": evidence["pae"]["state"],
    }


def draw_ca2(plt, dest: Path, runs: list[dict[str, Any]], case_id: str, lock: dict[str, Any]) -> dict[str, Any]:
    if lock.get("case_id") != case_id:
        raise SystemExit(f"source lock case_id is {lock.get('case_id')!r}, figure expects {case_id}")
    locked = {entry["entry_id"]: entry["sha256"] for entry in lock["structures"]}
    measured = []
    for run in runs:
        evidence = load_json(artifact_path(run, "structure_evidence"))
        analysis = load_json(artifact_path(run, "analysis_run"))
        metrics = coverage_from_evidence(evidence)
        entry_id = run["entry_id"]
        if metrics["structure_sha256"] != locked.get(entry_id):
            raise SystemExit(f"{entry_id} evidence checksum does not match the public source lock")
        if metrics["reference_id"] != lock["reference"]["id"]:
            raise SystemExit(f"{entry_id} reference is not {lock['reference']['id']}")
        if analysis["status"] != "scientifically_incomplete":
            raise SystemExit(f"{entry_id} run status is {analysis['status']!r}; this case supplies no external validation")
        metrics.update({
            "entry_id": entry_id,
            "analysis_id": run["analysis_id"],
            "run_id": run["run_id"],
            "run_status": analysis["status"],
        })
        measured.append(metrics)
    fig, ax = plt.subplots(num=5, figsize=(11, 4.8))
    labels = [f"{item['entry_id']}\nchain {item['chain']}" for item in measured]
    values = [float(item["coverage_fraction"]) for item in measured]
    ax.bar(labels, values, color=TEAL, width=0.55, zorder=2)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("Coverage of UniProt P00918")
    ax.set_title("Public CA II structures · sequence coverage", fontsize=16, fontweight="bold", color=NAVY, loc="left")
    ax.axhline(1, color="#d3d9de", lw=0.8, zorder=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=NAVY)
    for index, (item, value) in enumerate(zip(measured, values)):
        ax.text(
            index,
            value + 0.03,
            f"{item['mapped_residues']}/{item['reference_length']}",
            ha="center",
            va="bottom",
            color=NAVY,
            fontsize=11,
        )
    detail_lines = []
    for item in measured:
        detail_lines.append(
            f"{item['entry_id']} chain {item['chain']}: identity {item['identity_fraction']}, "
            f"missing {item['missing_label']}, median residue mean B {item['median_mean_b_factor']}"
        )
    limits = [
        "No external validation report was supplied. Both store runs are scientifically incomplete.",
        "Coverage describes these deposited coordinates. It is not a ligand, affinity, or activity result.",
        "B factors are deposited descriptors from separate experiments, not a shared quality score.",
        "Case " + case_id + ". Store runs: " + "; ".join(f"{item['analysis_id']} {item['run_id']}" for item in measured) + ".",
    ]
    fig.subplots_adjust(bottom=0.58)
    fig.text(0.125, 0.44, "\n".join(detail_lines), color=NAVY, fontsize=10, va="top")
    fig.text(0.125, 0.28, "\n".join(limits), color=GOLD, fontsize=9, va="top")
    save_figure(plt, fig, dest, "ca2-coverage")
    return {"case_id": case_id, "entries": measured, "limits": limits}


def public_artifact(artifact: dict[str, Any]) -> dict[str, str]:
    return {key: artifact[key] for key in ("role", "path", "sha256") if key in artifact} | (
        {"run_id": artifact["run_id"], "analysis_id": artifact["analysis_id"]} if "run_id" in artifact else {}
    )


def render(dest: Path, manifest: dict[str, Any], runs: dict[str, Path], inputs_path: Path) -> dict[str, Any]:
    plt, FancyBboxPatch = configure_matplotlib()
    plt.close("all")
    specs = manifest["figures"]
    resolved_runs: dict[str, list[dict[str, Any]]] = {}
    resolved_files: dict[str, list[dict[str, str]]] = {}
    for name in FIGURE_ORDER:
        spec = specs[name]
        if spec["input_class"] == "store_run":
            resolved_runs[name] = [resolve_run(run, runs) for run in spec["runs"]]
        elif spec["input_class"] == "historical_receipt":
            resolved_files[name] = [resolve_repo_artifact(artifact) for artifact in spec["artifacts"]]
        elif spec["input_class"] != "schematic":
            raise SystemExit(f"unknown input class for {name}: {spec['input_class']}")
    historical = {
        key: {
            **value,
            "artifacts": [resolve_repo_artifact(artifact) for artifact in value["artifacts"]],
        }
        for key, value in manifest["historical_receipts"].items()
    }

    plotted = {
        "architecture": draw_architecture(plt, FancyBboxPatch, dest),
        "worked-example": draw_worked_example(plt, dest, resolved_runs["worked-example"][0]),
        "qualification-history": draw_qualification(
            plt,
            dest,
            qualification_table(
                load_json(REPO / next(item["path"] for item in resolved_files["qualification-history"] if item["role"] == "execution_summary")),
                load_json(REPO / next(item["path"] for item in resolved_files["qualification-history"] if item["role"] == "panel_manifest")),
                load_json(REPO / next(item["path"] for item in resolved_files["qualification-history"] if item["role"] == "qualification_reconciliation")),
                load_json(REPO / next(item["path"] for item in resolved_files["qualification-history"] if item["role"] == "public_ci_observation")),
            ),
        ),
        "six-evidence-questions": draw_questions(plt, FancyBboxPatch, dest),
    }
    ca2_spec = specs["ca2-coverage"]
    lock_spec = resolve_repo_artifact(ca2_spec["source_lock"])
    plotted["ca2-coverage"] = draw_ca2(
        plt,
        dest,
        resolved_runs["ca2-coverage"],
        ca2_spec["case_id"],
        load_json(REPO / lock_spec["path"]),
    )

    figures: dict[str, Any] = {}
    for name in FIGURE_ORDER:
        spec = specs[name]
        runs = resolved_runs.get(name, [])
        artifacts = [public_artifact(artifact) for run in runs for artifact in run["artifacts"]]
        artifacts.extend(resolved_files.get(name, []))
        if name == "ca2-coverage":
            artifacts.append(lock_spec)
        note = spec["reason"]
        if name == "qualification-history":
            note = note + " " + plotted[name]["label_note"]
        if not cited_in_paper(name):
            note = note + " paper/paper.md does not reference this figure. This generator does not add or remove that reference."
        figures[name] = {
            "input_class": spec["input_class"],
            "run_ids": [run["run_id"] for run in runs],
            "cited_in_paper_md": cited_in_paper(name),
            "note": note,
            "artifacts": artifacts,
            "plotted_data": plotted[name],
            "files": {
                f"{name}.png": sha256(dest / f"{name}.png"),
                f"{name}.svg": sha256(dest / f"{name}.svg"),
            },
        }
    provenance = {
        "schema_version": "2.0",
        "figure_generator": {"path": "tools/build_scientific_figures.py", "sha256": sha256(GENERATOR)},
        "figure_inputs": {
            "path": inputs_path.resolve().relative_to(REPO).as_posix(),
            "sha256": sha256(inputs_path),
        },
        "matplotlib_determinism": {
            "backend": "Agg",
            "svg.fonttype": "none",
            "svg.hashsalt": SVG_HASHSALT,
            "SOURCE_DATE_EPOCH": SOURCE_DATE_EPOCH,
            "required_package": MATPLOTLIB_SPEC,
            "png_bytes": PNG_DIFFERENCE_REASON,
        },
        "historical_receipts": historical,
        "figures": figures,
        "limits": (
            "Worked example is the shipped two-residue synthetic fixture. The builder regenerates that "
            "store run and the public CA2 runs, then requires the run ids and artifact SHA-256s in "
            "FIGURE_INPUTS.json. Those run directories are not committed. The 2026-09-06 evidence files "
            "are retained as the historical receipt and are not a run id. CA2 coverage is the public case "
            "public-ca2-structure-coverage-v1; both runs are scientifically incomplete because no external "
            "validation report was supplied. The coordinate files are acquired at build time and are not "
            "committed. Qualification history presents three non-combinable records: adopted v2.11 panel "
            "counts, the retained single-machine summary, and public CI job metadata. That metadata is an "
            "actions run id, not a store run id, and it is not per-case evidence. Membrane orientation stays "
            "beta-barrel-only, research-only, and nonblocking. No unified accuracy score or "
            "changed-candidate qualification is claimed."
        ),
    }
    dump_json(dest / PROVENANCE_NAME, provenance)
    return provenance


def svg_canonical(text: str) -> str:
    text = re.sub(r"<dc:date>.*?</dc:date>", "", text)
    text = re.sub(r"<dc:title>.*?</dc:title>", "", text)
    text = re.sub(r"url\(#[^)]+\)", "url(#clip)", text)
    text = re.sub(r"\sid=\"[^\"]*\"", "", text)
    return "\n".join(line.rstrip() for line in text.splitlines()) + "\n"


def compare_provenance(rebuilt: dict[str, Any], committed: dict[str, Any]) -> list[str]:
    problems = []
    if rebuilt["schema_version"] != committed.get("schema_version"):
        problems.append(f"provenance schema {committed.get('schema_version')!r} is not {rebuilt['schema_version']}")
    if rebuilt["figure_generator"] != committed.get("figure_generator"):
        problems.append("FIGURE_PROVENANCE.json generator hash does not match tools/build_scientific_figures.py")
    if rebuilt["figure_inputs"] != committed.get("figure_inputs"):
        problems.append("FIGURE_PROVENANCE.json does not cite the current FIGURE_INPUTS.json checksum")
    if rebuilt["matplotlib_determinism"] != committed.get("matplotlib_determinism"):
        problems.append("matplotlib determinism record in FIGURE_PROVENANCE.json is stale")
    if rebuilt["historical_receipts"] != committed.get("historical_receipts"):
        problems.append("historical receipt checksums in FIGURE_PROVENANCE.json do not match the frozen inputs")
    if rebuilt["limits"] != committed.get("limits"):
        problems.append("FIGURE_PROVENANCE.json limits text is stale")
    committed_figures = committed.get("figures") or {}
    for name, figure in rebuilt["figures"].items():
        other = committed_figures.get(name)
        if other is None:
            problems.append(f"FIGURE_PROVENANCE.json has no entry for {name}")
            continue
        for field in ("input_class", "run_ids", "cited_in_paper_md", "note", "artifacts", "plotted_data"):
            if figure[field] != other.get(field):
                problems.append(f"{name} {field} in FIGURE_PROVENANCE.json does not match this rebuild")
    extra = set(committed_figures) - set(rebuilt["figures"])
    if extra:
        problems.append("FIGURE_PROVENANCE.json has unexpected figures: " + ", ".join(sorted(extra)))
    return problems


def compare_files(rebuilt_dir: Path, committed_dir: Path, rebuilt: dict[str, Any]) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    notes: list[str] = []
    for name in FIGURE_ORDER:
        for suffix in (".png", ".svg"):
            filename = name + suffix
            rebuilt_bytes = (rebuilt_dir / filename).read_bytes()
            committed_path = committed_dir / filename
            if not committed_path.is_file():
                failures.append(f"committed {filename} is missing")
                continue
            committed_bytes = committed_path.read_bytes()
            recorded_hash = rebuilt["figures"][name]["files"][filename]
            if hashlib.sha256(rebuilt_bytes).hexdigest() != recorded_hash:
                failures.append(f"rebuilt provenance hash does not match rebuilt {filename}")
            if rebuilt_bytes == committed_bytes:
                continue
            if suffix == ".png":
                svg_rebuilt = (rebuilt_dir / f"{name}.svg").read_text(encoding="utf-8")
                svg_committed = (committed_dir / f"{name}.svg").read_text(encoding="utf-8")
                if svg_rebuilt == svg_committed or svg_canonical(svg_rebuilt) == svg_canonical(svg_committed):
                    notes.append(f"{filename} bytes differ; SVG matches. {PNG_DIFFERENCE_REASON}")
                    continue
                notes.append(
                    f"{filename} and {name}.svg bytes differ. Extracted plotted data matches the "
                    f"committed provenance, so the check accepts the platform difference. {PNG_DIFFERENCE_REASON}"
                )
                continue
            svg_rebuilt = rebuilt_bytes.decode("utf-8")
            svg_committed = committed_bytes.decode("utf-8")
            if svg_canonical(svg_rebuilt) == svg_canonical(svg_committed):
                notes.append(f"{filename} differs only in SVG metadata or clip ids, which the fixed hashsalt and SOURCE_DATE_EPOCH are there to stop. Investigate if this persists on one machine.")
                continue
            notes.append(f"{filename} drawing bytes differ; extracted plotted data was compared separately.")
    return failures, notes


def check(manifest: dict[str, Any], runs: dict[str, Path]) -> int:
    first = Path(tempfile.mkdtemp(prefix="yauvi-figures-"))
    second = Path(tempfile.mkdtemp(prefix="yauvi-figures-"))
    try:
        rebuilt = render(first, manifest, runs, DEFAULT_INPUTS)
        render(second, manifest, runs, DEFAULT_INPUTS)
        drifted = []
        for path in sorted(first.rglob("*")):
            if not path.is_file():
                continue
            other = second / path.relative_to(first)
            if not other.is_file() or path.read_bytes() != other.read_bytes():
                drifted.append(path.name)
        if drifted:
            print("Two rebuilds on this machine were not byte-identical: " + ", ".join(drifted), file=sys.stderr)
            return 1
        committed = load_json(FIGURES_DIR / PROVENANCE_NAME)
        problems = compare_provenance(rebuilt, committed)
        file_failures, notes = compare_files(first, FIGURES_DIR, rebuilt)
        problems.extend(file_failures)
        if problems:
            print("\n".join(problems), file=sys.stderr)
            return 1
        if notes:
            print("Figure byte notes:")
            for note in notes:
                print(f"  {note}")
        else:
            print("Paper figures match the committed PNG and SVG bytes.")
        print("Extracted plotted data matches FIGURE_PROVENANCE.json.")
        return 0
    finally:
        shutil.rmtree(first, ignore_errors=True)
        shutil.rmtree(second, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inputs", type=Path, default=DEFAULT_INPUTS, help="Figure-inputs manifest. Default: paper/figures/FIGURE_INPUTS.json")
    parser.add_argument("--check", action="store_true", help="Regenerate the pinned runs, rebuild, and compare with the committed figures. Exit 1 if the plotted data differs.")
    args = parser.parse_args(argv)
    inputs = args.inputs if args.inputs.is_absolute() else REPO / args.inputs
    if not inputs.is_file():
        print(f"missing figure inputs: {inputs}", file=sys.stderr)
        return 2
    manifest = load_json(inputs)
    if manifest.get("schema_version") != "1":
        print(f"unsupported figure-inputs schema {manifest.get('schema_version')!r}", file=sys.stderr)
        return 2
    if args.check and inputs.resolve() != DEFAULT_INPUTS.resolve():
        print("--check compares the committed provenance, so it uses paper/figures/FIGURE_INPUTS.json.", file=sys.stderr)
        return 2
    scratch, runs = prepare_runs(manifest)
    try:
        if args.check:
            return check(manifest, runs)
        provenance = render(FIGURES_DIR, manifest, runs, inputs.resolve())
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    print(f"Wrote {FIGURES_DIR / PROVENANCE_NAME}")
    for name, figure in provenance["figures"].items():
        kind = figure["input_class"]
        cited = ", ".join(figure["run_ids"]) if figure["run_ids"] else "no store run id"
        print(f"  {name}: {kind}; {cited}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
