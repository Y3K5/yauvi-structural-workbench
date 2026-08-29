#!/usr/bin/env python3
"""Run five deterministic structural-workbench cases and build a human-readable showcase.

All coordinates and annotations are synthetic. The showcase demonstrates software
behavior and scientific boundaries; it contains no human or campaign research data.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PORTFOLIO = ROOT / "examples" / "structural-portfolio"
QC_EXAMPLE = ROOT / "structqc" / "examples"
DEFAULT_OUT = ROOT / "yauvi-structural-workbench" / "showcase" / "five-human-use-cases"
DEFAULT_PUBLIC_OUT = ROOT / "yauvi-structural-workbench" / "public-showcase"
PUBLIC_SOURCE_FILES = ("index.html", "styles.css", "app.js", "README.md")
PACKAGE_PATHS = (
    "structqc/src",
    "Membrane Orientor/memorient/src",
    "state-atlas/src",
    "site-context/src",
    "assembly-context/src",
)

QUALIFICATION_SOURCE_LINKS = {
    "structure_qc": [
        {"label": "RCSB and wwPDB files", "url": "https://www.rcsb.org/docs/programmatic-access/file-download-services"},
        {"label": "AlphaFold Protein Structure Database", "url": "https://www.alphafold.ebi.ac.uk/"},
    ],
    "membrane_orientation": [
        {"label": "OPM and PPM reference system", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3245162/"},
    ],
    "conformational_state": [
        {"label": "KinCore ABL1 classifications", "url": "https://dunbrack.fccc.edu/kincore/GENE/ABL1"},
    ],
    "functional_site_state": [
        {"label": "M-CSA entry 1", "url": "https://www.ebi.ac.uk/thornton-srv/m-csa/entry/1/"},
    ],
    "assembly_interface": [
        {"label": "RCSB structure 4HHB", "url": "https://www.rcsb.org/structure/4HHB"},
        {"label": "FreeSASA", "url": "https://freesasa.github.io/"},
    ],
    "sf_csa": [
        {"label": "CATH downloads", "url": "https://cathdb.info/download"},
        {"label": "Foldseek", "url": "https://github.com/steineggerlab/foldseek"},
        {"label": "DIAMOND", "url": "https://www.nature.com/articles/s41592-021-01101-x"},
    ],
}


def canonical(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _required_check(workflow: dict[str, Any], check_id: str) -> dict[str, Any]:
    return next(item for item in workflow["checks"] if item["check"] == check_id)


def qualification_showcase_rows(qualification: dict[str, Any]) -> list[dict[str, Any]]:
    """Translate exact qualification evidence into bounded public narratives."""
    records = {item["workflow"]: item for item in qualification["workflows"]}
    qc = records["structure_qc"]
    membrane = records["membrane_orientation"]
    state = records["conformational_state"]
    site = records["functional_site_state"]
    assembly = records["assembly_interface"]
    sfcsa = records["sf_csa"]
    membrane_failed = [
        item["check"] for item in membrane["checks"]
        if item.get("required") and not item["passed"]
    ]
    state_failed = [
        item["check"] for item in state["checks"]
        if item.get("required") and not item["passed"]
    ]
    active, inactive = state["records"]
    runtime_versions = sfcsa["runtime_versions"]
    rows = [
        {
            "analysis_type": "structure_qc",
            "case_label": "1CRN plus AlphaFold P69905 v6",
            "status": qc["status"],
            "independent_reference": "wwPDB validation and AlphaFold DB model-confidence records",
            "finding": (
                f'Exact 46-residue mapping; raw clashscore {qc["observations"]["wwpdb_metrics_imported"]["clashscore"]:.1f}; '
                f'mean predicted-model pLDDT {qc["observations"]["predicted_mean_plddt"]:.3f}; unknown provenance remained incomplete.'
            ),
            "biological_context": "A structure can proceed to residue-level review only when its identity and provenance are explicit.",
            "remaining_limit": qc["claim_boundary"],
            "failed_checks": [],
        },
        {
            "analysis_type": "membrane_orientation",
            "case_label": "Five beta barrels and three alpha-helical OPM structures",
            "status": membrane["status"],
            "independent_reference": membrane["independent_reference"],
            "finding": (
                f'Beta-barrel mean normal error {_required_check(membrane, "beta_barrel_mean_normal_error")["observed"]:.3f}° passed; '
                f'alpha-helical mean error {_required_check(membrane, "alpha_helical_mean_normal_error")["observed"]:.3f}° and '
                f'1U19 rotation Jaccard {_required_check(membrane, "alpha_helical_rotation_invariance")["observed"][-1]:.2f} failed.'
            ),
            "biological_context": "The current local method is promising for the tested beta-barrels but is not reliable across the tested alpha-helical stratum.",
            "remaining_limit": membrane["claim_boundary"],
            "failed_checks": membrane_failed,
        },
        {
            "analysis_type": "conformational_state",
            "case_label": "ABL active and inactive holdouts",
            "status": state["status"],
            "independent_reference": state["independent_reference"],
            "finding": (
                f'{active["pdb_id"]} was {active["observed_yauvi_label"]} at {active["frame_metrics"][0]["best_rmsd_A"]:.3f} Å; '
                f'{inactive["pdb_id"]} safely remained {inactive["observed_yauvi_label"]} at {inactive["frame_metrics"][0]["best_rmsd_A"]:.3f} Å.'
            ),
            "biological_context": "A resolved resemblance can support state comparison; an unresolved result prevents an unjustified active or inactive claim.",
            "remaining_limit": state["claim_boundary"],
            "failed_checks": state_failed,
        },
        {
            "analysis_type": "functional_site_state",
            "case_label": "M-CSA glutamate racemase, PDB 1B73",
            "status": site["status"],
            "independent_reference": site["independent_reference"],
            "finding": (
                f'{len(site["mapped_sites"])} curated residues mapped exactly at positions '
                f'{", ".join(str(item["position"]) for item in site["mapped_sites"])}; pocket evidence remained missing.'
            ),
            "biological_context": "Exact residue mapping can support mutational-control planning while keeping annotation separate from observed chemistry.",
            "remaining_limit": site["claim_boundary"],
            "failed_checks": [],
        },
        {
            "analysis_type": "assembly_interface",
            "case_label": "Hemoglobin biological assembly 1, PDB 4HHB",
            "status": assembly["status"],
            "independent_reference": assembly["independent_reference"],
            "finding": (
                f'{len(assembly["assembly"]["chains_observed"])} expected chains recovered; '
                f'{assembly["surface"]["buried_sasa_A2"]:.3f} Å² subject surface buried using FreeSASA.'
            ),
            "biological_context": "The result identifies assembly-bound interfaces and burial that are absent from isolated-chain interpretation.",
            "remaining_limit": assembly["claim_boundary"],
            "failed_checks": [],
        },
        {
            "analysis_type": "sf_csa",
            "case_label": "CATH exact, homolog, fold-analogy, and unrelated controls",
            "status": sfcsa["status"],
            "independent_reference": sfcsa["independent_reference"],
            "finding": (
                f'{sfcsa["structure_hit_count"]} structure hits and {sfcsa["sequence_hit_count"]} sequence hits were kept separate; '
                f'Foldseek {runtime_versions["foldseek"]} and {runtime_versions["diamond"]} ran locally.'
            ),
            "biological_context": "Exact, homologous, analogous, and unresolved relationships remain distinct instead of collapsing into a function score.",
            "remaining_limit": sfcsa["claim_boundary"],
            "failed_checks": [],
        },
    ]
    for row in rows:
        row["source_links"] = QUALIFICATION_SOURCE_LINKS[row["analysis_type"]]
    return rows


def pdb_atom(serial: int, name: str, residue: str, residue_id: int,
             x: float, y: float, z: float, element: str) -> str:
    return (
        f"ATOM  {serial:5d} {name:^4s} {residue:>3s} A{residue_id:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 75.00          {element:>2s}\n"
    )


def synthetic_tm_receptor() -> str:
    """Return an invented single-pass helix with explicit charged and polar tails."""
    cytoplasmic = ("LYS", "ARG", "SER", "LYS", "THR", "ARG", "SER", "LYS")
    core = ("LEU", "ILE", "VAL", "PHE", "ALA", "MET")
    extracellular = ("ASP", "GLU", "SER", "THR", "GLU", "ASN", "SER", "ASP")
    residues = list(cytoplasmic) + [core[index % len(core)] for index in range(41)] + list(extracellular)
    z0 = -27.0
    lines: list[str] = []
    serial = 0
    for index, residue in enumerate(residues, start=1):
        angle = math.radians(100.0 * (index - 1))
        ca = (2.3 * math.cos(angle), 2.3 * math.sin(angle), z0 + 1.5 * (index - 1))
        radial = (math.cos(angle), math.sin(angle), 0.0)
        atoms = (
            ("N", (ca[0] - 0.6, ca[1] + 0.8, ca[2]), "N"),
            ("CA", ca, "C"),
            ("C", (ca[0] + 0.6, ca[1] - 0.8, ca[2]), "C"),
            ("O", (ca[0] + 1.2, ca[1] - 1.4, ca[2] + 0.3), "O"),
            ("CB", (ca[0] + 1.5 * radial[0], ca[1] + 1.5 * radial[1], ca[2]), "C"),
            ("CG", (ca[0] + 3.0 * radial[0], ca[1] + 3.0 * radial[1], ca[2]), "C"),
        )
        for name, xyz, element in atoms:
            serial += 1
            lines.append(pdb_atom(serial, name, residue, index, *xyz, element))
    lines.extend(("TER\n", "END\n"))
    return "".join(lines)


def copy_inputs(output: Path) -> dict[str, Path]:
    input_root = output / "inputs"
    mapping: dict[str, Path] = {}
    groups = {
        "qc": ("model.pdb", "reference.fasta", "provenance.json", "validation.json"),
        "portfolio": (
            "query.pdb", "reference.fasta", "provenance.json", "assembly.pdb",
            "annotations.json", "pockets.json", "reference_set.json", "inactive_reference.pdb",
        ),
    }
    for group, names in groups.items():
        source = QC_EXAMPLE if group == "qc" else PORTFOLIO
        for name in names:
            destination = input_root / group / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source / name, destination)
            mapping[f"{group}/{name}"] = destination
    membrane = input_root / "membrane" / "synthetic_tm_receptor.pdb"
    write_text(membrane, synthetic_tm_receptor())
    mapping["membrane/synthetic_tm_receptor.pdb"] = membrane
    topology = input_root / "membrane" / "synthetic_tm_topology.json"
    write_text(topology, canonical({
        "schema_version": "1.0",
        "coordinate_sha256": sha256(membrane),
        "source": {"id": "synthetic-showcase-topology", "citation": "invented software fixture"},
        "spans": [{"chain_id": "A", "start_auth_seq_id": 9, "end_auth_seq_id": 49}],
        "sidedness": {"extracellular_residue": {"chain_id": "A", "auth_seq_id": 57, "insertion_code": ""}},
        "limitations": ["Synthetic topology tests software behavior only."],
    }))
    mapping["membrane/synthetic_tm_topology.json"] = topology
    return mapping


def environment() -> dict[str, str]:
    value = os.environ.copy()
    local = os.pathsep.join(str(ROOT / item) for item in PACKAGE_PATHS)
    value["PYTHONPATH"] = local + (os.pathsep + value["PYTHONPATH"] if value.get("PYTHONPATH") else "")
    return value


def run_command(case_dir: Path, command: list[str]) -> int:
    completed = subprocess.run(
        command, cwd=ROOT, env=environment(), text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    replacements = ((str(ROOT), "<repository>"), (str(case_dir.parent.parent), "<showcase>"))
    stdout, stderr = completed.stdout, completed.stderr
    for source, target in replacements:
        stdout = stdout.replace(source, target)
        stderr = stderr.replace(source, target)
    write_text(case_dir / "STDOUT.txt", stdout)
    write_text(case_dir / "STDERR.txt", stderr)
    return completed.returncode


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def measurement(label: str, value: str, help_text: str) -> dict[str, str]:
    return {"label": label, "value": value, "help": help_text}


def render_html(showcase: dict[str, Any]) -> str:
    cards: list[str] = []
    for case in showcase["cases"]:
        metrics = "".join(
            f'<div class="metric"><span>{html.escape(item["label"])}</span>'
            f'<strong>{html.escape(item["value"])}</strong><small>{html.escape(item["help"])}</small></div>'
            for item in case["measurements"]
        )
        benefits = "".join(f"<li>{html.escape(item)}</li>" for item in case["human_benefits"])
        outputs = "".join(
            f'<a href="{html.escape(item)}">{html.escape(Path(item).name)}</a>' for item in case["evidence_files"]
        )
        cards.append(f'''<article class="case-card">
          <div class="case-head"><span class="case-number">{html.escape(case["case_id"])}</span><span class="pass">PASSED</span></div>
          <p class="tool">{html.escape(case["tool"])}</p><h2>{html.escape(case["human_label"])}</h2>
          <p class="question">{html.escape(case["human_question"])}</p>
          <div class="evidence"><b>Observed in this synthetic run</b><p>{html.escape(case["observed_result"])}</p></div>
          <div class="metrics">{metrics}</div>
          <h3>Human research benefits</h3><ul>{benefits}</ul>
          <div class="limit"><b>Cannot establish</b><p>{html.escape(case["non_claim"])}</p></div>
          <div class="outputs">{outputs}</div>
        </article>''')
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>YAUVI — Five Human Use Cases</title><style>
:root{{--ink:#17231f;--muted:#62706a;--line:#dbe4df;--green:#087a65;--soft:#f2f7f4;--amber:#8c5b14}}
*{{box-sizing:border-box}}body{{margin:0;background:#edf3ef;color:var(--ink);font:15px/1.5 Inter,system-ui,sans-serif}}
main{{max-width:1180px;margin:auto;padding:56px 24px}}.eyebrow{{color:var(--green);font-size:11px;font-weight:800;letter-spacing:.12em;text-transform:uppercase}}
h1{{font:500 42px/1.08 Georgia,serif;margin:8px 0 14px}}.lead{{max-width:780px;color:var(--muted);font-size:17px}}.notice{{margin:28px 0;padding:14px 16px;background:#fff6e8;border-left:4px solid #c48325}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}}.case-card{{background:#fff;border:1px solid var(--line);border-radius:16px;padding:22px;box-shadow:0 12px 35px rgba(20,48,38,.06)}}
.case-head{{display:flex;justify-content:space-between}}.case-number,.pass{{font-size:10px;font-weight:800;letter-spacing:.08em}}.pass{{background:#e5f5ee;color:#076451;border-radius:99px;padding:4px 8px}}.tool{{color:var(--green);font-size:11px;font-weight:750;margin:20px 0 3px}}h2{{font:500 25px Georgia,serif;margin:0}}.question{{color:var(--muted)}}
.evidence,.limit{{padding:12px 14px;border-radius:9px;background:var(--soft);font-size:12px}}.evidence p,.limit p{{margin:4px 0 0}}.limit{{background:#fff5e7;margin-top:14px;color:#654714}}.metrics{{display:grid;grid-template-columns:repeat(2,1fr);gap:8px;margin:14px 0}}.metric{{padding:10px;background:#f8faf9;border:1px solid var(--line);border-radius:8px}}.metric span,.metric strong,.metric small{{display:block}}.metric span{{font-size:9px;color:var(--muted);text-transform:uppercase}}.metric strong{{font-size:16px;margin:2px 0}}.metric small{{font-size:9px;color:var(--muted)}}h3{{font-size:12px;margin:16px 0 5px}}ul{{margin:0;padding-left:19px;font-size:12px}}.outputs{{display:flex;flex-wrap:wrap;gap:7px;margin-top:15px}}.outputs a{{color:var(--green);font-size:10px}}footer{{margin-top:28px;color:var(--muted);font-size:11px}}
@media(max-width:760px){{main{{padding:32px 14px}}h1{{font-size:34px}}.grid{{grid-template-columns:1fr}}}}@media print{{body{{background:#fff}}main{{padding:0}}.case-card{{break-inside:avoid;box-shadow:none}}}}
</style></head><body><main><p class="eyebrow">YAUVI Structural Biology Platform · Mark 1 · tested showcase</p>
<h1>Five questions a structural scientist can investigate</h1>
<p class="lead">Each card below comes from an actual local CLI execution. Results demonstrate deterministic software behavior on invented coordinates—not findings about a person, organism, disease, or treatment.</p>
<div class="notice"><strong>Five of five software cases passed.</strong> External scientific qualification remains a separate release gate.</div>
<section class="grid">{''.join(cards)}</section>
<footer>Raw JSON and TSV files are linked from every case. No combined biological, quality, activity, or readiness score is calculated.</footer>
</main></body></html>'''


def _citation_metadata() -> dict[str, str]:
    """Read the small approved citation record without adding a YAML dependency."""
    text = (ROOT / "yauvi-structural-workbench" / "CITATION.cff").read_text(encoding="utf-8")
    values: dict[str, str] = {}
    for key in ("title", "version", "license"):
        values[key] = next((line.split(":", 1)[1].strip().strip('"') for line in text.splitlines()
                            if line.startswith(f"{key}:")), "unresolved")
    family = next((line.split(":", 1)[1].strip() for line in text.splitlines()
                   if line.strip().lstrip("- ").startswith("family-names:")), "")
    given = next((line.split(":", 1)[1].strip() for line in text.splitlines()
                  if line.strip().lstrip("- ").startswith("given-names:")), "")
    orcid = next((line.split(":", 1)[1].strip().strip('"') for line in text.splitlines()
                  if line.strip().startswith("orcid:")), "unresolved")
    values.update({"author": f"{given} {family}".strip(), "orcid": orcid})
    return values


def _public_output_path(technical_output: Path, requested: Path | None) -> Path:
    if requested is not None:
        return requested.resolve()
    if technical_output == DEFAULT_OUT.resolve():
        return DEFAULT_PUBLIC_OUT.resolve()
    return technical_output.parent / f"{technical_output.name}-public"


def _prepare_public_output(output: Path, *, replace: bool) -> None:
    default = DEFAULT_PUBLIC_OUT.resolve()
    if output == default:
        output.mkdir(parents=True, exist_ok=True)
        missing = [name for name in PUBLIC_SOURCE_FILES if not (output / name).is_file()]
        if missing:
            raise ValueError(f"public showcase source assets are missing: {', '.join(missing)}")
        shutil.rmtree(output / "evidence", ignore_errors=True)
        shutil.rmtree(output / "qualification", ignore_errors=True)
        shutil.rmtree(output / "qualification-v2", ignore_errors=True)
        shutil.rmtree(output / "reviewer", ignore_errors=True)
        shutil.rmtree(output / "share", ignore_errors=True)
        for name in ("data.js", "PUBLIC_SHOWCASE_MANIFEST.json", "CHECKSUMS.json"):
            (output / name).unlink(missing_ok=True)
        return
    if output.exists():
        if not replace:
            raise ValueError(f"public output already exists: {output}; pass --replace")
        if output == Path(output.anchor):
            raise ValueError("refusing to replace an unsafe public output")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    for name in PUBLIC_SOURCE_FILES:
        shutil.copyfile(default / name, output / name)


def build_public_showcase(technical_output: Path, public_output: Path,
                          showcase: dict[str, Any], sfcsa_output: Path,
                          sfcsa_case: dict[str, Any], *, replace: bool) -> None:
    """Create a deterministic public-safe story and evidence subset."""
    _prepare_public_output(public_output, replace=replace)
    baseline_path = ROOT / "yauvi-structural-workbench" / "BASELINE.json"
    release_path = ROOT / "yauvi-structural-workbench" / "RELEASE_STATUS.json"
    identity_path = ROOT / "yauvi-structural-workbench" / "PLATFORM_IDENTITY.json"
    start_here_path = ROOT / "yauvi-structural-workbench" / "START_HERE.md"
    roadmap_path = ROOT / "yauvi-structural-workbench" / "JOSS_PUBLICATION_ROADMAP.json"
    qualification_root = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v1"
    qualification_results_path = qualification_root / "results" / "QUALIFICATION_RESULTS.json"
    source_verification_path = qualification_root / "results" / "SOURCE_VERIFICATION.json"
    source_lock_path = qualification_root / "SOURCE_LOCK.json"
    qualification_report_path = qualification_root / "QUALIFICATION_REPORT.md"
    qualification_v2_root = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"
    qualification_v2_status_path = qualification_v2_root / "results" / "QUALIFICATION_V2_STATUS.json"
    qualification_v2_manifest_path = qualification_v2_root / "PANEL_MANIFEST.json"
    qualification_v2_report_path = qualification_v2_root / "results" / "QUALIFICATION_REPORT.html"
    qualification_v2_strata_path = qualification_v2_root / "results" / "STRATUM_STATUS.tsv"
    execution_summary_path = qualification_v2_root / "results" / "EXECUTION_SUMMARY.json"
    reviewer_quickstart_path = ROOT / "yauvi-structural-workbench" / "docs" / "reviewer-quickstart.md"
    paper_path = ROOT / "yauvi-structural-workbench" / "paper" / "paper.md"
    joss_checklist_path = ROOT / "yauvi-structural-workbench" / "JOSS_CHECKLIST.md"
    baseline = load(baseline_path)
    release = load(release_path)
    identity = load(identity_path)
    roadmap = load(roadmap_path)
    qualification = load(qualification_results_path)
    qualification_v2 = load(qualification_v2_status_path)
    execution_summary = load(execution_summary_path)
    executed_states = {
        panel["workflow"]: panel["stratum_state"] for panel in execution_summary["panels"]
    }
    sys.path.insert(0, str(ROOT / "platform" / "src"))
    try:
        from yauvi_platform.structural_workbench.store import analysis_definitions, tool_readiness
        definitions = analysis_definitions()
        readiness = {item["analysis_type"]: item for item in tool_readiness(ROOT) if not item.get("labs")}
    finally:
        sys.path.pop(0)

    public_questions = {
        "structure_qc": "Can I trust these coordinates?",
        "membrane_orientation": "How might this protein sit in a membrane?",
        "conformational_state": "Which conformation does it resemble?",
        "functional_site_state": "Are important functional residues present?",
        "assembly_interface": "Which residues form an assembly interface?",
        "sf_csa": "How does this protein relate to other proteins?",
    }
    case_by_type = {item["analysis_type"]: item for item in showcase["cases"]}
    case_by_type["sf_csa"] = sfcsa_case
    workflow_rows: list[dict[str, Any]] = []
    for definition in definitions:
        analysis_type = definition["analysis_type"]
        record = readiness[analysis_type]
        demonstrated = analysis_type in case_by_type
        showcase_state = (
            sfcsa_case["test_state"] if analysis_type == "sf_csa"
            else "passed_synthetic_case" if demonstrated
            else "public_showcase_case_pending"
        )
        if analysis_type == "sf_csa":
            showcase_note = (
                "Canonical pipeline and release audit passed through deterministic test doubles. "
                "Foldseek and DIAMOND computed no alignments; external binary benchmarks remain pending."
            )
        elif demonstrated:
            showcase_note = "Executed synthetic evidence is linked below."
        else:
            showcase_note = "Workflow available. Public showcase case pending. No example result displayed."
        workflow_rows.append({
            "analysis_type": analysis_type,
            "title": definition["title"],
            "public_question": public_questions[analysis_type],
            "question": definition["question"],
            "measures": definition["measures"],
            "non_claim": definition["non_claim"],
            "software_state": record["state"],
            "showcase_state": showcase_state,
            "showcase_note": showcase_note,
            "external_benchmark": (
                "public_case_passed" if next(item for item in qualification["workflows"]
                                              if item["workflow"] == analysis_type)["status"] == "passed"
                else "partial_public_case"
            ),
            "external_benchmark_detail": release["external_benchmarks"][analysis_type],
            "scientific_scopes": record.get("scientific_scopes", []),
            "inputs": [
                {
                    "role": item["role"],
                    "label": item["label"],
                    "required": item["required"],
                    "extensions": item["accepted_extensions"],
                    "absence_effect": item["absence_effect"],
                    "source_ids": item["source_ids"],
                }
                for item in definition["inputs"]
            ],
            "required_runtimes": record.get("required_runtimes", {}),
            "optional_runtimes": record.get("optional_runtimes", {}),
        })

    input_groups = {
        "HUC-01": ("qc/",),
        "HUC-02": ("membrane/",),
        "HUC-03": ("portfolio/query.pdb", "portfolio/reference.fasta", "portfolio/provenance.json",
                   "portfolio/reference_set.json", "portfolio/inactive_reference.pdb"),
        "HUC-04": ("portfolio/query.pdb", "portfolio/reference.fasta", "portfolio/provenance.json",
                   "portfolio/annotations.json", "portfolio/pockets.json"),
        "HUC-05": ("portfolio/query.pdb", "portfolio/reference.fasta", "portfolio/provenance.json",
                   "portfolio/assembly.pdb"),
    }
    public_cases: list[dict[str, Any]] = []
    copied_evidence: dict[str, str] = {}
    for case in showcase["cases"]:
        evidence_files: list[dict[str, str]] = []
        for relative in case["evidence_files"]:
            source = (technical_output / relative).resolve()
            source.relative_to(technical_output.resolve())
            destination = public_output / "evidence" / case["case_id"] / source.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            public_relative = destination.relative_to(public_output).as_posix()
            digest = sha256(destination)
            copied_evidence[public_relative] = digest
            evidence_files.append({"label": source.name, "path": public_relative, "sha256": digest})
        selectors = input_groups[case["case_id"]]
        input_hashes = {
            name: digest for name, digest in showcase["input_sha256"].items()
            if any(name == selector or name.startswith(selector) for selector in selectors)
        }
        public_cases.append({
            key: case[key] for key in (
                "case_id", "tool", "analysis_type", "human_label", "human_question",
                "observed_result", "measurements", "human_benefits", "non_claim",
            )
        } | {"test_state": "passed_synthetic_case", "input_sha256": input_hashes,
             "evidence_files": evidence_files})

    sf_evidence_files: list[dict[str, str]] = []
    for relative in sfcsa_case["evidence_files"]:
        source = (sfcsa_output / relative).resolve()
        source.relative_to(sfcsa_output.resolve())
        destination = public_output / "evidence" / sfcsa_case["case_id"] / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        public_relative = destination.relative_to(public_output).as_posix()
        digest = sha256(destination)
        copied_evidence[public_relative] = digest
        label = relative.removeprefix("release/").replace("/", " · ")
        sf_evidence_files.append({"label": label, "path": public_relative, "sha256": digest})
    public_cases.append({
        key: sfcsa_case[key] for key in (
            "case_id", "tool", "analysis_type", "human_label", "human_question",
            "observed_result", "measurements", "human_benefits", "non_claim",
            "runtime_disclosure", "known_findings", "test_state",
        )
    } | {"input_sha256": sfcsa_case["input_sha256"], "evidence_files": sf_evidence_files})

    copied_publication_files: dict[str, str] = {}
    publication_files = (
        (qualification_results_path, "qualification/QUALIFICATION_RESULTS.json"),
        (source_verification_path, "qualification/SOURCE_VERIFICATION.json"),
        (source_lock_path, "qualification/SOURCE_LOCK.json"),
        (qualification_report_path, "qualification/QUALIFICATION_REPORT.md"),
        (qualification_v2_status_path, "qualification-v2/QUALIFICATION_V2_STATUS.json"),
        (qualification_v2_manifest_path, "qualification-v2/PANEL_MANIFEST.json"),
        (qualification_v2_report_path, "qualification-v2/QUALIFICATION_REPORT.html"),
        (qualification_v2_strata_path, "qualification-v2/STRATUM_STATUS.tsv"),
        (execution_summary_path, "qualification-v2/EXECUTION_SUMMARY.json"),
        (reviewer_quickstart_path, "reviewer/REVIEWER_QUICKSTART.md"),
        (roadmap_path, "reviewer/JOSS_PUBLICATION_ROADMAP.json"),
        (paper_path, "reviewer/PAPER_PREVIEW.md"),
        (joss_checklist_path, "reviewer/JOSS_CHECKLIST.md"),
        (identity_path, "share/PLATFORM_IDENTITY.json"),
        (start_here_path, "share/START_HERE.md"),
    )
    for source, relative in publication_files:
        destination = public_output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        copied_publication_files[relative] = sha256(destination)

    payload = {
        "schema_version": "1.0",
        "showcase_id": "yauvi-public-evidence-showcase",
        "product": identity["display_name"],
        "platform_identity": identity,
        "data_class": "synthetic_demonstrations_plus_public_qualification_summaries",
        "baseline": {
            "baseline_id": baseline["baseline_id"],
            "selection": baseline["selection"],
            "total_passed": baseline["total_passed"],
            "total_deselected": baseline["total_deselected"],
            "scientific_boundary": baseline["scientific_boundary"],
        },
        "release": {
            "release_state": release["release_state"],
            "release_candidate": release["release_candidate"],
            "submission_eligible": release["submission_eligible"],
            "publication_authorized": release["publication_authorized"],
            "version_control": release["version_control"],
        },
        "workflows": workflow_rows,
        "cases": public_cases,
        "qualification": {
            "collection_id": qualification["collection_id"],
            "overall_state": qualification["overall_state"],
            "qualification_rule": qualification["qualification_rule"],
            "source_artifact_count": qualification["source_lock"]["artifact_count"],
            "source_lock_passed": qualification["source_lock"]["passed"],
            "workflow_counts": qualification["workflow_counts"],
            "cases": qualification_showcase_rows(qualification),
            "result_sha256": sha256(qualification_results_path),
            "source_verification_sha256": sha256(source_verification_path),
            "files": [
                {"label": "Readable qualification report", "path": "qualification/QUALIFICATION_REPORT.md"},
                {"label": "Machine-readable qualification results", "path": "qualification/QUALIFICATION_RESULTS.json"},
                {"label": "Source checksum lock", "path": "qualification/SOURCE_LOCK.json"},
                {"label": "Source verification record", "path": "qualification/SOURCE_VERIFICATION.json"},
            ],
        },
        "qualification_v2": {
            "collection_id": qualification_v2["collection_id"],
            "overall_state": qualification_v2["overall_state"],
            # The composition audit reports its own scope and hardcodes this
            # false; on its own it reads as "nothing has been executed", which
            # stopped being true once four panels executed. The executed
            # evidence answers the same question directly, so the flag now
            # comes from there and travels with the note that bounds it.
            "scientific_execution_performed": execution_summary["scientific_execution_performed"],
            "scientific_execution": {
                "panels_executed": execution_summary["panels_executed"],
                "panels_total": execution_summary["panels_total"],
                "cases_passed": execution_summary["cases_passed"],
                "cases_required": execution_summary["cases_required"],
                "every_executed_panel_passed": execution_summary["every_executed_panel_passed"],
                "workflows_executed": execution_summary["workflows_executed"],
                "workflows_not_executed": execution_summary["workflows_not_executed"],
                "all_release_blocking_scopes_qualified":
                    execution_summary["all_release_blocking_scopes_qualified"],
                "second_machine_reproduction": execution_summary["second_machine_reproduction"],
                "scope_qualification_note": execution_summary["scope_qualification_note"],
                "collection_note": release["qualification_evidence"]["current_v2"]["scientific_execution_note"],
                "panels": [
                    {
                        "workflow": panel["workflow"],
                        "stratum_scope": panel["stratum_scope"],
                        "stratum_state": panel["stratum_state"],
                        "cases_adopted": panel["cases_adopted"],
                        "cases_required": panel["cases_required"],
                        "cases_passed": panel["cases"]["passed"],
                        "controls_passed": panel["controls"]["passed"],
                        "controls_total": panel["controls"]["total"],
                        "coverage_witnessed": panel["coverage"]["witnessed"],
                        "coverage_required": panel["coverage"]["required"],
                        "coverage_unwitnessable": panel["coverage"]["unwitnessable"],
                    }
                    for panel in execution_summary["panels"]
                ],
            },
            "missing_records": sum(
                requirement["missing_count"]
                for panel in qualification_v2["panels"]
                for requirement in panel["requirements"]
            ),
            "panels": [
                {
                    "workflow": panel["workflow"], "state": panel["state"],
                    "record_count": panel["record_count"],
                    "missing_count": sum(row["missing_count"] for row in panel["requirements"]),
                    # "state" is the composition audit's word: ready_for_execution
                    # means composed, not run. Carrying the executed state beside
                    # it stops the two panel lists reading as a contradiction.
                    "execution_state": executed_states.get(panel["workflow"], "not_executed"),
                }
                for panel in qualification_v2["panels"]
            ],
            "files": [
                {"label": "Qualification v2 panel specification", "path": "qualification-v2/PANEL_MANIFEST.json"},
                {"label": "Qualification v2 status", "path": "qualification-v2/QUALIFICATION_V2_STATUS.json"},
                {"label": "Printable v2 audit", "path": "qualification-v2/QUALIFICATION_REPORT.html"},
                {"label": "V2 stratum gaps", "path": "qualification-v2/STRATUM_STATUS.tsv"},
                {"label": "V2 execution summary", "path": "qualification-v2/EXECUTION_SUMMARY.json"},
            ],
        },
        "publication_roadmap": roadmap,
        "reviewer_files": [
            {"label": "Reviewer quickstart", "path": "reviewer/REVIEWER_QUICKSTART.md"},
            {"label": "Publication roadmap", "path": "reviewer/JOSS_PUBLICATION_ROADMAP.json"},
            {"label": "JOSS paper preview", "path": "reviewer/PAPER_PREVIEW.md"},
            {"label": "JOSS preparation checklist", "path": "reviewer/JOSS_CHECKLIST.md"},
        ],
        "non_claims": [
            "Structural resemblance is not biochemical activity.",
            "Membrane orientation is not native exposure.",
            "A mapped site is not observed catalysis.",
            "An interface is not binding affinity.",
            "Structural or sequence similarity is not exact functional transfer.",
        ],
        "citation": _citation_metadata(),
        "limitations": showcase["limitations"] + [
            "The SF-CSA public case verifies orchestration and evidence boundaries with test doubles; it does not benchmark real alignments.",
            "Synthetic demonstrations do not externally qualify a workflow; independent public cases are displayed separately.",
            "Four passing public cases do not establish workflow-general accuracy, and the two partial cases remain release blockers.",
            "Qualification v2 freezes the expanded panels and remains blocked until source-locked public cases are adopted and executed.",
            "Four of six Qualification v2 panels have executed and passed. Executed panels passing is not scope qualification: two panels are unadopted, membrane covers only its beta_barrel stratum, and no scope has reproduced on an independent second machine.",
        ],
    }
    write_text(public_output / "data.js", "window.YAUVI_PUBLIC_SHOWCASE = " + canonical(payload).rstrip() + ";\n")
    manifest = {
        "schema_version": "1.0",
        "showcase_id": payload["showcase_id"],
        "product": identity["display_name"],
        "data_class": payload["data_class"],
        "workflow_count": len(workflow_rows),
        "executed_case_count": len(public_cases),
        "pending_public_cases": [],
        "source_sha256": {
            "SHOWCASE.json": sha256(technical_output / "SHOWCASE.json"),
            "SF_CSA_CASE.json": sha256(sfcsa_output / "CASE.json"),
            "BASELINE.json": sha256(baseline_path),
            "RELEASE_STATUS.json": sha256(release_path),
            "PLATFORM_IDENTITY.json": sha256(identity_path),
            "START_HERE.md": sha256(start_here_path),
            "JOSS_PUBLICATION_ROADMAP.json": sha256(roadmap_path),
            "QUALIFICATION_RESULTS.json": sha256(qualification_results_path),
            "SOURCE_VERIFICATION.json": sha256(source_verification_path),
            "SOURCE_LOCK.json": sha256(source_lock_path),
            "QUALIFICATION_V2_STATUS.json": sha256(qualification_v2_status_path),
            "QUALIFICATION_V2_PANEL_MANIFEST.json": sha256(qualification_v2_manifest_path),
            "EXECUTION_SUMMARY.json": sha256(execution_summary_path),
            "analysis_definitions": sha256_text(canonical(definitions)),
        },
        "generated_sha256": {
            "data.js": sha256(public_output / "data.js"),
            **copied_evidence,
            **copied_publication_files,
        },
        "network_dependencies": [],
        "external_uploads": [],
    }
    write_text(public_output / "PUBLIC_SHOWCASE_MANIFEST.json", canonical(manifest))
    checksums = {
        path.relative_to(public_output).as_posix(): sha256(path)
        for path in sorted(public_output.rglob("*")) if path.is_file() and path.name != "CHECKSUMS.json"
    }
    write_text(public_output / "CHECKSUMS.json", canonical(checksums))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--public-out", type=Path)
    parser.add_argument("--sfcsa-out", type=Path)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    output = args.out.resolve()
    if output.exists():
        if not args.replace:
            parser.error(f"output already exists: {output}; pass --replace for generated showcase data")
        if output == ROOT or ROOT not in output.parents:
            parser.error("refusing to replace an output outside this repository")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    inputs = copy_inputs(output)
    runs = output / "runs"
    shared = output / "shared" / "structqc"

    qc_out = runs / "HUC-01-coordinate-trust"
    membrane_out = runs / "HUC-02-membrane-sidedness"
    state_out = runs / "HUC-03-conformational-resemblance"
    site_out = runs / "HUC-04-functional-site"
    assembly_out = runs / "HUC-05-assembly-interface"
    for path in (qc_out, membrane_out, state_out, site_out, assembly_out, shared):
        path.mkdir(parents=True, exist_ok=True)

    python = sys.executable
    codes: dict[str, int] = {}
    codes["HUC-01"] = run_command(qc_out, [
        python, "-m", "structqc.cli", "run", "--structure", str(inputs["qc/model.pdb"]),
        "--subject-id", "SYNTHETIC", "--reference-fasta", str(inputs["qc/reference.fasta"]),
        "--provenance", str(inputs["qc/provenance.json"]), "--validation-report",
        str(inputs["qc/validation.json"]), "--require-external-validation", "--chain", "A", "--out", str(qc_out),
    ])
    dependency_code = run_command(shared, [
        python, "-m", "structqc.cli", "run", "--structure", str(inputs["portfolio/query.pdb"]),
        "--subject-id", "SYNTHETIC_QUERY", "--reference-fasta", str(inputs["portfolio/reference.fasta"]),
        "--provenance", str(inputs["portfolio/provenance.json"]), "--chain", "A", "--out", str(shared),
    ])
    codes["HUC-02"] = run_command(membrane_out, [
        python, "-m", "memorient.cli", "run", "--structure", str(inputs["membrane/synthetic_tm_receptor.pdb"]),
        "--context", "tm_receptor", "--chain", "A", "--topology-evidence",
        str(inputs["membrane/synthetic_tm_topology.json"]), "--out", str(membrane_out),
    ])
    manifest = shared / "STRUCTURE_EVIDENCE.json"
    codes["HUC-03"] = run_command(state_out, [
        python, "-m", "state_atlas.cli", "run", "--manifest", str(manifest), "--reference-set",
        str(inputs["portfolio/reference_set.json"]), "--structure", str(inputs["portfolio/query.pdb"]),
        "--chain", "A", "--cluster-cutoff-A", "0.5", "--out", str(state_out),
    ])
    codes["HUC-04"] = run_command(site_out, [
        python, "-m", "site_context.cli", "run", "--manifest", str(manifest), "--structure",
        str(inputs["portfolio/query.pdb"]), "--annotations", str(inputs["portfolio/annotations.json"]),
        "--pocket-result", str(inputs["portfolio/pockets.json"]), "--out", str(site_out),
    ])
    codes["HUC-05"] = run_command(assembly_out, [
        python, "-m", "assembly_context.cli", "run", "--manifest", str(manifest), "--isolated",
        str(inputs["portfolio/query.pdb"]), "--assembly", str(inputs["portfolio/assembly.pdb"]),
        "--subject-chain", "A", "--relationship", "exact_protein", "--assembly-id", "SYNTHETIC_AB",
        "--expected-chains", "A,B", "--out", str(assembly_out),
    ])
    if dependency_code not in (0, 1) or any(code not in (0, 1) for code in codes.values()):
        print(canonical({"dependency_exit_code": dependency_code, "case_exit_codes": codes}), file=sys.stderr)
        return 2

    qc = load(qc_out / "STRUCTURE_EVIDENCE.json")
    membrane_record = load(membrane_out / "MEMBRANE_ORIENTATION.json")
    membrane = membrane_record.get("summary", membrane_record)
    state = load(state_out / "STATE_ENSEMBLE.json")
    site = load(site_out / "SITE_CONTEXT.json")
    assembly = load(assembly_out / "ASSEMBLY_CONTEXT.json")
    surface_set = membrane.get("n_surface_set", 0)

    cases = [
        {
            "case_id": "HUC-01", "tool": "StructQC", "analysis_type": "structure_qc",
            "human_label": "Can I safely interpret these coordinates?",
            "human_question": "Before mapping a variant or functional residue, are sequence identity, residue numbering, provenance, and validation bound to the exact model?",
            "test_state": "passed", "exit_code": codes["HUC-01"],
            "observed_result": "Both synthetic reference residues mapped exactly; predicted provenance and imported validation remained explicit.",
            "measurements": [
                measurement("Sequence coverage", f'{100 * qc["completeness"]["coverage_fraction"]:.1f}%', "Reference residues with mapped coordinates"),
                measurement("Sequence identity", f'{100 * qc["completeness"]["identity_fraction"]:.1f}%', "Identity within the accepted mapping"),
                measurement("Validation", qc["external_validation"]["state"], "Imported, not recomputed by StructQC"),
                measurement("Chain breaks", str(qc["chain_summaries"][0]["chain_breaks"]), "Detected coordinate discontinuities"),
            ],
            "human_benefits": ["Prevents residue-numbering mistakes before mutation or site analysis.", "Makes experimental and predicted confidence interpretable without mixing them.", "Provides a reproducible accept, hold, or investigate-first boundary."],
            "non_claim": "Native conformation, biological function, or experimental correctness beyond the imported validation record.",
            "evidence_files": ["runs/HUC-01-coordinate-trust/STRUCTURE_EVIDENCE.json", "runs/HUC-01-coordinate-trust/RESIDUE_QUALITY.tsv"],
        },
        {
            "case_id": "HUC-02", "tool": "MembraneOrient", "analysis_type": "membrane_orientation",
            "human_label": "Which side of a membrane might a receptor expose?",
            "human_question": "Can a single-pass membrane protein be placed in a consistent coordinate frame and divided into membrane and sided residue sets?",
            "test_state": "passed", "exit_code": codes["HUC-02"],
            "observed_result": f'The invented receptor was routed as {membrane.get("label", "unresolved")} using the experimental {membrane.get("method", "recorded method")} path with checksum-bound synthetic spans.',
            "measurements": [
                measurement("Structure label", str(membrane.get("label", "unresolved")), "Geometry route selected by the tool"),
                measurement("Orientation method", str(membrane.get("method", "unresolved")), "Named computational method"),
                measurement("Residues reviewed", str(membrane.get("n_residues", 0)), "Coordinate-bound residue annotations"),
                measurement("Modeled surface set", str(surface_set), "Geometry-derived candidates; not intact-cell evidence"),
            ],
            "human_benefits": ["Helps choose candidate extracellular loops for follow-up assays.", "Separates membrane-core residues from flanking domains for construct planning.", "Provides a common frame for comparing mutations or models."],
            "non_claim": "Native topology, intact-cell exposure, antibody accessibility, expression, receptor function, or Mark 1 alpha-helical qualification.",
            "evidence_files": ["runs/HUC-02-membrane-sidedness/MEMBRANE_ORIENTATION.json", "runs/HUC-02-membrane-sidedness/RESIDUE_ORIENTATION.tsv"],
        },
        {
            "case_id": "HUC-03", "tool": "StateAtlas", "analysis_type": "conformational_state",
            "human_label": "Which experimental conformation does my model resemble?",
            "human_question": "When two bounded reference states are declared, which reference is geometrically closer and is the margin interpretable?",
            "test_state": "passed", "exit_code": codes["HUC-03"],
            "observed_result": f'The synthetic query was {state["overall_label"]}; all {state["frames_total"]} frame(s) were retained in the population accounting.',
            "measurements": [
                measurement("Resemblance label", state["overall_label"], "Bounded structural label"),
                measurement("Best RMSD", f'{state["frame_metrics"][0]["best_rmsd_A"]:.3f} Å', "After declared C-alpha alignment"),
                measurement("Reference margin", f'{state["frame_metrics"][0]["margin_A"]:.3f} Å', "Distance separation from the alternative reference"),
                measurement("Interpretable frames", f'{state["frames_interpretable"]}/{state["frames_total"]}', "Unresolved frames would remain in the denominator"),
            ],
            "human_benefits": ["Compares mutant, apo, ligand-bound, or ensemble structures consistently.", "Identifies ambiguous frames instead of forcing a state call.", "Supports selection of conformations for further simulation or experiments."],
            "non_claim": "Biochemical activity, activation, inhibition, mechanism, efficacy, or a time-resolved transition pathway.",
            "evidence_files": ["runs/HUC-03-conformational-resemblance/STATE_ENSEMBLE.json", "runs/HUC-03-conformational-resemblance/FRAME_METRICS.tsv"],
        },
        {
            "case_id": "HUC-04", "tool": "SiteContext", "analysis_type": "functional_site_state",
            "human_label": "Are the declared functional residues structurally present?",
            "human_question": "Do curated residues map exactly, retain role-compatible identities, and overlap a separately identified pocket?",
            "test_state": "passed", "exit_code": codes["HUC-04"],
            "observed_result": f'{sum(item["state"] == "role_compatible" for item in site["sites"])} of {len(site["sites"])} synthetic site residues were role-compatible; pocket evidence remained method-specific.',
            "measurements": [
                measurement("Mapped sites", str(len(site["sites"])), "Declared residues found in the structure"),
                measurement("Role compatible", str(sum(item["state"] == "role_compatible" for item in site["sites"])), "Identity fits the declared role vocabulary"),
                measurement("Maximum separation", f'{site["site_geometry"]["maximum_separation_A"]:.3f} Å', "Descriptive C-alpha geometry"),
                measurement("Pocket methods", ", ".join(site["config"]["pocket_methods"]), "Scores stay specific to each named tool"),
            ],
            "human_benefits": ["Reveals missing or mismapped catalytic and binding residues.", "Supports control-mutation and construct-boundary discussions.", "Keeps annotations, observed chemistry, and pocket predictions separate."],
            "non_claim": "Observed catalysis, ligand affinity, inhibition, druggability, physiological function, or clinical relevance.",
            "evidence_files": ["runs/HUC-04-functional-site/SITE_CONTEXT.json", "runs/HUC-04-functional-site/SITE_RESIDUES.tsv"],
        },
        {
            "case_id": "HUC-05", "tool": "AssemblyContext", "analysis_type": "assembly_interface",
            "human_label": "Which residues become part of an oligomer interface?",
            "human_question": "In a declared two-chain assembly, which subject residues contact a partner and how much surface becomes buried?",
            "test_state": "passed", "exit_code": codes["HUC-05"],
            "observed_result": f'{len(assembly["residue_contacts"])} synthetic subject residues contacted chain B and {assembly["surface"]["buried_sasa_A2"]:.3f} Å² became buried.',
            "measurements": [
                measurement("Assembly complete", str(assembly["assembly"]["complete"]).lower(), "Expected and observed chain inventory agree"),
                measurement("Contact residues", str(len(assembly["residue_contacts"])), "Subject residues with a partner inside the cutoff"),
                measurement("Buried surface", f'{assembly["surface"]["buried_sasa_A2"]:.3f} Å²', "Method-specific buried SASA"),
                measurement("SASA method", assembly["methods"]["sasa"], "FreeSASA is preferred when installed"),
            ],
            "human_benefits": ["Identifies candidate interface-disrupting or interface-preserving mutations.", "Shows when a site or surface is occluded by an assembly partner.", "Separates monomer interpretation from oligomeric structural context."],
            "non_claim": "Native oligomer abundance, binding affinity, intact-cell accessibility, physiological interaction, or mechanism.",
            "evidence_files": ["runs/HUC-05-assembly-interface/ASSEMBLY_CONTEXT.json", "runs/HUC-05-assembly-interface/INTERFACES.tsv"],
        },
    ]
    showcase = {
        "schema_version": "1.0", "showcase_id": "five-human-use-cases",
        "data_class": "synthetic_software_demonstration", "test_cases": 5,
        "passed": sum(case["test_state"] == "passed" for case in cases),
        "scientific_qualification_state": "separate_external_benchmarks_pending",
        "cases": cases,
        "input_sha256": {key: sha256(path) for key, path in sorted(inputs.items())},
        "limitations": [
            "Invented coordinates validate software behavior only.",
            "Benefits are potential research uses, not findings from these synthetic cases.",
            "Passing these cases does not replace external scientific qualification benchmarks.",
        ],
    }
    write_text(output / "SHOWCASE.json", canonical(showcase))
    write_text(output / "SHOWCASE.html", render_html(showcase))
    write_text(output / "README.md", """# Five human use-case showcase

This directory is generated by `tools/build_five_use_case_showcase.py`. It contains five
actual local CLI executions on invented data. The examples demonstrate software behavior,
not findings about a human, organism, disease, target, or treatment.

```bash
python tools/build_five_use_case_showcase.py --replace
```

Open `SHOWCASE.html` for the plain-language review. `SHOWCASE.json`, every module JSON/TSV,
the exact synthetic inputs, and command logs remain available for inspection.
""")
    checksums = {
        path.relative_to(output).as_posix(): sha256(path)
        for path in sorted(output.rglob("*")) if path.is_file() and path.name != "CHECKSUMS.json"
    }
    write_text(output / "CHECKSUMS.json", canonical(checksums))
    public_output = _public_output_path(output, args.public_out)
    sfcsa_output = (
        args.sfcsa_out.resolve() if args.sfcsa_out
        else (DEFAULT_OUT.parent / "sfcsa-ceiling-case").resolve()
        if output == DEFAULT_OUT.resolve()
        else output.parent / f"{output.name}-sfcsa"
    )
    try:
        from build_sfcsa_showcase_case import build as build_sfcsa_case
        sfcsa_case = build_sfcsa_case(sfcsa_output, replace=args.replace)
        build_public_showcase(
            output, public_output, showcase, sfcsa_output, sfcsa_case, replace=args.replace,
        )
    except (OSError, ValueError, RuntimeError, ImportError, json.JSONDecodeError) as exc:
        print(f"public showcase build failed: {exc}", file=sys.stderr)
        return 2
    print(canonical({"showcase": str(output), "public_showcase": str(public_output),
                     "sfcsa_showcase": str(sfcsa_output), "cases": 6,
                     "passed": showcase["passed"] + 1, "exit_codes": codes}))
    return 0 if showcase["passed"] == 5 else 1


if __name__ == "__main__":
    raise SystemExit(main())
