#!/usr/bin/env python3
"""Check completeness and honesty of the local pre-public JOSS materials."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import yaml

from yauvi_platform.structural_workbench import analysis_definitions


REQUIRED = (
    "README.md", "LICENSE", "NOTICE.md", "CITATION.cff", "CONTRIBUTING.md",
    "SUPPORT.md", "SECURITY.md", "GOVERNANCE.md", "CODE_OF_CONDUCT.md",
    "RELEASE_STATUS.json", "BASELINE.json", "JOSS_CHECKLIST.md",
    "JOSS_PUBLICATION_ROADMAP.json", "PLATFORM_IDENTITY.json", "START_HERE.md",
    "paper/paper.md", "paper/paper.bib",
    "docs/install.md", "docs/quickstart.md", "docs/workflows.md",
    "docs/files-and-sources.md", "docs/methods-and-limitations.md",
    "docs/benchmarks.md", "docs/reproducibility.md", "docs/reviewer-quickstart.md",
    "benchmarks/benchmark-manifest.yaml",
    "benchmarks/qualification-v2/PANEL_MANIFEST.json",
    "benchmarks/qualification-v2/SOURCE_LOCK.json",
    "benchmarks/qualification-v2/results/QUALIFICATION_V2_STATUS.json",
)
ALLOWED_RELEASE_STATES = {
    "pre_public_preparation", "local_release_candidate", "public_history_in_progress",
    "submission_eligible",
}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    staging = root / "yauvi-structural-workbench"
    problems = [f"missing {name}" for name in REQUIRED if not (staging / name).is_file()]
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    status = json.loads((staging / "RELEASE_STATUS.json").read_text(encoding="utf-8"))
    identity = json.loads((staging / "PLATFORM_IDENTITY.json").read_text(encoding="utf-8"))
    if identity.get("display_name") != "YAUVI Structural Biology Platform — Mark 1":
        problems.append("primary Mark 1 display identity is missing")
    if identity.get("release_state") != status.get("release_state"):
        problems.append("platform identity and release state differ")
    if identity.get("publication_authorized") or identity.get("edition") != "Mark 1":
        problems.append("platform identity overstates publication or has an invalid edition")
    state = status.get("release_state")
    if state not in ALLOWED_RELEASE_STATES:
        problems.append(f"unknown release state: {state}")
    gates = status.get("gates", {})
    if state == "submission_eligible" and not gates.get("all_six_mark_1_qualification_v2_scopes_passed"):
        problems.append("submission_eligible requires all six Mark 1 Qualification v2 scopes")
    if state == "submission_eligible" and not gates.get("public_history_requirement_satisfied"):
        problems.append("submission_eligible requires public-history evidence")

    license_text = (staging / "LICENSE").read_text(encoding="utf-8")
    if "TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION" not in license_text or len(license_text.splitlines()) < 180:
        problems.append("LICENSE is not the complete Apache-2.0 text")
    paper = (staging / "paper" / "paper.md").read_text(encoding="utf-8").lower()
    for phrase in (
        "statement of need", "state of the field", "software design",
        "research impact statement", "limitations", "ai usage disclosure",
    ):
        if phrase not in paper:
            problems.append(f"paper lacks {phrase}")
    benchmark = yaml.safe_load((staging / "benchmarks" / "benchmark-manifest.yaml").read_text(encoding="utf-8"))
    benchmark_ids = set(benchmark.get("workflows", {}))
    definitions = analysis_definitions()
    workflow_ids = {item["analysis_type"] for item in definitions}
    if benchmark_ids != workflow_ids:
        problems.append(f"benchmark coverage differs from workflows: {sorted(benchmark_ids ^ workflow_ids)}")
    v2 = json.loads((staging / "benchmarks" / "qualification-v2" / "results" / "QUALIFICATION_V2_STATUS.json").read_text(encoding="utf-8"))
    if v2.get("overall_state") != "blocked_panel_incomplete" or v2.get("scientific_execution_performed") is not False:
        problems.append("Qualification v2 must remain visibly blocked until public cases are adopted and executed")
    roadmap = json.loads((staging / "JOSS_PUBLICATION_ROADMAP.json").read_text(encoding="utf-8"))
    if roadmap.get("current_phase") != "local_hardening":
        problems.append("publication roadmap does not preserve the local-hardening state")
    roadmap_gates = {item.get("gate_id"): item for item in roadmap.get("gates", [])}
    if roadmap_gates.get("publication_approval", {}).get("state") != "blocked":
        problems.append("publication roadmap does not preserve the approval boundary")
    for definition in definitions:
        for field in ("use_when", "measures", "receives", "non_claim", "scientific_readiness"):
            if not definition.get(field):
                problems.append(f"{definition['analysis_type']} lacks {field}")
        for role in definition["inputs"]:
            for field in (
                "description", "why_needed", "absence_effect", "accepted_artifact_types",
                "accepted_extensions", "format_guide", "source_ids", "validator_id", "sensitivity",
            ):
                if field not in role:
                    problems.append(f"{definition['analysis_type']}:{role['role']} lacks {field}")
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print(f"verified pre-public JOSS staging: {len(definitions)} workflows; release_state={state}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
