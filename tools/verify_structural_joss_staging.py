#!/usr/bin/env python3
"""Check completeness and honesty of the local pre-public JOSS materials."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import yaml

from yauvi_platform.structural_workbench import analysis_definitions


# GitHub only detects community health files at the repository root, in
# .github/, or in docs/. They therefore live at the root, not in the staging
# directory, and are checked there.
REQUIRED_ROOT = (
    "README.md", "LICENSE", "CITATION.cff", "CONTRIBUTING.md", "SUPPORT.md",
    "SECURITY.md", "GOVERNANCE.md", "CODE_OF_CONDUCT.md", "NOTICE.md",
)
REQUIRED = (
    "README.md", "LICENSE", "CITATION.cff",
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
    problems += [f"missing root {name}" for name in REQUIRED_ROOT if not (root / name).is_file()]
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
    manifest = json.loads((staging / "benchmarks/qualification-v2/PANEL_MANIFEST.json").read_text())
    expected_workflows = {scope.split(":", 1)[0] for scope in manifest["release_blocking_scopes"]}
    if state == "submission_eligible":
        cross_path = staging / "benchmarks/qualification-v2/results/CROSS_MACHINE_REPRODUCTION.json"
        if not cross_path.is_file():
            problems.append("submission_eligible requires current cross-machine evidence")
        else:
            cross = json.loads(cross_path.read_text())
            reproduced = set(cross.get("release_blocking_panels_reproduced", []))
            if not cross.get("every_release_blocking_panel_reproduced") or reproduced != expected_workflows:
                problems.append("submission_eligible requires every manifest-defined blocking workflow")
        required_gates = ("independent_second_machine_reproduction_passed", "independent_research_use_recorded",
                          "license_and_third_party_audit_passed", "ai_tool_versions_fully_recovered",
                          "conflict_and_funding_statements_approved")
        for gate in required_gates:
            if gates.get(gate) is not True:
                problems.append(f"submission_eligible requires {gate}")
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
    if v2.get("scientific_execution_performed") is not False:
        problems.append("composition audit cannot claim scientific execution")
    if not v2.get("panel_composition_ready") and v2.get("overall_state") != "blocked_panel_incomplete":
        problems.append("incomplete composition must be reported as blocked")
    roadmap = json.loads((staging / "JOSS_PUBLICATION_ROADMAP.json").read_text(encoding="utf-8"))
    if not roadmap.get("current_phase"):
        problems.append("publication roadmap must identify its current phase")
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
