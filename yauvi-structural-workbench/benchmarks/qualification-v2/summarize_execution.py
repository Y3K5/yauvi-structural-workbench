#!/usr/bin/env python3
"""Collect what the executed panels actually measured into one readable file.

`run_qualification.py` audits panel *composition* and reports
`scientific_execution_performed: false` by construction -- correct for what it
does, and until now the only thing any reviewer-facing surface read. Every
`run_execution.py` result was written to `results/execution-*/` and read by
nothing, so four executed panels and sixty-four passing cases were invisible
outside the directory that held them.

This reads those results back and writes `results/EXECUTION_SUMMARY.json`, the
input the release status file and the public showcase need in order to state
what has been executed rather than assume it has not.

What this file may not do, deliberately:

* It never qualifies a scope. Qualification needs a *complete* panel, every
  case and control passing, reproduced on an independent second machine. Two
  panels are unadopted and no second-machine run is recorded, so
  `all_release_blocking_scopes_qualified` is false here by construction, the
  same way `scope_qualified` is false in every execution result it reads.
* It never fails the run on a non-blocking scope. Collection 2.4 moved
  membrane orientation to research-only, both strata, and a research-only
  stratum that misses its gate is a recorded finding, not a broken release.
  Which scopes block is read from the manifest, never hardcoded here.
* It records no absolute paths. Each case carries an `output_dir` naming the
  machine that ran it; that is provenance belonging to the run, not to a
  summary the showcase commits, and embedding it would leak a home directory
  into a public artifact.
* It *does* record the interpreter the counts came from. These totals are one
  machine's result, and the membrane stratum in particular reports different
  counts on different runners, so publishing a case total without naming the
  recorder would overstate what has been established -- and the recorder is an
  x86_64 build under Rosetta on arm64 hardware, which matches neither CI
  platform.

Usage:  python summarize_execution.py
Exit:   0 every executed *release-blocking* panel passed, 1 one did not or its
        evidence is unreadable. A failing non-blocking panel is reported in the
        summary and on stdout, and does not change the exit code.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from run_execution import drift_deltas

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MANIFEST = HERE / "PANEL_MANIFEST.json"


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")


def release_blocking(workflow: str, strata: list[str], non_blocking: frozenset[str]) -> bool:
    """Whether this panel's executed evidence gates the release.

    Read from the manifest's `non_blocking_scopes`, which collection 2.4 revised
    to carry both membrane strata. A panel is non-blocking only when *every*
    stratum it actually executed is named there: `conformational_state` has one
    blocking stratum and one non-blocking one, so executing a non-blocking
    stratum of an otherwise blocking workflow must not exempt the workflow.

    A panel that records no strata is treated as blocking. Failing closed is the
    right direction for a question about whether a result may be waived.
    """
    if not strata:
        return True
    return not all(f"{workflow}:{stratum}" in non_blocking for stratum in strata)


def panel_summary(
    status: Mapping[str, Any],
    workflow: str,
    required: int,
    non_blocking: frozenset[str],
) -> dict[str, Any]:
    """Reduce one execution result to the facts a surface can state."""
    cases = status.get("cases", [])
    controls = status.get("controls", [])
    coverage = status.get("coverage", {})
    drift = {
        name: {
            "cases": len(deltas),
            "max": deltas[-1],
            "median": deltas[len(deltas) // 2],
        }
        for name, deltas in drift_deltas(cases + controls).items()
    }
    return {
        "panel_id": status.get("collection_id"),
        "workflow": workflow,
        "stratum_scope": status.get("stratum"),
        "stratum_state": status.get("stratum_state"),
        "strata_executed": sorted(status.get("strata_executed", [])),
        # Whether a failure here blocks the release, from the manifest. The
        # summary states this per panel so a reader never has to infer it from
        # the exit code, which reports only the blocking set.
        "release_blocking": release_blocking(
            workflow, sorted(status.get("strata_executed", [])), non_blocking
        ),
        "cases": status.get("case_counts", {}),
        "controls": status.get("control_counts", {}),
        # Adopted cases against the panel's full requirement. A panel can pass
        # everything it has adopted while still covering only half its scope --
        # membrane is exactly that, with beta_barrel adopted and alpha_helical
        # deliberately deferred -- so both numbers have to travel together.
        "cases_adopted": len(cases),
        "cases_required": required,
        "coverage": {
            "required": len(coverage.get("required", [])),
            "witnessed": len(coverage.get("required", [])) - len(coverage.get("unmet", [])),
            "unmet": sorted(coverage.get("unmet", [])),
            "unwitnessable": sorted(coverage.get("unwitnessable", [])),
        },
        # Informational drift against the recorded run. Never a gate; see the
        # membrane note in run_execution.py. On the machine that recorded the
        # expectations these are zero, so a nonzero value here is cross-machine
        # drift measured against a true reference.
        "drift": drift,
    }


def main() -> int:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    panels = {
        panel["panel_id"]: (
            panel["workflow"],
            sum(requirement["count"] for requirement in panel["requirements"]),
        )
        for panel in manifest["panels"]
    }

    non_blocking = frozenset(manifest["non_blocking_scopes"])

    summaries: list[dict[str, Any]] = []
    for directory in sorted(RESULTS.glob("execution-*")):
        path = directory / "EXECUTION_STATUS.json"
        if not path.is_file():
            continue
        status = json.loads(path.read_text(encoding="utf-8"))
        panel_id = status.get("collection_id")
        if panel_id not in panels:
            print(f"unknown panel in {path.name}: {panel_id}")
            return 1
        workflow, required = panels[panel_id]
        summaries.append(panel_summary(status, workflow, required, non_blocking))

    # One entry per distinct interpreter across the results. Normally one; more
    # than one means the results were not all produced by the same run and the
    # totals should not be read as a single measurement.
    runtimes = sorted({
        json.dumps(json.loads((directory / "EXECUTION_STATUS.json").read_text(encoding="utf-8"))
                   .get("runtime", {}), sort_keys=True)
        for directory in sorted(RESULTS.glob("execution-*"))
        if (directory / "EXECUTION_STATUS.json").is_file()
    })

    executed = {summary["workflow"] for summary in summaries}
    blocking = {scope.split(":", 1)[0] for scope in manifest["release_blocking_scopes"]}
    cases_passed = sum(summary["cases"].get("passed", 0) for summary in summaries)
    cases_required = sum(required for _, required in panels.values())
    all_passed = all(summary["stratum_state"] == "passed" for summary in summaries)
    # Partitioned, because these two answer different questions. A blocking
    # panel that misses its gate stops the release; a non-blocking one records
    # a research finding. Collapsing them is what pinned every CI run red on
    # membrane orientation after collection 2.4 made that scope research-only.
    blocking_failed = sorted(
        summary["workflow"]
        for summary in summaries
        if summary["release_blocking"] and summary["stratum_state"] != "passed"
    )
    non_blocking_failed = sorted(
        summary["workflow"]
        for summary in summaries
        if not summary["release_blocking"] and summary["stratum_state"] != "passed"
    )

    result = {
        "schema_version": "1.0",
        "collection_id": manifest["collection_id"],
        "scientific_execution_performed": bool(summaries),
        "panels_executed": len(summaries),
        "panels_total": len(panels),
        "workflows_executed": sorted(executed),
        "workflows_not_executed": sorted(blocking - executed),
        "cases_passed": cases_passed,
        "cases_required": cases_required,
        "every_executed_panel_passed": all_passed,
        "every_executed_release_blocking_panel_passed": not blocking_failed,
        "release_blocking_panels_failed": blocking_failed,
        "non_blocking_panels_failed": non_blocking_failed,
        # The two gates this file must never be able to close.
        "all_release_blocking_scopes_qualified": False,
        "second_machine_reproduction": "not_recorded",
        "recorded_on": [json.loads(runtime) for runtime in runtimes],
        "counts_are_single_machine": len(runtimes) == 1,
        "scope_qualification_note": (
            "Executed panels passing is not scope qualification. A Mark 1 scope is qualified only "
            "when its panel composes in full, every case and control passes, and the result "
            "reproduces independently on a second machine. This summary reports execution only. "
            f"{len(panels) - len(summaries)} of {len(panels)} panels are unadopted and no "
            "second-machine reproduction is recorded, so no scope is qualified."
        ),
        "panels": summaries,
    }

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "EXECUTION_SUMMARY.json").write_bytes(canonical(result))

    print(f"{len(summaries)}/{len(panels)} panels executed, "
          f"{cases_passed}/{cases_required} cases passed"
          + (f" on {result['recorded_on'][0].get('platform')}/"
             f"{result['recorded_on'][0].get('machine')}"
             f" py{result['recorded_on'][0].get('python')}" if len(runtimes) == 1 else ""))
    for summary in summaries:
        counts = summary["cases"]
        print(f"  {summary['workflow']}: {counts.get('passed', 0)}/{counts.get('total', 0)} cases"
              f" [{summary['stratum_scope']}] {summary['stratum_state']}")
    if result["workflows_not_executed"]:
        print(f"not executed: {', '.join(result['workflows_not_executed'])}")
    if non_blocking_failed:
        print(f"non-blocking, does not gate the release: {', '.join(non_blocking_failed)}"
              " (see PANEL_MANIFEST.json non_blocking_scopes)")
    if blocking_failed:
        print(f"release-blocking failure: {', '.join(blocking_failed)}")
    print("all_release_blocking_scopes_qualified: false "
          "(panels unadopted; no second-machine reproduction recorded)")
    return 1 if blocking_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
