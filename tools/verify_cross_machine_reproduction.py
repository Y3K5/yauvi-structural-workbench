#!/usr/bin/env python3
"""Combine the per-runner execution evidence into one cross-machine verdict.

`summarize_execution.py` reduces *one* runner's results to one document, and
hardcodes `second_machine_reproduction: "not_recorded"` because from inside a
single runner that is the only honest answer. The qualification matrix then
produces six of those summaries, uploads six archives, and nothing has ever
combined them. So the evidence for cross-machine reproduction has been produced
weekly since the matrix existed and read by nobody, which is the same failure
the summarizer itself was written to fix one level down.

This reads the uploaded archives back and writes
`results/CROSS_MACHINE_REPRODUCTION.json`: which environments actually ran, and
whether each release-blocking panel produced the *same* outcome on more than one
of them.

What counts as a second machine
-------------------------------

Two runners are a different **environment** when they differ in `platform` or
`machine` -- a different operating system, or a different CPU architecture.
A different Python version on the same OS and architecture is recorded, and is
worth having, but does not on its own make a second machine: the interpreter
changed, the hardware and the system libraries did not, and the failures this
gate exists to catch (an endianness assumption, a libm difference, a BLAS
kernel, a compiler's floating-point contraction) do not vary with it.

That is the strict reading, and it is deliberate. `qualification-v2` currently
matrixes ubuntu-latest (x86_64) against macos-latest (arm64 since macos-14), so
under this definition the matrix already spans two environments on two
architectures. Whether *GitHub runners under one account* satisfy "independent"
in the sense the release gate means is a separate question this tool does not
answer and must not appear to: it reports environments and agreement, and
`independent_second_machine_reproduction_passed` stays a human decision recorded
in RELEASE_STATUS.json.

Agreement, not just success
---------------------------

Two runners both reporting "passed" is necessary but not sufficient. A panel
that passes 16/16 on one environment and 14/16 on another has *not* reproduced,
even if some looser gate calls both acceptable, so the comparison is over the
full outcome tuple -- case counts, control counts, stratum state -- and any
difference is reported as a disagreement rather than folded into a pass.

Disagreement is the finding worth having. The membrane panel is expected to
differ across environments (its drift deltas are why the matrix exists at all)
and is non-blocking by collection 2.4, so it is reported and never gates.

Reads only `EXECUTION_SUMMARY.json` from each archive, never the per-case
results: the summary is the surface that already excludes absolute paths, and
each case's `output_dir` names the machine that ran it.

Usage:
  verify_cross_machine_reproduction.py --evidence-dir DIR [--json-out PATH]
                                       [--min-environments N]
Exit:
  0  every executed release-blocking panel agreed across >= N environments
  1  a release-blocking panel disagreed between environments
  2  not enough evidence to say -- fewer than N environments reported

1 and 2 are separated deliberately, and CI must not collapse them. A
disagreement is a finding about the software. Missing evidence is usually a
transfer that GitHub could not complete: uploads here are `continue-on-error`
precisely because the artifact service is not evidence about the software, and
a dropped upload arriving as a red reproduction gate would say something false.
`verify_source_lock_health.py` draws the same line, exiting 0 when a provider is
merely unreachable.
"""
from __future__ import annotations

import argparse
import json
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
QUALIFICATION = ROOT / "yauvi-structural-workbench" / "benchmarks" / "qualification-v2"

# The outcome fields compared between environments. Anything that differs here
# means the panel did not reproduce, whatever the individual verdicts say.
OUTCOME_FIELDS = ("stratum_state", "cases_passed", "cases_failed", "cases_total",
                  "controls_passed", "controls_total")


def environment_key(runtime: dict[str, Any]) -> str:
    """The axis a second machine has to differ on. Python version is not it."""
    return f"{runtime.get('platform', '?')}/{runtime.get('machine', '?')}"


def runner_label(runtime: dict[str, Any], source: str) -> str:
    return (f"{environment_key(runtime)} py{runtime.get('python', '?')}"
            if runtime else f"unknown ({source})")


def outcome(panel: dict[str, Any]) -> dict[str, Any]:
    cases = panel.get("cases", {}) or {}
    controls = panel.get("controls", {}) or {}
    return {
        "stratum_state": panel.get("stratum_state"),
        "cases_passed": cases.get("passed", 0),
        "cases_failed": cases.get("failed", 0),
        "cases_total": cases.get("total", 0),
        "controls_passed": controls.get("passed", 0),
        "controls_total": controls.get("total", 0),
    }


def load_summaries(evidence_dir: Path) -> list[dict[str, Any]]:
    """Every EXECUTION_SUMMARY.json reachable under `evidence_dir`.

    Accepts the archives `actions/download-artifact` leaves behind, and also an
    already-extracted tree, so the same tool works in CI and by hand.
    """
    found: list[dict[str, Any]] = []

    for path in sorted(evidence_dir.rglob("EXECUTION_SUMMARY.json")):
        try:
            found.append({"source": str(path.relative_to(evidence_dir)),
                          "summary": json.loads(path.read_text(encoding="utf-8"))})
        except (OSError, ValueError) as exc:
            print(f"  unreadable: {path.name}: {exc}", file=sys.stderr)

    for archive in sorted(evidence_dir.rglob("*.tar.gz")):
        try:
            with tempfile.TemporaryDirectory() as tmp, tarfile.open(archive, "r:gz") as tar:
                members = [m for m in tar.getmembers()
                           if m.name.endswith("EXECUTION_SUMMARY.json")]
                if not members:
                    continue
                # filter="data" refuses absolute paths, traversal and device
                # nodes. Added in 3.12; on older interpreters the members are
                # selected by name above and extracted into a temp directory
                # that is discarded, so the blast radius is bounded either way.
                try:
                    tar.extractall(tmp, members=members, filter="data")
                except TypeError:
                    tar.extractall(tmp, members=members)
                for extracted in Path(tmp).rglob("EXECUTION_SUMMARY.json"):
                    found.append({"source": archive.name,
                                  "summary": json.loads(extracted.read_text(encoding="utf-8"))})
        except (OSError, tarfile.TarError, ValueError) as exc:
            print(f"  unreadable archive: {archive.name}: {exc}", file=sys.stderr)

    return found


def build_report(loaded: list[dict[str, Any]], min_environments: int) -> dict[str, Any]:
    runners: list[dict[str, Any]] = []
    for entry in loaded:
        summary = entry["summary"]
        recorded = summary.get("recorded_on") or [{}]
        # A summary carrying more than one runtime was not produced by a single
        # run; it is not a runner and cannot speak for an environment.
        if len(recorded) != 1:
            runners.append({"source": entry["source"], "usable": False,
                            "why": f"summary spans {len(recorded)} runtimes; not a single run"})
            continue
        runtime = recorded[0]
        runners.append({
            "source": entry["source"],
            "usable": True,
            "environment": environment_key(runtime),
            "label": runner_label(runtime, entry["source"]),
            "python": runtime.get("python"),
            "panels": {p.get("workflow"): {**outcome(p),
                                           "release_blocking": bool(p.get("release_blocking"))}
                       for p in summary.get("panels", [])},
        })

    usable = [r for r in runners if r["usable"]]
    environments = sorted({r["environment"] for r in usable})

    workflows = sorted({w for r in usable for w in r["panels"]})
    panels: list[dict[str, Any]] = []
    for workflow in workflows:
        by_environment: dict[str, list[dict[str, Any]]] = {}
        blocking = False
        for runner in usable:
            panel = runner["panels"].get(workflow)
            if panel is None:
                continue
            blocking = blocking or panel["release_blocking"]
            by_environment.setdefault(runner["environment"], []).append(
                {"label": runner["label"], **{k: panel[k] for k in OUTCOME_FIELDS}})

        # One representative outcome per environment, plus any environment that
        # disagreed with itself across Python versions -- which is its own
        # finding and must not be averaged away.
        per_environment: dict[str, Any] = {}
        internally_inconsistent: list[str] = []
        for env, observations in by_environment.items():
            distinct = {tuple(o[k] for k in OUTCOME_FIELDS) for o in observations}
            if len(distinct) > 1:
                internally_inconsistent.append(env)
            per_environment[env] = {
                "runners": [o["label"] for o in observations],
                "outcomes": [dict(zip(OUTCOME_FIELDS, d)) for d in sorted(distinct)],
            }

        distinct_outcomes = {tuple(o[k] for k in OUTCOME_FIELDS)
                             for observations in by_environment.values() for o in observations}
        agreed = len(distinct_outcomes) == 1 and not internally_inconsistent
        passed_everywhere = all(
            o["stratum_state"] == "passed"
            for observations in by_environment.values() for o in observations)

        panels.append({
            "workflow": workflow,
            "release_blocking": blocking,
            "environments": sorted(by_environment),
            "environment_count": len(by_environment),
            "agreed_across_environments": agreed,
            "passed_on_every_environment": passed_everywhere,
            "internally_inconsistent_environments": sorted(internally_inconsistent),
            "reproduced": bool(agreed and passed_everywhere
                               and len(by_environment) >= min_environments),
            "per_environment": per_environment,
        })

    blocking_panels = [p for p in panels if p["release_blocking"]]
    reproduced = [p for p in blocking_panels if p["reproduced"]]
    disagreed = [p for p in blocking_panels if not p["agreed_across_environments"]]
    thin = [p for p in blocking_panels
            if p["agreed_across_environments"] and p["environment_count"] < min_environments]

    return {
        "schema_version": "1.0",
        "min_environments": min_environments,
        "environment_definition": (
            "Two runners are different environments when they differ in platform or machine. "
            "A different Python version on the same platform and machine is recorded but is "
            "not a second machine."
        ),
        "runners_read": len(runners),
        "runners_usable": len(usable),
        "environments_observed": environments,
        "environment_count": len(environments),
        "panels": panels,
        "release_blocking_panels_reproduced": sorted(p["workflow"] for p in reproduced),
        "release_blocking_panels_disagreed": sorted(p["workflow"] for p in disagreed),
        "release_blocking_panels_under_minimum": sorted(p["workflow"] for p in thin),
        "every_release_blocking_panel_reproduced": bool(blocking_panels) and not disagreed and not thin,
        # This file reports agreement between runners. Whether runners rented
        # from one provider under one account are "independent" in the sense
        # the release gate means is a judgement, and belongs to a person.
        "independence_note": (
            "Agreement across environments is measured here. Whether those environments are "
            "independent in the sense the release gate requires is not a property of this "
            "evidence and is not decided by this tool."
        ),
        "unusable_runners": [r for r in runners if not r["usable"]],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--evidence-dir", type=Path, required=True,
                        help="Directory holding the downloaded per-runner archives.")
    parser.add_argument("--json-out", type=Path,
                        default=QUALIFICATION / "results" / "CROSS_MACHINE_REPRODUCTION.json")
    parser.add_argument("--min-environments", type=int, default=2,
                        help="Environments a release-blocking panel must agree across (default 2).")
    args = parser.parse_args(argv)

    if not args.evidence_dir.is_dir():
        print(f"missing evidence directory: {args.evidence_dir}", file=sys.stderr)
        return 1

    loaded = load_summaries(args.evidence_dir)
    if not loaded:
        print(f"no EXECUTION_SUMMARY.json found under {args.evidence_dir}", file=sys.stderr)
        return 1

    report = build_report(loaded, args.min_environments)

    print(f"read {report['runners_usable']} runner summaries across "
          f"{report['environment_count']} environment(s): "
          f"{', '.join(report['environments_observed'])}\n")
    for panel in report["panels"]:
        # Four distinct states, because collapsing them mislabels the panel the
        # matrix exists for: membrane agreeing on two environments and failing
        # its gate on both is neither a disagreement nor a thin sample.
        if panel["reproduced"]:
            mark = "reproduced"
        elif not panel["agreed_across_environments"]:
            mark = "DISAGREED"
        elif not panel["passed_on_every_environment"]:
            mark = "agreed, did not pass"
        else:
            mark = f"only {panel['environment_count']} environment(s)"
        scope = "blocking" if panel["release_blocking"] else "non-blocking"
        print(f"  {panel['workflow']:<22} {scope:<12} {panel['environment_count']} env  {mark}")
        if not panel["agreed_across_environments"]:
            for env, detail in sorted(panel["per_environment"].items()):
                for observed in detail["outcomes"]:
                    print(f"      {env:<18} {observed['stratum_state']:<8} "
                          f"cases {observed['cases_passed']}/{observed['cases_total']}"
                          f"  controls {observed['controls_passed']}/{observed['controls_total']}")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        print(f"\n  written to {args.json_out}")

    if report["release_blocking_panels_disagreed"]:
        print(f"\n{len(report['release_blocking_panels_disagreed'])} release-blocking panel(s) "
              f"disagreed between environments: "
              f"{', '.join(report['release_blocking_panels_disagreed'])}")
        print("A panel that produces different counts on different machines has not reproduced. "
              "This is the finding the matrix exists to surface; do not average it away.")
        return 1
    if not report["panels"]:
        print("\nNo panels found in the evidence.")
        return 2
    if report["release_blocking_panels_under_minimum"]:
        # Not a finding about the software. Every panel that reported agreed;
        # there were simply too few environments to call it reproduction, which
        # is almost always an upload that did not arrive.
        print(f"\nNo disagreement found. Only {report['environment_count']} environment(s) "
              f"reported, so {', '.join(report['release_blocking_panels_under_minimum'])} "
              f"cannot be called reproduced at --min-environments {args.min_environments}. "
              "That is missing evidence, not a failed reproduction.")
        return 2
    print(f"\nEvery release-blocking panel agreed across {report['environment_count']} "
          f"environments. Independence remains a human judgement; see independence_note.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
