# Reviewer quickstart

This path demonstrates software behavior and independent scientific case
evidence without requiring private data or a network connection. Passing
software tests and passing scientific qualification cases are separate claims.

## Install

From a clean checkout and Python 3.10-3.12 environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Run the release-scoped offline tests:

```bash
python -m pip install ".[dev]"
python tools/run_structural_workbench_tests.py
```

## Inspect the six workflows

```bash
structqc describe
memorient describe
state-atlas describe
site-context describe
assembly-context describe
sf-csa describe
```

Each command reports its inputs, outputs, scientific boundary, and optional
runtimes. The browser workbench generates manifests; reviewers do not need to
hand-author JSON.

## Run the synthetic showcase

```bash
python tools/build_five_use_case_showcase.py --replace
python tools/verify_five_use_case_showcase.py
python tools/verify_public_showcase.py
```

The synthetic cases demonstrate deterministic execution and fail-closed
boundaries. They do not establish accuracy against independent biology.

## Run the public qualification panel

After acquiring the exact artifacts in `benchmarks/qualification-v1/SOURCE_LOCK.json`:

```bash
python yauvi-structural-workbench/benchmarks/qualification-v1/run_qualification.py
```

The v1 runner performs no network access. Its current expected state is exit
`1`: four named cases pass and two remain partial. These are historical
qualification records, not the current Mark 1 release gate.

Audit the expanded Qualification v2 panel separately:

```bash
python yauvi-structural-workbench/benchmarks/qualification-v2/run_qualification.py
```

The expected state is `blocked_panel_incomplete`. No v2 scientific execution
has occurred and no missing case is counted as favorable evidence.

## Start the local interface

```bash
structqc describe
```

Open `/public-showcase/` for the narrative and `/#new` for the task-first
analysis builders. Nothing is uploaded or fetched unless a separately enabled,
registered public-source action is explicitly requested.

## Evidence to inspect

- `BASELINE.json`: recorded software regression selection.
- `RELEASE_STATUS.json`: authoritative readiness and approval state.
- `JOSS_PUBLICATION_ROADMAP.json`: remaining gates and relative sequence.
- `benchmarks/qualification-v1/QUALIFICATION_REPORT.md`: readable results.
- `benchmarks/qualification-v1/results/QUALIFICATION_RESULTS.json`: exact checks.
- `benchmarks/qualification-v1/SOURCE_LOCK.json`: public artifact identities.
- `benchmarks/qualification-v2/PANEL_MANIFEST.json`: frozen scope, strata,
  splits, evidence requirements, and tolerances.
- `benchmarks/qualification-v2/results/QUALIFICATION_V2_STATUS.json`: exact
  source-adoption and panel-completeness block.

No showcased result establishes biochemical activity, native membrane exposure,
binding affinity, observed catalysis, exact functional transfer, clinical
utility, or permission to publish.
