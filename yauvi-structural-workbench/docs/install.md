# Installation

## Supported environment

- macOS or Linux
- Python 3.10–3.12
- An installed wheel, or a local checkout for development

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python tools/run_structural_workbench_tests.py
```

Install `.[md]` for trajectory analysis and `.[source-fetch]` for explicitly
enabled public-reference acquisition. FreeSASA, Foldseek, DIAMOND, MolProbity,
Phenix, and mkdssp are external runtimes; the workbench reports their absence
instead of silently substituting another method.

Run the packaged local server with an explicit analysis workspace:

```bash
yauvi --workspace ./my-analysis workbench serve
```

The controller binds only to loopback and uses same-origin session tokens for
file ingestion, acquisition, analysis launch, and cancellation.

`python tools/run_structural_workbench_tests.py` is the authoritative local
reviewer check. It runs suites separately to avoid legacy test-module name
collisions and deliberately excludes private platform tests that are not present
in the reviewer wheel.

## Installed, offline example

The example is bundled in the wheel and works outside the checkout:

```bash
yauvi --workspace ./my-analysis example --analysis qc-example
yauvi --workspace ./my-analysis analysis run --analysis qc-example
yauvi --workspace ./my-analysis workbench serve
```

Open the loopback address printed by the command. The synthetic example establishes
software behavior only. To exercise missing evidence, create another example with
`--without-validation`; its run must remain scientifically incomplete.

`yauvi-fetch sources` uses the packaged structural registry. An explicit `--registry`
or `YAUVI_SOURCES_REGISTRY` overrides it; private parent catalogs are not discovered.
