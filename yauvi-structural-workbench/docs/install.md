# Installation

## Supported environment

- macOS or Linux
- Python 3.10–3.12
- A local checkout of the repository

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

Run the local server from the repository root:

```bash
structqc describe
```

The controller binds only to loopback and uses same-origin session tokens for
file ingestion, acquisition, analysis launch, and cancellation.

`python tools/run_structural_workbench_tests.py` is the authoritative local
reviewer check. It runs suites separately to avoid legacy test-module name
collisions and deliberately excludes private platform tests that are not present
in the reviewer wheel.
