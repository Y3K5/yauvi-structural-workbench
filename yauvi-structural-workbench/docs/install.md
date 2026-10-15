# Installation

## Supported environment

- macOS or Linux
- Python 3.10–3.13
- An installed wheel, or a local checkout for development

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python tools/run_structural_workbench_tests.py
```

## Pinned dependencies

The command above installs whatever currently satisfies the lower bounds in
`pyproject.toml`. Those bounds are unchanged by the files in this directory.
A constraint file pins one resolution for one CI-matrix environment so a later
install can repeat it. Pass it with `--constraint`, and install the checkout
or wheel separately:

```bash
python -m pip install --constraint yauvi-structural-workbench/docs/requirements-linux-x86_64-py312.txt -e ".[dev]"
```

Choose the file whose name matches the machine (`uname -s`, `uname -m`) and
the interpreter minor version:

- `requirements-linux-x86_64-py310.txt`, `py311`, `py312`, and `py313` — Linux x86_64
- `requirements-macos-arm64-py310.txt`, `py311`, and `py313` — macOS arm64
- `requirements-macos-arm64-py312.txt` — macOS arm64, Python 3.12. This is the previously verified environment (Python 3.12.0), not a fresh resolution
- `requirements-linux-aarch64-py312.txt` — Linux aarch64, the CI architecture job

Files other than the verified macOS arm64 / Python 3.12 pin were resolved from
PyPI on 2026-10-15 for that platform tag and CPython minor. For each minor,
the pin list is the same on every architecture that has a file. Linux aarch64
is recorded only for 3.12, which is the architecture job in CI. The files stay
separate so a later platform-specific wheel can be pinned without renaming the
series. On Linux x86_64, `pip install --constraint` of the Python 3.12 file
reproduced that file under `pip freeze`. The Python 3.10 and 3.11 files match
`pip freeze` from an install on those interpreters. CI does not install
through these files. It installs from the open bounds and records `pip freeze`
in the job log and a `resolved-deps-*` artifact.

Install `.[md]` for trajectory analysis and `.[source-fetch]` for explicitly
enabled public-reference acquisition. FreeSASA, Foldseek, DIAMOND, MolProbity,
Phenix, and mkdssp are external runtimes; the workbench reports their absence
instead of silently substituting another method.

Run the packaged local server with an explicit analysis workspace:

```bash
yauvi --workspace ./my-analysis workbench open
```

The controller binds only to loopback and uses same-origin session tokens for
file ingestion, acquisition, analysis launch, and cancellation.

The browser opens on port 8947 by default. Its footer identifies the exact build;
click it for version details. `workbench serve` starts without opening a browser.
See [Launch and versions](LAUNCH_AND_VERSIONS.md) for multiple local copies.

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
