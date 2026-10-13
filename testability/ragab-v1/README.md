# RagAB testability kit v1

Reproduce the retained RagAB structural-resemblance and protein-interface
measurements in a new output folder, with explicit failures and fixed comparison
rules. This kit is a candidate for independent use, not evidence that another
researcher or computer has reproduced the study.

## What goes in → what happens → what comes out

Six checksummed public structures and their validation reports, fixed chain
pairings and residue masks, three frozen engine snapshots and a declared runtime
go in. The kit re-prepares coordinates, runs the twelve original case controls,
then checks ten state calculations and twelve interface calculations. Every
worker must exit successfully and provide all required outputs. Forty-six
comparison groups and separate execution records come out.

The primary stopping rule requires Q1 open-like and Q2 closed-like before the
eight secondary queries run. A mismatch, missing file or failed worker exits
nonzero and preserves the failed attempt. Existing destinations are refused.

## Setup

The v1 profile is **Python 3.12**, the exact dependency versions in
`requirements.txt`, and a working native **FreeSASA 2.1.3** executable. Python,
FreeSASA and dependencies are not bundled. Other profiles must be separately
versioned before running; do not edit this profile to turn a failed attempt into
a pass. Only macOS arm64 has a local test route; Windows/Linux and independent
installation remain untested. The old Intel runtime remains a recorded failure
and is not replaced by this profile.

An operator who chooses to install dependencies can use a separate environment:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

The three engines are loaded from their frozen source copies in `engines/`;
there is no need to install the current Workbench package over them. Apache-2.0
licenses are retained with each source copy.

## Obtain and check the inputs

No source coordinates or validation-report bytes are shipped in the Git payload.
`SOURCES.json` carries their original provider URLs and SHA-256 values. An
operator explicitly choosing public downloads can run:

```sh
.venv/bin/python acquire_sources.py --download-public-sources
```

Provider bytes may have changed. Different bytes are refused, never adopted
automatically; use the exact historical source version from the provider or
register a different study. This acquisition route has not been tested live in
this preparation. Replay itself never downloads, installs or uploads anything.

```sh
.venv/bin/python -B replay.py verify-kit
.venv/bin/python -B replay.py check --freesasa /path/to/freesasa
.venv/bin/python -B test_replay.py
```

`verify-kit` checks shipped identities only. `check` also checks source bytes,
dependency versions and whether the chosen executable runs; it does no
scientific calculation. Run outputs contain local execution provenance and
must remain private.

## Run and inspect

Choose a new output directory; replace the executable placeholder with the
actual installed native FreeSASA path:

```sh
.venv/bin/python -B replay.py run --freesasa /path/to/freesasa --output runs/operator-001
.venv/bin/python -B replay.py compare --output runs/operator-001
```

Before calculation, `REGISTRATION.json` binds kit, sources, runtime and
comparison policy. `RESULT.json` reports the attempt status. `COMPARISON.json`
reports each of the ten state, twelve interface-summary and twenty-four raw
interface-record comparisons. The read-only `compare` command recalculates the
comparisons from retained measurements; it does not rerun the engines.

The surface-area absolute tolerance is 0.0001 square angstrom; other float
fields use 0.000001, with zero relative tolerance. Discrete fields are exact.
Only the raw input QC-manifest hash is omitted because regenerated QC provenance
depends on the execution context. Coordinate bytes must match their expected
hashes exactly. These are numerical limits, not biological uncertainty estimates.

## Results and scientific limits

Expected outputs are disclosed in `EXPECTED.json`; this is not a blinded test.
The retained primary contrast is 1,594.6596595 square angstrom more modeled
interface area for the closed-like pair. Six of the eight secondary pairs have
retained terminal-extent disagreements with the closed reference. Repeating
these calculations cannot establish transport activity, immune recognition or
biological validation. Symmetry-related copies are not independent observations.

The kit starts from deposited structural models, uses joint fitting that gives
RagA more weight, excludes the fixed signature positions from the scoring mask,
and omits bound peptides and nonprotein components from interface measurement.
Labels remain case-specific operational conventions, outside the qualified
StateAtlas exact-mapping scope. No Workbench qualification or submission gate
is cleared by a numerical replay.

## What another operator should record

Keep the run and its source registration locally. Record who ran it, whether
they helped develop the study, whether it is a different physical computer,
installation difficulties and their interpretation. The runner does not infer
independence from a successful exit. Agree the exact report and destination
before sharing run records.

## Files

- `replay.py`: new preflight, registration, worker checks, stopping and comparison.
- `test_replay.py`: synthetic failure checks; no scientific measurements.
- `frozen-adapter/`: unchanged original preparation, measurement and case controls.
- `engines/`: unchanged historical StructQC, StateAtlas and AssemblyContext code.
- `CONFIG.json`, `SOURCES.json`: fixed scientific choices and source identities.
- `EXPECTED.json`, `RUNTIME_PROFILE.json`, `ORIGIN.json`: fixed endpoints,
  runtime adoption and original protocol identities.
- `PAYLOAD.json`: exact outgoing member list and hashes.
- `sources/`, `runs/`, `.venv/`: local-only, ignored, never copied by the handoff.

Do not invoke the frozen adapter's old `run` or `replay` helpers directly. The
new runner deliberately calls and verifies each calculation leg independently.
