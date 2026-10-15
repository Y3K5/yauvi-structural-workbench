# Reproducibility

Each analysis stores input bytes by SHA-256 and emits:

- `REPORT_DATA.json`
- `REPORT.html`
- `RAW_EVIDENCE.zip`
- `CHECKSUMS.json`
- `RUN_MANIFEST.json`

The raw ZIP uses fixed entry order and normalized archive metadata. Scientific
JSON/TSV/HTML excludes absolute workspace paths. Logs scrub the workspace and
source-tree roots. External runtime versions, reference pack checksums,
parameters, missing evidence, and non-claims remain visible.

## Determinism

No step of the public case uses a random number generator, so no seed is
needed and none is recorded. Do not add a seed field to the case.

The case is `public-ca2-structure-coverage-v1` (PDB 2POW and 4RN4, chain A,
against UniProt P00918). `tools/check_determinism.py` runs it through
`StructuralAnalysisStore` twice, into two workspaces, and byte-compares every
emitted file. That includes `REPORT_DATA.json`, `CHECKSUMS.json`,
`RAW_EVIDENCE.zip`, `REPORT.html`, `RUN_MANIFEST.json`, and the StructQC
evidence written beside them. No file is left out of the compare. The script's
within-run exclusion list is empty: a wall-clock timestamp or an absolute path
that actually varied would have to be named there, with a reason, and a
difference that is not named fails the run.

The same command searches `software/**/src` and `software/**/examples` for
`random`, `np.random`, `default_rng`, `os.urandom`, and `secrets` token draws.
An unlisted hit fails the job. The public-case path has no generator. The hits
that do exist, and why they are not a seed for this case:

- MembraneOrient's rotation-invariance check calls `numpy.random.default_rng(0)`.
  The `0` is a constant written in that function, not an input. This case does
  not run MembraneOrient.
- `software/membrane-orientor/memorient/examples/p4_benchmark.py` is an example
  script with a local seed argument. It is not a store step.
- Ingest upload ids use `os.urandom`. `finalize_ingest` deletes the ingest
  directory and does not store the id. Repeated case runs are byte-identical.
- Browser tokens, job ids, import ids, and source-acquisition ids use
  `secrets`. This case does not read them.

Uses of the word "seed" in the barrel search and the assembly frame are
geometric starting vectors, not random draws. Test and fixture builders are
outside this search; they are not case steps.

CI runs that command on `ubuntu-latest` and `macos-latest` (arm64) for Python
3.10–3.13, uploads each workspace, and diffs them in a final job. A
cross-runner difference is accepted only when every changed value is one of
these fields, which the log prints with the values that were observed:

- `runtime_versions` — interpreter, OS, CPU, and the resolved numpy, scipy,
  gemmi, biopython, freesasa, foldseek, and diamond versions. The run id is a
  hash of an identity that includes this object.
- `run_id` — that hash, also used as the run directory name. Not a timestamp.
  A string whose only difference is the `run-<16 hex>` token is the same
  field. `runs[].path` (`runs/run-<id>/ANALYSIS_RUN.json`) is one of those.
- `latest_run_id` — the case manifest's copy of `run_id`.
- `revision_sha256` — SHA-256 of the case manifest with this field removed.
  Accepted only when every other difference in that document is one of the
  fields above. The content-addressed file
  `revisions/<revision_sha256>.json` is that hash, so the final revision's
  filename differs for the same reason. A hash that changes while the rest of
  the document is identical fails. A scientific field that also differs fails.
  Revisions written before the run id exists stay byte-identical.
- `gemmi_version` — the parser build string. Coordinate numbers stay rounded
  to 6 decimal places; a coordinate or B-factor change is not this field.
- `optional_runtimes.mkdssp` and `optional_runtimes.molprobity_or_phenix` —
  whether those optional executables are on `PATH`. This case does not invoke
  them. `optional_runtimes.gemmi` is not in this list.

`CHECKSUMS.json` and `RAW_EVIDENCE.zip` may differ only because they bind
files that differ in the fields above. Any other difference, including float
text that is not one of those fields, fails the job and names the file. The
store's structure QC adapter passes `--require-external-validation`. This case
does not supply a validation report, so both runs are scientifically
incomplete and must record the same status.

For a reviewer reproduction:

```bash
python -m pip install -e ".[dev]"
python tools/run_structural_workbench_tests.py
structqc describe
```

That install follows the lower bounds in `pyproject.toml`, so it is not pinned.
To repeat one recorded resolution, pass the constraint file for the machine as
documented in [Installation](install.md#pinned-dependencies). From the
repository root:

```bash
python -m pip install --constraint yauvi-structural-workbench/docs/requirements-linux-x86_64-py312.txt -e ".[dev]"
```

CI still installs from the open bounds and uploads the versions that install
actually resolved.

Network tests and live public-source smoke tests are opt-in. A genuine offline
demonstration must use bundled synthetic fixtures; cached or remote retrieval is
not treated as offline.

## Paper figures

Regenerate the JOSS figures with `python tools/build_scientific_figures.py`.
The command reruns the synthetic StructQC fixture and the public CA2 case,
then checks those runs against the ids and artifact SHA-256s in
`paper/figures/FIGURE_INPUTS.json`. The run directories are not committed.
The pins match Linux x86_64, CPython 3.12.3, and
`docs/requirements-linux-x86_64-py312.txt`. Qualification history is not a
store run. `evidence/implementation-evidence/2026-09-06/` stays as the worked
example's historical receipt. The command, the pin file, and the PNG
comparison rule are in [Paper figures](figures.md).

## Independent public qualification

The historical v1 scientific panel is separate from the synthetic software
tests. Its exact named-case result remains reproducible offline with:

```bash
python evidence/benchmarks/qualification-v1/run_qualification.py
```

The current Mark 1 release gate is Qualification v2. Audit its frozen strata,
splits, source adoption, and exact thresholds with:

```bash
python evidence/benchmarks/qualification-v2/run_qualification.py
```

The v2 runner makes no network requests and currently returns `1` with
`blocked_panel_incomplete`. No v2 scientific execution has occurred. Existing
v1 files remain candidate evidence only until each record receives its own v2
release, checksum, license, citation, split, expected result, mapping, and
exclusion rationale. Repeated v2 audits produce byte-identical JSON, TSV, HTML,
and checksum records.
