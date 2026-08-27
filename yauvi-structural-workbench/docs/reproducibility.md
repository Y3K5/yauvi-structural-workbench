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

For a reviewer reproduction:

```bash
python -m pip install -e ".[dev]"
python tools/run_structural_workbench_tests.py
structqc describe
```

Network tests and live public-source smoke tests are opt-in. A genuine offline
demonstration must use bundled synthetic fixtures; cached or remote retrieval is
not treated as offline.

## Independent public qualification

The historical v1 scientific panel is separate from the synthetic software
tests. Its exact named-case result remains reproducible offline with:

```bash
python yauvi-structural-workbench/benchmarks/qualification-v1/run_qualification.py
```

The current Mark 1 release gate is Qualification v2. Audit its frozen strata,
splits, source adoption, and exact thresholds with:

```bash
python yauvi-structural-workbench/benchmarks/qualification-v2/run_qualification.py
```

The v2 runner makes no network requests and currently returns `1` with
`blocked_panel_incomplete`. No v2 scientific execution has occurred. Existing
v1 files remain candidate evidence only until each record receives its own v2
release, checksum, license, citation, split, expected result, mapping, and
exclusion rationale. Repeated v2 audits produce byte-identical JSON, TSV, HTML,
and checksum records.
