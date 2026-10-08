# Region explorer software verification

The Region explorer is experimental, unqualified and nonblocking. It has no
qualification panel. Existing scientific panels and thresholds are unchanged.

The offline synthetic example and regression tests check exact source bindings,
selected conformers, mapping coverage, heavy-atom contacts, report/export inclusion,
separate review records and relocated region replay. Synthetic inputs carry no
biological claims. Missing provenance and experimental validation remain explicit.

From an installed checkout, reproduce the example with:

```bash
yauvi --workspace ./region-check example --analysis region-demo --regions
yauvi --workspace ./region-check analysis validate --analysis region-demo
yauvi --workspace ./region-check analysis run --analysis region-demo
python -m pytest software/platform/tests/test_region_explorer.py -q
```

The run is expected to return exit code 1 because StructQC is scientifically
incomplete; region processing returns 0. Inspect both recorded steps rather than
treating the complete pipeline as clean QC. Inputs, mapping, conformers and
represented contacts remain traceable in the report and evidence export.

Verification on 2026-10-08 used the installed distribution on macOS arm64 with
Python 3.12.14. The 31 synthetic region regression checks passed, and the installed
example retained incomplete QC with successful region processing. This is local
software verification, not second-machine or independent-researcher reproduction.

See [Regions, review and replay](REGION_EXPLORER.md) for inputs and interpretation
limits. Region replay recalculates mapping, atom selection and represented contact
measurements; other engine outputs are checksum verification only.
