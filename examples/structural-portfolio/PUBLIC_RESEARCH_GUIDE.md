# Reproduce the public CA II structure case

This case uses two public RCSB PDB entries and the UniProt P00918 human CA II
reference sequence. It analyzes chain A in PDB 2POW and 4RN4 with StructQC.
The protocol and fixed selections are in `PUBLIC_RESEARCH_PROTOCOL.md`; exact
source checksums and metadata are in `public_research_sources.json`.

The run acquires no private data and does not send analysis results anywhere.
The coordinates are downloaded at run time and are not duplicated in this
example directory. RCSB and UniProt are the primary data providers. If any
downloaded bytes differ from the lock, the acquisition step stops; update the
protocol and lock in a reviewed new case version rather than accepting the new
bytes silently.

## Recommended: run the installed wheel

From the workbench repository root, use the prepared wheel or another reviewed
wheel, and let `mktemp` choose a temporary directory outside the source tree:

```sh
CASE_TMP=$(mktemp -d)
python3.12 -m venv "$CASE_TMP/venv"
"$CASE_TMP/venv/bin/python" -m pip install build/joss-preparation-2026-09-20/yauvi_structural_workbench-0.1.0.dev0-py3-none-any.whl
"$CASE_TMP/venv/bin/python" examples/structural-portfolio/acquire_public_inputs.py \
  --out "$CASE_TMP/inputs"
"$CASE_TMP/venv/bin/python" examples/structural-portfolio/public_research_case.py \
  --inputs "$CASE_TMP/inputs" \
  --out "$CASE_TMP/run"
```

The prepared wheel is SHA-256
`b89402fa9a489455545f22bcf1a7c183b8a0285f23dc38bc48370fd1ee3e3007`.
The acquisition step needs access to the RCSB and UniProt public endpoints;
the analysis runs offline once the three locked input files are present. Keep
the temporary directory if you want to inspect its run bundle. The runner
refuses to overwrite an existing result.

## Alternate: run from the source tree

With the project dependencies available in the active Python environment, run
from the workbench repository root and point `PYTHONPATH` at the current
StructQC package:

```sh
CASE_TMP=$(mktemp -d)
python3 examples/structural-portfolio/acquire_public_inputs.py --out "$CASE_TMP/inputs"
PYTHONPATH=software/structqc/src python3 examples/structural-portfolio/public_research_case.py \
  --inputs "$CASE_TMP/inputs" \
  --out "$CASE_TMP/run"
```

Both input and output directories must remain outside the source tree. A new
run requires a new output path.

## Recorded run

The locked comparison ran locally on Python 3.12.0, Biopython 1.88, NumPy
2.5.3, and StructQC schema 1.1. The environment reported the installed root
distribution `yauvi-structural-workbench==0.1.0.dev0`; this run imported the
current pub-dev StructQC source module (`module_origin_state` is
`source_tree_or_editable`), whose SHA-256 is
`c15cb56344274806844919ba7be02b74a6ee5297741287871eb3bc9754d41ad4`. The
analysis driver SHA-256 is
`b8ce4d2a259e827f37d99225bc6b650d7c99a42fd19ec910432ba3e3c8400225` and the
acquisition script SHA-256 is
`a85a7ec8c3add651258b20e86af2c493d55ea98f6416aaee047f6c9d998d735c`. The full
output was written outside the source tree. The summary was:

| Entry | Resolved residues | Identity vs P00918 | P00918 coverage | Missing reference positions | Median residue mean B factor |
|---|---:|---:|---:|---|---:|
| 2POW, chain A | 257 | 1.000 | 0.988462 | 1–3 | 14.470 |
| 4RN4, chain A | 258 | 1.000 | 0.992308 | 1–2 | 14.197 |

These values describe the deposited coordinates and the mapping procedure. The
small coverage difference is not interpreted as a ligand effect. B factors
come from different crystal and refinement settings and are not treated as a
directly comparable quality score. StructQC imported no external geometry
validation report, and PAE was not applicable to these experimental structures.

The machine-readable run records exact coordinate, reference, protocol,
source-lock, and StructQC-core digests. It also writes the original StructQC
evidence bundle for each entry. No participant or independent researcher took
part; this demonstrates a reproducible example workflow, not independent
research use.

## Installed-wheel replay

The same locked case was replayed twice in a fresh isolated environment from
the wheel above, without `PYTHONPATH`. Both runs reported the installed
`yauvi-structural-workbench==0.1.0.dev0` distribution and
`module_origin_state=installed_distribution`; all 12 generated output files
were byte-identical across the two runs. The portable verification bundle is
[`PUBLIC_CASE_REPORT.md`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_REPORT.md),
with machine-readable results, input digests, runtime identity, and replay
verification in the neighboring `PUBLIC_CASE_*.json` files, including the
[focused safety-test record](../../evidence/preparation-2026-09-20/PUBLIC_CASE_SAFETY_TESTS.json).
This verifies a fresh installed-wheel execution; it is not independent human
research use.

## Scope limits

This case does not evaluate ligand contacts, binding affinity, inhibition,
catalysis, native-state structure, or biological function. It does not show
that one entry is better than the other. The pre-protocol 1CA2 smoke check is
excluded and not part of the source lock or reported analysis.
