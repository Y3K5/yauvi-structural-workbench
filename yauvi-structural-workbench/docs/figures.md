# Paper figures

The JOSS figures are rebuilt from a manifest of store run ids and historical
receipts. They are not read out of a frozen evidence folder.

```bash
python -m pip install "matplotlib==3.11.2"
python tools/build_scientific_figures.py
python tools/build_scientific_figures.py --check
```

Matplotlib is not part of the scientific install. `3.11.2` is the version the
committed figures were drawn with. `--check` is the CI gate. It rebuilds into
a temporary directory and byte-compares the PNG and SVG files with the
committed copies, then checks `paper/figures/FIGURE_PROVENANCE.json`.

## What each figure is built from

| Figure | Input | Run id |
|---|---|---|
| `architecture` | Schematic. No store run. | none |
| `six-evidence-questions` | Schematic. No store run. | none |
| `worked-example` | Recorded StructQC run of `software/structqc/examples` | `run-963b421095bb0b15` |
| `ca2-coverage` | Recorded runs of `public-ca2-structure-coverage-v1` (PDB 2POW and 4RN4, chain A, UniProt P00918) | `run-34f976fe741681be`, `run-cf802634e64a2e94` |
| `qualification-history` | Historical receipt. Not a store run. | none |

`paper/paper.md` cites `architecture`, `six-evidence-questions`, and
`qualification-history`. It does not cite `worked-example` or `ca2-coverage`.
The generator does not add or remove those references.

The manifest is
[`paper/figures/FIGURE_INPUTS.json`](../paper/figures/FIGURE_INPUTS.json).
Each store-run entry names the analysis id, the run id, and the SHA-256 of
every artifact the figure reads. Those artifacts are copied under
`yauvi-structural-workbench/paper/figures/recorded-runs/<analysis_id>/runs/<run_id>/`. The copy is the
StructQC evidence, the residue table, and `ANALYSIS_RUN.json`. It is not the
RCSB mmCIF. The public case still checks those outputs against
[`public_research_sources.json`](../../examples/structural-portfolio/public_research_sources.json).

## Historical receipt

`evidence/implementation-evidence/2026-09-06/` is the original worked-example
input. It is retained and it is not a run id. The residue table in that folder
matches the recorded synthetic run byte for byte. The evidence JSON there is
the earlier schema 1.1 export. Qualification history is a separate receipt:
the v2.11 panel manifest, the retained single-machine
`EXECUTION_SUMMARY.json`, the 2026-09-20 reconciliation, and
`CI_OBSERVATION.json`. Those records cannot be regenerated from the public
CA2 case. The figure cites the actions run id stored in the observation
(`35295451246`). That number is GitHub Actions job metadata, not a
structural-workbench run id.

`FIGURE_PROVENANCE.json` is schema 2.0. For every figure it records the input
class, the run ids (empty when there is no store run), the artifact SHA-256s,
the plotted data, and the generator script hash.

## A fresh store run

Run ids include the interpreter, OS, and library versions, so a new execution
of the same case gets a new run id. Point the builder at that workspace:

```bash
python tools/build_scientific_figures.py \
  --store /path/to/workspace \
  --run 2pow=run-0123abcd0123abcd \
  --run 4rn4=run-0123abcd0123abcd
```

Named analyses are read from the workspace
(`structural_analyses/cases/<analysis_id>/runs/<run_id>/`). Figures you do not
name keep the recorded runs. This command does not rewrite
`FIGURE_INPUTS.json`. `--check` without `--store` or `--run` uses the
manifest. To change the committed receipt, copy the cited artifacts into
`yauvi-structural-workbench/paper/figures/recorded-runs/`, put the new run ids and SHA-256s in the manifest, and rerun
the builder.

## Why PNG bytes can differ

Drawing uses the Agg backend. SVG keeps text as text (`svg.fonttype` is
`none`), clip-path ids use the fixed salt `yauvi-paper-figures`, and
`SOURCE_DATE_EPOCH=1700000000` fixes the embedded date. Two builds on one
machine are byte-identical.

PNG still rasterizes text with FreeType and compresses with libpng, so another
platform can produce different PNG bytes from the same drawing commands. When
that happens, `--check` compares the SVG and the extracted plotted data (the
sentences and numbers recorded in the provenance) and prints the reason. A
change in that plotted data fails the check.
