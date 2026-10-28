# Paper figures

The JOSS figures are rebuilt from public repository inputs and a small pin
file. The synthetic StructQC run and the public CA2 runs are regenerated into
a temporary store, then checked against the run ids and artifact SHA-256s in
the pin file. Those run directories are not committed.

The historical source revision is `21398f9a788870190385ee87397e1e488b38d09b`.
The paper-figures CI job checks out that exact revision into
`.paper-replay-source`, installs it in the pinned environment, and runs the
current figure checker against the current committed figures and pin file.
For a manual replay, make that historical checkout available at the same path:

```bash
python -m pip install --constraint yauvi-structural-workbench/docs/requirements-linux-x86_64-py312.txt -e "./.paper-replay-source[dev]" "matplotlib==3.11.2"
python tools/build_scientific_figures.py --check
```

Run ids also hash the scientific source and orchestration source. Adding a new
Workbench procedure changes that identity even when the plotted measurements
are identical. This job reproduces the historical figures with their recorded
source; it does not qualify the current application. The other reviewer,
installed-package and determinism jobs continue to test the current commit.
The original figure pins, provenance, figures and scientific receipts remain
unchanged. Do not update them merely to make a newer application match.

The pins were recorded on Linux x86_64, CPython 3.12.3, with
`yauvi-structural-workbench/docs/requirements-linux-x86_64-py312.txt`. A run id
hashes that interpreter, the CPU, and the numpy, scipy, gemmi, and biopython
versions, so another environment fails the pin check. `freesasa`, `foldseek`,
and `diamond` stay off `PATH`; the pins record them as unavailable.

Matplotlib is not part of the scientific install. `3.11.2` is the version the
committed figures were drawn with. `--check` is the CI gate. It regenerates
the runs, rebuilds into a temporary directory, and byte-compares the PNG and
SVG files with the committed copies, then checks
`paper/figures/FIGURE_PROVENANCE.json`. CA2 acquisition needs network access
to RCSB and UniProt. The downloaded coordinates are not written into the
repository.

## What each figure is built from

| Figure | Input | Run id |
|---|---|---|
| `architecture` | Schematic. No store run. | none |
| `six-evidence-questions` | Schematic. No store run. | none |
| `worked-example` | Regenerated StructQC run of `software/structqc/examples` | `run-963b421095bb0b15` |
| `ca2-coverage` | Regenerated runs of `public-ca2-structure-coverage-v1` (PDB 2POW and 4RN4, chain A, UniProt P00918) | `run-34f976fe741681be`, `run-cf802634e64a2e94` |
| `qualification-history` | Historical receipt. Not a store run. | none |

`paper/paper.md` cites `architecture`, `six-evidence-questions`, and
`qualification-history`. It does not cite `worked-example` or `ca2-coverage`.
The generator does not add or remove those references.

The pin file is
[`paper/figures/FIGURE_INPUTS.json`](../paper/figures/FIGURE_INPUTS.json).
Each store-run entry names the analysis id, the run id, and the SHA-256 of
`ANALYSIS_RUN.json`, `outputs/structqc/STRUCTURE_EVIDENCE.json`, and
`outputs/structqc/RESIDUE_QUALITY.tsv`. The builder reruns the case and
requires those digests. It does not keep the RCSB mmCIF or the evidence JSON
in the repository. The public case still checks the downloaded coordinates
against
[`public_research_sources.json`](../../examples/structural-portfolio/public_research_sources.json).

The synthetic case is analysis `synthetic-structqc`, question "What does this
two-residue fixture measure?", subject `SYNTHETIC`, parameters model `0` and
chain `A`, and the four files in `software/structqc/examples/`. The CA2 cases
use the same question, provenance document, and parameter rules as
`tools/check_determinism.py`.

## Historical receipt

`evidence/implementation-evidence/2026-09-06/` is the original worked-example
input. It is retained and it is not a run id. The residue table in that folder
matches the pinned synthetic run byte for byte. The evidence JSON there is
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

## Changing a pin

A new interpreter or library set produces a new run id. Record the new id and
the three artifact SHA-256s in `FIGURE_INPUTS.json` only after the regenerated
run is the one you intend to keep, then rerun the builder and commit the
figures and the pin file. Do not commit the run directory.

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
