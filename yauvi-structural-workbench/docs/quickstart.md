# Five-minute offline quickstart

Everything here runs offline from files bundled in the distribution. No network,
no accession fetch, no external database.

## 1. Install

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Nine console scripts land on `PATH`. Check one:

```bash
structqc describe
```

`describe` prints the module's contract as JSON: what it consumes, what it
emits, and the ceiling on what its output may be claimed to mean.

## 2. Run a complete StructQC analysis

```bash
structqc run \
  --structure structqc/examples/model.pdb \
  --reference-fasta structqc/examples/reference.fasta \
  --provenance structqc/examples/provenance.json \
  --validation-report structqc/examples/validation.json \
  --out qc-demo
```

Exit code `0` means the analysis completed. `1` means it ran but is
scientifically incomplete; `2` means the input or configuration was invalid.
An incomplete result is a real result — it is never upgraded to a favorable one.

## 3. Read what came out

`qc-demo/` contains four files:

| File | What it is |
|---|---|
| `STRUCTURE_EVIDENCE.json` | The evidence record: completeness, provenance class, imported validation, per-input digests |
| `RESIDUE_QUALITY.tsv` | Per-residue table backing the summary |
| `STRUCTURE_LAYER.json` | The typed layer other modules consume |
| `RUN_MANIFEST.json` | Inputs, parameters, runtime versions, optional runtimes, and anything that was missing |

Expected values for this input:

- `completeness.coverage_fraction` = `1.0`, `completeness.identity_fraction` = `1.0`
- `completeness.state` = `evaluated` — a reference sequence was supplied, so
  completeness could be assessed at all
- `external_validation.state` = `imported`, `clashscore` = `2.5` — **imported,
  not recomputed.** StructQC does not derive validation metrics
- `coordinate.parser.gemmi_version` records the exact parser build
- `coordinate.sha256` = `a598a5203e772bb96e4d2005f1f5cd5cb6d4d753511fce064025a258f34ac4b0`

That last digest is the anchor. It is the same coordinate file recorded in the
HUC-01 evidence record under `showcase/five-human-use-cases/`, so the run you
just did and the shipped evidence refer to provably identical input.

## 4. Confirm the output is deterministic

```bash
structqc run \
  --structure structqc/examples/model.pdb \
  --reference-fasta structqc/examples/reference.fasta \
  --provenance structqc/examples/provenance.json \
  --validation-report structqc/examples/validation.json \
  --out qc-demo-again

diff -r qc-demo qc-demo-again && echo "byte-identical"
```

The two directories are byte-identical. Determinism is what makes recording a
digest meaningful; without it every downstream checksum is noise.

## 5. See a fail-closed boundary

Ask for external validation that was never supplied:

```bash
structqc run \
  --structure structqc/examples/model.pdb \
  --require-external-validation \
  --out qc-incomplete
```

This exits `1`. `RUN_MANIFEST.json` names what was missing rather than scoring
around the gap. Missing evidence disables its leg; it never produces a
favorable value.

## Next

- [`cli-reference.md`](cli-reference.md) — every command in all nine CLIs
- [`workflows.md`](workflows.md) — the six analyses and their claim ceilings
- [`methods-and-limitations.md`](methods-and-limitations.md) — what each result may and may not be taken to mean
