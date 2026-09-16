# memorient

**Context-aware membrane orientation and per-residue accessibility for protein structures.**

Given a predicted or experimental structure, `memorient` estimates how the protein sits in a
membrane and labels every residue: which way is "out", where the bilayer is, and which residues
face lipid, lumen, or the context-declared outside region. These are geometric assignments,
not evidence of native antibody accessibility. It is
the geometry-and-accessibility stage of a reverse-vaccinology pipeline, generalized to stand
on its own.

The biological assumptions live in a **context object**, not in code branches, so one call
handles a gram-negative outer-membrane β-barrel, a eukaryotic single-pass receptor, or a
soluble secreted protein — each gets exactly the metrics that are meaningful for it.

## Install

```bash
pip install -e ".[compute]"     # numpy, scipy, biopython, matplotlib
```

`contexts.py` and the CLI's `contexts`/`describe` subcommands are stdlib-only; the numeric
pipeline needs the `compute` extra.

## Quickstart

```python
from memorient.geometry import load_structure
from memorient.contexts import get_context
from memorient.orientor import orient_structure

s = load_structure("1BXW.pdb", chain="A")
result = orient_structure(s, get_context("gram_negative_om"))

print(result.label)                     # 'barrel'
print(result.summary()["n_extracellular"])
print(result.extracellular_resids())    # antibody-facing loop residues
result.write_pdb("oriented.pdb")        # +Z is extracellular, core at origin
```

Command line:

```bash
python -m memorient contexts
python -m memorient describe gram_negative_om
python -m memorient orient 1BXW.pdb --context gram_negative_om --chain A \
    --out-pdb oriented.pdb --out-html view.html --out-viz viz.json --out-pymol view.pml
```

## What it produces

- **Orientation** — membrane normal, bilayer half-thickness, oriented coordinates (+Z out).
- **Per-residue labels** — zone (extracellular / interface / core / periplasmic), facing
  (lipid vs lumen), accessibility (antibody-accessible, LPS-shielded, lipid-embedded, …).
- **Geometric outside-facing surface set** — residues nominated for downstream
  review. This is not native exposure or antibody-accessibility evidence.
- **Rotation-invariance report** — Jaccard stability of the extracellular set under random
  input rotations (the honest correctness check).
- **3D exports** — 3Dmol.js HTML, a JSON descriptor, and a PyMOL script.

## Contexts

| name | membrane | orientation method |
|---|---|---|
| `gram_negative_om` *(default)* | asymmetric LPS | barrel normal |
| `eukaryotic_pm` | symmetric | TM helix-axis v2, experimental; checksum-bound spans required |
| `tm_receptor` | symmetric | TM helix-axis v2, experimental; sidedness may remain unresolved |
| `gram_positive_surface` | none | anchor-relative |
| `soluble_secreted` | none | SASA only |

## Worked examples

```bash
python examples/worked_examples.py     # OmpA β-barrel (1BXW) + glycophorin A TM helix (1AFO)
python examples/p4_benchmark.py        # correctness vs OPM: mean normal error 7.4°
```

Outputs land in `examples/out/` (oriented PDB, labels TSV, 3Dmol HTML, PyMOL script per protein).

## Correctness

In the historical checksum-locked v1 collection, five beta-barrel cases recover
the OPM/PPM reference orientation with a mean normal error of about 7.44 degrees.
Mark 1 therefore reports the beta-barrel scope as conditionally qualified while
independent reproduction remains pending. The v1 alpha-helical whole-structure
method failed its frozen normal-error gate; it has been replaced by the separate
experimental helix-axis path, whose Qualification v2 development and held-out
panel is not yet source-adopted or executed. OPM/PPM is an external computational
reference system, not direct experimental proof of membrane orientation.

## Documentation

`SPEC.md` is the module-by-module contract — signatures, guarantees, and design invariants —
suitable as an implementer hand-off.

## License / status

v0.2.0. Pure Python (NumPy/SciPy/BioPython); SASA computed in-package
(Shrake–Rupley). The clean reviewer selection currently reports 85 passing
offline tests and 4 network/adapter tests deselected.
