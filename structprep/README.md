# structprep

Deterministic structural editing with a derivation record.

Removes solvent and crystallisation additives, selects chains and alternate
locations, and flattens biological assemblies — then writes down exactly what it
removed and what it refused to.

## Why it derives instead of edits

`structqc` binds a SHA-256 to a set of coordinates, and every downstream module
refuses input that does not match. A tool that edits coordinates in place would
break that. So `structprep` never edits: it writes new coordinates plus a record
naming the parent checksum and every atom removed. Re-run `structqc` on the
child and the chain survives:

    structqc(raw) → structprep → structqc(prepared) → memorient / site-context / …
                        ↓
                 PREPARATION.json    parent sha256, policy, per-class counts

## Use

```bash
structprep classes                       # the component table, versioned
structprep validate --structure in.pdb   # what would happen, nothing written
structprep run --structure in.pdb --out prepared/
```

Defaults remove **solvent, cryoprotectant and buffer**. Ions and unclassified
components are reported, never dropped without being named: a magnesium may be a
crystallisation additive or the catalytic centre, and the residue code cannot
tell you which.

```bash
structprep run --structure in.pdb --out p/ --keep ZN --drop GOL --chains A
```

## What it refuses

- **Removing a classed component that touches the polymer** within 4 Å. Glycerol
  sitting in a pocket may be occupying a real site. Name it in `--drop` to
  override. Solvent is exempt — nearly every crystal water contacts protein, so
  the rule would refuse everything and mean nothing.
- **Silent alternate-location choice.** Keeps A by default and records that B
  existed.
- **Silent interface changes.** A water bridging two chains that gets removed
  earns a warning, because the interface geometry changes with it.

`run` exits **1** when anything was refused, so a caller who asked for a removal
that did not happen learns it from the exit code.

## What it does not do

No protonation, tautomer assignment, charge assignment or minimisation. Those are
force-field- and pH-dependent, and bundling them here would let a model-dependent
assumption travel under the authority of a deterministic step. **The output is
not docking-ready** until a step that declares its force field and pH has run.

Removing a component does not establish that it was not functional.

## Tests

    PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest tests/ -q

Two of them encode mistakes that reached working code before this module existed:
a haem nearly stripped with the waters, and a biological assembly whose duplicate
chain ids made a tetramer look like a dimer — 98 Å² of buried surface instead of
952.
