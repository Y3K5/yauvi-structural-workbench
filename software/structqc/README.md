# StructQC

StructQC binds a protein structure to an explicit identity, provenance class,
coordinate checksum, model, chains, and residue map before another module is
allowed to interpret it. It accepts PDB and mmCIF, an optional reference FASTA,
declared provenance JSON, and AlphaFold-style PAE JSON.

It never certifies experimental provenance from a filename or header. Missing
provenance is `unknown`; a recognizable predictor header is recorded only as a
warning. Experimental B-factors and predicted-model pLDDT are kept distinct.

```bash
structqc validate --structure model.pdb --provenance provenance.json
structqc run --structure model.pdb --provenance provenance.json --out out
structqc describe
structqc fetch --plan
```

Outputs are deterministic and contain no absolute paths:

- `STRUCTURE_EVIDENCE.json`
- `RESIDUE_QUALITY.tsv`
- `STRUCTURE_LAYER.json`
- `RUN_MANIFEST.json`

The module reports coordinate evidence. It does not establish native structure,
function, biological assembly, or experimental validity.
