# Public CA II structure coverage comparison

This report describes two hash-locked public coordinate files using StructQC. It does not test ligand effects, activity, affinity, native conformation, or which model is better.

## Results

| PDB entry | Context | Resolution (Å) | Resolved residues | Identity | Reference coverage | Missing P00918 positions | Median per-residue mean B factor |
|---|---|---:|---:|---:|---:|---|---:|
| [2POW](https://www.rcsb.org/structure/2POW) | I7C inhibitor; zinc is also present | 1.75 | 257 | 1.000 | 0.988 | 1, 2, 3 | 14.470 |
| [4RN4](https://www.rcsb.org/structure/4RN4) | acetazolamide derivative 3T7; zinc is also present | 1.53 | 258 | 1.000 | 0.992 | 1, 2 | 14.197 |

## Interpretation limits

Coverage and identity are computed against UniProt P00918 by the StructQC sequence mapping. B factors are shown only as deposited coordinate descriptors; different resolution, crystal environment, refinement, and ligand context prevent treating their difference as a controlled quality or ligand effect. No external geometry validation report or PAE was supplied, so these remain missing or unevaluated in the module evidence.

The exact input and software identities are in `INPUTS_CHECKED.json`, `RUN_ENVIRONMENT.json`, and each `structqc/<PDB>/RUN_MANIFEST.json`. The complete machine-readable comparison is `CASE_RESULTS.json`.

No participant, private use case, or independent researcher was involved. This is an example workflow, not evidence of independent use.
