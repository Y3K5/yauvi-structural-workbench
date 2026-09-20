# Sequence coverage in two human carbonic anhydrase II structures

**Case:** `public-ca2-structure-coverage-v1`  
**Analysis run:** 2026-09-20  
**Software:** `yauvi-structural-workbench==0.1.0.dev0` metadata was installed in the environment; the run imported the current source-tree StructQC module, schema 1.1, Biopython 1.88, NumPy 2.5.3, Python 3.12.0  
**Protocol SHA-256:** `39f9737f5615f9431183f0e9db1971dcdc78f8a838064c7d59a579e4bcba941d`  
**Source-lock SHA-256:** `0f341f810472544e8e205f8f9e190de2578aa0941ffa3313556a526fe3dd8a30`  
**StructQC core SHA-256:** `c15cb56344274806844919ba7be02b74a6ee5297741287871eb3bc9754d41ad4`  
**Analysis driver SHA-256:** `b8ce4d2a259e827f37d99225bc6b650d7c99a42fd19ec910432ba3e3c8400225`  
**Acquisition script SHA-256:** `a85a7ec8c3add651258b20e86af2c493d55ea98f6416aaee047f6c9d998d735c`

## Question

For chain A of two public, wild-type human carbonic anhydrase II (CA II)
crystal structures, what sequence identity and coverage does StructQC report
against UniProt P00918, and which canonical positions are absent from the
coordinate chain?

The fixed records are RCSB PDB [2POW](https://www.rcsb.org/structure/2POW),
with inhibitor I7C at 1.75 Å, and [4RN4](https://www.rcsb.org/structure/4RN4),
with acetazolamide derivative 3T7 at 1.53 Å. Both are X-ray diffraction
structures. These are separate experiments and differ in ligand, resolution,
crystal environment, and refinement. The deposited metadata motivates their
description; all numerical measurements below come from the hash-locked
coordinate bytes listed in `public_research_sources.json`.

## Results

| Entry | Chain | Resolved residues | Identity vs P00918 | Reference coverage | Missing reference positions | Median per-residue mean B factor |
|---|---|---:|---:|---:|---|---:|
| 2POW | A | 257 | 1.000 | 0.988462 | 1–3 | 14.470 |
| 4RN4 | A | 258 | 1.000 | 0.992308 | 1–2 | 14.197 |

StructQC maps each selected chain to the 260-residue P00918 sequence. Both
chains are fully identical at the mapped positions; the coordinate chain in
4RN4 maps one more reference position than the chain in 2POW. This is a
description of these two deposited records, not evidence that ligand identity
caused the difference.

## Evidence and limits

The locked inputs are UniProt P00918 FASTA (SHA-256
`74cc1e0de5c8488d747471c3b3f4d2c219bc0afbd3d324cc65e017a65209bf8b`), 2POW
mmCIF (`378221980ca5c658048f8cf1b88bad27ad0f7ad15673d6796a14a282329e1173`),
and 4RN4 mmCIF
(`5cd0cbc3e95161c58de03e27d8df772a388feaa6af14139128fbc44db4c04f43`). The
case runner verifies byte counts and digests before analysis, selects model 0
and author chain A, and stops if an input, model, or chain is missing or
changed. It writes the StructQC evidence bundles and deterministic
`CASE_RESULTS.json` to a new output path outside this source directory.

Per-residue B factors are shown only as deposited coordinate descriptors, not
as directly comparable quality scores. No wwPDB, MolProbity, or Phenix
validation report was supplied, so external geometry validation is missing;
PAE is unevaluated. This case does not test ligand contacts, binding affinity,
inhibition, catalytic activity, native-state structure, or biological
function. It does not claim that either structure is better. No participants,
private use cases, or independent researchers were involved.

The protocol was fixed before this comparison. An earlier exploratory 1CA2
StructQC smoke check occurred before protocol lock; 1CA2 and its measurements
are excluded. See `PUBLIC_RESEARCH_PROTOCOL.md` for the full protocol and
`PUBLIC_RESEARCH_GUIDE.md` for reproduction instructions.

## Installed-wheel replay

After the source-tree run above, the same hash-locked inputs and protocol were
replayed twice in a fresh isolated Python 3.12.0 environment from the prepared
`yauvi_structural_workbench-0.1.0.dev0` wheel, without `PYTHONPATH`. The wheel
SHA-256 was
`b89402fa9a489455545f22bcf1a7c183b8a0285f23dc38bc48370fd1ee3e3007`; runtime
metadata confirmed `module_origin_state=installed_distribution`. All 12
generated case-output files matched byte-for-byte between the two runs. The
portable replay records are in
[`PUBLIC_CASE_REPORT.md`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_REPORT.md),
[`PUBLIC_CASE_CASE_RESULTS.json`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_CASE_RESULTS.json),
[`PUBLIC_CASE_RUN_ENVIRONMENT.json`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_RUN_ENVIRONMENT.json),
[`PUBLIC_CASE_INPUTS_CHECKED.json`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_INPUTS_CHECKED.json),
[`PUBLIC_CASE_SAFETY_TESTS.json`](../../evidence/preparation-2026-09-20/PUBLIC_CASE_SAFETY_TESTS.json),
and [`CANDIDATE_VERIFICATION.json`](../../evidence/preparation-2026-09-20/CANDIDATE_VERIFICATION.json).
This replay verifies execution from an installed wheel; no independent human
research use is claimed.
