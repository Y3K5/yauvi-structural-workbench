# Public-data case: sequence coverage in two human carbonic anhydrase II structures

Protocol locked 2026-09-20 before the primary comparison was run. This is a
small, descriptive, reproducible software case study using public experimental
structures; it is not a biological validation or an independent-use claim.

## Question and scope

For chain A of two public human carbonic anhydrase II (CA II) crystal
structures, what sequence identity and canonical-sequence coverage does
StructQC report against UniProt P00918, and which canonical positions are
unrepresented in the coordinate chain? The fixed records are RCSB PDB 2POW
(human CA II with ligand I7C) and 4RN4 (human CA II with inhibitor 3T7).
RCSB describes both as wild type human CA II and reports chain A; the entries
are separate experiments with different ligands and crystallographic
conditions, not a controlled pair.

The predeclared measurements are: coordinate-residue count, reference length,
mapped-residue count, identity fraction, coverage fraction, missing reference
positions, and the minimum/median/maximum of StructQC's per-residue mean
B-factor values. B factors are reported as descriptive coordinate fields only;
they are not treated as a comparable model-quality score across separate
experiments. Missing external validation, PAE, or other evidence remains
missing/unevaluated.

## Sources and fixed selections

- Coordinates: `https://files.rcsb.org/download/2POW.cif` and
  `https://files.rcsb.org/download/4RN4.cif`; model index 0, author chain A.
- Reference: `https://rest.uniprot.org/uniprotkb/P00918.fasta` (canonical
  human CA II sequence, accession P00918).
- Entry metadata: [2POW](https://www.rcsb.org/structure/2POW) and
  [4RN4](https://www.rcsb.org/structure/4RN4). Metadata establishes record
  identity and experimental method only; the analysis uses the hash-locked
  coordinate bytes.
- The selected PDB identifiers were checked against the adopted Qualification
  v2 panel manifest and are not panel cases.

After acquisition and before analysis, `public_research_sources.json` records
the SHA-256 digest, byte count, URL, and selected chain for every input. A
digest mismatch, missing selected chain, invalid structure, or missing FASTA
stops the run. The workflow does not substitute another record or silently
drop an input.

## Analysis and interpretation

Run the supplied acquisition script to fetch the three public inputs, then run
`public_research_case.py` with those inputs and a new output directory outside
this source directory. The runner invokes the installed `structqc` command
once per entry, pins model 0 and chain A, and records tool/runtime versions,
input digests, command parameters, raw module outputs, and a deterministic
JSON/Markdown comparison. A run directory must not already exist.

Results describe what these deposited coordinate files encode relative to the
chosen UniProt sequence. They do not establish native conformation, activity,
ligand affinity, causal effects of ligand identity, or which structure is
“better.” The two structures differ in more than ligand identity. StructQC does
not calculate crystallographic validation metrics; without an imported
wwPDB/MolProbity/Phenix report, those metrics remain unavailable.

No biological or human-participant research is performed. The data are public
macromolecular structures and a public reference sequence. The output is an
example workflow, not evidence of independent research use.

## Development disclosure

An earlier exploratory StructQC smoke check used 1CA2 before this protocol was
locked. That record and its output are excluded from this case, its source lock,
and its reported measurements. The primary comparison is restricted to the
two records fixed above.
