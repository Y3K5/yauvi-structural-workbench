## 2026-09-11 — mmCIF is read and written natively, and a chain label is not an entity label

Capability change in `read_structure` and `write_outputs`, plus one mistake made and
fixed inside this change, recorded because it altered atom records while the
derivation record claimed it had not.

`read_structure` refused mmCIF outright — *"mmCIF is not supported yet; convert to
PDB first"* — which pushed every caller holding a cryo-EM entry into a conversion
step outside the evidence chain. That is the one thing this module exists to
prevent: the converted file is where `structqc` then binds its checksum, so the
parent sha256 names a re-rendered artifact rather than the deposited coordinates.
It also imposed PDB's ceilings on files that do not have them — one character of
chain id, 62 flattening slots, 99,999 atom serials.

The fix keeps the guarantee that makes the PDB path trustworthy rather than
reaching for a parser that would break it. `structprep` has no third-party
dependency and does not format coordinates: it slices the columns it needs out of
each source line and re-emits that line. The mmCIF path does the same thing with
tokens instead of columns — each `atom_site` row is tokenised with the character
offsets of every value, so a row is re-emitted exactly as written and only a
deliberately chosen token can be spliced. A `gemmi` or Biopython reader would have
re-rendered every coordinate and added the first numeric dependency to the one
module that is immune to a wheel or version mismatch.

**The mistake.** Flattening an assembly relabels chains. The first version rewrote
both `_atom_site.auth_asym_id` and `_atom_site.label_asym_id`, on the assumption
that they are two spellings of the chain. They are not: `label_asym_id` is an
entity label, and the two legitimately differ — in PDB 8DBK the ATP of auth chain A
carries `label_asym_id` G. Rewriting both changed 52 atom records in that entry and
reassigned their entity while `PREPARATION.json` reported a chain relabel and
"coordinates are unchanged". Worse, the splice ran on every row whether or not any
flattening had happened, so a single-model file came back modified for no reason.
Only `auth_asym_id` is rewritten now, only when a remap actually assigns a new
label, and `label_asym_id` is left as deposited.

Changes:

- `read_structure(path)` dispatches on content and suffix; `_read_pdb` and
  `_read_cif` return `(atoms, header, layout)` and `read_structure` keeps its
  two-value return for existing callers.
- `_locate_atom_site` finds the one `atom_site` loop and **refuses** a file with
  none or with more than one, rather than picking.
- `_cif_tokens` tokenises a data region into `(value, start, end)`, honouring
  single and double quotes, skipping `#` comments, and refusing a multi-line `;`
  value inside the loop rather than guessing where a row ends.
- `_read_cif` refuses a ragged row count, a missing `Cartn_x/y/z` column, and a
  non-numeric coordinate. `.` and `?` are read as absent. Values come from the
  `auth_*` columns first, so chain selection and residue numbering mean what they
  mean in the PDB path.
- Each row is stored as its whole physical line, so trailing whitespace and column
  alignment survive; a row wrapped across lines falls back to its token span.
- `_chain_ids` replaces the fixed 62-character alphabet: PDB keeps the
  single-character ceiling its format requires, mmCIF continues into two-character
  ids instead of refusing a large assembly.
- `write_outputs` emits `PREPARED.cif` for an mmCIF parent, carrying every category
  before and after the loop through byte for byte, and `PREPARED.pdb` as before.
- `PREPARATION.json` gains `coordinate_format`, a limitation recording that atom
  records are re-emitted rather than re-rendered, and a limitation stating that
  `label_asym_id` is left as deposited. `schema_version` 1.0 → 1.1.

Verified on PDB 8DBK, a 6-chain 2.10 Å cryo-EM hexamer: 14,612 atoms in, 2,472
kept, 12,140 removed, identical to the count the PDB-converted route reaches, and
all 2,472 kept rows plus all 4,855 non-atom lines byte-identical to the source. The
PDB path is unchanged — the same input produces the same `prepared_sha256` as
before this change. `structqc` re-binds `PREPARED.cif` at 316/318 mapped,
`identity_fraction` 1.0, no warnings, so the provenance chain now closes without
leaving mmCIF.

13 self-contained tests in `tests/test_structprep_mmcif.py`. They use synthetic
fixtures and always run, unlike `tests/test_structprep.py`, which skips unless
`STRUCTPREP_TEST_STRUCTURES` points at unshipped RCSB coordinates.
