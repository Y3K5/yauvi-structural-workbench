# Where the raw files come from

The machine-readable version is `src/sf_csa/sources.yaml`, resolved against the
workspace registry `catalogs/sources.yaml`.

    yauvi-fetch plan --for sf_csa
    yauvi-fetch get  --for sf_csa

## Required

### Protein Data Bank — the structural leg

- **Registry id:** `pdb` · **Licence:** CC0 1.0 (public domain)
- **Access:** `https://files.rcsb.org/download/<id>.cif`, or a Foldseek-formatted
  PDB100 database
- **Pinned:** the database manifest records the database checksum and its version
  string. Two releases built against different snapshots are **not comparable**,
  and `build-manifests` refuses a database with no recorded version.

Foldseek's own prebuilt databases: `foldseek databases PDB <out> tmp`. Record the
resulting `pdb.version` file next to the database — the builder reads it.

### UniProt reference proteomes — the sequence leg

- **Registry id:** `uniprot_proteomes` · **Licence:** CC BY 4.0
- **Access:** REST API, no registration
- **Retrieved by:** `yauvi-fetch get --for sf_csa --arg uniprot_proteomes=<UP-accession>`

UniProt records no release identifier inside the export; `yauvi-fetch` captures
the `X-UniProt-Release` header at retrieval time.

### DIAMOND — a runtime, not a file

- **Registry id:** `diamond` · **Licence:** GPL-3.0
- Install and put on PATH. Resolved by `shared/runtime-registry.yaml`, which is
  fail-closed: absent means the module blocks, never that the sequence leg is
  quietly skipped.

## Optional

### AlphaFold DB

- **Registry id:** `alphafold_db` · **Licence:** CC BY 4.0
- Predicted models for queries with no local structure. Absent models are
  reported as `candidate_missing_structure` — **never** as structural negatives.

### InterPro / Pfam

- **Registry id:** `interpro` · **Licence:** CC0 1.0 (data)
- Family assignment, for reading a hit in its mechanism class.

## Not a source: Foldseek

Foldseek is the structure aligner and is declared in
`shared/runtime-registry.yaml`, not in `catalogs/sources.yaml` — it is an
executable, not a data file. The database manifest pins the **required version**
(`required_foldseek_version`), and `core.run_pipeline` raises on a mismatch: an
aligner version change is a change to the results.

## What you must supply yourself

The **campaign spec** — the targets, their structures, and the ledgers that pin
their sequences. That is your experiment, not a download. `examples/` in this
package shows the shape; there is no default campaign, and the package ships no
target dictionaries.
