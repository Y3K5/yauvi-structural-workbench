# Where the raw files come from

Every input this module reads, with its licence and the exact way to obtain it.
The machine-readable version is `src/actstate/sources.yaml`, resolved against the
workspace registry `catalogs/sources.yaml`.

    yauvi-fetch plan --for actstate     # what is needed, and what is already here
    yauvi-fetch get  --for actstate     # retrieve what the licence permits

## Required

### UniProtKB annotation export

The annotation table, carrying the sequence features the classifier reads.

- **Registry id:** `uniprot_proteomes`
- **Licence:** CC BY 4.0
- **Access:** REST API, no registration
- **Retrieved by:** `yauvi-fetch get --for actstate --arg uniprot_proteomes=<UP-accession>`

Manual equivalent — note the four feature fields, which most existing exports omit:

```
https://rest.uniprot.org/uniprotkb/stream
  ?query=proteome:UP000005640
  &format=tsv
  &fields=accession,protein_name,gene_names,protein_existence,cc_function,
          cc_subcellular_location,go_id,keyword,xref_pfam,xref_interpro,
          lit_pubmed_id,ec,ft_act_site,ft_binding,ft_site,cc_cofactor,
          cc_activity_regulation
```

**Versioning caveat.** UniProt records no release identifier inside the export.
`yauvi-fetch` captures the `X-UniProt-Release` response header at retrieval time
and stores it in the cache manifest; a file obtained any other way is pinned only
by its digest and its retrieval date.

## Optional

### AlphaFold Protein Structure Database

Predicted monomer models, keyed by UniProt accession.

- **Registry id:** `alphafold_db`
- **Licence:** CC BY 4.0
- **Access:** public API — resolve `pdbUrl` from
  `https://alphafold.ebi.ac.uk/api/prediction/<accession>` rather than guessing a
  `model_v4` path; model versions advance
- **Effect if absent:** `geometry` and `occupancy` report `unevaluated`
- **Effect if present:** enables geometry, but the label is capped at
  `probable_active` — a prediction cannot establish a functional state

### Protein Data Bank

Experimental coordinates.

- **Registry id:** `pdb`
- **Licence:** CC0 1.0 (public domain)
- **Access:** `https://files.rcsb.org/download/<id>.cif`
- **Effect if present:** the only source that can support `active_state_supported`,
  and the only real evidence of cofactor occupancy

### InterPro / Pfam

Family and domain assignment, for reading a site in its mechanism class.

- **Registry id:** `interpro`
- **Licence:** CC0 1.0 for the data; InterProScan is Apache-2.0
- **Access:** API, or InterProScan locally
- **Effect if absent:** nothing is blocked; family context is not used to reach a label

## Not required, and not shipped

### Reference structures of known state

The `conformation` signal compares coordinates to references whose active or
inactive state is established. **No curated set ships with this module**, and none
is declared in the registry, so the signal reports `unavailable` on every run
until one is supplied via `--reference-comparison`.

This is stated here rather than left to be discovered because it is the single
largest limitation of the module: without it, no protein can reach
`active_state_supported` on the strength of conformation alone.

Building one means choosing, per mechanism family, PDB entries whose state is
documented in the literature — a curation task, not a download.
