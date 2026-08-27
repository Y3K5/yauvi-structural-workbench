# SiteContext

SiteContext maps declared functional residues and cofactors onto checksum-bound
coordinates. Catalytic roles, observed ligands, residue geometry, and predicted
pockets remain separate evidence legs.

```bash
site-context run --manifest STRUCTURE_EVIDENCE.json --structure model.pdb \
  --annotations sites.json --pocket-result fpocket.json --out out
```

Annotations use reference-sequence positions. A site with no unique StructQC
mapping remains unresolved. Cofactor identities match exact component or
declared ChEBI mappings; unknown synonyms are never silently accepted.
