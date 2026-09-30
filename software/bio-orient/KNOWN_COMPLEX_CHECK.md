# AQP1 tetramer: known-complex check

Checked locally on 2026-09-30 using the public bovine aquaporin-1 structure
[PDB 1J4N](https://www.rcsb.org/structure/1J4N), author-assigned biological assembly
**1** and deposited model **1**. Assembly 2 was not substituted. The primary
[structure paper](https://doi.org/10.1038/414872a) describes extracellular and
cytoplasmic vestibules connected through the water-selectivity pore.

The selected assembly contains four protein copies, twelve nonyl
beta-D-glucopyranoside (BNG) detergent molecules and 456 waters: **8,116 atoms**.
BNG is a sugar-containing detergent; its presence does not establish native
glycosylation or a membrane lipid. See [CCD BNG](https://www.rcsb.org/ligand/BNG).

## Inputs and exact scope

- Coordinates: `1J4N.cif`, SHA-256
  `dedbebe3321bf70ed404115b47ff5b78c13e6a3df65b275c5d49acb5842c6ae9`.
- Topology: [UniProt P47865](https://www.uniprot.org/uniprotkb/P47865/entry),
  JSON SHA-256 `3ba996f827e6077bc1f5182c21c50f6a997fb524c4199394d356a9180ebe0c99`.
  All 249 observed protein residues per copy matched the accession at their
  deposited author numbers, with zero substitutions. The 271-residue reference
  contains unobserved C-terminal positions; no completeness claim is made.
- Six mapped transmembrane spans supported the existing MemOrient placement.
  These UniProt spans carry similarity evidence ECO:0000250 from P29972;
  topological domains carry curated ECO:0000305. They are not newly measured
  experimental annotations. Residues 45 and 163 supplied cited compartment markers,
  explicitly remapped to copies `A:1` through `A:4`.
- Twenty-two local CCD definitions were checksum locked. Manifest SHA-256:
  `6cdf8112a60de7ca5fe6f6f38a89fdbb7c8e112ee0e8b161e2740d3983eb74be`.
- Experimental provenance and the producing wwPDB validation XML were retained.
  Imported validation includes 2.2 Å resolution, R-free 0.308 and clashscore 33.93;
  those are imported measurements, not recalculated quality findings.

## Graph findings

The final workbench run completed in approximately **24 seconds** with nine
declared/evidence-mapped patches, 473 represented/context objects and **4,256
relationships**. Run identity: `run-2830abf97adc7b0b`. Selection identity and
coordinate hash remain visible in both viewers.

| Selection | Recorded interpretation |
|---|---|
| Extracellular loop A, residues 38–46, copies 1 and 2 | All selected atoms lie on the supported extracellular face; sampled outward direction is unresolved |
| Cytoplasmic loop D, residues 156–165, copies 1 and 2 | Curated cytoplasmic role retained; whole patch intersects the modeled core, so orientational availability stays unresolved |
| Cytoplasmic tip, residue 164, copies 1 and 2 | All atoms lie on the supported cytoplasmic face; no coherent single outward normal |
| His182 selectivity site, copy 1 | Source-annotated selectivity role retained; core location and zero sampled exposure do not establish current transport activity |
| 73-residue measured tetramer interface, copy 1 | Measured contacts and reduced exposure; no coherent outward normal and functional availability unknown |
| BNG801, copy `B:1` | Detergent component selected separately; all five evaluated sugar centers match the CCD reference |

The membrane normal is approximately `(0.01142, 0.00856, -0.99990)`, with center
`(34.10989, 31.15686, 25.24883)` Å and a modeled half-thickness of 14 Å. Its
positive face is mapped to **cytoplasm** and negative face to **extracellular**.
The labels come from mapped evidence, not from the deposited +Z axis.

For the interface patch, removing neighbor copy `A:3` recovered **1,178.32 Å²**
of sampled SASA, with **675** heavy-atom pairs within 5 Å and a minimum distance
of **2.51 Å**. This is measured occlusion in this represented assembly, not physical
contact area or a functional-availability score. Every side and edge keeps
functional availability unknown; static coordinates generated no probabilities.

Multisided, homosided and heterosided descriptors coexist. Equivalence groups
refer to explicitly declared homologous patch positions across copies; they
do not establish equivalent functional activity.

## Direction and chemistry controls

All **1,108 evaluated centers** matched their reference geometries, including
**60 BNG sugar centers** across twelve molecules. Proper rotation/translation
preserved all findings. An artificial in-memory mirror produced 1,108 reference
mismatches, and the biology-preserving rotation path rejected a reflection.
The perturbed controls were never attached as biological evidence to the real view.

Chemistry coverage remains **partial**: 100 glycine component rows did not pass
whole-component reference coverage, while 456 waters have no declared tetrahedral
centers. No absolute configurations were invented for these components. Matching
evaluated centers is not a blanket chemical-validation claim.

Excluding the two marker positions, 432 topology-annotated C-alpha positions were
checked across the four copies. **308** lay outside the modeled slab and all
agreed with their recorded compartment labels; **124** lay in the modeled core
and stayed unassigned. This is a geometry sanity check against curated annotations,
not an independent biological qualification benchmark.

The installed CLI reproduced the graph byte for byte, SHA-256:
`3256211e505506b9a5f70c6bb95070c45d486fa72a95e24237a5a57c2ab11942`.

## Inspect locally

The loopback test workspace is separate from the software distribution. Its
browser entry is `AQP1 1J4N · membrane sidedness and chemistry`. In Findings,
select `extracellular-loop-1`, `cytoplasmic-tip-1` or `tetramer-interface`, then
inspect a neighbor. Select `BNG · B:1/801/BNG` and atom `C1` to inspect sugar
chirality. Evidence supplies the exact source and full ZIP. The connected viewer
uses the same graph. Reset, zoom, focus and light/dark background controls affect
presentation only; labeled outlined planes show their evidence state.

See [delivery checks](VERIFICATION.md) for environment and broader test limits.
