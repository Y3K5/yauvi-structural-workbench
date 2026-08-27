# Files and official sources

Every input card includes its purpose, accepted extension, minimum format,
missing-evidence effect, official sources, and an optional checked template.

| Artifact | Provider | Identifier | Workbench behavior |
|---|---|---|---|
| Experimental coordinates | RCSB PDB | PDB ID | Verified mmCIF import when enabled |
| Biological assembly | RCSB PDB | `PDB:assembly` | Verified assembly mmCIF import when enabled |
| Validation XML | wwPDB/RCSB | PDB ID | Bounded download, safe decompression, checksum import |
| Predicted model and PAE | AlphaFold DB | UniProt accession | API-resolved model/PAE import when enabled |
| Sequence and annotations | UniProtKB | UniProt accession | FASTA or feature TSV import when enabled |
| Proteome FASTA | UniProt | UP identifier | Frozen local comparison universe |
| Residue mappings | SIFTS | PDB ID | Official link and manual normalized import |
| Catalytic sites | M-CSA | Entry/PDB ID | Official link and checked annotation template |
| Components/cofactors | PDB CCD and ChEBI | Exact identifier | Official link and checked mapping template |
| Membrane benchmark | OPM/PPM | PDB ID | Curator-frozen external comparison pack |

Default mode performs zero acquisition requests. `--allow-reference-fetch`
enables only registered artifact/identifier builders. It does not enable
arbitrary URLs, database-wide browsing, private uploads, or target-derived
queries. Cache refresh and analysis adoption are separate actions.

Trajectories, pocket-tool results, SF-CSA interpretation tables, and frozen
Foldseek databases must be generated or staged locally using their displayed
format guide.
