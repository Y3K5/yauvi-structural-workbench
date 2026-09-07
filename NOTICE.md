# Notice and third-party boundary

YAUVI Structural Workbench invokes or links to third-party software and public
data but does not claim ownership of them. FreeSASA, Foldseek, DIAMOND,
MDAnalysis, Gemmi, Biopython, RCSB PDB, wwPDB validation, AlphaFold DB, UniProt,
SIFTS, M-CSA, PDB CCD, ChEBI, CATH, and OPM/PPM retain their own licenses, terms,
and citation requirements.

The public release gate must verify every bundled file. External binaries,
database snapshots, downloaded structures, validation reports, and model weights
are excluded unless a reviewed redistribution record explicitly permits them.

## Reference-data licences, as reviewed on 7 September 2026

Every provider recorded in a `SOURCE_LOCK.json` was checked against its own
published terms. All are public-domain or attribution licences; none restricts
redistribution, and none is copyleft. This repository redistributes none of the
data regardless — the source lock ships so a reviewer re-acquires the identical
checksummed bytes themselves.

| Provider | Licence |
|---|---|
| wwPDB / RCSB PDB, incl. validation reports and the PDB CCD | CC0 1.0 public domain dedication |
| UniProt / UniProtKB | CC BY 4.0 |
| AlphaFold DB | CC BY 4.0, plus EMBL-EBI Terms of Use |
| M-CSA | CC BY 4.0 |
| CATH | CC BY 4.0, attribution to Sillitoe, Dawson, Lewis, Lee, Lees and Orengo |
| SIFTS (PDBe-KB) | CC BY 4.0 |
| ChEBI | CC BY 4.0, plus EMBL-EBI Terms of Use |
| OPM / PPM | reported as CC BY 3.0 by secondary sources; **not stated on the OPM site itself**, so treat the attribution requirement as binding and the version as unconfirmed |

The CC BY entries carry an attribution obligation on this project, not only on
its users: where the benchmark collections and showcases present data derived
from these resources, the provider is named in the record that presents it.
