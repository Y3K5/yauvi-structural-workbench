# Public structural qualification v1

This directory is an external scientific qualification panel, not a unit-test
fixture. It asks whether the current YAUVI outputs agree with independent public
records under predeclared gates. A passed case does not by itself qualify a
workflow across protein classes.

## Reproduce locally

The runner performs no network access. Public artifacts must already match
[`SOURCE_LOCK.json`](SOURCE_LOCK.json).

```bash
python yauvi-structural-workbench/benchmarks/qualification-v1/run_qualification.py
```

Exit code `0` means every required check for all six workflows passed. Exit code
`1` means one or more scientific checks are partial or failed. Exit code `2`
means source drift blocked the run.

Primary outputs:

- [`results/QUALIFICATION_RESULTS.json`](results/QUALIFICATION_RESULTS.json)
- [`results/SOURCE_VERIFICATION.json`](results/SOURCE_VERIFICATION.json)
- [`QUALIFICATION_REPORT.md`](QUALIFICATION_REPORT.md)

`build/` contains reproducible derived inputs and raw workflow outputs. It can be
deleted and recreated by the runner. It is not the source lock.

## Independent authorities

| Evidence leg | Independent authority | How it is used |
|---|---|---|
| Coordinates, assemblies, validation | [wwPDB/RCSB file services](https://www.rcsb.org/docs/programmatic-access/file-download-services) | Exact mmCIF, biological assembly, and wwPDB validation XML artifacts |
| Predicted confidence | [AlphaFold Protein Structure Database](https://www.alphafold.ebi.ac.uk/) | Versioned model, declared pLDDT encoding, and PAE matrix |
| Membrane placement | [OPM/PPM](https://pmc.ncbi.nlm.nih.gov/articles/PMC3245162/) | OPM-oriented coordinates define the Z membrane normal and report half-thickness |
| Kinase state labels | [KinCore ABL1](https://dunbrack.fccc.edu/kincore/GENE/ABL1) | Experimental-chain active/inactive labels independent of StateAtlas RMSD |
| Catalytic residues | [M-CSA entry 1](https://www.ebi.ac.uk/thornton-srv/m-csa/entry/1/) | Glutamate-racemase residue identities, positions, and chemical-role text |
| SASA | [FreeSASA](https://freesasa.github.io/) | Lee-Richards solvent-accessible surface area through a pinned local CLI |
| Structural relationships | [CATH downloads](https://cathdb.info/download) | Frozen v4.3.0 exact, homologous-superfamily, topology-analogy, and unrelated controls |
| Structure search | [Foldseek](https://github.com/steineggerlab/foldseek) | Real local structure searches against the CATH-labeled mini-database |
| Sequence search | [DIAMOND](https://www.nature.com/articles/s41592-021-01101-x) | Real local protein-sequence searches kept separate from Foldseek results |

## Frozen case gates

- StructQC: exact author/label identity, exact reference mapping, correct raw
  wwPDB geometry values, prediction-only pLDDT interpretation, valid PAE, and
  fail-closed unknown provenance.
- MembraneOrient: mean normal error at most 15 degrees, mean half-thickness
  error at most 2.5 A, and every repeated-orientation residue-set Jaccard at
  least 0.95, evaluated separately for beta-barrel and alpha-helical strata.
- StateAtlas: experimental active and inactive ABL references, maximum RMSD
  2.5 A, minimum margin 0.25 A, and no confident opposite-state held-out call.
- SiteContext: the six M-CSA residues for 1B73 must map exactly; absent pocket
  evidence remains incomplete and no catalytic-activity claim is allowed.
- AssemblyContext: exact 4HHB tetramer stoichiometry, 5 A contacts, positive
  burial, Gemmi assembly handling, and real FreeSASA invocation.
- SF-CSA: real Foldseek and DIAMOND executions, exact and CATH-homolog controls,
  no promotion of analogy/unrelated controls, separate evidence tables, and a
  checksum spanning every Foldseek database sidecar.

## Current interpretation

Four cases pass their predeclared gates. Membrane orientation is partial because
the alpha-helical stratum fails normal-error and rotation-invariance gates.
Conformational state is partial because held-out 8SSN chain A remains unresolved
under the frozen two-reference full-domain CA comparison. These are release
blockers, not favorable or silently omitted results.

The current qualification does not establish clinical utility, biochemical
activity, native membrane exposure, binding affinity, catalytic activity, exact
functional transfer, or general accuracy outside the named public cases.

## Redistribution boundary

The artifacts were acquired for local qualification. Before a public repository
or release includes source files, each provider's current redistribution and
attribution terms must be reviewed. If redistribution is not approved, publish
only the lock/acquisition manifest and require reviewers to acquire the exact
checksummed artifacts themselves.
