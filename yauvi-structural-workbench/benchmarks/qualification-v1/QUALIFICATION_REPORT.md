# YAUVI public qualification report v1

Overall state: **incomplete or failed gate**. Four public cases pass and two are
partial. Therefore the six-workflow release gate is **not satisfied**.

| Workflow | Case result | Independent result |
|---|---|---|
| StructQC | Passed | 1CRN author/label mapping and wwPDB metrics imported correctly; AlphaFold P69905 pLDDT/PAE handled only under declared prediction provenance; unknown provenance returned incomplete |
| MembraneOrient | Partial | Five beta barrels passed with 7.441992 degree mean normal error and 0.743508 A mean half-thickness error; the alpha-helical stratum failed with 31.60872 degree mean normal error and a 1U19 rotation Jaccard of 0.88 |
| StateAtlas | Partial | Held-out 2V7A chain A matched the KinCore active label; held-out inactive 8SSN chain A remained unresolved at 6.23899 A from the closest inactive reference; no confident opposite-state call occurred |
| SiteContext | Passed | M-CSA positions 7, 8, 70, 147, 178, and 180 mapped exactly on 1B73; missing pocket evidence remained explicit |
| AssemblyContext | Passed | 4HHB assembly 1 retained A/B/C/D stoichiometry, 201 residue-pair contacts, and 1760.002685 A2 subject-chain burial by FreeSASA (version not captured by the runner) |
| SF-CSA | Passed | Foldseek 10.941cd33 and DIAMOND 2.1.11 ran on four CATH-labeled controls; exact and homolog controls were recovered, analogy/unrelated controls were not promoted, and 14 Foldseek database files were checksum-bound |

## Defects found and fixed during independent execution

1. StructQC confused a wwPDB clashscore percentile with the raw clashscore and
   missed hyphenated Ramachandran/rotamer fields. The parser now normalizes names,
   prefers the shortest exact metric suffix, and has a realistic XML regression
   test.
2. AssemblyContext used obsolete FreeSASA arguments and expected a different JSON
   key. It now uses `--depth=chain`, accepts FreeSASA 2.1.x `structure` output,
   and uses mmCIF for expanded multi-character chain identifiers.
3. SF-CSA checksummed only the main Foldseek prefix. It now records and verifies
   a composite checksum over all 14 prefix files and fails closed on sidecar
   drift.

Focused regression after these changes: **123 tests passed** across StructQC,
AssemblyContext, and SF-CSA. The complete recorded offline workbench selection
then passed **508 reviewer-scope tests with 5 network/adapter tests deselected
and 1 optional MDAnalysis test skipped**. This software
regression result is separate from the six scientific case gates.

The complete qualification runner was executed twice from the frozen source
lock. `QUALIFICATION_RESULTS.json` and `SOURCE_VERIFICATION.json` were
byte-identical across the two runs.

## Remaining scientific work

- Restrict MembraneOrient's qualified scope to the passing beta-barrel class or
  improve and rebenchmark the alpha-helical placement method without changing
  the frozen tolerances after seeing results.
- Replace whole-domain RMSD as the only StateAtlas discriminator with a curated
  kinase-state selection or state-defining collective variables, then rerun a
  larger KinCore-held-out panel. Until then, unresolved is the correct safe call.
- Expand each passing case to a stratified multi-protein benchmark and complete
  redistribution/license review before claiming workflow-level qualification.
- Rerun the complete macOS/Linux Python 3.10-3.12 matrix. This local run used the
  current Python 3.12 environment only.

## Evidence identity

- Qualification JSON SHA-256:
  `06c7e84e10e45335d506554d0587ea4979ca8c0b113c4c4e2d0a9c015d3f8eed`
- Source-verification JSON SHA-256:
  `48091cfb8fc1385d751d3226c503f96b7312f3cfcfbbae596724629dd561f3d1`

This report does not claim biochemical activity, clinical usefulness, native
surface exposure, affinity, catalytic activity, therapeutic effect, or exact
functional transfer.
