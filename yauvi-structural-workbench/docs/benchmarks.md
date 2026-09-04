# Benchmarks and independent qualification

Software tests and scientific qualification are independent gates.

| Workflow | Public qualification collection | Current result |
|---|---|---|
| StructQC | 1CRN wwPDB validation, AlphaFold P69905 v6 model/PAE, and unknown-provenance control | Public case passed; workflow-general qualification not established |
| MembraneOrient | Five beta-barrel and three alpha-helical OPM structures with rotation invariance | Partial: beta-barrel stratum passed; alpha-helical normal error and 1U19 rotation invariance failed |
| StateAtlas | KinCore-labeled two-sided ABL references with active and inactive holdouts | Partial: active holdout passed; inactive holdout remained unresolved; no opposite-state false call |
| Functional site | M-CSA glutamate racemase/1B73 | Public case passed; absent pocket evidence remained incomplete |
| AssemblyContext | 4HHB biological assembly, stoichiometry, contacts, and FreeSASA burial | Public case passed with FreeSASA (version not captured) |
| SF-CSA | CATH-labeled exact, homolog, topology-analogy, and unrelated controls | Public mini-case passed with real Foldseek 10.941cd33 and DIAMOND 2.1.11 runs |

Public artifacts are acquired from versioned checksum manifests unless a
reviewed redistribution record permits bundling. Stochastic or platform-specific
external tools compare scientific invariants and tolerances rather than claiming
cross-hardware byte identity.

The historical lock, measurements, and offline reproduction command remain in
[`benchmarks/qualification-v1/`](../benchmarks/qualification-v1/README.md).
Those named cases are retained; their thresholds are not rewritten.

The current release gate is
[`benchmarks/qualification-v2/`](../benchmarks/qualification-v2/README.md).
Qualification v2 freezes scope-specific readiness, the full stratified panel,
development and held-out splits, exact ABL mappings, and unchanged numerical
gates. Its current audit state is `blocked_panel_incomplete`: four of six panels
are adopted and executed, and two -- ABL StateAtlas and SF-CSA, both
release-blocking -- are not. That visible block prevents historical
demonstrations or absent data from becoming a favorable release result, and it
holds regardless of how the executed panels performed.

All five Mark 1 release-blocking scopes must pass v2 and reproduce on a second
machine before `local_release_candidate`. Alpha-helical membrane orientation is
an experimental, non-blocking scope and cannot appear qualified.
