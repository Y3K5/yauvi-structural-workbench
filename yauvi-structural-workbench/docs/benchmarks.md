# Benchmarks and independent qualification

Software tests and scientific qualification are independent gates.

| Workflow | First public case set (qualification v1) | Result in that historical case set |
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
[`qualification-v1/`](../../evidence/benchmarks/qualification-v1/README.md).
Those named cases are retained; their thresholds are not rewritten.

The retained Qualification v2 source manifest is collection 2.11. It adopts 94
of 110 required records across six panels: ABL StateAtlas (14), StructQC (16),
functional-site context (16), AssemblyContext (16), SF-CSA (16), and beta-barrel
MembraneOrient (16). The alpha-helical membrane stratum remains unadopted and
non-blocking.

Public workflow run 35295451246 on commit
`7981148e70c5eaa6424e3608f09cd65bcfb35834` recorded passing results for the five
release-blocking panels on seven CI runners. This historical evidence belongs to
that public commit. The checked-in aggregate `EXECUTION_SUMMARY.json` is older:
it totals 67/110 across five executed panels and predates the later SF-CSA
result. It is neither the collection's adopted-record count nor a result for
the changed local candidate. Consult the [JOSS research-use and human-review
handoff](JOSS_RESEARCH_HANDOFF.md) for evidence boundaries and remaining work.

The collection and CI outcomes describe curated cases and recorded software
behavior only. They do not establish biological validation or qualify source
changes made after the cited public commit. Alpha-helical membrane orientation
remains experimental; beta-barrel orientation failed its objective accuracy
gate and carries no Mark 1 accuracy claim.
