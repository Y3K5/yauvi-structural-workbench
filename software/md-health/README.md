# MD Health: GROMACS preparation and early signals

Supplied GROMACS files -> a separately recorded test copy -> explicit preprocessing and a bounded pilot -> inspectable early-failure evidence.

This is a local prototype. It is included in root Workbench package discovery and has its own `md-health` command/report. It is not yet an analysis type in the Workbench browser. No real GROMACS execution has been verified here; current tests use invented data and mocked engine calls.

## What it does

- Checks a GRO file's atom count, finite coordinates/velocities and positive-volume box. The authoritative topology/coordinate correspondence check still belongs to GROMACS.
- Inventories a supplied MDP and the relative topology/include graph. Missing dependencies, escaping paths, duplicate settings and unsupported external include paths are refused. Conditional includes are conservatively inventoried on both branches; no complete topology preprocessor is reimplemented.
- Copies only the declared dependency set into a new directory. Originals remain untouched. A separate `pilot.mdp` changes only the explicitly requested `nsteps`; its derivation and source fingerprints are recorded.
- Offers an explicit single-process GROMACS pilot with strict `grompp -maxwarn 0`, processed topology/resolved MDP outputs, a fresh TPR and bounded `mdrun`. Commands, runtime identity, stdout, stderr, timeouts and exit codes are preserved. No warning override, automatic parameter repair or old TPR reuse.
- Inspects recorded log/XVG copies for fatal errors, constraint warnings, non-finite values, non-increasing times, truncated series, absent required metrics and excursions outside **user-declared** bands.
- Runs an optional recorded `gmx check -f ... -s1 ...` on supplied trajectory/TPR copies. Engine exit zero requires review of retained diagnostics; it is not automatically labeled valid.

## What it preserves

Force field, protonation assumptions, water/ion/ligand choices, integration settings, coupling and existing coordinates come from the supplied prepared bundle. MD Health does not choose or repair them. It does not turn a raw PDB into a fully parameterized system, generate a force field, solvate, assign charges or rebuild missing residues. StructPrep still performs deterministic coordinate selection; StateAtlas still interprets mapped trajectory resemblance.

Changing test length is a new run, not an extension/restart with guaranteed continuity. Include any supplied checkpoint, restraint coordinates, index or energy state using its explicit manifest role. Random velocity generation and runtime differences may change trajectories; the actual resolved MDP and runtime outputs must be reviewed.

## Input contract

Start with [input-template.json](examples/input-template.json). Set the actual stage and supplied force field, then name `coordinates`, `topology` and `mdp` relative to the bundle directory. Supported optional file roles are `restraint_coordinates`, `index`, `checkpoint` and `energy`. Include all referenced `.itp`/force-field files locally. v1 requires GRO coordinates and refuses external MDP include-search paths. It does not claim full GROMACS syntax/physics validation.

The pilot step cap (1..1,000,000) and command wall-clock limit (1..3,600 seconds) are software execution bounds, not recommended physical durations. Set them for the selected stage and available runtime.

For diagnosis, export one named GROMACS metric dataset to XVG with `@ sN legend` labels and `Time (ps)` units. [policy-template.json](examples/policy-template.json) intentionally contains placeholders and refuses execution until supplied with finite stage/system-specific limits. State required metrics, expected recorded endpoint and output-cadence tolerance. Pressure fluctuations, RMSD changes and energy changes do not have universal failure thresholds. This first diagnostic does not estimate statistical convergence, drift significance or equilibrium.

## Commands

After a separately approved installation of the local candidate, `md-health` is the entry point. Without installation, from `pub-dev` prefix commands with `PYTHONPATH=software/md-health/src python -m md_health.cli`:

```bash
md-health describe
md-health validate --bundle /path/to/existing-bundle --manifest /path/to/input.json
md-health prepare --bundle /path/to/existing-bundle --manifest /path/to/input.json --pilot-steps YOUR_DECLARED_STEPS --out /path/to/new-prepared-copy
md-health pilot --prepared /path/to/new-prepared-copy --out /path/to/new-pilot --gmx /path/to/gmx --timeout-seconds YOUR_DECLARED_LIMIT
md-health check-trajectory --trajectory /path/to/run.xtc --tpr /path/to/run.tpr --out /path/to/new-integrity-check
md-health assess --log /path/to/captured-run.log --series /path/to/metrics.xvg --policy /path/to/declared-policy.json --run-receipt /path/to/MDRUN_RECEIPT.json --out /path/to/new-health-report
```

The pilot receipt binds its exit result to the captured log. Historical logs without that receipt can still reveal failures, but cannot receive the fully evidenced `no_failure_signal_observed` status. Changing files during assessment are refused: take a local captured copy for each early checkpoint. This is snapshot diagnosis, not a continuously running monitor. The prototype reads named XVG exports; it does not directly decode EDR, parse all logs or automatically run `gmx energy` selections.

`REPORT.html` explains observed signals, missing evidence, series ranges and source locations. `REPORT.json` preserves exact policy and fingerprints; `CHECKSUMS.json` identifies both reports. Existing output directories are never overwritten.

## Status meanings

| Result | Meaning |
|---|---|
| `failure_signal_observed` | A concrete fatal/failed-exit/non-finite signal was found in supplied evidence |
| `review_required` | Warning, declared-band departure or partial/time-order evidence needs inspection |
| `insufficient_evidence` | Required evidence is absent, even if no failure pattern was found |
| `no_failure_signal_observed` | No configured signal found in the supplied finite, labeled and sufficiently evidenced window |

The last status is not proof of stable production dynamics, equilibrium, convergence or biological validity. Warning detection is deliberately bounded; unfamiliar messages can be missed. A failing run is preserved for diagnosis, not silently adjusted until it passes.

`assess` exits 0 only for its last status, 1 for review/insufficient evidence and 2 for failure/refusal. `pilot` exit zero establishes only the recorded engine-exit condition. `check-trajectory` exits 1 even when GROMACS exits zero, to require diagnostic review.

## Checked software and open work

See [VERIFICATION.json](checks/VERIFICATION.json) and [the session report](SESSION_REPORT.md) for exact current checks. The diagnostics demonstration uses invented fixtures only: [open its local page](checks/demo-001/index.html). No private molecular data is included in that demonstration.

Real GROMACS preprocessing/pilot/check, actual input bundles, version/OS compatibility, installed runtime performance, interrupted-process descendants and full UI integration remain unverified. Only command-line/report integration is implemented. The complete prior Workbench qualification scope is unchanged.

## Official references

Command contracts and conservative warning handling were checked against official GROMACS 2026.4 documentation:

- [gmx grompp](https://manual.gromacs.org/documentation/2026.4/onlinehelp/gmx-grompp.html): preprocessing, include behavior, checkpoints/restraints and warnings.
- [gmx mdrun](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html): explicit run input, output and step/thread options.
- [gmx check](https://manual.gromacs.org/current/onlinehelp/gmx-check.html): trajectory/TPR and bond diagnostics.
- [Numerical instability terminology](https://manual.gromacs.org/documentation/2026.4/user-guide/terminology.html): constraint warnings can be symptoms of broader instability, rather than proof that constraints caused it.

## Privacy

Raw bundles and generated pilots/reports can contain unpublished structures and local provenance. Keep them outside publication discovery. `checks/` is excluded from this module's Git discovery. No dataset, installation, external execution, commit, push or publication is authorized by this README.
