# memorient documentation

- [SPEC.md](SPEC.md) — module-by-module contract: signatures, guarantees, design invariants.
- [../README.md](../README.md) — install, quickstart, contexts, correctness summary.
- [../examples/worked_examples.py](../examples/worked_examples.py) — OmpA β-barrel + glycophorin A TM helix, end to end.
- [../examples/p4_benchmark.py](../examples/p4_benchmark.py) — correctness benchmark vs OPM (mean normal error 7.4°).

## Pipeline at a glance

```
load_structure ─▶ canonicalize ─▶ compute_sasa ─┐
                                                 ▼
   context.orientation_method ──▶ fit_membrane ──▶ classify ──▶ call_extracellular_side
                                                 ▼
                        reframe (+Z = extracellular) ──▶ project_membrane ──▶ label_residues
                                                 ▼
                        context_metrics + five_fold_validate ──▶ OrientationResult
                                                 ▼
                        viz: display_oriented / write_3dmol_html / write_pymol_script
```
