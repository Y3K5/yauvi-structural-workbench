# StateAtlas

StateAtlas compares static structures, multi-model ensembles, or MD trajectories
to a two-sided reference set of experimentally supported conformations. Its
labels are deliberately bounded to structural resemblance:
`active_like`, `inactive_like`, `mixed`, or `unresolved`.

```bash
state-atlas run --manifest STRUCTURE_EVIDENCE.json --structure ensemble.pdb \
  --reference-set references.json --out out

# MD trajectory support requires: pip install 'yauvi-state-atlas[md]'
state-atlas run --manifest STRUCTURE_EVIDENCE.json --topology top.pdb \
  --trajectory run.xtc --pbc none --stride 10 \
  --reference-set references.json --out out
```

Reference sets must contain experimental active and inactive structures, state
evidence, a maximum interpretable RMSD, and a between-state margin. StateAtlas
reports unresolved frames in every denominator; it never hides them.
