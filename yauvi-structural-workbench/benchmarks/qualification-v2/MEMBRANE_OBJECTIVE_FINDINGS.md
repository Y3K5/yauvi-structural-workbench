# Membrane orientation: current state of knowledge

Written 2026-08-31. `PANEL_MANIFEST.json` carries the full provenance —
`threshold_revisions` for gate changes, `record_corrections` for claims made and
withdrawn during the investigation. That trail is the evidence; this file is the
synthesis, so a reader does not have to reconstruct the conclusions from it.

Nothing here is inference. Each statement names the measurement behind it.

---

## 1. What the panel measures, and what it used to

Until collection 2.3 every normal gate in this panel measured **rotational
self-consistency** — whether the fit reproduces when the input is rotated.
`orientor.py` says so in as many words: *"Self-consistency, not correctness — a
stable wrong answer still passes."*

`half_thickness_error_A` was the only comparison against OPM, and it constrains
slab thickness, not direction. There was no orientation-accuracy gate, and the
fitted normal was not even recorded in `MEMBRANE_ORIENTATION.json` — the
quantity to be judged was never written down.

Collection 2.3 added `opm_normal_error_deg` (bound 2.0°, reference `[0,0,1]`) and
`summary.normal`. The stratum went from 14/16 to 5/16. The software did not
change; what the panel can ask did.

## 2. Accuracy, measured

`opm_normal_error_deg`, undirected angle to OPM's deposited z-axis:

| | |
|---|---|
| within 0.14° | 1PRN 0.0010, 2OMF 0.0067, 2POR 0.0511, 1E54 0.0803, 2MPR 0.1384 |
| 2.5°–17.4° | 1QD6, 1BXW, 1FEP, 2F1C, 1TLY, 1K24, 1P4T, 2ERV, 1T16, 1UYN, 2F1T |

Bimodal, with a 2.341° empty band. The fit either recovers the orientation
essentially exactly or misses by degrees. Nine of the eleven failing values are
identical across Rosetta x86_64, macOS arm64 and ubuntu x64; 1BXW, 1QD6, 2F1T,
2F1C and 2MPR vary by environment.

## 3. The correct answer always exists

Polishing from OPM's z with `fatol=1e-8` stays at z for **all 16** — drift below
0.13° for fifteen, 0.70° for 2F1C. The objective has a local optimum at the
correct orientation in every case. The problem is ranking, not representation.

An earlier claim that five structures admitted no optimum near OPM was an
artifact of an 80-point grid at ~20° spacing with loose polish, and is retracted.

## 4. Which term is responsible

Decomposing the ranking gap across all ten mis-ranked cases:

| dominant term | cases |
|---|---|
| `delta_kd` | 8 — 1BXW, 1P4T, 1K24, 1PRN, 1UYN, 2F1C, 2ERV, 2MPR |
| `belt` | 2 — 2F1T (dkd −0.512, belt **+0.431**), 1T16 (dkd −0.136, belt **+0.152**) |

Both hydrophobicity terms can prefer a tilted slab; which one does is
structure-dependent. An earlier claim that `delta_kd` was solely responsible came
from two cases and is withdrawn.

`girdle` — the only term encoding barrel-specific interface structure rather than
generic hydrophobic contrast — favours the correct orientation in 8 of 10 but
contributes at most 0.058 against gaps of 0.065 to 0.786. At weight 0.5 it is
effectively inert.

## 5. Re-weighting cannot fix it — proven

`J` is linear in its weights, so requiring the reference to outscore every
alternative optimum is a linear program:

```
maximise t  subject to  w · (c_z − c_alt) ≥ t   for every structure, every alternative
```

1303 constraints over 16 structures. Best margin across eight normalisations
(each weight pinned to ±1, excluding the trivial `w = 0`): **−0.158825.
Infeasible.** Leave-one-out remains infeasible for every case.

This is a proof, not a search result. Alternatives were located under the current
weights; re-optimising under any other weighting can only surface more optima,
adding constraints and never relaxing them. Infeasibility here implies
infeasibility on the true optimum set.

Every structure is individually satisfiable — they cannot be satisfied jointly.
A greedy maximal feasible subset reaches 6 of 16, and requires `girdle = −18.7`,
inverting the aromatic girdle term. That is the LP solving a problem that is not
numerical.

**Consequence: repairing this scope requires a new term, or a change to the
functional form of the existing ones. Not a re-weighting.**

## 6. Eliminated approaches

- **Re-weighting** — proven infeasible, §5.
- **Deterministic basin selection** (best-of-all-starts) — would return the wrong
  orientation for the 10 mis-ranked cases, since the objective ranks them first.
- **Barrel-axis alignment term** — for a β-barrel the membrane normal is the
  barrel axis, so this looks compelling. Principal axes of Cα track OPM to within
  0.16° for the five disc-like trimers, which already pass, and err by 6.1°–26.3°
  for the monomeric barrels that fail. PCA is contaminated by extramembrane loops
  and plugs. Accurate only where it is not needed.
- **Finer search** — `n_scan` 80 → 400 makes rotation drift worse, not better.

## 7. Open

What physical or geometric signal distinguishes a correctly placed bilayer for
monomeric β-barrels, which the four present terms do not capture. This is a
question about membrane protein biology.

A note on how to read §4 and §6: seven claims were proposed and withdrawn during
this investigation, each after a measurement that agreed with the explanation
then current, each corrected when a second measurement disagreed. The surviving
statements are the ones checked against source or against an independent
measurement. `record_corrections` lists what was withdrawn and why.

## 8. Related defect

`fit_membrane` is not rotation-equivariant: rotating the input moves the
back-rotated normal by up to 16°, and rotated fits often score higher than the
unrotated one. Pinned by `test_fit_membrane_is_rotation_equivariant` (strict
xfail) in `memorient/tests/test_barrel.py`.

It is visible in the production path, not only under synthetic rotation:
`fit_membrane` on 1T16 gives 11.08°, while the pipeline — which fits a
canonicalised structure and maps back — gives 16.30°. Same file, same 427
residues.
