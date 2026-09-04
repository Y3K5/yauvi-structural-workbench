# ABL StateAtlas record selection — pool, evidence, and a feasibility problem

Screened 2026-09-01, before any record is frozen. Companion to
`ABL_STATE_ATLAS_PANEL_SCOPE.draft.md`. Nothing here is curated yet, and §5 is the
reason.

## 1. The pool

All 85 experimental human ABL1 (UniProt P00519) entries were enumerated from the
RCSB search API. Murine ABL1 (P00520) is excluded throughout — 1IEP, 1FPU, 1OPJ,
3KF4, 3KFA, 1M52 and 3OXZ are murine, which removes several textbook DFG-out
imatinib complexes; 2HYY is the human imatinib complex that remains.

    85 human P00519 entries
    47 whose mapping spans the frozen mask (UniProt 242-495)
    47 pass coverage >= 0.90 and sequence identity >= 0.95

All 47 qualify. Chain-mapped by global alignment to the P00519 canonical, not by
PDBe residue-level SIFTS, which returns `author_residue_number: null` for most of
these entries and silently drops them.

## 2. Three independent evidence axes

State is never assigned from the geometry the module itself measures — that is
circular, and it is how a false-positive gate passes for the wrong reason. Each
record carries up to three axes that can be checked separately:

1. **Backbone clustering** — hierarchical clustering of the 47x47 CA-RMSD matrix.
2. **Ligand identity** — the chemistry of what is bound, resolved from the PDBe
   compound API. Type II (DFG-out) versus type I (DFG-in) binding mode.
3. **Ligand position** — every entry superposed into the 1OPL frame; the ATP-site
   ligand centroid separates the modes with a clean gap (type I 0.2-3.9 A, type II
   6.0-10.3 A from reference). The ~4 A offset is the DFG-out back pocket.

A fourth marker is decisive where present: **phosphotyrosine (PTR) on the
activation loop**. Ten entries carry it, all in one cluster and none in the other,
which is what orients the clustering rather than any assumption of mine.

## 3. Pool by state

    active_like   29    DFG-in cluster; type I ligands, ADP, or pTyr
    inactive      14    DFG-out cluster; type II ligands
    assembled      4    SH3-SH2 clamped with the myristoyl groove filled
                        (1OPL, 2FO0, 5MO4, 8SSN — the third state of scope draft section 8)

## 4. Discordances, recorded not resolved

Four entries where the axes disagree. None should be curated until resolved.

| entry | cluster | ligand mode | note |
|---|---|---|---|
| 2HZ0 | active-like | type II (6.5 A) | the only clear cross-axis conflict |
| 2G2F | active-like | between (5.3 A) | ATP-gamma-S conjugate, as 2G1T |
| 2G1T | inactive | between (5.1 A) | outlier in the matrix, 3.1-4.1 A from all inactive entries |
| 2HIW | inactive | between (5.8 A) | |

Also noted: **7N9G** registers myristoyl-groove occupancy with only 9 N-terminal
residues, so a ligand sits in the lipid pocket without the SH3-SH2 clamp. 3PYY,
6NPE, 6NPU and 6NPV do the same with imatinib in the ATP cleft — dual occupancy,
consistent and expected for allosteric myristoyl-site binders, but it means
"myristoyl pocket filled" and "assembled" are two different properties and must be
recorded as two fields.

## 5. The feasibility problem

The panel requires **12 inactive records** (4 reference, 8 held-out) and 6 active.
The pool has 14 inactive entries, so the requirement is satisfiable **by count**.
It is not satisfiable by independent evidence.

Pairwise RMSD within the inactive side collapses it into roughly seven distinct
conformations:

    {2E2B, 2HYY, 3CS9}          0.26 - 0.44 A apart
    {3PYY, 6NPE, 6NPU, 6NPV}    0.14 - 0.26 A apart
    {3QRI, 3QRJ, 3QRK}          0.51 - 0.68 A apart
    2HIW, 5HU9, 6XRG, 2G1T      distinct

Nine near-duplicate pairs under 0.50 A, involving 7 of the 14 entries. **Twelve
inactive records would represent about seven independent conformations, four of
them copies separated by 0.14-0.26 A** — below the panel's own 0.25 A margin
threshold, which is to say indistinguishable by the very quantity being gated.

A panel reporting 12/12 would be reporting far less evidence than the number
implies. That is a count overstating its own strength, which is the failure this
programme exists to catch.

Only 11 of the 14 carry concordant ligand-mode support; 2G1T, 2HIW and 6XRG do
not. So the strict pool is 11 records for a 12-record requirement.

The active side has room — 29 entries — but 16 near-duplicate pairs sit inside it,
so the same rule must apply there.

## 6. What this needs from Yuvraj

The splits are frozen: `immutability_policy` says splits are never edited in
place, and a scientific change requires a new collection version. So this is a
decision, not an edit.

1. **Add a minimum within-stratum separation to the coverage rules** — records in
   one stratum must differ by more than the gated margin, otherwise the panel
   counts copies as evidence. Recommended, and it is the rule that makes §5
   visible in the manifest rather than only in this document.
2. **Then rebalance the split.** With a separation rule the inactive side supports
   roughly 7 independent records, not 12. Either the 4/8 inactive requirement
   drops, or the panel accepts fewer records per stratum and says so, or the
   inactive side is broadened beyond ABL1 — which changes the scope name.
3. **Assembled structures**: 1OPL, 2FO0, 5MO4 and 8SSN are a third state and must
   not be curated as `inactive`. Excluded by a stated rule, or curated by kinase
   conformation with assembly context recorded.
4. Resolve or exclude the four discordances in §4.

Nothing should be curated before 1 and 2 are settled, because both change how many
records exist and what they are allowed to be.

---

## 7. Found by running the module, 2026-09-01

Everything above was reached by screening coordinates. These were reached by
running `state-atlas` itself, and two of them contradict what was written earlier
in this session.

1. **`state-atlas` was not installed.** Five of the six engine CLIs were present in
   the working venv; this one was not, and `import state_atlas` failed. An
   `ENGINES` entry naming it would have died at "state-atlas is not on PATH" the
   first time the panel ran. Installed from local source with `--no-deps`.
2. **Reference Set v2 requires at least two references per state.** The 2 + 2
   reference requirement in collection 2.5 is exactly the module's floor. There is
   no slack: losing one reference on either side makes the set invalid.
3. **Modified residues are invisible to the module.** StateAtlas reads coordinates
   through Bio.PDB, which skips HETATM records, so phosphotyrosine is in the file
   and absent from the model. An alignment map declaring UniProt 393 exact for
   2GQG or 7W7X is rejected outright. This lands on exactly the entries proposed
   as independent active-state evidence: **their maps must declare 253 of 254
   positions, not 254.** Coverage stays well above the 0.90 rule, but the map must
   be built from standard residues only.
4. **The path works end to end.** References 2GQG and 7W7X (active), 2HYY and 3QRI
   (inactive), query 3CS9. Result `inactive_like`, best reference 2HYY at
   0.440 A, margin 4.534 A against a 0.25 A threshold. The reported RMSDs agree
   with the independent screen in `abl_state_separability_screen.py` to about
   0.1 A, which is two routes to the same number rather than one.
5. **It fails closed.** Corrupting a reference `structure_sha256` is rejected with
   a checksum mismatch rather than run.
6. **The runner still cannot execute this panel.** `MEASURERS` in
   `run_execution.py` has no `conformational_state` entry, so the command builds
   and the result cannot be read back. A measurement function and a gate function
   remain, which is what the session primer meant by machinery plus curation.

The smoke-test reference set and alignment map are not committed here: they were
built to answer whether the path works, not to be the panel's records, and the
panel's references must be chosen under the §6 rules rather than for convenience.

---

## 8. Records curated, 2026-09-01 — `ADOPTION_DRAFT_ABL.json`

Fourteen records and one control, source-locked, expectations recorded by
execution rather than transcribed. **Not merged into `PANEL_MANIFEST.json`:**
adoption needs independent reproduction on a second machine (protocol rule 7),
which has not happened.

### Selection rules, in the order they were applied

1. Human P00519, discordant entries (2HZ0, 2G2F, 2G1T, 2HIW) and the four
   quaternary-autoinhibited entries excluded.
2. State from ligand binding mode and activation-loop phosphorylation, never from
   the RMSD the module computes.
3. **References are the two most representative members of a stratum that are
   more than 1.0 A apart.** The first pass chose 6NPU and 3PYY, two members of
   the same near-duplicate cluster 0.26 A apart — a reference set pinned to one
   crystal form. An order of magnitude above the 0.25 A gate is the rule that
   prevents it.
4. Held-out records by farthest-point spread, each more than 0.25 A from every
   other record in its stratum.

    active    references 4XEY, 7W7X (1.12 A apart)   held-out 2GQG, 7N9G, 6BL8, 4TWP, 2G2I
              stratum min separation 0.924 A, max spread 1.692 A, 2 of 7 phosphorylated
    inactive  references 2E2B, 3QRK (1.09 A apart)   held-out 6NPE, 3QRJ, 3QRI, 2HYY, 3CS9
              stratum min separation 0.261 A, max spread 3.420 A, all type II

Only 2 of 7 active records carry phosphotyrosine, deliberately: a stratum where
every member is phosphorylated would conflate DFG-in with phosphorylation.

### Result

    14/14 cases passed   1/1 controls passed   6/6 coverage features witnessed

Both references self-match at 0.000 A with margins of 3.3 to 5.0 A. Every held-out
call is on the correct side. Zero confident opposite-state calls.

### 5HU9: a case that became a control

The farthest-point rule selected 5HU9 for the inactive stratum, and on execution
it returned **unresolved at 2.883 A** — past the 2.5 A interpretability bound.
That is the module behaving correctly and the selection rule being wrong:
**maximising spread and respecting the interpretability bound pull in opposite
directions**, and nothing in section 6 said so.

5HU9 is now a fail-closed control with that as its declared property: a structure
the reference set cannot speak for must be reported unresolved with
`interpretable_state_frames` withheld, rather than forced onto the nearer side.
3CS9 takes the fifth inactive held-out slot. The control was kept because of what
it demonstrates, not discarded because it failed.

### Protocol rule 2 — two independent falsifications

- **Mis-curated stratum.** 2HYY relabelled `active` on a copy: fails on
  `confident_opposite_state_call`, and only that. This is the 1OPL / 2FO0 error
  from scope draft section 8, caught by the panel.
- **Tampered expectation.** 3QRI's recorded `margin_A` set to 99.0 on a copy:
  fails on `margin_A_matches_recorded_expectation`.

The second failed to fail on the first attempt. Only the *call* was compared
against the recorded expectation, so the recorded distances were inert — a run
could drift by an angstrom, keep the same side, and pass. `best_rmsd_A`,
`margin_A` and `best_reference` are now all compared at a declared 0.001 A
tolerance, which is the quantity the second-machine gate will have to reproduce.

### Two harness facts found only by executing

- Every run withholds `qualification_v2_held_out_panel` — the module declining to
  claim a qualification this panel has not earned. `expected_missing_by_stratum`
  must list it; the empty list first written there would have failed all fourteen.
- `missing_evidence_behavior` and `coverage_requirements` are read from the top
  level of `gate_semantics`, beside the workflow key, not from inside it. Nested,
  they raised `KeyError` on the first real execution.

### What adoption still needs

Independent reproduction on a second machine, and review of the state assignments
by someone other than their author.
