# ABL StateAtlas panel: what it measures — DRAFT, not adopted

Written before any record is curated and before the engine is registered, per the
adoption protocol's closing requirement: state plainly what each gate measures and
confirm it is what the scope claims. The membrane panel satisfied every procedural
rule for its entire life while measuring rotational self-consistency and appearing
to measure orientation accuracy. This document exists to make that failure visible
here before it can happen.

Scope name: `conformational_state:abl_family`. Panel id
`qualification-v2-abl-state-atlas`. 18 records: 6 reference (2 active, 4 inactive)
and 12 held-out (4 active, 8 inactive).

Everything numeric below was measured, not reasoned about. Reproduce with:

    python3 abl_state_separability_screen.py

---

## 1. What the panel measures

**Whether StateAtlas assigns a kinase-domain conformation to the correct side of a
two-sided experimental reference set, and — above all — whether it ever makes a
confident call to the opposite state.**

The four declared gates split into three kinds:

| gate | value | kind | what it asserts |
|---|---|---|---|
| `max_best_reference_rmsd_A` | 2.5 | interpretability bound | the case resembles *some* reference closely enough to be interpretable at all |
| `minimum_opposite_state_margin_A` | 0.25 | separation | the two sides are far enough apart to distinguish |
| `confident_opposite_state_calls_max` | 0 | **false-positive** | a confident call is never to the wrong state |
| `correct_interpretable_coverage_per_state_min` | 0.8 | recall | at least 80 percent of each state's cases are both interpretable and correct |

The third is the scientifically load-bearing one and is an absolute bound: a single
confident opposite-state call fails the panel. The first is not a result about
conformation — it is the condition under which a result is allowed to exist — and
should be labelled an interpretability bound so a green panel is not read as four
scientific results.

The method is global CA RMSD after superposition over UniProt P00519 242-495,
against active and inactive reference structures. Labels are bounded to structural
resemblance: `active_like`, `inactive_like`, `mixed`, `unresolved`.

## 2. Pre-adoption measurement, and one hypothesis withdrawn

Sixteen human ABL1 entries, chain-mapped to P00519 by global sequence alignment
(coverage 0.921-1.000, identity >= 0.990), all 120 pairwise RMSDs computed.

**Withdrawn before it was written into anything.** The concern raised at the start
of this work was that global CA RMSD over a 254-residue mask would be dominated by
the invariant bulk of the kinase domain, since the active/inactive difference is
local — and would therefore measure overall structural similarity while appearing
to measure conformational state. That is the membrane failure mode, and it is
**not what the data shows.**

- Two clusters separate cleanly: within-cluster median **1.30 A**, between-cluster
  median **5.09 A**.
- The difference is localised where the biology says it should be. Mean per-residue
  deviation between clusters, after superposition:

      UniProt 382-397   16 residues, mean 16.4 A   activation loop, immediately C-terminal to DFG (D381-F382-G383)
      UniProt 274-279    6 residues, mean  5.2 A   beta3-alphaC

  Only 54 of 254 positions deviate by less than 1.0 A, and the top 15 positions are
  all in 384-397. The metric is reading the activation loop and the alphaC region —
  the canonical conformational switches — not incidental structure.

The metric measures what the scope claims. That is a positive result and it is
recorded as such.

## 3. The real finding: margin headroom is a property of the reference set

The frozen margin is 0.25 A. On this pool that threshold is not comfortably far
from the data — and *how* far depends entirely on which structures are chosen as
references, not on the held-out cases.

With all sixteen available as references, **5MO4 is the nearest opposite-group
reference for every member of cluster A**, and every cluster-A margin collapses
into 0.22-0.48 A. Two structures fall below the threshold and go `unresolved`:

    1OPL  margin 0.22 A
    5MO4  margin 0.20 A

Withhold 5MO4 from the reference set and nothing else changes, but the same
cluster-A margins become **2.03-2.71 A** — an order of magnitude of headroom, and
only one structure below threshold.

**Why this is a trap and not a tuning knob.** A curator choosing references without
this measurement in hand could produce a comfortably passing panel or a marginal
one, and the choice would not be visible anywhere as a scientific decision. That is
the membrane failure mode relocated from the gate to the curation step: the panel
would be measuring a reference-set choice while appearing to measure the module.

The rule that follows: **the reference set is frozen by a criterion stated in
advance, never by observed difficulty.** Excluding a structure *because* it makes
the panel harder is exactly the move the protocol forbids. If 5MO4 is excluded
from the references it must be for a declared property — it is an allosteric
myristoyl-site (asciminib) complex whose kinase domain sits ~1 A from both clusters
— and it should then appear as a held-out case, where being hard is the point.

Note also that **zero opposite-state calls occurred** on this pool under either
reference set. The load-bearing gate holds; it is the margin gate that is close to
the data.

## 4. What the panel does NOT measure

- **Not the autoinhibited assembly state.** The mask is 242-495, the kinase domain
  alone. 1OPL (UniProt span 6-512), 5MO4 (27-515), 4XEY (119-515) and 2FO0 (38-512)
  are multi-domain constructs whose SH3 and SH2 domains lie entirely outside the
  measured region. A structure whose inactive character is carried by the
  autoinhibitory assembly rather than by the kinase domain cannot be recognised as
  such by any method restricted to this mask. The panel must say so.
- **Not function, activity or drug response.** Structural resemblance to a reference
  is not catalytic state.
- **Not a claim about murine ABL1.** Seven canonical entries — 1IEP, 1FPU, 1OPJ,
  3KF4, 3KFA, 1M52, 3OXZ — map to P00520 (ABL1_MOUSE) and are excluded, because
  mixing accessions breaks the residue-equivalence basis the mask depends on. This
  removes several textbook DFG-out imatinib complexes from the pool; 2HYY is the
  human imatinib complex that remains.
- **Not conformational classification in the Dunbrack/Kincore sense.** This is a
  two-sided RMSD comparison, not a DFG/alphaC spatial-label assignment.

## 5. Freeze items (protocol rule 1)

- The **reference set**: six entries, chains, and the criterion by which they were
  chosen — fixed before any margin is observed on the held-out cases.
- The **mask**: UniProt P00519 242-495 at >= 0.90 coverage. Already frozen in
  `state_atlas.core` as `ABL_DOMAIN_START/END` with a validator that rejects any
  other mask for Mark 1.
- The **thresholds** `max_rmsd_A = 2.5` and `min_margin_A = 0.25`, already frozen in
  code as `ABL_MAX_RMSD_A` / `ABL_MIN_MARGIN_A` and enforced by
  `validate_reference_set`. Changing either is a code change and a new collection
  version, not a manifest edit.
- The **residue-equivalence route**. PDBe residue-level SIFTS returns
  `author_residue_number: null` for twelve of these sixteen entries, which silently
  drops them. The screen uses global sequence alignment to the P00519 canonical,
  matching what StateAtlas does internally. Whichever route the panel adopts must
  be declared, because the coverage rule permits either.
- The **state evidence source** — see below.

## 6. Open questions for Yuvraj

1. **Which cluster is active and which is inactive.** The screen deliberately does
   not say. Its two groups come from clustering the RMSD matrix, which establishes
   separability, not identity. Coverage rule 3 requires independent state evidence
   per reference. Kincore was unreachable from here (`dunbrack.fccc.edu` redirects,
   `kincore.research.fchampalimaud.org` does not resolve). The alternatives are
   per-structure literature assignment recorded with citations, or a reachable
   mirror. **This must be sourced, not inferred from the clustering** — deriving
   states from the same geometry the module uses is circular, and it is exactly how
   a false-positive gate passes for the wrong reason.
2. **5MO4's role**: reference or held-out case, decided by a stated property rather
   than by its effect on the margins. Recommendation above is held-out.
3. **The record count.** 18 records are required and 16 human entries were found.
   Either the pool grows, or some entries contribute more than one record via
   distinct chains — in which case the rule against a single entry appearing in
   both the reference and held-out splits needs writing down, because that is
   leakage.
4. **The 2:1 state imbalance.** Requirements ask for 2 active and 4 inactive
   references, 4 active and 8 inactive held-out. The empirical clusters are 8 and 8.
   Whichever cluster turns out to be active is short of what the inactive side
   needs, or over-supplied — check against the pool once states are assigned.
5. **Gate labelling**: whether `max_best_reference_rmsd_A` is recorded as an
   interpretability bound rather than a scientific gate, as §1 proposes.

## 7. Blocking work, in order

1. Decide §6.1 — the state evidence source. Nothing else can proceed honestly.
2. Decide §6.2 and §6.5.
3. Register `conformational_state` in the runner's `ENGINES`. It is absent; the
   `state-atlas` CLI is `state-atlas run --manifest ... --structure ...
   --reference-set ... --out ...`, and the evidence document name must be declared.
4. Write `gate_semantics` for the panel. It is currently empty and the runner
   refuses to run without it — correctly, since it will not assume gate meanings.
5. Curate 6 reference and 12 held-out records against the frozen criterion.
6. Prove the gate can fail (protocol rule 2): corrupt a case on a copy and confirm
   rejection, and confirm an opposite-state call is actually caught.
7. Cross-machine reproduction before adoption (protocol rule 7).

---

## 8. Quaternary context — a third state the two-sided design cannot hold

Added 2026-09-01 after Yuvraj proposed using assembled complexes and quaternary
structure. The suggestion turned out to identify a panel-defining problem, so it is
recorded here in full rather than as a note.

### 8.1 Ligand chemistry is independent state evidence, and it confirms the geometry

Bound ligands were resolved from the PDBe compound API and are evidence of
conformation that does not come from the coordinates the module measures — which
is what protocol rule 4 asks for. Assignments below use the standard type I
(DFG-in) / type II (DFG-out) binding-mode distinction.

| cluster | entry | ligand evidence | implied kinase-domain conformation |
|---|---|---|---|
| B | 2HYY | STI, imatinib | DFG-out |
| B | 3CS9 | NIL, nilotinib | DFG-out |
| B | 3QRI, 3QRJ | 919, type II pyrazole-quinoline | DFG-out |
| B | 2E2B | 406, bipyrimidine type II | DFG-out |
| B | 5MO4 | AY7 asciminib (allosteric myristoyl site) + NIL nilotinib | DFG-out, allosterically occupied |
| A | 2GQG | 1N1 dasatinib **+ PTR, phosphotyrosine** | DFG-in, activation loop phosphorylated |
| A | 2G2I | ADP | DFG-in |
| A | 4WA9 | AXI, axitinib | DFG-in |
| A | 4XEY, 2F4J | 1N1 dasatinib / VX6 type I | DFG-in |

Cluster B is the type-II / DFG-out set and cluster A is the type-I / DFG-in set.
The empirical clustering of §2 and the ligand chemistry agree, by two independent
routes. **§6.1 is answered for these twelve entries**: cluster A is active-like,
cluster B is inactive, and each assignment is citable per structure rather than
inferred from the RMSD matrix.

### 8.2 Where it breaks, and why that matters more than where it agrees

**1OPL and 2FO0 carry myristic acid (MYR) and 180 / 196 N-terminal residues in
contact with the kinase domain — the assembled, SH3-SH2-clamped, myristoylated
autoinhibited form. They are biologically inactive. They sit in cluster A,
0.78-1.7 A from the phosphorylated active structure and ~5.1 A from imatinib-bound
2HYY.**

That is not an error in the module. The autoinhibited kinase domain is DFG-in; what
holds it inactive is the quaternary assembly, and the assembly lies entirely
outside the 242-495 mask. The module reports the kinase-domain conformation
correctly and the biological state incorrectly, and it cannot do otherwise.

The consequence is concrete and would have surfaced only after curation: **a
curator labelling 1OPL or 2FO0 `inactive` on functional grounds produces a
confident opposite-state call, and `confident_opposite_state_calls_max = 0` fails
the panel — for a structurally correct answer.** The panel would have reported a
module defect that is really a scope error in the record.

So ABL1 presents at least three distinguishable states where the panel design
allows two:

1. **DFG-out**, type II inhibitor bound — inactive, kinase domain alone.
2. **DFG-in**, type I bound or phosphorylated or ADP-loaded — active-like.
3. **Assembled autoinhibited** — DFG-in kinase domain, SH3-SH2 docked, myristate
   in the C-lobe pocket. Inactive as a molecule, cluster A as a domain.

5MO4 is the fourth combination and explains its position between the clusters: an
assembled construct (178 N-terminal residues) with the allosteric myristoyl site
occupied by asciminib *and* nilotinib in the ATP site. Its 0.20 A margin is not
noise. It is a structure that genuinely belongs to both sides of a two-sided model.

### 8.3 What this requires of the panel

- **Declare the state axis as kinase-domain conformation**, DFG-in-like versus
  DFG-out-like, not biological activity. The gates then mean what they say.
- **Record quaternary context per record** as a declared covariate, never as a
  gate: N-terminal residue span, SH3-SH2 contact count, and myristoyl-site
  occupancy. This is Yuvraj's assembled-complex proposal in the form the panel can
  use, and it is what makes the 1OPL / 2FO0 / 5MO4 cases legible instead of
  anomalous.
- **Do not curate 1OPL, 2FO0 or 5MO4 as `inactive`.** Either exclude them by a
  stated rule, or curate them by conformation with the assembly context recorded.
- The platform already has an adopted `assembly_interface:deposited_biological_assembly`
  scope at 16/16 cross-machine. If quaternary context is to be recorded, it should
  come from there rather than from a new mechanism.

### 8.4 Recorded, not acted on

2G1T carries 112, an ATP-gamma-S peptide conjugate — a substrate-mimicking
nucleotide, which would suggest an active-like domain — yet it clusters with the
DFG-out set and is 3.1-4.4 A from everything, forming its own cluster at k=3. Its
margin is 1.75 A, so it would not fail the panel, but it does not fit either
account above. Worth resolving before it is curated, and worth not explaining away.

## 9. The grooves, measured — and why the mask cannot see the energy

Three binding grooves resolve from the sixteen structures with no prior
assignment, by superposing every entry onto 1OPL over UniProt 242-495 and
carrying each ligand through the same transform.

**Myristoyl pocket (C-lobe).** Three occupants within 2 A of each other:
MYR in 1OPL (reference), MYR in 2FO0 (0.9 A), and **AY7 asciminib in 5MO4
(1.9 A)**. The allosteric drug occupies the natural myristate groove — measured
here, not assumed.

**ATP cleft.** Every ATP-site ligand lands 0.3-7.5 A from the reference, and the
spread inside the cleft separates the two binding modes cleanly, with no overlap:

    type I  (DFG-in)   P16 0.3, AXI 1.5, VX6 2.4, 1N1 2.6/2.8, ADP 3.9
    -- gap --
    type II (DFG-out)  NIL 6.0/6.9, STI 7.0, 919 7.2/7.5, 406 7.4

The ~4 A centroid shift is the **DFG-out back pocket**: type II compounds extend
past the gatekeeper into space that only exists when DFG is out. This is a third
independent confirmation of the §8.1 state assignment, from ligand geometry rather
than from backbone RMSD or from compound identity.

**2G1T sits in the gap at 5.1 A**, between the two modes — consistent with its
intermediate position in the RMSD matrix and its own cluster at k=3. The anomaly
of §8.4 is internally consistent on both axes rather than noise.

**SH3-SH2 interface.** Buried surface area per side, Shrake-Rupley:

    5MO4  178 N-terminal residues   1355 A^2   7.6 per residue
    1OPL  180                       1335 A^2   7.4
    2FO0  195                       1276 A^2   6.5
    4XEY  105 (SH2 only)             892 A^2   8.5

~1300 A^2 per side is a substantial interface, in the range of a genuine
protein-protein interaction rather than a crystal contact.

### What this means for the panel

The energy holding 1OPL, 2FO0 and 5MO4 inactive is carried by two mechanisms, and
**a CA-RMSD comparison over 242-495 can see neither**:

- The **SH3-SH2 clamp** lies outside the mask entirely.
- The **myristoyl groove is inside the mask**, but occupancy is a ligand, and a
  backbone RMSD cannot register whether a pocket is filled.

So §8.3's covariate is better recorded as **buried interface area and myristoyl-site
occupancy** than as the raw contact counts first proposed. Both are cheap, both are
reproducible, and together they state exactly what the panel is blind to.

## 10. Machinery complete, 2026-09-01 — and it can fail

The four pieces adoption needs are in place: engine registration and command
builder, a measurement function, a gate function, and a coverage witness
(`MEASURERS`, `GATE_CHECKS`, `WITNESSES` in `run_execution.py`).

Demonstrated against real coordinates rather than a fixture. References 2GQG and
7W7X (active), 2HYY and 3QRI (inactive); query 3CS9, nilotinib-bound, DFG-out:

    call inactive_like   best reference INACTIVE_2HYY at 0.440 A   margin 4.534 A
    rmsd to each reference: 2GQG 5.081, 7W7X 4.974, 2HYY 0.440, 3QRI 1.029

**Curated as an inactive record, all six per-case checks pass. Curated as an
active record, exactly one fails — `confident_opposite_state_call`.** That is
protocol rule 2 satisfied on the gate the panel exists for, and the failure it
catches is precisely the one section 8 predicted: a structure whose biological
inactivity is carried by something the mask cannot see, curated by activity
rather than by conformation.

Three design points worth stating, because each was a choice:

- **The margin gate binds only on a confident call.** An unresolved frame has no
  margin claim to make, and requiring one would turn "I cannot tell" into a
  failure — when reporting it is the behaviour the module promises.
- **Both sides must have been compared.** A margin computed against one side is
  not a margin, and the module would still emit a call. Without this check a
  panel could pass while measuring nothing about separation, so it is also a
  required coverage feature.
- **The alignment map checksum is gated.** A run that lost it still produces
  numbers, from a different set of residues.

`correct_interpretable_coverage_per_state_min` is the one gate that is not per
case. It is a fraction over a stratum, and the runner's `stratum_state` already
requires every case to pass, which is stricter. It is declared and **not yet
enforced anywhere**; enforcing it needs an aggregate step in
`run_execution.main`, and it only starts to bind once records are permitted to be
unresolved and still pass.

What remains before adoption: curate the records under the section 6 rules,
source-lock every structure (`execution_policy.network_access` is `forbidden`, so
nothing may be fetched at run time), and reproduce cross-machine.
