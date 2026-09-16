# sf-csa pre-adoption findings — two corrections to the panel scope draft

Written before curation, before any code change. Both findings were reached by
running the module, not by reading it, and both revise
`SF_CSA_PANEL_SCOPE.draft.md`. Neither is a result about protein biology; both
are results about what the declared gates would actually measure.

Verification environment: `$TMPDIR/joss-venv`, `sf_csa` editable from the working
copy. Runtimes verified present and matching the v1 pins: `foldseek 10.941cd33`,
`diamond 2.1.11`.

---

## Finding 1 — `probable_same_function` is reachable, and reachable the wrong way

**The draft says the label "cannot be emitted". That is withdrawn.**

The draft's two stated causes are both real: the RBH computation runs after
classification (`core.py:591`) and writes to `h["orthology_status"]` on sequence
rows (`core.py:598`), not to `target_meta`. They correctly explain why the
*computed* reciprocal-best-hit never reaches the label.

They do not make the label unreachable. A third path exists:

| step | location | effect |
|---|---|---|
| query records validated | `core.py:517` | `row = {**q, ...}` — every curator key preserved verbatim |
| target metadata built | `core.py:533` | `target_meta = {q["accession"]: q for q in validated}` |
| classification reads it | `core.py:348` | `target_meta.get("rbh")` |
| manifest validation | `manifests.py:89` | checks only for *missing* required fields; never rejects unknown keys |

So `"rbh": true` written into a query record in the queries manifest passes
validation, survives into `target_meta`, and makes the module emit
`probable_same_function` with the stated basis *"RBH plus compatible whole
architecture and mechanism annotation"* — for a pair on which no reciprocal best
hit was ever computed.

Confirmed by direct call (`classify_hit` with `{"mechanism_group": ..., "rbh": True}`,
category `whole_architecture_match`) returning `probable_same_function`.

**Why this matters more than a dead path.** A dead path fails safe: the label is
never emitted, and the panel's false-positive bound holds trivially. This path
fails open. The strongest non-exact functional label in the closed vocabulary can
be switched on by a manifest field, carrying a basis string that is false, and
`analogy_or_unrelated_promoted_to_function_max = 0` — the panel's one
scientifically load-bearing gate — no longer requires any computed evidence to
be satisfied or violated.

**Consequence for the draft's recommendation.** Option B was "include a control
that *should* produce `probable_same_function`; it will fail; the panel correctly
reports a known-broken path." The natural way to curate that control is to assert
the relationship in the manifest — which is exactly the field above. The control
would then **pass**, and the panel would report that the computed RBH path works
when it has never run. That is the membrane failure mode restated: a gate that
appears to measure one thing while measuring curator input.

Option B is therefore not safe as written. See the decision note below.

## Finding 2 — with default tables, the load-bearing gate cannot fail

`DEFAULT_MECHANISM_FAMILIES` (`core.py:295`) is periodontal / Bacteroidetes outer
membrane biology: `omp85_bama`, `susc_raga_importer`, `t9ss_porg`, `msp_contested`,
`generic_om_barrel`. For any header outside that biology, `classify_title` returns
`"unknown"`.

`classify_hit` promotes only under `qgroup == tgroup and qgroup != "unknown"`
(`core.py:344`). Measured, with default tables and non-periodontal families:

| case | target group | label |
|---|---|---|
| periodontal homolog | `omp85_bama` | `same_mechanism_class` |
| **true homolog** (serine protease) | `unknown` | `structural_analogy_only` |
| **true homolog** (globin) | `unknown` | `structural_analogy_only` |
| **unrelated pair** | `unknown` | `structural_analogy_only` |
| exact self-match | `unknown` | `exact_function_supported` |

A true homolog and an unrelated pair are indistinguishable. The consequences for
the four declared gates are not symmetric:

- `exact_controls_recovered` — **still passes.** Exact self-match is tested first,
  by accession, before any mechanism-group logic.
- `structural_homolog_controls_recovered = all` — **fails on every case.**
  `same_mechanism_class` is unreachable when `tgroup` is `unknown`.
- `analogy_or_unrelated_promoted_to_function_max = 0` — **passes vacuously.**
  Nothing can be promoted at all, so the bound cannot be violated. By protocol
  rule 2, this is not a gate.

So an unoverridden default table does not merely use "the wrong biology". It
produces a panel that fails its recall gate 16/16 while its false-positive gate
reports success without measuring anything.

**Where discrimination does exist.** `target_meta` contains campaign query
accessions only. For `experimental_pdb` hits, `meta` is always `None` and the
periodontal regex is the sole source of `tgroup`. For `campaign_models` hits the
target *is* another campaign query, so `tgroup` is the curator-supplied
`mechanism_group` and classification separates correctly:

    same curated group      -> same_mechanism_class
    different curated group -> structural_analogy_only
    below threshold         -> unresolved_or_conflicted

The panel has discriminating power along the query-vs-query axis, and none along
the query-vs-experimental-PDB axis unless the tables are overridden.

Secondary observation, recorded not acted on: a `domain_or_partial_match` with a
matching curated group also returns `same_mechanism_class`. Partial and whole
architecture matches carry the same label, so the `homologous_superfamily`
stratum cannot distinguish them.

---

## What these findings change

The draft's freeze list already required overriding the mechanism tables. Finding
2 raises that from a correctness precaution to a precondition: without it, two of
four gates are meaningless in opposite directions.

Finding 1 changes the shape of the §4 decision. The choice is no longer A-or-B:

- **A — exclude the label.** Map `homologous_superfamily` to
  `same_mechanism_class`; record the exclusion and why. Defensible, but leaves the
  fail-open manifest path unmeasured and undocumented in the panel.
- **B — gate it as drafted.** Now known unsafe: the control most likely passes via
  the manifest field and certifies a path that never ran.
- **C — gate the defect as a defect.** Two controls, not one:
  a *positive* control with `rbh` absent, expected to fail, proving the computed
  path does not reach the label; and a *negative* control with `"rbh": true` on a
  pair with no reciprocal best hit, expected to be **rejected**, proving the panel
  distinguishes computed evidence from asserted evidence. This is the only option
  under which the false-positive gate is a real gate.

C is the recommendation. It costs two record slots and requires deciding whether
the negative control's expected rejection is enforced by the panel or by a fix to
the module.

## Withdrawn from the draft

- "`probable_same_function` cannot be emitted" — withdrawn; reachable via the
  queries manifest (Finding 1).
- "the panel must override [the tables] or it is silently using the wrong table" —
  understated rather than wrong; the failure is not silent misuse but a vacuous
  gate plus a total recall failure (Finding 2).

Both original statements are left in place in the draft, with a pointer here.

---

## Finding 3 — the flag is not pairwise, so a naive fix produces a working defect

Found 2026-09-01 while implementing the decision below, by reading the call site
rather than the classifier.

`target_meta = {q["accession"]: q for q in validated}` (`core.py:533`) is keyed by
**target accession alone**, and `classify_hit` read `target_meta.get("rbh")` off
that row. Reciprocal-best-hit is a relation between one query and one target: the
same target may be reciprocal for query A and not for query B.

So an `rbh` flag stored on a target row is query-independent by construction. The
repair implied by Findings 1 and 2 — move the RBH computation earlier and write
it into `target_meta` — would have produced a path that runs, returns the label,
and is still wrong: the first query to establish RBH against a target would flag
that target for every other query hitting it.

This is why "wrong order, wrong destination" was an incomplete diagnosis. There
were three defects, and the third is the only one that survives a fix aimed at
the first two.

## Decisions — 2026-09-01, Yuvraj

**§4: fix the module first.** Not A, B or C as drafted. With the module repaired,
the question the three options were dividing no longer exists: the label becomes
reachable by measurement and unreachable by assertion, so `homologous_superfamily`
can map to it on evidence, and the panel needs no control to route around a
defect that is gone. Option A's objection — that it "leaves the fail-open manifest
path unmeasured and undocumented" — is answered by the manifest gate below rather
than by a curated record.

**The gate `structural_and_sequence_outputs = separate` is labelled a contract
check**, not a scientific gate. A passing panel reads as three scientific results
plus one contract check. See `SF_CSA_PANEL_SCOPE.draft.md` §1, which already said
this and is now binding.

**The periodontal table control is a control, not a family.** See
`SF_CSA_RECORD_SELECTION.md` §1. The four families stay uniformly SCOP-derived.

### What was changed in the module

Three changes, `sf-csa`, all covered by `tests/test_rbh_provenance.py` (8 tests,
written red before the fix and failing for the three predicted reasons):

1. `classify_hit` takes `rbh` as an explicit keyword-only argument and no longer
   reads it from `target_meta`. Absence defaults to unpromoted.
2. The sequence leg — proteome universe, DIAMOND, reciprocal best hits — now runs
   *before* structural classification, and RBH is carried as
   `rbh_targets[query_accession] -> {target accessions}`. Pairwise by
   construction; the shared target row can no longer carry it.
3. `reject_reserved_fields` refuses any curator record declaring `rbh`,
   `orthology_status`, `identity_status`, `sequence_length`,
   `structure_residue_count`, `geometry` or `uniprot_annotation`. Applied to
   query records in `run_pipeline` and to campaign spec targets in
   `read_campaign_spec`. Assertion and measurement are now separable at read
   time, which is the property the false-positive gate depends on.

### A finding about the test suite, recorded rather than dropped

Two existing tests asserted the fail-open behaviour as correct
(`test_rbh_plus_whole_architecture_reaches_probable_same_function`,
`test_the_title_trap_is_a_release_audit_not_a_classifier_guard`), and both went
red on the fix. Ten further call sites passed `"rbh"` inside `target_meta` while
asserting *non*-promotion; those stayed green and would have continued to pass
while exercising nothing, since the classifier no longer reads the key. All
twelve were migrated to the real parameter. The suite is 119 passing.

This is the membrane failure mode inside the test suite: assertions that look
like coverage of a path and do not touch it. Worth stating because the panel's
own gates are only as good as the same distinction.

---

## Finding 4 — two search filters act in series, and only one is visible

Found 2026-09-01 by running the campaign, not by reading the manifest.

Under collection 2.5 settings, **eight of the sixteen declared judgments produced
no Foldseek row at all** — every `fold_analogy` and every `unrelated` pair. The
panel's own missing-evidence rule says a pair below the structural threshold still
yields a row saying so, and that an absent row means the release never compared
the pair. That distinction was unobservable, because absence covered both.

The cause is not the e-value. Foldseek prefilters before scoring, and the
prefilter ignores the e-value; the e-value then applies to whatever survives.
Measured for P00198 against the twelve-entry campaign database:

    -e 0.01                            2 rows
    -e 10000                           2 rows     prefilter, e-value irrelevant
    -e 0.01   --exhaustive-search      2 rows     no prefilter, e-value binds
    -e 10000  --exhaustive-search     12 rows

Neither setting alone reports a distant pair. Collection 2.6 sets both, and moves
filtering to `same_fold_tm` (0.5) and `whole_architecture_coverage` (0.7) — the
thresholds this panel actually declares — rather than leaving it to a search
heuristic that was never intended as a scientific gate. Both scientific thresholds
are unchanged.

### Result, against the frozen expectation

All sixteen judgments now produce a row, and all sixteen pass:

    exact                    4/4   exact_function_supported recovered
    homologous_superfamily   4/4   same_mechanism_class recovered
    fold_analogy             4/4   no promotion to a function claim
    unrelated                4/4   no promotion to a function claim
    analogy or unrelated promoted to a function claim: 0   (the bound allows 0)

A fresh re-run reproduces the previous run exactly: 144 rows compared, zero
TM-score differences. That is same-machine determinism only, and says nothing
about the cross-machine question.

Six of the sixteen are `unresolved_or_conflicted` rather than
`structural_analogy_only`, and that is the honest answer: 1CSP→1E6A at TM 0.314
and 1CGN→1C02 at TM 0.393 are genuinely below the declared `same_fold_tm` of 0.5.
SCOP calls those pairs same-fold; the geometry does not support it at this
threshold, and the module says so instead of splitting the difference.

## Withdrawn — a mis-scoring of the run above

**Claimed:** ten of sixteen judgments matched their stratum and six were
mismatches.

**Actual:** sixteen of sixteen pass. The comparison was made against an
expected-label table written from memory, which demanded `structural_analogy_only`
for `fold_analogy` and `unrelated`. The panel's frozen
`stratum_expected_label` requires **"any label outside `function_claim_labels`"**
for both — a false-positive bound, not an exact-label match — and
`unresolved_or_conflicted` satisfies it exactly as `structural_analogy_only` does.

The module was correct throughout; the scoring was not. Recorded because it is the
failure mode this programme exists to catch: a measurement built to fit the story
already in hand, and reported before it was checked against the frozen artifact.
It was caught by Yuvraj, not by a second measurement, which is the weaker of the
two ways to find it.

---

## Finding 5 — five queries can never have a sequence leg. The panel is blocked.

Found 2026-09-01 while recording expectations. The records are built and execute,
and they must not be frozen.

**Symptom.** Thirteen of sixteen records report `sequence_comparison` in
`missing_evidence`, against a declared expectation of nothing missing in any
stratum. Two of the thirteen are `exact` self-matches — a protein compared with
itself cannot lack sequence similarity, which is what showed this was a data
question rather than a biological one.

**Cause.** Six of the twelve queries do not appear in the proteome they declare:

    P00193  UP000192368   0 occurrences      P26394  UP000002695   0 occurrences
    P23370  UP000000532   0 occurrences      P00138  UP000595916   0 occurrences
    P45850  UP000011116   0 occurrences      P00147  UP000002361   0 occurrences

The files exist and load; they simply do not contain the protein. `run_pipeline`
raises only when the file is absent, so this passed silently.

**Why it cannot be fixed by pointing at a different file.** Checked against
UniProt: **five of the six belong to no reference proteome at all.**

    P00193  Peptoniphilus asaccharolyticus   none
    P23370  Thermus thermophilus             none
    P45850  Hordeum vulgare                  none
    P00138  Alcaligenes xylosoxydans         none
    P00147  Rhodobacter capsulatus           none
    P26394  Salmonella typhimurium LT2       UP000001014  -- wrong file declared

Only P26394 is a pointer error, and the correct proteome is the one
`SF_CSA_RECORD_SELECTION.md` already names for 1DZR. The other five are classic
Swiss-Prot entries whose organisms or strains carry no reference proteome, so no
file exists that would give them a sequence leg.

**What the selection gate actually checked.** The fourth gate was written as "the
organism must have a reference proteome", and it removed 2HMZ, 2MHR and 1JUH on
that basis. It verified that a proteome *exists*, never that the query is *in* it.
Five entries passed a gate that was measuring the wrong thing — the same shape of
error as the membrane panel, in the curation step rather than the gate.

**Scope of the damage.** Two of the four family anchor queries are affected:
P45850 (1FI2, ds beta-helix) and P00138 (1CGN, four-helical bundle). Each anchors
four judgments, so **eight of the sixteen records belong to families that cannot
produce the panel's second evidence leg.** The panel's own coverage rule requires
that "Foldseek and DIAMOND execute as separate checksum-bound evidence legs";
for half the panel, DIAMOND contributes nothing.

**This is a curation decision, not an engineering one.** Three ways out, and the
choice belongs to Yuvraj:

1. **Re-select the two affected families** with proteins that are in a reference
   proteome, adding "the query appears in its declared proteome" as a fifth,
   stated selection gate. Preserves the four-family rule and both legs.
2. **Accept single-leg records** for those families and declare it. Cheapest, and
   it weakens the claim the panel exists to support.
3. **Reduce the family count** — but the coverage rule requires exactly four
   unrelated families, so this is a rule change and a new collection version.

Recommendation is 1. Nothing should be frozen until it is settled: the recorded
expectations in `ADOPTION_DRAFT_SFCSA.json` are real measurements of a panel that
is half-evidenced, and adopting them would freeze that in.

---

## Finding 6 — the reselection fixed the data, and exposed a semantics error

The twelve queries were reselected on 2026-09-02 with the missing fifth gate
applied: the accession must appear in the proteome it declares, not merely belong
to an organism that has one. All twelve pass all five gates, and the structural
side is clean:

    located 16/16 | passing the frozen expectation 16/16 | promotions 0 (bound 0)

**The data-integrity failure is gone.** Every query now finds itself in its own
proteome as a reciprocal best hit, and the queries that previously produced
nothing now produce real homology: P00198 finds a ferredoxin at 55.9 percent
identity, P26394 finds three homologs between 57 and 65 percent.

**Eleven of sixteen judgments still carry no sequence-leg row, and that is a
different fact from before.** It is now the declared search reporting no homology,
not a broken input. Measured per pair:

- `fold_analogy` and `unrelated` pairs have no sequence homology. That is the
  point of those strata; a row would be the surprise.
- Of the three `homologous_superfamily` pairs, two are below the declared
  thresholds and are correctly absent: 1DZR→1J1L at 24.5 percent identity over
  29.0 percent of the query, and 1C02→2OOC at 34.8 percent over 13.8 percent,
  against declared minima of 20 percent identity and 40 percent coverage.
- 1FDN→1BLU is the interesting one. A local alignment scores it 50.9 percent
  identity at 98.2 percent coverage, which looks like it should be found — and
  **DIAMOND does not find it at any e-value.** Rerun at `--evalue 10` against the
  same database it still returns two hits, neither of them P00208. The seed
  filter does not reach this pair. The local-alignment number was the misleading
  one, not the tool.

### What this means for the panel, and it is a decision

`expected_missing_by_stratum` declares that nothing may be missing in any
stratum. That is false as written, and it is the sequence-side twin of Finding 4:
the module treats an absent row as missing evidence, which conflates *searched
and found nothing* with *never searched*. On the structural side that was fixed by
making the search exhaustive, so every declared pair yields a row saying where it
fell. **There is no equivalent for DIAMOND** — a sequence aligner does not emit
non-hits, and forcing it to would not mean anything.

So the fix has to be in the semantics rather than the search:

1. **Declare the sequence leg legitimately absent for `fold_analogy` and
   `unrelated`.** No sequence homology is the expected finding for those strata.
2. **Record it per record for `homologous_superfamily`**, with the measured
   identity and coverage, so a reader sees *why* a homolog produced no sequence
   row rather than inferring a fault.
3. Keep the coverage rule that both legs execute — DIAMOND does execute for every
   query. What varies is whether it finds anything, and that is a result.

This changes what a green panel asserts, so it belongs to Yuvraj rather than to
whoever is next at the keyboard. Nothing is frozen until it is settled.

---

## Records built and falsified — 2026-09-02, collection 2.7

`ADOPTION_DRAFT_SFCSA.json`: sixteen records, four families by four strata,
expectations recorded by execution. **Not merged into `PANEL_MANIFEST.json`** —
adoption needs independent reproduction on a second machine.

    16/16 cases passed | 8/8 coverage features witnessed | promotions 0 (bound 0)
    strata executed: exact, fold_analogy, homologous_superfamily, unrelated

### Protocol rule 2, twice

**A — an exception stripped of its justification.** A record legitimately
declaring `expected_missing: [sequence_comparison]` had its measurement deleted
on a copy:

    15/16 passed; v2-sfcsa-homologous_superfamily-ferredoxin_like-1FDN-1BLU failed
      record_expected_missing_is_justified            observed None
      missing_evidence_matches_stratum_expectation    expected [], observed [sequence_comparison]

It fails **twice**: once for claiming an exception without evidence, and again
because the strict stratum expectation snaps back when the justification is gone.
That is the property the mechanism needed. Had it failed only the first check, or
passed, the per-record exception added the same day would have been a hole in the
gate rather than a way to state a measured fact.

**B — an unrelated pair relabelled a homolog.** 1CSP→1DZR, which has no homology,
recurated as `homologous_superfamily` on a copy:

    15/16 passed; v2-sfcsa-unrelated-ob_fold-1CSP-1DZR failed
      homologous_superfamily_control_recovered   expected same_mechanism_class,
                                                 observed unresolved_or_conflicted

The recall gate catches a curator claiming ancestry the evidence does not support,
which is the direction that would flatter the panel.

### What a green panel here does and does not assert

It asserts that four exact self-matches are recovered, that four homologous
superfamily pairs are recovered as `same_mechanism_class`, and that **no fold
analogy and no unrelated pair was promoted to a functional claim** — the absolute
bound, held at zero.

It does not assert that both evidence legs support every judgment. No
`homologous_superfamily` record in this selection carries a DIAMOND row, so that
stratum rests on structural evidence alone, and `probable_same_function` is
unreachable in this campaign at all. Both are recorded in
`gate_semantics.sf_csa.labels_not_exercised` and
`.sequence_leg_measurements` rather than left for a reader to infer from a
passing summary line.

---

## Finding 7 — the panel's load-bearing gate cannot fail

Found 2026-09-02 by running a held-out test on proteins the module had never seen,
at Yuvraj's instruction. The test passed. What it exposed is the finding.

### The held-out test

SCOP fold `d.15` (beta-Grasp), untouched by the four panel families. Four
judgments, predictions written down in `SF_CSA_HELDOUT_PREREGISTRATION.md` before
the data was acquired, compared once:

    exact                   1DOI->1DOI  exact_function_supported    PASS
    homologous_superfamily  1DOI->1DOX  same_mechanism_class        PASS
    fold_analogy            1DOI->1EO6  not a function claim        PASS
    unrelated               1DOI->1CSP  not a function claim        PASS

    4/4 predictions met, 0 rows missing, 0 promotions

The recall gates generalise: an unseen homolog pair in an unseen fold is recovered
as `same_mechanism_class` at 0.915 query coverage, and an unseen self-match is
recovered exactly.

### What passing revealed

Predictions 3 and 4 passed, and **they could not have failed.**

`classify_hit` promotes only when `qgroup == tgroup`. The curated mechanism group
*is* the SCOP superfamily — that is the design decision recorded in
`sf_csa_build_campaign.py`, made so that homologs share a group and analogues do
not. But `fold_analogy` and `unrelated` are *defined* as different superfamily.
So for every record in those strata `qgroup != tgroup` holds by construction, and
promotion is unreachable at any structural score.

Demonstrated directly rather than argued:

    analogy pair, perfect whole_architecture_match, rbh=False -> structural_analogy_only
    analogy pair, perfect whole_architecture_match, rbh=True  -> structural_analogy_only
    same-group pair, identical structural input,    rbh=True  -> probable_same_function

The structure is as promotable as it can be made. The label does not move.

### Consequence

`analogy_or_unrelated_promoted_to_function_max = 0` is described in the panel's
own gate semantics as **"the scientifically load-bearing one … an absolute bound,
not a tolerance"**. It is satisfied by the definition of the strata, not by
anything the module computed. By protocol rule 2 — *a gate that cannot fail is not
a gate* — it is not currently a gate.

This is Finding 2's shape returning at the centre of the panel. There, an
unoverridden mechanism table made the same bound pass vacuously while recall
failed 16/16. The tables were overridden and recall was fixed; the bound stayed
vacuous for a different reason, and the fix hid it rather than removing it.

### What would make it a real gate

Not a threshold change. The bound needs a case where the module could plausibly
promote and must not, which means a pair whose *target* group is inferred rather
than curated. That is the query-versus-experimental-PDB axis, which
`gate_semantics.sf_csa.comparison_axis` currently declares out of scope: off the
campaign axis the target's group comes from the mechanism-family regex over a PDB
title, and a title that matches the query's family would produce `qgroup ==
tgroup` on a pair with no curated relationship. That is exactly the promotion the
bound claims to forbid, and it is the only route by which it can occur.

Whether to bring that axis into scope is a scientific decision and Yuvraj's.
Until then the gate should be **relabelled from scientific to definitional**, so a
green panel is not read as three scientific results plus a contract check when one
of the three cannot fail.

---

## Finding 8 — a sha256 on a live query is not a lock

Found 2026-09-04 by pushing the sf-csa CI enablement and reading why all six
runners failed. The failing step was `Acquire the locked artifacts`, before any
panel executed.

### What broke

`SOURCE_LOCK.json` gained 48 sf-csa entries. Ten of them are reference proteomes
recorded as UniProt *stream queries*:

    https://rest.uniprot.org/uniprotkb/stream?query=%28proteome%3AUP000005640%29&format=fasta

That endpoint has no version parameter. It returns whatever the current UniProt
release holds. Measured against the locked digest on the day of the push:

    locked   eb9a0f1cd363c8e7afa20517c46c7ae517ff80564cbfbe4d73db956f78118f7d   144,818 sequences
    fetched  4477c2a858646a958fc84fca351b21a4f479ce4744cf3c054bc889fa69759369   147,520 sequences

`acquire_sources.py` retries "until its digest matches", which for this endpoint
is never, so every job spent its retries and failed. Fifteen minutes, six
runners, no science.

### Why it is a finding and not a CI defect

`SOURCE_LOCK.json`'s stated purpose is that a reviewer "acquires the exact same
bytes themselves". For these ten entries no reviewer ever can, including the
author on a second machine. The lock recorded a digest for a file the URL cannot
return again.

That is worse than an unreproducible build step, because the entry's own note
says the proteome set defines the measurement:

> Part of the declared search space: coverage_boundary states the search is
> exhaustive only within these proteomes, so a different file is a different
> measurement.

Taken together: the sf-csa sequence leg is computed against a search space that
silently changes with every UniProt release, and nothing in the panel could
detect it. The digest was the thing that was supposed to detect it.

### What was done, and what was not

The ten proteome entries are **withdrawn from the lock** (276 -> 266) and the
sf-csa CI steps are removed, restoring the workflow byte-identically to the last
green run. The other 38 sf-csa entries — AlphaFold models, per-accession UniProt
FASTAs, RCSB coordinate files — are addressed by stable URLs, verified reachable,
and stay locked.

Nothing was repaired, because the repair is a scientific choice and is Yuvraj's:

1.  **Pin to a UniProt release.** `previous_releases/release-2026_02/` and its
    siblings are versioned and stable. This makes the lock a real lock, but the
    files it pins are a *different* release from the ones the sixteen sf-csa
    expectations were recorded against, so those expectations must be re-recorded
    and the sequence leg re-measured.
2.  **Declare the search space by accession list** rather than by live query, and
    lock the accessions instead of the bytes.
3.  **Narrow the search space** to proteomes small enough to redistribute, which
    changes what "exhaustive" covers.

Option 1 is the smallest change that makes the claim true. It is not free: 153 MB
per runner, and re-recorded expectations.

### One thing the failed run did establish

Every step before acquisition passed on all six runners, including the pinned
install and the version assertion. Foldseek 10.941cd33 and DIAMOND 2.1.11 install
reproducibly from conda-forge/bioconda on ubuntu-latest and macos-latest across
Python 3.10, 3.11 and 3.12. That part of the enablement is measured and can be
re-landed as-is.
