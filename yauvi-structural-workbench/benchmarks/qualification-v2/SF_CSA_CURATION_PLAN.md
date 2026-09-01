# sf-csa curation plan — structure frozen, biology not yet

Written after the three §6 decisions were settled and before any record exists.
What is frozen here is the *shape* of the 18 records and the rule each must
satisfy. The specific protein pairs are candidates, not selections: assigning
homology is a literature judgment, and freezing it from memory is the failure
this project exists to avoid.

Decisions in force (2026-08-31):

1. `probable_same_function` — **Option C**, gated as a defect with two controls.
2. Families — **three classic plus one periodontal table control.**
3. `structural_and_sequence_outputs` — **labelled a contract check.**

All three, with their rationale, are frozen in the panel's `gate_semantics`.

---

## 1. Record inventory: 18, not 16

| kind | stratum | count | notes |
|---|---|---|---|
| `relationship_judgment` | `exact` | 4 | one per family |
| `relationship_judgment` | `homologous_superfamily` | 4 | one per family |
| `relationship_judgment` | `fold_analogy` | 4 | one per family |
| `relationship_judgment` | `unrelated` | 4 | one per family |
| `control_case` | `homologous_superfamily` | 1 | `rbh_computed_path` |
| `control_case` | `homologous_superfamily` | 1 | `rbh_asserted_rejected` |

**Corrected 2026-08-31.** An earlier draft of this plan said the two controls
needed their own `requirements` entries. They do not. Controls live in the
panel's separate `controls` key, not in `records`: the runner executes them in
their own pass, counts them under `control_counts`, and its exit code requires
`controls_passed == len(controls)` independently of the case count. The
`requirements` block governs `records` only. The precedent is
`ADOPTION_DRAFT_XRAY.json`, which carries two controls (`fail_closed` and
`coverage`) declared exactly this way and no control requirement entries.

One interaction to check at curation rather than assume. Beyond the sf-csa gates,
`run_case` applies a shared control rule: `control_purpose == "fail_closed"`
demands a non-empty missing-evidence list, and **any other purpose** is held to
`expected_missing_by_stratum`. Both sf-csa controls use neither of the two
established purposes, so they fall into the second branch and must report no
missing evidence — which is right, since both are meant to produce complete
structural and sequence legs and differ only in what the classifier does with
them. That path is not exercised by `sf_csa_gate_falsification.py`, which tests
the sf-csa gate function directly and does not reach the shared block.

## 2. Every pair is query-to-query

Fixed in `gate_semantics.comparison_axis`, and it constrains the campaign spec:
every judged pair is between two campaign queries, so both sides carry a curated
`mechanism_group`. Off that axis the target's group comes from a regex over the
PDB title and the panel would be gating free text.

Consequence: the campaign needs roughly 10–12 queries to yield 16 pairs — each
family contributes a reference member plus the partners its four strata need,
with `unrelated` partners drawn from the other families rather than added.

## 3. Families — revised against sources, 2026-08-31

The first draft of this plan proposed serine proteases, TIM barrels and
immunoglobulin-like sandwiches as the three classic families. **Two of those three
are withdrawn.** They were proposed from recall and the literature does not
support them for the `fold_analogy` stratum:

- **TIM barrels — withdrawn.** Copley and Bork found statistically reliable
  sequence evidence that at least 12 of 23 SCOP (βα)8 superfamilies share a
  common origin, covering all but one of the barrels in central metabolism.
  Nagano and colleagues describe 21 homologous superfamilies whose active sites
  all cluster at one end of the barrel. Two TIM barrels are therefore a probable
  distant-homology pair, not an analogy pair. Curating them as `fold_analogy`
  would put homologs in the stratum whose whole purpose is to contain
  non-homologs — and the false-positive gate would then pass for the wrong
  reason, which is the failure this panel exists to detect.
- **Immunoglobulin-like sandwich — withdrawn.** The recent Ig-fold literature
  argues the opposite of independent origin: a shared "irreducible structural
  signature" across all Ig-fold variants, and unsuspected structural homologies
  between folds classified separately. The independent-origin claim rests on a
  single low-citation hypothesis paper.
- **Serine proteases — retained, but not as an analogy source.** Chymotrypsin and
  subtilisin are the textbook convergence case at the *catalytic triad*, not at
  the fold: α/β versus β-barrel. The pair sits below the structural threshold and
  belongs in `unrelated`.

### The operational criterion

Fixed, and taken from the source rather than from judgement:

> **analogues** — same SCOP fold, different SCOP superfamily: similar
> three-dimensional structure with little evidence of a common ancestor.
> **homologues** — same superfamily.

That is Russell, Saqi, Sayle, Bates and Sternberg (J Mol Biol, 1997), and it is
the definition the `fold_analogy` and `homologous_superfamily` strata now use. It
is citable, checkable per record, and it does not require the curator to settle
what actually happened in evolution — only what the classification asserts.

A record additionally fails curation if there is published sequence-level
evidence of common origin for the specific pair, which is what disqualified the
TIM barrels above.

### Revised families

Russell and colleagues (J Mol Biol, 1998) identify superfolds that contain
genuinely analogous members and — usefully — exclude the TIM barrel, whose
supersite they place among *homologous* proteins. The first three families are
taken from that list:

| # | family | fold_analogy source | role |
|---|---|---|---|
| 1 | ferredoxin-like fold | named superfold with analogous members | classic |
| 2 | four-helical bundle | named superfold with analogous members | classic |
| 3 | double-stranded β-helix | named superfold with analogous members | classic |
| 4 | `omp85_bama` | n/a | table control |

Family 4 is unchanged and still earns its slot: it is the one family the
periodontal defaults already classify, so it must produce the same judgment with
the defaults and with the override. That is the evidence the override changed the
biology it was meant to change and nothing else.

## 4. How rare analogy is, and what that costs this panel

Wright (Genome Biol Evol, 2025) is directly on point, because it asks this
question with the tool this panel pins. Across Foldseek clusters, only about
2.6% lack sequence-level support for homology, and only about 1% of strong
matches at TM-score ≥ 0.5. The paper's conclusion is that caution is warranted
when inferring homology from structural resemblance alone.

Two consequences, both binding:

1. **A casually chosen "fold analogy" pair is far more likely to be a distant
   homolog than a true analogue.** Four such records cannot be assembled by
   picking structurally similar proteins with different names.
2. **Homology must never be inferred from the TM-score the panel is testing.**
   Using structural similarity to assign the stratum, then testing whether the
   module agrees, is circular and passes by construction. The stratum comes from
   the SCOP classification and the literature; the TM-score is what is under
   test.

Each record carries its citation and the boundary of what that source claims. A
source asserting "no detectable sequence similarity" is not a source asserting
non-homology — that distinction is exactly what the Wright result turns on.

## 5. Data acquisition

Reachable from here, verified: `files.rcsb.org` and `rest.uniprot.org` both
return 200. Needed and not yet present under `sources/`:

- per-family structures (queries are predicted or experimental monomers whose
  sequence must match their FASTA exactly — `run_pipeline` refuses otherwise),
- a source proteome FASTA per query, for the RBH leg,
- a Foldseek-formatted structure database, plus its version and checksum, which
  `build-manifests` requires and refuses to fabricate.

`sources/` currently holds wwpdb (119), uniprot (39), alphafold (8), opm (23) and
no proteomes and no structure database. Acquisition is a real step, not a
formality, and `acquire_sources.py` is where it belongs.

## 6. Order of work

1. Verify candidate pairs against sources; freeze the 18 records with citations.
2. Add the two `control_case` requirement entries.
3. Acquire sources; extend `SOURCE_LOCK.json`.
4. Build the campaign spec, run `build-manifests`, and freeze
   `gate_semantics.sf_csa.database_manifest_sha256` — the gate fails closed
   until this is set, deliberately.
5. Run the table control both ways (defaults vs override) and record both.
6. Execute the panel; record expectations by measurement, never by hand.
7. Reproduce on a second machine before calling the scope adopted.

Steps 1 and 3 are the long ones. Nothing before step 7 qualifies the scope.
