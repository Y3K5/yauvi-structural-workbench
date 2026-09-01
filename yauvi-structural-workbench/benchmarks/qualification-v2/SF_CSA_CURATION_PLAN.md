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

The frozen `requirements` block covers only the 16 judgments. Composition counts
requirements by exact equality on `(record_kind, stratum, split)`, so the two
controls need their own requirement entries or they will execute without being
counted. Add them in the same edit that adds the records.

## 2. Every pair is query-to-query

Fixed in `gate_semantics.comparison_axis`, and it constrains the campaign spec:
every judged pair is between two campaign queries, so both sides carry a curated
`mechanism_group`. Off that axis the target's group comes from a regex over the
PDB title and the panel would be gating free text.

Consequence: the campaign needs roughly 10–12 queries to yield 16 pairs — each
family contributes a reference member plus the partners its four strata need,
with `unrelated` partners drawn from the other families rather than added.

## 3. Families

| # | family | default tables cover it | role |
|---|---|---|---|
| 1 | serine proteases (chymotrypsin clan) | no | classic |
| 2 | TIM-barrel enzymes | no | classic |
| 3 | immunoglobulin-like β-sandwich | no | classic |
| 4 | `omp85_bama` | **yes** | table control |

Family 4 earns its slot: it is the one family the periodontal defaults already
classify. It must produce the same judgment with the defaults and with the
override. That is the evidence the override changed the biology it was meant to
change and nothing else — otherwise "we overrode the tables" is an assertion
about a file, not a measured fact.

## 4. What is NOT frozen here, and why

**The specific protein pairs.** Each `fold_analogy` record requires a pair whose
folds match above the release threshold and which the literature holds to be
non-homologous. That is the hardest judgment in the panel and the easiest to get
wrong in a way that flatters the result: picking a "fold analogy" that is
actually a distant homolog turns a passing false-positive gate into a
meaningless one, and nothing downstream would catch it.

Two rules, fixed now:

- **Homology comes from a cited source, never from the TM-score the panel is
  testing.** Using the structural similarity to decide the stratum, then testing
  whether the module agrees, is circular and would pass by construction.
- **Each record carries its citation and the boundary of what that source
  claims.** A source asserting "no detectable sequence similarity" is not a
  source asserting non-homology.

Candidate directions, all requiring verification before selection: the
`(β/α)8` barrel is the standard non-homologous-same-fold case and is where the
`fold_analogy` stratum is most defensible; the immunoglobulin-like sandwich is
shared across domains with independent origins argued in the literature. The
chymotrypsin/subtilisin pair, the textbook convergence example, is convergent at
the *catalytic triad* and not at the fold — the two are α/β and β-barrel — so it
would sit below the structural threshold and belongs in `unrelated`, not
`fold_analogy`. That distinction is exactly the kind of thing to check against
sources rather than recall.

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
