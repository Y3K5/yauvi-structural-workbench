
## 2026-09-01 — reciprocal best hit is evidence, not a manifest field

Behaviour change in `classify_hit`, made before sf-csa panel curation and
recorded here because it changes what a release means, not only how it is built.

`probable_same_function` — the strongest non-exact label in the closed vocabulary
— could be switched on by a curator writing `"rbh": true` into a query record.
The key survived manifest validation, survived `row = {**q, ...}` into
`target_meta`, and reached the classifier, which emitted the label with the basis
string "RBH plus compatible whole architecture and mechanism annotation" for a
pair on which no reciprocal best hit had been computed. The computed RBH ran
after classification and wrote elsewhere, so the label was never reached by
measurement. And `target_meta` is keyed by target accession alone, so the flag
was query-independent even in principle.

Three changes:

- `classify_hit(..., *, rbh: bool = False)` — explicit, keyword-only, defaulting
  to unpromoted. It no longer reads `rbh` from `target_meta`.
- The sequence leg (proteome universe, DIAMOND, reciprocal best hits) runs before
  structural classification. RBH is carried as `rbh_targets[query] -> {targets}`,
  pairwise by construction.
- `reject_reserved_fields` refuses curator records declaring any field the
  pipeline computes (`rbh`, `orthology_status`, `identity_status`,
  `sequence_length`, `structure_residue_count`, `geometry`,
  `uniprot_annotation`), in both `run_pipeline` and `read_campaign_spec`.

Covered by `tests/test_rbh_provenance.py`, written red first. Two existing tests
asserted the old fail-open behaviour and were corrected; ten further call sites
were passing `rbh` through `target_meta` while asserting non-promotion and would
have gone on passing without exercising the path. All migrated. Suite: 119 passing.

Full record, including the two withdrawn statements this supersedes:
`yauvi-structural-workbench/benchmarks/qualification-v2/SF_CSA_PREADOPTION_FINDINGS.md`.
