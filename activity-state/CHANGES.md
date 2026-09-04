## 2026-09-02 — a broad residue set may cap a claim, not make one

Behaviour change in `completeness_signal` and `assign_label`, recorded here
because it changes what a label means, not only how it is computed.

`active_site_disrupted` — the strongest negative label in the closed vocabulary
— was emitted whenever a residue at an annotated `ACT_SITE` position fell
outside `CATALYTICALLY_COMPETENT`, a 13-letter set covering nucleophiles,
acid/base pairs, metal ligands and the two residues that act through backbone
geometry. The set knows nothing about the role the position plays and nothing
about which residue the annotation expects there, so the screen was wrong in
both directions:

- **It over-claimed.** Seven residues fall outside the set. Any of them at an
  annotated position produced the disruption label, without a role, without an
  expected residue, and without comparison against a validated ortholog. A
  sequence that is the wrong isoform, or numbered against a different entry, is
  indistinguishable from a pseudoenzyme by that test.
- **It was silent on the commoner degradation.** A catalytic Cys to Ser, or His
  to Asn, stays inside the set, so a genuinely dead site read as `supported`.

Four changes:

- `completeness_signal(record, features, *, expected_residues=None)` — the
  expectation is explicit and keyword-only, for the reason `sf_csa.classify_hit`
  takes `rbh` that way: a field read out of a record the caller also controls is
  not provenance. It maps an annotated position to the residue an
  experimentally validated reference carries there.
- The position-specific comparison runs first and is decisive. It is the only
  evidence that can establish disruption, and the only one that catches a
  within-set substitution. A position whose expectation matches is settled
  whatever the competence set thinks of that residue — the reference is the
  authority on its own site.
- Without an expectation, a non-competent residue is still reported
  `contradicted`, and `assign_label` caps at `indeterminate` rather than
  inverting: the observation is kept in the signal, the rationale names the
  positions and what would raise it, and the branch returns immediately so no
  downstream signal can lift it into a positive claim. Capping is not
  fail-open — `active_state_supported`, `probable_active` and
  `apo_but_competent` all stay unreachable on that path.
- `normalize_expected_residues` rejects an expectation the pipeline cannot
  check: a position carrying no `ACT_SITE` annotation, or anything that is not a
  standard one-letter residue code. Rejections are named in the signal, in the
  run config as `rejected_expectations`, and printed — a mistyped entry doing
  nothing silently is how a curator concludes it was applied. `assess` stays
  total; validation never raises.

`actstate run` gains `--expected-residues`, and the bundled `P_DISRUPTED`
fixture now ships `examples/expected_residues.json`, so the fixtures still
exercise every label in the vocabulary and the example demonstrates the path
that reaches the strong one.

Covered by `tests/test_site_disruption_provenance.py`, written red first.
Four existing tests asserted the old behaviour and were migrated: two in
`test_core.py` and `test_properties.py` now pass an expectation and are joined
by their capped counterparts, the golden fixture run supplies the sidecar, and
property P3 is restated — a degraded position is still decisive, and which label
it is decisive *for* now depends on how the degradation was established. Suite:
146 passing, up from 127.

This closes the ActState half of the two interpretation defects the pre-public
audit recorded; the SF-CSA half was closed on 2026-09-01. The occupancy caveat
is untouched: a non-solvent heteroatom is detected, but its identity is still
not proven against the declared cofactor.
