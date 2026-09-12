## 2026-09-11 — absent density may not improve a geometry call

Behaviour change in `geometry_signal` and `assign_label`, recorded here because it
changes what a label means, not only how it is computed. It continues the
2026-09-02 entry below: that one stopped `completeness_signal` claiming more than a
broad residue set can support, and this one stops `geometry_signal` claiming
anything at all about a site it only partly observed.

`geometry_signal` bailed out only when fewer than two declared catalytic positions
were resolved. With three of four present it measured the cluster over those three
and returned `supported`. The screen was wrong in both directions:

- **It over-claimed, and did so more the less it saw.** Dropping a declared position
  that has no density also drops it from the clustering test, so a site scored
  *better* the more of it was missing. On human PRPS1 entries 8DBG, 8DBI and 8DBL
  the declared set is 130, 196, 220 and 221; R196 has no density in any chain, and
  in every entry that does resolve it, it sits about 20 A from D221. Removing it
  turned a 20.64 A failure into a 10.48 A pass, and the resulting rationale read
  "the site is intact" about coordinates with nothing at the catalytic residue. For
  a protein whose signature is a mobile catalytic loop that is usually unmodelled,
  that inverts the answer: every entry that resolved the loop was labelled
  `inactive_conformation` for showing it retracted, and the three that resolved
  nothing were labelled `probable_active`.
- **The cap that limited the damage was an accident.** `probable_active` was not a
  ceiling those three had earned; it was where they landed because `conformation`,
  `occupancy` and `assembly` all happened to be unavailable. Supplying a
  conformational reference set removes one of those reasons, so the next caller to
  build one would have lifted all three to `active_state_supported`, the strongest
  label in the vocabulary, on coordinates with no density at the catalytic residue.

A third gap sat alongside it and is now closed by the same coverage field:
`completeness_signal` reads the **sequence** and reported "all 4 annotated catalytic
position(s) hold competent residues (130H, 196R, 220D, 221D)", while
`geometry_signal` reads the **coordinates**. Nothing reconciled them, so one signal
could assert presence while the other was measuring absence.

Changes:

- `SIGNAL_STATES` gains **`unevaluable`**. It is distinct from `unavailable` on
  purpose: both mean "not judged", but `unavailable` means no structure was
  supplied, and `unevaluable` means the structure was supplied and the declared set
  was incomplete. Collapsing them would hide exactly the case this entry is about.
- `geometry_signal` returns `unevaluable` whenever **any** declared catalytic
  position is unresolved, rather than only when fewer than two survive. It no longer
  computes a separation over a partial set.
- Every geometry signal now carries `declared_set_coverage` as `{declared, resolved}`
  alongside `missing_from_structure`, so a reader can see how much of the site was
  observed without reconstructing it from two other fields.
- `assign_label` returns `indeterminate` on an `unevaluable` geometry, and returns
  rather than collecting a cap, so no downstream signal can lift it to a positive
  label. The `supported` detail string no longer carries an "N position(s)
  unresolved" clause, because that branch is now only reachable with a complete set.
- `SCHEMA_VERSION` in `io.py` 1.0 -> 1.1: the state vocabulary is part of the
  published contract.

Verified on the eight-state human PRPS1 ladder that exposed it. 8DBG, 8DBI and 8DBL
move `probable_active` -> **`indeterminate`** with geometry `unevaluable` at 3/4
coverage; 2H06, 8DBE, 8DBK, 8DBN and 8DBO are unchanged at `inactive_conformation`
with geometry `contradicted` at 4/4. No other label moved.

Tests: 146 -> 154. One migration in `test_core.py`, where
`test_unresolved_positions_block_geometry_rather_than_biasing_it` already held this
principle for the two-of-three case and now expects `unevaluable`. Eight new tests in
`test_geometry_declared_set_coverage.py`, including the one that matters — a pair of
resolved residues that cluster while a third declared position is simply absent must
not produce a verdict — and a guard that a predicted-active reference set plus an
active fold state still cannot lift partial coverage off `indeterminate`.

No qualification panel runs `actstate`, so no adopted baseline changes. This does
advance the "resolve the ActState occupancy caveat" item on the pre-release checklist.

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
