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
