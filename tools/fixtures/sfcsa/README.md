# SF-CSA offline golden fixture

An end-to-end exercise of the SF-CSA pipeline with no foldseek, no diamond, no
reference PDB database and no network. It produces two releases: one that passes its own audit and one that fails it. The
driver takes about a second; the pytest form (21 tests, which run each scenario
more than once) takes about three.

```
PY=/path/to/python ./run_fixture.sh                 # run both scenarios, diff against golden/
PY=/path/to/python ./run_fixture.sh --update-golden # accept the current output as golden
/path/to/python -m pytest tools/fixtures/sfcsa/ -q  # the same thing as assertions
```

`PY` must point at an interpreter with `sf_csa` importable; it defaults to
`python3`. From the repo root, `PY=$PWD/.venv/bin/python` works.

## Why this exists

Before this fixture, `sf-csa`'s 110 tests all called functions. Nothing ran
`run_pipeline`, because running it meant installing two structural-biology
binaries and a copy of the PDB. The consequence was not a coverage gap in the
abstract: the pipeline's shape — which files a release contains, whether the
audit can fail, whether the classification of a hit survives the trip through
Foldseek's positional output columns and back — was untested. Three of the
findings in `REVIEW.md` are things you cannot see by reading `classify_hit` in
isolation.

## How it stands in for the tools

The pipeline shells out. `foldseek_search` and `diamond_search` build argv
lists and hand them to `run_cmd`, which is `subprocess.run`, and the
executables are found with `shutil.which`. So the fixture substitutes at the
process boundary: `stub_bin/foldseek` and `stub_bin/diamond` are executable
Python scripts, and `run_fixture.sh` prepends `stub_bin/` to `PATH`.

Substituting there rather than monkey-patching Python internals means the
argv construction, the `--format-output` field ordering, the TSV parsing, the
version gate and the temp-directory handling are all still under test. A
monkey-patched `foldseek_search` would skip every one of them.

The stubs compute nothing. Each number in `stub_bin/hits.json` was written by
hand to sit on a chosen side of a threshold declared in the fixture's
`database_manifest.json`. What the fixture tests is the pipeline's
*interpretation* of alignment output — not the alignment.

## The two scenarios

| | `main` | `trap` |
|---|---|---|
| manifests | `query_manifest.json`, `database_manifest.json` | `query_manifest_trap.json`, `database_manifest_trap.json` |
| queries | QRY_A, QRY_B | QRY_A, QRY_B (QRY_B declares `rbh: true`) |
| audit | must **pass** | must **fail** |
| golden | `golden/main/` | `golden/trap/` |

A fixture whose audit always returns clean does not demonstrate that the audit
works. `trap` exists so that `verify_release` has something to catch, and
`run_fixture.sh` treats a clean `trap` audit as a broken fixture, not a pass.

Both scenarios read the same `hits.json`. A hit carrying `"_scenarios":
["trap"]` belongs to that scenario alone; a hit with no `_scenarios` key
belongs to both. The stub reads `SFCSA_FIXTURE_SCENARIO` from the environment.

## What each hit demonstrates

Thresholds in force: `same_fold_tm` 0.5, `whole_architecture_coverage` 0.8,
`sequence_min_identity` 30, `sequence_min_query_coverage` 60.

**QRY_A vs QRY_A** (`campaign_models`) — the self-match control. TM 1.00.
`classify_hit` short-circuits an exact accession match to
`exact_function_supported` before any threshold is consulted, and records
`"self-match is a control, not an independent function experiment"`. This is
the only route to that label in the fixture.

**QRY_A vs SYN_WHOLE** (`experimental_pdb`) — TM 0.86, coverage 0.91/0.89:
above `same_fold_tm` and above `whole_architecture_coverage` on both sides, so
`whole_architecture_match`. Classification stops at `same_mechanism_class`.

**QRY_A vs SYN_PART** — TM 0.62, coverage 0.44/0.51: above the fold threshold,
below the coverage threshold, so `domain_or_partial_match` →
`structural_analogy_only`. Pushing the TM below 0.5 moves this hit to a third
category, which is how the fixture pins that boundary from both sides.

**QRY_A vs SYN_TRAP** — TM 0.88, coverage 0.93/0.90, title `"... toluene
transporter, unrelated substrate"`. The structurally strongest experimental hit
in the fixture, and it carries the trap substring — yet the `main` scenario
still audits clean, because the hit only reaches `same_mechanism_class`, which
is not in the trap's `must_not_promote_to` list. That is the correct behaviour
and worth having a case for: the trap fires on the *label*, not on the title. A
trap that failed any release mentioning "toluene" would be unusable.

**QRY_B vs SYN_BELOW** — TM 0.31, below `same_fold_tm`, so
`below_structural_similarity_threshold` → `unresolved_or_conflicted` with
`"no same-fold interpretation allowed"`. QRY_B's only experimental hit, so
QRY_B's dossier is the "nothing can be said" case, which is the one most worth
having a golden copy of.

**QRY_A vs QRY_B** (`campaign_models`) — TM 0.83, coverage 0.87/0.94, same
`mechanism_group`. In `main` this reaches `same_mechanism_class`. In `trap`,
where QRY_B declares `rbh: true`, the identical geometry reaches
`probable_same_function` — and because the trap variant's title carries
"toluene", `verify_release` catches it. The pair does double duty: it is both
the audit's demonstration and the fixture's record of the defect below.

Sequence side: `QRY_A → SYN011` and `QRY_B → SYN029` are reciprocal;
`QRY_A → SYN041` is deliberately not (`reverse_best_hit` points it at
`SOMETHING_ELSE`), so both orthology statuses appear in the golden output.

## The RBH promotion gate is dead

Found while building this fixture, and the reason the `trap` scenario has to
declare `rbh` by hand.

`classify_hit` gates promotion to `probable_same_function` on a truthy `rbh`
key in the target's metadata:

```python
if target_meta and target_meta.get("rbh") and category == "whole_architecture_match":
    return ("probable_same_function", ...)
```

Nothing in the pipeline ever writes that key. `run_pipeline` builds target
metadata as `target_meta = {q["accession"]: q for q in validated}` — purely
from validated query manifest entries. `reciprocal_best_hits` does compute
reciprocity, and its result is written into the *sequence* rows as
`orthology_status`, but it is never fed back to the structural
classification. So the second-strongest label in a six-label vocabulary is
unreachable on any real input: a hit can only be promoted if a human hand-wrote
`"rbh": true` into the query manifest.

This is the same failure family as the assembly-signal defect in `actstate`
(`REVIEW.md` §4) — a signal computed and then not consulted — but it fails in
the opposite direction. There it inflated a label; here it suppresses one. A
campaign comparing two of its own structures, whole-architecture match, same
mechanism group, reciprocal best hits in both proteomes, still reports
`same_mechanism_class`.

`test_the_main_scenario_cannot_reach_probable_same_function` asserts the
current behaviour. When the defect is fixed that test fails, which is the
intended signal — the fix is a behaviour change and the golden trees change
with it.

## Golden comparison and what varies

A release is not byte-comparable across machines. Exactly three things vary,
and `canonicalise.py` handles them rather than leaving the golden trees to rot:

1. `SF_CSA_RELEASE_MANIFEST.json` records `release_id`, taken from the basename
   of `--output`. Replaced with `<RELEASE_ID>`.
2. `proteome_denominator.json` records each proteome's resolved absolute path
   (`build_proteome_universe` stores `str(path)` post-`resolve()`). The fixture
   root is replaced with `<FIXTURE_ROOT>`.
3. `CHECKSUMS.json` digests both of the above, so it inherits both. It is
   **recomputed** over the canonical bytes rather than blanked, so the digests
   stay load-bearing: a content change in any release file shows up twice, once
   in the file and once in the checksum.

`work/` is excluded — pipeline scratch, not part of the release contract, and
its `rbh_source_<hash>.dmnd` filenames hash the absolute source path.

Everything else is byte-identical run to run, which is worth knowing on its
own: it means the release format has no embedded timestamps, no iteration
order dependence, and no nondeterministic ordering in the TSVs.

## Files

| File | What it is |
|---|---|
| `run_fixture.sh` | driver: regenerate inputs, run both scenarios, diff against golden |
| `build_inputs.py` | generates `inputs/`; `--check` proves the generation is reproducible |
| `canonicalise.py` | rewrites a release into machine-independent form for diffing |
| `test_fixture.py` | the same checks as pytest assertions, plus stub-level tests |
| `stub_bin/foldseek` | stand-in structural search; emits canned rows from `hits.json` |
| `stub_bin/diamond` | stand-in sequence search; handles `makedb`, forward and reverse `blastp` |
| `stub_bin/hits.json` | every canned number, with per-hit rationale |
| `inputs/` | generated (queries, proteomes, annotations, four manifests) |
| `golden/main/`, `golden/trap/` | canonical expected releases |

`inputs/` is generated rather than committed because the query manifest carries
SHA-256 checksums of both the FASTA sequence and the PDB file, and the pipeline
refuses to run on a mismatch. Hand-maintained checksums are a guaranteed source
of "the fixture is broken and nobody knows why". `build_inputs.py --check`
rebuilds into a scratch directory and compares per-file digests, so a
non-reproducible generator is caught as itself rather than as a pipeline
failure.

## What this fixture does not do

It does not test foldseek or diamond, or any structural or sequence alignment.
It does not validate that the thresholds are scientifically well-chosen — only
that the code applies the ones it is given. The synthetic sequences are
generated from an index formula and the synthetic structures are ideal
CA-only helices; they are not protein-like and are not intended to be. Nothing
in `golden/` is a biological result.
