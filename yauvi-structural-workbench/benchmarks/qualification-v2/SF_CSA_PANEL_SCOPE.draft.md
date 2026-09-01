# sf-csa panel: what it measures — DRAFT, not adopted

Written before any record is curated or any code changed, per the adoption
protocol's closing requirement: state plainly what each gate measures and confirm
it is what the scope claims. The membrane panel satisfied every procedural rule
for its entire life while measuring rotational self-consistency and appearing to
measure orientation accuracy. This document exists to make that failure visible
here before it can happen.

Scope name: `sf_csa:curated_structure_sequence_comparison`.

---

## 1. What the panel measures

**Whether the module's closed-vocabulary classification correctly separates four
curated relationship strata, and — above all — whether it ever promotes a weak
relationship to a functional claim.**

The four declared gates split into two kinds:

| gate | kind | what it asserts |
|---|---|---|
| `exact_controls_recovered = all` | recall | genuinely equivalent pairs are recognised |
| `structural_homolog_controls_recovered = all` | recall | superfamily homologs are recognised |
| `analogy_or_unrelated_promoted_to_function_max = 0` | **false-positive** | analogy and unrelated pairs are *never* promoted to a functional claim |
| `structural_and_sequence_outputs = separate` | contract | the two evidence legs are not merged |

The third is the scientifically load-bearing one. It is an absolute bound, not a
tolerance: a single promotion of a fold analogy to a functional claim fails the
panel. That matches the module's own declared limitation that structural and
sequence similarity "must not be merged into a single similarity claim."

The fourth is a contract check, not a scientific one, and should be labelled as
such so it is not read as evidence of accuracy.

## 2. What the panel does NOT measure

Stated explicitly, because each of these could be mistaken for something the
panel establishes:

- **Not the accuracy of Foldseek or DIAMOND.** Those are pinned runtimes. The
  panel measures how their output is *interpreted*, not whether their alignments
  are correct.
- **Not function.** The module's own contract says substrate and native activity
  "still require direct validation." A passing panel establishes classification
  behaviour, never biochemical function.
- **Not cross-manifest comparability.** Structural-category thresholds live in
  the database manifest, and the module states two releases built against
  different manifests "are not comparable on category alone." The panel is valid
  only against its own frozen manifest.
- **Not structure quality.** Query structures are predicted monomers, not
  experimental assemblies or active poses.

## 3. Freeze items (protocol rule 1)

Fixed before the first run, and checksummed:

- The **database manifest**, because it carries the structural-category
  thresholds. Without pinning it, category results are not comparable between
  runs and the panel measures nothing stable.
- The **mechanism-family, contested-group and divergence tables**
  (**see `SF_CSA_PREADOPTION_FINDINGS.md`, Finding 2 — overriding these is a
  precondition, not a precaution; without it two of four gates are meaningless
  in opposite directions**). Their defaults
  are periodontal-pathogen biology. A four-family fold-relationship benchmark is
  almost certainly not periodontal, so **the panel must override them or it is
  silently using the wrong table.** The module records which table a release
  used, so this is checkable rather than a matter of trust.
- **Runtime versions**: foldseek `10.941cd33`, diamond `2.1.11` — verified
  present and matching the v1 pins.
- The **stratum → label mapping** (see §4), which is a scientific decision and
  not an implementation detail.

## 4. The unreachable label — a decision, not a detail

> **SUPERSEDED IN PART — see `SF_CSA_PREADOPTION_FINDINGS.md`, Finding 1.**
> The claim below that the label "cannot be emitted" is withdrawn: it is
> reachable through the queries manifest, and Option B is unsafe as written.
> The text is left intact as the superseded record.

The module declares a closed six-label vocabulary. **One label,
`probable_same_function`, cannot be emitted.** `classify_hit` reaches it only via

```python
if target_meta and target_meta.get("rbh") and category == "whole_architecture_match":
```

and `target_meta[...]["rbh"]` is never assigned anywhere in the package. The
reciprocal-best-hit computation runs *after* classification and writes its result
to `h["orthology_status"]` on sequence-hit rows — a different structure. Two
independent causes, either of which alone would suffice: wrong order, wrong
destination. The repository records this as an ordering defect, which understates
it; fixing the order alone would not make the label reachable.

This forces a decision before curation, because it determines what the
`homologous_superfamily` stratum expects:

**Option A — exclude it.** Map `homologous_superfamily` to `same_mechanism_class`,
record `probable_same_function` as excluded and why. The panel passes on a
module with a dead code path.

**Option B — gate it.** Include a control case that *should* produce
`probable_same_function`. It will fail, and the panel will correctly report a
known-broken path. This matches the panel's existing fail-closed control pattern
and means the benchmark detects the defect rather than routing around it.

**Recommendation: B.** A qualification panel that avoids a known-broken path is
measuring the path it chose, not the module. Option A is defensible only if the
defect is fixed first, in which case the question disappears.

Either way the choice is recorded before records are curated, not discovered
afterwards.

## 5. Blocking work, in order

1. Decide §4.
2. Register `sf_csa` in the runner's `ENGINES`. It is absent, along with
   `conformational_state`. Invocation differs from the four existing workflows:
   `--queries` and `--databases` manifests, not `--structure`. Evidence document
   is `SF_CSA_RELEASE_MANIFEST.json`.
3. Write `gate_semantics` for the panel. It is currently empty and the runner
   refuses to run without it — correctly, since it will not assume gate meanings.
4. Curate 16 records: four families × four strata, per the coverage rule.
5. Cross-machine reproduction before adoption (protocol rule 7).

## 6. Open questions for Yuvraj

- §4, Option A or B.
- Which four protein families. The coverage rule requires exactly four unrelated
  families with one judgment per stratum each; the choice determines whether the
  periodontal default tables are appropriate or must be overridden.
- Whether `structural_and_sequence_outputs = separate` should be labelled a
  contract check rather than a scientific gate, so a passing panel is not read as
  four scientific results when it is three plus one.
