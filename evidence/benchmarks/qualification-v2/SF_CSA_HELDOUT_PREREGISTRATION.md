# Held-out test of the sf-csa panel — preregistered 2026-09-02

Written **before the data was acquired and before anything was run.** Nothing
below is adjusted afterwards; the comparison is made once.

## Why

Every record in `ADOPTION_DRAFT_SFCSA.json` is `curator_frozen`. There is no
held-out split, so 16/16 is consistent with two different worlds: the module
classifies fold relationships correctly, or the selection was shaped until it
appeared to. Those are distinguishable only on proteins the module has not seen.

## The held-out set

SCOP fold **d.15, beta-Grasp (ubiquitin-like)** — not used by any of the four
panel families, and containing three superfamilies whose entries pass all five
selection gates including proteome membership.

| stratum | query | target | basis |
|---|---|---|---|
| `exact` | 1DOI | 1DOI | accession identity (P00217) |
| `homologous_superfamily` | 1DOI | 1DOX | both `d.15.4.1`, 2Fe-2S ferredoxin-like; P00217 vs P27320 |
| `fold_analogy` | 1DOI | 1EO6 | fold `d.15`; `d.15.4` vs `d.15.1` ubiquitin-like |
| `unrelated` | 1DOI | 1CSP | folds `d.15` vs `b.40` |

Mechanism groups are the SCOP superfamily, as in the panel:
`scop_d15_4_2fe2s_ferredoxin_like`, `scop_d15_1_ubiquitin_like`,
`scop_b40_4_nucleic_acid_binding`.

## Predictions

1. `1DOI -> 1DOI` is **`exact_function_supported`**.
2. `1DOI -> 1DOX` is **`same_mechanism_class`**.
3. `1DOI -> 1EO6` is **not** in `{exact_function_supported, probable_same_function}`.
4. `1DOI -> 1CSP` is **not** in `{exact_function_supported, probable_same_function}`.
5. Promotions of an analogy or unrelated pair to a functional claim: **zero**.

Predictions 3 and 4 are deliberately bounds rather than exact labels, matching the
panel's own `stratum_expected_label`, which asks only that these strata stay
outside the function-claim labels. Committing to `structural_analogy_only`
specifically would be predicting a threshold outcome, not a classification — the
mistake made on 2026-09-01 and withdrawn.

## Kill criterion

**Any promotion of the analogy or the unrelated pair to a function claim falsifies
the panel's central bound on unseen data**, and the 16/16 result should not be
trusted until it is explained. A wrong label that is not a promotion — an exact
match not recovered, say — is a weaker failure but still a failure of prediction
1 or 2, and is reported as one.

A pair producing no structural row at all is **not** counted as a pass. Under
collection 2.6 the search is exhaustive, so every declared pair must yield a row
stating where it fell.

## What this cannot show

Four judgments on one fold. A pass is evidence the machinery generalises beyond
the frozen twelve; it is not a second qualification panel, and it is not
independent — it runs on the same machine, by the same author, against the same
tool versions.
