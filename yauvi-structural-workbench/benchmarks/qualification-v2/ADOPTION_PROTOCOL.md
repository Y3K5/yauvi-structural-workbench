# Adoption protocol

Binding for the two remaining release-blocking panels — `conformational_state`
(ABL StateAtlas, 18 records) and `sf_csa` (16 records).

Every rule below exists because its absence caused a specific failure during the
membrane investigation of 2026-08-30/31. The incidents are named so the rule is
not mistaken for ceremony. Full provenance is in `PANEL_MANIFEST.json`
(`threshold_revisions`, `record_corrections`) and
`MEMBRANE_OBJECTIVE_FINDINGS.md`.

---

## 1. Freeze the inputs before running anything

Fix the exact records, the fields each gate reads, the identity of the code under
test, the mapping rules, and the expected record counts.

> **Caught:** the panel spent an entire run executing unmodified `memorient` —
> a stale copy in `site-packages` shadowed the editable install, so a rotation
> change had no effect and the result looked exactly as predicted. Separately, a
> gate was named `normal_drift_across_20_rotations_deg_max` while the validation
> performed five.

## 2. Prove the gate can fail

Corrupt a case deliberately and confirm the gate rejects it. Do this on a copy,
never on the adopted evidence.

> **Caught:** the digest chain in `RELEASE_STATUS.json` recorded a sha256 beside
> every evidence document and nothing compared them; a stale membrane digest sat
> in the public repository unnoticed. Separately, the membrane accuracy bound of
> 15° sat inside the failure mode and could not separate a converged fit from a
> failed one — it passed an 8.39° error.

## 3. Confirm every record was consumed, before reading any verdict

Check the count, and check that each gated field is present in every record.
Then read the result.

> **Caught:** a broken virtualenv made all 16 cases fail at the CLI, and because
> the runner clears its output directory first, the committed evidence was
> replaced by nothing. The summary line still parsed. On a later run only 1 of 16
> evidence files carried the field a new gate depended on.

## 4. Cross-check against an independently computed value

Compute at least one gated quantity a second way. Disagreement is the signal;
agreement proves nothing on its own.

> **Caught:** an accuracy gate read the fitted normal from the *canonicalised*
> frame rather than the input frame — 29.54° against a true 16.30° on 1T16. Both
> the gate and the independent measurement were wrong, and their disagreement is
> what forced a read of the source.

## 5. Fix thresholds and withdrawal criteria before results exist

Write the pass threshold, the verdict bands and the kill condition into the
script or the manifest before the first value is computed.

> **Caught:** the circumferential-coherence test returned 7 of 16 against a
> preregistered withdrawal threshold of 10 and was withdrawn intact. Adjusting
> the percentile, flank width or sector minimum until it passed would have been
> fitting, and would have looked like a discovery.

## 6. Preserve exclusions, failures and withdrawn claims beside the survivors

Record what was excluded and why, what failed, and every claim that was made and
later withdrawn.

> **Caught:** eight claims were withdrawn across the membrane investigation. Each
> was reported after a measurement that agreed with the explanation then current,
> and corrected only when a second measurement disagreed. The surviving results
> are trustworthy in proportion to that trail being visible.

## 7. Reproduce cross-machine before calling a scope adopted

A stratum that passes only where its expectations were recorded is not adopted.

> **Caught:** the membrane stratum passed 16/16 on the recording machine and
> 9–10/16 on CI. The recording machine was an x86_64 build under Rosetta on arm64
> hardware, matching neither CI platform.

---

## Panel-specific notes

**`sf_csa` (16).** Needs Foldseek and DIAMOND, version-pinned in v1 as
`foldseek 10.941cd33` and `diamond 2.1.11`. Read how `sf_csa` invokes them before
writing code — the module is at `sf-csa/src/sf_csa/`, at the repository root, not
under `yauvi-structural-workbench/`. It carries a known defect: reciprocal-best-hit
is computed *after* structural classification, so `probable_same_function` is
unreachable end to end.

**`conformational_state` (18).** `execution_policy` sets
`predicted_structures_allowed_for_state_atlas: false`, so every ABL case must be
experimental. Tighter curation than any panel adopted so far.

## What this protocol does not cover

It governs adoption, not scientific validity. A panel can satisfy all seven rules
and still gate the wrong quantity — the membrane panel did exactly that for its
entire life, measuring rotational self-consistency while appearing to measure
orientation accuracy. Before adopting, state plainly what each gate measures and
confirm it is what the scope claims.
