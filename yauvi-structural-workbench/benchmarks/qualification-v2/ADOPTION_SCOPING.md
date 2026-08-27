# Qualification v2 adoption scoping — StructQC panel

Scoping date: 2026-08-26
Panel: `qualification-v2-structqc` (16 of the 114 required cases)
Current state: `blocked_panel_incomplete`, 0 records adopted, 0 sources locked

This document turns "114 cases" into a concrete work plan, using StructQC as the
worked example. It adopts nothing and changes no threshold. Every claim about
existing artifacts below was verified by reading the files.

---

## 1. What adoption actually requires

The runner (`run_qualification.py`) enforces a strict contract. Each record needs
**13 fields**:

`record_id`, `record_kind`, `stratum`, `split`, `expected_result`,
`source_release`, `checksum`, `license`, `citation`, `artifact`,
`exclusion_rationale`, `pdb_entry_id`, `homolog_group`

And four structural rules:

1. **The artifact must exist on disk and match its recorded SHA-256.** A missing
   or drifted file is an error, not a skip.
2. **The artifact path must resolve inside `qualification-v2/`.** Absolute paths
   and `..` are rejected.
3. **No PDB entry may appear in both `development` and `held_out`.**
4. **No `homolog_group` may straddle the two splits.** This is the stricter rule
   and the one most likely to be tripped by accident — see §3.

`SOURCE_LOCK.json` separately records `artifact`, `provider`, `source_id`, `url`,
`sha256` per source, following the v1 shape.

The StructQC panel requires **4 strata × 2 splits × 2 cases = 16**:

| Stratum | development | held_out |
|---|---:|---:|
| `x_ray` | 2 | 2 |
| `cryo_em` | 2 | 2 |
| `nmr` | 2 | 2 |
| `alphafold` | 2 | 2 |

Gates: `residue_identity: exact`, `official_metric_import: exact`,
`missing_evidence_behavior: fail_closed`.

Coverage the panel must exercise: insertion codes, alternate locations, missing
residues, nonstandard residues, multichain coordinates, pLDDT, PAE, and
experiment-specific wwPDB validation.

---

## 2. What Qualification v1 can contribute — verified, not assumed

v1 locked 27 artifacts. I read `_exptl.method` from all 11 wwPDB mmCIF files:

| Entry | `_exptl.method` |
|---|---|
| 1B73, 1CRN, 1OAI, 2G1T, 2HZ4, 2JFN, 2V7A, 4EA9, 4HHB, 4HHB-assembly1, 8SSN | **X-RAY DIFFRACTION** (all 11) |

**Every v1 experimental entry is X-ray. There is no cryo-EM and no NMR anywhere
in the v1 collection.** Plus one AlphaFold model set (AF-P69905-F1 v6: model,
PAE, API metadata).

> **Curation note.** Determine the stratum from the `_exptl.method` key only.
> Grepping the whole mmCIF for `SOLUTION NMR` or `ELECTRON MICROSCOPY` produces a
> false positive on 2V7A, whose citation block references the NMR entry 1AWO.
> Likewise, v1 appears to hold two validation files but they are `1crn_validation.xml`
> and its `.gz` — one report, one entry.

Resulting coverage against the StructQC panel:

| Stratum | Needed | v1 candidates | Gap |
|---|---:|---:|---:|
| `x_ray` | 4 | 11 | **0** — ample choice |
| `alphafold` | 4 | 1 | **3** |
| `cryo_em` | 4 | 0 | **4** |
| `nmr` | 4 | 0 | **4** |
| **Total** | **16** | **5** | **11** |

`candidate_migrations` already records v1 as `not_adopted` because its artifacts
"do not satisfy v2 split, mapping, stratum, or expected-result metadata by
themselves." That is accurate: re-adopting an X-ray entry still requires adding
`stratum`, `split`, `expected_result`, `license`, `citation`, `homolog_group`,
and `exclusion_rationale`. The coordinates carry over; the curation does not.

### The validation-report gap is larger than the coordinate gap

The `official_metric_import: exact` gate needs an official validation report per
experimental case. v1 locked **exactly one**: `PDB:1CRN:validation_xml`.

Twelve experimental cases (4 x-ray + 4 cryo-EM + 4 NMR) each need one, so
**11 additional wwPDB validation XMLs** must be acquired even for the four x-ray
cases drawn from entries already on hand.

---

## 3. A leakage trap already present in the candidate set

4HHB is human deoxyhaemoglobin. AF-P69905-F1 is the AlphaFold model of UniProt
P69905 — haemoglobin subunit alpha. **They are the same protein.** If 4HHB lands
in `development` and the AlphaFold case in `held_out` (or vice versa), the runner
fails the panel on homolog-group leakage — correctly.

This has to be settled by assigning `homolog_group` deliberately across all four
strata before any record is written, not discovered case by case. The same
question applies to 2V7A and 8SSN, which are both ABL-family entries used by the
StateAtlas panel.

---

## 4. Two design questions the frozen manifest does not answer

Both need a scientific decision before records can be written. Neither is a
threshold change, so neither breaks the immutability policy.

**(a) `official_metric_import: exact` means something different per stratum.**
X-ray validation reports carry clashscore, Ramachandran and rotamer outliers,
plus resolution and R-factors. NMR reports have no resolution or R-factors and
are multi-model. Cryo-EM reports substitute map-model metrics. AlphaFold has no
wwPDB validation report at all — its equivalents are pLDDT (in the B-factor
column) and the PAE JSON. The gate needs a per-stratum definition of which
metric set must import exactly, or it is untestable for three of four strata.

**(b) NMR ensembles force a model-selection decision.** NMR entries deposit many
models. StructQC has `--model` (default 0). The panel must state whether the
expected result is defined against model 1 only, or whether ensemble handling is
itself part of what the stratum tests. This interacts with the "multichain
coordinates" and "alternate locations" coverage rules.

---

## 4b. What drafting the x-ray four actually uncovered

Two findings emerged from reading the candidate coordinates, and both change the
work rather than merely describing it.

### No candidate on hand has an insertion code

The coverage rules require the panel to exercise insertion codes. Across all ten
distinct v1 entries, the insertion-code count is **zero**. So even the "free"
x-ray stratum cannot fully satisfy its coverage rule from material already
acquired; at least one entry with author-assigned insertion codes (antibody
structures under Kabat/Chothia numbering are the usual source) must be curated
in. This is recorded honestly in the draft's `coverage_achieved` as
`insertion_codes: NOT COVERED`.

### `residue_identity: exact` collides with real deposited entries

Reading `_struct_ref_seq_dif` across the candidates:

| Entry | Declared differences vs reference | Kind |
|---|---:|---|
| 1OAI | 0 | none |
| 8SSN | 2 | expression tag |
| 4EA9 | 6 | expression tag |
| 2V7A | 2 | **engineered mutation T315I** |
| 1CRN (the v1 passing case) | **0** | none |

The single public case that has ever passed this gate is the one entry with no
sequence differences at all. That is worth stating plainly: the `exact` identity
gate has only ever been exercised on the easiest possible input.

How much this matters depends on how StructQC computes identity, so I checked.
`structqc/src/structqc/core.py` runs a **global Biopython pairwise alignment** and
reports `identity_fraction = identity / len(mapping)` — identity over aligned
positions, not over reference length. The consequence:

- **Expression tags are harmless.** Extra residues absent from the reference are
  gapped by the global alignment and do not depress identity. 8SSN and 4EA9 can
  still reach 1.0.
- **Engineered point mutations are not.** T315I aligns and mismatches, so 2V7A
  lands at roughly `1 - 1/mapped` ≈ 0.998, never 1.0.

**So the gate needs one decision, and it is not a threshold change:** does
`exact` mean *identity must equal 1.0*, or *the computed identity must exactly
match the case's declared expectation*? The first reading disqualifies every
mutant structure in the PDB and makes the panel far harder to populate. The
second is stronger science — it turns 2V7A into a positive test that a declared
engineered mutation is **detected rather than silently absorbed**, which is
precisely what a coordinate-trust tool should do.

The draft is written for the second reading. If the first is chosen instead,
2V7A must be replaced and the nonstandard-residue coverage re-sourced.

## 4c. Executing the x-ray four — what running it actually showed

The four x-ray artifacts were acquired and StructQC was run on each. All 11
sources are locked and checksum-verified, and the x-ray half of the panel now
composes (2/2 development, 2/2 held_out, no leakage, no missing fields).

| Entry | chain | exit | identity | coverage | mapped/ref | breaks |
|---|---|---:|---:|---:|---|---:|
| 4EA9 | — | 1 | 1.000000 | 0.962791 | 207/215 | 1 |
| 1OAI | — | 0 | 1.000000 | 0.095315 | 59/619 | 0 |
| 2V7A | A | 1 | 0.996296 | 0.238938 | 270/1130 | 1 |
| 8SSN | A | 1 | 0.997619 | 0.371681 | 420/1130 | 5 |

**Source stability, verified.** All four freshly downloaded coordinate files were
byte-identical to the SHA-256 values recorded in the Qualification v1 source
lock. The upstream artifacts have not drifted, which is what makes a checksum
lock worth keeping.

### Three behaviours the panel must encode

**(a) Chain selection is mandatory, and getting it wrong is silent.** 2V7A and
8SSN are homodimers — two copies of the same entity. Run without `--chain`,
StructQC concatenates both copies and aligns them against a single reference,
giving identity **0.649** for 2V7A and **0.656** for 8SSN. With `--chain A` the
same entries give **0.996296** and **0.997619**. Nothing warns you; the low
number just looks like a bad structure. Every multi-copy case in all six panels
must record its selected chain as part of the expected result.

**(b) `2V7A` at 0.996296 is exactly `1 - 1/270` — the single declared T315I
mismatch, and nothing else.** The engineered mutation is detected rather than
absorbed. This settles the §4b question in favour of the second reading: `exact`
must mean *the computed identity matches the case's declared expectation*, not
*identity equals 1.0*. Under the first reading this correct result would be
recorded as a failure. Expression tags behave as predicted and cost nothing —
4EA9 carries six and still reaches 1.0.

**(c) Exit code 1 is the normal result for real crystal structures.** Three of
four cases exit 1, driven by chain breaks in the selected chain — disordered
loops, which nearly every deposited structure has. Only 1OAI, with zero breaks,
exits 0. The panel must treat exit 1 as an expected scientifically-incomplete
result carrying its reason, not as failure. An executor that requires exit 0
would reject almost the entire PDB.

### A coverage rule that cannot currently be met

2V7A was selected to carry the panel's **nonstandard-residue** coverage: gemmi
sees `PTR393` (phosphotyrosine) in both chains. StructQC reports
`nonstandard_residues: 0` and `residues: 270` where gemmi counts 271, and the
string `PTR` appears nowhere in the evidence record, the residue table, or the
run manifest. The counter is `nonstandard += aa == "X"`, so a modified residue
that maps to its unmodified parent — phosphotyrosine to tyrosine — is normalised
away silently.

The direction is conservative: it omits information rather than inventing a
favourable claim. But for a coordinate-trust tool, a post-translational
modification is precisely the thing a user needs surfaced before mapping a
functional residue, and the panel's `nonstandard_residues` coverage rule cannot
be satisfied by any entry whose modification maps to a standard parent.

**Two of the panel's coverage rules are therefore unmet by the current x-ray
four**: insertion codes (no candidate has any) and nonstandard residues
(reported as zero). Both need a decision before the x-ray stratum is adopted
into the frozen manifest.

## 5. Adoption is not execution — the executor does not exist

`run_qualification.py` is a **composition audit only**: it validates that records
are present, checksummed, non-leaking, and completely specified. It sets
`scientific_execution_performed: false` unconditionally and never invokes an
analysis engine. The manifest says so: "a later immutable executor version must
run the canonical engines and evaluate every gate."

So there are two distinct pieces of work, and finishing the first changes the
state from `blocked_panel_incomplete` to composed-but-unexecuted:

1. **Adopt** — curate entries, acquire artifacts, record 13 fields per case, lock
   sources. → panel composes, runner exits 0.
2. **Execute** — write the executor that runs StructQC per case, compares against
   `expected_result`, evaluates the three gates, and records pass/fail. → the
   scope can become qualified.

v1's runner is the working reference for step 2: 800 lines that ran six workflows
and produced the check structure already in use (`experimental_cli_completed`,
`gemmi_coordinate_validation`, `author_label_identity_present`,
`reference_mapping_exact`, `wwpdb_validation_imported`,
`wwpdb_geometry_metrics_imported`). Those checks map cleanly onto the v2 gates.

**`execution_policy.network_access` is `forbidden`**, so acquisition is strictly a
prior step; execution reads only locked local files. `independent_second_machine_required`
is `true`, so a passing local run is not the end of the gate either.

---

## 6. What already works, and does not need building

The acquisition layer is done. `yauvi-fetch plan --for structqc` resolves all
four sources StructQC needs, each with a fetch policy:

| Source | Role |
|---|---|
| `pdb` | experimental coordinates and declared experimental metadata |
| `alphafold_db` | predicted coordinates, PAE, prediction confidence |
| `uniprot_proteomes` | reference sequence for mapping and completeness |
| `wwpdb_validation` | community geometry and experiment-specific validation |

So the gap is **curation and the executor**, not plumbing.

---

## 7. Effort shape

For the StructQC panel:

| Step | Work | Blocked on |
|---|---|---|
| Resolve §4(a) and §4(b) | Scientific decision, per-stratum metric definitions | Judgement |
| Select 16 entries across 4 strata, assign homolog groups and splits | Curation; the x-ray four can come from v1 | §4 decisions |
| Acquire 11 coordinate sets + 11 validation XMLs + reference sequences | `yauvi-fetch get`, then checksum | Network, one session |
| Write 16 records and the source lock | Mechanical once curation is settled | Above |
| Run the composition audit | Already built; exits 0 when complete | Above |
| Build the executor and define expected results per case | New code, modelled on the v1 runner | All of the above |
| Second-machine reproduction | Policy requirement | Executor |

**Extrapolating to all six panels** (114 cases): the other five need 98 more
cases. Membrane is the largest at 32 and its α-helical half is explicitly
non-blocking. StateAtlas at 18 is constrained further —
`predicted_structures_allowed_for_state_atlas` is `false`, so every ABL case must
be experimental. The per-panel executor logic differs, but the adoption contract
and acquisition path are shared, so the second panel should cost materially less
than the first.

---

## 8. Recommended first move

Adopt the **four x-ray cases** end to end: entries already on hand, one stratum,
one metric definition, no new scientific decisions beyond homolog grouping. That
exercises the entire pipeline — curation, acquisition of the missing validation
XMLs, record writing, checksum locking, composition audit — on the cheapest
possible slice, and produces a template the remaining 110 cases follow.

It will also flush out §4(a) concretely, because the x-ray metric set is the one
the existing v1 checks already validate.
