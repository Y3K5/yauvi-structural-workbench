# actstate

Is this protein in a functionally competent state?

Five signals are computed and **reported separately** — active-site
completeness, geometry, cofactor occupancy, conformation against references of
known state, and assembly dependence. They are never summed, averaged, or
collapsed into one score, for the same reason `sf-csa` keeps structural and
sequence similarity apart: they answer different questions, fail in different
ways, and a reader told only the total cannot tell which one carried it.

```bash
actstate run --in examples --out results/
actstate describe                 # machine-readable IO contract
actstate validate --in examples   # check inputs without running
actstate fetch --plan             # what raw files are needed, and where from
```

No workspace, no project, no campaign. It runs on files you already have.

## The labels

A closed six-value vocabulary. Nothing else is ever emitted.

| label | meaning |
|---|---|
| `active_state_supported` | every signal evaluated and consistent, on experimental coordinates |
| `probable_active` | residues intact, but a needed signal is unavailable or the structure is predicted |
| `apo_but_competent` | site intact and unoccupied; a declared cofactor is absent from the coordinates |
| `inactive_conformation` | residues present but not mutually positioned as a site |
| `active_site_disrupted` | an annotated position does not hold the residue an experimentally validated reference carries there. Requires `--expected-residues`; see rule 4 |
| `indeterminate` | not enough is annotated to make any claim |

## Four rules

1. **No annotated active site means `indeterminate`, never `inactive`.** Absence
   of annotation is not evidence of absence of function. Most proteins in a
   proteome carry no `ACT_SITE` line and are not thereby pseudoenzymes.
2. **A predicted model alone can never yield `active_state_supported`.** A
   predictor reproduces the fold it was trained to reproduce; that is not an
   observation of a functional state. The reader detects AlphaFold and similar
   headers and the label is capped accordingly.
3. **An unavailable signal is recorded as unavailable.** It never becomes a
   neutral value, and it never quietly drops out of the summary. `probable_active`
   most often means *a signal could not be evaluated*, not *the evidence was weak* —
   the `rationale` field says which.
4. **`active_site_disrupted` requires a position-specific expected residue.**
   Supply it with `--expected-residues`, a JSON map keyed by accession then by
   annotated position, holding the residue an experimentally validated reference
   carries there. Without one, a residue outside the competence set is reported
   as contradicting evidence and the label caps at `indeterminate` — the
   observation is kept and the `rationale` names the positions that would need
   an expectation. Membership in a broad residue set is not a chemistry test for
   a position.

## The signals

| signal | needs | says |
|---|---|---|
| `completeness` | annotation + sequence | do the annotated catalytic positions hold residues that can do chemistry? |
| `geometry` | annotation + structure | are those residues clustered as one site? |
| `occupancy` | annotation + structure | apo or holo — is a declared cofactor actually present? |
| `conformation` | structural aligner + curated references | does this resemble a known active or inactive state? |
| `assembly` | `fold_state` output | is this the isolated fold or the working assembly? |

`completeness` compares against an expected residue where one is supplied, and
otherwise against a generic 13-letter competence set. The generic screen is not
role-aware, and it is wrong in both directions on its own: it cannot separate a
real substitution from a sequence numbered against a different entry, and it is
blind to the commoner degradations that stay inside the set — catalytic Cys to
Ser, His to Asn. So it may cap a claim and may not make one. Only the
position-specific comparison reaches `active_site_disrupted`, and only it
detects a within-set substitution. An expectation the pipeline cannot check —
a position carrying no `ACT_SITE` annotation, or anything that is not a standard
one-letter code — is rejected, named in `rejected_expectations`, and has no
effect on the label.

`occupancy` ignores a maintained set of common solvent, buffer, and
cryoprotectant components. It presently detects candidate non-solvent
heteroatoms but does not prove that the observed component identity matches the
declared cofactor; mismatches require expert review.

`conformation` reports `unavailable` out of the box. It needs a curated set of
reference structures of known state, and **no such set ships with this module** —
so the honest thing is to say so rather than to quietly stop asking the question.
Supply one with `--reference-comparison`.

## Inputs

Only the annotation table is required.

```
examples/
  annotations.tsv              UniProt export (Entry, ft_act_site, ft_binding, cc_cofactor, ec)
  sequences.fasta              optional; accepts sp|ACC|NAME or bare accession headers
  structures/                  optional; PDB or mmCIF, named AF-<acc>-F1.pdb or <acc>.pdb
  fold_state.json              optional; keyed by accession
  reference_comparison.json    optional; keyed by accession
```

Column names are matched against both UniProt's field names (`ft_act_site`) and
its display headers (`Active site`), so an export from either the API or the web
UI works unchanged.

`yauvi-fetch get --for actstate` retrieves the annotation export with the right
columns. The four features this module needs — `ft_act_site`, `ft_binding`,
`ft_site`, `cc_cofactor` — are absent from most existing exports, so an older
annotation TSV will produce `indeterminate` for everything. `actstate validate`
says so explicitly rather than letting you discover it from the results.

## Outputs

`ACTIVITY_STATE.json` — the label, the rationale, and all five signals with their
own reasons, per protein.
`ACTIVITY_STATE.tsv` — one row per protein, for reading alongside other channels.

Both are **byte-deterministic**: sorted keys, no timestamps, no absolute paths.
Two runs over the same inputs produce identical bytes, which is what makes
recording a digest meaningful.

## Dependencies

None. Standard library only, following the house rule for these packages: heavy
tools are consumed through their outputs rather than imported. The structure
reader and the geometry test are written against the handful of fields they
need, so neither numpy nor biopython is required.

    pip install actstate                 # the module
    pip install 'actstate[sources]'      # adds yauvi-fetch integration

## Tests

    pytest                                  # offline, no external runtime
    pytest -m "adapter or network"          # adds the external-tool checks

The bundled fixtures are synthetic and exercise every label in the vocabulary
exactly once, so the golden tests need no download and no licence.

## Inside the platform

Registered as `catalogs/modules/activity_state.yaml`, adapter at
`platform/src/yauvi_platform/modules/native/activity_state.py`. Only
`active_state_supported` maps to a passing `EvidenceState`; `probable_active`
maps to `indeterminate`, because a signal that could not be evaluated is not a
favourable one.

It composes with `fold_state` rather than duplicating it: that module says
*which* fold the evidence describes, this one asks whether that fold is competent.
