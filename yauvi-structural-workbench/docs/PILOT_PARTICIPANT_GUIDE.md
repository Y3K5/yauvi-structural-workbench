# Pilot participant guide

This is the sheet a pilot participant receives. It implements
[PILOT_PROTOCOL.md](PILOT_PROTOCOL.md), which the facilitator follows. No
participants are enrolled yet.

## Before the session (facilitator)

- Prepare the exact wheel under test and write its SHA-256 into the record
  (`artifact_sha256`).
- Build the showcase: `python3 tools/build_3d_showcase.py --workbench`, run
  from the repository root after `evidence/benchmarks/qualification-v2/acquire_sources.py`.
  It writes a standalone page, `build/showcase-3d/index.html`, and a copy the
  workbench can serve, `build/showcase-3d/workbench/`.
- Assign one workflow and one reviewed public or synthetic input bundle. Do not
  ask for confidential research inputs.
- Copy [PILOT_RECORD_TEMPLATE.json](PILOT_RECORD_TEMPLATE.json) to a private
  location outside this repository. Filled records are never committed here.

## Step 0: orientation (unscored, about 10 minutes)

Open the workbench with the showcase as its landing page:

```bash
YAUVI_SHOWCASE_DIR=build/showcase-3d/workbench yauvi --workspace ./pilot-case workbench open
```

On the welcome screen, choose **See what each tool does**. Each panel shows one
public example from the qualification panel, what the tool reported for it,
what the tool does *not* claim, and the exact commands that produced it. If the
workbench cannot be installed yet, open the standalone `index.html` in a browser
instead; it needs no server or network.

Questions are welcome during this step. Coaching stops when step 1 begins.

### Conformational state: see the structural difference

Open **Conformational state** and choose the activation-loop or DFG close-up.
The aligned endpoints are chain A of 2GQG (active-like, dasatinib-bound) and
2HYY (inactive-like, imatinib-bound). Select a residue to see its Cα displacement,
endpoint side chains, backbone φ/ψ torsions and N–Cα–C bond angle. Asp/Phe
residues also report χ1 when its four atoms are unambiguous.

Use the slider or **Play illustrative motion** to track corresponding Cα
positions. The amber trace is linear interpolation; intermediate shapes are
not experimentally observed, and the connecting lines are not chemical bonds.
The angle table always describes the two actual endpoints. Backbone torsions,
side-chain torsions and bond angles are different measurements.

The comparison fits all 252 shared, same-residue Cα positions in the frozen
P00519 242–495 maps. Positions 275 and 393 are excluded and the trace leaves
gaps. Residue 393 is phosphotyrosine in 2GQG, outside the standard-residue map.
This is a comparison between differently prepared, inhibitor-bound structures,
not a time-resolved transition or a measurement of cellular activity.

**Download geometry and provenance** exports the fit, source hashes, exclusions,
per-residue displacements and endpoint angles. The build also writes
`build/showcase-3d/ABL-endpoint-geometry.json`. These new illustrative measurements
do not replace frozen StateAtlas classifications. Building requires NumPy and
Biopython; the generated page remains self-contained and offline.

### Cross-species comparison: two guided examples

Open **Cross-proteome comparison**. Its four models form two pairs:

| Fold space | Cross-species pair | Research-supported shared role |
|---|---|---|
| Ferredoxin-like (`d.58`) | *G. acidurici* FdxA P00198 / *A. vinosum* ferredoxin P00208 | Two-cluster iron-sulfur electron carriers |
| OB-fold (`b.40`) | *B. subtilis* CspB P32081 / *E. coli* CspA P0A9X9 | Single-stranded nucleic-acid binding in the cold-shock family |

Read the linked primary research beside each pair. The structures displayed are
checksum-verified AlphaFold models; the named PDB entries identify experimental
references. Each model rotates independently; the views are not a superposition.
The colours identify the query and comparison protein, not measured activity.

Follow the evidence in this order:

1. Identify the species, accession, experimental reference and fold family.
2. Read the recorded structural coverage and RMSD. These values come from the
   frozen qualification manifest, not a new analysis run during this session.
3. Read the separate sequence result. The ferredoxin pair has no retained
   sequence-comparison row. The cold-shock pair has a match classified as
   `paralog_or_secondary_candidate`, not a reciprocal best hit.
4. Explain what the cited experiments support independently of the geometry.
   Both pairs carry `same_mechanism_class`, not proof of identical function.
5. Read the same-fold counterexample. Acylphosphatase does not inherit
   electron-carrier function, and pyrophosphatase does not inherit RNA-binding
   function merely from a shared broad fold assignment.

“Similar function” here means the stated shared role within each homologous
family. It does not establish equal substrates, rates, partners or physiological
effects. The separate ferredoxin/flavodoxin literature example explains functional
analogy across different folds; it is not one of the two executed comparisons.

The references and case identities are maintained in
[SF_CSA_SHOWCASE_CASES.json](SF_CSA_SHOWCASE_CASES.json). SF-CSA reproduction
requires the locked source collection plus Foldseek and DIAMOND. The displayed
commands document the qualification setup; use a fresh output directory outside
the source or publication tree for a learner's run and preserve the frozen inputs.

## Scored task (no coaching)

1. Install the wheel you were given.
2. Create a case for your assigned workflow.
3. Check whether it is ready to run, and resolve anything it asks for.
4. Run it.
5. Find one measurement and show which input it came from.
6. Explain any evidence the run reports as missing.
7. Export the case, then reopen the report.
8. In one or two sentences, state what this result does **not** show.

A run can end three ways. The exit code or status tells you which:

| code | meaning |
|---|---|
| 0 | complete |
| 1 | ran, but some evidence is unresolved. This is reported, not hidden |
| 2 | stopped before producing a result |

An exit of 1 is not a failure of the task. Explaining what was unresolved is
part of step 6.

## After the session (facilitator)

Fill in every field of the record. Separate usability failures from missing
scientific prerequisites. Record any reading that turns geometry into
biological proof under `misunderstandings`, because the protocol treats it as a
blocker. Record a quote only when `consent_to_quote` is true.
