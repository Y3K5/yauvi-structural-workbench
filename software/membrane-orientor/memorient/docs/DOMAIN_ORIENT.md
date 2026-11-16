# memorient.domain_orient — domain-aware orientation from C-alpha coordinates

A small, offline, numpy-only routine that places a multi-domain membrane protein in a flat bilayer when
only C-alpha coordinates are trustworthy (e.g. predicted models). It targets the cases where a
whole-protein hydrophobic fit fails: outer-membrane β-barrels with large soluble domains
(POTRA domains, plugs, lipoprotein anchors) that drag a global fit sideways.

It is PPM-like in spirit, but much simpler. The result is an **estimate**, not a measured pose.

## Method

1. **Secondary structure** — P-SEA-style assignment from C-alpha distances, plus C-alpha ladder strands
   (pairs at < 5.5 Å whose i±2 neighbours are also paired). The ladders are robust on predicted models.
2. **β-barrel core** — strands from both sources are paired by C-alpha contacts. A closed barrel is the
   2-core of the strand graph (≥ 8 strands). If there is none, an open ring is accepted and flagged.
   The axis is fitted to the strand ring.
3. **Lipid-face check** — the mean Kyte-Doolittle value of outward-facing core residues in the slab must
   be > 0. Soluble barrels (GFP, streptavidin, lipocalins) fail this check and are rejected.
4. **Pose scan** — tilt within a 25° cone and depth within ±14 Å. Score = outward hydropathy in the slab
   (default half-thickness 12 Å) + aromatic-girdle reward − penalty for non-core residues in the slab.
5. **Side (sign)** — votes, each recorded in `sign_evidence`:
   - N-terminal soluble domains (POTRA) far from the slab face the periplasm. Plugs within half+8 Å are ignored.
   - Longer loops face outside.
   - The first strand starts on the periplasmic side.
   - A lipid anchor sits in the periplasmic leaflet.
6. **Helical core** — apolar helices (mean KD ≥ 1.0, ≥ 14 residues, ≥ 18 Å) after `skip_n` residues, so a
   signal-peptide helix is not counted. The side comes from the positive-inside rule. A single helix
   is reported as low confidence.
7. **No core** — `method: none`. Nothing is placed.

Output: `normal` (+normal = extracellular / outside), `origin`, `half`, `tilt`, `method`, `n_strands`,
`confidence` (`high` / `medium` / `low`), `sign_evidence`, and `rejected_barrel` when applicable.

```bash
python -m memorient.domain_orient model.pdb --chain A --skip-n 25 --json pose.json
```

```python
from memorient.domain_orient import read_pdb_ca, orient
ca, seq, num = read_pdb_ca("4c00.pdb", "A")
pose = orient(ca, seq, num)
```

## Command and Workbench

`memorient domain-orient` runs the method on a local PDB file and writes four files into a new folder:

```bash
memorient domain-orient model.pdb --chain A --skip-n 25 --out pose/
```

- `DOMAIN_ORIENTATION.json`: the result above, the input file name and SHA-256, the settings, the core
  residues, zone counts, and the claim limit.
- `ORIENTED_STRUCTURE.pdb`: the first model's atoms moved into the estimated frame (slab centre at the
  origin, membrane normal along +Z, +Z = estimated outside).
- `MEMBRANE_LAYER.json`: slab bounds and residue zone colours for a 3D viewer.
- `RESIDUE_ZONES.tsv`: per residue, height Z (Å) and zone (`membrane_core` within ±half, `interface` up to
  5 Å beyond, then `positive_side` / `negative_side`). Zones use the rounded normal and origin, so the
  core count can differ by one or two residues from `in_slab`.

Exit codes: `0` core found and pose written; `1` no membrane core (only the JSON is written); `2` input
error (nothing written). The same output is byte-identical when rerun on the same file.

In the Structural Workbench this is the analysis type `domain_membrane_orientation`
(*Membrane orientation · domain-aware estimate*): the run writes these files under
`outputs/memorient/`, the report lists the residue-zone table and the recorded gaps, and the viewer
offers the oriented model with the membrane slab. Tests: `tests/test_domain_report.py` (command) and
`software/structural-workbench/tests/test_domain_orientation_workflow.py` (Workbench). Both use invented
coordinates and are software checks only.

## Validation (public structures)

The structures are superposed onto the OPM-oriented file, and the predicted normal and slab centre are
compared with OPM's z axis. OPM places the periplasmic side of outer-membrane proteins at −z.

| PDB | protein | angle to OPM normal | depth offset | side | confidence |
|---|---|---|---|---|---|
| 4C00 | TamA | 9.8° | −1.0 Å | agrees | high |
| 4K3B | BamA | 2.8° | −0.1 Å | agrees | high |
| 5D0O | BamA (open) | 12.4° | 0.8 Å | agrees | low (open ring) |
| 1BXW | OmpA | 7.4° | 1.4 Å | agrees | high |
| 2OMF | OmpF | 8.5° | −1.7 Å | agrees | high |
| 1QD6 | OMPLA | 10.8° | 0.2 Å | agrees | high |
| 1C3W | bacteriorhodopsin | 6.9° | −4.5 Å | agrees | medium |

The soluble barrels 1EMA (GFP), 1STP (streptavidin) and 1BBP (lipocalin) are rejected: their lipid-face
KD is ≤ −0.26, versus ≥ 1.27 for the outer-membrane proteins above.

Tests: `tests/test_domain_orient.py`. The offline tests use invented coordinates only. To run the public
OPM comparison, set `MEMORIENT_NETWORK_TESTS=1` (downloads from RCSB/OPM) or set `MEMORIENT_PDB_DIR`
to a directory of `rcsb_<id>.pdb` and `opm_<id>.pdb` files.

## Limits

- C-alpha only. There are no side chains, no electrostatics and no implicit-solvent energy, so this is
  not PPM.
- A flat slab with a fixed thickness.
- The sides of single-helix proteins and lipoproteins are low confidence.
- Open or partial barrels are flagged as low confidence.
- The side vote can disagree with a full-atom orientation for small barrels without soluble domains. Treat
  `sign_evidence` as the explanation, and treat disagreement as a reason to check.
