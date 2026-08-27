#!/usr/bin/env python
"""Two worked examples end-to-end, on real structures fetched from RCSB.

1. **Beta-barrel OMP** — OmpA transmembrane domain (PDB 1BXW), gram_negative_om context.
   Barrel-normal fit, extracellular-loop labelling, LPS-shielded band, antibody-accessible
   surface set.

2. **Single-pass TM helix** — glycophorin A TM domain (PDB 1AFO), tm_receptor context.
   Hydrophobic-belt orientation with the positive-inside rule breaking the sign: the
   C-terminal Arg/Lys cluster (…RRLIKK) is placed cytoplasmic, the N-terminal ectodomain
   extracellular — glycophorin A's known topology, recovered from geometry + charge alone.

Each example prints a summary + a per-residue label table and writes, into examples/out/:
  <name>_oriented.pdb   — coordinates in the membrane frame (+Z extracellular, core at origin)
  <name>_viz.json       — 3Dmol.js descriptor (per-residue colours + membrane slab)
  <name>.pml            — PyMOL script (colour by accessibility, epitope surface as sticks)
  <name>_labels.tsv     — full per-residue label table

Run:  python examples/worked_examples.py     (needs network + compute extras)
"""

from __future__ import annotations

import json
import os
import urllib.request

from memorient.contexts import get_context
from memorient.geometry import load_structure
from memorient.orientor import orient_structure
from memorient.viz import display_oriented, write_3dmol_html, write_pymol_script

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "out")
RCSB = "https://files.rcsb.org/download/{pdb}.pdb"

EXAMPLES = [
    {"pdb": "1BXW", "chain": "A", "context": "gram_negative_om", "name": "ompa_barrel"},
    {"pdb": "1AFO", "chain": "A", "context": "tm_receptor", "name": "glycophorin_tm_helix"},
]


def _fetch(pdb: str) -> str:
    path = os.path.join("/tmp", f"{pdb}.pdb")
    if not os.path.exists(path):
        urllib.request.urlretrieve(RCSB.format(pdb=pdb), path)
    return path


def run_one(ex: dict) -> dict:
    os.makedirs(OUT, exist_ok=True)
    s = load_structure(_fetch(ex["pdb"]), chain=ex["chain"])
    ctx = get_context(ex["context"])
    result = orient_structure(s, ctx, n_points=240, validate=True)

    base = os.path.join(OUT, ex["name"])
    result.write_pdb(base + "_oriented.pdb")
    with open(base + "_viz.json", "w") as fh:
        json.dump(display_oriented(result), fh, indent=2)
    write_pymol_script(result, base + ".pml")
    write_3dmol_html(result, base + ".html")
    rows = result.residue_table()
    hdr = ["resid", "resname", "chain", "zone", "facing", "accessibility", "extracellular", "rsa"]
    with open(base + "_labels.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in hdr) + "\n")

    sm = result.summary()
    print(f"\n=== {ex['pdb']} ({ex['context']}) -> {result.label} ===")
    for k in ("method", "confidence", "n_residues", "n_extracellular", "n_surface_set",
              "host_antibody_accessible", "half_thickness", "delta_kd", "rotation_invariant",
              "mean_jaccard", "ec_sign_confidence"):
        if k in sm:
            print(f"  {k}: {sm[k]}")
    print("  metrics: " + ", ".join(f"{k[7:]}={v}" for k, v in sm.items() if k.startswith("metric.")))
    return {"pdb": ex["pdb"], "context": ex["context"], "summary": sm,
            "surface_set": sorted(result.labels.surface_set)}


if __name__ == "__main__":
    out = [run_one(ex) for ex in EXAMPLES]
    with open(os.path.join(OUT, "worked_examples_summary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote outputs to {OUT}/")
