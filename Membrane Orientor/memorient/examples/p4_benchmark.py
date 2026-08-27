#!/usr/bin/env python
"""P4 correctness benchmark — accuracy, not just self-consistency.

We take gram-negative outer-membrane proteins from OPM (Orientations of Proteins in
Membranes), where the deposited coordinates are already oriented so the membrane normal is
+Z. We then **un-orient** each structure with a set of random rotations, re-fit it with
memorient, and measure the angle between the recovered membrane normal and the known OPM
normal. This tests correctness against an external experimental reference — distinct from the
rotation-invariance (self-consistency) check in the test suite.

Also compares memorient's fitted bilayer half-thickness against OPM's reported value.

Run:  python examples/p4_benchmark.py            # regenerate examples/p4_results.json + figure
Needs network (downloads from OPM) and the compute extras (numpy/scipy/biopython/matplotlib).
"""

from __future__ import annotations

import json
import os
import re
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tests"))

from memorient.barrel import fit_membrane
from memorient.contexts import get_context
from memorient.geometry import load_structure
from memorient.orientor import five_fold_validate

# PDB id -> common name. All are gram-negative OMP beta-barrels in OPM.
PANEL = {
    "1bxw": "OmpA",
    "2por": "Porin",
    "1qd6": "OmpF-family",
    "2f1t": "NspA",
    "1p4t": "OmpLA/PldA",
}
OPM_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb"
HERE = os.path.dirname(__file__)


def _random_rotation(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def _opm_half_thickness(path: str):
    for ln in open(path):
        if "1/2 of bilayer thickness" in ln:
            return float(re.search(r"([\d.]+)", ln.split(":")[1]).group(1))
    return None


def run(cache_dir: str = "/tmp") -> list:
    ctx = get_context("gram_negative_om")
    rows = []
    for pdb, name in PANEL.items():
        path = os.path.join(cache_dir, f"opm_{pdb}.pdb")
        if not os.path.exists(path):
            urllib.request.urlretrieve(OPM_URL.format(pdb=pdb), path)
        s = load_structure(path)                 # hardened loader tolerates OPM headers
        errs, ds = [], []
        for seed in range(4):
            R = _random_rotation(seed)
            fit = fit_membrane(s.transformed(R), ctx)
            true_n = R @ np.array([0.0, 0.0, 1.0])
            cos = min(abs(float(np.dot(fit.normal, true_n))), 1.0)
            errs.append(np.degrees(np.arccos(cos)))
            ds.append(fit.half_thickness)
        v = five_fold_validate(s, ctx, n_points=200)
        rows.append({
            "pdb": pdb.upper(), "name": name,
            "opm_half": _opm_half_thickness(path),
            "fit_half": round(float(np.mean(ds)), 2),
            "angle_err": round(float(np.mean(errs)), 2),
            "angle_sd": round(float(np.std(errs)), 2),
            "mean_jaccard": v["mean_jaccard"], "n_res": len(s),
        })
    return rows


def make_figure(rows: list, out_png: str) -> None:
    import matplotlib.pyplot as plt

    labels = [f"{r['pdb']}\n{r['name']}" for r in rows]
    ang = [r["angle_err"] for r in rows]
    angsd = [r["angle_sd"] for r in rows]
    opm = [r["opm_half"] for r in rows]
    fit = [r["fit_half"] for r in rows]
    x = np.arange(len(rows))
    mean_ang = float(np.mean(ang))
    FOCAL = "#2166ac"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ax1.errorbar(x, ang, yerr=angsd, fmt="o", color=FOCAL, ms=7, capsize=3, lw=1.2, zorder=3)
    ax1.axhline(15.0, ls="--", lw=1.0, color="#999999", zorder=1)
    ax1.text(len(rows) - 1, 15.6, "15\u00b0 \u2014 well-oriented", ha="right", va="bottom", fontsize=6, color="#666666")
    ax1.axhline(mean_ang, ls=":", lw=1.0, color=FOCAL, zorder=1)
    ax1.text(0.02, mean_ang + 0.4, f"mean {mean_ang:.1f}\u00b0", ha="left", va="bottom", fontsize=6, color=FOCAL)
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=6)
    ax1.set_ylabel("membrane-normal error vs OPM (\u00b0)")
    ax1.set_title("Recovered normal matches the experimental orientation", fontsize=8)
    ax1.set_ylim(-1, 22)
    ax1.text(0.99, 0.02, "lower = better", transform=ax1.transAxes, ha="right", va="bottom", fontsize=6, color="#666666")

    ax2.plot([9, 15], [9, 15], ls="--", lw=1.0, color="#999999", zorder=1)
    ax2.scatter(opm, fit, color=FOCAL, s=42, zorder=3)
    for xi, yi, r in zip(opm, fit, rows):
        ax2.annotate(r["pdb"], (xi, yi), textcoords="offset points", xytext=(5, 3), fontsize=6, color="#333333")
    ax2.set_xlabel("OPM half-thickness (\u00c5)")
    ax2.set_ylabel("memorient fitted half-thickness (\u00c5)")
    ax2.set_title("Fitted bilayer thickness tracks OPM", fontsize=8)
    ax2.set_xlim(10.5, 14); ax2.set_ylim(10.5, 14); ax2.set_aspect("equal")

    fig.suptitle("P4 correctness benchmark \u2014 5 gram-negative OMPs, coordinates un-oriented then re-fit", fontsize=8.5, y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    rows = run()
    with open(os.path.join(HERE, "p4_results.json"), "w") as fh:
        json.dump(rows, fh, indent=2)
    make_figure(rows, os.path.join(HERE, "p4_benchmark.png"))
    print(f"mean normal angular error: {np.mean([r['angle_err'] for r in rows]):.2f} deg")
    for r in rows:
        print(f"  {r['pdb']:6} {r['name']:14} err={r['angle_err']:5.1f}deg  "
              f"fit_d={r['fit_half']:.1f}  OPM_d={r['opm_half']}  J={r['mean_jaccard']}")
