#!/usr/bin/env python3
"""Pre-adoption screen for the ABL StateAtlas panel — separability and margin headroom.

Run before curating any record. Answers three questions with measurement rather
than argument:

  1. Does global CA RMSD over the frozen mask (UniProt P00519 242-495) separate
     ABL1 kinase-domain conformations at all?
  2. If it does, is the separation driven by the activation loop and the alphaC
     region — the biological switch — or by something incidental?
  3. Is the frozen margin (0.25 A) far from or close to the data, and what
     determines that?

State labels are deliberately NOT assumed. The two groups below come from
hierarchical clustering of the RMSD matrix, so this establishes separability and
threshold feasibility only. Which group is active and which inactive requires
independent state evidence per coverage rule 3.

Requires network (RCSB, PDBe, UniProt). Writes nothing into the panel manifest.
"""
from __future__ import annotations

import itertools
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
from Bio import Align
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

# Human ABL1 only. 1IEP, 1FPU, 1OPJ, 3KF4, 3KFA, 1M52 and 3OXZ are murine
# (P00520) and are excluded: mixing accessions breaks the residue-equivalence
# basis the mask depends on.
PDBS = ["1opl", "2hyy", "2gqg", "2g2i", "2g1t", "2e2b", "3qri", "3qrj",
        "4wa9", "6xr6", "6xrg", "2f4j", "3cs9", "5mo4", "4xey", "2fo0"]
ACCESSION = "P00519"
LO, HI = 242, 495          # frozen in state_atlas.core as ABL_DOMAIN_START/END
MAX_RMSD_A = 2.5           # frozen as ABL_MAX_RMSD_A
MIN_MARGIN_A = 0.25        # frozen as ABL_MIN_MARGIN_A


def fetch(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    for p in PDBS:
        cif, sifts = work / f"{p}.cif", work / f"{p}.sifts.json"
        if not cif.exists() or cif.stat().st_size < 1000:
            subprocess.run(["curl", "-sS", "--retry", "3", "--max-time", "90", "-o", str(cif),
                            f"https://files.rcsb.org/download/{p.upper()}.cif"], check=True)
        if not sifts.exists() or sifts.stat().st_size < 100:
            subprocess.run(["curl", "-sS", "--retry", "3", "--max-time", "60", "-o", str(sifts),
                            f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{p}"], check=True)
    fasta = work / f"{ACCESSION}.fasta"
    if not fasta.exists():
        subprocess.run(["curl", "-sS", "--retry", "3", "--max-time", "60", "-o", str(fasta),
                        f"https://rest.uniprot.org/uniprotkb/{ACCESSION}.fasta"], check=True)


def unp_sequence(work: Path) -> str:
    return "".join(l.strip() for l in (work / f"{ACCESSION}.fasta").read_text().splitlines()
                   if not l.startswith(">"))


def one_letter(res) -> str | None:
    info = gemmi.find_tabulated_residue(res.name)
    return info.one_letter_code.upper() if info and info.is_amino_acid() else None


def chain_coords(work: Path, pdb: str, unp: str) -> tuple[str, dict[int, np.ndarray], float]:
    """CA coordinates keyed by UniProt number, from the best-covered chain.

    Mapping is by global pairwise alignment of the observed sequence to the
    canonical UniProt sequence, matching what StateAtlas itself does. The PDBe
    residue-level SIFTS route was tried first and abandoned: `author_residue_number`
    is null for twelve of these sixteen entries, which silently drops them.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode, aligner.match_score, aligner.mismatch_score = "global", 2, -1
    aligner.open_gap_score, aligner.extend_gap_score = -5, -0.5
    st = gemmi.read_structure(str(work / f"{pdb}.cif"))
    st.setup_entities()
    best = None
    for ch in st[0]:
        obs = [(r, one_letter(r)) for r in ch
               if r.find_atom("CA", "*") is not None and one_letter(r)]
        if len(obs) < 100:
            continue
        seq = "".join(c for _, c in obs)
        aln = aligner.align(seq, unp)[0]
        pairs: dict[int, int] = {}
        for (qs, qe), (rs, _re) in zip(aln.aligned[0], aln.aligned[1]):
            for k in range(qe - qs):
                pairs[qs + k] = rs + k + 1          # UniProt numbering is 1-based
        ident = sum(1 for i, u in pairs.items() if seq[i] == unp[u - 1]) / max(1, len(pairs))
        coords = {u: np.array(r.find_atom("CA", "*").pos.tolist())
                  for i, (r, _) in enumerate(obs)
                  if (u := pairs.get(i)) and LO <= u <= HI}
        if best is None or len(coords) > len(best[1]):
            best = (ch.name, coords, ident)
    if best is None:
        raise SystemExit(f"no usable chain in {pdb}")
    return best


def kabsch(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ac, bc = a - a.mean(0), b - b.mean(0)
    v, _s, wt = np.linalg.svd(ac.T @ bc)
    d = np.diag([1, 1, np.sign(np.linalg.det(v @ wt))])
    return ac @ (v @ d @ wt), bc


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    ar, br = kabsch(a, b)
    return float(np.sqrt(np.mean(np.sum((ar - br) ** 2, axis=1))))


def main() -> int:
    work = Path(sys.argv[1] if len(sys.argv) > 1 else "abl_screen_work")
    fetch(work)
    unp = unp_sequence(work)
    span = HI - LO + 1

    data: dict[str, dict[int, np.ndarray]] = {}
    print(f"{'pdb':7}{'chain':6}{'mapped':>7}{'cover':>8}{'ident':>7}")
    for p in PDBS:
        ch, coords, ident = chain_coords(work, p, unp)
        data[p] = coords
        note = "" if len(coords) / span >= 0.90 else "   below the 0.90 coverage rule"
        print(f"{p:7}{ch:6}{len(coords):7}{len(coords)/span:8.3f}{ident:7.3f}{note}")

    order = sorted(data)
    n = len(order)
    m = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        common = sorted(set(data[order[i]]) & set(data[order[j]]))
        a = np.array([data[order[i]][u] for u in common])
        b = np.array([data[order[j]][u] for u in common])
        m[i, j] = m[j, i] = rmsd(a, b)

    print("\npairwise CA RMSD (A)\n")
    print("       " + "".join(f"{p:>7}" for p in order))
    for i, p in enumerate(order):
        print(f"{p:6} " + "".join("      ." if i == j else f"{m[i, j]:7.2f}" for j in range(n)))

    groups = fcluster(linkage(squareform(m, checks=False), method="average"), 2, "maxclust")
    g = {c: [order[i] for i in range(n) if groups[i] == c] for c in (1, 2)}
    print(f"\ntwo clusters (empirical, NOT state labels):\n  A {g[1]}\n  B {g[2]}")
    within = [m[i, j] for i, j in itertools.combinations(range(n), 2) if groups[i] == groups[j]]
    between = [m[i, j] for i, j in itertools.combinations(range(n), 2) if groups[i] != groups[j]]
    print(f"  within  median {np.median(within):.2f} A   between median {np.median(between):.2f} A")

    print("\nper-residue contribution to the between-cluster difference")
    acc: dict[int, list[float]] = {}
    for a_id in g[1]:
        for b_id in g[2]:
            common = sorted(set(data[a_id]) & set(data[b_id]))
            ar, br = kabsch(np.array([data[a_id][u] for u in common]),
                            np.array([data[b_id][u] for u in common]))
            for u, x in zip(common, np.linalg.norm(ar - br, axis=1)):
                acc.setdefault(u, []).append(float(x))
    prof = {u: float(np.mean(v)) for u, v in acc.items()}
    runs, cur = [], []
    for u in sorted(u for u, x in prof.items() if x >= 4.0):
        if cur and u == cur[-1] + 1:
            cur.append(u)
        else:
            if len(cur) >= 3:
                runs.append(cur)
            cur = [u]
    if len(cur) >= 3:
        runs.append(cur)
    for r in runs:
        print(f"  UniProt {r[0]}-{r[-1]}  {len(r)} residues, mean {np.mean([prof[u] for u in r]):.1f} A")

    idx = {p: i for i, p in enumerate(order)}
    for withheld in (None, "5mo4"):
        refs_a = [p for p in g[1] if p != withheld]
        refs_b = [p for p in g[2] if p != withheld]
        label = "all sixteen usable as references" if withheld is None else f"{withheld.upper()} withheld from references"
        print(f"\nmargin simulation — {label}")
        below = 0
        for p in order:
            a = min((m[idx[p], idx[q]] for q in refs_a if q != p), default=float("inf"))
            b = min((m[idx[p], idx[q]] for q in refs_b if q != p), default=float("inf"))
            margin = abs(a - b)
            if margin < MIN_MARGIN_A:
                below += 1
                print(f"  {p:7} margin {margin:.2f} A  UNRESOLVED")
        print(f"  below the {MIN_MARGIN_A} A margin: {below}/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
