"""memorient.sasa — in-package Shrake-Rupley solvent-accessible surface area and RSA.

No external SASA binary. The Shrake-Rupley (1973) algorithm rolls a probe sphere over the
van der Waals surface: each atom is inflated by its VdW radius + the probe radius, a
Fibonacci sphere of test points is placed on that inflated surface, and a test point counts
as accessible if it lies outside every *other* inflated atom. Atom SASA is the accessible
fraction times the inflated-sphere area; residue SASA is the sum over its atoms.

Relative SASA (RSA) divides residue SASA by the residue's theoretical maximum accessible
area (Tien et al. 2013, *PLoS ONE*), giving a 0..1 exposure that is comparable across
residue types. RSA >= ~0.20 is the usual "exposed" threshold.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

# Bondi (1964) van der Waals radii (Angstrom) by element.
_BONDI = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "SE": 1.90, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
}
_DEFAULT_RADIUS = 1.70
PROBE_RADIUS = 1.40  # water

# Tien et al. 2013 theoretical maximum accessible surface area (Angstrom^2), per residue.
MAXASA_TIEN = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


def _fibonacci_sphere(n: int) -> np.ndarray:
    """``n`` roughly-uniform unit vectors on the sphere (golden-spiral construction)."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    golden = np.pi * (1.0 + 5.0 ** 0.5)
    theta = golden * i
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack([x, y, z])


def _atom_radius(element: str, name: str) -> float:
    el = (element or "").strip().upper()
    if el in _BONDI:
        return _BONDI[el]
    # fall back to first alpha char of the atom name
    for ch in name.upper():
        if ch.isalpha():
            return _BONDI.get(ch, _DEFAULT_RADIUS)
    return _DEFAULT_RADIUS


def atom_sasa(
    coords: np.ndarray,
    radii: np.ndarray,
    n_points: int = 240,
    probe: float = PROBE_RADIUS,
) -> np.ndarray:
    """Per-atom SASA (Angstrom^2) via Shrake-Rupley.

    Parameters
    ----------
    coords : (M,3) atom coordinates
    radii  : (M,) van der Waals radii (probe NOT yet added)
    n_points : test points per atom (240 is a good accuracy/speed tradeoff)
    """
    coords = np.asarray(coords, dtype=float).reshape(-1, 3)
    radii = np.asarray(radii, dtype=float).reshape(-1)
    M = len(coords)
    sphere = _fibonacci_sphere(n_points)
    inflated = radii + probe
    out = np.zeros(M)

    # neighbour cutoff: two atoms can occlude each other only if closer than the sum of
    # their inflated radii. Use a simple spatial hash on a grid of the max inflated diameter.
    max_r = float(inflated.max()) if M else 0.0
    cell = 2.0 * max_r + 1e-6
    grid: Dict[tuple, list] = {}
    keys = np.floor(coords / cell).astype(int)
    for idx in range(M):
        grid.setdefault(tuple(keys[idx]), []).append(idx)

    area_per_point = 4.0 * np.pi / n_points

    for i in range(M):
        ri = inflated[i]
        ci = coords[i]
        # candidate neighbours from the 27 surrounding cells
        base = keys[i]
        neigh = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neigh.extend(grid.get((base[0] + dx, base[1] + dy, base[2] + dz), []))
        neigh = [j for j in neigh if j != i]
        if neigh:
            nc = coords[neigh]
            nr = inflated[neigh]
            d = np.linalg.norm(nc - ci, axis=1)
            keep = d < (ri + nr)
            nc = nc[keep]
            nr = nr[keep]
        else:
            nc = np.empty((0, 3))
            nr = np.empty(0)
        pts = ci + ri * sphere  # (n_points, 3)
        if len(nc):
            # a point is buried if within any neighbour's inflated radius
            # dist2 (n_points, n_neigh)
            diff = pts[:, None, :] - nc[None, :, :]
            dist2 = np.einsum("ijk,ijk->ij", diff, diff)
            buried = np.any(dist2 < (nr[None, :] ** 2), axis=1)
            accessible = np.count_nonzero(~buried)
        else:
            accessible = n_points
        out[i] = area_per_point * accessible * ri * ri
    return out


def compute_sasa(structure, n_points: int = 240, probe: float = PROBE_RADIUS,
                 heavy_only: bool = True) -> Dict[str, np.ndarray]:
    """Per-residue SASA and RSA for a :class:`memorient.geometry.Structure`.

    Returns a dict with parallel arrays indexed like ``structure``:

    * ``sasa``      (N,) residue SASA (Angstrom^2)
    * ``rsa``       (N,) relative SASA in [0, ~1] via Tien 2013 MaxASA
    * ``atom_sasa`` (M,) per-atom SASA (all atoms in ``structure.atoms``)

    Uses all atoms in ``structure.atoms``. Hydrogens are dropped by default (``heavy_only``).
    """
    atoms = structure.atoms
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)
    elems = [a.get("element", "") for a in atoms]
    names = [a.get("name", "") for a in atoms]
    if heavy_only:
        mask = np.array([
            (e.strip().upper() not in ("H", "D")) for e in elems
        ])
    else:
        mask = np.ones(len(atoms), dtype=bool)
    radii = np.array([_atom_radius(elems[i], names[i]) for i in range(len(atoms))])

    asa_all = np.zeros(len(atoms))
    sub_idx = np.where(mask)[0]
    if len(sub_idx):
        asa_sub = atom_sasa(coords[sub_idx], radii[sub_idx], n_points=n_points, probe=probe)
        asa_all[sub_idx] = asa_sub

    # aggregate to residues, keyed by (chain, resseq, icode)
    res_sasa: Dict[tuple, float] = {}
    for a, sa in zip(atoms, asa_all):
        key = (a["chain"], a["resseq"], a.get("icode", ""))
        res_sasa[key] = res_sasa.get(key, 0.0) + sa

    N = len(structure)
    sasa = np.zeros(N)
    rsa = np.zeros(N)
    for i in range(N):
        # match residue arrays back to atom keys; icode not tracked in per-res arrays, so
        # match on (chain, resid) and take the residue's atom-summed SASA
        chain = structure.chains[i]
        resid = int(structure.resids[i])
        # sum across any icode variants for this (chain, resid)
        total = 0.0
        for (c, r, ic), val in res_sasa.items():
            if c == chain and r == resid:
                total += val
        sasa[i] = total
        maxasa = MAXASA_TIEN.get(str(structure.resnames[i]), 200.0)
        rsa[i] = total / maxasa if maxasa > 0 else 0.0

    return {"sasa": sasa, "rsa": rsa, "atom_sasa": asa_all}


def total_sasa(structure, **kw) -> float:
    """Total molecular SASA (Angstrom^2)."""
    return float(compute_sasa(structure, **kw)["sasa"].sum())

