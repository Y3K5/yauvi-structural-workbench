"""memorient.barrel — membrane fit (find the bilayer) + barrel/surface classifier.

The membrane placement — unit normal ``n``, centre offset ``c`` along the normal, and
half-thickness ``d`` — is found by **maximizing the physical membrane signature**, not by
taking a principal axis. A plain "first SVD component = barrel axis" heuristic collapses
whenever a large extramembrane domain (a BamA POTRA arm, a TonB plug) dominates the CA
variance; the global scan below is immune to that.

Objective (maximized)::

    J(n, c, d) =  1.00 * dKD(lipid - pore)     # hydrophobic outside, polar lumen
                + 0.70 * belt_contrast          # embedded residues more hydrophobic than flanks
                + 0.50 * aromatic_girdle         # Trp/Tyr enrichment at the two interfaces
                - 0.25 * thickness_prior(d)      # soft half-thickness prior (never a hard clamp)

Search
------
1. **Global scan** of ~80 candidate normals on a Fibonacci hemisphere, plus the three PCA
   axes and a hydrophobic-moment seed.
2. Per normal, a cheap 1-D inner search over ``(c, d)``.
3. **Nelder-Mead polish** of the top-6 seeds (scipy) over the full ``(theta, phi, c, d)``.
4. A final **embedded-residues-only re-fit**: recompute the normal from just the residues
   inside the fitted slab (a PCA of the membrane-embedded CA), so the axis is derived from
   the membrane strands alone rather than the whole protein.

The context supplies the thickness prior and *which metric terms are active* (a single-pass
helix has no lumen, so ``lipid_pore_gap`` is dropped for it — see :mod:`memorient.contexts`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

from .contexts import MembraneContext, Metric

# Kyte-Doolittle hydropathy (higher = more hydrophobic).
KD = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3,
    "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
    "LYS": -3.9, "ARG": -4.5,
}
AROMATIC_INTERFACE = {"TRP", "TYR"}  # the aromatic girdle residues


def _kd(resnames: np.ndarray) -> np.ndarray:
    return np.array([KD.get(str(r), 0.0) for r in resnames])


def _fibonacci_hemisphere(n: int) -> np.ndarray:
    """``n`` roughly-uniform unit vectors on the upper hemisphere (z >= 0).

    A membrane normal and its negation define the same plane, so we only need a hemisphere.
    """
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - i / n)  # z from 1 down to ~0
    golden = np.pi * (1.0 + 5.0 ** 0.5)
    theta = golden * i
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack([x, y, z])


@dataclass
class MembraneFit:
    """Result of fitting a bilayer slab to a structure."""

    normal: np.ndarray          # unit membrane normal
    center: float               # slab centre offset along the normal (from CA centroid)
    half_thickness: float       # slab half-thickness (Angstrom)
    centroid: np.ndarray        # CA centroid the offset is measured from
    score: float                # objective J at the optimum
    components: Dict[str, float]
    embedded_mask: np.ndarray   # (N,) residues inside the slab
    n_embedded: int
    inner_frac: float           # fraction of embedded residues near the axis (hollowness probe)
    delta_kd: float             # lipid-minus-pore hydrophobicity gap

    def depth(self, ca: np.ndarray) -> np.ndarray:
        """Signed distance of each CA from the slab centre along the normal."""
        return (ca - self.centroid) @ self.normal - self.center


def _project(ca: np.ndarray, centroid: np.ndarray, normal: np.ndarray) -> np.ndarray:
    return (ca - centroid) @ normal


def _radial_distances(ca: np.ndarray, centroid: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Perpendicular distance of each CA from the membrane axis (through centroid)."""
    rel = ca - centroid
    along = (rel @ normal)[:, None] * normal[None, :]
    perp = rel - along
    return np.linalg.norm(perp, axis=1)


def _score_placement(
    ca: np.ndarray, sc_vec: np.ndarray, kd: np.ndarray, arom: np.ndarray,
    centroid: np.ndarray, normal: np.ndarray, center: float, half_thick: float,
    ctx: MembraneContext,
) -> Tuple[float, Dict[str, float], np.ndarray]:
    """Evaluate objective J for a candidate placement. Returns (J, components, embedded_mask)."""
    proj = _project(ca, centroid, normal) - center
    embedded = np.abs(proj) <= half_thick
    n_emb = int(embedded.count_nonzero() if hasattr(embedded, "count_nonzero") else embedded.sum())
    comps = {"delta_kd": 0.0, "belt": 0.0, "girdle": 0.0, "thickness": 0.0}
    if n_emb < 6:
        return -10.0, comps, embedded

    radial = _radial_distances(ca, centroid, normal)
    # lipid- vs pore-facing among embedded residues: side-chain vector's radial component.
    rel = ca - centroid
    along = (rel @ normal)[:, None] * normal[None, :]
    perp = rel - along
    perp_norm = perp / (np.linalg.norm(perp, axis=1, keepdims=True) + 1e-9)
    facing_out = np.einsum("ij,ij->i", sc_vec, perp_norm)  # >0 lipid-facing, <0 pore-facing

    emb = embedded
    lipid = emb & (facing_out > 0.15)
    pore = emb & (facing_out < -0.15)

    # -- delta KD (lipid minus pore) -------------------------------------------------
    if Metric.LIPID_PORE_GAP in ctx.metrics and lipid.sum() >= 3 and pore.sum() >= 3:
        dkd = float(kd[lipid].mean() - kd[pore].mean())
    else:
        dkd = 0.0
    comps["delta_kd"] = dkd

    # -- hydrophobic belt contrast: embedded vs flanking -----------------------------
    flank = (~emb) & (np.abs(proj) <= half_thick + 8.0)
    if Metric.HYDROPHOBIC_BELT in ctx.metrics and emb.sum() >= 3 and flank.sum() >= 3:
        belt = float(kd[emb].mean() - kd[flank].mean())
    else:
        belt = 0.0
    comps["belt"] = belt

    # -- aromatic girdle: Trp/Tyr enrichment at the two interfaces -------------------
    if Metric.AROMATIC_GIRDLE in ctx.metrics and emb.sum() >= 3:
        interface = emb & (np.abs(proj) > half_thick - 3.5)
        core = emb & (np.abs(proj) <= half_thick - 3.5)
        f_int = arom[interface].mean() if interface.sum() else 0.0
        f_core = arom[core].mean() if core.sum() else 0.0
        girdle = float(f_int - f_core)
    else:
        girdle = 0.0
    comps["girdle"] = girdle

    # -- soft thickness prior --------------------------------------------------------
    tp = ctx.thickness_prior
    thick_pen = tp.penalty(half_thick) if tp is not None else 0.0
    comps["thickness"] = thick_pen

    J = 1.0 * dkd + 0.7 * belt + 0.5 * girdle - 0.25 * thick_pen
    return J, comps, embedded


def _inner_search_cd(
    ca, sc_vec, kd, arom, centroid, normal, ctx, half_range=(8.0, 20.0),
) -> Tuple[float, float, float]:
    """Cheap grid search over centre offset and half-thickness for a fixed normal."""
    proj = _project(ca, centroid, normal)
    best = (-1e9, 0.0, 13.0)
    tp = ctx.thickness_prior
    d0 = tp.mean if tp is not None else 13.0
    centers = np.linspace(proj.min() * 0.5, proj.max() * 0.5, 9)
    # allow centre near the projection median too
    centers = np.unique(np.concatenate([centers, [np.median(proj)]]))
    halfs = np.linspace(half_range[0], half_range[1], 7)
    for c in centers:
        for d in halfs:
            J, _, _ = _score_placement(ca, sc_vec, kd, arom, centroid, normal, c, d, ctx)
            if J > best[0]:
                best = (J, float(c), float(d))
    return best


def fit_membrane(structure, ctx: MembraneContext, n_scan: int = 80,
                 polish: bool = True) -> MembraneFit:
    """Fit a bilayer slab by maximizing the membrane signature (see module docstring)."""
    ca = structure.ca
    sc_vec = structure.sc_vec
    kd = _kd(structure.resnames)
    arom = np.array([1.0 if str(r) in AROMATIC_INTERFACE else 0.0 for r in structure.resnames])
    centroid = ca.mean(axis=0)

    # candidate normals: fibonacci hemisphere + PCA axes + hydrophobic-moment seed
    normals = list(_fibonacci_hemisphere(n_scan))
    from .geometry import principal_axes
    _, vecs, _ = principal_axes(ca)
    for i in range(3):
        normals.append(vecs[:, i])
    # hydrophobic-moment seed: vector from polar-centroid to hydrophobic-centroid
    w = kd - kd.min()
    if w.sum() > 0:
        hydro_center = (ca * w[:, None]).sum(axis=0) / w.sum()
        seed = hydro_center - centroid
        if np.linalg.norm(seed) > 1e-6:
            normals.append(seed / np.linalg.norm(seed))

    scored = []
    for n in normals:
        n = np.asarray(n, dtype=float)
        n = n / (np.linalg.norm(n) + 1e-12)
        J, c, d = _inner_search_cd(ca, sc_vec, kd, arom, centroid, n, ctx)
        scored.append((J, n, c, d))
    scored.sort(key=lambda t: t[0], reverse=True)

    best = scored[0]
    if polish:
        best = _polish(ca, sc_vec, kd, arom, centroid, ctx, scored[:6])

    J, n, c, d = best
    # embedded-residues-only re-fit: re-derive the normal from the membrane strands alone
    n, c, d, J = _embedded_refit(ca, sc_vec, kd, arom, centroid, ctx, n, c, d)

    Jfinal, comps, embedded = _score_placement(ca, sc_vec, kd, arom, centroid, n, c, d, ctx)
    # hollowness: fraction of embedded residues within half the barrel radius of the axis
    radial = _radial_distances(ca, centroid, n)
    emb_radial = radial[embedded]
    if len(emb_radial):
        rmax = np.percentile(emb_radial, 90)
        inner_frac = float(np.mean(emb_radial < 0.5 * rmax)) if rmax > 0 else 1.0
    else:
        inner_frac = 1.0
    return MembraneFit(
        normal=n, center=c, half_thickness=d, centroid=centroid, score=Jfinal,
        components=comps, embedded_mask=embedded, n_embedded=int(embedded.sum()),
        inner_frac=inner_frac, delta_kd=comps["delta_kd"],
    )


def fit_membrane_on_normal(
    structure,
    ctx: MembraneContext,
    normal: np.ndarray,
    *,
    half_range: tuple[float, float] = (8.0, 20.0),
) -> MembraneFit:
    """Fit slab centre and thickness while keeping a declared normal fixed.

    Alpha-helical membrane placement uses a normal derived from explicitly mapped
    transmembrane helices.  Reusing :func:`fit_membrane` would allow the barrel-oriented
    whole-structure search to replace that evidence, so this public helper optimizes only
    the remaining centre/thickness dimensions.
    """
    ca = structure.ca
    normal = np.asarray(normal, dtype=float).reshape(3)
    norm = float(np.linalg.norm(normal))
    if norm <= 1e-12:
        raise ValueError("membrane normal cannot be zero")
    normal = normal / norm
    centroid = ca.mean(axis=0)
    kd = _kd(structure.resnames)
    arom = np.array([1.0 if str(r) in AROMATIC_INTERFACE else 0.0 for r in structure.resnames])
    score, center, half_thickness = _inner_search_cd(
        ca, structure.sc_vec, kd, arom, centroid, normal, ctx, half_range=half_range,
    )
    final, components, embedded = _score_placement(
        ca, structure.sc_vec, kd, arom, centroid, normal, center, half_thickness, ctx,
    )
    radial = _radial_distances(ca, centroid, normal)
    embedded_radial = radial[embedded]
    if len(embedded_radial):
        rmax = float(np.percentile(embedded_radial, 90))
        inner_frac = float(np.mean(embedded_radial < 0.5 * rmax)) if rmax > 0 else 1.0
    else:
        inner_frac = 1.0
    return MembraneFit(
        normal=normal,
        center=float(center),
        half_thickness=float(half_thickness),
        centroid=centroid,
        score=float(final if np.isfinite(final) else score),
        components=components,
        embedded_mask=embedded,
        n_embedded=int(embedded.sum()),
        inner_frac=inner_frac,
        delta_kd=float(components["delta_kd"]),
    )


def _polish(ca, sc_vec, kd, arom, centroid, ctx, seeds):
    """Nelder-Mead polish over (theta, phi, c, d) for each seed; keep the best."""
    from scipy.optimize import minimize

    def unpack(x):
        theta, phi, c, d = x
        n = np.array([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])
        return n, c, max(d, 4.0)

    best = seeds[0]
    for J0, n0, c0, d0 in seeds:
        theta0 = np.arctan2(n0[1], n0[0])
        phi0 = np.arccos(np.clip(n0[2], -1, 1))
        x0 = np.array([theta0, phi0, c0, d0])

        def neg(x):
            n, c, d = unpack(x)
            nn = n / (np.linalg.norm(n) + 1e-12)
            J, _, _ = _score_placement(ca, sc_vec, kd, arom, centroid, nn, c, d, ctx)
            return -J

        res = minimize(neg, x0, method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 400})
        n, c, d = unpack(res.x)
        n = n / (np.linalg.norm(n) + 1e-12)
        J = -res.fun
        if J > best[0]:
            best = (J, n, c, d)
    return best


def _embedded_refit(ca, sc_vec, kd, arom, centroid, ctx, n, c, d, iters: int = 3):
    """Re-derive the normal from the membrane-embedded residues alone, a few times."""
    from .geometry import principal_axes

    best = (n, c, d, _score_placement(ca, sc_vec, kd, arom, centroid, n, c, d, ctx)[0])
    for _ in range(iters):
        proj = _project(ca, centroid, n) - c
        embedded = np.abs(proj) <= d
        if embedded.sum() < 8:
            break
        # The barrel strands run along the normal, so the *smallest*-variance in-plane
        # direction of the embedded shell is noise; the normal is the axis along which the
        # embedded CAs are most extended relative to their radial spread. Use the embedded
        # covariance and pick the eigenvector best aligned with the current normal.
        _, vecs, vals = principal_axes(ca[embedded])
        # choose eigenvector most aligned with current normal (sign-agnostic)
        aligns = np.abs(vecs.T @ n)
        pick = vecs[:, int(np.argmax(aligns))]
        if np.dot(pick, n) < 0:
            pick = -pick
        # blend to avoid oscillation, renormalize
        n_new = pick / (np.linalg.norm(pick) + 1e-12)
        J_new, c_new, d_new = _inner_search_cd(ca, sc_vec, kd, arom, centroid, n_new, ctx)
        if J_new > best[3]:
            best = (n_new, c_new, d_new, J_new)
            n, c, d = n_new, c_new, d_new
        else:
            break
    return best[0], best[1], best[2], best[3]


# --------------------------------------------------------------------------------------
# Classifier
# --------------------------------------------------------------------------------------


@dataclass
class MembraneClass:
    label: str            # "barrel" | "surface" | "uncertain"
    confidence: float
    reasons: Dict[str, float]
    fit: Optional[MembraneFit]


def classify_membrane_protein(
    structure, ctx: MembraneContext, fit: Optional[MembraneFit] = None,
    dkd_barrel: float = 1.0, dkd_surface: float = 0.6,
    inner_frac_max: float = 0.25, n_embedded_min: int = 40,
) -> MembraneClass:
    """Classify barrel vs surface from the fit's own metrics.

    * ``delta_kd > dkd_barrel`` AND hollow cross-section (``inner_frac < inner_frac_max``)
      AND enough embedded residues (``n_embedded >= n_embedded_min``) -> **barrel**
    * ``delta_kd < dkd_surface`` OR a filled centre OR too few embedded -> **surface**
    * otherwise **uncertain**
    """
    if fit is None:
        fit = fit_membrane(structure, ctx)
    dkd = fit.delta_kd
    inner = fit.inner_frac
    n_emb = fit.n_embedded
    reasons = {"delta_kd": dkd, "inner_frac": inner, "n_embedded": float(n_emb)}

    hollow = inner < inner_frac_max
    enough = n_emb >= n_embedded_min

    if dkd > dkd_barrel and hollow and enough:
        conf = min(1.0, 0.5 + 0.5 * (dkd - dkd_barrel) / max(dkd_barrel, 1e-6))
        return MembraneClass("barrel", conf, reasons, fit)
    if dkd < dkd_surface or not hollow or not enough:
        # confidence higher the more clearly it fails the barrel test
        conf = 0.6
        if dkd < dkd_surface:
            conf = min(1.0, 0.6 + 0.4 * (dkd_surface - dkd))
        return MembraneClass("surface", conf, reasons, fit)
    return MembraneClass("uncertain", 0.4, reasons, fit)
