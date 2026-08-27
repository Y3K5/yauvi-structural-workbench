"""memorient.membrane — context-selected bilayer zones, facing, and validation metrics.

Given an oriented structure and its :class:`~memorient.barrel.MembraneFit`, this module places
each residue in a **bilayer zone** along the membrane normal, decides whether an embedded
residue faces **lipid or the pore/water**, and folds zone + facing + SASA into a single
**membrane_accessibility** category that epitope mining reads.

Everything here is gated by the :class:`~memorient.contexts.MembraneContext`:

* the **bilayer model** is asymmetric (LPS outer leaflet + PL inner leaflet, with an LPS
  buffer that shields the proximal extracellular surface) for ``gram_negative_om`` and
  symmetric (two equivalent PL leaflets, no shield) for ``eukaryotic_pm`` / ``tm_receptor``;
* a **non-membrane context computes no zones at all** — a secreted protein is never assigned
  a membrane band or drawn inside a slab it does not have;
* :func:`context_metrics` evaluates **only the metrics the context declares active**, so the
  aromatic girdle / lipid-pore gap are never scored on a single TM helix that has no lumen.

The signed membrane coordinate is ``ec_depth = ec_sign * (proj - center)`` where ``ec_sign``
puts the extracellular side at positive depth. The extracellular side itself is decided by
:mod:`memorient.labeler`; here we take it as input.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

from .barrel import KD, MembraneFit
from .contexts import MembraneContext, MembraneModel, Metric

# Zone names along the membrane normal (extracellular at the top).
ZONE_EXTRACELLULAR = "extracellular"
ZONE_EC_INTERFACE = "extracellular_interface"   # headgroup / LPS band on the EC side
ZONE_CORE = "hydrophobic_core"
ZONE_PERI_INTERFACE = "periplasmic_interface"    # headgroup band on the periplasmic/cyto side
ZONE_PERIPLASMIC = "periplasmic"                 # or cytoplasmic for a plasma membrane
ZONES = (ZONE_EXTRACELLULAR, ZONE_EC_INTERFACE, ZONE_CORE, ZONE_PERI_INTERFACE, ZONE_PERIPLASMIC)

# Accessibility categories (what epitope mining is allowed to draw from).
ACC_ANTIBODY = "antibody_accessible"     # the only publishable category
ACC_LPS_SHIELDED = "lps_shielded"
ACC_LIPID = "lipid_embedded"
ACC_PORE = "pore_lumen_facing"
ACC_PERIPLASMIC = "periplasmic"
ACC_BURIED = "buried_interior"

INTERFACE_WIDTH = 4.0   # Angstrom band on each side of the core for headgroups
LPS_BUFFER = 6.0        # Angstrom above the EC interface that LPS shields (gram-neg OM only)
RSA_EXPOSED = 0.20      # relative SASA above which a residue is "exposed"


@dataclass
class MembraneProjection:
    """Per-residue membrane placement for an oriented structure."""

    ec_depth: np.ndarray        # (N,) signed depth; + = extracellular side
    zone: np.ndarray            # (N,) zone label (object)
    facing: np.ndarray          # (N,) "lipid" | "pore" | "" (object)
    accessibility: np.ndarray   # (N,) accessibility category (object)
    ec_sign: int                # +1 or -1: sign that maps proj-center to ec_depth

    def in_zone(self, zone: str) -> np.ndarray:
        return self.zone == zone


def _radial_unit(ca: np.ndarray, centroid: np.ndarray, normal: np.ndarray) -> np.ndarray:
    rel = ca - centroid
    along = (rel @ normal)[:, None] * normal[None, :]
    perp = rel - along
    return perp / (np.linalg.norm(perp, axis=1, keepdims=True) + 1e-9)


def project_membrane(
    structure, fit: MembraneFit, ctx: MembraneContext, ec_sign: int,
    rsa: Optional[np.ndarray] = None,
) -> MembraneProjection:
    """Assign zone, facing, and accessibility to every residue.

    Parameters
    ----------
    ec_sign : int
        +1 if the extracellular side is at positive ``(proj - center)``, else -1. Supplied by
        the labeler's extracellular-side call.
    rsa : (N,) array, optional
        Relative SASA; residues below :data:`RSA_EXPOSED` are treated as buried.
    """
    if not ctx.has_bilayer:
        raise ValueError(f"context {ctx.name!r} has no bilayer; do not project membrane zones")

    ca = structure.ca
    N = len(structure)
    proj = (ca - fit.centroid) @ fit.normal - fit.center
    ec_depth = ec_sign * proj
    d = fit.half_thickness

    zone = np.empty(N, dtype=object)
    facing = np.array([""] * N, dtype=object)
    acc = np.empty(N, dtype=object)

    core = np.abs(proj) <= d
    ec_side = ec_depth > 0
    # zones
    zone[core] = ZONE_CORE
    ec_iface = (~core) & (ec_depth > 0) & (ec_depth <= d + INTERFACE_WIDTH)
    peri_iface = (~core) & (ec_depth < 0) & (ec_depth >= -(d + INTERFACE_WIDTH))
    zone[ec_iface] = ZONE_EC_INTERFACE
    zone[peri_iface] = ZONE_PERI_INTERFACE
    zone[(~core) & (ec_depth > d + INTERFACE_WIDTH)] = ZONE_EXTRACELLULAR
    zone[(~core) & (ec_depth < -(d + INTERFACE_WIDTH))] = ZONE_PERIPLASMIC

    # facing for core residues: side-chain radial component
    perp_unit = _radial_unit(ca, fit.centroid, fit.normal)
    facing_out = np.einsum("ij,ij->i", structure.sc_vec, perp_unit)
    facing[core & (facing_out > 0.15)] = "lipid"
    facing[core & (facing_out < -0.15)] = "pore"

    # accessibility
    if rsa is None:
        rsa = np.ones(N)
    exposed = rsa >= RSA_EXPOSED

    for i in range(N):
        z = zone[i]
        if z == ZONE_CORE:
            acc[i] = ACC_PORE if facing[i] == "pore" else ACC_LIPID
        elif z in (ZONE_PERIPLASMIC, ZONE_PERI_INTERFACE):
            acc[i] = ACC_PERIPLASMIC
        elif z == ZONE_EC_INTERFACE:
            # LPS band shields this on a gram-negative OM
            if ctx.lps_shielding:
                acc[i] = ACC_LPS_SHIELDED
            else:
                acc[i] = ACC_ANTIBODY if exposed[i] else ACC_BURIED
        else:  # extracellular
            # LPS extends a short buffer above the interface on gram-neg OM
            if ctx.lps_shielding and ec_depth[i] <= d + INTERFACE_WIDTH + LPS_BUFFER:
                acc[i] = ACC_LPS_SHIELDED
            elif exposed[i]:
                acc[i] = ACC_ANTIBODY
            else:
                acc[i] = ACC_BURIED

    return MembraneProjection(ec_depth=ec_depth, zone=zone, facing=facing,
                              accessibility=acc, ec_sign=ec_sign)


# --------------------------------------------------------------------------------------
# Context-selected validation metrics
# --------------------------------------------------------------------------------------


def context_metrics(
    structure, fit: MembraneFit, ctx: MembraneContext, proj: MembraneProjection,
) -> Dict[str, float]:
    """Evaluate ONLY the metrics the context declares active (the gating contract).

    A metric that is not in ``ctx.active_metrics()`` is absent from the returned dict — never
    computed, never reported. Running an inactive metric would measure noise.
    """
    ca = structure.ca
    kd = np.array([KD.get(str(r), 0.0) for r in structure.resnames])
    d = fit.half_thickness
    depth = proj.ec_depth
    core = np.abs((ca - fit.centroid) @ fit.normal - fit.center) <= d
    out: Dict[str, float] = {}

    if Metric.LIPID_PORE_GAP in ctx.metrics:
        lipid = core & (proj.facing == "lipid")
        pore = core & (proj.facing == "pore")
        if lipid.sum() >= 3 and pore.sum() >= 3:
            out[Metric.LIPID_PORE_GAP] = float(kd[lipid].mean() - kd[pore].mean())
        else:
            out[Metric.LIPID_PORE_GAP] = float("nan")

    if Metric.HYDROPHOBIC_BELT in ctx.metrics:
        flank = (~core) & (np.abs(depth) <= d + 8.0)
        if core.sum() >= 3 and flank.sum() >= 3:
            out[Metric.HYDROPHOBIC_BELT] = float(kd[core].mean() - kd[flank].mean())
        else:
            out[Metric.HYDROPHOBIC_BELT] = float("nan")

    if Metric.AROMATIC_GIRDLE in ctx.metrics:
        arom = np.array([1.0 if str(r) in ("TRP", "TYR") else 0.0 for r in structure.resnames])
        interface = core & (np.abs((ca - fit.centroid) @ fit.normal - fit.center) > d - 3.5)
        deep = core & ~interface
        f_int = arom[interface].mean() if interface.sum() else 0.0
        f_deep = arom[deep].mean() if deep.sum() else 0.0
        out[Metric.AROMATIC_GIRDLE] = float(f_int - f_deep)

    if Metric.POSITIVE_INSIDE in ctx.metrics:
        # cytoplasmic (periplasmic-side here) loops enriched in Lys/Arg vs extracellular loops
        pos = np.array([1.0 if str(r) in ("LYS", "ARG") else 0.0 for r in structure.resnames])
        ec_loops = depth > d + 3.0
        cyto_loops = depth < -(d + 3.0)
        f_cyto = pos[cyto_loops].mean() if cyto_loops.sum() else 0.0
        f_ec = pos[ec_loops].mean() if ec_loops.sum() else 0.0
        # positive value => more Lys/Arg on the cytoplasmic side (correct sign)
        out[Metric.POSITIVE_INSIDE] = float(f_cyto - f_ec)

    if Metric.LPS_SHIELDING in ctx.metrics:
        # fraction of would-be extracellular-exposed residues that LPS actually shields
        shielded = proj.accessibility == ACC_LPS_SHIELDED
        ec_exposed = depth > d
        denom = int(ec_exposed.sum())
        out[Metric.LPS_SHIELDING] = float(shielded.sum() / denom) if denom else 0.0

    return out

