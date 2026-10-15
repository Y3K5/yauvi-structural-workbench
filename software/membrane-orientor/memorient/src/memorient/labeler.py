"""memorient.labeler — extracellular-side call, per-residue topology, and the surface set.

Two jobs:

1. **Which membrane face is extracellular.** For a gram-negative OMP barrel the extracellular
   loops are long and irregular while the periplasmic turns are short. We score each face by
   the **median connecting-segment length + the long-loop fraction** (deliberately NOT the
   mean, so one giant soluble domain — a POTRA arm, a TonB plug — cannot invert the call),
   and take orthogonal evidence from **per-strand up/down connectivity** (a proper barrel
   alternates the direction of successive membrane-spanning segments; the face where more
   segments *terminate* carries the loops). These two signals **vote**; the positive-inside
   charge bias (fewer Arg/Lys on the extracellular side) is used **only as a weak tiebreaker**
   when the primary signals are close. The winning ``ec_sign`` and the per-signal votes are
   returned so a reviewer can see which evidence agreed.

2. **Per-residue topology + the antibody-accessible surface set.** Combines the membrane
   projection (zone / facing / accessibility) with SASA to emit, for every residue, a label
   and an ``extracellular`` boolean; the residues that are ``antibody_accessible`` form the
   surface set that epitope mining is allowed to draw from.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .barrel import MembraneFit
from .contexts import MembraneContext
from .membrane import (
    ACC_ANTIBODY,
    MembraneProjection,
    RSA_EXPOSED,
)

LONG_LOOP = 6.0  # residues; a connecting segment longer than this is a "long loop"


@dataclass
class SideCall:
    """Result of the extracellular-side vote."""

    ec_sign: int                      # +1 or -1: sign mapping (proj - center) -> ec_depth
    votes: Dict[str, int]             # per-signal vote (+1 / -1 / 0)
    scores: Dict[str, float]          # raw signal values
    agreement: float                  # fraction of non-abstaining signals that agreed
    confidence: float


def _segment_membrane_crossings(proj_along: np.ndarray, half_thick: float, center: float):
    """Identify membrane-spanning segments and the connecting loops between them.

    Returns (crossings, loops) where each crossing is (start_idx, end_idx, direction) with
    direction +1 for a segment running toward +proj and -1 toward -proj, and each loop is
    (start_idx, end_idx, side) with side = +1 if the loop sits on the +proj face.
    """
    z = proj_along - center
    inside = np.abs(z) <= half_thick
    crossings = []
    i = 0
    n = len(z)
    while i < n:
        if inside[i]:
            j = i
            while j < n and inside[j]:
                j += 1
            # segment i..j-1
            direction = 1 if z[min(j - 1, n - 1)] >= z[i] else -1
            crossings.append((i, j - 1, direction))
            i = j
        else:
            i += 1
    # loops = contiguous outside runs bounded by membrane on both sides (or termini)
    loops = []
    i = 0
    while i < n:
        if not inside[i]:
            j = i
            while j < n and not inside[j]:
                j += 1
            side = 1 if z[i:j].mean() > 0 else -1
            loops.append((i, j - 1, side, j - i))
            i = j
        else:
            i += 1
    return crossings, loops


def call_extracellular_side(
    structure, fit: MembraneFit, ctx: MembraneContext,
) -> SideCall:
    """Decide which membrane face is extracellular by voting loop architecture + connectivity.

    ``ec_sign = +1`` means the extracellular side is at positive ``(proj - center)``.
    """
    proj_along = (structure.ca - fit.centroid) @ fit.normal
    crossings, loops = _segment_membrane_crossings(proj_along, fit.half_thickness, fit.center)

    votes: Dict[str, int] = {}
    scores: Dict[str, float] = {}

    # -- signal 1: loop architecture (median length + long-loop fraction), per face --------
    plus_lengths = [L for (_, _, side, L) in loops if side > 0]
    minus_lengths = [L for (_, _, side, L) in loops if side < 0]

    def face_score(lengths: List[int]) -> float:
        if not lengths:
            return 0.0
        arr = np.array(lengths, dtype=float)
        median = float(np.median(arr))
        long_frac = float(np.mean(arr > LONG_LOOP))
        return median + 8.0 * long_frac  # long-loop fraction weighted to matter

    s_plus = face_score(plus_lengths)
    s_minus = face_score(minus_lengths)
    scores["loop_plus"] = s_plus
    scores["loop_minus"] = s_minus
    if abs(s_plus - s_minus) < 1e-6:
        votes["loop_architecture"] = 0
    else:
        votes["loop_architecture"] = 1 if s_plus > s_minus else -1

    # -- signal 2: terminus topology — classic OMP has both termini periplasmic ------------
    # A membrane-spanning barrel with an even strand count places both the N- and C-terminus
    # on the periplasmic face. The face the termini sit on is therefore periplasmic, so the
    # OTHER face is extracellular. Orthogonal to loop length; can be fooled by a large
    # terminal extension, so it is voted (weight 0.7), never decisive.
    z_all = proj_along - fit.center
    term_z = np.array([z_all[0], z_all[-1]])
    # only informative if the termini are actually off the core (out of the membrane)
    off = np.abs(term_z) > fit.half_thickness
    if off.any():
        term_side = np.sign(term_z[off].mean())  # +1 if termini on +proj face
        scores["terminus_side"] = float(term_side)
        # termini periplasmic -> extracellular is the opposite face
        votes["terminus"] = int(-term_side) if term_side != 0 else 0
    else:
        scores["terminus_side"] = 0.0
        votes["terminus"] = 0

    # -- signal 3 (weak tiebreaker): positive-inside — fewer Arg/Lys on extracellular side --
    pos = np.array([1.0 if str(r) in ("LYS", "ARG") else 0.0 for r in structure.resnames])
    z = proj_along - fit.center
    plus_mask = z > fit.half_thickness
    minus_mask = z < -fit.half_thickness
    f_plus = pos[plus_mask].mean() if plus_mask.sum() else 0.0
    f_minus = pos[minus_mask].mean() if minus_mask.sum() else 0.0
    scores["poscharge_plus"] = float(f_plus)
    scores["poscharge_minus"] = float(f_minus)
    # extracellular side has FEWER positive charges -> vote toward the lower-charge face
    if abs(f_plus - f_minus) < 1e-6:
        votes["positive_inside"] = 0
    else:
        votes["positive_inside"] = 1 if f_plus < f_minus else -1

    # -- combine: loop architecture is primary (weight 1.0); terminus topology is strong ---
    #    orthogonal evidence (0.7); positive-inside is only a weak tiebreaker (0.25).
    w = {"loop_architecture": 1.0, "terminus": 0.7, "positive_inside": 0.25}
    tally = sum(w[k] * v for k, v in votes.items())
    ec_sign = 1 if tally >= 0 else -1
    # degenerate all-zero -> default +1 but flag low confidence
    nonzero = [k for k, v in votes.items() if v != 0]
    if nonzero:
        agree = np.mean([1.0 if np.sign(votes[k]) == ec_sign else 0.0 for k in nonzero])
    else:
        agree = 0.0
        ec_sign = 1
    # confidence: how decisive the primary signals were
    primary = abs(votes.get("loop_architecture", 0)) + abs(votes.get("terminus", 0))
    conf = 0.4 + 0.3 * min(primary, 2) / 2 + 0.3 * agree
    return SideCall(ec_sign=ec_sign, votes=votes, scores=scores,
                    agreement=float(agree), confidence=float(min(conf, 1.0)))


# --------------------------------------------------------------------------------------
# Per-residue topology + surface set
# --------------------------------------------------------------------------------------


@dataclass
class ResidueLabel:
    resid: int
    resname: str
    chain: str
    zone: str
    facing: str
    accessibility: str
    extracellular: bool
    rsa: float
    ec_depth: float
    insertion_code: str = ""


@dataclass
class LabelSet:
    labels: List[ResidueLabel]
    surface_set: List[int] = field(default_factory=list)  # resids that are antibody_accessible

    def extracellular_resids(self) -> List[int]:
        return [l.resid for l in self.labels if l.extracellular]

    def to_rows(self) -> List[dict]:
        return [
            {
                "resid": l.resid, "insertion_code": l.insertion_code,
                "resname": l.resname, "chain": l.chain, "zone": l.zone,
                "facing": l.facing, "accessibility": l.accessibility,
                "extracellular": l.extracellular, "rsa": round(l.rsa, 3),
                "ec_depth": round(l.ec_depth, 2),
            }
            for l in self.labels
        ]


def label_residues(
    structure, proj: MembraneProjection, rsa: np.ndarray, ctx: MembraneContext,
    fit: Optional[MembraneFit] = None,
) -> LabelSet:
    """Emit per-residue labels + the antibody-accessible surface set.

    A residue is ``extracellular`` when it is **confidently** on the extracellular side and
    exposed: zone ``extracellular`` (not the headgroup interface band, which is buried in the
    LPS/lipid headgroups and not an antibody target) with an RSA margin above the exposure
    threshold. The margin matters because residues sitting exactly on the RSA cutoff or a zone
    boundary would otherwise flicker in/out under Shrake-Rupley sampling noise, breaking the
    rotation-invariance guarantee for no biological reason. The surface set is the residues
    whose accessibility is ``antibody_accessible`` — the only category epitope mining may
    draw from.
    """
    N = len(structure)
    labels: List[ResidueLabel] = []
    surface: List[int] = []
    from .membrane import ZONE_EXTRACELLULAR
    rsa_margin = 0.02  # hysteresis band around the exposure cutoff
    # A residue is confidently extracellular when it lies in the extracellular zone proper —
    # NOT the headgroup interface band, which is LPS/lipid-buried and not an antibody target —
    # and is exposed with a small RSA margin. Both criteria use the discrete zone label (whose
    # boundary the projection places from the fit) plus an RSA hysteresis band, so residues
    # clustered near the exposure cutoff don't flicker in and out under Shrake-Rupley sampling
    # noise. Zone membership tracks the fit, so a couple-Angstrom half-thickness wobble moves
    # the interface boundary and the residues' zone assignment together, keeping the set stable
    # across input frames.
    for i in range(N):
        zone = str(proj.zone[i])
        acc = str(proj.accessibility[i])
        exposed = rsa[i] >= RSA_EXPOSED
        # confident extracellular: clearly exposed AND in the extracellular zone proper
        is_ec = (rsa[i] >= RSA_EXPOSED + rsa_margin) and (zone == ZONE_EXTRACELLULAR)
        lbl = ResidueLabel(
            resid=int(structure.resids[i]), resname=str(structure.resnames[i]),
            chain=str(structure.chains[i]), zone=zone, facing=str(proj.facing[i]),
            accessibility=acc, extracellular=bool(is_ec), rsa=float(rsa[i]),
            ec_depth=float(proj.ec_depth[i]),
            insertion_code=str(structure.icodes[i]),
        )
        labels.append(lbl)
        if acc == ACC_ANTIBODY:
            surface.append(lbl.resid)
    return LabelSet(labels=labels, surface_set=surface)
