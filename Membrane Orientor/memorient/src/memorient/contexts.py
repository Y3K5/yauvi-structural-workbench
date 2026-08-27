"""Membrane-context registry — the "metrics are not a blanket" core of memorient.

A beta-barrel in a gram-negative outer membrane, a single-pass receptor in a plasma
membrane, and a secreted S-layer antigen are three different physics problems. Running one
membrane model and one set of validation metrics on all of them produces confident nonsense.

Each :class:`MembraneContext` declares:

* the **membrane model** (asymmetric LPS / symmetric phospholipid / none),
* the **orientation method** the orientor should route to,
* a soft **thickness prior** (half-thickness mean/sd, Angstrom) — ``None`` when there is no bilayer,
* the set of **validation metrics that are meaningful** for that physics.

The orientor, labeler, and viz layers consult the context instead of running everything on
everything. This module is **stdlib-only** so the registry + CLI install and run with zero
third-party dependencies; the numeric layers import numpy/biopython lazily.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

__all__ = [
    "MembraneModel",
    "OrientationMethod",
    "Metric",
    "ThicknessPrior",
    "MembraneContext",
    "REGISTRY",
    "get_context",
    "list_contexts",
    "default_context",
]

# --------------------------------------------------------------------------------------
# Controlled vocabularies. Plain string constants (not Enums) so the registry stays
# trivially serializable to JSON for the CLI and cross-tool handoff.
# --------------------------------------------------------------------------------------


class MembraneModel:
    """How the bilayer (if any) is modelled along the membrane normal."""

    ASYMMETRIC_LPS = "asymmetric_lps"          # gram-negative OM: LPS outer + PL inner leaflet
    SYMMETRIC_PHOSPHOLIPID = "symmetric_pl"    # plasma / inner membrane: two equivalent PL leaflets
    NONE = "none"                              # no bilayer (cell-wall-anchored or secreted)

    ALL = (ASYMMETRIC_LPS, SYMMETRIC_PHOSPHOLIPID, NONE)


class OrientationMethod:
    """Which orientation algorithm the orientor routes to for this context."""

    BARREL_NORMAL = "barrel_normal"        # maximize the beta-barrel signature (fit_membrane)
    TM_HELIX_BELT = "tm_helix_belt"        # hydrophobic-belt slab + positive-inside sign (alpha-helical)
    ANCHOR_RELATIVE = "anchor_relative"    # SVD axis + N-terminal membrane-proximal anchor sign
    SASA_ONLY = "sasa_only"                # no membrane frame; solvent accessibility only

    ALL = (BARREL_NORMAL, TM_HELIX_BELT, ANCHOR_RELATIVE, SASA_ONLY)


class Metric:
    """Self-validation / scoring metrics. A context lists only the ones that are meaningful.

    Running a metric outside the physics it was designed for measures noise, so the orientor
    calls :meth:`MembraneContext.active_metrics` and evaluates only those.
    """

    AROMATIC_GIRDLE = "aromatic_girdle"        # Trp/Tyr enrichment at the two interfaces (barrels & helices)
    LIPID_PORE_GAP = "lipid_pore_gap"          # hydrophobic lipid face vs polar lumen (barrels: has a lumen)
    HYDROPHOBIC_BELT = "hydrophobic_belt"      # embedded residues more hydrophobic than flanking (any TM protein)
    POSITIVE_INSIDE = "positive_inside"        # cytoplasmic loops enriched in Lys/Arg (alpha-helical sign break)
    LPS_SHIELDING = "lps_shielding"            # LPS leaflet hides the proximal extracellular surface (gram-neg OM only)
    ROTATION_INVARIANCE = "rotation_invariance"  # answer independent of input frame (ALWAYS meaningful)

    ALL = (
        AROMATIC_GIRDLE,
        LIPID_PORE_GAP,
        HYDROPHOBIC_BELT,
        POSITIVE_INSIDE,
        LPS_SHIELDING,
        ROTATION_INVARIANCE,
    )


@dataclass(frozen=True)
class ThicknessPrior:
    """Soft prior on the bilayer **half-thickness** (Angstrom).

    Used as a gentle penalty ``((d - mean) / sd)^2`` in the membrane fit, never a hard clamp,
    so an unusual bilayer can still win if the data strongly support it. ``None`` for
    non-membrane contexts (see :attr:`MembraneContext.thickness_prior`).
    """

    mean: float
    sd: float

    def penalty(self, half_thickness: float) -> float:
        """Dimensionless soft penalty for a candidate half-thickness."""
        z = (float(half_thickness) - self.mean) / self.sd
        return z * z


@dataclass(frozen=True)
class MembraneContext:
    """A biological setting that fixes the membrane physics and the meaningful metrics."""

    name: str
    description: str
    membrane_model: str
    orientation_method: str
    thickness_prior: Optional[ThicknessPrior]
    metrics: Tuple[str, ...]
    # Whether the context has a defined extracellular / periplasmic (or in/out) side at all.
    has_membrane_sides: bool = True
    # Whether an LPS shielding band applies to the proximal outer-leaflet surface.
    lps_shielding: bool = False

    def __post_init__(self) -> None:
        if self.membrane_model not in MembraneModel.ALL:
            raise ValueError(f"{self.name}: unknown membrane_model {self.membrane_model!r}")
        if self.orientation_method not in OrientationMethod.ALL:
            raise ValueError(f"{self.name}: unknown orientation_method {self.orientation_method!r}")
        bad = tuple(m for m in self.metrics if m not in Metric.ALL)
        if bad:
            raise ValueError(f"{self.name}: unknown metric(s) {bad}")
        if Metric.ROTATION_INVARIANCE not in self.metrics:
            raise ValueError(f"{self.name}: rotation_invariance is always required")
        if self.membrane_model == MembraneModel.NONE and self.thickness_prior is not None:
            raise ValueError(f"{self.name}: non-membrane context must not carry a thickness prior")
        if self.membrane_model != MembraneModel.NONE and self.thickness_prior is None:
            raise ValueError(f"{self.name}: membrane context requires a thickness prior")
        if self.lps_shielding and self.membrane_model != MembraneModel.ASYMMETRIC_LPS:
            raise ValueError(f"{self.name}: lps_shielding only valid with asymmetric_lps model")

    # -- queries the numeric layers use -------------------------------------------------

    def active_metrics(self) -> Tuple[str, ...]:
        """The metrics that are meaningful in this context (the gating contract)."""
        return self.metrics

    def is_metric_active(self, metric: str) -> bool:
        return metric in self.metrics

    @property
    def has_bilayer(self) -> bool:
        return self.membrane_model != MembraneModel.NONE

    @property
    def is_asymmetric(self) -> bool:
        return self.membrane_model == MembraneModel.ASYMMETRIC_LPS

    def to_dict(self) -> Dict[str, object]:
        """JSON-serializable view (used by the CLI and cross-tool handoff)."""
        tp = None
        if self.thickness_prior is not None:
            tp = {"mean": self.thickness_prior.mean, "sd": self.thickness_prior.sd}
        return {
            "name": self.name,
            "description": self.description,
            "membrane_model": self.membrane_model,
            "orientation_method": self.orientation_method,
            "thickness_prior": tp,
            "metrics": list(self.metrics),
            "has_membrane_sides": self.has_membrane_sides,
            "lps_shielding": self.lps_shielding,
        }


# --------------------------------------------------------------------------------------
# The registry. Five contexts spanning the physics memorient must handle.
# --------------------------------------------------------------------------------------

_CONTEXTS = (
    MembraneContext(
        name="gram_negative_om",
        description=(
            "Gram-negative bacterial outer membrane. Beta-barrel OMPs sit in an asymmetric "
            "bilayer: an LPS outer leaflet that shields the proximal extracellular surface and "
            "a glycerophospholipid inner leaflet. Extracellular loops are antibody-reachable; "
            "periplasmic turns and lipid-facing strands are not."
        ),
        membrane_model=MembraneModel.ASYMMETRIC_LPS,
        orientation_method=OrientationMethod.BARREL_NORMAL,
        thickness_prior=ThicknessPrior(mean=13.0, sd=2.0),  # ~11-16 A half-thickness
        metrics=(
            Metric.AROMATIC_GIRDLE,
            Metric.LIPID_PORE_GAP,
            Metric.HYDROPHOBIC_BELT,
            Metric.LPS_SHIELDING,
            Metric.ROTATION_INVARIANCE,
        ),
        has_membrane_sides=True,
        lps_shielding=True,
    ),
    MembraneContext(
        name="eukaryotic_pm",
        description=(
            "Eukaryotic plasma membrane. Symmetric phospholipid bilayer, no LPS leaflet. "
            "Alpha-helical single- or multi-pass proteins; oriented by hydrophobic belt with the "
            "positive-inside rule breaking the extracellular/cytoplasmic sign."
        ),
        membrane_model=MembraneModel.SYMMETRIC_PHOSPHOLIPID,
        orientation_method=OrientationMethod.TM_HELIX_BELT,
        thickness_prior=ThicknessPrior(mean=15.0, sd=2.0),
        metrics=(
            Metric.AROMATIC_GIRDLE,
            Metric.POSITIVE_INSIDE,
            Metric.HYDROPHOBIC_BELT,
            Metric.LIPID_PORE_GAP,
            Metric.ROTATION_INVARIANCE,
        ),
        has_membrane_sides=True,
        lps_shielding=False,
    ),
    MembraneContext(
        name="tm_receptor",
        description=(
            "Single-pass transmembrane receptor in a symmetric phospholipid membrane. Like "
            "eukaryotic_pm but a single TM helix means no water-filled lumen, so the lipid/pore "
            "gap is not meaningful."
        ),
        membrane_model=MembraneModel.SYMMETRIC_PHOSPHOLIPID,
        orientation_method=OrientationMethod.TM_HELIX_BELT,
        thickness_prior=ThicknessPrior(mean=15.0, sd=2.0),
        metrics=(
            Metric.AROMATIC_GIRDLE,
            Metric.POSITIVE_INSIDE,
            Metric.HYDROPHOBIC_BELT,
            Metric.ROTATION_INVARIANCE,
        ),
        has_membrane_sides=True,
        lps_shielding=False,
    ),
    MembraneContext(
        name="gram_positive_surface",
        description=(
            "Gram-positive cell-wall-anchored surface protein (e.g. LPXTG-sortase, S-layer, LRR "
            "adhesin). No bilayer to fit; the outward axis is defined intrinsically by the "
            "principal axis with the membrane-proximal N-terminal anchor breaking the sign."
        ),
        membrane_model=MembraneModel.NONE,
        orientation_method=OrientationMethod.ANCHOR_RELATIVE,
        thickness_prior=None,
        metrics=(Metric.ROTATION_INVARIANCE,),
        has_membrane_sides=False,
        lps_shielding=False,
    ),
    MembraneContext(
        name="soluble_secreted",
        description=(
            "Secreted / periplasmic soluble antigen (protease, sheath subunit). No membrane at "
            "all: no bilayer fit, no membrane zones, no membrane slab in the viewer. Only "
            "solvent accessibility (SASA) is computed."
        ),
        membrane_model=MembraneModel.NONE,
        orientation_method=OrientationMethod.SASA_ONLY,
        thickness_prior=None,
        metrics=(Metric.ROTATION_INVARIANCE,),
        has_membrane_sides=False,
        lps_shielding=False,
    ),
)

REGISTRY: Dict[str, MembraneContext] = {c.name: c for c in _CONTEXTS}

_DEFAULT = "gram_negative_om"


def get_context(name: str) -> MembraneContext:
    """Look up a context by name, with a helpful error listing valid names."""
    try:
        return REGISTRY[name]
    except KeyError:
        raise KeyError(
            f"unknown membrane context {name!r}; choose one of: {', '.join(sorted(REGISTRY))}"
        ) from None


def list_contexts() -> Tuple[MembraneContext, ...]:
    """All registered contexts, in declaration order."""
    return _CONTEXTS


def default_context() -> MembraneContext:
    """The default context (gram-negative outer membrane — memorient's origin domain)."""
    return REGISTRY[_DEFAULT]
