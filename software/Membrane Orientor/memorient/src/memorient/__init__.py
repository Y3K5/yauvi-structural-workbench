"""memorient — context-aware membrane orientor, structural-biology labeler, and 3D viz.

Public API is deliberately small. The membrane-context registry is always importable
(stdlib-only); the numeric layers (geometry/sasa/barrel/membrane/orientor/labeler/viz)
import numpy + biopython lazily, so ``import memorient`` never pulls heavy deps until you
call into them.
"""

from __future__ import annotations

from .contexts import (
    MembraneContext,
    MembraneModel,
    Metric,
    OrientationMethod,
    ThicknessPrior,
    REGISTRY,
    default_context,
    get_context,
    list_contexts,
)

__version__ = "0.3.0"

__all__ = [
    "MembraneContext",
    "MembraneModel",
    "Metric",
    "OrientationMethod",
    "ThicknessPrior",
    "REGISTRY",
    "default_context",
    "get_context",
    "list_contexts",
    "__version__",
]
