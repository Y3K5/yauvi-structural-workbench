"""actstate — is this protein in a functionally competent state?

Five independent signals — active-site completeness, geometry, cofactor
occupancy, conformation against references of known state, and assembly
dependence — are computed and reported separately, then reduced to one label
from a closed six-value vocabulary.

The module is deliberately conservative about what it will conclude. An entry
with no annotated active site is `indeterminate`, not `inactive`. A predicted
model alone can never yield `active_state_supported`. A signal that could not be
evaluated is recorded as unavailable rather than quietly dropped.
"""
from .core import (
    ActivityAssessment,
    LABELS,
    ProteinRecord,
    SIGNAL_STATES,
    Signal,
    assess,
    assign_label,
)
from .features import Feature, FeatureSet, parse_cofactors, parse_features
from .structure import Residue, Structure, StructureError, read_structure

__version__ = "0.1.0"

__all__ = [
    "ActivityAssessment",
    "Feature",
    "FeatureSet",
    "LABELS",
    "ProteinRecord",
    "Residue",
    "SIGNAL_STATES",
    "Signal",
    "Structure",
    "StructureError",
    "__version__",
    "assess",
    "assign_label",
    "parse_cofactors",
    "parse_features",
    "read_structure",
]
