"""sf_csa — Structure-Function Comparative Species Analysis.

Compares a checksum-pinned protein model against a frozen experimental structure
database and a local proteome universe, keeping structural similarity and
sequence homology as **separate outputs**. A fold match can nominate an
architecture; it can never, by itself, become a function claim.

Interpretation is bounded to a closed six-label vocabulary, and every
campaign-specific judgement — mechanism families, contested groups, divergence
sets, title traps — is a manifest entry rather than a literal in the code.
"""
from .core import (
    CLASSIFICATION_VOCABULARY,
    DEFAULT_CONTESTED_GROUPS,
    DEFAULT_DIVERGENCE_SETS,
    DEFAULT_MECHANISM_FAMILIES,
    DEFAULT_TITLE_TRAPS,
    SFCSError,
    classify_hit,
    classify_title,
    run_pipeline,
    structural_category,
    verify_release,
)

__version__ = "1.1.0"

__all__ = [
    "CLASSIFICATION_VOCABULARY",
    "DEFAULT_CONTESTED_GROUPS",
    "DEFAULT_DIVERGENCE_SETS",
    "DEFAULT_MECHANISM_FAMILIES",
    "DEFAULT_TITLE_TRAPS",
    "SFCSError",
    "__version__",
    "classify_hit",
    "classify_title",
    "run_pipeline",
    "structural_category",
    "verify_release",
]
