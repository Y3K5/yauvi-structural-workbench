"""Reference-bounded conformational state and MD ensemble analysis."""

from .core import InputError, analyze, validate_reference_set, write_outputs

__all__ = ["InputError", "analyze", "validate_reference_set", "write_outputs"]
__version__ = "0.2.0"
