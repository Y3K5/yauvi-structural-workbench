"""yauvi-sources — the executable half of the evidence-source registry.

`catalogs/sources.yaml` declares every external database, reference panel, and
predictor the platform can compare a protein against. This package reads that
declaration and acts on it: it reports what a module needs, retrieves what the
licence permits, refuses what it does not, and records the digest of everything
that lands on disk.

It never invents a substitute for a source it could not obtain.
"""
from .cache import CacheEntry, CacheError, SourceCache, default_cache_dir, sha256_file
from .manifest import ManifestError, ModuleManifest, Requirement, load_manifest, resolve_manifest
from .planner import Plan, PlanItem, build_plan, render_plan
from .policy import FetchClass, PolicyError, classify, instructions_for
from .registry import RegistryError, Source, SourceRegistry
from .refresh import ReferenceRefreshManager

__version__ = "0.1.0"

__all__ = [
    "CacheEntry",
    "CacheError",
    "FetchClass",
    "ManifestError",
    "ModuleManifest",
    "Plan",
    "PlanItem",
    "PolicyError",
    "RegistryError",
    "Requirement",
    "ReferenceRefreshManager",
    "Source",
    "SourceCache",
    "SourceRegistry",
    "__version__",
    "build_plan",
    "classify",
    "default_cache_dir",
    "instructions_for",
    "load_manifest",
    "render_plan",
    "resolve_manifest",
    "sha256_file",
]
