"""Load `catalogs/sources.yaml` as typed records.

The registry file is the identity/provenance/meaning layer for every external
database, reference panel, and predictor the platform can compare a protein
against. It was written to be read by humans; this module makes it executable
without changing its meaning.

Two rules are enforced here rather than left to callers:

* **A declaration is not a download.** Loading a source tells you what it is,
  how it is reached, and what its licence permits. It never implies the file is
  present. `Source.access` is the only thing that decides whether code may
  fetch it, and that decision lives in `policy.py`.
* **An unknown source is an error, not a default.** Asking for a `source_id`
  that the registry does not declare raises, so a typo cannot silently become
  an unfetched-but-unreported input.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import yaml


class RegistryError(RuntimeError):
    """The registry file is malformed, or a requested source is not declared."""


# Fields every source entry must carry. The registry's own header documents
# these; we refuse to load a file where one is missing rather than filling in a
# default, because every one of them is provenance.
REQUIRED_FIELDS = (
    "source_id",
    "display_name",
    "kind",
    "channel",
    "status",
    "access",
)

# The `status` vocabulary, copied from the registry header. Anything else means
# the file and this loader have drifted apart.
STATUS_VOCABULARY = frozenset(
    {
        "wired",
        "configured_optional",
        "heuristic_stub",
        "table_only",
        "dead_config",
        "proposed",
    }
)


@dataclass(frozen=True)
class Citation:
    pmid: str = ""
    ref: str = ""


@dataclass(frozen=True)
class Source:
    """One declared evidence source, exactly as the registry describes it."""

    source_id: str
    display_name: str
    kind: str
    channel: str
    status: str
    access: str
    what_it_is: str = ""
    what_we_check: str = ""
    hit_means: str = ""
    stage: str = ""
    operation: str = "NONE"
    config_key: str = ""
    code_ref: str = ""
    versioning: str = ""
    license_note: str = ""
    thresholds: Mapping[str, Any] = field(default_factory=dict)
    limits: Sequence[str] = ()
    citations: Sequence[Citation] = ()
    # Optional retrieval hints. Absent for most entries today; the fetchers fall
    # back to per-kind defaults when they are missing.
    url: str = ""
    homepage: str = ""
    documentation_url: str = ""
    artifact_types: Sequence[str] = ()
    identifier_pattern: str = ""
    format_guides: Mapping[str, str] = field(default_factory=dict)
    manual_instructions: str = ""

    @property
    def is_internal(self) -> bool:
        """True when the 'source' is an in-code heuristic with nothing to fetch."""
        return self.access == "internal"


def _coerce_citations(raw: Any) -> tuple[Citation, ...]:
    if not raw:
        return ()
    out: list[Citation] = []
    for item in raw:
        if isinstance(item, Mapping):
            out.append(Citation(pmid=str(item.get("pmid", "")), ref=str(item.get("ref", ""))))
        else:
            out.append(Citation(ref=str(item)))
    return tuple(out)


def _build_source(raw: Mapping[str, Any], *, index: int) -> Source:
    missing = [key for key in REQUIRED_FIELDS if not raw.get(key)]
    if missing:
        label = raw.get("source_id") or f"entry #{index}"
        raise RegistryError(f"source {label} is missing required field(s): {', '.join(missing)}")

    status = str(raw["status"])
    if status not in STATUS_VOCABULARY:
        raise RegistryError(
            f"source {raw['source_id']} declares unknown status {status!r}; "
            f"expected one of {sorted(STATUS_VOCABULARY)}"
        )

    return Source(
        source_id=str(raw["source_id"]),
        display_name=str(raw["display_name"]),
        kind=str(raw["kind"]),
        channel=str(raw["channel"]),
        status=status,
        access=str(raw["access"]),
        what_it_is=str(raw.get("what_it_is", "")).strip(),
        what_we_check=str(raw.get("what_we_check", "")).strip(),
        hit_means=str(raw.get("hit_means", "")).strip(),
        stage=str(raw.get("stage", "")),
        operation=str(raw.get("operation", "NONE")),
        config_key=str(raw.get("config_key", "")),
        code_ref=str(raw.get("code_ref", "")),
        versioning=str(raw.get("versioning", "")).strip(),
        license_note=str(raw.get("license_note", "")).strip(),
        thresholds=dict(raw.get("thresholds") or {}),
        limits=tuple(str(item) for item in (raw.get("limits") or ())),
        citations=_coerce_citations(raw.get("citations")),
        url=str(raw.get("url", "")),
        homepage=str(raw.get("homepage", "")),
        documentation_url=str(raw.get("documentation_url", "")),
        artifact_types=tuple(str(item) for item in (raw.get("artifact_types") or ())),
        identifier_pattern=str(raw.get("identifier_pattern", "")),
        format_guides={str(key): str(value) for key, value in dict(raw.get("format_guides") or {}).items()},
        manual_instructions=str(raw.get("manual_instructions", "")).strip(),
    )


class SourceRegistry:
    """The declared evidence sources, keyed by `source_id`."""

    def __init__(self, sources: Iterable[Source], *, catalog_id: str = "", updated_at: str = ""):
        self._sources: dict[str, Source] = {}
        for source in sources:
            if source.source_id in self._sources:
                raise RegistryError(f"duplicate source_id in registry: {source.source_id}")
            self._sources[source.source_id] = source
        self.catalog_id = catalog_id
        self.updated_at = updated_at

    # -- construction -----------------------------------------------------

    @classmethod
    def load(cls, path: str | Path) -> "SourceRegistry":
        path = Path(path)
        if not path.is_file():
            raise RegistryError(f"source registry not found: {path}")
        try:
            document = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            raise RegistryError(f"source registry is not valid YAML: {exc}") from exc
        if not isinstance(document, Mapping):
            raise RegistryError(f"source registry must be a mapping at the top level: {path}")

        raw_sources = document.get("sources")
        if not isinstance(raw_sources, list) or not raw_sources:
            raise RegistryError(f"source registry declares no sources: {path}")

        sources = [_build_source(item, index=i) for i, item in enumerate(raw_sources)]
        return cls(
            sources,
            catalog_id=str(document.get("catalog_id", "")),
            updated_at=str(document.get("updated_at", "")),
        )

    # -- lookup -----------------------------------------------------------

    def __contains__(self, source_id: object) -> bool:
        return source_id in self._sources

    def __len__(self) -> int:
        return len(self._sources)

    def __iter__(self):
        return iter(self._sources.values())

    def get(self, source_id: str) -> Source:
        """Return one declared source, or raise. There is no default."""
        try:
            return self._sources[source_id]
        except KeyError:
            raise RegistryError(
                f"source {source_id!r} is not declared in the registry. "
                f"Add it to catalogs/sources.yaml before any code depends on it."
            ) from None

    def resolve_many(self, source_ids: Iterable[str]) -> list[Source]:
        """Resolve a module's declared requirements, reporting every miss at once."""
        wanted = list(source_ids)
        unknown = [sid for sid in wanted if sid not in self._sources]
        if unknown:
            raise RegistryError(
                "module requires source(s) not declared in the registry: "
                + ", ".join(sorted(unknown))
            )
        return [self._sources[sid] for sid in wanted]

    def ids(self) -> list[str]:
        return sorted(self._sources)
