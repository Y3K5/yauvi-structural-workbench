"""Join a module's requirements to the registry, the policy, and the cache.

This is the part with the answer in it: for each source a module declares, is it
here, may we go and get it, and if not, what must the human do? `plan` reports
that without touching the network; `get` acts on it.

The plan is also the honest-accounting layer. A source that is required, absent,
and unfetchable makes the plan *not satisfied*, and that fact is surfaced rather
than left for a stage to discover halfway through a run.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .cache import CacheEntry, SourceCache
from .manifest import ModuleManifest, Requirement
from .policy import FetchClass, classify, instructions_for
from .registry import Source, SourceRegistry


@dataclass(frozen=True)
class PlanItem:
    """One declared requirement, resolved against everything we know."""

    requirement: Requirement
    source: Source
    fetch_class: FetchClass
    cached: CacheEntry | None

    @property
    def source_id(self) -> str:
        return self.source.source_id

    @property
    def present(self) -> bool:
        return self.cached is not None

    @property
    def can_fetch(self) -> bool:
        """True when this layer is permitted to retrieve it automatically."""
        return self.fetch_class.may_download

    @property
    def needs_human(self) -> bool:
        """Absent, required in practice, and not something code may go and get."""
        if self.present:
            return False
        return self.fetch_class in (FetchClass.LICENSE_GATED, FetchClass.TABLE_ONLY)

    @property
    def status(self) -> str:
        """One word for the report."""
        if self.fetch_class is FetchClass.INTERNAL:
            return "internal"
        if self.fetch_class is FetchClass.RUNTIME:
            return "runtime"
        if self.present:
            return "present"
        if self.can_fetch:
            return "fetchable"
        return "manual"

    @property
    def blocks(self) -> bool:
        """Would a run be unable to proceed because of this item?"""
        if not self.requirement.required:
            return False
        if self.fetch_class in (FetchClass.INTERNAL, FetchClass.RUNTIME):
            return False
        return not self.present


@dataclass(frozen=True)
class Plan:
    module_id: str
    items: Sequence[PlanItem]

    @property
    def satisfied(self) -> bool:
        """True when nothing required is missing. Fetchable-but-absent still blocks."""
        return not any(item.blocks for item in self.items)

    def blocking(self) -> list[PlanItem]:
        return [item for item in self.items if item.blocks]

    def fetchable(self) -> list[PlanItem]:
        return [item for item in self.items if item.can_fetch and not item.present]

    def manual(self) -> list[PlanItem]:
        return [item for item in self.items if item.needs_human]


def build_plan(
    manifest: ModuleManifest,
    registry: SourceRegistry,
    cache: SourceCache,
) -> Plan:
    """Resolve every requirement. Raises if the module names an undeclared source."""
    # resolve_many reports all unknown ids at once rather than failing on the first.
    registry.resolve_many(manifest.source_ids())

    items = []
    for requirement in manifest.requires:
        source = registry.get(requirement.source_id)
        items.append(
            PlanItem(
                requirement=requirement,
                source=source,
                fetch_class=classify(source),
                cached=cache.latest(source.source_id),
            )
        )
    return Plan(module_id=manifest.module_id, items=tuple(items))


# -- rendering ------------------------------------------------------------

_STATUS_GLYPH = {
    "present": "ok  ",
    "fetchable": "get ",
    "manual": "MANUAL",
    "runtime": "rt  ",
    "internal": "--  ",
}


def render_plan(plan: Plan, *, verbose: bool = False) -> str:
    """Human-readable plan. Never prints anything it has not actually checked."""
    lines: list[str] = []
    lines.append(f"module: {plan.module_id}")
    lines.append(f"sources declared: {len(plan.items)}")
    lines.append("")

    width = max((len(i.source_id) for i in plan.items), default=10)
    for item in plan.items:
        glyph = _STATUS_GLYPH.get(item.status, item.status)
        optional = "" if item.requirement.required else "  (optional)"
        lines.append(f"  [{glyph:<6}] {item.source_id:<{width}}  {item.source.display_name}{optional}")
        if item.requirement.role:
            lines.append(f"           {' ' * width}  role: {item.requirement.role}")
        if item.cached:
            lines.append(
                f"           {' ' * width}  cached {item.cached.sha256[:12]} "
                f"({item.cached.bytes} bytes, {item.cached.retrieved_at})"
            )
            if item.cached.version:
                lines.append(f"           {' ' * width}  version: {item.cached.version}")
        if verbose and item.source.license_note:
            lines.append(f"           {' ' * width}  licence: {item.source.license_note}")

    manual = plan.manual()
    if manual:
        lines.append("")
        lines.append("Sources this tool will not retrieve for you:")
        for item in manual:
            lines.append("")
            lines.append(f"  {item.source_id} — {item.source.display_name}")
            for line in instructions_for(item.source).splitlines():
                lines.append(f"    {line}")

    lines.append("")
    blocking = plan.blocking()
    if blocking:
        lines.append(
            f"NOT SATISFIED — {len(blocking)} required source(s) absent: "
            + ", ".join(i.source_id for i in blocking)
        )
    else:
        lines.append("SATISFIED — every required source is present.")
    return "\n".join(lines)
