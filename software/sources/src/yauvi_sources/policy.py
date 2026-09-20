"""What code is allowed to do with a declared source.

The registry's `access` field already records how each source is reached. This
module turns that description into a decision, so the rule lives in one place
instead of being re-derived at every call site.

The governing principle is the one `shared/runtime-registry.yaml` states for
runtimes and this layer inherits for data: **fail closed**. A source that cannot
be fetched is reported as unfetched. It is never approximated, never substituted,
and never silently skipped.

Four classes:

``OPEN_FETCHABLE``  a public endpoint we may call; download and hash it.
``LICENSE_GATED``   redistribution or automated retrieval is not ours to perform.
                    Print what the human must do; download nothing; exit non-zero.
``TABLE_ONLY``      the platform consumes someone else's export. We validate the
                    shape of a user-supplied file; we never produce it.
``RUNTIME``         not a file at all — an executable, resolved by the runtime
                    registry's preflight, not by this layer.
``INTERNAL``        an in-code heuristic. Nothing exists to fetch.
"""
from __future__ import annotations

from enum import Enum
from typing import Mapping

from .registry import Source


class FetchClass(str, Enum):
    OPEN_FETCHABLE = "open_fetchable"
    LICENSE_GATED = "license_gated"
    TABLE_ONLY = "table_only"
    RUNTIME = "runtime"
    INTERNAL = "internal"

    @property
    def may_download(self) -> bool:
        return self is FetchClass.OPEN_FETCHABLE


# Every `access` value that appears in a registry file, mapped to what code
# may do with it. Keeping this exhaustive is deliberate: an unmapped access value
# raises rather than defaulting, so adding a source with a new access mode forces
# an explicit decision about whether it may be downloaded.
ACCESS_TO_CLASS: Mapping[str, FetchClass] = {
    # public endpoints
    "network_rest_api": FetchClass.OPEN_FETCHABLE,
    "network_api": FetchClass.OPEN_FETCHABLE,
    "network_download": FetchClass.OPEN_FETCHABLE,
    "network_api_and_download": FetchClass.OPEN_FETCHABLE,
    "ftp_download": FetchClass.OPEN_FETCHABLE,
    # a human must accept a licence or complete a registration first
    "manual_download": FetchClass.LICENSE_GATED,
    "manual_download_licensed": FetchClass.LICENSE_GATED,
    # someone else's tool produced a table; we only read it
    "manual_export": FetchClass.TABLE_ONLY,
    "web_server_manual": FetchClass.TABLE_ONLY,
    # executables, not files
    "local_binary": FetchClass.RUNTIME,
    "local_binary_or_docker": FetchClass.RUNTIME,
    "local_binary_or_biolib": FetchClass.RUNTIME,
    "local_binary_academic": FetchClass.RUNTIME,
    "local_binary_or_api": FetchClass.RUNTIME,
    "local_docker": FetchClass.RUNTIME,
    # computed in-process
    "internal": FetchClass.INTERNAL,
}


class PolicyError(RuntimeError):
    """An access mode that no policy decision has been made about."""


def classify(source: Source) -> FetchClass:
    """Decide what may be done with a source.

    `status` is consulted before `access`, because the two answer different
    questions and the registry is explicit about which one governs. `access`
    describes the transport a source *has*; `status` records the relationship the
    platform has *chosen* to it. The registry defines `table_only` as:

        "core accepts its export as a table; the platform never fetches or runs
         it — the human does, out of band"

    IEDB is the case that makes this matter: it is reachable over an API
    (`access: network_api_and_download`), so transport alone would mark it
    downloadable, but its status is `table_only` and the registry states plainly
    that "no code path fetches it". Honouring access over status would have this
    layer start retrieving a source the platform deliberately does not automate.
    """
    if source.status == "table_only":
        return FetchClass.TABLE_ONLY
    try:
        return ACCESS_TO_CLASS[source.access]
    except KeyError:
        raise PolicyError(
            f"source {source.source_id!r} declares access mode {source.access!r}, "
            f"which has no fetch policy. Add it to ACCESS_TO_CLASS with an explicit "
            f"decision — an unmapped access mode must never default to downloadable."
        ) from None


# Acquisition instructions and table shapes used to live here, keyed by
# source_id. They no longer do. A description of how to obtain a source belongs
# to the entry that declares the source, not to the code that reads it: keeping
# it here meant this package carried the source list of every catalogue that had
# ever been pointed at it, including entries no public registry declares and no
# public code path can reach. `Source.manual_instructions` and
# `Source.table_expectation` carry that text now, so a catalogue describes its
# own sources and this module only decides what kind of instruction to print.
#
# `manual_instructions` was already parsed into `Source` and then ignored, which
# is how the two drifted apart in the first place.


def instructions_for(source: Source) -> str:
    """Human-facing acquisition text for a source code may not fetch."""
    fetch_class = classify(source)
    if fetch_class is FetchClass.LICENSE_GATED:
        base = source.manual_instructions or (
            f"{source.display_name} must be obtained manually; the registry records\n"
            f"its access mode as {source.access!r}. Stage it with:\n"
            f"  yauvi-fetch stage {source.source_id} <path>"
        )
        if source.license_note:
            base += f"\n  Licence: {source.license_note}"
        return base
    if fetch_class is FetchClass.TABLE_ONLY:
        expectation = (
            source.table_expectation or "a tabular export keyed by protein identifier"
        )
        return (
            f"{source.display_name} is consumed as a table the platform never produces.\n"
            f"  Run it yourself, then stage the export:\n"
            f"    yauvi-fetch stage {source.source_id} <path>\n"
            f"  Expected shape: {expectation}"
        )
    if fetch_class is FetchClass.RUNTIME:
        return (
            f"{source.display_name} is an executable, not a data file.\n"
            f"  Install it and put it on PATH; preflight resolves it through\n"
            f"  shared/runtime-registry.yaml, which is fail-closed on absence."
        )
    if fetch_class is FetchClass.INTERNAL:
        return f"{source.display_name} is computed in-process. Nothing to acquire."
    return f"{source.display_name} is fetched automatically by `yauvi-fetch get`."
