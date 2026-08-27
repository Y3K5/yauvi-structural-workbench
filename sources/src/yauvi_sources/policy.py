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


# Every `access` value that appears in catalogs/sources.yaml, mapped to what code
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


# Acquisition instructions for the license-gated sources, so `plan` can tell a
# user exactly what to do instead of just refusing. Keyed by source_id; the
# registry's own `license_note` and `versioning` are printed alongside these.
MANUAL_INSTRUCTIONS: Mapping[str, str] = {
    "deg": (
        "DEG (Database of Essential Genes) requires accepting its academic terms.\n"
        "  1. Register at http://origin.tubic.org/deg/public/index.php\n"
        "  2. Download the bacterial protein set (DEG10 / deg_bacteria.faa)\n"
        "  3. Stage it and record its digest:\n"
        "       yauvi-fetch stage deg <path-to-deg_bacteria.faa>"
    ),
    "drugbank": (
        "DrugBank requires a licence for the full data download; the academic\n"
        "licence is free but must be requested per user.\n"
        "  1. Request access at https://go.drugbank.com/releases/latest\n"
        "  2. Download the protein identifiers / target sequences FASTA\n"
        "  3. yauvi-fetch stage drugbank <path>"
    ),
    "ogee": (
        "OGEE bulk downloads are served from https://v3.ogee.info/#/downloads\n"
        "after accepting the site terms.\n"
        "  yauvi-fetch stage ogee <path>"
    ),
    "depmap": (
        "DepMap releases are versioned and require accepting the DepMap terms.\n"
        "  1. Choose a release at https://depmap.org/portal/data_page/\n"
        "  2. Download the gene-effect matrix\n"
        "  3. yauvi-fetch stage depmap <path>   (record the release name)"
    ),
    "human_protein_atlas": (
        "The Human Protein Atlas bulk TSV is at https://www.proteinatlas.org/about/download\n"
        "and is CC BY-SA 3.0 — attribution is required wherever its values are shown.\n"
        "  yauvi-fetch stage human_protein_atlas <path>"
    ),
}


# Sources whose files we never hold at all: the human runs a web server or a
# licensed binary and gives us the export. `plan` prints the expected shape.
TABLE_EXPECTATIONS: Mapping[str, str] = {
    "vaxijen": "TSV/CSV export with columns: protein_id, score, prediction (threshold declared in the registry).",
    "allertop": "Export with columns: protein_id, prediction (allergen | non-allergen).",
    "toxinpred": "Export with columns: peptide_or_protein_id, score, prediction.",
    "deeploc2": "DeepLoc 2 CSV: protein_id, localization, per-compartment probabilities.",
    "cello": "CELLO output table: protein_id, localization, reliability.",
    "netmhcpan": "NetMHCpan output: peptide, allele, rank, affinity.",
    "iedb": "IEDB export (epitope table): epitope, source antigen, assay, MHC restriction.",
    "chembl": "ChEMBL export: target_chembl_id, uniprot accession, activity summary.",
}


def instructions_for(source: Source) -> str:
    """Human-facing acquisition text for a source code may not fetch."""
    fetch_class = classify(source)
    if fetch_class is FetchClass.LICENSE_GATED:
        base = MANUAL_INSTRUCTIONS.get(
            source.source_id,
            f"{source.display_name} must be obtained manually; the registry records\n"
            f"its access mode as {source.access!r}. Stage it with:\n"
            f"  yauvi-fetch stage {source.source_id} <path>",
        )
        if source.license_note:
            base += f"\n  Licence: {source.license_note}"
        return base
    if fetch_class is FetchClass.TABLE_ONLY:
        expectation = TABLE_EXPECTATIONS.get(
            source.source_id, "a tabular export keyed by protein identifier"
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
