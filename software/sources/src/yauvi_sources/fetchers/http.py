"""Network retrieval for the sources whose access mode permits it.

Every function here follows the contract set by the bridge's AlphaFold fetcher,
which is the one retrieval path in the tree that was already written correctly:
it resolves the real download URL from the provider's API rather than guessing a
version-shaped path, it caches, and it returns a result carrying a
machine-readable reason instead of raising on a network miss. The orchestrator
decides what a miss means; the fetcher only reports it.

`requests` is imported lazily so that `yauvi-fetch plan` — which never touches
the network — works in an install without the `fetch` extra.
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
import re
import socket
import threading
import time
from typing import Any, Callable, Mapping
import urllib.request
from urllib.parse import urlparse

# (connect, read). Split deliberately: a host that is unreachable should be
# reported in seconds, while a legitimate multi-hundred-megabyte proteome stream
# needs a long read budget. A single scalar cannot serve both, and using one
# makes an offline machine wait the full download timeout on every source before
# reporting what it could already have known immediately.
DEFAULT_TIMEOUT: tuple[int, int] = (10, 120)

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_PROTEOMES = "https://rest.uniprot.org/proteomes"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
RCSB_FILE = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_ASSEMBLY = "https://files.rcsb.org/download/{pdb_id}-assembly{assembly_id}.cif"
WWPDB_VALIDATION = "https://files.rcsb.org/pub/pdb/validation_reports/{middle}/{pdb_id}/{pdb_id}_validation.xml.gz"
UNIPROT_ENTRY_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
HPA_BULK = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
RHEA_BULK = "https://ftp.expasy.org/databases/rhea/tsv/rhea-tsv.tar.gz"
PROTEOMEXCHANGE_DATASET = "https://proteomecentral.proteomexchange.org/cgi/GetDataset"
PUBLIC_MAX_BYTES = 128 * 1024 * 1024
PUBLIC_ALLOWED_HOSTS = {
    "www.proteinatlas.org", "ftp.expasy.org", "proteomecentral.proteomexchange.org",
    "v31a.homd.org", "homd.org", "www.homd.org",
    "files.rcsb.org", "rest.uniprot.org", "alphafold.ebi.ac.uk",
}
_PROVIDER_LOCKS = {host: threading.Lock() for host in PUBLIC_ALLOWED_HOSTS}

# The annotation columns the platform reads. The first twelve are what the
# existing per-taxon annotation TSVs already carry; the last five are added for
# the activity-state module, which needs catalytic and cofactor features that no
# current export includes.
UNIPROT_ANNOTATION_FIELDS = (
    "accession",
    "protein_name",
    "gene_names",
    "protein_existence",
    "cc_function",
    "cc_subcellular_location",
    "go_id",
    "keyword",
    "xref_pfam",
    "xref_interpro",
    "lit_pubmed_id",
    "ec",
    # --- activity-state inputs ---
    "ft_act_site",
    "ft_binding",
    "ft_site",
    "cc_cofactor",
    "cc_activity_regulation",
)


class FetchError(RuntimeError):
    """Raised only for programming errors, never for a network miss."""


@dataclass(frozen=True)
class FetchOutcome:
    """The result of one retrieval attempt. `ok=False` never means an exception."""

    ok: bool
    payload: bytes = b""
    filename: str = ""
    origin: str = ""
    version: str = ""
    reason: str = ""


def _requests():
    try:
        import requests  # noqa: PLC0415 — lazy so `plan` works without the extra
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise FetchError(
            "network retrieval needs the 'fetch' extra: pip install 'yauvi-sources[fetch]'"
        ) from exc
    return requests


def _get(url: str, *, params: Mapping[str, Any] | None = None, timeout=DEFAULT_TIMEOUT):
    """GET a URL, converting transport failures into a reason string."""
    requests = _requests()
    try:
        response = requests.get(url, params=params, timeout=timeout)
    except requests.exceptions.ConnectTimeout:
        return None, "network_unreachable:connect_timeout"
    except requests.exceptions.ConnectionError:
        return None, "network_unreachable:connection_error"
    except requests.exceptions.RequestException as exc:
        return None, f"network_error:{type(exc).__name__}"
    if response.status_code == 404:
        return None, "not_found"
    if response.status_code != 200:
        return None, f"http_{response.status_code}"
    return response, ""


def _uniprot_release(response: Any) -> str:
    """The UniProt release the response was served from.

    The registry records this gap honestly: "the release is whatever UniProt
    serves that day." UniProt does report the release in a response header, so
    the fix is to capture it at retrieval time rather than leaving the run
    pinned only by date.
    """
    for header in ("X-UniProt-Release", "x-uniprot-release"):
        value = response.headers.get(header)
        if value:
            date = response.headers.get("X-UniProt-Release-Date", "")
            return f"{value}" + (f" ({date})" if date else "")
    return ""


# -- reachability ---------------------------------------------------------

# Name resolution is not covered by a socket timeout: `getaddrinfo` is a
# blocking call into the system resolver, so on a machine with no route out a
# request can sit far longer than any `timeout=` argument suggests. Probing once,
# with a bound enforced from outside the call, means an offline run reports that
# it is offline in seconds instead of stalling on every source in turn.

PROBE_TIMEOUT = 3.0


def proxy_configured() -> bool:
    """Is an HTTP(S) proxy in effect for outbound requests?"""
    proxies = urllib.request.getproxies()
    return bool(proxies.get("http") or proxies.get("https"))


def host_reachable(host: str, *, timeout: float = PROBE_TIMEOUT) -> bool:
    """Can this host be resolved within `timeout` seconds?

    **Returns True unconditionally when a proxy is configured.** Behind a proxy
    the client never resolves the target host — the proxy does — so a direct
    `getaddrinfo` says nothing about reachability and fails even when requests
    succeed. Trusting it there would refuse to fetch on exactly the networks
    where fetching works: CI runners, and anything corporate. The retrieval
    itself still reports its own failure, so skipping the probe costs a slower
    error message and nothing else.

    Otherwise a daemon thread bounds the lookup, because that is the only
    mechanism that bounds `getaddrinfo` on every platform. If the probe has not
    finished when the budget expires the thread is abandoned; it holds no
    resource the caller needs.
    """
    if proxy_configured():
        return True

    import threading  # noqa: PLC0415

    result: list[bool] = []

    def resolve() -> None:
        try:
            socket.getaddrinfo(host, None)
            result.append(True)
        except OSError:
            result.append(False)

    worker = threading.Thread(target=resolve, daemon=True)
    worker.start()
    worker.join(timeout)
    return bool(result) and result[0]


def endpoint_reachable(url: str, *, timeout: float = PROBE_TIMEOUT) -> bool:
    """Reachability of the host named in a URL."""
    host = urlparse(url).hostname
    return host_reachable(host, timeout=timeout) if host else False


# The hosts this module retrieves from, for a single up-front probe.
FETCH_HOSTS = (
    "rest.uniprot.org",
    "alphafold.ebi.ac.uk",
    "files.rcsb.org",
    "www.proteinatlas.org",
    "ftp.expasy.org",
    "proteomecentral.proteomexchange.org",
    "v31a.homd.org",
)


# -- individual sources ---------------------------------------------------


def fetch_uniprot_proteome(proteome_id: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Stream one reference proteome as FASTA, recording the UniProt release."""
    params = {
        "query": f"proteome:{proteome_id}",
        "format": "fasta",
        "includeIsoform": "false",
    }
    response, reason = _public_get(UNIPROT_STREAM, params=params, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    payload = response.content
    if not payload.startswith(b">"):
        return FetchOutcome(ok=False, reason="not_fasta")
    return FetchOutcome(
        ok=True,
        payload=payload,
        filename=f"{proteome_id}.fasta",
        origin=response.url,
        version=_uniprot_release(response),
    )


def fetch_uniprot_annotation(
    query: str, *, filename: str, timeout=DEFAULT_TIMEOUT
) -> FetchOutcome:
    """Stream the annotation TSV, including the catalytic/cofactor features.

    `query` is a UniProt query string, e.g. `taxonomy_id:837` or
    `proteome:UP000005640`.
    """
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(UNIPROT_ANNOTATION_FIELDS),
    }
    response, reason = _public_get(UNIPROT_STREAM, params=params, timeout=timeout, max_bytes=32 * 1024 * 1024)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    payload = response.content
    if not payload or b"\t" not in payload.split(b"\n", 1)[0]:
        return FetchOutcome(ok=False, reason="not_tsv")
    return FetchOutcome(
        ok=True,
        payload=payload,
        filename=filename,
        origin=response.url,
        version=_uniprot_release(response),
    )


def fetch_alphafold_model(accession: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Fetch the current AlphaFold model for a UniProt accession.

    The static `AF-<acc>-F1-model_v4.pdb` path is not stable — model versions
    advance — so the download URL is resolved from the prediction API, exactly as
    `bridge/fetch.py` does.
    """
    record, reason = _alphafold_record(accession, timeout=timeout)
    if record is None:
        return FetchOutcome(ok=False, reason=reason)
    url = record.get("pdbUrl")
    if not url:
        return FetchOutcome(ok=False, reason="no_model")

    model, model_reason = _public_get(str(url), timeout=timeout)
    if model is None:
        return FetchOutcome(ok=False, reason=model_reason)
    payload = model.content
    if not payload.startswith((b"HEADER", b"ATOM", b"CRYST", b"REMARK", b"TITLE", b"data_")):
        return FetchOutcome(ok=False, reason="not_a_structure")
    return FetchOutcome(
        ok=True,
        payload=payload,
        filename=f"AF-{accession}-F1.pdb",
        origin=url,
        version=str(record.get("latestVersion", "")),
    )


def fetch_pdb_structure(pdb_id: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Fetch one experimental structure from RCSB as mmCIF."""
    identifier = pdb_id.lower().strip()
    url = RCSB_FILE.format(pdb_id=identifier)
    response, reason = _public_get(url, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    payload = response.content
    if not payload.startswith(b"data_"):
        return FetchOutcome(ok=False, reason="not_mmcif")
    return FetchOutcome(
        ok=True, payload=payload, filename=f"{identifier}.cif", origin=url
    )


def _pdb_identifier(value: str) -> str | None:
    identifier = value.strip().lower()
    return identifier if re.fullmatch(r"[a-z0-9]{4}", identifier) else None


def _uniprot_identifier(value: str) -> str | None:
    identifier = value.strip().upper()
    return identifier if re.fullmatch(r"[A-Z0-9][A-Z0-9_-]{5,19}", identifier) else None


def fetch_pdb_assembly(value: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Fetch an explicitly selected deposited biological assembly as mmCIF."""
    pdb_value, separator, assembly_value = value.partition(":")
    identifier = _pdb_identifier(pdb_value)
    assembly_id = assembly_value if separator else "1"
    if identifier is None or not re.fullmatch(r"[1-9][0-9]{0,2}", assembly_id):
        return FetchOutcome(ok=False, reason="invalid_pdb_assembly_identifier")
    url = RCSB_ASSEMBLY.format(pdb_id=identifier, assembly_id=assembly_id)
    response, reason = _public_get(url, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    payload = response.content
    if not payload.lstrip().startswith(b"data_"):
        return FetchOutcome(ok=False, reason="not_mmcif")
    return FetchOutcome(True, payload, f"{identifier}-assembly{assembly_id}.cif", response.url,
                        str(response.headers.get("Last-Modified", "")))


def fetch_wwpdb_validation(value: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Fetch and safely decompress the wwPDB validation XML for one PDB ID."""
    identifier = _pdb_identifier(value)
    if identifier is None:
        return FetchOutcome(ok=False, reason="invalid_pdb_id")
    url = WWPDB_VALIDATION.format(middle=identifier[1:3], pdb_id=identifier)
    response, reason = _public_get(url, timeout=timeout, max_bytes=64 * 1024 * 1024)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    try:
        payload = gzip.decompress(response.content)
    except (OSError, EOFError):
        return FetchOutcome(ok=False, reason="not_gzip_xml")
    if len(payload) > PUBLIC_MAX_BYTES or not payload.lstrip().startswith(b"<"):
        return FetchOutcome(ok=False, reason="not_validation_xml")
    return FetchOutcome(True, payload, f"{identifier}_validation.xml", response.url,
                        str(response.headers.get("Last-Modified", "")))


def fetch_uniprot_sequence(value: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    identifier = _uniprot_identifier(value)
    if identifier is None:
        return FetchOutcome(ok=False, reason="invalid_uniprot_accession")
    url = UNIPROT_ENTRY_FASTA.format(accession=identifier)
    response, reason = _public_get(url, timeout=timeout, max_bytes=8 * 1024 * 1024)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    if not response.content.lstrip().startswith(b">"):
        return FetchOutcome(ok=False, reason="not_fasta")
    return FetchOutcome(True, response.content, f"{identifier}.fasta", response.url, _uniprot_release(response))


def _alphafold_record(accession: str, *, timeout=DEFAULT_TIMEOUT) -> tuple[Mapping[str, Any] | None, str]:
    identifier = _uniprot_identifier(accession)
    if identifier is None:
        return None, "invalid_uniprot_accession"
    response, reason = _public_get(AFDB_API.format(accession=identifier), timeout=timeout, max_bytes=2 * 1024 * 1024)
    if response is None:
        return None, reason
    try:
        entries = response.json()
    except ValueError:
        return None, "bad_api_json"
    if not entries or not isinstance(entries, list) or not isinstance(entries[0], Mapping):
        return None, "no_model"
    return entries[0], ""


def fetch_alphafold_pae(value: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    record, reason = _alphafold_record(value, timeout=timeout)
    if record is None:
        return FetchOutcome(ok=False, reason=reason)
    url = str(record.get("paeDocUrl") or "")
    if not url:
        return FetchOutcome(ok=False, reason="no_pae")
    response, reason = _public_get(url, timeout=timeout, max_bytes=64 * 1024 * 1024)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    try:
        value_json = response.json()
    except ValueError:
        return FetchOutcome(ok=False, reason="not_pae_json")
    if not isinstance(value_json, (list, dict)):
        return FetchOutcome(ok=False, reason="not_pae_json")
    identifier = str(value).strip().upper()
    return FetchOutcome(True, response.content, f"AF-{identifier}-F1-pae.json", response.url,
                        str(record.get("latestVersion", "")))


def fetch_structural_artifact(artifact_type: str, identifier: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Acquire one allowlisted structural artifact selected by type and accession."""
    if artifact_type == "pdb.coordinates":
        pdb_id = _pdb_identifier(identifier)
        return fetch_pdb_structure(pdb_id or "", timeout=timeout) if pdb_id else FetchOutcome(ok=False, reason="invalid_pdb_id")
    if artifact_type == "pdb.biological_assembly":
        return fetch_pdb_assembly(identifier, timeout=timeout)
    if artifact_type == "wwpdb.validation":
        return fetch_wwpdb_validation(identifier, timeout=timeout)
    if artifact_type == "alphafold.model":
        accession = _uniprot_identifier(identifier)
        return fetch_alphafold_model(accession or "", timeout=timeout) if accession else FetchOutcome(ok=False, reason="invalid_uniprot_accession")
    if artifact_type == "alphafold.pae":
        return fetch_alphafold_pae(identifier, timeout=timeout)
    if artifact_type == "uniprot.sequence":
        return fetch_uniprot_sequence(identifier, timeout=timeout)
    if artifact_type == "uniprot.annotations":
        accession = _uniprot_identifier(identifier)
        if accession is None:
            return FetchOutcome(ok=False, reason="invalid_uniprot_accession")
        return fetch_uniprot_annotation(f"accession:{accession}", filename=f"{accession}-annotations.tsv", timeout=timeout)
    if artifact_type == "uniprot.proteome":
        value = identifier.strip().upper()
        if not re.fullmatch(r"UP[0-9]{9}", value):
            return FetchOutcome(ok=False, reason="invalid_uniprot_proteome_id")
        return fetch_uniprot_proteome(value, timeout=timeout)
    return FetchOutcome(ok=False, reason="unsupported_structural_artifact")


def fetch_url(url: str, *, filename: str = "", timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Retrieve a source whose registry entry names a direct download URL."""
    response, reason = _get(url, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    if not response.content:
        return FetchOutcome(ok=False, reason="empty_response")
    name = filename or url.rstrip("/").rsplit("/", 1)[-1] or "download"
    return FetchOutcome(ok=True, payload=response.content, filename=name, origin=url)


def _public_get(url: str, *, params: Mapping[str, Any] | None = None,
                headers: Mapping[str, str] | None = None, timeout=DEFAULT_TIMEOUT,
                max_bytes: int = PUBLIC_MAX_BYTES) -> tuple[Any | None, str]:
    """Bounded public-reference request with strict redirect and provider policy."""
    requests = _requests()
    current = url
    for redirect_count in range(4):
        host = (urlparse(current).hostname or "").lower()
        if host not in PUBLIC_ALLOWED_HOSTS:
            return None, "unsafe_host"
        lock = _PROVIDER_LOCKS[host]
        for attempt in range(3):
            try:
                with lock:
                    response = requests.get(current, params=params, headers=dict(headers or {}),
                                            timeout=timeout, allow_redirects=False)
            except requests.exceptions.ConnectTimeout:
                return None, "network_unreachable:connect_timeout"
            except requests.exceptions.ConnectionError:
                return None, "network_unreachable:connection_error"
            except requests.exceptions.RequestException as exc:
                return None, f"network_error:{type(exc).__name__}"
            if response.status_code not in {429, 500, 502, 503, 504}:
                break
            if attempt == 2:
                return None, f"http_{response.status_code}"
            retry_after = str(getattr(response, "headers", {}).get("Retry-After", "0"))
            try:
                delay = min(max(float(retry_after), 0.0), 5.0)
            except ValueError:
                delay = 0.0
            if delay:
                time.sleep(delay)
        if response.status_code in {301, 302, 303, 307, 308}:
            location = str(getattr(response, "headers", {}).get("Location", ""))
            if not location or redirect_count == 3:
                return None, "unsafe_redirect"
            from urllib.parse import urljoin  # noqa: PLC0415
            current = urljoin(current, location)
            params = None
            continue
        if response.status_code == 304:
            return response, "not_modified"
        if response.status_code == 404:
            return None, "not_found"
        if response.status_code != 200:
            return None, f"http_{response.status_code}"
        content = getattr(response, "content", b"")
        declared = str(getattr(response, "headers", {}).get("Content-Length", ""))
        if declared.isdigit() and int(declared) > max_bytes:
            return None, "response_too_large"
        if len(content) > max_bytes:
            return None, "response_too_large"
        return response, ""
    return None, "unsafe_redirect"


def fetch_hpa_salivary_gland(_: str = "current", *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    response, reason = _public_get(HPA_BULK, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    if not response.content.startswith(b"PK"):
        return FetchOutcome(ok=False, reason="not_zip")
    return FetchOutcome(True, response.content, "proteinatlas.tsv.zip", response.url,
                        str(response.headers.get("Last-Modified", "")))


def fetch_rhea_release(_: str = "current", *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    response, reason = _public_get(RHEA_BULK, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    if not response.content.startswith(b"\x1f\x8b"):
        return FetchOutcome(ok=False, reason="not_gzip")
    return FetchOutcome(True, response.content, "rhea-tsv.tar.gz", response.url,
                        str(response.headers.get("Last-Modified", "")))


def fetch_proteomexchange_metadata(dataset_id: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    identifier = dataset_id.strip().upper()
    if not re.fullmatch(r"PXD[0-9]{6,}", identifier):
        return FetchOutcome(ok=False, reason="invalid_dataset_id")
    response, reason = _public_get(PROTEOMEXCHANGE_DATASET, params={"ID": identifier, "outputMode": "JSON"}, timeout=timeout,
                                   max_bytes=8 * 1024 * 1024)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    payload = response.content
    if identifier.encode() not in payload:
        return FetchOutcome(ok=False, reason="dataset_identity_mismatch")
    return FetchOutcome(True, payload, f"{identifier}.json", response.url,
                        str(response.headers.get("Last-Modified", "")))


def fetch_homd_release(download_url: str, *, timeout=DEFAULT_TIMEOUT) -> FetchOutcome:
    """Fetch an explicitly selected HOMD/eHOMD release URL from an allowlisted host."""
    response, reason = _public_get(download_url, timeout=timeout)
    if response is None:
        return FetchOutcome(ok=False, reason=reason)
    payload = response.content
    if not payload or payload.lstrip().lower().startswith(b"<html"):
        return FetchOutcome(ok=False, reason="not_release_data")
    name = urlparse(response.url).path.rsplit("/", 1)[-1] or "ehomd-release"
    return FetchOutcome(True, payload, name, response.url, str(response.headers.get("Last-Modified", "")))


# Retrieval strategies keyed by source_id, for the sources that need a bespoke
# call rather than a plain URL GET. Anything not listed here falls back to
# `fetch_url` using the registry entry's `url` field.
NAMED_FETCHERS: Mapping[str, Callable[..., FetchOutcome]] = {
    "uniprot_proteomes": fetch_uniprot_proteome,
    "human_proteome": fetch_uniprot_proteome,
    "commensal_panel": fetch_uniprot_proteome,
    "strain_panel": fetch_uniprot_proteome,
    "alphafold_db": fetch_alphafold_model,
    "pdb": fetch_pdb_structure,
    "human_protein_atlas": fetch_hpa_salivary_gland,
    "rhea": fetch_rhea_release,
    "proteomexchange": fetch_proteomexchange_metadata,
    "homd": fetch_homd_release,
}
