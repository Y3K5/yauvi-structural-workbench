"""Retrieval logic, exercised without a network.

A fake `requests` module is injected so the parsing, validation, and
miss-reporting paths are all covered offline. Tests that need a real endpoint are
marked `network` and deselected by default.
"""
from __future__ import annotations

import sys
import types
import gzip

import pytest

from yauvi_sources.fetchers import (
    UNIPROT_ANNOTATION_FIELDS,
    proxy_configured,
    fetch_alphafold_model,
    fetch_alphafold_pae,
    fetch_pdb_assembly,
    fetch_pdb_structure,
    fetch_structural_artifact,
    fetch_uniprot_annotation,
    fetch_uniprot_proteome,
    fetch_uniprot_sequence,
    fetch_wwpdb_validation,
    fetch_url,
    fetch_homd_release,
    fetch_hpa_salivary_gland,
    fetch_proteomexchange_metadata,
    fetch_rhea_release,
    host_reachable,
)


class _Response:
    def __init__(self, *, status=200, content=b"", headers=None, json_data=None, url="https://x"):
        self.status_code = status
        self.content = content
        self.headers = headers or {}
        self._json = json_data
        self.url = url

    def json(self):
        if self._json is None:
            raise ValueError("not json")
        return self._json


class _Exceptions:
    class RequestException(Exception):
        pass

    class ConnectTimeout(RequestException):
        pass

    class ConnectionError(RequestException):
        pass


@pytest.fixture
def fake_requests(monkeypatch):
    """Install a stub `requests` that returns queued responses."""
    module = types.ModuleType("requests")
    module.exceptions = _Exceptions
    queue: list = []
    calls: list = []

    def get(url, params=None, timeout=None, **kwargs):
        calls.append({"url": url, "params": params or {}, "timeout": timeout, **kwargs})
        item = queue.pop(0)
        if isinstance(item, Exception):
            raise item
        return item

    module.get = get
    monkeypatch.setitem(sys.modules, "requests", module)
    module.queue = queue
    module.calls = calls
    return module


# --- UniProt -------------------------------------------------------------


def test_proteome_fetch_records_the_uniprot_release(fake_requests):
    fake_requests.queue.append(
        _Response(
            content=b">sp|P1|X\nMKV\n",
            headers={"X-UniProt-Release": "2026_03", "X-UniProt-Release-Date": "2026-06-18"},
            url="https://rest.uniprot.org/uniprotkb/stream?x",
        )
    )
    outcome = fetch_uniprot_proteome("UP000005640")
    assert outcome.ok
    assert outcome.filename == "UP000005640.fasta"
    # The registry records that UniProt has no pinned version. Capturing the
    # release header at retrieval time is what closes that gap.
    assert outcome.version == "2026_03 (2026-06-18)"


def test_proteome_fetch_rejects_a_non_fasta_body(fake_requests):
    fake_requests.queue.append(_Response(content=b"<html>error</html>"))
    outcome = fetch_uniprot_proteome("UP000005640")
    assert not outcome.ok and outcome.reason == "not_fasta"


def test_proteome_fetch_reports_a_404_without_raising(fake_requests):
    fake_requests.queue.append(_Response(status=404))
    outcome = fetch_uniprot_proteome("UP0000000")
    assert not outcome.ok and outcome.reason == "not_found"


def test_proteome_fetch_reports_an_unexpected_status(fake_requests):
    fake_requests.queue.extend([_Response(status=503), _Response(status=503), _Response(status=503)])
    assert fetch_uniprot_proteome("UP000005640").reason == "http_503"


def test_connection_failure_is_reported_not_raised(fake_requests):
    fake_requests.queue.append(_Exceptions.ConnectionError("no route"))
    outcome = fetch_uniprot_proteome("UP000005640")
    assert not outcome.ok
    assert outcome.reason == "network_unreachable:connection_error"


def test_annotation_fetch_requests_the_activity_state_fields(fake_requests):
    fake_requests.queue.append(_Response(content=b"Entry\tEC number\nP1\t1.1.1.1\n"))
    outcome = fetch_uniprot_annotation("taxonomy_id:837", filename="837.tsv")
    assert outcome.ok
    requested = fake_requests.calls[0]["params"]["fields"].split(",")
    # These four are what the activity-state module needs and no existing
    # annotation export in the tree contains.
    for field in ("ft_act_site", "ft_binding", "ft_site", "cc_cofactor"):
        assert field in requested, f"{field} not requested from UniProt"


def test_annotation_fetch_rejects_a_body_that_is_not_tsv(fake_requests):
    fake_requests.queue.append(_Response(content=b"not a table"))
    assert fetch_uniprot_annotation("q", filename="x.tsv").reason == "not_tsv"


# --- AlphaFold -----------------------------------------------------------


def test_alphafold_resolves_the_url_from_the_api(fake_requests):
    fake_requests.queue.append(
        _Response(json_data=[{"pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P1-F1-model_v6.pdb", "latestVersion": 6}])
    )
    fake_requests.queue.append(_Response(content=b"ATOM      1  N   MET A   1\n"))
    outcome = fetch_alphafold_model("P12345")
    assert outcome.ok
    assert outcome.filename == "AF-P12345-F1.pdb"
    assert outcome.version == "6"
    # The version-shaped path is never guessed; it comes from the API response.
    assert fake_requests.calls[1]["url"] == "https://alphafold.ebi.ac.uk/files/AF-P1-F1-model_v6.pdb"


def test_alphafold_reports_an_accession_with_no_model(fake_requests):
    fake_requests.queue.append(_Response(json_data=[]))
    assert fetch_alphafold_model("P12345").reason == "no_model"


def test_alphafold_rejects_a_body_that_is_not_a_structure(fake_requests):
    fake_requests.queue.append(_Response(json_data=[{"pdbUrl": "https://alphafold.ebi.ac.uk/files/x.pdb"}]))
    fake_requests.queue.append(_Response(content=b"<html>404</html>"))
    assert fetch_alphafold_model("P12345").reason == "not_a_structure"


def test_alphafold_handles_a_non_json_api_response(fake_requests):
    fake_requests.queue.append(_Response(content=b"oops"))
    assert fetch_alphafold_model("P12345").reason == "bad_api_json"


# --- PDB and generic -----------------------------------------------------


def test_pdb_fetch_normalises_the_identifier(fake_requests):
    fake_requests.queue.append(_Response(content=b"data_1ABC\n_entry.id 1ABC\n"))
    outcome = fetch_pdb_structure("1ABC")
    assert outcome.ok and outcome.filename == "1abc.cif"
    assert "1abc.cif" in fake_requests.calls[0]["url"]


def test_pdb_fetch_rejects_a_non_mmcif_body(fake_requests):
    fake_requests.queue.append(_Response(content=b"HEADER not cif"))
    assert fetch_pdb_structure("1abc").reason == "not_mmcif"


def test_structural_fetchers_are_artifact_specific_and_identifier_bounded(fake_requests):
    fake_requests.queue.append(_Response(content=b"data_4hhb\n_entry.id 4HHB\n", url="https://files.rcsb.org/download/4hhb-assembly1.cif"))
    assembly = fetch_pdb_assembly("4HHB:1")
    assert assembly.ok and assembly.filename == "4hhb-assembly1.cif"
    assert fetch_pdb_assembly("../../etc/passwd").reason == "invalid_pdb_assembly_identifier"
    assert fetch_structural_artifact("unknown.kind", "1ABC").reason == "unsupported_structural_artifact"


def test_wwpdb_validation_is_decompressed_and_checked(fake_requests):
    payload = b"<ValidationReport><Entry id='1crn'/></ValidationReport>"
    fake_requests.queue.append(_Response(content=gzip.compress(payload), url="https://files.rcsb.org/validation.xml.gz"))
    outcome = fetch_wwpdb_validation("1CRN")
    assert outcome.ok and outcome.payload == payload and outcome.filename == "1crn_validation.xml"


def test_uniprot_sequence_and_alphafold_pae_are_checked(fake_requests):
    fake_requests.queue.append(_Response(content=b">sp|P69905|HBA_HUMAN\nMVLSPADK\n", url="https://rest.uniprot.org/uniprotkb/P69905.fasta"))
    assert fetch_uniprot_sequence("P69905").ok
    fake_requests.queue.append(_Response(json_data=[{"paeDocUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-pae.json", "latestVersion": 6}]))
    fake_requests.queue.append(_Response(content=b'[{"predicted_aligned_error":[[0.0]]}]', json_data=[{"predicted_aligned_error": [[0.0]]}], url="https://alphafold.ebi.ac.uk/files/AF-P69905-F1-pae.json"))
    pae = fetch_alphafold_pae("P69905")
    assert pae.ok and pae.filename == "AF-P69905-F1-pae.json"


def test_alphafold_model_refuses_an_api_redirect_to_an_unregistered_host(fake_requests):
    fake_requests.queue.append(_Response(json_data=[{"pdbUrl": "https://attacker.invalid/model.pdb"}]))
    assert fetch_alphafold_model("P12345").reason == "unsafe_host"


def test_generic_url_fetch_derives_a_filename(fake_requests):
    fake_requests.queue.append(_Response(content=b"payload"))
    outcome = fetch_url("https://example.invalid/data/VFDB_setA.fas.gz")
    assert outcome.ok and outcome.filename == "VFDB_setA.fas.gz"


def test_generic_url_fetch_reports_an_empty_body(fake_requests):
    fake_requests.queue.append(_Response(content=b""))
    assert fetch_url("https://example.invalid/x").reason == "empty_response"


# --- Oral ecosystem public reference sources ----------------------------


def test_hpa_and_rhea_fetchers_check_archive_signatures(fake_requests):
    fake_requests.queue.append(_Response(content=b"PK\x03\x04data", url="https://www.proteinatlas.org/download/proteinatlas.tsv.zip"))
    fake_requests.queue.append(_Response(content=b"\x1f\x8bdata", url="https://ftp.expasy.org/databases/rhea/tsv/rhea-tsv.tar.gz"))
    assert fetch_hpa_salivary_gland("current").ok
    assert fetch_rhea_release("current").ok


def test_oral_public_fetchers_reject_unsafe_hosts_and_wrong_dataset_identity(fake_requests):
    assert fetch_homd_release("https://attacker.invalid/release.tsv").reason == "unsafe_host"
    fake_requests.queue.append(_Response(content=b'{"dataset":"PXD999999"}', url="https://proteomecentral.proteomexchange.org/cgi/GetDataset"))
    assert fetch_proteomexchange_metadata("PXD006367").reason == "dataset_identity_mismatch"


def test_homd_redirect_cannot_escape_allowlist(fake_requests):
    fake_requests.queue.append(_Response(status=302, headers={"Location": "https://attacker.invalid/file.tsv"}, url="https://v31a.homd.org/download"))
    assert fetch_homd_release("https://v31a.homd.org/download").reason == "unsafe_host"


def test_public_fetcher_bounds_declared_response_size(fake_requests):
    fake_requests.queue.append(_Response(content=b"PK", headers={"Content-Length": str(200 * 1024 * 1024)},
                                           url="https://www.proteinatlas.org/download/proteinatlas.tsv.zip"))
    assert fetch_hpa_salivary_gland("current").reason == "response_too_large"


# --- reachability --------------------------------------------------------


def test_unresolvable_host_is_reported_within_the_budget(monkeypatch):
    """The probe must bound name resolution, which no socket timeout covers."""
    import time

    monkeypatch.setattr("yauvi_sources.fetchers.http.proxy_configured", lambda: False)
    started = time.monotonic()
    assert host_reachable("no-such-host.invalid", timeout=2.0) is False
    assert time.monotonic() - started < 5.0


def test_the_probe_defers_to_the_proxy_when_one_is_configured(monkeypatch):
    """Behind a proxy the client never resolves the target host -- the proxy does.

    A direct getaddrinfo then fails even though requests succeed, so trusting it
    would refuse to fetch on exactly the networks where fetching works: CI
    runners, and anything corporate.
    """
    monkeypatch.setattr("yauvi_sources.fetchers.http.proxy_configured", lambda: True)
    assert host_reachable("no-such-host.invalid", timeout=2.0) is True


def test_proxy_configured_reads_the_environment(monkeypatch):
    monkeypatch.setattr("urllib.request.getproxies", lambda: {})
    assert proxy_configured() is False
    monkeypatch.setattr("urllib.request.getproxies", lambda: {"https": "http://p:8080"})
    assert proxy_configured() is True


def test_annotation_field_list_covers_the_existing_export():
    """The columns the current per-taxon TSVs carry must not be dropped."""
    for field in ("accession", "protein_name", "gene_names", "xref_interpro", "xref_pfam", "ec"):
        assert field in UNIPROT_ANNOTATION_FIELDS


@pytest.mark.network
def test_real_uniprot_release_header_is_present():
    outcome = fetch_uniprot_proteome("UP000000625")
    assert outcome.ok, outcome.reason
    assert outcome.version, "UniProt stopped reporting its release header"
