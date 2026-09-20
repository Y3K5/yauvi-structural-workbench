"""Typed public reference discovery; provider assertions stay distinct from validation."""
from __future__ import annotations

import hashlib
import io
import re
from urllib.parse import urlencode, urlsplit

PROTEOME = re.compile(r'UP[0-9]{9}')


def proteome_identifier(value):
    """Return a proteome ID, or None so the strict protein parser can handle it."""
    if not isinstance(value, str) or len(value) > 500:
        return None
    text = value.strip()
    if '://' in text:
        url = urlsplit(text)
        if (url.scheme != 'https' or url.hostname not in {'uniprot.org', 'www.uniprot.org', 'rest.uniprot.org'}
                or url.username or url.password or url.port not in (None, 443) or url.query or url.fragment):
            return None
        match = re.fullmatch(r'/proteomes/(UP[0-9]{9})(?:/entry|\.json)?/?', url.path, re.I)
        return match[1].upper() if match else None
    return text.upper() if PROTEOME.fullmatch(text.upper()) else None


def resource(kind, identifier, label, *, context=None, url=None):
    return {'resource_id': kind + ':' + identifier, 'artifact_type': kind,
            'identifier': identifier, 'label': label, 'context': context or {}, 'url': url}


def protein_context(entry, entry_sha256):
    organism = entry.get('organism', {})
    context = {'taxon_id': organism.get('taxonId'), 'organism': organism.get('scientificName'),
               'lineage': organism.get('lineage', []), 'entry_type': entry.get('entryType'),
               'entry_audit': entry.get('entryAudit', {}), 'genes': entry.get('genes', []),
               'sequence_sha256': hashlib.sha256(entry['sequence']['value'].encode()).hexdigest(),
               'evidence_basis': 'UniProt provider assertions', 'source_sha256': entry_sha256}
    resources = []
    seen = set()
    for ref in entry.get('uniProtKBCrossReferences', []):
        identifier = str(ref.get('id', '')).upper()
        key = (ref.get('database'), identifier)
        if key in seen:
            continue
        seen.add(key)
        properties = {p['key']: p.get('value') for p in ref.get('properties', []) if 'key' in p}
        if ref.get('database') == 'Proteomes' and PROTEOME.fullmatch(identifier):
            resources.append(resource('uniprot.proteome', identifier, identifier + ' · whole proteome FASTA',
                context={'component': properties.get('Component'), 'taxon_id': organism.get('taxonId'),
                         'relationship': 'UniProt entry cross-reference; proteome metadata checked on retrieval',
                         'source_sha256': entry_sha256}, url=f'https://www.uniprot.org/proteomes/{identifier}'))
        if ref.get('database') == 'PDB' and re.fullmatch(r'[0-9][A-Z0-9]{3}', identifier):
            detail = {'method': properties.get('Method'), 'resolution': properties.get('Resolution'),
                      'chain_coverage': properties.get('Chains'), 'source_sha256': entry_sha256,
                      'relationship': 'UniProt PDB cross-reference; exact chain sequence and strain not verified',
                      'assembly': 'deposited entry / asymmetric unit; biological assembly not selected'}
            resources.append(resource('pdb.coordinates', identifier, identifier + ' · experimental coordinates (mmCIF)',
                context=detail, url=f'https://www.rcsb.org/structure/{identifier}'))
            resources.append(resource('pdb.legacy', identifier, identifier + ' · legacy PDB coordinates',
                context=detail, url=f'https://www.rcsb.org/structure/{identifier}'))
    return context, resources


def proteome_context(doc, identifier, expected_taxon=None):
    if doc.get('id') != identifier:
        raise ValueError('Proteome metadata returned a different identifier.')
    taxonomy = doc.get('taxonomy', {})
    taxon = taxonomy.get('taxonId')
    if not isinstance(taxon, int) or not taxonomy.get('scientificName'):
        raise ValueError('Proteome taxonomy is unavailable; context cannot be checked.')
    if expected_taxon is not None and taxon != expected_taxon:
        raise ValueError('Proteome taxon differs from the selected protein; review the organism/strain relationship separately.')
    count = doc.get('proteinCount')
    if type(count) is not int or count < 1:
        raise ValueError('Proteome protein count is unavailable.')
    return {'proteome_id': identifier, 'organism': taxonomy['scientificName'], 'taxon_id': taxon,
            'strain': doc.get('strain'), 'proteome_type': doc.get('proteomeType', 'unspecified'),
            'protein_count': count, 'genome_assembly': doc.get('genomeAssembly', {}),
            'components': doc.get('components', []), 'evidence_basis': 'UniProt proteome metadata',
            'isoform_policy': 'canonical entries only; includeIsoform=false'}


def proteome_url(identifier):
    if not PROTEOME.fullmatch(identifier):
        raise ValueError('Invalid proteome identifier.')
    return 'https://rest.uniprot.org/uniprotkb/stream?' + urlencode(
        {'query': f'proteome:{identifier}', 'format': 'fasta', 'includeIsoform': 'false'})


def validate_proteome(body, context, accession=None, sequence=None):
    """Check every FASTA record, count, taxon and optional exact protein membership."""
    from Bio import SeqIO
    text = body.decode('utf-8')
    if not text.startswith('>'):
        raise ValueError('Expected proteome FASTA, not a page or error response.')
    ids = set()
    found = False
    for record in SeqIO.parse(io.StringIO(text), 'fasta'):
        match = re.fullmatch(r'(?:sp|tr)\|([^|]+)\|[^|]+', record.id)
        seq = str(record.seq)
        taxon = re.search(r'\bOX=(\d+)\b', record.description)
        if not match or not seq or not re.fullmatch('[A-Z]+', seq) or not taxon:
            raise ValueError('Proteome contains an invalid record or missing taxon identity.')
        rid = match[1]
        if rid in ids or '-' in rid or int(taxon[1]) != context['taxon_id']:
            raise ValueError('Proteome contains duplicate/isoform entries or a different taxon.')
        ids.add(rid)
        if rid == accession:
            if seq != sequence:
                raise ValueError('Proteome sequence differs from the locked protein entry.')
            found = True
    if len(ids) != context['protein_count']:
        raise ValueError('Proteome FASTA count differs from the metadata; incomplete or changed release.')
    if accession and not found:
        raise ValueError('Selected protein is absent from its linked proteome.')
    return {'record_count': len(ids), 'taxon_check': 'every FASTA record matched',
            'protein_membership': 'exact accession and sequence matched' if accession else 'not applicable'}


def validate_pdb(body, identifier, legacy=False):
    """Check deposited file identity; this deliberately does not infer chain equivalence."""
    if legacy:
        text = body.decode('utf-8')
        header = next((line for line in text.splitlines() if line.startswith('HEADER')), '')
        if header[62:66].strip().upper() != identifier or not any(line.startswith('ATOM  ') for line in text.splitlines()):
            raise ValueError('PDB file identity or protein coordinates are missing.')
        return {'entry_id': identifier, 'format': 'legacy PDB', 'chain_mapping': 'requires review'}
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    doc = MMCIF2Dict(io.StringIO(body.decode('utf-8')))
    if doc.get('_entry.id', [''])[0].upper() != identifier or not doc.get('_atom_site.Cartn_x'):
        raise ValueError('mmCIF entry identity or coordinates do not match.')
    return {'entry_id': identifier, 'format': 'PDBx/mmCIF', 'methods': doc.get('_exptl.method', []),
            'source_taxa': doc.get('_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id', [])
                + doc.get('_entity_src_nat.pdbx_ncbi_taxonomy_id', []),
            'revision_dates': doc.get('_pdbx_audit_revision_history.revision_date', []),
            'chain_mapping': 'requires review', 'assembly_selection': 'not performed'}
