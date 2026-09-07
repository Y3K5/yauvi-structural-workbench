"""Explicit public-accession acquisition, followed by reviewed case adoption.

Only a UniProt accession leaves the machine. Provider responses and derived
inputs are retained by hash. No annotations are promoted to measured biology.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import re
import secrets
from pathlib import Path
from urllib.parse import urlsplit

from .sources import StructuralSourceStore, StructuralSourceError
from .store import AnalysisError, _write_json

ACCESSION = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})')
IMPORT_ID = re.compile(r'protein_[0-9a-f]{24}')


def parse_uniprot_link(value):
    if not isinstance(value, str) or len(value) > 500:
        raise ValueError('Enter a UniProt entry link or accession, such as P02931.')
    value = value.strip()
    if '://' in value:
        url = urlsplit(value)
        if url.scheme != 'https' or url.hostname not in {'uniprot.org', 'www.uniprot.org', 'rest.uniprot.org'} or url.username or url.password or url.port not in {None, 443}:
            raise ValueError('Use an https UniProt entry link; other hosts are not accepted.')
        match = re.fullmatch(r'/(?:uniprotkb|uniprot)/([A-Za-z0-9]+)(?:/entry|\.fasta|\.json)?/?', url.path)
        if not match:
            raise ValueError('Use a single UniProt entry link, not a search or proteome page.')
        value = match[1]
    value = value.upper()
    if not ACCESSION.fullmatch(value):
        raise ValueError('Enter a canonical UniProt accession. Isoforms need explicit sequence and structure selection.')
    return value


class ProteinImportStore:
    def __init__(self, store, *, get=None):
        from yauvi_sources.fetchers.http import _public_get
        self.store = store
        self.sources = StructuralSourceStore(store.workspace)
        self.root = self.sources.root / 'protein-imports'
        self.root.mkdir(exist_ok=True)
        self.get = get or _public_get

    def load(self, identifier):
        if not isinstance(identifier, str) or not IMPORT_ID.fullmatch(identifier):
            raise ValueError('Invalid protein import identifier.')
        return json.loads((self.root / identifier / 'IMPORT.json').read_text())

    def save(self, record):
        _write_json(self.root / record['import_id'] / 'IMPORT.json', record)
        return record

    def download(self, kind, accession, url, name, release='', limit=8*1024*1024):
        from yauvi_sources.fetchers.http import FetchOutcome
        response, reason = self.get(url, timeout=(8, 30), max_bytes=limit)
        if response is None:
            raise StructuralSourceError(f'{kind}: {reason}')
        return self.sources.acquire(kind, accession, prepared_outcome=FetchOutcome(
            True, response.content, name, response.url, release or response.headers.get('X-UniProt-Release', 'unknown')))

    def document(self, acquisition):
        return json.loads(self.sources.artifact_path(acquisition['acquisition_id']).read_text())

    def lookup(self, link):
        accession = parse_uniprot_link(link)
        entry_artifact = self.download('uniprot.entry', accession,
            f'https://rest.uniprot.org/uniprotkb/{accession}.json', f'{accession}-entry.json')
        entry = self.document(entry_artifact)
        if entry.get('primaryAccession') != accession:
            raise ValueError('UniProt returned a different or redirected accession. Use its current canonical entry explicitly.')
        sequence = entry.get('sequence', {}).get('value', '')
        if not re.fullmatch('[A-Z]+', sequence) or len(sequence) > 20000:
            raise ValueError('The UniProt sequence is missing or exceeds the supported import size.')
        description = entry.get('proteinDescription', {})
        recommended = description.get('recommendedName') or next(iter(description.get('submissionNames', [])), {})
        record = {'schema_version':'1.0', 'import_id':'protein_'+secrets.token_hex(12),
            'accession':accession, 'name':recommended.get('fullName', {}).get('value', accession),
            'organism':entry.get('organism', {}).get('scientificName', ''), 'length':len(sequence),
            'sequence':sequence, 'state':'identified', 'metadata':[entry_artifact], 'files':[], 'warnings':[],
            'pdb_ids':[r['id'] for r in entry.get('uniProtKBCrossReferences', []) if r.get('database') == 'PDB' and re.fullmatch('[A-Za-z0-9]{4}', r.get('id',''))],
            'limitations':['AlphaFold coordinates are predicted; confidence is not experimental validation.',
                'Assemblies, state references, curated site mappings, validation and benchmark universes require separate evidence.']}
        return self.save(record)

    def prepare(self, identifier):
        from yauvi_sources.fetchers.http import FetchOutcome
        record = self.load(identifier)
        if record['state'] == 'prepared':
            return record
        if record['state'] != 'identified':
            raise ValueError('This import was interrupted; start a new lookup to preserve its record.')
        record['state'] = 'retrieving'
        self.save(record)
        accession, sequence = record['accession'], record['sequence']
        def add(kind, outcome, role):
            acquired = self.sources.acquire(kind, accession, prepared_outcome=outcome)
            record['files'].append({'role':role, 'acquisition':acquired})
            return acquired
        entry = record['metadata'][0]['artifact']
        # Derive FASTA from the exact locked UniProt JSON, avoiding release drift.
        fasta = f'>{accession}\n'+'\n'.join(sequence[i:i+60] for i in range(0,len(sequence),60))+'\n'
        add('uniprot.sequence', FetchOutcome(True, fasta.encode(), f'{accession}.fasta',
            entry['origin']+'#sequence.value', entry['release']), 'reference_fasta')
        try:
            metadata = self.download('alphafold.metadata', accession,
                f'https://alphafold.ebi.ac.uk/api/prediction/{accession}', f'{accession}-alphafold.json')
            record['metadata'].append(metadata)
            entries = self.document(metadata)
            candidates = [r for r in entries if isinstance(r,dict) and r.get('uniprotAccession') == accession
                and r.get('sequence', r.get('uniprotSequence')) == sequence
                and r.get('sequenceStart',r.get('uniprotStart')) == 1
                and r.get('sequenceEnd',r.get('uniprotEnd')) == len(sequence)] if isinstance(entries,list) else []
            if len(candidates) != 1:
                raise ValueError('No single full-length AlphaFold model matches this UniProt sequence. Fragmented or ambiguous models require manual selection.')
            af = candidates[0]
            entity = af.get('modelEntityId',af.get('entryId',''))
            version = str(af.get('latestVersion',''))
            if entity != f'AF-{accession}-F1' or not version.isdigit():
                raise ValueError('The provider model identity/version needs manual review.')
            model_url = self.model_url(af.get('pdbUrl'), entity, version, 'model', 'pdb')
            model = self.download('alphafold.model',accession,model_url,f'{entity}-model_v{version}.pdb',version,32*1024*1024)
            model_path = self.sources.artifact_path(model['acquisition_id'])
            chain = self.verify_structure(model_path.read_bytes(),sequence)
            record['files'].append({'role':'structure','acquisition':model})
            record['chain'] = chain
            provenance = {'class':'predicted','method':'AlphaFold DB','source_id':entity,
                'confidence_encoding':'plddt_in_bfactor', 'coordinate_sha256':model['artifact']['sha256'],
                'source_url':model_url,'release':version,'uniprot_accession':accession,
                'uniprot_entry_sha256':entry['sha256'],'model_metadata_sha256':metadata['artifact']['sha256'],
                'limitations':record['limitations']}
            path = self.root/identifier/f'{accession}-provenance.json'
            _write_json(path,provenance)
            record['provenance']={'file_name':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
            try:
                pae_url = self.model_url(af.get('paeDocUrl'), entity, version, 'predicted_aligned_error', 'json')
                pae = self.download('alphafold.pae',accession,pae_url,f'{entity}-pae_v{version}.json',version,64*1024*1024)
                doc = self.document(pae)
                matrix = (doc[0] if isinstance(doc,list) and len(doc)==1 else doc).get('predicted_aligned_error')
                if not isinstance(matrix,list) or len(matrix)!=len(sequence) or any(not isinstance(row,list) or len(row)!=len(sequence) or any(type(v) not in (int,float) or not math.isfinite(v) or v<0 for v in row) for row in matrix):
                    raise ValueError('PAE dimensions or values do not match this model.')
                record['files'].append({'role':'pae','acquisition':pae})
            except (StructuralSourceError,ValueError,AttributeError,TypeError) as exc:
                record['warnings'].append('PAE unavailable: '+str(exc))
        except (StructuralSourceError,ValueError,TypeError,KeyError) as exc:
            record['warnings'].append('Structure unavailable: '+str(exc))
        record['state']='prepared'
        return self.save(record)

    @staticmethod
    def model_url(url, entity, version, kind, suffix):
        expected=f'https://alphafold.ebi.ac.uk/files/{entity}-{kind}_v{version}.{suffix}'
        if url != expected:
            raise ValueError('Provider file URL does not match the selected model identity and version.')
        return url

    @staticmethod
    def verify_structure(body, sequence):
        from Bio.PDB import PDBParser
        from Bio.SeqUtils import seq1
        structure=PDBParser(QUIET=True).get_structure('import',io.StringIO(body.decode()))
        models=list(structure)
        if len(models)!=1 or len(list(models[0]))!=1:
            raise ValueError('Automatic import requires one model and one chain.')
        chain=next(iter(models[0]));residues=list(chain)
        if ''.join(seq1(r.resname) for r in residues)!=sequence or any(r.id != (' ',i,' ') for i,r in enumerate(residues,1)):
            raise ValueError('Downloaded coordinates do not exactly match the UniProt sequence and numbering.')
        return chain.id

    def adopt(self, identifier, analysis_id, revision):
        record=self.load(identifier)
        if record['state']!='prepared':
            raise ValueError('Retrieve the protein files before using them.')
        with self.store._lock:
            case=self.store.load(analysis_id)
            if case['revision_sha256']!=revision:
                raise AnalysisError('Inputs changed during retrieval. Review the current case before trying again.')
            definition=self.store._definitions[case['analysis_type']]
            roles={r['role'] for r in definition['inputs']}
            sf=case['analysis_type']=='sf_csa'
            mapped=lambda role: {'structure':'query_structure','reference_fasta':'query_fasta'}.get(role,role) if sf else role
            files=[(mapped(row['role']),row['acquisition']) for row in record['files'] if mapped(row['role']) in roles]
            existing={r['role']:r for r in case['inputs']}
            # Never bind new confidence/provenance to unrelated or stale coordinates.
            coordinate=next((a for role,a in files if role in {'structure','query_structure'}),None)
            for role,acquisition in files:
                if role in existing and existing[role]['sha256']!=acquisition['artifact']['sha256']:
                    raise AnalysisError('An attached input differs from this protein bundle. Create a new analysis to keep the original evidence intact.')
                self.sources.artifact_path(acquisition['acquisition_id'])
            if not coordinate and any(r in existing for r in ('structure','query_structure')):
                raise AnalysisError('No matching downloaded model is available to verify your existing structure. Use a new analysis or review files manually.')
            provenance=record.get('provenance')
            if provenance:
                path=self.root/identifier/provenance['file_name']
                if hashlib.sha256(path.read_bytes()).hexdigest()!=provenance['sha256']:
                    raise AnalysisError('Provenance checksum mismatch.')
                if 'provenance' in existing and existing['provenance']['sha256']!=provenance['sha256']:
                    raise AnalysisError('This analysis already has a different provenance declaration. Use a new analysis.')
            attached=[]
            for role,acquisition in files:
                if role not in existing:
                    self.store.adopt_source_artifact(analysis_id,role=role,artifact_manifest=acquisition['artifact'],path=self.sources.artifact_path(acquisition['acquisition_id']))
                    attached.append(role)
            if provenance and 'provenance' in roles and 'provenance' not in existing:
                self.store.add_file(analysis_id,role='provenance',path=path);attached.append('provenance')
            parameters=dict(case['parameters']);allowed={p['name'] for p in definition.get('parameters',[])}
            for key,value in {'chain':record.get('chain'),'model':0 if coordinate else None,'accession':record['accession'],'organism':record['organism']}.items():
                if key in allowed and value is not None and parameters.get(key) in (None,''):
                    parameters[key]=value
            if parameters!=case['parameters']:
                self.store.update_parameters(analysis_id,parameters)
            latest=self.store.load(analysis_id)
            latest['protein_import']={k:record[k] for k in ('import_id','accession','name','organism','metadata','warnings','limitations','pdb_ids')}
            latest['revision']+=1
            self.store._commit(self.store._case_dir(analysis_id),latest)
            missing=[r['label'] for r in definition['inputs'] if not any(i['role']==r['role'] for i in latest['inputs'])]
            return {'attached':attached,'missing':missing,'warnings':record['warnings'],'analysis':self.store.load(analysis_id)}
