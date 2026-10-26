"""Offline source-bound organisms, coordinates and schematic envelope layers."""
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
from .surfaces import Surface

HOME=Path(__file__).with_name('catalogue')
CACHE=Path(os.environ['MEMORIENT_CACHE']).expanduser().resolve() if os.environ.get('MEMORIENT_CACHE') else None


def configure(cache=None):
    """Select a caller-owned coordinate cache; never fetch or write on import."""
    global CACHE
    value=cache if cache is not None else os.environ.get('MEMORIENT_CACHE')
    CACHE=Path(value).expanduser().resolve() if value is not None else None
    load_catalogue.cache_clear()

@lru_cache(maxsize=1)
def load_catalogue():
    if CACHE is not None or not (HOME/'starter.json').is_file():
        document=json.loads((HOME/'metadata.json').read_text())
        for protein in document['proteins']:
            protein['availability']='missing_local_coordinates'
            if CACHE is None: continue
            path=CACHE/(protein['view_pdb_sha256']+'.pdb')
            if not path.is_file(): continue
            raw=path.read_bytes()
            if hashlib.sha256(raw).hexdigest()!=protein['view_pdb_sha256']:
                protein['availability']='checksum_mismatch'
                continue
            protein['pdb']=raw.decode('utf-8')
            protein['availability']='available'
        return document
    document=json.loads((HOME/'starter.json').read_text())
    for source in document['sources']:
        relative=Path(source['file'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Catalogue source path escapes its home')
        if hashlib.sha256((HOME/relative).read_bytes()).hexdigest()!=source['sha256']:
            raise ValueError(f'Catalogue source identity changed: {source["file"]}')
    for protein in document['proteins']:
        if hashlib.sha256(protein['pdb'].encode()).hexdigest()!=protein['view_pdb_sha256']:
            raise ValueError('Catalogue view coordinate identity mismatch')
        protein['availability']='available'
    return document


def catalogue_index():
    document=load_catalogue()
    return {**document,'proteins':[{k:v for k,v in p.items() if k!='pdb'} for p in document['proteins']]}


def protein_record(identifier):
    for record in load_catalogue()['proteins']:
        if record['id']==identifier:
            if 'pdb' not in record:
                raise ValueError('Coordinates unavailable: supply the checksum-matched view PDB in the local cache, or import your own PDB')
            return deepcopy(record)
    raise KeyError('Unknown starter protein')


def organism_overview(identifier):
    for organism in load_catalogue()['organisms']:
        if organism['id']!=identifier:continue
        layers=[]
        for layer in organism['layers']:
            surface={**organism['overview_surface']}
            surface['size']=[x+layer['offset_angstrom'] for x in surface['size']]
            layers.append({**layer,'surface':Surface(**surface).sample(24,13)})
        return {'organism':organism,'layers':layers,'sampling':{'nu':24,'nv':13},
                'units':'angstrom','evidence_state':organism['overview_evidence_state'],
                'protein_locations':'unknown; no locations or abundance inferred'}
    raise KeyError('Unknown starter organism')
