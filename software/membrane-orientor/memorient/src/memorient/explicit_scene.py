"""Read-only access to the checked, prebuilt K-12 explicit coordinate model."""
from functools import lru_cache
import hashlib
import json
from .explicit_outer_membrane import HOME

MODEL=HOME/'model-v2'
FILES={'ECOLI_K12_OM.cif':'coordinate_sha256','LIPID_BONDS.tsv':'lipid_bonds_sha256','ATOM_MAP.tsv':'atom_map_sha256'}

@lru_cache(maxsize=1)
def load_explicit_scene():
    scene=json.loads((MODEL/'scene.json').read_text())
    receipt=json.loads((MODEL/'BUILD_RECEIPT.json').read_text())
    if any(scene.get(key)!=value for key,value in receipt.items()):
        raise ValueError('Explicit scene metadata disagrees with its receipt')
    for name,key in FILES.items():
        if hashlib.sha256((MODEL/name).read_bytes()).hexdigest()!=receipt[key]:
            raise ValueError('Explicit model output identity mismatch')
    for name,expected in receipt['source_sha256'].items():
        if name not in ('K12-AA-rep1.pdb','K12-AA-rep1.psf','zenodo-18009953.json'):
            raise ValueError('Unexpected source in explicit model receipt')
        if hashlib.sha256((HOME/'sources'/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Explicit model source identity mismatch')
    for part in scene['parts']:
        if hashlib.sha256(part['pdb'].encode()).hexdigest()!=receipt['view_part_sha256'][part['id']]:
            raise ValueError('Explicit model display identity mismatch')
    return scene


def model_file(name):
    if name not in FILES:raise KeyError('Unknown explicit model output')
    load_explicit_scene()
    return (MODEL/name).read_bytes()
