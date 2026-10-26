"""Read-only MARSE ecosystem frames summary; no ecological-to-molecular inference."""
import hashlib
import json
import math


def read_marse_context(text):
    if not isinstance(text, str) or len(text.encode()) > 8_000_000:
        raise ValueError('MARSE frames JSON must be text under 8 MB')
    document = json.loads(text)
    if not isinstance(document, dict):
        raise ValueError('Expected MARSE ecosystem frames object')
    species = document.get('species')
    frames = document.get('frames')
    if not isinstance(species, list) or not species or not all(isinstance(s,str) and s.strip() for s in species):
        raise ValueError('Expected named species list')
    if len(species) != len(set(species)):
        raise ValueError('Duplicate species names')
    if not isinstance(frames,list) or not frames or not isinstance(frames[-1],dict):
        raise ValueError('Expected nonempty ecosystem frames')
    def number(x):
        if isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) or x < 0:
            raise ValueError('MARSE context values must be finite nonnegative numbers')
        return x
    last = frames[-1]
    totals = last['statistics']['species_total_biomass']
    if not isinstance(totals,dict) or set(totals) != set(species):
        raise ValueError('Species totals do not match species list')
    cell_size = number(document['cell_size_um'])
    if not cell_size:
        raise ValueError('cell_size_um must be positive')
    return {'schema_version':'organism_orient.marse_context.v1',
            'source_sha256':hashlib.sha256(text.encode()).hexdigest(),
            'experiment_id':str(document['experiment_id']), 'species':species,
            'frame_count':len(frames), 'last_time_h':number(last['time_h']),
            'last_species_total_biomass':{s:number(totals[s]) for s in species},
            'grid_cell_size_um':cell_size,
            'integration_state':'read_only_context_import',
            'limitations':['MARSE biomass is not a receptor count', 'Grid cell size is not an organism radius',
                           'No membrane composition or molecular coordinates inferred',
                           'No simulation execution or feedback into MARSE']}
