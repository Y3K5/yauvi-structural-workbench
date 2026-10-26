"""Explicit-lipid preparation contract. Never packs lipids or runs a builder."""
import math
from .surfaces import Surface


def plan_lipids(document):
    if not isinstance(document, dict):
        raise ValueError('Lipid recipe must be an object')
    if document.get('schema_version') != 'organism_orient.lipid_recipe.v1':
        raise ValueError('Expected organism_orient.lipid_recipe.v1')
    Surface(**document['surface'])
    missing = []
    for key in ('representation', 'force_field', 'composition_source', 'builder', 'builder_version', 'template_manifest_sha256'):
        if not isinstance(document.get(key),str) or not document[key].strip():
            missing.append(key)
    for side in ('positive','negative'):
        entries = document.get('leaflets',{}).get(side,[])
        if not isinstance(entries,list):
            raise ValueError('Leaflet compositions must be lists')
        if not entries:
            missing.append(f'leaflets.{side}')
            continue
        ids=[];total=0.
        for entry in entries:
            lipid=entry['lipid_id'];fraction=entry['mole_fraction']
            if not isinstance(lipid,str) or not lipid.strip() or lipid in ids:
                raise ValueError('Lipid IDs must be nonempty and unique within each leaflet')
            if isinstance(fraction,bool) or not isinstance(fraction,(int,float)) or not math.isfinite(fraction) or not 0 < fraction <= 1:
                raise ValueError('Mole fractions must be finite and in (0,1]')
            ids.append(lipid);total+=fraction
        if not math.isclose(total,1.,abs_tol=1e-6):
            raise ValueError(f'{side} leaflet mole fractions must sum to 1')
    return {'schema_version':'organism_orient.lipid_plan.v1', 'recipe':document,
            'state':'held_no_packing_adapter', 'missing_inputs':missing,
            'next_checks':['Review organism/strain/environment-specific composition evidence',
                           'Verify lipid templates and force-field compatibility against pinned builder',
                           'Implement and test a local packing adapter for the selected morphology',
                           'Check overlaps, topology, charge, solvation and minimization before MD'],
            'limitations':['No molecular lipids generated','Composition completeness is not biological validity',
                           'Geometric interfaces do not establish bilayer equilibrium thickness',
                           'No archaeal spanning-monolayer or cell-wall representation']}
