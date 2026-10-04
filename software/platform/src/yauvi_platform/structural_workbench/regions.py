"""Versioned, identity-bound local regions and observed heavy-atom contacts.

Annotations, measured proximity and literature context remain separate. No
pocket discovery, interaction chemistry, dynamics or functional score is inferred.
"""
from __future__ import annotations

import hashlib
import json
import math

import numpy as np
from scipy.spatial import cKDTree

CONTACT_CUTOFF_A = 5.0  # AssemblyContext's existing heavy-atom contact convention.
AA = dict(zip(
    ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'],
    'ARNDCQEGHILKMFPSTWYV'))
METHOD = 'coordinate_bound_regions_and_heavy_atom_contacts_v1'
BASES = {'mapped_annotation', 'measured_contact', 'literature_context', 'declared_region'}


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False,
                       allow_nan=False) + '\n').encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def view_binding(view):
    """A selection binds atom identities/coordinates as well as its header."""
    if view.get('atom_identity_sha256') != digest([atom_identity(a) for a in view['atoms']]):
        raise ValueError('View atom identities or coordinates differ from their recorded binding')
    return {key: view[key] for key in ('structure_id','coordinate_sha256','sequence_sha256',
            'case_manifest_sha256','source_identities_sha256',
            'reference_record_sha256','mapping_sha256','method_sha256','model_id','assembly_id',
            'frame','conformer_selection','atom_identity_sha256')}


def select_region(view, region_id, binding):
    if binding != view_binding(view):
        raise ValueError('Region binding differs from selected coordinates/model/assembly/conformer/frame')
    matches = [region for region in view['regions'] if region['id'] == region_id]
    if len(matches) != 1 or matches[0]['binding'] != binding:
        raise ValueError('Region does not belong to this exact view')
    return matches[0]


def atom_identity(atom):
    return {key: atom.get(key) for key in ('source_atom_id','chain','resn','atom','elem','resi',
            'icode','label_seq_id','label_alt_id','occupancy','x','y','z')}


def measure_distance(view, first, second, binding):
    if binding != view_binding(view):
        raise ValueError('Measurement binding differs from selected view')
    if any(type(index) is not int or not 0 <= index < len(view['atoms']) for index in (first, second)):
        raise ValueError('Select represented atoms from this view')
    atoms = [view['atoms'][index] for index in (first, second)]
    return {'schema_version':'1.0', 'binding':binding, 'method':'euclidean_atom_distance_v1',
            'atom_indices':[first,second], 'atoms':[atom_identity(atom) for atom in atoms],
            'distance_A':math.dist([atoms[0][k] for k in 'xyz'], [atoms[1][k] for k in 'xyz']),
            'interpretation_limit':'Geometric distance only; no bond, affinity or biological interaction is assigned.'}


def _rows(doc, prefix):
    fields = [key[len(prefix):] for key in doc if key.startswith(prefix)]
    columns = [doc[prefix + field] for field in fields]
    if not columns:
        return []
    if len({len(column) for column in columns}) != 1:
        raise ValueError('Incomplete deposited category: ' + prefix)
    return [dict(zip(fields, row)) for row in zip(*columns)]


def _connectivity(view, doc):
    """Only deposited component bonds and short consecutive peptide links."""
    components = [r for c in view['chains'] for r in c['residues']]
    bonds = {}
    for row in _rows(doc, '_chem_comp_bond.'):
        bonds.setdefault(row['comp_id'], []).append(row)
    def connect(i, j, order):
        if i == j or j in view['atoms'][i]['bonds']:
            return
        for left, right in ((i,j),(j,i)):
            view['atoms'][left]['bonds'].append(right)
            view['atoms'][left]['bondOrder'].append(order)
    by_residue = {}
    for residue in components:
        names = {view['atoms'][i]['atom']: i for i in residue['atom_indices']}
        if len(names) != len(residue['atom_indices']):
            raise ValueError('Ambiguous atom identity within selected component conformer')
        by_residue[residue['id']] = names
        for row in bonds.get(residue['comp_id'], []):
            a, b = row['atom_id_1'], row['atom_id_2']
            if a in names and b in names:
                connect(names[a], names[b], {'SING':1,'DOUB':2,'TRIP':3}.get(row.get('value_order'), 1))
    for chain in view['chains']:
        polymer = {r['label_seq_id']:r for r in chain['residues'] if r['label_seq_id'] is not None}
        for position, residue in polymer.items():
            following = polymer.get(position+1)
            if not following:
                continue
            i = by_residue[residue['id']].get('C')
            j = by_residue[following['id']].get('N')
            if i is not None and j is not None:
                if math.dist([view['atoms'][i][k] for k in 'xyz'], [view['atoms'][j][k] for k in 'xyz']) < 2:
                    connect(i,j,1)
    # Keep reacted/adduct links explicitly declared in deposited coordinates.
    for row in _rows(doc, '_struct_conn.'):
        if row.get('conn_type_id') != 'covale' or any(
                row.get(f'ptnr{n}_symmetry', '1_555') not in ('1_555','?','.') for n in (1,2)):
            continue
        for operator in sorted({chain['operator_ids'][0] for chain in view['chains']}):
            partners=[]
            for n in (1,2):
                partners.append([a['index'] for a in view['atoms'] if
                    a['chain'] == f"{row.get(f'ptnr{n}_label_asym_id')}:{operator}" and
                    a['atom'] == row.get(f'ptnr{n}_label_atom_id') and
                    a['resn'] == row.get(f'ptnr{n}_label_comp_id') and
                    str(a['resi']) == row.get(f'ptnr{n}_auth_seq_id') and
                    a['icode'] == ('' if row.get(f'pdbx_ptnr{n}_PDB_ins_code', '?') in ('?','.') else row.get(f'pdbx_ptnr{n}_PDB_ins_code'))])
            if len(partners[0]) == len(partners[1]) == 1:
                connect(partners[0][0], partners[1][0], 1)


def _region(binding, region_id, title, basis, source_id, **fields):
    return {'id':region_id, 'binding':binding, 'title':title, 'basis':basis, 'source_id':source_id,
            'atom_indices':[], 'ligand_atom_indices':[], 'residue_ids':[],
            'mapping_state':'unmapped', 'coverage':{'observed':0,'expected':None},
            'limitations':[], **fields}


def build_regions(case, view, doc):
    """Called only after case identity, sequence, assembly and conformer checks."""
    view['schema_version']='1.0'
    view['frame']='deposited_coordinates'
    view['reference_record_sha256']=case.sources[case.protein['source_id']]['sha256']
    view['case_manifest_sha256']=hashlib.sha256((case.directory/'case.json').read_bytes()).hexdigest()
    view['source_identities_sha256']=digest({sid:source['sha256'] for sid,source in case.sources.items()})
    view['mapping_sha256']=digest(case.structures[view['structure_id']][0].get('sequence_mappings',[]))
    from pathlib import Path
    method_sources={'regions_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    'adapter_sha256':hashlib.sha256((Path(__file__).parent/'biological_case.py').read_bytes()).hexdigest()}
    view['method_sha256']=digest(method_sources)
    view['atom_identity_sha256']=digest([atom_identity(atom) for atom in view['atoms']])
    binding = view_binding(view)
    all_residues=[r for chain in view['chains'] for r in chain['residues']]
    regions=[]
    coordinate_source=case.structures[view['structure_id']][0]['source_id']
    reference=json.loads(case.payloads[case.protein['source_id']])
    annotations=[]
    for number, feature in enumerate(reference.get('features', [])):
        if feature.get('type') not in {'Active site','Binding site','Site','Domain','Region','Transmembrane','Topological domain'}:
            continue
        start=feature.get('location',{}).get('start',{})
        end=feature.get('location',{}).get('end',{})
        if start.get('modifier','EXACT') != 'EXACT' or end.get('modifier','EXACT') != 'EXACT':
            annotations.append({'id':f'feature-{number}', 'title':feature.get('type','Annotation'),
                'source_id':case.protein['source_id'], 'interval':None, 'feature_pointer':f'/features/{number}',
                'basis':'mapped_annotation', 'reason':'Source interval is not exact; no atoms are highlighted.'})
            continue
        lo, hi=start.get('value'),end.get('value')
        if type(lo) is not int or type(hi) is not int or not 1 <= lo <= hi <= len(case.protein['sequence']):
            raise ValueError('Reference feature interval differs from reference sequence')
        annotations.append({'id':f'feature-{number}', 'title':feature['type']+' · '+(
            feature.get('ligand',{}).get('name') or feature.get('description') or 'source annotation'),
            'source_id':case.protein['source_id'], 'interval':(lo,hi), 'feature_pointer':f'/features/{number}',
            'basis':'mapped_annotation', 'source_evidence':feature.get('evidences',[])})
    for assertion in case.document.get('assertions', []):
        scope=assertion['scope']
        if scope.get('frame','deposited_coordinates') != view['frame']:
            continue
        if scope.get('structure_id') and any(scope.get(k) != view[k] for k in
                ('structure_id','coordinate_sha256','model_id','assembly_id')):
            continue
        if scope.get('conformer_selection') is not None and scope['conformer_selection'] != view['conformer_selection']:
            continue
        annotations.append({'id':'assertion-'+assertion['id'], 'title':assertion['property'].replace('_',' '),
            'source_id':assertion['source_id'], 'interval':(scope['start'],scope['end']) if 'start' in scope else None,
            'basis':'mapped_annotation' if 'start' in scope else 'literature_context',
            'recorded_assertion':assertion, 'chain_copy_id':scope.get('chain_copy_id')})
    for annotation in annotations:
        interval=annotation.pop('interval')
        selected=[r for r in all_residues if interval and r['sequence_position'] is not None
                  and interval[0] <= r['sequence_position'] <= interval[1]
                  and (not annotation.get('chain_copy_id') or r['id'].split('/')[0] == annotation['chain_copy_id'])
                  and AA.get(r['comp_id']) == case.protein['sequence'][r['sequence_position']-1]]
        observed=[r for r in selected if r['observed']]
        present=sorted({r['sequence_position'] for r in observed})
        expected=list(range(interval[0],interval[1]+1)) if interval else []
        regions.append(_region(binding, annotation['id'], annotation['title'],annotation['basis'],annotation['source_id'],
            **{k:v for k,v in annotation.items() if k not in ('id','title','basis','source_id')},
            reference_interval=list(interval) if interval else None, observed_positions=present,
            missing_positions=sorted(set(expected)-set(present)),
            atom_indices=[i for r in observed for i in r['atom_indices']], residue_ids=[r['id'] for r in observed],
            mapping_state='mapped' if observed else 'unmapped', coverage={'observed':len(present),'expected':len(expected) if interval else None},
            limitations=['Annotation mapping does not establish function or functional availability.',
                        'Missing, ambiguous and chemically modified positions remain unhighlighted.']))
    expected_atoms={}
    for row in _rows(doc, '_chem_comp_atom.'):
        if row.get('type_symbol') not in ('H','D'):
            expected_atoms.setdefault(row['comp_id'],set()).add(row['atom_id'])
    polymer_positions={}
    for r in all_residues:
        if r['label_seq_id'] is not None:
            polymer_positions.setdefault(r['id'].split('/')[0],set()).add(r['label_seq_id'])
    for r in all_residues:
        expected=expected_atoms.get(r['comp_id'])
        represented={view['atoms'][i]['atom'] for i in r['atom_indices'] if view['atoms'][i]['elem'] not in ('H','D')}
        # The terminal OXT of a free component is not required for an internal peptide.
        if expected is not None:
            expected=set(expected)
            if r['label_seq_id'] is not None and r['label_seq_id']+1 in polymer_positions.get(r['id'].split('/')[0],set()):
                expected.discard('OXT')
        r['missing_heavy_atom_names']=sorted(expected-represented) if expected is not None else None
        r['atom_coverage_basis']='Deposited component atom names; not density or independent chemistry validation.'
    _connectivity(view,doc)
    protein=[a for a in view['atoms'] if a['target_protein'] and a['label_seq_id'] is not None
             and a['elem'] not in ('H','D') and a['occupancy'] > 0]
    tree=cKDTree([[a[k] for k in 'xyz'] for a in protein]) if protein else None
    for component in all_residues:
        if component['label_seq_id'] is not None or not component['observed'] or component['comp_id'] in ('HOH','DOD','WAT'):
            continue
        ligand=[view['atoms'][i] for i in component['atom_indices'] if view['atoms'][i]['elem'] not in ('H','D')
                and view['atoms'][i]['occupancy'] > 0]
        distances={}
        minimum=None
        closest=[]
        if tree is not None and ligand:
            for atom in ligand:
                xyz=[atom[k] for k in 'xyz']
                distance, near=tree.query(xyz)
                minimum=min(minimum if minimum is not None else float(distance),float(distance))
                for neighbor in tree.query_ball_point(xyz,CONTACT_CUTOFF_A):
                    target=protein[neighbor]
                    d=math.dist(xyz,[target[k] for k in 'xyz'])
                    rid=target['properties']['residue_id']
                    pair=(d,target['index'],atom['index'])
                    if rid not in distances or pair < distances[rid]:distances[rid]=pair
            closest=[{'residue_id':rid,'distance_A':pair[0],'protein_atom_index':pair[1],
                      'ligand_atom_index':pair[2]} for rid,pair in sorted(distances.items())]
        partner=[r for r in all_residues if r['id'] in distances]
        missing=component['missing_heavy_atom_names']
        regions.append(_region(binding,'contact-'+component['id'],component['comp_id']+' · observed heavy-atom contacts',
            'measured_contact',coordinate_source, component_id=component['id'],
            atom_indices=[i for r in partner for i in r['atom_indices']], ligand_atom_indices=component['atom_indices'],
            residue_ids=[r['id'] for r in partner], mapping_state='measured' if tree is not None and ligand else 'unevaluated',
            coverage={'observed':len(ligand),'expected':len(ligand)+len(missing) if missing is not None else None},
            missing_heavy_atom_names=missing, measurement={'method':'heavy_atom_min_distance_v1',
                'cutoff_A':CONTACT_CUTOFF_A,'minimum_distance_A':minimum,'contact_residue_count':len(partner),
                'closest_pairs':closest}, limitations=['Proximity does not establish bonding, affinity, importance or efficacy.',
                    'Contacts describe represented heavy atoms only; missing chemical atoms remain missing.',
                    'Component names do not classify natural substrates, antibiotics or reacted forms without annotations.']))
    view['regions']=regions
    view['region_method']={'id':METHOD,'contact_cutoff_A':CONTACT_CUTOFF_A,
                          'source_sha256':method_sources['regions_sha256'],'adapter_sha256':method_sources['adapter_sha256'],
                          'method_sha256':view['method_sha256']}
    view['recorded_states']=[a for a in case.document.get('assertions',[]) if a.get('property') in
        ('recorded_state','catalytic_state','conformation') and a.get('scope',{}).get('structure_id') == view['structure_id']
        and all(a['scope'].get(k) == view[k] for k in ('coordinate_sha256','model_id','assembly_id'))
        and a['scope'].get('conformer_selection') == view['conformer_selection']
        and a['scope'].get('frame') == view['frame']]
    view['limitations']=['Static observations do not imply chronology, transition probabilities or energies.',
                        'Predicted coordinates do not supply an observed catalytic state.',
                        'Membrane overlays require separately compatible recorded placement and sidedness evidence.']
    return view
