"""Case-local preparation only; installed engine sources remain unchanged."""
from __future__ import annotations
import copy, hashlib, json, shutil, subprocess
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, MMCIFIO, Structure, Model, Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from structqc import core as qc
from state_atlas import core as state
from assembly_context import core as assembly

ROOT = Path(__file__).resolve().parents[1]
NUMERICAL_TOLERANCE_A = 1e-4

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text())

def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists(): raise FileExistsError(f'Refusing to overwrite {path}')
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')

def load(path):
    structure = MMCIFParser(QUIET=True).get_structure('ragab', str(path))
    if len(structure) != 1: raise ValueError('Expected exactly one deposited model')
    return structure[0]

def protein_residues(chain):
    result = {}
    for r in chain:
        if r.id[0] != ' ': continue
        if r.resname not in state.AA3: raise ValueError(f'Unrecognized protein residue: {r.resname}')
        key = (int(r.id[1]), r.id[2].strip())
        if key in result: raise ValueError('Duplicate author residue key')
        result[key] = r
    return result

def pair(source, spec):
    """Use frozen chain choices; reject rather than optimize a wrong choice."""
    model = load(source)
    a, b = spec['ragA_auth'], spec['ragB_auth']
    metadata = MMCIF2Dict(str(source))
    entity_desc = dict(zip(metadata['_entity.id'], metadata['_entity.pdbx_description']))
    identities = {}
    label_to_auth = {}
    for auth,label,entity,group in zip(metadata['_atom_site.auth_asym_id'], metadata['_atom_site.label_asym_id'], metadata['_atom_site.label_entity_id'], metadata['_atom_site.group_PDB']):
        if group == 'ATOM':
            identities.setdefault(auth,set()).add(entity_desc[entity]); label_to_auth[label]=auth
    if identities.get(a) != {'RagA protein'} or identities.get(b) != {'Lipoprotein RagB'}:
        raise ValueError('Frozen pair does not have unique RagA/RagB polymer identities')
    recipes = list(zip(metadata['_pdbx_struct_assembly_gen.assembly_id'], metadata['_pdbx_struct_assembly_gen.oper_expression'],metadata['_pdbx_struct_assembly_gen.asym_id_list']))
    matching = [r for r in recipes if r[0]==spec['assembly_id'] and {a,b} <= {label_to_auth.get(x) for x in r[2].split(',')}]
    if len(matching)!=1: raise ValueError('Pair is not jointly in exactly one declared assembly recipe')
    anchor = model[b][20]['CA'].coord.astype(float)
    distances = {}
    for candidate in spec['candidate_ragA_auth']:
        coords=np.array([r['CA'].coord for r in protein_residues(model[candidate]).values() if 'CA' in r],dtype=float)
        distances[candidate]=float(np.linalg.norm(coords-anchor,axis=1).min())
    ordered=sorted(distances,key=distances.get)
    if ordered[0]!=a or (len(ordered)>1 and distances[ordered[1]]-distances[a]<1e-3):
        raise ValueError('Frozen pairing contradicts unique N-terminal adjacency')
    s=Structure.Structure('protein_pair'); m=Model.Model(0); s.add(m)
    altloc=[]
    for canonical,original in [('A',a),('B',b)]:
        ch=Chain.Chain(canonical); m.add(ch)
        for key,r in sorted(protein_residues(model[original]).items()):
            nr=r.copy(); nr.detach_parent()
            for atom in list(nr): nr.detach_child(atom.id)
            by_name={}
            for atom in r.get_unpacked_list(): by_name.setdefault(atom.name,[]).append(atom)
            for name,choices in sorted(by_name.items()):
                chosen=sorted(choices,key=lambda x: (-(x.occupancy if x.occupancy is not None else 0.0),x.altloc.strip()))[0]
                if len(choices)>1: altloc.append({'auth_chain':original,'residue':list(key),'atom':name,'selected':chosen.altloc,'options':[x.altloc for x in choices]})
                if chosen.element in ('H','D'): continue
                atom=chosen.copy(); atom.detach_parent(); atom.set_altloc(' '); atom.disordered_flag=0
                nr.add(atom)
            nr.disordered=0; ch.add(nr)
    return s, {'pairing_anchor_distances_A':distances,'recipe':matching[0],'altloc_choices':altloc,'canonical_chain_roles':{'A':'RagA','B':'RagB'},'original_chain_ids':{'A':a,'B':b}}

def save_structure(path, structure):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists(): raise FileExistsError(path)
    writer=MMCIFIO(); writer.set_structure(structure); writer.save(str(path))

def residue_dict(structure):
    model=next(structure.get_models()) if structure.level=='S' else structure
    return {(ch.id,*key):r for ch in model for key,r in protein_residues(ch).items()}

def ca_mask(open_pair, closed_pair):
    a,b=residue_dict(open_pair),residue_dict(closed_pair)
    return sorted(k for k in a.keys()&b.keys() if 'CA'in a[k] and 'CA'in b[k]
                  and not(k[0]=='A' and (103<=k[1]<=114 or 838<=k[1]<=841)))

def selected(structure, keys, ca_only=True, atom_keys=None):
    wanted=set(keys); result=structure.copy()
    model=next(result.get_models())
    for ch in model:
        for r in list(ch):
            key=(ch.id,int(r.id[1]),r.id[2].strip())
            if key not in wanted: ch.detach_child(r.id); continue
            for atom in list(r):
                if (ca_only and atom.name!='CA') or (atom_keys is not None and (*key,atom.name) not in atom_keys): r.detach_child(atom.id)
            if len(r)==0: ch.detach_child(r.id)
    return result

def verify_mapping(paths, keys):
    expected=[tuple(k) for k in keys]; arrays=[]; sequences=[]
    for p in paths:
        frames,actual,sequence=state._bio_models(Path(p))
        if actual != expected: raise ValueError('Missing or reordered mapped coordinates')
        arrays.append(frames[0]); sequences.append(sequence)
    for seq in sequences[1:]:
        q,r=state._sequence_indices(sequences[0],seq)
        expected_indices=np.arange(len(expected))
        if sequences[0]!=seq or not np.array_equal(q,expected_indices) or not np.array_equal(r,expected_indices):
            raise ValueError('Engine alignment differs from exact chain-aware correspondence')
    return arrays

def manifest(path,out,source,chain=None,validation=None):
    result=qc.analyze(path,subject_id=Path(path).stem,chain=chain,
       provenance={'class':'experimental','method':source['method'],'source_id':source['pdb_entry_id'],'resolution_angstrom':source['resolution_A']},
       validation_report=qc.read_validation_report(validation) if validation else None)
    qc.write_outputs(out,result)
    return result

def ensure_freesasa():
    exe=shutil.which('freesasa')
    if not exe: raise RuntimeError('Pinned FreeSASA required; fallback forbidden')
    expected=read(ROOT/'environment'/'ENVIRONMENT.json')['freesasa']
    if digest(exe)!=expected['sha256']: raise RuntimeError('FreeSASA executable changed')
    return exe

def burial(pair_path, out, source):
    ensure_freesasa(); out=Path(out); out.mkdir(parents=True,exist_ok=False)
    components={}
    structure=load(pair_path)
    for ch in ('A','B'):
        isolated=Structure.Structure('isolated'); model=Model.Model(0); isolated.add(model)
        chain=structure[ch].copy(); chain.detach_parent(); model.add(chain)
        isolated_path=out/f'isolated_{ch}.cif'; save_structure(isolated_path,isolated)
        evidence=manifest(isolated_path,out/f'qc_{ch}',source)
        result=assembly.analyze(evidence,isolated_path,pair_path,subject_chain=ch,relationship='exact_protein',expected_chains=['A','B'])
        if not result['methods']['sasa'].startswith('freesasa_'): raise RuntimeError('Unexpected SASA backend')
        assembly.write_outputs(out/f'engine_{ch}',result)
        components[ch]=result['surface']['buried_sasa_A2']
    result={'ragA_burial_A2':components['A'],'ragB_burial_A2':components['B'],
            'interface_area_A2':(components['A']+components['B'])/2,
            'convention':'half sum of the two protein burial values; pair-only modeled heavy atoms'}
    write(out/'summary.json',result); return result

def horn_rmsd(moving,fixed):
    """Independent quaternion eigenvalue solution; does not call engine SVD."""
    x=np.asarray(moving,dtype=float); y=np.asarray(fixed,dtype=float)
    x=x-x.mean(0); y=y-y.mean(0); M=x.T@y
    xx,xy,xz,yx,yy,yz,zx,zy,zz=M.ravel()
    K=np.array([[xx+yy+zz,yz-zy,zx-xz,xy-yx],
      [yz-zy,xx-yy-zz,xy+yx,zx+xz],
      [zx-xz,xy+yx,-xx+yy-zz,yz+zy],
      [xy-yx,zx+xz,yz+zy,-xx-yy+zz]])
    e=float(np.linalg.eigvalsh(K)[-1])
    return float(np.sqrt(max(0.0,(np.sum(x*x)+np.sum(y*y)-2*e)/len(x))))
