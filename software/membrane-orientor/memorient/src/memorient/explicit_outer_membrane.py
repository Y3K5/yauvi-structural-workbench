"""Source-bound finite E. coli K-12 atomic lipid patch; no MD execution.

Tiles an author-deposited AA membrane snapshot using PSF covalent components,
unwraps covalent molecules, and removes whole lipids colliding with OmpF/FhuA.
This is a coordinate assembly, never a complete or equilibrated MD system.
"""
from __future__ import annotations
from collections import Counter, defaultdict, deque
import hashlib
import json
from pathlib import Path
import math
import numpy as np
from scipy.spatial import cKDTree
from .starter_catalogue import protein_record

HOME=Path(__file__).with_name('ecoli_om')
EXPECTED_MD5={'K12-AA-rep1.pdb':'71d5cc6c3c316f38a39ebeec1d77cb2a','K12-AA-rep1.psf':'9fae2cc690e5f37060732c1f850af807'}


def load_donor(home=HOME):
    source=Path(home)/'sources'
    for name,expected in EXPECTED_MD5.items():
        if hashlib.md5((source/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('Donor source identity mismatch: '+name)
    lines=(source/'K12-AA-rep1.psf').read_text().splitlines()
    start=next(i for i,s in enumerate(lines) if '!NATOM' in s)
    count=int(lines[start].split()[0])
    atoms=[]
    for line in lines[start+1:start+1+count]:
        f=line.split();mass=float(f[7])
        element='H' if mass<2.1 else 'C' if mass<13 else 'N' if mass<15 else 'O' if mass<17 else 'P' if 30<mass<32 else 'S' if 32<=mass<33 else 'MG' if 24<mass<25 else 'CA' if 39<mass<41 else 'K' if 38<mass<40 else 'CL' if 35<mass<36 else 'X'
        element={'POT':'K','CLA':'CL','CAL':'CA','MG':'MG'}.get(f[3],element)
        atoms.append({'source_index':int(f[0])-1,'segment':f[1],'resid':f[2], 'resname':f[3], 'name':f[4], 'type':f[5], 'charge':float(f[6]),'mass':mass,'element':element})
    p=(source/'K12-AA-rep1.pdb').read_text().splitlines()
    box=np.array([float(p[0][a:b]) for a,b in ((6,15),(15,24),(24,33))])
    angles=[float(p[0][a:b]) for a,b in ((33,40),(40,47),(47,54))]
    if angles!=[90.,90.,90.]:raise ValueError('Only the pinned orthorhombic donor is supported')
    plines=[s for s in p if s[:6] in ('ATOM  ','HETATM')]
    if len(plines)!=count:raise ValueError('Donor PDB/PSF atom counts differ')
    # Source PDB truncates long residue/segment names; PSF is authoritative.
    for i,line in enumerate(plines):
        if line[12:16].strip()!=atoms[i]['name'][:4]:raise ValueError('Donor PDB/PSF atom order mismatch')
    xyz=np.array([[float(s[a:b]) for a,b in ((30,38),(38,46),(46,54))] for s in plines])
    start=next(i for i,s in enumerate(lines) if '!NBOND' in s);nbonds=int(lines[start].split()[0]);numbers=[]
    for line in lines[start+1:]:
        if len(numbers)>=nbonds*2:break
        numbers.extend(int(x)-1 for x in line.split())
    bonds=np.asarray(numbers[:nbonds*2],dtype=int).reshape(-1,2)
    adjacency=[[] for _ in atoms]
    for a,b in bonds:adjacency[a].append(b);adjacency[b].append(a)
    allowed={i for i,a in enumerate(atoms) if a['segment'].startswith('LPSA') or a['segment']=='MEMB'}
    seen=set();molecules=[]
    for root in sorted(allowed):
        if root in seen:continue
        indices=[];queue=[root];seen.add(root)
        while queue:
            i=queue.pop();indices.append(i)
            for j in adjacency[i]:
                if j not in allowed:raise ValueError('Unexpected lipid bond to excluded source atom')
                if j not in seen:seen.add(j);queue.append(j)
        indices.sort();idx=np.array(indices);local={old:new for new,old in enumerate(indices)}
        unwrapped={root:xyz[root].copy()};queue=deque([root])
        while queue:
            i=queue.popleft()
            for j in adjacency[i]:
                delta=xyz[j]-xyz[i];delta-=box*np.round(delta/box)
                candidate=unwrapped[i]+delta
                if j in unwrapped:
                    if np.linalg.norm(candidate-unwrapped[j])>1e-5:raise ValueError('Inconsistent periodic covalent unwrap')
                else:unwrapped[j]=candidate;queue.append(j)
        coords=np.array([unwrapped[i] for i in indices])
        # Membrane crosses the z box boundary in this source snapshot.
        center=coords.mean(axis=0);target=np.array([0.,0.,box[2]/2])
        coords-=box*np.round((center-target)/box);coords[:,2]-=box[2]/2
        group_atoms=[atoms[i] for i in indices]
        names={a['resname'] for a in group_atoms}
        kind='LPS' if 'ECLIPA' in names else group_atoms[0]['resname']
        if kind not in {'LPS','PVPE','PVPG','PVCL2'}:raise ValueError('Unclassified donor molecule')
        group_bonds=[(local[i],local[j]) for i in indices for j in adjacency[i] if i<j]
        molecules.append({'kind':kind,'atoms':group_atoms,'xyz':coords,'bonds':group_bonds,'source_root':root})
    # Midpoint of the two leaflet phosphorus distributions; not core thickness.
    upper=[m['xyz'][i,2] for m in molecules if m['kind']=='LPS' for i,a in enumerate(m['atoms']) if a['element']=='P' and a['resname']=='ECLIPA']
    lower=[m['xyz'][i,2] for m in molecules if m['kind']!='LPS' for i,a in enumerate(m['atoms']) if a['element']=='P']
    center_z=(np.mean(upper)+np.mean(lower))/2
    for m in molecules:m['xyz'][:,2]-=center_z
    return atoms,xyz,box,molecules,float(center_z),float(np.mean(upper)-np.mean(lower))


def protein_part(identifier,translation):
    record=protein_record(identifier)
    atoms=[];coords=[]
    for line in record['pdb'].splitlines():
        if line[:6]!='ATOM  ':continue
        atoms.append({'source_index':int(line[6:11])-1,'segment':identifier+line[21], 'resid':line[22:27].strip(), 'resname':line[17:20].strip(), 'name':line[12:16].strip(), 'element':line[76:78].strip() or line[12:16].strip()[0], 'charge':None})
        coords.append([float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))])
    original=np.array(coords);coords=original+translation
    return {'kind':identifier,'atoms':atoms,'xyz':coords,'bonds':[], 'source_sha256':record['coordinate_source_sha256'],'view_sha256':record['view_pdb_sha256'],'translation':list(translation),'removed_heteroatoms':record['view_atoms']-len(atoms)}


def view_pdb(molecules,chain):
    """Heavy-atom display derivative, explicit lipid CONECT, independent serials."""
    rows=[];bondrows=[];serial=0;residue=0;mapping={};groupmap=[]
    for molecule_id,m in enumerate(molecules):
        atommap={};old_residue=None;group_start=serial+1
        for i,(a,p) in enumerate(zip(m['atoms'],m['xyz'])):
            if a['element']=='H':continue
            key=(a['segment'],a['resid'])
            if key!=old_residue:residue+=1;old_residue=key
            serial+=1;atommap[i]=serial
            if serial>99999 or residue>9999:raise ValueError('Display PDB field limit exceeded')
            name=a['name'][:4];resname=a['resname'][:3]
            row_chain=a['segment'][-1] if m['kind'] in ('2OMF','2FCP') else chain
            row_residue=int(a['resid']) if m['kind'] in ('2OMF','2FCP') else residue
            rows.append(f"{'ATOM' if m['kind'] in ('2OMF','2FCP') else 'HETATM':<6}{serial:5d} {name:>4} {resname:>3} {row_chain}{row_residue:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {a['element']:>2}")
        for a,b in m['bonds']:
            if a in atommap and b in atommap:bondrows.append(f'CONECT{atommap[a]:5d}{atommap[b]:5d}')
        groupmap.append({'molecule_index':molecule_id,'kind':m['kind'],'serial_start':group_start,'serial_end':serial,'source_root_index':m.get('source_root'),'tile':m.get('tile'),'source_segment':m['atoms'][0]['segment'],'source_residue_names':sorted(set(a['resname'] for a in m['atoms']))})
    return '\n'.join(rows+bondrows+['END'])+'\n',groupmap,serial


def write_cif(path,molecules):
    rows=['data_ecoli_k12_outer_membrane','#','_audit_creation_method \'MemOrient source-derived finite coordinate assembly; unrelaxed\'','#','loop_','_atom_site.group_PDB','_atom_site.id','_atom_site.type_symbol','_atom_site.label_atom_id','_atom_site.label_alt_id','_atom_site.label_comp_id','_atom_site.label_asym_id','_atom_site.label_entity_id','_atom_site.label_seq_id','_atom_site.pdbx_PDB_ins_code','_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z','_atom_site.occupancy','_atom_site.B_iso_or_equiv','_atom_site.auth_atom_id','_atom_site.auth_seq_id','_atom_site.auth_asym_id','_atom_site.pdbx_PDB_model_num']
    serial=0;mapping=[]
    for mid,m in enumerate(molecules,1):
        ids=[];names_seen=Counter()
        for a,p in zip(m['atoms'],m['xyz']):
            serial+=1;ids.append(serial)
            name_key=(a['segment'],a['resid'],a['name']);names_seen[name_key]+=1
            label_name=a['name'] if names_seen[name_key]==1 else a['name']+'__'+str(names_seen[name_key])
            # Original full CHARMM residue/atom names are retained in mmCIF.
            rows.append(f"{'ATOM' if m['kind'] in ('2OMF','2FCP') else 'HETATM'} {serial} {a['element']} {json.dumps(label_name)} . {a['resname']} M{mid}_{a['segment']} {mid} . ? {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} 1.0 0.0 {json.dumps(a['name'])} {json.dumps(a['resid'])} M{mid}_{a['segment']} 1")
        mapping.append(ids)
    rows.append('#');Path(path).write_text('\n'.join(rows)+'\n')
    return mapping


def build_patch(*,tiles=(3,2),protein_positions=((-48.,0.,0.),(48.,0.,0.)),clash_distance=2.6):
    if tiles!=(3,2) or clash_distance!=2.6:raise ValueError('This version supports only its recorded 3x2 finite-patch recipe')
    atoms,source_xyz,box,donor,center_z,p_spacing=load_donor()
    proteins=[protein_part(id,np.array(pos)) for id,pos in zip(('2OMF','2FCP'),protein_positions)]
    pxyz=np.concatenate([p['xyz'] for p in proteins]);tree=cKDTree(pxyz)
    min_protein=float(cKDTree(proteins[0]['xyz']).query(proteins[1]['xyz'])[0].min())
    if min_protein<clash_distance:raise ValueError('OmpF/FhuA placement collides')
    accepted=[];removed=Counter();candidates=Counter();min_contact=float('inf')
    for ix in range(tiles[0]):
        for iy in range(tiles[1]):
            shift=np.array([(ix-(tiles[0]-1)/2)*box[0],(iy-(tiles[1]-1)/2)*box[1],0.])
            for m in donor:
                coords=m['xyz']+shift;heavy=np.array([a['element']!='H' for a in m['atoms']])
                distance=float(tree.query(coords[heavy])[0].min());candidates[m['kind']]+=1
                if distance<clash_distance:removed[m['kind']]+=1;continue
                min_contact=min(min_contact,distance)
                accepted.append({**m,'xyz':coords,'tile':[ix,iy]})
    # Include source divalent ions only near retained LPS and away from proteins.
    lps_xyz=np.concatenate([m['xyz'][[a['element']!='H' for a in m['atoms']]] for m in accepted if m['kind']=='LPS'])
    lpstree=cKDTree(lps_xyz);ions=[]
    for ix in range(tiles[0]):
        for iy in range(tiles[1]):
            shift=np.array([(ix-1)*box[0],(iy-.5)*box[1],0.])
            for i,a in enumerate(atoms):
                if a['resname'] not in {'MG','CAL'}:continue
                coord=source_xyz[i].copy();coord[2]-=box[2]/2;coord-=box*np.round(coord/box);coord[2]-=center_z;coord+=shift
                if lpstree.query(coord)[0]>6.0 or tree.query(coord)[0]<clash_distance:continue
                ions.append({'kind':a['resname'],'atoms':[a],'xyz':coord[None,:],'bonds':[]})
    all_molecules=proteins+accepted+ions
    # Only severe new lipid-lipid overlaps are checked, not force-field energies.
    heavyxyz=np.concatenate([m['xyz'][[a['element']!='H' for a in m['atoms']]] for m in accepted])
    molids=np.concatenate([np.full(sum(a['element']!='H' for a in m['atoms']),i) for i,m in enumerate(accepted)])
    pairs=cKDTree(heavyxyz).query_pairs(1.5,output_type='ndarray');severe=sum(molids[a]!=molids[b] for a,b in pairs)
    if severe:raise ValueError(f'{severe} severe intermolecular lipid overlaps after tiling')
    ion_xyz=np.concatenate([m['xyz'] for m in ions]) if ions else np.empty((0,3))
    ion_lipid_severe=int(sum(len(x) for x in cKDTree(heavyxyz).query_ball_point(ion_xyz,1.5))) if len(ion_xyz) else 0
    ion_ion_severe=len(cKDTree(ion_xyz).query_pairs(1.5)) if len(ion_xyz) else 0
    if ion_lipid_severe or ion_ion_severe:raise ValueError('Severe ion overlap in assembled patch')
    groups={'ompf':[proteins[0]],'fhua':[proteins[1]],'outer_lps':[m for m in accepted if m['kind']=='LPS'],'inner_pl':[m for m in accepted if m['kind']!='LPS'],'divalent_ions':ions}
    colors={'ompf':'#70c9fa','fhua':'#ff8fa7','outer_lps':'#72d9ae','inner_pl':'#e6bb68','divalent_ions':'#c8a5ff'}
    parts=[]
    for chain,(kind,molecules) in zip('ABCDI',groups.items()):
        pdb,groups_map,count=view_pdb(molecules,chain)
        parts.append({'id':kind,'pdb':pdb,'color':colors[kind],'view_atoms':count,'molecules':groups_map,'full_atom_count':sum(len(m['atoms']) for m in molecules)})
    counts=Counter(m['kind'] for m in accepted);ion_counts=Counter(m['kind'] for m in ions)
    recipe={'schema_version':'organism_orient.explicit_om.v1','organism':'Escherichia coli K-12','state':'explicit_atomic_coordinate_assembly_unrelaxed','representation':'Atomistic source lipids (hydrogens retained in mmCIF); heavy-atom experimental proteins. Viewer omits hydrogens.',
            'donor':'Zenodo 18009953 K12-AA-rep1 PDB + PSF','donor_doi':'https://doi.org/10.5281/zenodo.18009953','paper':'https://doi.org/10.1021/acs.jpcb.6c00050','source_license':'CC-BY-4.0',
            'source_composition':dict(Counter(m['kind'] for m in donor)), 'candidate_lipids':dict(candidates),'removed_whole_lipids':dict(removed),'retained_lipids':dict(counts),'retained_divalent_ions':dict(ion_counts),
            'composition_discrepancy':'Paper describes 75% POPE /25% POPG; deposited coordinates + PSF contain 90 PVPE /5 PVPG /5 PVCL2 per donor patch. Deposited identities retained; not relabelled as POPE/POPG.',
            'donor_box_angstrom':box.tolist(),'tiles':list(tiles),'nominal_patch_xy_angstrom':(box[:2]*tiles).tolist(),'pbc_handling':'PSF covalent graph unwrap, molecule-wise recenter, whole-molecule 3x2 tile; final assembly is a finite nonperiodic patch',
            'donor_phosphorus_leaflet_spacing_angstrom':p_spacing,'source_membrane_center_offset_angstrom':center_z,
            'protein_placements':[{'id':p['kind'],'translation_angstrom':p['translation'],'source_sha256':p['source_sha256'],'view_sha256':p['view_sha256'],'protein_atoms':len(p['atoms']),'excluded_heteroatoms':p['removed_heteroatoms']} for p in proteins],
            'checks':{'protein_protein_min_heavy_distance_angstrom':min_protein,'retained_protein_lipid_min_heavy_distance_angstrom':min_contact,'protein_lipid_removal_cutoff_angstrom':clash_distance,'intermolecular_lipid_pairs_below_1_5_angstrom':int(severe),'ion_lipid_pairs_below_1_5_angstrom':ion_lipid_severe,'ion_ion_pairs_below_1_5_angstrom':ion_ion_severe},
            'lipid_partial_charge_e':round(sum(float(a['charge']) for m in accepted for a in m['atoms']),6),'retained_ion_charge_e':round(sum(float(a['charge']) for m in ions for a in m['atoms']),6),
            'source_atom_name_note':'Donor PSF repeats some hydrogen atom names within a residue. Export label_atom_id adds occurrence suffixes; original names remain in auth_atom_id and ATOM_MAP.tsv; atom-order identity and coordinates are preserved.',
            'protein_charge':'not assigned; no complete force-field protein topology','full_atoms':sum(len(m['atoms']) for m in all_molecules),'view_atoms':sum(p['view_atoms'] for p in parts),
            'limitations':['New assembly has not been minimized or equilibrated; not MD ready','No native OmpF/FhuA relative placement or copy-number inference','Whole-lipid exclusion can leave gaps requiring repacking/relaxation','Water and monovalent ions omitted; charge neutralization is not established','Source Ca/Mg kept only within 6 A of retained LPS; not a new ion-equilibration calculation','No periodic combined protein-system topology or force-field parameters supplied','Heavy-atom overlap checks do not assess energies, all vdW clashes or ring penetration','This is a modeled molecular patch, not a measured native membrane or whole cell'],
            'parts':parts}
    return recipe,all_molecules


def save_patch(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing to overwrite an existing model directory')
    scene,molecules=build_patch()
    destination.mkdir(parents=True)
    mapping=write_cif(destination/'ECOLI_K12_OM.cif',molecules)
    with (destination/'LIPID_BONDS.tsv').open('w') as stream:
        stream.write('atom_id_1\tatom_id_2\n')
        for m,ids in zip(molecules,mapping):
            for a,b in m['bonds']:stream.write(f'{ids[a]}\t{ids[b]}\n')
    atom_id=0;mapping_rows=['atom_id\tmolecule_id\tkind\tsource_atom_index_1based\tsource_segment\tsource_residue_id\tsource_residue_name\tsource_atom_name\ttile\n']
    for mid,m in enumerate(molecules,1):
        for a in m['atoms']:
            atom_id+=1
            mapping_rows.append(f"{atom_id}\t{mid}\t{m['kind']}\t{a['source_index']+1}\t{a['segment']}\t{a['resid']}\t{a['resname']}\t{a['name']}\t{m.get('tile',[])}\n")
    (destination/'ATOM_MAP.tsv').write_text(''.join(mapping_rows))
    scene['atom_map_file']='ATOM_MAP.tsv'
    scene['atom_map_sha256']=hashlib.sha256((destination/'ATOM_MAP.tsv').read_bytes()).hexdigest()
    scene['view_part_sha256']={p['id']:hashlib.sha256(p['pdb'].encode()).hexdigest() for p in scene['parts']}
    scene['coordinate_file']='ECOLI_K12_OM.cif';scene['lipid_bond_file']='LIPID_BONDS.tsv'
    scene['coordinate_sha256']=hashlib.sha256((destination/'ECOLI_K12_OM.cif').read_bytes()).hexdigest()
    scene['lipid_bonds_sha256']=hashlib.sha256((destination/'LIPID_BONDS.tsv').read_bytes()).hexdigest()
    scene['source_sha256']={name:hashlib.sha256((HOME/'sources'/name).read_bytes()).hexdigest() for name in ('zenodo-18009953.json',*EXPECTED_MD5)}
    (destination/'scene.json').write_text(json.dumps(scene,separators=(',',':'))+'\n')
    receipt={k:v for k,v in scene.items() if k!='parts'}
    (destination/'BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt
