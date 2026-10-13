#!/usr/bin/env python3
"""Frozen case-study workflow. All engine calls are local, from installed wheels."""
from __future__ import annotations
import argparse, copy, datetime, json, os, subprocess, sys, traceback
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from adapter import (ROOT, read, write, digest, pair, ca_mask, residue_dict, selected,
    save_structure, verify_mapping, manifest, burial, load, state, horn_rmsd)

def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def config(): return read(ROOT/'CONFIG.json')
def sources(): return {x['pdb_entry_id']:x for x in read(ROOT/'SOURCES.json')['sources']}

def check_lock(name):
    lock=read(ROOT/name)
    for rel,expected in lock['files'].items():
        if digest(ROOT/rel)!=expected: raise RuntimeError(f'Frozen artifact changed: {rel}')
    return lock

def lock(name, paths):
    write(ROOT/name,{'created_utc':utc(),'files':{str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths)) if p.is_file()}})

def prepare():
    check_lock('PROTOCOL_LOCK.json')
    target=ROOT/'prepared'; target.mkdir(exist_ok=False)
    src=sources(); cfg=config(); errors={}; pairs={}; ledger={}; triage={}
    for pdb,item in src.items():
        path=ROOT/item['artifact']
        try:
            if digest(path)!=item['sha256']: raise ValueError('Source checksum mismatch')
            validation=ROOT/item['validation']['xml_artifact'] if item['validation']['identity_revision_match'] else None
            ev=manifest(path,target/'triage'/pdb,item,validation=validation)
            triage[pdb]={'status':'computed','external_validation':ev['external_validation'],'limitations':ev.get('limitations',[])}
        except Exception as e:
            triage[pdb]={'status':'failed','error':str(e)}; errors[pdb]=str(e)
    write(target/'TRIAGE.json',triage)
    for qid,spec in cfg['protomers'].items():
        try:
            if spec['pdb'] in errors: raise ValueError(errors[spec['pdb']])
            original=ROOT/src[spec['pdb']]['artifact']
            structure,record=pair(original,spec)
            pairs[qid]=structure
            path=target/qid/'pair.cif'; save_structure(path,structure)
            manifest(path,target/qid/'pair-qc',src[spec['pdb']])
            record.update(spec); record['source_sha256']=digest(original); record['derived_pair_sha256']=digest(path)
            d=MMCIF2Dict(str(original)); atom_rows={}
            for auth,label,entity,authnum,icode,labelnum,group in zip(d['_atom_site.auth_asym_id'],d['_atom_site.label_asym_id'],d['_atom_site.label_entity_id'],d['_atom_site.auth_seq_id'],d['_atom_site.pdbx_PDB_ins_code'],d['_atom_site.label_seq_id'],d['_atom_site.group_PDB']):
                if group=='ATOM': atom_rows[(auth,int(authnum),'' if icode in ('.','?') else icode)]=(label,entity,labelnum)
            record['residues']=[]
            for key,r in residue_dict(structure).items():
                original_chain=spec['ragA_auth'] if key[0]=='A' else spec['ragB_auth']
                label,entity,labelnum=atom_rows[(original_chain,key[1],key[2])]
                record['residues'].append({'role':'RagA' if key[0]=='A' else 'RagB','canonical_chain':key[0],'auth_chain':original_chain,'auth_seq_id':key[1],'insertion_code':key[2],'label_chain':label,'label_seq_id':labelnum,'entity_id':entity,'resname':r.resname})
            ledger[qid]=record
        except Exception as e: errors[qid]=str(e)
    write(target/'CHAIN_RESIDUE_LEDGER.json',ledger)
    for qid in ('REF_OPEN','REF_CLOSED','Q1','Q2'):
        if qid in errors or qid not in pairs:
            write(target/'PREPARATION.json',{'status':'primary_blocked','errors':errors}); return
    refs=[pairs['REF_OPEN'],pairs['REF_CLOSED']]; mask=ca_mask(*refs)
    baseline={ch:sum(k[0]==ch for k in mask) for ch in ('A','B')}
    write(target/'REFERENCE_MASK.json',{'keys':mask,'baseline_counts':baseline,'excluded_RagA_positions':list(range(103,115))+list(range(838,842))})
    for qid,structure in pairs.items():
        try:
            rd=residue_dict(structure); keys=[k for k in mask if k in rd and 'CA'in rd[k]]
            coverage={ch:sum(k[0]==ch for k in keys)/baseline[ch] for ch in baseline}
            if min(coverage.values())<cfg['minimum_coverage_per_protein']: raise ValueError(f'Insufficient reference-mask coverage: {coverage}')
            for ref in refs:
                rrd=residue_dict(ref)
                if any(rd[k].resname!=rrd[k].resname for k in keys): raise ValueError('Residue identity mismatch on declared author numbering')
            p=target/qid/'state-input'; paths=[]
            for label,s in [('query',structure),('open',refs[0]),('closed',refs[1])]:
                path=p/f'{label}.cif'; save_structure(path,selected(s,keys)); paths.append(path)
                manifest(path,p/f'{label}-qc',src[cfg['protomers'][qid if label=='query' else 'REF_'+label.upper()]['pdb']])
            verify_mapping(paths,keys)
            write(p/'MAPPING.json',{'keys':keys,'coverage':coverage,'baseline_counts':baseline,'identical_mask_for_both_references':True,'engine_alignment_exactly_verified':True})
            reference_set={'schema_version':'1.0','reference_set_id':'ragab-frozen-v1','qualification_scope':'ragab_case_study_unqualified','decision_rules':cfg['decision_rules'],
              'references':[{'reference_id':label,'state':alias,'structure':f'{label}.cif','provenance':{'class':'experimental','method':'ELECTRON MICROSCOPY'},'state_evidence':{'basis':f'Depositor assigned {label}-{label} assembly; compatibility alias only'}} for label,alias in [('open','active'),('closed','inactive')]]}
            write(p/'references.json',reference_set)
        except Exception as e: errors[qid]=str(e)
    # Common heavy-atom inventory for the paired burial sensitivity, defined before scoring.
    def atom_keys(s):
        return {(*k,a.name) for k,r in residue_dict(s).items() for a in r}
    common=atom_keys(pairs['Q1'])&atom_keys(pairs['Q2'])
    for qid in ('Q1','Q2'):
        save_structure(target/qid/'pair-common-atoms.cif',selected(pairs[qid],{k[:3] for k in common},ca_only=False,atom_keys=common))
        manifest(target/qid/'pair-common-atoms.cif',target/qid/'common-atoms-qc',src[cfg['protomers'][qid]['pdb']])
    write(target/'COMMON_ATOMS.json',{'keys':sorted(common),'count':len(common),'caveat':'Common-inventory truncation is a sensitivity analysis, not a repaired structure.'})
    write(target/'PREPARATION.json',{'status':'primary_blocked' if any(q in errors for q in ('REF_OPEN','REF_CLOSED','Q1','Q2')) else 'ready','errors':errors})

def controls():
    check_lock('PROTOCOL_LOCK.json')
    out=ROOT/'controls'; out.mkdir(exist_ok=False)
    cmd=[sys.executable,'-m','pytest',str(ROOT/'scripts'/'test_adapter.py'),'-q','-p','no:cacheprovider']
    p=subprocess.run(cmd,capture_output=True,text=True,cwd=ROOT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
    (out/'pytest.log').write_text(p.stdout+p.stderr)
    write(out/'RESULT.json',{'utc':utc(),'command':cmd,'exit_code':p.returncode})
    if p.returncode: raise RuntimeError('Controls failed; queries must not run')
    pth=ROOT/'prepared'/'REF_OPEN'/'state-input'
    keys=read(pth/'MAPPING.json')['keys']; arrays=verify_mapping([pth/'query.cif',pth/'open.cif',pth/'closed.cif'],keys)
    cross=state._superpose(arrays[1],arrays[2])[1]
    write(out/'REFERENCE_COMPARISON.json',{'rmsd_open_to_closed_A':cross,'independent_quaternion_rmsd_A':horn_rmsd(arrays[1],arrays[2]),'numerical_tolerance_A':1e-4,'note':'Numerical control only; not an uncertainty estimate.'})

def freeze():
    check_lock('PROTOCOL_LOCK.json')
    prep=read(ROOT/'prepared'/'PREPARATION.json')
    if prep['status']=='ready' and read(ROOT/'controls'/'RESULT.json')['exit_code']!=0: raise RuntimeError('Controls are not passing')
    paths=[]
    for name in ['prepared','controls']:
        paths.extend(p for p in (ROOT/name).rglob('*') if p.is_file())
    paths.extend(ROOT/x for x in read(ROOT/'PROTOCOL_LOCK.json')['files'])
    paths.append(ROOT/'PROTOCOL_LOCK.json')
    lock('EXECUTION_LOCK.json',paths)

def leg(kind,qid,out):
    check_lock('EXECUTION_LOCK.json')
    cfg=config(); source=sources()[cfg['protomers'][qid]['pdb']]
    inp=ROOT/'prepared'/qid; out=Path(out)
    if kind=='state':
        out.mkdir(parents=True,exist_ok=False); inp=inp/'state-input'
        m=read(inp/'MAPPING.json'); verify_mapping([inp/'query.cif',inp/'open.cif',inp/'closed.cif'],m['keys'])
        evidence=read(inp/'query-qc'/'STRUCTURE_EVIDENCE.json')
        result=state.analyze(evidence,read(inp/'references.json'),reference_base=inp,structure_path=inp/'query.cif')
        state.write_outputs(out/'engine',result)
        row=result['frame_metrics'][0]
        aliases={'active_like':'open_like','inactive_like':'closed_like','unresolved':'unresolved'}
        write(out/'summary.json',{'call':aliases[row['call']],'engine_call':row['call'],'rmsd_open_A':row['rmsd_open_A'],'rmsd_closed_A':row['rmsd_closed_A'],'signed_margin_closed_minus_open_A':row['rmsd_closed_A']-row['rmsd_open_A'],'coverage':m['coverage'],'operational_thresholds_not_validated_for_RagAB':True})
    else:
        burial(inp/('pair-common-atoms.cif' if kind=='common' else 'pair.cif'),out,source)

def invoke(kind,qid,base):
    out=base/qid/kind; out.parent.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,str(Path(__file__).resolve()),'leg',kind,qid,str(out)]
    p=subprocess.run(cmd,capture_output=True,text=True,cwd=ROOT)
    (out.parent/f'{kind}.log').write_text(p.stdout+p.stderr)
    write(out.parent/f'{kind}.invocation.json',{'command':cmd,'exit_code':p.returncode,'utc':utc()})
    return read(out/'summary.json') if p.returncode==0 else {'status':'failed','error':(p.stderr or p.stdout)[-3000:]}

def run():
    check_lock('EXECUTION_LOCK.json')
    base=ROOT/'runs'/'execution-001'; base.mkdir(parents=True,exist_ok=False)
    write(base/'START.json',{'utc':utc(),'execution_lock_sha256':digest(ROOT/'EXECUTION_LOCK.json')})
    prep=read(ROOT/'prepared'/'PREPARATION.json'); cfg=config()
    rows={f'Q{i}':{'query':f'Q{i}','status':'not_run_after_stopping'} for i in range(1,11)}
    primary_pass=False
    if prep['status']=='ready':
        # Each worker reads only frozen inputs, never another measurement's output.
        for q in ('Q1','Q2'):
            rows[q]={'query':q,'state':invoke('state',q,base),'burial':invoke('burial',q,base),'common_atoms_burial':invoke('common',q,base)}
            rows[q]['status']='failed' if rows[q]['state'].get('status')=='failed' else ('unresolved' if rows[q]['state']['call']=='unresolved' else 'resolved')
        primary_pass=rows['Q1']['state'].get('call')=='open_like' and rows['Q2']['state'].get('call')=='closed_like'
        write(base/'PRIMARY_DECISION.json',{'utc':utc(),'passed':primary_pass,'criterion':'Q1 open_like AND Q2 closed_like under frozen thresholds','rows':[rows['Q1'],rows['Q2']]})
        if primary_pass:
            for i in range(3,11):
                q=f'Q{i}'
                if q in prep['errors']:
                    rows[q]={'query':q,'status':'excluded','reason':prep['errors'][q]}; continue
                rows[q]={'query':q,'state':invoke('state',q,base),'burial':invoke('burial',q,base)}
                rows[q]['status']='failed' if rows[q]['state'].get('status')=='failed' else ('unresolved' if rows[q]['state']['call']=='unresolved' else 'resolved')
    else:
        for q in ('Q1','Q2'): rows[q]={'query':q,'status':'failed','reason':'primary preparation blocked','preparation_errors':prep['errors']}
    write(base/'CALLS_BEFORE_SIGNATURES.json',{'utc':utc(),'primary_pass':primary_pass,'rows':list(rows.values())})
    assemblies=[]
    for spec in cfg['assemblies']:
        members=spec['query_members']; labels=[rows[q].get('state',{}).get('call') for q in members]
        state_name='unresolved'
        if primary_pass and all(x in ('open_like','closed_like') for x in labels):
            state_name='-'.join(sorted([x.replace('_like','') for x in labels],reverse=True))
        assemblies.append({**spec,'composed_state':state_name,'symmetry_copies_are_not_independent':len(set(members))<len(members)})
    signatures=[]
    if primary_pass:
        for i in range(3,11):
            q=f'Q{i}'; spec=cfg['protomers'][q]; model=load(ROOT/sources()[spec['pdb']]['artifact'])
            pos={int(r.id[1]) for r in model[spec['ragA_auth']] if r.id[0]==' '}
            signatures.append({'query':q,'positions_103_114_present':sorted(pos&set(range(103,115))),'positions_838_841_present':sorted(pos&set(range(838,842))),'interpretation':'model-derived consistency check; counts/spans already exposed before execution'})
    write(base/'ASSEMBLIES_AND_SIGNATURES.json',{'utc':utc(),'assemblies':assemblies,'signatures':signatures,'calls_sha256':digest(base/'CALLS_BEFORE_SIGNATURES.json')})
    write(base/'END.json',{'utc':utc(),'primary_pass':primary_pass})
    print(json.dumps({'primary_pass':primary_pass,'queries':[{k:v for k,v in row.items() if k not in ('burial','common_atoms_burial')} for row in rows.values()]},indent=2))

def verify():
    check_lock('PROTOCOL_LOCK.json'); check_lock('EXECUTION_LOCK.json')
    if (ROOT/'OUTPUT_MANIFEST.json').exists(): check_lock('OUTPUT_MANIFEST.json')
    data=read(ROOT/'runs'/'execution-001'/'CALLS_BEFORE_SIGNATURES.json')
    assert [x['query'] for x in data['rows']]==[f'Q{i}' for i in range(1,11)]
    if not data['primary_pass']:
        assert all(x['status']=='not_run_after_stopping' for x in data['rows'][2:])
        assert not any((ROOT/'runs'/'execution-001'/f'Q{i}'/'state').exists() for i in range(3,11))
    print('All frozen hashes and all ten query records verified; stopping rule enforced.')

if __name__=='__main__':
    command=sys.argv[1]
    if command=='leg': leg(*sys.argv[2:])
    else: {'prepare':prepare,'controls':controls,'freeze':freeze,'run':run,'verify':verify}[command]()
