"""Dedicated MemOrient adapters over the existing custody and job machinery."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

TYPES={'membrane_surface_placement','membrane_preparation','membrane_result_import'}
CLAIM='Experimental geometry or restrained minimization only; no native exposure, equilibration or biological validation.'


def definitions():
    def role(name,label,ext):return {'role':name,'label':label,'required':True,'multiple':False,'extensions':ext}
    common={'readiness':'experimental','module_ids':['membrane_orientation'],'claim_ceiling':CLAIM,
            'use_when':'Explore geometric membranes or a separately reviewed molecular system.',
            'measures':'Declared placement, source identity and bounded numerical acceptance.',
            'receives':['Portable membrane result','Exact source records','Processing stage and limitations'],
            'non_claim':CLAIM}
    return [{**common,'analysis_type':'membrane_surface_placement','title':'Organism Orient · surface placement',
        'question':'How does this receptor sit on the declared membrane geometry?',
        'inputs':[role('structure','Local single-model receptor PDB',['.pdb','.ent']),
                  role('placement_spec','Saved surface and placement settings',['.json']),
                  {'role':'validation_report','label':'Optional coordinate validation','required':False,'multiple':False,'extensions':['.json','.xml']}],'parameters':[]},
        {**common,'analysis_type':'membrane_preparation','title':'Organism Orient · prepare and minimize',
        'question':'Can these reviewed molecular inputs be hydrated, neutralized and minimized?',
        'inputs':[role('preparation_profile','Reviewed CHARMM36 planar profile',['.json']),
                  role('molecular_topology','Complete standalone topology',['.top']),
                  role('molecular_coordinates','Checksum-matched solute coordinates',['.gro'])],
        'parameters':[{'name':'threads','label':'CPU threads (1–16)','type':'integer','default':4}]},
        {**common,'analysis_type':'membrane_result_import','title':'Organism Orient · import result',
        'question':'Which processing stage and source identities does this local result preserve?',
        'inputs':[role('membrane_bundle','Portable local result ZIP',['.zip'])],'parameters':[]}]


def runtime():
    python=os.environ.get('YAUVI_MEMORIENT_PYTHON')
    gmx=os.environ.get('YAUVI_MEMORIENT_GMX','gmx')
    if not python or not shutil.which(python) or not shutil.which(gmx):raise ValueError('Preparation held: configure the optional preparation Python and GROMACS; nothing is installed automatically')
    probe=subprocess.run([python,'-c','import numpy, scipy; import memorient.md_preparation; print("ready")'],
                         capture_output=True,text=True,timeout=15)
    if probe.returncode:raise ValueError('Preparation held: selected Python cannot import the MemOrient preparation runtime')
    return python,gmx


def readiness():
    from .store import _installed_package_root, _tree_sha
    package=_installed_package_root('memorient')
    digest=_tree_sha(package) if package else None
    rows=[]
    for definition in definitions():
        detail=CLAIM;state='experimental'
        if not package:state='blocked_missing_runtime';detail='MemOrient is unavailable'
        rows.append({'contract_id':'tool_readiness_record','analysis_type':definition['analysis_type'],
            'title':definition['title'],'state':state,'package_source_sha256':digest,
            'required_runtimes':{},'optional_runtimes':{},'scientific_scopes':[],
            'external_benchmark':{'state':'not_evaluated'},'limitations':[detail]})
    return rows


def bind_profile(profile,topology,coordinates,out):
    from memorient.md_preparation import load_profile,sha256
    doc=json.loads(Path(profile).read_text());out=Path(out)
    for key,path in [('topology',topology),('coordinates',coordinates)]:
        if sha256(path)!=doc['inputs'][key]['sha256']:raise ValueError(key+' does not match reviewed profile')
        target=out/('input.top' if key=='topology' else 'input.gro')
        shutil.copyfile(path,target)
        doc['inputs'][key]['path']=target.name
    target=out/'BOUND_PROFILE.json';target.write_text(json.dumps(doc,indent=2)+'\n')
    load_profile(target)
    return target


def unpack(source,destination):
    try:return _unpack(source,destination)
    except (zipfile.BadZipFile,UnicodeDecodeError,NotImplementedError) as exc:
        raise ValueError('Invalid or unsupported membrane result archive') from exc


def _unpack(source,destination):
    from memorient.run_bundle import artifact_path,load_manifest
    destination=Path(destination)
    with zipfile.ZipFile(source) as archive:
        rows=archive.infolist()
        names=[r.filename for r in rows]
        if len(rows)>64 or len(set(names))!=len(names) or sum(r.file_size for r in rows)>512*1024*1024:
            raise ValueError('Result archive exceeds 64 files / 512 MiB or has duplicate entries')
        if 'MEMBRANE_RUN.json' not in names:raise ValueError('Missing membrane result manifest')
        doc=json.loads(archive.read('MEMBRANE_RUN.json'))
        allowed={'MEMBRANE_RUN.json'}|{r['path'] for r in doc['artifacts'].values()}
        if set(names)!=allowed:raise ValueError('Archive must contain exactly the declared artifacts and manifest')
        for row in rows:
            path=artifact_path(destination,row.filename)
            if (row.external_attr>>16)&0o170000==0o120000:raise ValueError('Archive shortcuts are forbidden')
            path.parent.mkdir(parents=True,exist_ok=True)
            with archive.open(row) as stream,path.open('xb') as out:shutil.copyfileobj(stream,out)
    return load_manifest(destination/'MEMBRANE_RUN.json')


def preflight(store,manifest,inputs):
    kind=manifest['analysis_type']
    with tempfile.TemporaryDirectory(prefix='memorient-preflight-') as temporary:
        root=Path(temporary)
        if kind=='membrane_surface_placement':
            from memorient.pdb_surface import build_scene
            spec=json.loads(store._one(inputs,'placement_spec').read_text())
            spec['pdb']=store._one(inputs,'structure').read_text();build_scene(spec)
        elif kind=='membrane_preparation':
            threads=manifest.get('parameters',{}).get('threads',4)
            if type(threads) is not int or not 1<=threads<=16:raise ValueError('Use 1–16 CPU threads')
            bind_profile(store._one(inputs,'preparation_profile'),store._one(inputs,'molecular_topology'),
                         store._one(inputs,'molecular_coordinates'),root)
            runtime()
        else:unpack(store._one(inputs,'membrane_bundle'),root)


def execute(store,manifest,inputs,run_dir,cancel_event=None,on_process=None,on_progress=None):
    kind=manifest['analysis_type'];output=run_dir/'outputs/membrane'
    command=[sys.executable,'-m','yauvi_platform.structural_workbench.membrane_workflow','--kind',kind,'--out',str(output)]
    for role in {r['role'] for d in definitions() if d['analysis_type']==kind for r in d['inputs'] if r['required']}:
        command.extend(['--'+role.replace('_','-'),str(store._one(inputs,role))])
    if kind=='membrane_preparation':
        python,gmx=runtime()
        command.extend(['--python',python,'--gmx',gmx,'--threads',str(manifest.get('parameters',{}).get('threads',4))])
    last=None
    def tick():
        nonlocal last
        receipt=output/'system/PREPARATION_RECEIPT.json'
        if on_progress is None or not receipt.is_file():return
        try:
            record=json.loads(receipt.read_text());stage=record.get('commands',[{}])[-1].get('stage',record.get('state','preparing'))
        except (ValueError,OSError,IndexError):return
        if stage!=last:
            last=stage;on_progress('membrane_preparation','Local stage: '+stage.replace('-',' '))
    code=store._execute(command,cwd=store.source_root,log_path=run_dir/'logs/membrane.log',
                        cancel_event=cancel_event,on_process=on_process,on_tick=tick)
    code=0 if code==0 else 2
    return [{'module_id':kind,'exit_code':code}],code


def recorded_manifest(store,analysis,run):
    from memorient.run_bundle import load_manifest,digest
    case=store.load(analysis)
    if run not in {r['run_id'] for r in case.get('runs',[])}:raise ValueError('Unknown recorded membrane run')
    path=store.artifact_path(analysis,run,'outputs/membrane/MEMBRANE_RUN.json')
    checks=json.loads(store.artifact_path(analysis,run,'CHECKSUMS.json').read_text())['files']
    if checks.get('outputs/membrane/MEMBRANE_RUN.json')!=digest(path):raise ValueError('Recorded membrane manifest changed')
    return path,load_manifest(path)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind',choices=sorted(TYPES),required=True);parser.add_argument('--out',required=True)
    for role in ['structure','placement-spec','preparation-profile','molecular-topology','molecular-coordinates','membrane-bundle']:
        parser.add_argument('--'+role)
    parser.add_argument('--python',default=sys.executable);parser.add_argument('--gmx',default='gmx');parser.add_argument('--threads',type=int,default=4)
    args=parser.parse_args(argv);output=Path(args.out);output.mkdir(parents=True,exist_ok=False)
    if args.kind=='membrane_surface_placement':
        from memorient.pdb_surface import build_scene
        from memorient.run_bundle import write_manifest
        spec=json.loads(Path(args.placement_spec).read_text());spec['pdb']=Path(args.structure).read_text()
        result=build_scene(spec)
        (output/'scene.json').write_text(json.dumps(result,allow_nan=False)+'\n')
        (output/'configuration.json').write_text(json.dumps({k:v for k,v in spec.items() if k!='pdb'},allow_nan=False)+'\n')
        shutil.copyfile(args.structure,output/'source.pdb');(output/'placed.pdb').write_text(result['pdb'])
        write_manifest(output,stage='synthetic' if spec.get('synthetic') is True else 'unrelaxed',
            artifacts={'scene':'scene.json','configuration':'configuration.json','source':'source.pdb','coordinates':'placed.pdb'},
            membrane_geometry=spec['surface'],organism=spec.get('organism','unknown'),
            sources=spec.get('sources',[]),assumptions=['User-supplied axis and placement; native sidedness unknown.'])
    elif args.kind=='membrane_result_import':unpack(args.membrane_bundle,output)
    else:
        bound=bind_profile(args.preparation_profile,args.molecular_topology,args.molecular_coordinates,output)
        code=subprocess.call([args.python,'-m','memorient.cli','md-prepare','--profile',str(bound),'--out',str(output/'system'),
                              '--gmx',args.gmx,'--threads',str(args.threads)])
        if code:return code
        from memorient.run_bundle import load_manifest
        doc=load_manifest(output/'system/MEMBRANE_RUN.json')
        for ref in doc['artifacts'].values():ref['path']='system/'+ref['path']
        (output/'MEMBRANE_RUN.json').write_text(json.dumps(doc,indent=2)+'\n')
    return 0


if __name__=='__main__':raise SystemExit(main())
