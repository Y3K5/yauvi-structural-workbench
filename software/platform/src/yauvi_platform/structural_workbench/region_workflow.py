"""Build a self-contained, source-bound region case from local inputs."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import gemmi

from .biological_case import BiologicalCase
from .regions import canonical, view_binding


def auto_mappings(doc, sequence, accession=None):
    """Unique exact canonical correspondence; ambiguous matches stay unmapped.

    Modified component positions are separately withheld from annotation mapping
    by the region adapter. This correspondence does not establish strain identity.
    """
    rows=[]
    entities=doc.get('_entity_poly.entity_id',[])
    polymers=doc.get('_entity_poly.pdbx_seq_one_letter_code_can',[])
    if len(entities) != len(polymers):raise ValueError('Incomplete polymer identity')
    from .regions import AA, _rows
    scheme=_rows(doc,'_pdbx_poly_seq_scheme.')
    differences=_rows(doc,'_struct_ref_seq_dif.')
    for entity, raw in zip(entities,polymers):
        polymer=''.join(raw.split())
        if not polymer:
            continue
        excluded={}
        if accession:
            for difference in differences:
                if difference.get('pdbx_seq_db_accession_code') != accession:continue
                label=difference.get('seq_num')
                matching=[r for r in scheme if r.get('entity_id') == entity and r.get('seq_id') == label and
                    r.get('pdb_strand_id') == difference.get('pdbx_pdb_strand_id') and r.get('mon_id') == difference.get('mon_id')]
                if matching and label and label.isdecimal() and AA.get(difference.get('db_mon_id')):
                    excluded[int(label)]=AA[difference['db_mon_id']]
        if len(polymer)<30 and polymer != sequence and not excluded:continue
        offsets=[i for i in range(len(sequence)-len(polymer)+1) if all(
            sequence[i+j] == (excluded[j+1] if j+1 in excluded else aa) for j,aa in enumerate(polymer))]
        if len(offsets)==1:
            # Split exact segments around deposited differences. No changed
            # component is normalized into an annotated reference residue.
            start=None
            for position in range(1,len(polymer)+2):
                available=position<=len(polymer) and position not in excluded
                if available and start is None:start=position
                if not available and start is not None:
                    rows.append({'entity_id':entity,'label_start':start,'reference_start':offsets[0]+start,'length':position-start})
                    start=None
    return rows


def build_case(coordinates, reference_record, output, *, declaration=None, model=0, assembly='asu', conformer='auto'):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    import io
    output=Path(output)
    if output.exists() and any(output.iterdir()):raise ValueError('Region output directory must be empty')
    coordinate_bytes=Path(coordinates).read_bytes()
    reference_bytes=Path(reference_record).read_bytes()
    doc=MMCIF2Dict(io.StringIO(coordinate_bytes.decode()))
    if not doc.get('_entry.id'):raise ValueError('Region explorer requires mmCIF with deposited entry identity')
    reference=json.loads(reference_bytes)
    sequence=reference.get('sequence',{}).get('value')
    accession=reference.get('primaryAccession')
    organism=reference.get('organism',{}).get('scientificName')
    if not sequence or not accession or not organism:raise ValueError('Complete reference identity is required')
    sequence_sha=hashlib.sha256(sequence.encode('ascii')).hexdigest()
    coordinate_sha=hashlib.sha256(coordinate_bytes).hexdigest()
    mappings=auto_mappings(doc,sequence,accession)
    assertions=[]
    declared_bytes=None
    if declaration:
        declared_bytes=Path(declaration).read_bytes()
        declared=json.loads(declared_bytes)
        if declared.get('schema_version') != '1.0' or declared.get('binding') != {
                'coordinate_sha256':coordinate_sha,'sequence_sha256':sequence_sha}:
            raise ValueError('Region declaration coordinate/reference identity mismatch')
        mappings=declared.get('sequence_mappings',mappings)
        assertions=declared.get('assertions',[])
    gemmi.make_structure_from_block(gemmi.cif.read_string(coordinate_bytes.decode()).sole_block())
    models=sorted(set(doc.get('_atom_site.pdbx_PDB_model_num',[])),key=int)
    if type(model) is not int or not 0 <= model < len(models):raise ValueError('Select a deposited model index')
    output.mkdir(parents=True,exist_ok=True)
    sources=[]
    for sid,name,raw in [('coordinates','coordinates.cif',coordinate_bytes),('reference','reference.json',reference_bytes)]:
        (output/name).write_bytes(raw)
        sources.append({'id':sid,'path':name,'sha256':hashlib.sha256(raw).hexdigest(),
                        'label':'User-supplied local '+sid,'retrieved_at':'supplied_local; acquisition_date_not_declared'})
    if declared_bytes is not None:
        (output/'declaration.json').write_bytes(declared_bytes)
        sources.append({'id':'declaration','path':'declaration.json','sha256':hashlib.sha256(declared_bytes).hexdigest(),
                        'label':'User-supplied local region declaration','retrieved_at':'supplied_local; acquisition_date_not_declared'})
    record={'schema_version':'1.2','id':'recorded-regions','title':accession+' · local region explorer',
        'protein':{'accession':accession,'organism':organism,'sequence':sequence,'sequence_sha256':sequence_sha,'source_id':'reference'},
        'sources':sources,'structures':[{'id':'coordinates','entry_id':doc['_entry.id'][0],
            'source_id':'coordinates','sequence_mappings':mappings}],
        'assertions':assertions,'evidence':[],'membranes':[],
        'limitations':['Coordinate-to-reference correspondence does not establish native strain, function or observed state.',
                        'Unmapped and modified positions remain unhighlighted. No membrane placement is inferred.']}
    (output/'case.json').write_bytes(canonical(record))
    case=BiologicalCase(output)
    view=case.view('coordinates',models[model],assembly,conformer)
    if not mappings:view['warnings'].append('No unique exact sequence correspondence; annotation regions stay unmapped.')
    return case,view


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--structure',required=True)
    parser.add_argument('--reference-record',required=True)
    parser.add_argument('--declaration')
    parser.add_argument('--out',required=True)
    parser.add_argument('--model',type=int,default=0)
    parser.add_argument('--assembly',default='asu')
    parser.add_argument('--conformer',default='auto')
    parser.add_argument('--qc-manifest',help='Separately recorded StructQC output from the same coordinates and model.')
    args=parser.parse_args(argv)
    try:
        output=Path(args.out)
        output.mkdir(parents=True,exist_ok=True)
        case,view=build_case(args.structure,args.reference_record,output/'CASE',declaration=args.declaration,
            model=args.model,assembly=args.assembly,conformer=args.conformer)
        if args.qc_manifest and args.assembly == 'asu' and args.conformer == 'auto':
            # Retain the original QC result and a separate method/binding receipt.
            # It remains checksum verification only in region replay.
            raw=Path(args.qc_manifest).read_bytes();qc=json.loads(raw)
            if qc.get('coordinate',{}).get('sha256') != view['coordinate_sha256'] or qc['coordinate'].get('selected_model') != args.model:
                raise ValueError('StructQC coordinate/model identity mismatch')
            binding={k:view[k] for k in ('structure_id','coordinate_sha256','sequence_sha256','model_id','assembly_id')}
            import structqc.chirality
            run={'binding':binding,'output_sha256':hashlib.sha256(raw).hexdigest(),'method':'structqc.structure_quality',
                'method_source_sha256':hashlib.sha256(Path(structqc.chirality.__file__).read_bytes()).hexdigest()}
            for sid,name,value in [('structqc','STRUCTURE_MANIFEST.json',raw),('structqc-run','STRUCTQC_RECEIPT.json',canonical(run))]:
                (case.directory/name).write_bytes(value)
                case.document['sources'].append({'id':sid,'path':name,'label':'Recorded '+sid,'retrieved_at':'local_run',
                    'sha256':hashlib.sha256(value).hexdigest()})
            case.document['evidence'].append({'id':'structure-quality','engine':'structqc','source_id':'structqc','run_source_id':'structqc-run','binding':binding})
            (case.directory/'case.json').write_bytes(canonical(case.document))
            case=BiologicalCase(case.directory);view=case.view('coordinates',view['model_id'],args.assembly,args.conformer)
        (output/'REGION_VIEW.json').write_bytes(canonical(view))
        (output/'REGION_MANIFEST.json').write_bytes(canonical({'schema_version':'1.0','binding':view_binding(view),
            'method':view['region_method'],'sources':case.summary()['sources'],
            'limitations':view['limitations']}))
        import csv
        with (output/'REGIONS.tsv').open('w',newline='') as f:
            writer=csv.writer(f,delimiter='\t')
            writer.writerow(['region_id','title','basis','mapping_state','observed','expected','source_id'])
            for region in view['regions']:
                writer.writerow([region['id'],region['title'],region['basis'],region['mapping_state'],
                    region['coverage']['observed'],region['coverage']['expected'],region['source_id']])
        print(json.dumps({'state':'recorded_regions','binding':view_binding(view),'region_count':len(view['regions'])}))
        return 0
    except (OSError,ValueError,KeyError,TypeError,RuntimeError) as exc:
        print('Region explorer blocked: '+str(exc))
        return 2


if __name__=='__main__':raise SystemExit(main())
