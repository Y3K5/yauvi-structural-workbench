"""Bundled synthetic inputs demonstrate execution, not external validation."""
from importlib.resources import files, as_file
from yauvi_platform.structural_workbench import StructuralAnalysisStore

def create_example(workspace,analysis_id='qc-example',*,without_validation=False):
    store=StructuralAnalysisStore(workspace)
    store.create(analysis_id,analysis_type='structure_qc',question='Inspect a synthetic coordinate model and its declared evidence.',subject_id='synthetic-example')
    for role,name in [('structure','model.pdb'),('provenance','provenance.json'),('reference_fasta','reference.fasta'),('validation_report','validation.json')]:
        if without_validation and role=='validation_report':continue
        with as_file(files('yauvi_structural_workbench').joinpath('examples',name)) as path:
            store.add_file(analysis_id,role=role,path=path)
    return store.snapshot(analysis_id)

def create_region_example(workspace,analysis_id='region-example'):
    store=StructuralAnalysisStore(workspace)
    store.create(analysis_id,analysis_type='region_explorer',question='Inspect synthetic regions, alternate conformers and missing coverage.',subject_id='synthetic-regions')
    for role,name in [('structure','regions.cif'),('reference_record','regions-reference.json')]:
        with as_file(files('yauvi_structural_workbench').joinpath('examples',name)) as path:
            store.add_file(analysis_id,role=role,path=path)
    return store.snapshot(analysis_id)
