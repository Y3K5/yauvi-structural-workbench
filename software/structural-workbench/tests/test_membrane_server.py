"""Loopback HTTP integration of the installed-style membrane routes and jobs."""
import hashlib
import io
import json
import threading
import time
from urllib.request import Request,urlopen
from urllib.error import HTTPError
import zipfile
import pytest
from yauvi_structural_workbench.server import WorkbenchServer


def test_native_membrane_flow_and_guard(tmp_path):
    server=WorkbenchServer(('127.0.0.1',0),tmp_path/'workspace')
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    base=f'http://127.0.0.1:{server.server_port}'
    def request(path,data=None,token=None):
        headers={}
        if token:headers['X-Yauvi-Token']=token
        if data is not None:headers['Content-Type']='application/json'
        with urlopen(Request(base+path,data=None if data is None else json.dumps(data).encode(),headers=headers),timeout=30) as response:return response.read()
    try:
        assert b'Organism Orient' in request('/organism-orient/')
        assert b'wbSave' in request('/organism-orient/workbench.js')
        assert json.loads(request('/organism-orient/connection'))['workbench']
        with pytest.raises(HTTPError) as error:request('/organism-orient/place',{})
        assert error.value.code==403
        token=json.loads(request('/api/session'))['token']
        spec={'pdb':'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\nEND\n',
              'surface':{'kind':'sphere','size':[100],'thickness':30},'placement':{'axis':[0,0,1],'u':0,'v':0},'synthetic':True}
        assert json.loads(request('/organism-orient/place',spec,token))['atom_count']==1
        store=server.store;store.create('membrane-http',analysis_type='membrane_surface_placement',question='Invented HTTP check')
        for role,name,raw in [('structure','fixture.pdb',spec['pdb'].encode()),('placement_spec','spec.json',json.dumps({k:v for k,v in spec.items() if k!='pdb'}).encode())]:
            upload=store.begin_ingest('membrane-http',role=role,file_name=name,size=len(raw),expected_sha256=hashlib.sha256(raw).hexdigest())
            store.ingest_chunk(upload['upload_id'],0,raw);store.finalize_ingest(upload['upload_id'])
        with pytest.raises(HTTPError):request('/api/analyses/membrane-http/run',{},token)
        job=json.loads(request('/api/analyses/membrane-http/run',{'experimental_acknowledged':True},token))
        deadline=time.monotonic()+40
        while time.monotonic()<deadline:
            row=next(r for r in server.jobs.list() if r['job_id']==job['job_id'])
            if row['state'] not in {'queued','running'}:break
            time.sleep(.05)
        assert row['state']=='completed',row
        path=f'/api/analyses/membrane-http/membrane-result/{row["run_id"]}'
        doc=json.loads(request(path+'/manifest'));assert doc['stage']=='synthetic'
        raw=request(path+'/bundle')
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            assert json.loads(archive.read('MEMBRANE_RUN.json'))==doc
            for ref in doc['artifacts'].values():assert hashlib.sha256(archive.read(ref['path'])).hexdigest()==ref['sha256']
        source=store.artifact_path('membrane-http',row['run_id'],'outputs/membrane/scene.json');source.write_text('{}')
        with pytest.raises(HTTPError):request(path+'/scene')
    finally:server.shutdown();server.server_close();thread.join(timeout=3)


@pytest.mark.parametrize('kind',['membrane_surface_placement','membrane_preparation','membrane_result_import'])
def test_all_membrane_jobs_require_experimental_ack(store,manager,kind):
    store.add_case('invented',analysis_type=kind)
    with pytest.raises(Exception,match='experimental'):manager.submit('invented')
