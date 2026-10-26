import json
import threading
import http.client
import numpy as np
import pytest
from http.server import ThreadingHTTPServer
from memorient.pdb_surface import read_pdb, build_scene
from memorient.surface_viewer import ViewerHandler
from memorient.cli import main

PDB = '\n'.join(f'ATOM  {i+1:5d}  CA  ALA A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  ' for i,(x,y,z) in enumerate([(0,0,-15),(2,0,0),(0,0,15)]))+'\nEND\n'
def spec():
    return {'pdb':PDB,'surface':{'kind':'sphere','size':[100],'thickness':30},'placement':{'axis':[0,0,1],'u':0,'v':0}}

def test_pdb_roundtrip_and_records():
    s=spec();s['pdb']+='ANISOU    1\n'
    result=build_scene(s)
    xyz=read_pdb(result['pdb'])[2]
    exact=np.array(result['coordinates'])
    np.testing.assert_allclose(xyz,exact,atol=.00051)
    assert result['atom_count']==3 and len(result['ca_atoms'])==3
    assert 'ANISOU' not in result['pdb']
    assert len(result['input_sha256'])==64
    assert result['pdb'].splitlines()[0][:30]==PDB.splitlines()[0][:30]

@pytest.mark.parametrize('text',['','ATOM      1  CA  ALA A   1       nan     0.0     0.0','MODEL        1\n'+PDB+'ENDMDL\nMODEL        2\n'+PDB,'ATOM      1'])
def test_invalid_pdb(text):
    with pytest.raises(ValueError): read_pdb(text)

def test_cli_pdb(tmp_path):
    p=tmp_path/'receptor.pdb';p.write_text(PDB)
    s=spec();del s['pdb'];j=tmp_path/'spec.json';j.write_text(json.dumps(s));out=tmp_path/'scene.json'
    args=['place-pdb','--pdb',str(p),'--spec',str(j),'--out',str(out)]
    assert main(args)==0
    original=out.read_bytes()
    assert main(args)==2 and out.read_bytes()==original and p.read_text()==PDB

def test_local_http():
    server=ThreadingHTTPServer(('127.0.0.1',0),ViewerHandler)
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    host=f'127.0.0.1:{server.server_port}'
    def request(method,path,body=None,headers=None):
        c=http.client.HTTPConnection(host);c.request(method,path,body,headers or {});r=c.getresponse();data=r.read();c.close();return r,data
    try:
        r,html=request('GET','/');assert r.status==200 and b'src="3Dmol-min.js"' in html
        assert "connect-src 'self'" in r.getheader('Content-Security-Policy')
        r,data=request('GET','/catalogue');assert r.status==200 and len(json.loads(data)['organisms'])==3
        r,data=request('GET','/catalogue/protein/2OMF');assert r.status in (200,503)
        if r.status==200:assert json.loads(data)['view_atoms']==8481
        r,data=request('GET','/catalogue/organism/bsubtilis');assert r.status==200 and len(json.loads(data)['layers'])==2
        r,data=request('GET','/explicit/ecoli-k12');assert r.status in (200,503)
        if r.status==200:
            assert json.loads(data)['full_atoms']==150032
            r,data=request('GET','/explicit/ecoli-k12.cif');assert r.status==200 and data.startswith(b'data_ecoli_k12_outer_membrane')
        assert request('GET','/explicit/../../sources/K12-AA-rep1.psf')[0].status==404
        r,data=request('GET','/md/ecoli-k12')
        assert r.status==200 and json.loads(data)['equilibration'].startswith('not_started')
        assert request('GET','/md/../../runs/system.top')[0].status==404
        assert request('POST','/md/ecoli-k12','{}',{'Content-Type':'application/json'})[0].status==404
        assert request('GET','/md/ecoli-k12',headers={'Origin':'https://evil.example'})[0].status==403
        assert request('GET','/catalogue/protein/unknown')[0].status==404
        assert request('GET','/../../etc/passwd')[0].status==404
        assert request('GET','/',headers={'Host':'evil.example'})[0].status==403
        assert request('POST','/place',json.dumps(spec()),{'Content-Type':'application/json','Origin':'https://evil.example'})[0].status==403
        r,data=request('POST','/place',json.dumps(spec()),{'Content-Type':'application/json'})
        assert r.status==200 and json.loads(data)['atom_count']==3
        assert request('POST','/place','{}',{'Content-Type':'application/json'})[0].status==400
        context={'experiment_id':'toy','species':['A'],'cell_size_um':1,'frames':[{'time_h':0,'statistics':{'species_total_biomass':{'A':1}}}]}
        r,data=request('POST','/marse-context',json.dumps({'text':json.dumps(context)}),{'Content-Type':'application/json'})
        assert r.status==200 and json.loads(data)['integration_state']=='read_only_context_import'
    finally:
        server.shutdown();server.server_close();thread.join()


def test_cli_retains_crlf_identity(tmp_path):
    import hashlib
    raw=PDB.replace('\n','\r\n').encode()
    p=tmp_path/'crlf.pdb';p.write_bytes(raw)
    j=tmp_path/'spec.json';s=spec();del s['pdb'];j.write_text(json.dumps(s));out=tmp_path/'out.json'
    assert main(['place-pdb','--pdb',str(p),'--spec',str(j),'--out',str(out)])==0
    assert json.loads(out.read_text())['input_sha256']==hashlib.sha256(raw).hexdigest()
