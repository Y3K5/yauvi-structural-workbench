import json
import numpy as np
import pytest
from memorient.surfaces import Surface, place_receptor
from memorient.cli import main

@pytest.mark.parametrize('kind,size,u,v', [('plane',(100,80),10,5),('sphere',(100,),.8,.3),('cylinder',(100,200),.8,30),('ellipsoid',(100,120,150),.8,.3)])
@pytest.mark.parametrize('axis', [(0,0,1),(0,0,-1),(1,2,3)])
@pytest.mark.parametrize('side', [-1,1])
def test_rigid_placement(kind,size,u,v,axis,side):
    surface=Surface(kind,size,30,center=(4,5,6))
    xyz=np.array([[1,2,3],[3,7,1],[-2,1,4],[4,-1,2.]])
    result=place_receptor(xyz,surface,anchor=xyz[0],axis=axis,u=u,v=v,side=side,depth=2,spin=.7)
    out=np.array(result['coordinates']); r=np.array(result['rotation'])
    point,normal=surface.frame(u,v)
    np.testing.assert_allclose(out[0],point+2*side*normal)
    np.testing.assert_allclose(r @ (np.array(axis)/np.linalg.norm(axis)),side*normal,atol=1e-12)
    np.testing.assert_allclose(r.T@r,np.eye(3),atol=1e-12)
    assert np.linalg.det(r)==pytest.approx(1)
    np.testing.assert_allclose(np.linalg.norm(out[:,None]-out,axis=-1),np.linalg.norm(xyz[:,None]-xyz,axis=-1),atol=1e-12)

@pytest.mark.parametrize('kind,size',[('plane',(100,80)),('sphere',(100,)),('cylinder',(100,200)),('ellipsoid',(100,120,150))])
def test_interfaces(kind,size):
    data=Surface(kind,size,30).sample(8,5)
    np.testing.assert_allclose(np.linalg.norm(np.array(data['positive_interface'])-data['negative_interface'],axis=1),30)
    np.testing.assert_allclose(np.linalg.norm(data['normals'],axis=1),1)

@pytest.mark.parametrize('kind,size,thickness',[('unknown',(1,),1),('sphere',(0,),1),('sphere',(10,),20),('ellipsoid',(10,20,100),3),('plane',(10,10),float('nan'))])
def test_bad_surface(kind,size,thickness):
    with pytest.raises(ValueError): Surface(kind,size,thickness)

def test_invalid_placement():
    s=Surface('sphere',(100,),30)
    for update in [{'axis':[0,0,0]},{'coordinates':[[float('nan'),0,0]]},{'side':0},{'v':2}]:
        args=dict(coordinates=[[0,0,0]],anchor=[0,0,0],axis=[0,0,1],u=0,v=0)
        args.update(update)
        with pytest.raises(ValueError): place_receptor(surface=s,**args)

def test_cli(tmp_path):
    source=tmp_path/'spec.json'; output=tmp_path/'result.json'
    source.write_text(json.dumps({'schema_version':'memorient.surface.v1','surface':{'kind':'sphere','size':[100],'thickness':30}}))
    original=source.read_bytes()
    assert main(['surface','--spec',str(source),'--out',str(output)])==0
    result=output.read_bytes()
    assert json.loads(result)['scientific_state']=='geometric_model_only'
    assert main(['surface','--spec',str(source),'--out',str(output)])==2
    assert output.read_bytes()==result and source.read_bytes()==original
