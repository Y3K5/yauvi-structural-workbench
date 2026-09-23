"""The displayed membrane must be one complete fit in the input coordinate frame."""
import numpy as np
import pytest
from memorient.contexts import get_context
from memorient.orientor import orient_structure
from synthetic import make_tm_helix, make_soluble_blob, random_rotation


def test_input_plane_transforms_with_coordinates_and_keeps_unsigned_side():
    structure=make_tm_helix(seed=4)
    topology={'spans':[{'chain_id':'A','start_auth_seq_id':11,'end_auth_seq_id':31}],
              'source':{'id':'synthetic','citation':'Synthetic geometry'}}
    result=orient_structure(structure,get_context('tm_receptor'),validate=False,n_points=80,topology_evidence=topology)
    pose=result.to_dict()['input_coordinate_membrane']
    assert pose['frame']=='input_coordinates' and pose['sidedness']=='unknown'
    rotation=random_rotation(44);translation=np.array([17.,-8.,4.])
    moved=structure.transformed(rotation,t=translation)
    changed=orient_structure(moved,get_context('tm_receptor'),validate=False,n_points=80,topology_evidence=topology).to_dict()['input_coordinate_membrane']
    normal=np.array(pose['normal']);center=np.array(pose['center'])
    moved_normal=np.array(changed['normal']);moved_center=np.array(changed['center'])
    assert abs(np.dot(rotation@normal,moved_normal))==pytest.approx(1,abs=1e-9)
    # A plane's center may differ by an in-plane vector; signed distance is the invariant.
    assert abs(np.dot(moved_center-(rotation@center+translation),moved_normal))<1e-8
    assert changed['half_thickness']==pytest.approx(pose['half_thickness'],abs=1e-8)
    np.testing.assert_allclose(np.abs((structure.ca-center)@normal),np.abs((moved.ca-moved_center)@moved_normal),atol=1e-8)


def test_soluble_route_has_no_input_membrane():
    result=orient_structure(make_soluble_blob(seed=3),get_context('soluble_secreted'),validate=False,n_points=80)
    assert result.to_dict()['input_coordinate_membrane']=={}
