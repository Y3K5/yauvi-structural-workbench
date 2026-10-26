import json
import pytest
from memorient.organism_context import read_marse_context
from memorient.lipid_plan import plan_lipids


def test_marse_context():
    raw=json.dumps({'experiment_id':'invented','species':['Toy A'], 'cell_size_um':2,
                    'frames':[{'time_h':1,'statistics':{'species_total_biomass':{'Toy A':3}}}]})
    result=read_marse_context(raw)
    assert result['species']==['Toy A'] and result['last_species_total_biomass']=={'Toy A':3}
    assert result['integration_state']=='read_only_context_import'
    assert len(result['source_sha256'])==64
    assert 'surface' not in result

@pytest.mark.parametrize('update',[{'species':['Toy A','Toy A']},{'frames':[]},{'cell_size_um':0},{'cell_size_um':float('nan')}])
def test_bad_marse(update):
    raw={'experiment_id':'invented','species':['Toy A'], 'cell_size_um':2,'frames':[{'time_h':1,'statistics':{'species_total_biomass':{'Toy A':3}}}]}
    raw.update(update)
    with pytest.raises(ValueError):read_marse_context(json.dumps(raw))

def recipe():
    return {'schema_version':'organism_orient.lipid_recipe.v1','surface':{'kind':'sphere','size':[100],'thickness':30}}

def test_lipid_draft_held():
    result=plan_lipids(recipe())
    assert result['state']=='held_no_packing_adapter'
    assert 'force_field' in result['missing_inputs'] and 'leaflets.positive' in result['missing_inputs']

@pytest.mark.parametrize('fraction',[0,-1,.8,float('nan'),True])
def test_bad_lipid_fraction(fraction):
    r=recipe();r['leaflets']={'positive':[{'lipid_id':'invented','mole_fraction':fraction}]}
    with pytest.raises(ValueError):plan_lipids(r)

def test_complete_recipe_still_held():
    r=recipe();r.update({key:'invented' for key in ('representation','force_field','composition_source','builder','builder_version','template_manifest_sha256')})
    r['leaflets']={side:[{'lipid_id':'toy','mole_fraction':1}] for side in ('positive','negative')}
    result=plan_lipids(r)
    assert result['missing_inputs']==[] and result['state']=='held_no_packing_adapter'


def test_recipe_object_required():
    with pytest.raises(ValueError):plan_lipids([])
