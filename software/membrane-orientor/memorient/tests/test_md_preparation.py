"""Failure boundaries for the real preparation runner, using tiny invented files."""
import json
from pathlib import Path

import pytest

from memorient.md_preparation import (add_backbone_restraints, filter_core_water,
    load_profile, minimization_result, sha256, topology_summary,require_explicit_pair_parameters)


TOP = '''[ atomtypes ]
C 6 12 0 A .35 .2
O 8 16 0 A .31 .6
[ moleculetype ]
; name nrexcl
Protein 3
[ atoms ]
1 C 1 ALA CA 1 0 12
2 C 1 ALA CB 1 0 12
[ bonds ]
1 2 1 .15 1000
[ moleculetype ]
SOL 2
[ atoms ]
1 O 1 SOL OW 1 0 16
[ system ]
invented test
[ molecules ]
Protein 1
SOL 1
'''


def row(res, name, xyz, index=1):
    return f'{1:5d}{res:<5}{name:>5}{index:5d}' + ''.join(f'{v:8.3f}' for v in xyz)


def profile(tmp_path):
    top = tmp_path/'input.top'; top.write_text(TOP)
    gro = tmp_path/'input.gro'; gro.write_text('invented\n3\n'+ '\n'.join([
        row('ALA','CA',(1,1,1)),row('ALA','CB',(1,1,1)),row('SOL','OW',(1,1,1))])+'\n5 5 5\n')
    p = dict(id='invented', organism='invented organism', strain='test', membrane_layer='test',
             force_field='CHARMM36 invented test settings', sources=['test source'], assumptions=['invented'],
             water_policy={'center_z_nm':2.5,'exclude_half_width_nm':1.4},
             protein_molecule_types=['Protein'], admission='parameterized_inputs_reviewed',
             stop_after='minimization', equilibration_authorized=False,
             minimization_protocol={'contact_iterations':250},
             inputs={k:dict(path=f.name,sha256=sha256(f)) for k,f in [('topology',top),('coordinates',gro)]})
    path = tmp_path/'profile.json'; path.write_text(json.dumps(p)); return path,p


def test_tampering_and_equilibration_are_held(tmp_path):
    path,p = profile(tmp_path)
    assert load_profile(path)['organism'] == 'invented organism'
    p['equilibration_authorized'] = True; path.write_text(json.dumps(p))
    with pytest.raises(ValueError, match='before equilibration'): load_profile(path)
    p['equilibration_authorized'] = False; path.write_text(json.dumps(p))
    (tmp_path/'input.top').write_text(TOP+';tamper\n')
    with pytest.raises(ValueError, match='checksum'): load_profile(path)
    p['force_field']='unsupported force field'; path.write_text(json.dumps(p))
    with pytest.raises(ValueError,match='CHARMM36-family'):load_profile(path)


def test_unknown_chemistry_and_order_are_held(tmp_path):
    path,p = profile(tmp_path)
    p['admission'] = 'held_missing_lipid_parameters'; path.write_text(json.dumps(p))
    with pytest.raises(ValueError, match='held'): load_profile(path)
    path,p = profile(tmp_path)
    gro = tmp_path/'input.gro'; gro.write_text(gro.read_text().replace('   CA','   XX'))
    p['inputs']['coordinates']['sha256'] = sha256(gro); path.write_text(json.dumps(p))
    with pytest.raises(ValueError, match='order'): load_profile(path)


def test_negative_lj_is_rejected_and_restraints_use_local_indices():
    with pytest.raises(ValueError, match='Negative'): topology_summary(TOP.replace('.35 .2','.35 -.2'))
    output, counts = add_backbone_restraints(TOP,['Protein'])
    assert counts == {'Protein':1}
    protein, water = output.split('[ moleculetype ]\nSOL')
    assert '[ position_restraints ]' in protein and '1 1 1000 1000 1000' in protein
    assert '[ position_restraints ]' not in water
    assert topology_summary(output)['atoms'] == 3


def test_charmm_special_pairs_cannot_fall_back_to_generic_mixing():
    # Invented two-atom pair checks parameter coverage only, not a valid molecule.
    text=TOP.replace('[ bonds ]','[ pairs ]\n1 2 1\n[ bonds ]')
    with pytest.raises(ValueError,match='special 1-4'):require_explicit_pair_parameters(text)
    text=text.replace('[ moleculetype ]','[ pairtypes ]\nC C 1 .28 .1\n[ moleculetype ]',1)
    assert require_explicit_pair_parameters(text)==1


def test_core_water_policy_keeps_pores_and_complete_waters(tmp_path):
    waters=[]
    for xyz in [(1,1,2.5),(3,3,2.5),(3,3,4.5)]:
        waters += [row('SOL',name,xyz) for name in ['OW','HW1','HW2']]
    source=tmp_path/'in.gro'; dest=tmp_path/'out.gro'
    source.write_text('invented\n9\n'+'\n'.join(waters)+'\n5 5 5\n')
    result=filter_core_water(source,dest,dict(center_z_nm=2.5,exclude_half_width_nm=1.4,
        pore_cylinders_nm=[dict(xy=[1,1],radius=.5)]),0)
    assert result['removed_core_waters'] == 1
    assert result['added_retained_waters'] == 2
    assert int(dest.read_text().splitlines()[1]) == 6


def test_step_limit_is_not_convergence():
    log='Maximum force = 1.2345e+04 on atom 10\nPotential Energy = -9.01e+06\n'
    result=minimization_result(log)
    assert not result['converged'] and result['maximum_force_kj_mol_nm'] == 12345
    assert not minimization_result('Steepest Descents converged to Fmax < 1000 in 500 steps\n'+log)['converged']
    clean=log.replace('1.2345e+04','9.9e+02')
    assert minimization_result('Steepest Descents converged to Fmax < 1000 in 500 steps\n'+clean)['converged']
    assert not minimization_result('Steepest Descents converged to Fmax < 1000\n')['converged']
    result=minimization_result('Steepest Descents converged to Fmax < 1000\n'+log+'LINCS WARNING\n')
    assert result['constraint_or_fatal_messages'] == ['LINCS WARNING']
