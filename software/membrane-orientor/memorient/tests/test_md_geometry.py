"""Invented three-site-water geometry acceptance; no dynamics or external files."""
import pytest
from memorient.md_geometry import check_geometry


def fixture(root,bad=False):
    (root/'system.top').write_text('''[ atomtypes ]
OW 8 16.0 0 A 0.3 0.5
HW 1 1.0 0 A 0 0
[ moleculetype ]
SOL 2
[ atoms ]
1 OW 1 SOL OW 1 -0.834 16
2 HW 1 SOL HW1 1 0.417 1
3 HW 1 SOL HW2 1 0.417 1
[ system ]
Invented water geometry
[ molecules ]
SOL 1
''')
    xyz=[(.2,.2,.2),(.4 if bad else .296,.2,.2),(.176,.293,.2)]
    rows=['Invented software fixture','3']
    for i,(name,p) in enumerate(zip(['OW','HW1','HW2'],xyz),1):
        rows.append(f'{1:5d}{"SOL":<5}{name:>5}{i:5d}'+''.join(f'{v:8.3f}' for v in p))
    (root/'em.gro').write_text('\n'.join(rows+['2.0 2.0 2.0'])+'\n')


def test_water_acceptance_and_bad_geometry_hold(tmp_path):
    fixture(tmp_path);assert check_geometry(tmp_path)['passed']
    fixture(tmp_path,bad=True);assert not check_geometry(tmp_path)['passed']


def test_atom_order_and_box_are_checked(tmp_path):
    fixture(tmp_path)
    p=tmp_path/'em.gro';p.write_text(p.read_text().replace('  HW1','  BAD'))
    with pytest.raises(ValueError,match='order'):check_geometry(tmp_path)
    fixture(tmp_path);p.write_text(p.read_text().replace('2.0 2.0 2.0','0.0 2.0 2.0'))
    with pytest.raises(ValueError,match='boxes'):check_geometry(tmp_path)
