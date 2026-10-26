"""Read-only geometry checks for admitted, orthorhombic CHARMM systems."""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .md_preparation import topology_summary, sha256


def check_geometry(run):
    run=Path(run);summary=topology_summary((run/'system.top').read_text())
    rows=(run/'em.gro').read_text().splitlines();n=int(rows[1]);atoms=rows[2:2+n]
    xyz=np.array([[float(line[a:b]) for a,b in ((20,28),(28,36),(36,44))] for line in atoms])
    box=np.array([float(x) for x in rows[-1].split()])
    if len(box)!=3 or not np.isfinite(box).all() or np.any(box<=0):
        raise ValueError('Only finite positive orthorhombic boxes are supported')
    expected=[];masses=[];molecules=[];waters=[];offset=0;mid=0
    for name,count in summary['sequence']:
        for _ in range(count):
            mid+=1;aa=summary['types'][name]
            expected.extend(a['name'] for a in aa);masses.extend(a['mass'] for a in aa)
            molecules.extend([mid]*len(aa))
            if name=='SOL':
                if len(aa)!=3:raise ValueError('Only three-site waters are supported')
                waters.append(offset)
            offset+=len(aa)
    if (n!=summary['atoms'] or len(atoms)!=n or [s[10:15].strip() for s in atoms]!=expected
        or not np.isfinite(xyz).all() or abs(summary['charge_e'])>1e-4):
        raise ValueError('Final coordinate count, order, finiteness or charge check failed')
    heavy=np.array(masses)>2;mol=np.array(molecules)[heavy]
    pairs=cKDTree(xyz[heavy]%box,boxsize=box).query_pairs(.15,output_type='ndarray')
    clashes=int(np.count_nonzero(mol[pairs[:,0]]!=mol[pairs[:,1]]))
    oh_error=hh_error=0.
    if waters:
        water=xyz[np.array(waters)[:,None]+np.arange(3)[None,:]]
        oh=water[:,1:]-water[:,:1];oh-=box*np.round(oh/box)
        hh=water[:,2]-water[:,1];hh-=box*np.round(hh/box)
        oh_error=float(np.max(abs(np.linalg.norm(oh,axis=2)-.09572)))
        hh_error=float(np.max(abs(np.linalg.norm(hh,axis=1)-.15139)))
    passed=clashes==0 and oh_error<.002 and hh_error<.002
    return {'passed':passed,'atoms':n,'atom_order_matches_topology':True,'charge_e':summary['charge_e'],
        'coordinate_sha256':sha256(run/'em.gro'),'topology_sha256':sha256(run/'system.top'),
        'box_nm':box.tolist(),'intermolecular_heavy_pairs_below_1_5_angstrom':clashes,
        'waters_checked':len(waters),'max_water_OH_error_nm':oh_error,'max_water_HH_error_nm':hh_error,
        'claim_limit':'Numerical geometry only; no ring-penetration, stability or biological validation'}
