"""Tests for memorient.domain_orient. Offline tests use invented coordinates only.
The public OPM comparison runs only with MEMORIENT_NETWORK_TESTS=1 (downloads public RCSB and OPM files)
or MEMORIENT_PDB_DIR=<dir with rcsb_<id>.pdb and opm_<id>.pdb>."""
import os, sys, math
import numpy as np
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from memorient.domain_orient import orient, read_pdb_ca, find_barrel, sse_psea


def invented_barrel(n=12, polar_outside=False, with_n_domain=True):
    """An idealised n-strand antiparallel barrel along +z (9 residues per strand, 3.2 A rise, radial zig-zag),
    long loops on +z (outside), 2-residue turns on -z, and an optional 60-residue polar helical domain below."""
    R = 4.8 / (2 * math.sin(math.pi / n)); ca, seq = [], []
    if with_n_domain:
        for i in range(60):
            a = math.radians(100 * i); ca.append([2.3 * math.cos(a) + 3, 2.3 * math.sin(a), -45 + 1.5 * i * .4]); seq.append('EKSQ'[i % 4])
    for k in range(n):
        th = 2 * math.pi * k / n; up = k % 2 == 0
        zs = [-12.8 + 3.2 * j for j in range(9)]
        if not up: zs = zs[::-1]
        for j, z in enumerate(zs):
            out = (j % 2 == 0)
            r = R + (0.9 if out else -0.9)
            ca.append([r * math.cos(th), r * math.sin(th), z])
            face = 'S' if out == polar_outside else 'L'
            seq.append('Y' if out and abs(z) > 9 and not polar_outside else face)
        th2 = 2 * math.pi * (k + .5) / n
        if up:     # long extracellular loop
            for m in range(10):
                h = 8 * math.sin(math.pi * (m + 1) / 11); ca.append([(R + 2) * math.cos(th2), (R + 2) * math.sin(th2), 14 + h]); seq.append('GDNT'[m % 4])
        else:      # short periplasmic turn
            for m in range(2):
                ca.append([(R + 1) * math.cos(th2), (R + 1) * math.sin(th2), -15 - m]); seq.append('GN'[m])
    return np.array(ca, float), ''.join(seq)


def rot(ax, deg):
    ax = np.asarray(ax, float) / np.linalg.norm(ax); a = math.radians(deg); K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + math.sin(a) * K + (1 - math.cos(a)) * K @ K


def test_invented_barrel_axis_side_and_depth():
    ca, seq = invented_barrel()
    r = orient(ca, seq)
    assert r['method'].startswith('beta_barrel_core')
    n = np.array(r['normal'])
    assert math.degrees(math.acos(min(1, n[2]))) < 15          # +normal = outside = +z (long loops, N-terminal domain below)
    assert abs(r['origin'][2]) < 3.0                           # slab centred on the barrel
    assert r['confidence'] in ('high', 'medium')


def test_rotation_equivariance():
    ca, seq = invented_barrel()
    M = rot([1, 2, 3], 57)
    r = orient(ca @ M.T + np.array([10, -4, 7]), seq)
    n = np.array(r['normal'])
    assert math.degrees(math.acos(min(1, n @ (M @ np.array([0, 0, 1.0]))))) < 15


def test_soluble_barrel_rejected():
    ca, seq = invented_barrel(polar_outside=True, with_n_domain=False)
    r = orient(ca, seq)
    assert not r['method'].startswith('beta_barrel_core')
    assert 'rejected_barrel' in r


def test_soluble_helix_bundle_has_no_core():
    ca, seq = [], ''
    for h in range(4):
        for i in range(20):
            a = math.radians(100 * i); ca.append([2.3 * math.cos(a) + 10 * (h % 2), 2.3 * math.sin(a) + 10 * (h // 2), 1.5 * i]); seq += 'EKAQ'[i % 4]
    r = orient(np.array(ca), seq)
    assert r['method'] == 'none'


OPM = [('4c00', 'A'), ('4k3b', 'A'), ('1bxw', 'A'), ('2omf', 'A'), ('1qd6', 'C')]


def _public_files(tmp_path, pid):
    """Public RCSB + OPM files: from MEMORIENT_PDB_DIR (rcsb_<id>.pdb, opm_<id>.pdb) or downloaded."""
    d = os.environ.get('MEMORIENT_PDB_DIR')
    if d:
        return os.path.join(d, f'rcsb_{pid}.pdb'), os.path.join(d, f'opm_{pid}.pdb')
    import ssl, urllib.request
    try:
        import certifi; ctx = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        ctx = ssl.create_default_context()
    out = []
    for name, url in (('rcsb', f'https://files.rcsb.org/download/{pid.upper()}.pdb'), ('opm', f'https://opm-assets.storage.googleapis.com/pdb/{pid}.pdb')):
        f = tmp_path / f'{name}.pdb'; f.write_bytes(urllib.request.urlopen(url, context=ctx, timeout=60).read()); out.append(str(f))
    return out


@pytest.mark.skipif(os.environ.get('MEMORIENT_NETWORK_TESTS') != '1' and not os.environ.get('MEMORIENT_PDB_DIR'),
                    reason='set MEMORIENT_NETWORK_TESTS=1 (download) or MEMORIENT_PDB_DIR (local public files)')
@pytest.mark.parametrize('pid,chain', OPM)
def test_public_opm_agreement(tmp_path, pid, chain):
    f1, f2 = _public_files(tmp_path, pid)
    ca, seq, num = read_pdb_ca(f1, chain); oca, _, onum = read_pdb_ca(f2, chain)
    common = sorted(set(num) & set(onum)); A = ca[[list(num).index(k) for k in common]]; B = oca[[list(onum).index(k) for k in common]]
    pa, pb = A.mean(0), B.mean(0); U, S, Vt = np.linalg.svd((A - pa).T @ (B - pb)); d = np.sign(np.linalg.det(Vt.T @ U.T)); Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    r = orient(ca, seq, num); n = Rm @ np.array(r['normal']); o = Rm @ (np.array(r['origin']) - pa) + pb
    assert math.degrees(math.acos(min(1, n[2]))) < 15     # OPM puts the periplasmic side at -z for outer-membrane proteins
    assert abs(o[2]) < 3.0
