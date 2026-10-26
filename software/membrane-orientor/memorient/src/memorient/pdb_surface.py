"""Strict local PDB coordinate import and shape placement. No network access."""
import hashlib
import numpy as np
from .surfaces import Surface, place_receptor


def read_pdb(text):
    if not isinstance(text, str) or len(text.encode('utf-8')) > 8_000_000:
        raise ValueError('PDB must be text under 8 MB')
    lines = text.splitlines()
    if sum(line.startswith('MODEL ') for line in lines) > 1:
        raise ValueError('Multiple models: supply a single-model PDB explicitly')
    indices, xyz, ca = [], [], []
    for i, line in enumerate(lines):
        if line[:6] not in ('ATOM  ', 'HETATM'):
            continue
        try:
            point = [float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))]
        except ValueError as exc:
            raise ValueError(f'Invalid PDB coordinates on line {i+1}') from exc
        if not np.isfinite(point).all():
            raise ValueError(f'Nonfinite PDB coordinate on line {i+1}')
        indices.append(i); xyz.append(point)
        if line[12:16].strip() == 'CA' and line[:6] == 'ATOM  ':
            ca.append({'atom_index':len(xyz)-1, 'label':f'{line[21:22] or "_"}:{line[22:27].strip()} {line[17:20].strip()}', 'coordinate':point})
    if not xyz:
        raise ValueError('No ATOM/HETATM coordinates found')
    if len(xyz) > 60000:
        raise ValueError('Viewer limit: 60000 atoms')
    return lines, indices, np.asarray(xyz), ca


def place_pdb(text, surface, placement):
    lines, indices, xyz, ca = read_pdb(text)
    options = dict(placement)
    options.setdefault('anchor', xyz.mean(axis=0).tolist())
    result = place_receptor(xyz, surface, **options)
    for index, point in zip(indices, result['coordinates']):
        fields = [f'{v:8.3f}' for v in point]
        if any(len(f) != 8 for f in fields):
            raise ValueError('Placed coordinates exceed PDB fixed-width limits')
        line = lines[index]
        lines[index] = line[:30] + ''.join(fields) + line[54:]
    # ANISOU/SIGATM/SIGUIJ refer to the original coordinate frame; do not export
    # these untransformed records. Other source metadata is not an output claim.
    omitted = ('ANISOU', 'SIGATM', 'SIGUIJ')
    output = '\n'.join(line for line in lines if line[:6] not in omitted) + '\n'
    return {**result, 'pdb':output, 'input_sha256':hashlib.sha256(text.encode('utf-8')).hexdigest(),
            'atom_count':len(xyz), 'ca_atoms':ca, 'source_centroid':xyz.mean(axis=0).tolist(),
            'import_notes':['All ATOM/HETATM records and alternate locations retained; multiple models rejected',
                            'Output coordinates rounded to PDB precision (0.001 Angstrom)',
                            'ANISOU/SIGATM/SIGUIJ records omitted because their frame was not transformed',
                            'Header assembly and symmetry metadata describe the source frame, not the placed model']}


def build_scene(spec):
    if not isinstance(spec, dict):
        raise ValueError('Request must be an object')
    surface = Surface(**spec['surface'])
    result = place_pdb(spec['pdb'], surface, spec['placement'])
    return {**result, 'surface_spec':spec['surface'], 'placement_spec':spec['placement'],
            'surface':surface.sample(24, 13), 'sampling':{'nu':24,'nv':13},
            'units':'angstrom', 'angle_units':'radian'}
