"""Local, source-bound GROMACS preparation. No equilibration entry point.

Accepts an admitted standalone topology, not an organism name or guessed chemistry.
Each execution gets a new directory and preserves command logs and stage receipts.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def topology_summary(text):
    """Inspect standalone topology atom tables and molecule sequence, without includes."""
    section = ""; current = None; types = {}; sequence = []; atomtypes = {}
    for raw in text.splitlines():
        line = raw.split(';', 1)[0].strip()
        if not line or line.startswith('#'): continue
        if line.startswith('['):
            section = line.strip('[] ').strip(); continue
        f = line.split()
        if section == 'moleculetype':
            current = f[0]; types[current] = []
        elif section == 'atoms':
            types[current].append({'index': int(f[0]), 'type': f[1], 'resid':int(f[2]), 'residue': f[3],
                                   'name': f[4], 'charge': float(f[6]), 'mass': float(f[7])})
        elif section == 'molecules': sequence.append((f[0], int(f[1])))
        elif section == 'atomtypes': atomtypes[f[0]] = [float(f[-2]), float(f[-1])]
    if not sequence or not atomtypes: raise ValueError('Standalone topology is incomplete')
    if any(eps < 0 or sig < 0 for sig, eps in atomtypes.values()):
        raise ValueError('Negative GROMACS Lennard-Jones parameters')
    atoms = sum(len(types[name]) * n for name, n in sequence)
    charge = sum(sum(a['charge'] for a in types[name]) * n for name, n in sequence)
    return dict(atoms=atoms, charge_e=charge, types=types, sequence=sequence, atomtypes=atomtypes)


def load_profile(path):
    path = Path(path).resolve(); p = json.loads(path.read_text())
    for key in ('id', 'organism', 'strain', 'membrane_layer', 'force_field', 'sources',
                'assumptions', 'water_policy', 'protein_molecule_types'):
        if not p.get(key): raise ValueError(f'Missing reviewed profile field: {key}')
    if p.get('admission') != 'parameterized_inputs_reviewed':
        raise ValueError('Profile is held: parameterized inputs have not been reviewed')
    if p.get('stop_after') != 'minimization' or p.get('equilibration_authorized') is not False:
        raise ValueError('This preparation command stops before equilibration')
    if 'CHARMM36' not in str(p['force_field']).upper():
        raise ValueError('This workflow supports reviewed all-atom CHARMM36-family inputs only')
    if p['water_policy'].get('mode','planar_slab_with_pores') != 'planar_slab_with_pores':
        raise ValueError('Only a reviewed planar slab hydration policy is implemented')
    contact_steps=p.get('minimization_protocol',{}).get('contact_iterations')
    if not isinstance(contact_steps,int) or not 1<=contact_steps<=10000:
        raise ValueError('A reviewed contact_iterations count between 1 and 10000 is required')
    protocol=p['minimization_protocol']
    if protocol.get('final_block_iterations',10)!=10 or protocol.get('final_force_target_kj_mol_nm',1000)!=1000:
        raise ValueError('Final minimization uses ten-iteration blocks and force target1000')
    if not isinstance(protocol.get('max_final_blocks',30),int) or not 1<=protocol.get('max_final_blocks',30)<=30:
        raise ValueError('Use between1 and30 final minimization blocks')
    resolved = {}
    for kind in ('topology', 'coordinates'):
        ref = p['inputs'][kind]; f = (path.parent / ref['path']).resolve()
        if not f.is_file() or sha256(f) != ref['sha256']:
            raise ValueError(f'{kind} source checksum mismatch')
        resolved[kind] = f
    top = resolved['topology'].read_text()
    if re.search(r'^\s*#\s*include', top, re.M):
        raise ValueError('Use a complete standalone topology; unpinned includes are held')
    summary = topology_summary(top)
    if 'CHARMM' in str(p['force_field']).upper():
        require_explicit_pair_parameters(top,summary)
    if int(resolved['coordinates'].read_text().splitlines()[1]) != summary['atoms']:
        raise ValueError('Topology and coordinate atom counts differ')
    names = [line[10:15].strip() for line in resolved['coordinates'].read_text().splitlines()[2:-1]]
    expected = [a['name'] for name, n in summary['sequence'] for _ in range(n) for a in summary['types'][name]]
    if names != expected: raise ValueError('Topology and coordinate atom order differs')
    if summary['sequence'][-1][0] != 'SOL':
        raise ValueError('Seed SOL must be last, to keep all solvent contiguous for solvation/genion')
    for name in p['protein_molecule_types']:
        if name not in summary['types']: raise ValueError('Unknown restrained protein molecule type')
    p['_resolved'] = resolved
    return p


def require_explicit_pair_parameters(text, summary=None):
    """CHARMM special 1-4 terms must not silently use generic LJ mixing."""
    summary = summary or topology_summary(text)
    section=''; current=None; pairtypes=set(); pending=[]
    for raw in text.splitlines():
        line=raw.split(';',1)[0].strip()
        if not line or line.startswith('#'): continue
        if line.startswith('['): section=line.strip('[] '); continue
        f=line.split()
        if section=='moleculetype': current=f[0]
        elif section=='pairtypes': pairtypes.add(tuple(sorted(f[:2])))
        elif section=='pairs':
            pair=tuple(sorted(summary['types'][current][int(f[i])-1]['type'] for i in (0,1)))
            if len(f)<5: pending.append(pair)
    missing=set(pending)-pairtypes
    if missing: raise ValueError('CHARMM special 1-4 parameters are missing; generic mixing is held')
    return len(set(pending))


def add_backbone_restraints(text, protein_types, force=1000.):
    """Use molecule-local indices; preserve every parameter and molecule count."""
    summary = topology_summary(text)
    chunks = re.split(r'(?=^\[ moleculetype \])', text, flags=re.M)
    restrained = {}
    for i, chunk in enumerate(chunks):
        match = re.match(r'\[ moleculetype \]\s*\n(?:;[^\n]*\n)*([^\s;]+)', chunk)
        if not match or match[1] not in protein_types: continue
        name = match[1]
        ids = [a['index'] for a in summary['types'][name] if a['name'] in {'N','CA','C','O'}]
        if not ids: raise ValueError('No backbone atoms in protein molecule')
        # The final molecule includes system/molecules, so insert before those sections.
        pos = re.search(r'^\[ (system|molecules) \]', chunk, re.M)
        offset = pos.start() if pos else len(chunk)
        block = '\n#ifdef POSRES_PROTEIN_BACKBONE\n[ position_restraints ]\n; molecule-local backbone indices\n'
        block += ''.join(f'{idx} 1 {force:g} {force:g} {force:g}\n' for idx in ids)
        block += '#endif\n\n'
        chunks[i] = chunk[:offset] + block + chunk[offset:]; restrained[name] = len(ids)
    if set(restrained) != set(protein_types): raise ValueError('Failed to insert every protein restraint')
    return ''.join(chunks), restrained


def set_water_count(text, count):
    before, table = text.rsplit('[ molecules ]', 1)
    table, n = re.subn(r'^SOL\s+\d+\s*$', f'SOL {count}', table, flags=re.M)
    if n != 1: raise ValueError('Expected one SOL molecule-count entry')
    return before + '[ molecules ]' + table


def filter_core_water(source, destination, policy, solute_atoms):
    """Keep complete waters; explicitly declared pore cylinders override core exclusion."""
    import numpy as np
    rows = Path(source).read_text().splitlines(); n = int(rows[1]); atoms = rows[2:2+n]
    kept = list(atoms[:solute_atoms]); removed = 0
    remainder = atoms[solute_atoms:]
    if len(remainder) % 3: raise ValueError('Expected three-site added solvent')
    center = float(policy['center_z_nm']); half = float(policy['exclude_half_width_nm'])
    if not 0 < half < 5: raise ValueError('Explicit core exclusion half-width required')
    pores = policy.get('pore_cylinders_nm', [])
    box = np.array([float(x) for x in rows[-1].split()[:3]])
    for i in range(0, len(remainder), 3):
        water = remainder[i:i+3]
        if any(line[5:10].strip() != 'SOL' for line in water): raise ValueError('Non-SOL added solvent')
        xyz = np.array([float(water[0][a:b]) for a,b in ((20,28),(28,36),(36,44))])
        inside_pore = False
        for pore in pores:
            delta = xyz[:2] - np.array(pore['xy']); delta -= box[:2] * np.round(delta/box[:2])
            inside_pore |= float(delta @ delta) < float(pore['radius'])**2
        if abs(xyz[2]-center) < half and not inside_pore: removed += 1
        else: kept.extend(water)
    Path(destination).write_text(rows[0]+'\n'+f'{len(kept):5d}\n'+'\n'.join(kept)+'\n'+rows[-1]+'\n')
    return {'removed_core_waters': removed, 'added_retained_waters': (len(kept)-solute_atoms)//3,
            'policy': policy}


def minimization_result(log):
    force = re.search(r'Maximum force\s*=\s*([\d.eE+\-]+)', log)
    energy = re.search(r'Potential Energy\s*=\s*([\d.eE+\-]+)', log)
    maximum_force = float(force[1]) if force else None
    potential = float(energy[1]) if energy else None
    tolerance = re.search(r'converged to Fmax <\s*([\d.eE+\-]+)', log)
    converged = bool(tolerance)
    converged = (converged and maximum_force is not None and potential is not None
                 and math.isfinite(maximum_force) and math.isfinite(potential)
                 and maximum_force < float(tolerance[1]))
    warnings = [line.strip() for line in log.splitlines() if any(term in line for term in
                ('can not be settled', 'Constraint error', 'LINCS WARNING', 'Fatal error'))]
    return {'converged': converged, 'constraint_or_fatal_messages': warnings,
            'maximum_force_kj_mol_nm': maximum_force,
            'potential_energy_kj_mol': potential}


def read_minimization(run, stem='em'):
    run = Path(run)
    command_log = run/('minimize.command.log' if stem == 'em' else f'{stem}.command.log')
    return minimization_result((run/f'{stem}.log').read_text() + '\n' +
                               (command_log.read_text() if command_log.exists() else ''))


EM_MDP = '''; Local restrained minimization only; no temperature or pressure coupling.
integrator = steep
nsteps = 10000
emtol = 1000
emstep = 0.005
define = -DPOSRES_PROTEIN_BACKBONE
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw-switch = 1.0
rvdw = 1.2
DispCorr = no
constraints = h-bonds
pbc = xyz
nstenergy = 100
nstlog = 100
nstxout = 250
'''

# Flexible-water steepest-descent cleanup avoids forcing distorted starting
# contacts through SETTLE/LINCS. The final stage restores the rigid CHARMM water
# and hydrogen-bond constraints, and must independently converge without warnings.
CLEANUP_MDP = (EM_MDP.replace('nsteps = 10000','nsteps = 250').replace('emtol = 1000','emtol = 2000')
               .replace('emstep = 0.005','emstep = 0.001')
               .replace('define = -DPOSRES_PROTEIN_BACKBONE','define = -DPOSRES_PROTEIN_BACKBONE -DFLEXIBLE')
               .replace('constraints = h-bonds','constraints = none'))


def final_minimization_blocks(run,record,command,threads,max_blocks=30):
    """Bound adaptive step growth; retain G96 precision and one fixed restraint reference."""
    run=Path(run); start='cleanup.gro'; record['final_blocks']=[]
    for index in range(max_blocks):
        prefix=f'final-block-{index:03d}'; folder=run/prefix; folder.mkdir()
        mdp=(EM_MDP.replace('nsteps = 10000','nsteps = 10')
             .replace('emstep = 0.005','emstep = 0.001').replace('nstxout = 250','nstxout = 0'))
        (folder/'em.mdp').write_text(mdp)
        command(prefix+'/em-input',['grompp','-f',prefix+'/em.mdp','-c',start,'-r','neutralized.gro',
                                    '-p','system.top','-o',prefix+'/em.tpr','-po',prefix+'/em-resolved.mdp'])
        command(prefix+'/minimize',['mdrun','-deffnm',prefix+'/em','-c',prefix+'/em.g96',
                                    '-ntmpi','1','-ntomp',str(threads),'-pin','off','-v'])
        result=read_minimization(folder)
        record['final_blocks'].append({'directory':prefix,**result})
        if result['constraint_or_fatal_messages'] or result['converged']: break
        if result['maximum_force_kj_mol_nm'] is None or result['potential_energy_kj_mol'] is None:
            raise ValueError('Final minimization did not report finite force/energy')
        start=prefix+'/em.g96'
    for name in ['em.mdp','em-resolved.mdp','em.tpr','em.g96','em.log','em.edr',
                 'em-input.command.log','minimize.command.log']:
        shutil.copy2(folder/name,run/name)
    command('export-final-coordinates',['editconf','-f','em.g96','-o','em.gro'])
    return result


def prepare_and_minimize(profile_path, destination, gmx='gmx', threads=4):
    """Run only the explicitly authorized stages, with no warning override flags."""
    p = load_profile(profile_path)
    if not isinstance(threads, int) or not 1 <= threads <= 16: raise ValueError('Use 1–16 CPU threads')
    run = Path(destination).resolve(); run.mkdir(parents=True, exist_ok=False)
    record = {k:v for k,v in p.items() if not k.startswith('_')}
    record.update(state='preparing', commands=[], started_utc=datetime.now(timezone.utc).isoformat())
    def save(): (run/'PREPARATION_RECEIPT.json').write_text(json.dumps(record, indent=2)+'\n')
    def command(stage, argv, stdin=None):
        log = run/f'{stage}.command.log'
        record['commands'].append({'stage':stage, 'argv':argv, 'log':str(log.relative_to(run))}); save()
        with log.open('w') as out:
            result = subprocess.run([gmx, *argv], cwd=run, input=stdin, text=True,
                                    stdout=out, stderr=subprocess.STDOUT)
        if result.returncode:
            record.update(state='failed', failed_stage=stage); save()
            raise RuntimeError(f'{stage} failed; see {log}')
    try:
        text = p['_resolved']['topology'].read_text()
        text, restraints = add_backbone_restraints(text, p['protein_molecule_types'])
        (run/'system.top').write_text(text); shutil.copy2(p['_resolved']['coordinates'], run/'solute.gro')
        record['protein_backbone_restraints_per_molecule'] = restraints
        command('solvate', ['solvate','-cp','solute.gro','-cs','spc216.gro','-o','solvated.gro','-p','system.top'])
        summary = topology_summary(text)
        water = filter_core_water(run/'solvated.gro',run/'hydrated.gro',p['water_policy'],summary['atoms'])
        # One source seed water was already present in the admitted solute.
        initial_waters = sum(n for name,n in summary['sequence'] if name == 'SOL')
        (run/'system.top').write_text(set_water_count(text, water['added_retained_waters']+initial_waters))
        record['solvation'] = water; record['state'] = 'solvated'; save()
        # This cutoff TPR is used only to select replacement waters, never for dynamics.
        (run/'ions.mdp').write_text('integrator = steep\nnsteps = 0\ncutoff-scheme = Verlet\ncoulombtype = Cut-off\nrcoulomb = 1.2\nrvdw = 1.2\npbc = xyz\n')
        command('ions-input', ['grompp','-f','ions.mdp','-c','hydrated.gro','-p','system.top','-o','ions.tpr'])
        command('neutralize', ['genion','-s','ions.tpr','-o','neutralized.gro','-p','system.top',
                              '-pname','SOD','-nname','CLA','-neutral','-conc',str(p.get('salt_molar',.15)),
                              '-seed','20261026'], 'SOL\n')
        neutral = topology_summary((run/'system.top').read_text())
        if abs(neutral['charge_e']) > 1e-4: raise ValueError('Final system is not neutral')
        if neutral['atoms'] != int((run/'neutralized.gro').read_text().splitlines()[1]):
            raise ValueError('Neutralized coordinate/topology count mismatch')
        record.update(state='neutralized', final_atoms=neutral['atoms'], charge_e=neutral['charge_e'],
                      molecule_counts=dict(Counter({name:sum(n for k,n in neutral['sequence'] if k==name)
                                                   for name,_ in neutral['sequence']})))
        save()
        record['state'] = 'minimizing'; save()
        (run/'cleanup.mdp').write_text(CLEANUP_MDP.replace('nsteps = 250',f"nsteps = {p['minimization_protocol']['contact_iterations']}"))
        command('cleanup-input', ['grompp','-f','cleanup.mdp','-c','neutralized.gro','-r','neutralized.gro',
                                 '-p','system.top','-o','cleanup.tpr','-po','cleanup-resolved.mdp'])
        command('cleanup', ['mdrun','-deffnm','cleanup','-ntmpi','1','-ntomp',str(threads),'-pin','off','-v'])
        record['flexible_cleanup'] = read_minimization(run,'cleanup'); save()
        if record['flexible_cleanup']['constraint_or_fatal_messages']:
            raise ValueError('Cleanup has unresolved constraint/fatal messages')
        record['minimization'] = final_minimization_blocks(run,record,command,threads,
                                    p['minimization_protocol'].get('max_final_blocks',30))
        result = record['minimization']
        record['state'] = ('minimized' if result['converged'] else 'minimization_not_converged')
        if result['constraint_or_fatal_messages']: record['state'] = 'minimization_needs_review'
        if record['state']=='minimized':
            from .md_geometry import check_geometry
            geometry=check_geometry(run)
            (run/'FINAL_GEOMETRY_CHECK.json').write_text(json.dumps(geometry,indent=2)+'\n')
            record['geometry_checks']=geometry
            if not geometry['passed']:record['state']='minimization_needs_review'
        record['equilibration'] = 'not_started_user_requested_pause'
        record['outputs_sha256'] = {n:sha256(run/n) for n in ['system.top','neutralized.gro','em.mdp','em.tpr','em.gro','em.g96','em.log']}
        save()
        if record['state']=='minimized':
            from .run_bundle import write_manifest
            from .molecular_scene import build_display
            artifacts={name:name for name in ['em.gro','system.top','FINAL_GEOMETRY_CHECK.json','PREPARATION_RECEIPT.json']}
            artifacts['minimization_log']='em.log'
            artifacts['coordinates']=artifacts.pop('em.gro');artifacts['topology']=artifacts.pop('system.top')
            try:
                build_display(run,p)
                artifacts.update(scene='scene.json',atom_map='ATOM_MAP.tsv')
            except ValueError as exc:
                record['display_hold']=str(exc);save()
            write_manifest(run,stage='minimized',artifacts=artifacts,units='nm',profile_id=p['id'],
                organism=p['organism'],strain=p['strain'],membrane_geometry={'kind':'plane','periodic':True},
                sources=p['sources'],assumptions=p['assumptions'],charge_e=record['charge_e'],
                minimization=result,geometry_checks=record['geometry_checks'])
        return record
    except Exception as exc:
        if record['state'] != 'failed': record.update(state='failed', error=str(exc)); save()
        raise
