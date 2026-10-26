"""Verified legacy MD views; portable callers supply explicit run/display roots."""
import json
import os
import re
from pathlib import Path

from .md_preparation import read_minimization, sha256

RUN = Path(os.environ['MEMORIENT_LEGACY_RUN']) if os.environ.get('MEMORIENT_LEGACY_RUN') else None
DISPLAY = Path(os.environ['MEMORIENT_LEGACY_DISPLAY']) if os.environ.get('MEMORIENT_LEGACY_DISPLAY') else None


def preparation_status():
    if RUN is None or DISPLAY is None:
        return {'state':'not_prepared','display_ready':False,'equilibration':'not_started'}
    path = RUN/'MINIMIZATION_CONTINUATION.json'
    if not path.exists(): path = RUN/'PREPARATION_RECEIPT.json'
    if not path.exists():
        return {'state':'not_prepared','display_ready':False,'equilibration':'not_started'}
    record = json.loads(path.read_text())
    result = read_minimization(RUN) if (RUN/'em.log').exists() else {}
    ready = False
    if (DISPLAY/'VIEW_RECEIPT.json').is_file():
        try:
            minimized_scene(); ready = True
        except (ValueError,OSError,KeyError):
            pass
    state = record['state']
    stage = 'preparation'
    iteration = None
    if state == 'minimizing':
        state = 'minimizing'
        filename = 'em.log' if (RUN/'em.log').exists() else 'cleanup.log'
        stage = 'final rigid-water minimum' if filename == 'em.log' else 'flexible-water contact cleanup'
        if (RUN/filename).exists():
            steps = re.findall(r'Step\s+Time\s*\n\s*(\d+)\s+', (RUN/filename).read_text())
            if steps: iteration = int(steps[-1])
    return {'state':state, 'full_atoms':record.get('final_atoms'),
            'stage':stage, 'minimizer_iteration':iteration,
            'molecule_counts':record.get('molecule_counts',{}), 'charge_e':record.get('charge_e'),
            'minimization':result, 'display_ready':ready, 'equilibration':'not_started_user_requested_pause',
            'force_field':record['force_field'], 'profile_id':record['id'],
            'assumptions':record['assumptions'], 'sources':record['sources']}


def minimized_scene(*,run=None,display=None):
    run=Path(run) if run is not None else RUN
    display=Path(display) if display is not None else DISPLAY
    if run is None or display is None:raise ValueError('Select an explicit result manifest or legacy run/display directories')
    return _minimized_scene(run,display)


def _minimized_scene(RUN,DISPLAY):
    receipt = json.loads((DISPLAY/'VIEW_RECEIPT.json').read_text())
    if receipt['state'] != 'minimized' or not receipt['minimization']['converged']:
        raise ValueError('Minimized display is not accepted')
    if receipt['minimization'].get('constraint_or_fatal_messages'):
        raise ValueError('Minimization has unresolved constraint messages')
    if sha256(RUN/'em.gro') != receipt['coordinate_sha256']:
        raise ValueError('Minimized coordinate identity mismatch')
    if sha256(RUN/'system.top') != receipt['topology_sha256']:
        raise ValueError('Minimized topology identity mismatch')
    if sha256(DISPLAY/'scene.json') != receipt['scene_sha256']:
        raise ValueError('Minimized display identity mismatch')
    if read_minimization(RUN) != receipt['minimization']:
        raise ValueError('Minimization log/result identity mismatch')
    return json.loads((DISPLAY/'scene.json').read_text())


def minimized_coordinates():
    minimized_scene()
    return (RUN/'em.gro').read_bytes()


def minimized_topology():
    minimized_scene()
    return (RUN/'system.top').read_bytes()


def minimized_atom_map():
    minimized_scene()
    receipt=json.loads((DISPLAY/'VIEW_RECEIPT.json').read_text())
    if sha256(DISPLAY/'ATOM_MAP.tsv') != receipt['atom_map_sha256']:
        raise ValueError('Minimized atom-map identity mismatch')
    return (DISPLAY/'ATOM_MAP.tsv').read_bytes()
