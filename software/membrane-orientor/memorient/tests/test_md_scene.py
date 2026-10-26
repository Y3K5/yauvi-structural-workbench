import json

import pytest

from memorient import md_scene
from memorient.md_preparation import sha256


def test_missing_and_failed_preparations_never_look_complete(tmp_path,monkeypatch):
    run=tmp_path/'run';run.mkdir();display=tmp_path/'view';display.mkdir()
    monkeypatch.setattr(md_scene,'RUN',run);monkeypatch.setattr(md_scene,'DISPLAY',display)
    assert md_scene.preparation_status()['state']=='not_prepared'
    (run/'PREPARATION_RECEIPT.json').write_text(json.dumps(dict(state='failed',force_field='test',
        id='invented',assumptions=['test'],sources=['test'])))
    (run/'em.log').write_text('unfinished minimization')
    status=md_scene.preparation_status()
    assert status['state']=='failed' and not status['display_ready']


def test_unconverged_and_tampered_displays_are_held(tmp_path,monkeypatch):
    run=tmp_path/'run';run.mkdir();display=tmp_path/'view';display.mkdir()
    monkeypatch.setattr(md_scene,'RUN',run);monkeypatch.setattr(md_scene,'DISPLAY',display)
    receipt=dict(state='minimized',minimization={'converged':False})
    (display/'VIEW_RECEIPT.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='not accepted'):md_scene.minimized_scene()
    (run/'em.gro').write_text('invented coordinate bytes')
    receipt.update(minimization={'converged':True},coordinate_sha256=sha256(run/'em.gro'))
    (display/'VIEW_RECEIPT.json').write_text(json.dumps(receipt))
    (run/'em.gro').write_text('changed coordinate bytes')
    with pytest.raises(ValueError,match='coordinate identity'):md_scene.minimized_scene()


def test_constraint_messages_prevent_minimized_display(tmp_path,monkeypatch):
    monkeypatch.setattr(md_scene,'RUN',tmp_path)
    monkeypatch.setattr(md_scene,'DISPLAY',tmp_path)
    (tmp_path/'VIEW_RECEIPT.json').write_text(json.dumps(dict(state='minimized',
        minimization={'converged':True,'constraint_or_fatal_messages':['LINCS WARNING']})))
    with pytest.raises(ValueError,match='unresolved'):md_scene.minimized_scene()
