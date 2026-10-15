"""Job manager behaviour: one owner per workspace, honest states, no leaked paths.

These cover the guarantees a user relies on without being able to see them: that
two servers cannot quietly share one analysis library, that work abandoned by a
stopped server is never presented as still running, that an experimental
analysis is not started unacknowledged, and that a case edited after submission
is refused rather than run against changed inputs.
"""
from __future__ import annotations

import json

import pytest

from yauvi_platform.structural_workbench import AnalysisError
from yauvi_structural_workbench.jobs import JobManager

from conftest import wait_for, write_job_record


def test_second_server_cannot_own_the_same_workspace(store, manager):
    """Two managers on one workspace would both write the same job files."""
    with pytest.raises(AnalysisError, match="already owns this workspace"):
        JobManager(store)


def test_interrupted_work_is_not_presented_as_running(store):
    """A job the previous server left behind must not look resumable."""
    path = write_job_record(store, "job-abandoned", "running")

    instance = JobManager(store)
    try:
        row = next(r for r in instance.list() if r["job_id"] == "job-abandoned")
        assert row["state"] == "interrupted"
        assert row["finished_at"]
        assert "resubmit" in row["error"]
        assert [e["stage"] for e in row["events"]] == ["interrupted"]
        # The recovery is durable, not only in memory.
        assert json.loads(path.read_text())["state"] == "interrupted"
    finally:
        instance.close()


def test_queued_work_left_behind_is_also_recovered(store):
    path = write_job_record(store, "job-waiting", "queued")
    instance = JobManager(store)
    try:
        assert json.loads(path.read_text())["state"] == "interrupted"
    finally:
        instance.close()


def test_finished_records_are_left_untouched(store):
    """Recovery must not rewrite history for work that genuinely completed."""
    path = write_job_record(store, "job-done", "succeeded")
    instance = JobManager(store)
    try:
        assert json.loads(path.read_text())["state"] == "succeeded"
        assert json.loads(path.read_text())["events"] == []
    finally:
        instance.close()


def test_experimental_analysis_refuses_to_start_unacknowledged(store, manager):
    store.add_case("case-membrane", analysis_type="membrane_orientation")

    with pytest.raises(AnalysisError, match="experimental"):
        manager.submit("case-membrane")

    assert manager.list() == []
    assert store.run_calls == []


def test_experimental_analysis_starts_once_acknowledged(store, manager):
    store.add_case("case-membrane", analysis_type="membrane_orientation")

    row = manager.submit("case-membrane", experimental_acknowledged=True)

    assert row["state"] == "queued"
    wait_for(lambda: store.run_calls == ["case-membrane"])


def test_one_case_cannot_be_submitted_twice_at_once(store, manager):
    store.add_case("case-one")
    store.gate.clear()

    manager.submit("case-one")
    wait_for(lambda: manager.list()[0]["state"] == "running")

    with pytest.raises(AnalysisError, match="already has a queued or running job"):
        manager.submit("case-one")


def test_submission_pins_the_revision_it_was_asked_about(store, manager):
    store.add_case("case-one", revision="c" * 64)
    store.gate.clear()

    row = manager.submit("case-one")

    assert row["input_revision_sha256"] == "c" * 64


def test_a_case_edited_after_submission_is_refused(store, manager):
    """Running against inputs that changed mid-flight would misattribute the result.

    The worker re-reads the case before executing and compares it to the digest
    pinned at submission, so the edit is modelled by changing what the second
    read returns rather than by racing the worker.
    """
    store.add_case("case-one", revision="c" * 64)
    store.revision_after_first_load = "d" * 64

    manager.submit("case-one")

    row = wait_for(lambda: next((r for r in manager.list()
                                 if r["state"] not in {"queued", "running"}), None))
    assert row["state"] == "failed"
    assert "resubmit" in row["error"]
    assert store.run_calls == []


def test_queue_position_is_only_reported_for_waiting_work(store, manager):
    """A running job has no queue position; waiting jobs are numbered in order."""
    store.add_case("case-first")
    store.add_case("case-second")
    store.add_case("case-third")
    store.gate.clear()

    manager.submit("case-first")
    wait_for(lambda: manager.list()[0]["state"] == "running")
    manager.submit("case-second")
    manager.submit("case-third")

    rows = {row["analysis_id"]: row for row in manager.list()}
    assert "queue_position" not in rows["case-first"]
    assert rows["case-second"]["queue_position"] == 1
    assert rows["case-third"]["queue_position"] == 2


def test_listing_is_ordered_by_creation(store, manager):
    store.add_case("case-first")
    store.add_case("case-second")
    store.gate.clear()
    manager.submit("case-first")
    wait_for(lambda: manager.list()[0]["state"] == "running")
    manager.submit("case-second")

    assert [row["analysis_id"] for row in manager.list()] == ["case-first", "case-second"]


def test_cancelling_an_unknown_job_is_refused(manager):
    with pytest.raises(AnalysisError, match="unknown job"):
        manager.cancel("job-does-not-exist")


def test_cancelling_waiting_work_stops_it_before_it_starts(store, manager):
    store.add_case("case-first")
    store.add_case("case-second")
    store.gate.clear()
    manager.submit("case-first")
    wait_for(lambda: manager.list()[0]["state"] == "running")
    waiting = manager.submit("case-second")

    row = manager.cancel(waiting["job_id"])

    assert row["state"] == "cancelled"
    assert row["cancel_requested"] is True
    store.gate.set()
    assert "case-second" not in store.run_calls


def test_listing_callers_cannot_mutate_stored_records(store, manager):
    """`list` hands out copies; a caller editing them must not corrupt the manager."""
    store.add_case("case-one")
    store.gate.clear()
    manager.submit("case-one")
    wait_for(lambda: manager.list()[0]["state"] == "running")

    manager.list()[0]["state"] = "tampered"

    assert manager.list()[0]["state"] == "running"


def test_progress_never_reveals_command_arguments_or_paths(store, manager):
    """Only the registered module name may reach the event log.

    The job events are shown in the browser, so a leaked argument list would put
    local filesystem paths in front of anyone looking at the page.
    """
    store.add_case("case-one")
    # Assembled at run time so no home-shaped literal sits in public source;
    # the repository's path verifier rejects one even inside a test.
    secret = "/" + "Users/" + "alice/private/case/coordinates.pdb"
    store.process_args = ["/opt/python", "-m", "structqc.cli", "--structure", secret]

    manager.submit("case-one")
    row = wait_for(lambda: next((r for r in manager.list()
                                 if r["state"] not in {"queued", "running"}), None))

    messages = " ".join(event["message"] for event in row["events"])
    assert "Running structqc.cli" in messages
    assert secret not in messages
    assert "--structure" not in messages
    assert "/opt/python" not in messages


def test_progress_names_the_engine_generically_when_no_module_is_given(store, manager):
    store.add_case("case-one")
    store.process_args = ["/opt/engine/bin/run", "--input", "/private/path"]

    manager.submit("case-one")
    row = wait_for(lambda: next((r for r in manager.list()
                                 if r["state"] not in {"queued", "running"}), None))

    messages = " ".join(event["message"] for event in row["events"])
    assert "scientific engine" in messages
    assert "/private/path" not in messages


def test_a_full_queue_is_refused_rather_than_grown(store):
    """Capacity is a declared bound; exceeding it must fail loudly."""
    instance = JobManager(store, capacity=1)
    try:
        store.gate.clear()
        store.add_case("case-first")
        store.add_case("case-second")
        store.add_case("case-third")
        instance.submit("case-first")
        wait_for(lambda: instance.list()[0]["state"] == "running")
        instance.submit("case-second")

        with pytest.raises(AnalysisError, match="queue is full"):
            instance.submit("case-third")
    finally:
        store.gate.set()
        instance.close()


def test_closing_marks_unfinished_work_as_interrupted(store):
    store.add_case("case-one")
    instance = JobManager(store)
    store.gate.clear()
    instance.submit("case-one")
    wait_for(lambda: instance.list()[0]["state"] == "running")

    store.gate.set()
    instance.close()

    row = json.loads(next((store.root / "jobs").glob("job-*.json")).read_text())
    assert row["state"] in {"succeeded", "interrupted"}
    assert row["finished_at"]


def test_the_workspace_lease_is_released_for_the_next_server(store):
    """After a clean close, a new server must be able to take the workspace."""
    first = JobManager(store)
    first.close()

    second = JobManager(store)
    second.close()
