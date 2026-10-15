"""Fixtures for the local application tests.

Nothing here touches a real protein, a real workspace or a scientific engine.
The job manager is exercised against a stand-in store whose `run` is a gate the
test opens when it wants the worker to finish, so queued and running states can
be observed deliberately instead of raced against.
"""
from __future__ import annotations

import json
import threading
import time
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path

import pytest


def wait_for(predicate, timeout: float = 5.0, interval: float = 0.01):
    """Wait for a background worker to reach a state, or fail with a clear message.

    Used instead of a bare sleep so a slow machine does not turn into a flaky
    test and a genuine hang does not turn into a silent pass.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = predicate()
        if value:
            return value
        time.sleep(interval)
    raise AssertionError(f"condition not reached within {timeout}s")


class FakeProcess:
    """Stands in for a launched engine process, carrying only its argument list."""

    def __init__(self, args):
        self.args = list(args)


class FakeStore:
    """A stand-in for the analysis store, with no engine and no real case.

    `run` blocks on `gate` until a test opens it. That keeps a submitted job in
    the running state for as long as the test needs, and lets a second job sit
    genuinely queued behind it.
    """

    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.manifests: dict[str, dict] = {}
        self.gate = threading.Event()
        self.gate.set()
        self.result = {"status": "succeeded", "run_id": "run-synthetic", "error": ""}
        self.process_args: list[str] | None = None
        self.run_calls: list[str] = []
        # Set to a digest to model a case edited between submission and execution.
        # The submit path reads the manifest once; the worker reads it again, so
        # changing the answer from the second read on needs no timing games.
        self.revision_after_first_load: str | None = None
        self._loads = 0

    def add_case(self, analysis_id: str, analysis_type: str = "structure_quality",
                 revision: str = "a" * 64) -> dict:
        self.manifests[analysis_id] = {
            "analysis_id": analysis_id,
            "analysis_type": analysis_type,
            "revision_sha256": revision,
        }
        return self.manifests[analysis_id]

    def load(self, analysis_id: str) -> dict:
        if analysis_id not in self.manifests:
            raise KeyError(analysis_id)
        self._loads += 1
        manifest = deepcopy(self.manifests[analysis_id])
        if self.revision_after_first_load is not None and self._loads > 1:
            manifest["revision_sha256"] = self.revision_after_first_load
        return manifest

    @contextmanager
    def _transaction(self):
        yield

    def run(self, analysis_id, *, cancel_event=None, on_process=None, on_progress=None):
        self.run_calls.append(analysis_id)
        if on_process is not None and self.process_args is not None:
            on_process(FakeProcess(self.process_args))
        if not self.gate.wait(timeout=5):
            raise AssertionError("store.run gate was never opened")
        return dict(self.result)


@pytest.fixture
def store(tmp_path: Path) -> FakeStore:
    return FakeStore(tmp_path / "workspace")


@pytest.fixture
def manager(store: FakeStore):
    """A running job manager that is always shut down, even when a test fails."""
    from yauvi_structural_workbench.jobs import JobManager

    instance = JobManager(store)
    try:
        yield instance
    finally:
        store.gate.set()
        instance.close()


def write_job_record(store: FakeStore, job_id: str, state: str) -> Path:
    """Leave a job file behind as a stopped server would have left it."""
    jobs = store.root / "jobs"
    jobs.mkdir(parents=True, exist_ok=True)
    path = jobs / f"{job_id}.json"
    path.write_text(json.dumps({
        "schema_version": "1.0",
        "job_id": job_id,
        "analysis_id": "case-earlier",
        "input_revision_sha256": "b" * 64,
        "state": state,
        "run_id": None,
        "error": "",
        "created_at": "2026-01-01T00:00:00+00:00",
        "started_at": None,
        "finished_at": None,
        "events": [],
    }), encoding="utf-8")
    return path
