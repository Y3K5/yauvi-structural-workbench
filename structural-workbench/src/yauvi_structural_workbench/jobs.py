"""One bounded worker with durable jobs and explicit interruption recovery."""
from __future__ import annotations
from datetime import datetime, timezone
from copy import deepcopy
import fcntl
import json
import queue
import secrets
import threading
from pathlib import Path
from yauvi_platform.structural_workbench import AnalysisError
from yauvi_platform.structural_workbench.store import _write_json

def _now():
    return datetime.now(timezone.utc).isoformat()

class JobManager:
    def __init__(self, store, capacity=8):
        self.store = store
        self.root = store.root / 'jobs'; self.root.mkdir(exist_ok=True)
        self._lease = (self.root / '.server.lock').open('a')
        try: fcntl.flock(self._lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            self._lease.close()
            raise AnalysisError('a workbench server already owns this workspace') from exc
        self._lock = threading.RLock(); self._stop = threading.Event()
        self._queue = queue.Queue(maxsize=capacity); self._events = {}; self._records = {}
        for path in sorted(self.root.glob('job-*.json')):
            row = json.loads(path.read_text())
            if row['state'] in {'queued','running'}:
                row.update(state='interrupted', finished_at=_now(), error='Previous server stopped; resubmit deliberately.')
                row.setdefault('events', []).append({'at': _now(), 'stage': 'interrupted', 'message': 'Server restarted; previous work was interrupted'})
                _write_json(path,row)
            self._records[row['job_id']] = row
        self._worker = threading.Thread(target=self._work, name='yauvi-analysis-worker', daemon=True)
        self._worker.start()

    def _save(self,row): _write_json(self.root / (row['job_id']+'.json'),row)

    def _progress(self, row, stage, message):
        with self._lock:
            row['updated_at'] = _now()
            row['stage'] = stage
            row.setdefault('events', []).append({'at': row['updated_at'], 'stage': stage, 'message': message})
            self._save(row)

    def list(self):
        with self._lock:
            rows = deepcopy(sorted(self._records.values(), key=lambda row: row.get("created_at", "")))
            position = 0
            for row in rows:
                if row['state'] == 'queued':
                    position += 1
                    row['queue_position'] = position
            return rows

    def submit(self,analysis_id, *, experimental_acknowledged=False):
        manifest=self.store.load(analysis_id)
        if manifest['analysis_type']=='membrane_orientation' and not experimental_acknowledged:
            raise AnalysisError('membrane orientation is experimental; acknowledge its limits before running')
        with self._lock:
            if self._stop.is_set(): raise AnalysisError('server is stopping')
            if any(r['analysis_id']==analysis_id and r['state'] in {'queued','running'} for r in self._records.values()):
                raise AnalysisError('this analysis already has a queued or running job')
            if self._queue.full(): raise AnalysisError('analysis queue is full; retry after a job finishes')
            job_id='job-'+secrets.token_hex(12)
            row={'schema_version':'1.0','job_id':job_id,'analysis_id':analysis_id,
                 'input_revision_sha256':manifest['revision_sha256'],'state':'queued','run_id':None,'error':'',
                 'created_at':_now(),'started_at':None,'finished_at':None,'events':[]}
            self._records[job_id]=row;self._events[job_id]=threading.Event();self._progress(row,'queued','Waiting for the local worker')
            self._queue.put_nowait(job_id)
            return dict(row)

    def cancel(self,job_id):
        with self._lock:
            if job_id not in self._records: raise AnalysisError('unknown job')
            row=self._records[job_id]
            if row['state'] in {'queued','running'}:
                self._events[job_id].set()
                if row['state']=='queued':row.update(state='cancelled',finished_at=_now())
                row['cancel_requested'] = True
                self._progress(row,'cancelling','Cancellation requested; stopping owned processes')
            return dict(row)

    def _work(self):
        try:
            while not self._stop.is_set():
                try: job_id=self._queue.get(timeout=.1)
                except queue.Empty: continue
                with self._lock:
                    row=self._records[job_id];event=self._events[job_id]
                    if event.is_set(): self._queue.task_done();continue
                    row.update(state='running',started_at=_now())
                    self._progress(row,'starting','Worker started')
                try:
                    with self.store._transaction():
                        if self.store.load(row['analysis_id'])['revision_sha256']!=row['input_revision_sha256']:
                            raise AnalysisError('case changed after submission; review and resubmit')
                        def process_started(process):
                            # Expose only the registered module name, never command arguments or paths.
                            args = process.args
                            module = args[args.index('-m') + 1] if '-m' in args else 'scientific engine'
                            self._progress(row, 'engine', 'Running ' + module)
                        result=self.store.run(row['analysis_id'],cancel_event=event,
                            on_process=process_started,
                            on_progress=lambda stage, message: self._progress(row, stage, message))
                    with self._lock:
                        row.update(state=result['status'],run_id=result['run_id'],error=result.get('error',''))
                except Exception as exc:
                    with self._lock:row.update(state='failed',error=f'{type(exc).__name__}: {exc}')
                finally:
                    with self._lock:
                        row['finished_at'] = _now()
                        self._progress(row,row['state'],'Execution ended: ' + row['state'].replace('_',' '))
                    self._queue.task_done()
        finally:
            with self._lock:
                for row in self._records.values():
                    if row['state'] in {'queued','running'}:
                        row.update(state='interrupted',finished_at=_now(),error='Server stopped; resubmit deliberately.');self._save(row)
            fcntl.flock(self._lease,fcntl.LOCK_UN);self._lease.close()

    def close(self):
        self._stop.set()
        with self._lock:
            for event in self._events.values():event.set()
        self._worker.join(timeout=7)
