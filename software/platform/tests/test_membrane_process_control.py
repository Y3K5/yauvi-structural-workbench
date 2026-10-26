"""Cancel a real nested local process tree without a molecular runtime."""
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
import pytest
from yauvi_platform.structural_workbench.store import StructuralAnalysisStore


@pytest.mark.skipif(os.name!='posix',reason='Preparation process control supports macOS/Linux first')
def test_nested_preparation_worker_is_cancelled(tmp_path):
    store=StructuralAnalysisStore(tmp_path/'workspace');cancel=threading.Event();pid_file=tmp_path/'child.pid'
    child='import time; time.sleep(30)'
    script=('import subprocess,sys,time; from pathlib import Path; '
            f'p=subprocess.Popen([sys.executable,"-c",{child!r}]); '
            f'Path({str(pid_file)!r}).write_text(str(p.pid)); time.sleep(30)')
    def request_cancel():
        deadline=time.monotonic()+3
        while not pid_file.is_file() and time.monotonic()<deadline:time.sleep(.01)
        cancel.set()
    thread=threading.Thread(target=request_cancel);thread.start()
    try:
        with pytest.raises(InterruptedError):
            store._execute([sys.executable,'-c',script],cwd=tmp_path,log_path=tmp_path/'worker.log',cancel_event=cancel)
    finally:cancel.set();thread.join(timeout=4)
    assert pid_file.is_file()
    pid=int(pid_file.read_text())
    state=subprocess.run(['ps','-p',str(pid),'-o','stat='],capture_output=True,text=True,check=False).stdout.strip()
    assert not state or state.startswith('Z'),state
    assert (tmp_path/'worker.log').is_file()
