"""Failure controls for the new runner; no scientific calculation or network."""
import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('ragab_replay', Path(__file__).with_name('replay.py'))
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class FailureControls(unittest.TestCase):
    def worker(self, text, outputs=()):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            replay.checked_worker([sys.executable, '-c', text], base, base / 'worker.log',
                                  [base / name for name in outputs])

    def test_nonzero_worker_refused(self):
        with self.assertRaisesRegex(RuntimeError, 'Worker failed'):
            self.worker('raise SystemExit(7)')

    def test_missing_output_refused(self):
        with self.assertRaisesRegex(RuntimeError, 'Required output missing'):
            self.worker('pass', ['summary.json'])

    def test_nonfinite_json_refused(self):
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            self.worker("from pathlib import Path; Path('summary.json').write_text('{\"area\":NaN}')", ['summary.json'])

    def test_empty_json_refused(self):
        with self.assertRaises(ValueError):
            self.worker("from pathlib import Path; Path('summary.json').write_text('')", ['summary.json'])

    def test_zero_exit_complete_output_accepted(self):
        self.worker("from pathlib import Path; Path('summary.json').write_text('{\"area\":1.0}')", ['summary.json'])

    def test_existing_output_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(RuntimeError, 'already exists'):
                replay.execute(Path(folder), Path(folder), 'missing')

    def test_existing_symlink_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'link'
            path.symlink_to(Path(folder) / 'absent')
            with self.assertRaisesRegex(RuntimeError, 'already exists'):
                replay.execute(path, Path(folder), 'missing')

    def test_checksum_tamper_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'member').write_text('altered')
            (root / 'PAYLOAD.json').write_text('{"files":{"member":"' + '0' * 64 + '"}}')
            old = replay.ROOT
            try:
                replay.ROOT = root
                with self.assertRaisesRegex(RuntimeError, 'Kit member changed'):
                    replay.verify_payload()
            finally:
                replay.ROOT = old

    def test_missing_source_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(FileNotFoundError):
                replay.verify_sources(Path(folder))

    def test_comparison_ceilings(self):
        policy = replay.read(replay.ROOT / 'EXPECTED.json')['policy']
        self.assertTrue(replay.differences({'call': 'open_like'}, {'call': 'closed_like'}, policy))
        self.assertTrue(replay.differences({'x': 1.0}, {'x': math.nan}, policy))
        self.assertTrue(replay.differences({'area_A2': 1.0}, {'area_A2': 1.001}, policy))
        self.assertFalse(replay.differences({'area_A2': 1.0}, {'area_A2': 1.00001}, policy))
        self.assertTrue(replay.differences({'count': 1}, {'count': True}, policy))


if __name__ == '__main__':
    unittest.main()
