"""Cancellation tests use only a synthetic sleeping worker, never scientific tools."""
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil
import pytest
from fastapi.testclient import TestClient

from .. import execution_control as control, server
from ..configuration import get_config, configuration_scope


def synthetic_job(mode, output):
    from ..run_records import RunRecord
    directory = Path(get_config().data.base_dir) / 'upstream_20260929T000000Z_12345678'
    record = RunRecord(directory, 'UPSTREAM')
    record.data.update(running=True, execution_profile='upstream_analysis')
    record.event('synthetic', 'RUNNING')
    server.run_state['report_run_id'] = directory.name
    if mode == 'sleep':
        child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(90)'], start_new_session=True)
        Path(output).write_text(str(child.pid))
        child.wait()
    terminal = {'type': 'done', 'execution_status': 'COMPLETED', 'report_run_id': directory.name}
    server.run_state.update(running=False, result=terminal)
    server.run_state['history'].append(terminal)


def wait_for(predicate, timeout=20):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(.05)
    pytest.fail('Timed out waiting for worker')


@pytest.mark.external_tool
def test_cancel_stops_separate_session_child_preserves_report_and_allows_retry(tmp_path, monkeypatch):
    monkeypatch.setattr(server, 'run_state', dict(running=False, run_id=0, result=None, history=[]))
    monkeypatch.setattr(server, 'REPORT_RUNS_ROOT', tmp_path)
    config = get_config().model_copy(deep=True)
    config.data.base_dir = str(tmp_path)
    output = tmp_path / 'child.pid'
    with TestClient(server.app) as client, configuration_scope(config):
        first = server.begin_pipeline_run('synthetic cancellation')
        control.start_worker(synthetic_job, 'sleep', str(output))
        try:
            wait_for(output.is_file)
            pid = int(output.read_text())
            assert psutil.pid_exists(pid)
            assert client.post('/api/stop', json={'run_id': first + 1}).status_code == 409
            assert not control.cancel_event.is_set()
            assert client.post('/api/stop', json={'run_id': first}).json()['status'] == 'stopping'
            assert server.begin_pipeline_run('must remain reserved') is None
            wait_for(lambda: not server.run_state['running'])
            assert not psutil.pid_exists(pid) or psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
            assert server.run_state['result']['execution_status'] == 'CANCELLED'
            summary = tmp_path / server.run_state['report_run_id'] / 'run_summary.json'
            assert json.loads(summary.read_text())['execution_status'] == 'CANCELLED'
            assert client.post('/api/stop', json={'run_id': first}).json()['status'] == 'not_running'
            # A fresh directory for the retry keeps the first partial report intact.
            config.data.base_dir = str(tmp_path / 'retry')
            assert server.begin_pipeline_run('retry') == first + 1
            assert not control.cancel_event.is_set()
            control.start_worker(synthetic_job, 'done', str(output))
            wait_for(lambda: not server.run_state['running'])
            assert server.run_state['result']['execution_status'] == 'COMPLETED'
        finally:
            if server.run_state['running']:
                control.cancel_event.set()
                wait_for(lambda: not server.run_state['running'])
