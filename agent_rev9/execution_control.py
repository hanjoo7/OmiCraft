"""Supervise web jobs in separate processes so cancellation stops their compute tree."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time

import psutil

cancel_event = threading.Event()


def stop_tree(process):
    """Freeze parents before enumerating children, including separate process sessions."""
    pending = [psutil.Process(process.pid)] if psutil.pid_exists(process.pid) else []
    frozen = []
    while pending:
        item = pending.pop()
        try:
            item.suspend()
            frozen.append(item)
            pending.extend(item.children())
        except psutil.NoSuchProcess:
            pass
    for item in reversed(frozen):
        try:
            item.terminate()
            item.resume()
        except psutil.NoSuchProcess:
            pass
    _, alive = psutil.wait_procs(frozen, timeout=3)
    for item in alive:
        try:
            item.kill()
        except psutil.NoSuchProcess:
            pass
    psutil.wait_procs(alive, timeout=3)
    process.wait(timeout=5)


def save_interruption(server, status, message):
    """Keep partial artifacts and mark only the interrupted run and its active children."""
    root = server.REPORT_RUNS_ROOT
    identifier = server.run_state.get('report_run_id')
    if not identifier:
        return
    records = {}
    for path in root.glob('*/run_summary.json'):
        try:
            records[path.parent.name] = (path, json.loads(path.read_text()))
        except (OSError, ValueError):
            continue
    selected = {identifier}
    while True:
        children = {key for key, (_, data) in records.items()
                    if data.get('parent_run_id') in selected or data.get('batch_run_id') in selected}
        if children <= selected:
            break
        selected |= children
    for key in selected:
        if key not in records:
            continue
        path, data = records[key]
        if key != identifier and not data.get('running') and data.get('execution_status') not in {'RUNNING', 'QUEUED', 'PARTIAL', 'NOT_RUN'}:
            continue
        data.update(running=False, execution_status=status, workflow_status='FINISHED',
                    validation_decision='NOT_EVALUATED', scientific_validation_status='NOT_EVALUATED')
        data.setdefault('blockers', []).append(message)
        for stage in data.get('stages', []):
            if stage.get('execution_status') in {'RUNNING', 'QUEUED'}:
                stage.update(execution_status=status, status=status.lower())
        tmp = path.with_suffix('.tmp')
        tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2))
        tmp.replace(path)
        try:
            server.write_competition_report(data, path.parent)
        except Exception:
            # The summary remains authoritative even if partial artifacts cannot render.
            pass


def start_worker(function, *args, **kwargs):
    from . import server
    from .configuration import get_config
    with server.run_lock:
        initial = copy.deepcopy(server.run_state)
    def encode(value):
        return {'__omicraft_config__': value.model_dump()} if hasattr(value, 'model_dump') else value
    job = {'module': function.__module__.removeprefix(__package__ + '.'), 'function': function.__name__,
           'args': [encode(value) for value in args], 'kwargs': kwargs, 'state': initial,
           'config': get_config().model_dump()}
    thread = threading.Thread(target=supervise, args=(job,), daemon=True)
    thread.start()


def supervise(job):
    from . import server
    process = None
    error = None
    cancelled = False
    last = None
    try:
        with tempfile.TemporaryDirectory(prefix='omicraft-worker-') as temporary:
            folder = Path(temporary)
            request = folder / 'input.json'
            request.write_text(json.dumps(job, ensure_ascii=False))
            request.chmod(0o600)
            with (folder / 'worker.log').open('w+') as log:
                process = subprocess.Popen([sys.executable, str(Path(__file__).parent / 'scripts/execution_worker.py'), str(request)],
                                           stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                while True:
                    snapshot = folder / 'state.json'
                    if snapshot.is_file():
                        last = json.loads(snapshot.read_text())
                        with server.run_lock:
                            # Keep the reservation until the worker actually exits.
                            server.run_state.update({key: value for key, value in last.items()
                                                     if key not in {'running', 'cancelling', 'result'}})
                    if cancel_event.is_set() and process.poll() is None:
                        cancelled = True
                        stop_tree(process)
                        break
                    if process.poll() is not None:
                        # The final atomic snapshot is written before normal worker exit.
                        if snapshot.is_file():
                            last = json.loads(snapshot.read_text())
                        if process.returncode:
                            log.seek(0)
                            error = 'Worker exited: ' + str(process.returncode) + '\n' + log.read()[-2000:]
                        break
                    time.sleep(.15)
    except Exception as exc:
        error = str(exc)
    finally:
        if process is not None and process.poll() is None:
            stop_tree(process)
        if cancelled or error:
            status = 'CANCELLED' if cancelled else 'FAILED'
            message = '사용자가 실험을 중단했습니다. 완료된 파일은 보존됩니다.' if cancelled else error
            terminal = {'type': 'done' if cancelled else 'error', 'execution_status': status,
                        'message': message, 'errors': [] if cancelled else [message],
                        'blockers': [message], 'scientific_validation_status': 'NOT_EVALUATED',
                        'report_run_id': server.run_state.get('report_run_id'),
                        'total_time': round(time.time() - job['state']['started_at'], 1)}
            try:
                save_interruption(server, status, message)
                if job['function'] == 'run_lab_docking':
                    from .web_execution import lab_directory
                    (lab_directory(job['args'][0]) / 'status.json').write_text(json.dumps(
                        {'session': job['args'][0], 'status': 'error', 'execution_status': status, 'error': message}))
            except Exception as exc:
                terminal.setdefault('errors', []).append('Partial report: ' + str(exc))
            from .workflow_monitor import finish_interrupted
            finish_interrupted(server, terminal)
            with server.run_lock:
                server.run_state['history'].append(terminal)
        else:
            terminal = (last or {}).get('result') or {'type': 'error', 'execution_status': 'FAILED', 'message': 'Worker returned no result'}
            with server.run_lock:
                if last:
                    server.run_state.update(last)
        with server.run_lock:
            server.run_state.update(running=False, cancelling=False, result=terminal)
            server.event_queue.append(terminal)
