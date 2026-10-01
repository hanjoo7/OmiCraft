"""Bound design jobs to an isolated process and persist a report on every exit."""
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from .small_molecule_io import write_json


def stop_process_group(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def interrupted_report(directory, body, reason):
    from . import server
    path = directory / 'run_summary.json'
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        data = {'run_id':directory.name, 'modality':body['modality'], 'stages':[], 'artifacts':[], 'is_mock':False}
    for stage in data.get('stages', []):
        if stage.get('execution_status') in {'RUNNING','QUEUED'}:
            stage.update(execution_status='INTERRUPTED', status='interrupted')
    data.update(execution_profile='therapeutic_design', gene=body['gene'],
                parent_run_id=body['parent_run_id'], batch_run_id=body.get('batch_run_id'),
                execution_status='FAILED', running=False, validation_decision='NOT_EVALUATED',
                blockers=list(dict.fromkeys(data.get('blockers', [])+[reason])),
                critic={'decision':'NOT_EVALUATED'}, workflow_status='FINISHED', report_status='READY')
    if data.get('modality') == 'DEGRADER':
        summary = data.get('summary', {})
        prediction = summary.get('ternary_prediction', {})
        if prediction.get('status') == 'running':
            prediction.update(status='interrupted', error_message=reason)
            summary['ternary_status'] = 'INTERRUPTED'
    registered = {a['path'] for a in data.get('artifacts',[]) if isinstance(a,dict)}
    for file in directory.rglob('*'):
        if file.is_file() and file.suffix in {'.log','.csv','.sdf','.cif','.pdb','.png'} and str(file.resolve()) not in registered:
            data['artifacts'].append({'path':str(file.resolve()),'label':file.name})
    write_json(path,data)
    server.write_competition_report(data,directory)
    return {'type':'done','execution_status':'FAILED','scientific_validation_status':'NOT_EVALUATED',
            'report_run_id':directory.name,'errors':[reason]}



def execution_budget(body, config):
    default = 2100
    if body.get('modality') == 'DEGRADER' and config.get('run_ternary'):
        from .small_molecule_config import LigandAF3Config
        model = LigandAF3Config(**config.get('af3', {}))
        default = model.data_pipeline_timeout_seconds + (model.inference_timeout_seconds or model.timeout_seconds) + 600
    budget = float(config.get('execution_timeout_seconds', default))
    if not math.isfinite(budget) or budget <= 0:
        raise ValueError('execution_timeout_seconds must be finite and positive')
    return min(14400, budget)


def run_bounded_design(identifier, body, data, config, qualification, app_config, *, finalize=False):
    from . import server
    directory = Path(app_config.data.base_dir) / identifier
    job = directory/'worker_input.json'
    start = time.monotonic()
    events = directory/'worker_events.jsonl'
    cursor = 0
    terminal = None
    process = None
    try:
        directory.mkdir(parents=True, exist_ok=True)
        request = {'identifier':identifier,'body':body,'data':data,'config':config,
                   'qualification':qualification,'app_config':app_config.model_dump()}
        write_json(job,request)
        os.chmod(job,0o600)
        timeout = execution_budget(body, config)
        script = Path(__file__).parent/'scripts/design_worker.py'
        with (directory/'worker.log').open('w') as log:
            process = subprocess.Popen([sys.executable,str(script),str(job)],stdout=log,stderr=subprocess.STDOUT,
                                       start_new_session=True)
            while True:
                if events.exists():
                    with events.open() as stream:
                        stream.seek(cursor)
                        while line := stream.readline():
                            if not line.endswith('\n'):
                                break
                            try:
                                event = json.loads(line)
                            except ValueError:
                                break
                            cursor = stream.tell()
                            with server.run_lock:
                                server.run_state['history'].append(event)
                                server.event_queue.append(event)
                if process.poll() is not None:
                    break
                if time.monotonic()-start >= timeout:
                    stop_process_group(process)
                    terminal = interrupted_report(directory,body,f'EXECUTION_TIMEOUT: {timeout:g}s limit reached; completed files retained')
                    break
                time.sleep(.2)
        if terminal is None:
            output = directory/'worker_result.json'
            if process.returncode==0 and output.is_file():
                terminal=json.loads(output.read_text())
            else:
                terminal=interrupted_report(directory,body,f'DESIGN_WORKER_EXIT: {process.returncode}; see worker.log')
    except Exception as exc:
        if process is not None and process.poll() is None:
            stop_process_group(process)
        terminal=interrupted_report(directory,body,str(exc))
    finally:
        if process is not None:
            stop_process_group(process)
        job.unlink(missing_ok=True)
        if terminal is not None:
            terminal['total_time']=round(time.monotonic()-start,1)
            if finalize:
                with server.run_lock:
                    server.run_state.update(running=False,result=terminal)
                    server.run_state['history'].append(terminal)
                    server.event_queue.append(terminal)
    return terminal
