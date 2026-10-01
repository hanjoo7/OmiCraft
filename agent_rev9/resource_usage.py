"""Execution-local, durable accounting. No prompts, credentials or inferred tokens."""
from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps
from pathlib import Path
import json
import threading
import time
import uuid

_DIRECTORY = ContextVar('resource_directory', default=None)
_LOCK = threading.RLock()


@contextmanager
def resource_scope(directory):
    if not directory:
        yield
        return
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    token = _DIRECTORY.set(directory)
    try:
        with _LOCK:
            if not (directory / 'resource_events.jsonl').exists():
                record_usage('started')
        yield
    finally:
        _DIRECTORY.reset(token)


def record_usage(kind, **values):
    directory = _DIRECTORY.get()
    if directory is None:
        return
    with _LOCK:
        with (directory / 'resource_events.jsonl').open('a') as stream:
            stream.write(json.dumps({'kind': kind, 'timestamp': time.time(), **values}, allow_nan=False)+'\n')
        usage = read_usage(directory)
        temporary = directory / 'resource_usage.tmp'
        temporary.write_text(json.dumps(usage, indent=2))
        temporary.replace(directory / 'resource_usage.json')


def read_usage(directory):
    path = Path(directory) / 'resource_events.jsonl'
    if not path.is_file():
        return None
    requests, responses, leases = set(), {}, {}
    tools = hits = 0
    gpu_seconds = 0.0
    for line in path.read_text().splitlines():
        try:
            row = json.loads(line)
        except ValueError:
            continue  # A worker may have been stopped during its final write.
        kind = row.get('kind')
        if kind == 'llm_request':
            requests.add(row['id'])
        elif kind == 'llm_response':
            responses[row['id']] = row.get('usage') or {}
        elif kind == 'tool_call':
            tools += 1
        elif kind == 'cache_hit':
            hits += 1
        elif kind == 'gpu_start':
            leases[row['id']] = row
        elif kind == 'gpu_end' and row['id'] in leases:
            gpu_seconds += max(0, row['seconds'])
            leases.pop(row['id'])
    tokens = 0
    missing = 0
    for request in requests:
        usage = responses.get(request, {})
        total = usage.get('total_tokens')
        if total is None and all(type(usage.get(k)) is int for k in ('input_tokens', 'output_tokens')):
            total = usage['input_tokens'] + usage['output_tokens']
        if type(total) is int and total >= 0:
            tokens += total
        else:
            missing += 1
    return {'schema_version': 1, 'total_tokens': tokens if not missing else None,
            'recorded_tokens': tokens, 'llm_calls': len(requests), 'llm_usage_missing': missing,
            'tool_calls': tools, 'cache_hits': hits,
            'gpu_hours': gpu_seconds / 3600 if not leases else None,
            'gpu_recorded_hours': gpu_seconds / 3600, 'gpu_unfinished_leases': len(leases),
            'gpu_measurement': 'allocated_device_wall_time',
            'tool_measurement': 'instrumented_analysis_and_design_calls',
            'scope': 'current_execution_only'}


def measured_tool(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        if not kwargs.get('dry_run') and not kwargs.get('is_mock'):
            record_usage('tool_call', tool=function.__module__+'.'+function.__qualname__)
        return function(*args, **kwargs)
    return wrapped


def measured_design(function):
    @wraps(function)
    def wrapped(identifier, body, data, config, qualification, app_config, **kwargs):
        with resource_scope(Path(app_config.data.base_dir) / identifier):
            return function(identifier, body, data, config, qualification, app_config, **kwargs)
    return wrapped


@contextmanager
def gpu_allocation(device):
    identity = uuid.uuid4().hex
    started = time.monotonic()
    record_usage('gpu_start', id=identity, device=device)
    try:
        yield
    finally:
        record_usage('gpu_end', id=identity, device=device, seconds=time.monotonic()-started)


def report_usage(directory, children=()):
    """Aggregate only original child executions, never previous/follow-up reports."""
    directory = Path(directory).resolve()
    own = read_usage(directory)
    if own is None:
        return None
    result = dict(own)
    missing_children = []
    for child in sorted(set(children)):
        path = (directory.parent / child).resolve()
        if path.parent != directory.parent or not child.startswith('design_'):
            continue
        usage = read_usage(path)
        if usage is None:
            try:
                old = json.loads((path/'run_summary.json').read_text())
            except (OSError, ValueError):
                old = {}
            if old.get('execution_mode') != 'NOT_RUN':
                missing_children.append(child)
            continue
        for key in ('recorded_tokens','llm_calls','llm_usage_missing','tool_calls','cache_hits','gpu_recorded_hours','gpu_unfinished_leases'):
            result[key] += usage[key]
        for key in ('total_tokens','gpu_hours'):
            result[key] = result[key]+usage[key] if result[key] is not None and usage[key] is not None else None
    result['unmetered_children'] = missing_children
    if missing_children:
        result.update(total_tokens=None, gpu_hours=None, tool_calls=None)
    return result
