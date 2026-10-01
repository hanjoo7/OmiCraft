"""Read-only LangGraph callback observer. Never changes analysis inputs or outputs."""
import copy
import json
import math
import re
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    from langchain_core.callbacks import BaseCallbackHandler
except ImportError:
    class BaseCallbackHandler:
        pass


def now():
    return datetime.now(timezone.utc).isoformat()


def preview(value):
    """Bound traversal as well as serialized size; do not repr arbitrary objects."""
    budget = [6000]
    def visit(item, depth=0):
        if budget[0] <= 0:
            return '[preview limit]'
        if isinstance(item, float) and not math.isfinite(item):
            return str(item)
        if item is None or isinstance(item, (bool, int, float)):
            return item
        if isinstance(item, str):
            limit = min(800, budget[0]); budget[0] -= min(len(item), limit)
            return item[:limit] + ('…' if len(item) > limit else '')
        if depth >= 4:
            return {'type': type(item).__name__, 'summary': 'Nested value omitted'}
        if isinstance(item, dict):
            out = {}
            for index, (key, val) in enumerate(item.items()):
                if index >= 30 or budget[0] <= 0:
                    out['_preview'] = f'{len(item)} keys; remaining omitted'; break
                key = str(key)[:120]; budget[0] -= len(key)
                if re.search(r'secret|password|api.?key|authorization|access.?token|credential', key, re.I):
                    out[key] = '[redacted]'
                elif callable(val):
                    out[key] = '[callback omitted]'
                else:
                    out[key] = visit(val, depth + 1)
            return out
        if isinstance(item, (list, tuple)):
            return {'type': type(item).__name__, 'count': len(item),
                    'preview': [visit(v, depth + 1) for v in item[:8]]}
        shape = getattr(item, 'shape', None)
        return {'type': type(item).__name__, 'shape': list(shape) if isinstance(shape, tuple) else None}
    return visit(value)


def graph_definition(graph):
    if not hasattr(graph, 'get_graph'):
        return {'nodes': [], 'edges': [], 'available': False, 'reason': 'LangGraph topology unavailable'}
    raw = graph.get_graph().to_json()
    return {'available': True, 'nodes': [{'id': n['id']} for n in raw['nodes']],
            'edges': [{'source': e['source'], 'target': e['target'],
                       'conditional': bool(e.get('conditional')), 'label': str(e.get('data') or '')}
                      for e in raw['edges']]}


def persist(snapshot, directory):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / 'workflow.json'
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(snapshot, ensure_ascii=False, allow_nan=False))
    temporary.replace(path)


class WorkflowObserver(BaseCallbackHandler):
    run_inline = True

    def __init__(self, graph, run_id, directory, publish):
        self.lock = threading.RLock()
        self.directory = Path(directory)
        self.publish = publish
        self.tasks = {}
        self.clock = time.monotonic()
        definition = graph_definition(graph)
        self.data = {'schema_version': 1, 'run_id': run_id, 'graph': definition,
                     'status': 'RUNNING', 'started_at': now(), 'ended_at': None,
                     'elapsed_ms': 0, 'version': 0, 'events': [], 'nodes': {
                         n['id']: {'node': n['id'], 'status': 'WAITING', 'attempt': 0,
                                   'input': None, 'output': None, 'state': None, 'metadata': {}, 'error': None}
                         for n in definition['nodes'] if not n['id'].startswith('__')}}
        self.emit('run_start')

    def emit(self, kind, node=None, **fields):
        self.data['version'] += 1
        self.data['elapsed_ms'] = round((time.monotonic() - self.clock) * 1000)
        event = {'type': kind, 'run_id': self.data['run_id'], 'node': node, 'timestamp': now(),
                 'sequence': self.data['version'], **fields}
        self.data['events'].append(event)
        self.data['events'] = self.data['events'][-300:]
        # Monitoring failures never prevent scientific execution.
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            with (self.directory / 'workflow_events.jsonl').open('a') as stream:
                stream.write(json.dumps(event, ensure_ascii=False) + '\n')
            persist(self.data, self.directory)
        except (OSError, ValueError, TypeError) as exc:
            self.data['persistence_error'] = str(exc)[:300]
        self.publish(copy.deepcopy(self.data))

    def on_chain_start(self, serialized, inputs, *, run_id, metadata=None, name=None, **kwargs):
        metadata = metadata or {}
        if name not in self.data['nodes'] or metadata.get('langgraph_node') != name:
            return
        with self.lock:
            node = self.data['nodes'][name]
            self.tasks[str(run_id)] = (name, time.monotonic())
            node.update(status='RUNNING', started_at=now(), ended_at=None, duration_ms=None,
                        attempt=node['attempt'] + 1, input=preview(inputs), state=preview(inputs),
                        output=None, error=None, metadata={**preview(metadata), 'state_scope': 'node_entry'})
            self.emit('node_start', name, metadata={'attempt': node['attempt']})

    def on_chain_end(self, outputs, *, run_id, **kwargs):
        with self.lock:
            task = self.tasks.pop(str(run_id), None)
            if not task:
                return
            name, start = task
            output = outputs if isinstance(outputs, dict) else {}
            execution = output.get('execution_status')
            failed = execution == 'FAILED' or (isinstance(output.get('r_analysis_result'), dict)
                     and output['r_analysis_result'].get('success') is False and execution != 'NOT_RUN')
            status = 'FAILED' if failed else 'SKIPPED' if execution in {'NOT_RUN','SKIPPED','NOT_REQUESTED'} else 'COMPLETED'
            error = str(output.get('errors') or (output.get('r_analysis_result') or {}).get('error') or execution)[:1000] if failed else None
            node = self.data['nodes'][name]
            node.update(status=status, ended_at=now(), duration_ms=round((time.monotonic()-start)*1000),
                        output=preview(outputs), error=error)
            node['metadata']['execution_status'] = execution
            for message in (output.get('messages') or [])[-12:]:
                self.emit('log', name, message=str(message)[:1000])
            self.emit('node_error' if failed else 'node_complete', name,
                      duration_ms=node['duration_ms'], status=status, error=error)

    def on_chain_error(self, error, *, run_id, **kwargs):
        with self.lock:
            task = self.tasks.pop(str(run_id), None)
            if task:
                name, start = task
                node = self.data['nodes'][name]
                node.update(status='FAILED', ended_at=now(), error=str(error)[:1000],
                            duration_ms=round((time.monotonic()-start)*1000))
                self.emit('node_error', name, error=node['error'], duration_ms=node['duration_ms'])

    def log(self, message):
        with self.lock:
            active = [n for n, row in self.data['nodes'].items() if row['status']=='RUNNING']
            self.emit('log', active[-1] if len(active)==1 else None, message=str(message)[:1000])

    def finish(self, terminal, pending=()):
        with self.lock:
            failed = terminal.get('execution_status') == 'FAILED' or any(n['status']=='FAILED' for n in self.data['nodes'].values())
            status = 'FAILED' if failed else 'PAUSED' if pending else 'COMPLETED'
            if terminal.get('execution_status') == 'CANCELLED':
                status = 'CANCELLED'
            for node in self.data['nodes'].values():
                if node['status']=='RUNNING':
                    node.update(status='CANCELLED' if status=='CANCELLED' else 'FAILED', error=terminal.get('message','Run ended before node completion'), ended_at=now())
                elif node['status']=='WAITING' and status!='PAUSED':
                    node['status']='SKIPPED'
            self.data.update(status=status, ended_at=now(), execution_status=terminal.get('execution_status'),
                             pending_nodes=list(pending), error=terminal.get('message'))
            self.emit('run_complete', status=status, error=terminal.get('message'))


def finish_interrupted(server, terminal):
    """Called by supervisor after a killed worker can no longer emit callbacks."""
    with server.run_lock:
        snapshot = copy.deepcopy(server.run_state.get('workflow'))
    if not snapshot or snapshot['status'] != 'RUNNING':
        return
    status = terminal['execution_status']
    snapshot.update(status=status, ended_at=now(), version=snapshot['version']+1,
                    execution_status=status, error=terminal.get('message'))
    for node in snapshot['nodes'].values():
        if node['status']=='RUNNING':
            node.update(status=status, ended_at=now(), error=terminal.get('message'))
        elif node['status']=='WAITING':
            node['status']='SKIPPED'
    event = dict(type='run_complete', run_id=snapshot['run_id'], node=None, timestamp=now(),
                 sequence=snapshot['version'], status=status, error=terminal.get('message'))
    snapshot['events']=(snapshot['events']+[event])[-300:]
    directory = server.REPORT_RUNS_ROOT / snapshot['run_id']
    try:
        persist(snapshot, directory)
        with (directory/'workflow_events.jsonl').open('a') as stream:
            stream.write(json.dumps(event, ensure_ascii=False)+'\n')
    except (OSError, ValueError, TypeError) as exc:
        snapshot['persistence_error'] = str(exc)[:300]
    with server.run_lock:
        server.run_state['workflow']=snapshot
