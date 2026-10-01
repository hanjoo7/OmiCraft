"""Auditable local result reuse for debugging; never turn a failed check into PASS."""
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import asdict, is_dataclass
from functools import wraps
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile

_POLICY = ContextVar('result_reuse_policy', default=None)
# Scheduling and destination changes do not change a scientific calculation.
_EPHEMERAL = {'output_dir', 'base_dir', 'gpu_device', 'device', 'timeout_seconds',
              'execution_timeout_seconds', 'gpu_wait_timeout_seconds', 'cache_dir',
              'target_msa_cache_dir', 'reuse_results', 'dry_run', 'execution_mode',
              'elapsed_seconds', 'source_output_dir', 'source_run_id', 'cache_key', 'reused_stages'}
_PACKAGE = Path(__file__).parent


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def file_stamp(path):
    path = Path(path)
    st = path.stat()
    if path.is_file() and st.st_size <= 8_000_000:
        return {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'size': st.st_size}
    return {'size': st.st_size, 'mtime_ns': st.st_mtime_ns}


def normalize(value):
    if is_dataclass(value):
        value = asdict(value)
    elif hasattr(value, 'model_dump'):
        value = value.model_dump()
    if isinstance(value, dict):
        return {k: normalize(v) for k, v in value.items()
                if k not in _EPHEMERAL and not k.startswith('_') and not callable(v)}
    if isinstance(value, (list, tuple)):
        return [normalize(v) for v in value]
    if isinstance(value, (str, Path)) and str(value).startswith('/'):
        p = Path(value)
        if p.is_file():
            return {'file': file_stamp(p)}
        if p.is_dir():
            return {'directory': str(p.resolve()), 'entries': {
                x.name: file_stamp(x) for x in sorted(p.iterdir())}}
    return value


def code_stamp(names):
    return {name: file_stamp(_PACKAGE / name) for name in names if (_PACKAGE / name).is_file()}


def artifacts(value):
    """Record referenced files, including AF3 per-sample confidence files."""
    found = {}
    def walk(item):
        if isinstance(item, dict):
            for k, v in item.items():
                if k not in {'command', 'environment', 'provenance', 'configuration', 'agent_notes'}:
                    walk(v)
        elif isinstance(item, (tuple, list)):
            for v in item:
                walk(v)
        elif isinstance(item, str) and item.startswith('/') and '\n' not in item:
            p = Path(item)
            if p.is_file():
                found[str(p.resolve())] = file_stamp(p)
    walk(value)
    return found


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=path.parent, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as out:
            json.dump(value, out, ensure_ascii=False, default=str)
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


class ResultCache:
    def __init__(self, root, namespace):
        self.directory = Path(root) / 'cache' / 'results' / namespace

    def read(self, key):
        try:
            entry = json.loads((self.directory / (key + '.json')).read_text())
            if entry['key'] != key or entry['result_hash'] != digest(entry['result']):
                return None
            if any(not Path(p).is_file() or file_stamp(p) != stamp
                   for p, stamp in entry['artifacts'].items()):
                return None
            return entry
        except (OSError, ValueError, TypeError, KeyError, AttributeError):
            return None

    def write(self, key, result, source_run_id='', **extra):
        entry = dict(key=key, result=result, result_hash=digest(result), artifacts=artifacts(result),
                     source_run_id=source_run_id, **extra)
        atomic_json(self.directory / (key + '.json'), entry)
        return entry


@contextmanager
def reuse_scope(root, enabled=True, progress=None):
    token = _POLICY.set((Path(root), enabled, progress))
    try:
        yield
    finally:
        _POLICY.reset(token)


def cached_tool(namespace, dependencies):
    """Cache only successful real tool executions, independently of later gates."""
    def decorate(function):
        @wraps(function)
        def wrapped(self, *args, **kwargs):
            from .resource_usage import measured_tool, record_usage
            measured = measured_tool(function)
            policy = _POLICY.get()
            if not policy or kwargs.get('dry_run'):
                return measured(self, *args, **kwargs)
            cache = ResultCache(policy[0], namespace)
            key = digest({'tool': normalize(vars(self)), 'args': normalize(args),
                          'kwargs': normalize(kwargs), 'code': code_stamp(dependencies)})
            entry = cache.read(key) if policy[1] else None
            if entry:
                record_usage('cache_hit', tool=namespace)
                if policy[2]:
                    policy[2]('result_reuse', 'COMPLETED', namespace + ' · ' + entry['source_run_id'])
                value = copy.deepcopy(entry['result'])
                metadata_key = ('metadata' if entry.get('dataclass') == 'RFdiffusion3Output' else
                                'runtime_metadata' if entry.get('dataclass') else 'provenance')
                value.setdefault(metadata_key, {}).update(
                    execution_mode='REUSED_RESULT', source_run_id=entry['source_run_id'], cache_key=key)
                if entry.get('dataclass'):
                    # The return class lives beside the decorated method.
                    cls = function.__globals__[entry['dataclass']]
                    return cls(**value)
                return value
            output = measured(self, *args, **kwargs)
            value = asdict(output) if is_dataclass(output) else output
            if (isinstance(value, dict) and value.get('status') == 'success'
                    and not value.get('is_mock') and not value.get('metadata', {}).get('is_mock')
                    and not value.get('metrics', {}).get('is_mock')):
                output_dir = kwargs.get('output_dir', '') or (getattr(args[0], 'output_dir', '') if args else '')
                source = next((p.name for p in Path(output_dir).parents if p.name.startswith('design_')), '')
                try:
                    cache.write(key, value, source, dataclass=type(output).__name__ if is_dataclass(output) else None)
                except OSError:
                    pass  # A cache write must not discard completed scientific work.
            return output
        return wrapped
    return decorate


_ROUTE_CODE = ['modality_dispatch.py', 'screening_agent.py', 'screening_tools.py',
    'protein_design_route.py', 'rfdiffusion3_stage.py', 'protein_mpnn_tool.py',
    'validation_agent.py', 'critic_agent.py', 'modality_validation.py',
    'small_molecule_workflow.py', 'af3_ligand_tool.py', 'af3_ternary_tool.py',
    'degrader_execution.py', 'antibody_mapping.py', 'antibody_numbering.py',
    'configuration.py', 'small_molecule_config.py', 'screening_gates.py',
    'binder_analysis.py', 'binder_redesign.py', 'model_runtime.py']


def route_key(body, data, config):
    return digest({'gene': body['gene'], 'modality': body['modality'],
                   'origin': 'reference' if body['parent_run_id'].startswith('reference_') else 'analysis',
                   'input': normalize(data), 'config': normalize(config), 'code': code_stamp(_ROUTE_CODE)})


def completed_real(result):
    return (result.get('execution_status') == 'COMPLETED' and result.get('is_mock') is False
            and not result.get('running') and not result.get('blockers'))


def find_design(root, body, data, config):
    """Import a matching historical result once; subsequent reads check its manifest."""
    if body.get('reuse_results', True) is False:
        return None
    cache = ResultCache(root, 'design')
    key = route_key(body, data, config)
    entry = cache.read(key)
    if entry:
        return entry
    for folder in sorted(Path(root).glob('design_*'), reverse=True):
        try:
            result = json.loads((folder / 'run_summary.json').read_text())
            if not completed_real(result) or result.get('reuse'):
                continue
            old_body = json.loads((folder / 'request.json').read_text())
            old_data = json.loads((folder / 'input.json').read_text())
            old_config = json.loads((folder / 'config.json').read_text())
            if (folder / 'reuse_manifest.json').is_file():
                manifest = json.loads((folder / 'reuse_manifest.json').read_text())
                if manifest.get('key') != key:
                    continue
                if any(not Path(p).is_file() or file_stamp(p) != stamp
                       for p, stamp in manifest.get('artifacts', {}).items()):
                    continue
            else:
                # Historical runs predate fingerprints. Import only exact recorded
                # inputs/settings, mark the provenance limitation, then pin files.
                if route_key(old_body, old_data, old_config) != key:
                    continue
            paths = [a.get('path') for a in result.get('artifacts', []) if isinstance(a, dict)]
            if not paths or any(not p or not Path(p).is_file() for p in paths):
                continue
            entry = cache.write(key, result, folder.name,
                                legacy_import=not (folder / 'reuse_manifest.json').exists())
            atomic_json(folder / 'reuse_manifest.json', {'key': key, 'artifacts': entry['artifacts']})
            return entry
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return None


def remember_design(root, body, data, config, result):
    if not completed_real(result) or result.get('reuse'):
        return
    key = route_key(body, data, config)
    entry = ResultCache(root, 'design').write(key, result, result['run_id'])
    atomic_json(Path(root) / result['run_id'] / 'reuse_manifest.json',
                {'key': key, 'artifacts': entry['artifacts']})


def reused_design(entry, config=None):
    from .validation_agent import validate_modality
    from .critic_agent import review_modality
    output = copy.deepcopy(entry['result'])
    for key in ('run_id', 'parent_run_id', 'selection_event', 'batch_run_id', 'gpu_device'):
        output.pop(key, None)
    output['stages'] = []
    output['reuse'] = {'source_run_id': entry['source_run_id'], 'cache_key': entry['key'],
                       'source_report_url': '/report?run_id=' + entry['source_run_id'],
                       'legacy_import': entry.get('legacy_import', False)}
    output['execution_mode'] = 'REUSED_RESULT'
    if output.get('modality') == 'DE_NOVO_BINDER':
        from .protein_design_route import AlphaFold3BinderAdapter
        from .screening_tools import StructuralScreeningTool, ScreeningOptions
        options = ScreeningOptions.from_state({"screening_options": (config or {}).get("screening_options", {})})
        summary = output['summary']
        contexts = {t['target_id']: t.get('design_conditions', {}) for t in summary.get('targets', [])}
        for candidate in summary.get('candidates', []):
            previous = candidate.get('af3_validation', {})
            if previous.get('status') != 'success':
                continue
            metrics = AlphaFold3BinderAdapter.parse_output(Path(previous['output_dir']))
            af3 = {**previous, 'metrics': metrics}
            checked = StructuralScreeningTool(options.thresholds).run(af3, contexts.get(candidate.get('target_id'), {}),
                                                    is_mock=candidate.get('is_mock', True))
            candidate.update(checked, final_verdict=checked['verdict'], af3_validation=af3)
        for label, verdict in [('pass_count', 'PASS'), ('review_count', 'REVIEW'), ('fail_count', 'FAIL')]:
            summary[label] = sum(c.get('final_verdict') == verdict for c in summary.get('candidates', []))
        summary['accepted_candidate_ids'] = [c['candidate_id'] for c in summary.get('candidates', [])
                                             if c.get('final_verdict') == 'PASS']
        summary['summary'] = {label: summary[label + '_count'] for label in ('pass', 'review', 'fail')}
        summary['final_binder_hits'] = summary['pass_count']
        for branch in summary.get('branch_results', {}).values():
            if branch.get('modality') == 'DE_NOVO_BINDER':
                for key in ('candidates', 'summary', 'final_binder_hits', 'accepted_candidate_ids',
                            'pass_count', 'review_count', 'fail_count'):
                    branch[key] = copy.deepcopy(summary[key])
    summary = output.get('summary', {})
    for view in [summary, *summary.get('branch_results', {}).values()]:
        view['source_tool_calls'] = {k: view[k] for k in ('mpnn_tool_calls', 'af3_tool_calls', 'protenix_tool_calls') if k in view}
        for key in ('mpnn_tool_calls', 'af3_tool_calls', 'protenix_tool_calls', 'mpnn_attempted', 'af3_attempted'):
            if key in view:
                view[key] = 0
        view['execution_mode'] = 'REUSED_RESULT'
    output['validation'] = validate_modality(output)
    output['validation_decision'] = output['validation']['decision']
    output['critic'] = review_modality(output)
    return output


def cached_node(name, function, state):
    # Cache expensive analysis nodes, never the current permission gate or critic.
    directory = state.get('_run_directory')
    if not directory or state.get('dry_run'):
        return function(state)
    from .configuration import get_config
    ignored = {'run_id', 'messages', 'errors', 'agent_notes', 'review_issues', 'reused_stages',
               'user_selection_events', 'design_outcome', 'automatic_design_results'}
    inputs = {k: v for k, v in state.items() if k not in ignored and not k.startswith('_')}
    config = get_config().model_dump()
    deps = [p.name for p in _PACKAGE.glob('*.py') if p.name not in
            {'server.py', 'web_execution.py', 'automatic_design.py', 'result_reuse.py'}]
    # Data settings and model choices are fingerprinted; credentials are not saved.
    environment = {k: v for k, v in os.environ.items() if
                   (k.startswith('HR_') and k != 'HR_OUTPUT_DIR') or
                   k.endswith(('_MODEL', '_MODEL_NAME', '_BASE_URL', '_CONFIG'))}
    key = digest({'state': normalize(inputs), 'config': normalize(config),
                  'environment': normalize(environment), 'code': code_stamp(deps)})
    cache = ResultCache(Path(directory).parent, 'node_' + name)
    entry = cache.read(key) if state.get('reuse_results', True) else None
    if entry:
        from .resource_usage import record_usage
        record_usage('cache_hit', tool=name)
        output = copy.deepcopy(entry['result'])
        output.setdefault('messages', []).append('[Reuse] ' + name + ' · ' + entry['source_run_id'])
        output['reused_stages'] = [{'stage': name, 'source_run_id': entry['source_run_id']}]
        if name == 'r_analysis' and isinstance(output.get('r_analysis_result'), dict):
            output['r_analysis_result']['execution_mode'] = 'REUSED_RESULT'
        return output
    output = function(state)
    if (isinstance(output, dict) and not output.get('errors')
            and output.get('execution_status') not in {'FAILED', 'BLOCKED', 'CANCELLED'}
            and (name != 'r_analysis' or output.get('r_analysis_result', {}).get('success'))):
        try:
            cache.write(key, output, Path(directory).name)
        except OSError:
            pass
    return output


def resume_partial(root, body, data, config, qualification):
    if body.get('reuse_results', True) is False or body['modality'] != 'DE_NOVO_BINDER':
        return None
    from .web_execution import resume_binder_input
    key = route_key(body, data, config)
    for folder in sorted(Path(root).glob('design_*'), reverse=True):
        try:
            result = json.loads((folder / 'run_summary.json').read_text())
            if result.get('running') or result.get('execution_status') not in {'FAILED', 'CANCELLED', 'PARTIAL'}:
                continue
            old_body = json.loads((folder / 'request.json').read_text())
            manifest = folder / 'reuse_request.json'
            if manifest.is_file():
                if json.loads(manifest.read_text()).get('key') != key:
                    continue
            elif route_key(old_body, json.loads((folder / 'input.json').read_text()),
                           json.loads((folder / 'config.json').read_text())) != key:
                continue
            restored, _, _ = resume_binder_input({**body, 'parent_run_id': old_body['parent_run_id']},
                                                  folder.name, qualification)
            return restored, {'source_run_id': folder.name,
                              'source_report_url': '/report?run_id=' + folder.name,
                              'candidate_count': len(restored['protein_design_result']['candidates'])}
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return None
