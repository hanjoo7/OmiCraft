import json
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from ..af3_msa_cache import TargetMSACache
from ..protein_design_route import AlphaFold3BinderAdapter
from .. import gpu_scheduler as scheduler


def cache_for(tmp_path):
    script = tmp_path / 'af3.py'
    script.write_text('# fixture')
    db = tmp_path / 'db'
    db.mkdir()
    (db / 'uniref90_2022_05.fa').write_text('fixture DB')
    return TargetMSACache(str(tmp_path / 'cache'), str(script), str(db))


def processed(path, sequence='AAAA'):
    path.mkdir(parents=True, exist_ok=True)
    protein = {'id': 'A', 'sequence': sequence, 'unpairedMsa': '>query\n' + sequence + '\n', 'pairedMsa': ''}
    (path / 'test_data.json').write_text(json.dumps({'dialect': 'alphafold3', 'sequences': [{'protein': protein}]}))


def test_cache_rejects_corruption_sequence_and_provenance_changes(tmp_path):
    cache = cache_for(tmp_path)
    processed(tmp_path / 'predictions')
    assert cache.capture(tmp_path / 'predictions', 'AAAA', 'A')
    assert cache.load('AAAA')['unpairedMsa'] == '>query\nAAAA\n'
    assert cache.load('CCCC') is None
    path = cache.path('AAAA')
    path.write_text('{broken')
    assert cache.load('AAAA') is None
    assert cache.capture(tmp_path / 'predictions', 'AAAA', 'A')
    cache.script.write_text('# changed search defaults')
    assert cache.load('AAAA') is None
    assert cache.capture(tmp_path / 'predictions', 'AAAA', 'A')
    (cache.database / 'uniref90_2022_05.fa').write_text('new DB release')
    assert cache.load('AAAA') is None


def test_af3_second_run_reuses_only_target_msa(tmp_path, monkeypatch):
    cache = cache_for(tmp_path)
    model = tmp_path / 'models';model.mkdir()
    adapter = AlphaFold3BinderAdapter(python_path=sys.executable, script_path=str(cache.script),
        model_dir=str(model), database_dir=str(cache.database), target_msa_cache_dir=str(cache.directory),
        binder_msa_mode='single_sequence', use_templates=False)
    inputs = []
    def execute(command, **kwargs):
        args = dict(item[2:].split('=', 1) for item in command if item.startswith('--'))
        payload = json.loads(Path(args['json_path']).read_text())
        inputs.append(payload)
        processed(Path(args['output_dir']))
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(subprocess, 'run', execute)
    first = adapter.run(target_sequence='AAAA', binder_sequence='GGGG', output_dir=str(tmp_path / 'first'))
    second = adapter.run(target_sequence='AAAA', binder_sequence='CCCC', output_dir=str(tmp_path / 'second'))
    assert first.runtime_metadata['target_msa_cache']['stored']
    assert second.runtime_metadata['target_msa_cache']['status'] == 'HIT'
    assert 'unpairedMsa' not in inputs[0]['sequences'][0]['protein']
    assert inputs[1]['sequences'][0]['protein']['unpairedMsa'] == '>query\nAAAA\n'
    assert inputs[1]['sequences'][1]['protein']['sequence'] == 'CCCC'
    assert len(inputs) == 2  # New candidate still requires its own inference.


def test_explicit_msa_and_dry_run_do_not_populate_cache(tmp_path):
    cache = cache_for(tmp_path)
    for explicit, dry_run in [(True, False), (False, True)]:
        directory = tmp_path / str((explicit, dry_run))
        predictions = directory / 'predictions'
        with cache.execution('AAAA', directory, predictions, 'A', explicit=explicit, dry_run=dry_run) as audit:
            processed(predictions)
        assert not audit['stored']
    assert not cache.directory.exists()
    with cache.execution('AAAA', tmp_path / 'failed', tmp_path / 'failed/predictions', 'A') as audit:
        pass
    assert not audit['stored'] and cache.load('AAAA') is None


def test_cold_cache_search_is_shared_between_workers(tmp_path):
    cache = cache_for(tmp_path)
    entered = threading.Event();release = threading.Event();waiting = threading.Event()
    def first():
        with cache.execution('AAAA', tmp_path / 'a', tmp_path / 'a/predictions', 'A'):
            entered.set();assert release.wait(3)
            processed(tmp_path / 'a/predictions')
    def second():
        assert entered.wait(3);waiting.set()
        with cache.execution('AAAA', tmp_path / 'b', tmp_path / 'b/predictions', 'A') as audit:
            return audit['status']
    with ThreadPoolExecutor(max_workers=2) as pool:
        a = pool.submit(first);b = pool.submit(second)
        assert waiting.wait(3);release.set()
        a.result();assert b.result() == 'HIT'


def test_four_workers_receive_distinct_gpus_and_release_on_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(scheduler, 'lock_directory', lambda: tmp_path)
    devices = [scheduler.GPUDevice(i, i, 100000) for i in (4, 5, 6, 7)]
    monkeypatch.setattr(scheduler, 'available_devices', lambda: devices)
    barrier = threading.Barrier(4)
    def job():
        with scheduler.reserve_gpu(timeout=2, poll_seconds=.01) as gpu:
            barrier.wait(timeout=2)
            return gpu.physical
    with ThreadPoolExecutor(max_workers=4) as pool:
        assert sorted(pool.map(lambda _: job(), range(4))) == [4, 5, 6, 7]
    with pytest.raises(RuntimeError):
        with scheduler.reserve_gpu(timeout=.01):
            raise RuntimeError('failed model')
    with scheduler.reserve_gpu(timeout=.01) as gpu:
        assert gpu.physical == 4


def test_gpu_queue_timeout_and_physical_mask(tmp_path, monkeypatch):
    monkeypatch.setattr(scheduler, 'lock_directory', lambda: tmp_path)
    device = scheduler.GPUDevice(7, 0, 100000)
    monkeypatch.setattr(scheduler, 'available_devices', lambda: [device])
    with scheduler.reserve_gpu(timeout=.01):
        with pytest.raises(ValueError, match='대기'):
            with scheduler.reserve_gpu(timeout=.02, poll_seconds=.005):
                pytest.fail('GPU lease was shared')


def test_ligand_consumes_binder_target_cache_without_reusing_binder(tmp_path):
    from ..af3_ligand_tool import AF3LigandTool
    from ..small_molecule_config import LigandAF3Config
    cache = cache_for(tmp_path)
    processed(tmp_path / 'source')
    assert cache.capture(tmp_path / 'source', 'AAAA', 'A')
    config = LigandAF3Config(script_path=str(cache.script), database_dir=str(cache.database),
        target_msa_cache_dir=str(cache.directory), protein_chain_id='C', ligand_chain_id='D')
    tool = AF3LigandTool(config)
    result = tool.run(candidate={'canonical_smiles': 'CCO', 'candidate_id': 'ligand'},
                      context={'target_sequence': 'AAAA'}, output_dir=str(tmp_path / 'ligand'), dry_run=True)
    payload = json.loads(Path(result['input_json']).read_text())
    assert payload['sequences'][0]['protein']['id'] == 'C'
    assert payload['sequences'][0]['protein']['unpairedMsa'] == '>query\nAAAA\n'
    assert 'ligand' in payload['sequences'][1]


def test_design_workers_assign_gpu_in_execute_and_keep_parent_running(tmp_path, monkeypatch):
    from collections import deque
    from .. import web_execution as web, server, modality_dispatch
    from ..configuration import OmiCraftConfig, get_config
    monkeypatch.setattr(scheduler, 'lock_directory', lambda: tmp_path / 'locks')
    monkeypatch.setattr(scheduler, 'available_devices',
                        lambda: [scheduler.GPUDevice(i, i, 100000) for i in (4, 5, 6, 7)])
    monkeypatch.setattr(server, 'run_state', {'running': True, 'history': []})
    monkeypatch.setattr(server, 'event_queue', deque())
    monkeypatch.setattr(server, 'write_competition_report', lambda *a: None)
    monkeypatch.setattr(web, 'register_demo_artifact', lambda *a, **k: None)
    config = OmiCraftConfig();config.data.base_dir = str(tmp_path);config.agent_logging = False
    barrier = threading.Barrier(4)
    allocated = []
    def execute(data, options):
        assert get_config().data.base_dir == str(tmp_path)
        allocated.append(options['screening_options']['gpu_device'])
        barrier.wait(timeout=3)
        return {'stage': 'binder', 'execution_status': 'COMPLETED', 'validation_decision': 'NOT_EVALUATED'}
    monkeypatch.setattr(modality_dispatch, 'execute', execute)
    body = {'gene': 'ERBB2', 'modality': 'DE_NOVO_BINDER', 'parent_run_id': 'reference_erbb2'}
    with ThreadPoolExecutor(max_workers=4) as pool:
        tasks = [pool.submit(web.run_design, 'design_' + str(i), body, {}, {}, None, config, finalize=False)
                 for i in range(4)]
        assert all(f.result()['execution_status'] == 'COMPLETED' for f in tasks)
    assert sorted(allocated) == [4, 5, 6, 7]
    assert server.run_state['running']
    assert not any(e['type'] == 'done' for e in server.event_queue)
    saved = [json.loads(p.read_text())['gpu_device'] for p in tmp_path.glob('design_*/run_summary.json')]
    assert sorted(saved) == [4, 5, 6, 7]
