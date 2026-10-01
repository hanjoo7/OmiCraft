import io
import json
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
import pytest
from ..resource_usage import resource_scope, record_usage, read_usage, report_usage, gpu_allocation


def test_llm_usage_is_durable_and_missing_response_is_not_zero(tmp_path, monkeypatch):
    from .. import llm_client
    from .test_agent_narration import settings
    settings(tmp_path, monkeypatch)
    class Opener:
        def open(self, *args, **kwargs):
            return io.BytesIO(json.dumps({'usage': {'input_tokens': 10, 'output_tokens': 4, 'total_tokens': 14},
                'output': [{'type': 'message', 'content': [{'type': 'output_text', 'text': 'ok'}]}]}).encode())
    monkeypatch.setattr(llm_client.urllib.request, 'build_opener', lambda *args: Opener())
    with resource_scope(tmp_path):
        llm_client.complete('test', {})
        assert json.loads((tmp_path/'resource_usage.json').read_text())['total_tokens'] == 14
        record_usage('llm_request', id='failed')
    usage = read_usage(tmp_path)
    assert usage['total_tokens'] is None
    assert usage['recorded_tokens'] == 14 and usage['llm_usage_missing'] == 1
    assert 'private-test-key' not in (tmp_path/'resource_events.jsonl').read_text()


def test_parallel_tools_and_gpu_lease_on_failure(tmp_path, monkeypatch):
    from .. import resource_usage as module
    clock = iter([100.0, 1900.0])
    monkeypatch.setattr(module.time, 'monotonic', lambda: next(clock))
    with resource_scope(tmp_path):
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(copy_context().run, record_usage, 'tool_call', tool='test') for _ in range(20)]
            for future in futures: future.result()
        with pytest.raises(ValueError):
            with gpu_allocation(4):
                raise ValueError('tool failed')
    assert read_usage(tmp_path)['tool_calls'] == 20
    assert read_usage(tmp_path)['gpu_hours'] == .5


def test_original_children_only_and_no_reused_cost(tmp_path):
    parent = tmp_path/'upstream_fixture'
    with resource_scope(parent):
        record_usage('tool_call', tool='R')
    with resource_scope(tmp_path/'design_original'):
        record_usage('cache_hit', tool='design')
    with resource_scope(tmp_path/'design_followup'):
        record_usage('tool_call', tool='unrelated-to-original')
    usage = report_usage(parent, ['design_original', 'design_original'])
    assert usage['tool_calls'] == 1 and usage['cache_hits'] == 1
    assert usage['total_tokens'] == 0 and usage['gpu_hours'] == 0
    assert read_usage(tmp_path/'old_execution') is None


def test_unfinished_gpu_and_child_usage_remain_unknown(tmp_path):
    parent = tmp_path/'upstream_fixture'
    with resource_scope(parent):
        record_usage('gpu_start', id='interrupted', device=4)
    assert read_usage(parent)['gpu_hours'] is None
    usage = report_usage(parent, ['design_missing'])
    assert usage['total_tokens'] is None and usage['tool_calls'] is None


def test_graph_meter_does_not_enter_checkpoint(tmp_path, monkeypatch):
    from .. import graph, planner_agent, orchestrator
    def planner(state):
        record_usage('llm_request', id='planner')
        record_usage('llm_response', id='planner', usage={'total_tokens': 12})
        return {}
    monkeypatch.setattr(planner_agent, 'planner_node', planner)
    monkeypatch.setattr(orchestrator, 'r_analysis_node', lambda state: {'r_analysis_result': {'success': False}})
    runner = graph.build_graph()
    config = {'configurable': {'thread_id': 'meter'}}
    list(runner.stream({'_run_directory': str(tmp_path), 'reuse_results': False}, config=config))
    assert read_usage(tmp_path)['total_tokens'] == 12
    assert 'resource_usage' not in runner.get_state(config).values


def test_r_reuse_counts_no_new_tool_and_corrupt_cache_recomputes(tmp_path, monkeypatch):
    import subprocess
    from pathlib import Path
    from .. import r_tool
    from ..configuration import get_config
    from .test_upstream_execution import rds_inputs
    inputs = rds_inputs(tmp_path)
    get_config().upstream.cache_dir = str(tmp_path/'cache')
    monkeypatch.setattr(r_tool, 'find_rscript', lambda: '/fixture/Rscript')
    calls = []
    def execute(command, **kwargs):
        calls.append(command)
        out = Path(kwargs['env']['HR_OUTPUT_DIR'])
        (out/'03_de').mkdir()
        (out/'03_de/DE_all_genes.tsv').write_text('gene_symbol\nERBB2\n')
        (out/'report_status.json').write_text('{"status":"SUCCEEDED"}')
        return subprocess.CompletedProcess(command, 0, stdout='', stderr='')
    monkeypatch.setattr(subprocess, 'run', execute)
    with resource_scope(tmp_path/'fresh'):
        first = r_tool.run_r_analysis(**inputs, output_dir=str(tmp_path/'fresh/r'))
    with resource_scope(tmp_path/'reused'):
        second = r_tool.run_r_analysis(**inputs, output_dir=str(tmp_path/'reused/r'))
    assert second['execution_mode'] == 'VERIFIED_CACHE_REUSE' and len(calls) == 1
    assert read_usage(tmp_path/'fresh')['tool_calls'] == 1
    assert read_usage(tmp_path/'reused')['tool_calls'] == 0
    assert read_usage(tmp_path/'reused')['cache_hits'] == 1
    (tmp_path/'cache'/first['cache_key']/'cache_manifest.json').write_text('{partial')
    r_tool.run_r_analysis(**inputs, output_dir=str(tmp_path/'recompute'))
    assert len(calls) == 2
    r_tool.run_r_analysis(**inputs, output_dir=str(tmp_path/'force'), reuse_results=False)
    assert len(calls) == 3


def test_report_uses_metered_original_children(tmp_path):
    from ..pipeline_report import assemble_pipeline_report
    parent = tmp_path/'upstream_fixture'
    with resource_scope(parent):
        record_usage('llm_request', id='a')
        record_usage('llm_response', id='a', usage={'total_tokens': 18})
    child = tmp_path/'design_fixture'
    with resource_scope(child):
        record_usage('tool_call', tool='design')
    report = assemble_pipeline_report({'run_id':parent.name, 'summary':{
        'modalities':[{'gene':'MMP7','modality':'DE_NOVO_BINDER','run_id':child.name}]}}, tmp_path)
    assert report['summary']['resource_usage']['total_tokens'] == 18
    assert report['summary']['resource_usage']['tool_calls'] == 1
