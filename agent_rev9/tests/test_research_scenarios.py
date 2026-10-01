import asyncio
import json
from collections import deque
from pathlib import Path

import pytest
from fastapi import HTTPException

from .. import server, web_execution as web, r_tool, orchestrator, critic_agent, modality_dispatch
from ..configuration import OmiCraftConfig, configuration_scope
from ..research_scenarios import RESEARCH_SCENARIOS


@pytest.mark.parametrize('scenario', RESEARCH_SCENARIOS)
def test_fixed_question_launches_the_correct_real_worker(scenario, monkeypatch):
    jobs = []
    monkeypatch.setattr(server, 'run_state', {'running': False, 'run_id': 0, 'history': []})
    monkeypatch.setattr(server, 'event_queue', deque())
    monkeypatch.setattr(server, 'start_worker', lambda fn, *args, **kw: jobs.append((fn, args, kw)))
    monkeypatch.setattr(server, '_validate_research_question', lambda q: pytest.fail('Fixed question must not require LLM validation'))
    response = asyncio.run(server.api_run({'scenario': scenario, 'question': 'ignored', 'reuse_results': False}))
    assert response.status_code == 200
    assert server.run_state['research_question'] == RESEARCH_SCENARIOS[scenario]
    fn, args, kw = jobs[0]
    if scenario == 'tnbc_discovery':
        assert fn is server.run_pipeline_sync
        assert args[0] == RESEARCH_SCENARIOS[scenario]
        assert args[1]['design_mode'] == 'all_eligible'
        assert args[1]['reuse_results'] is False
        assert not args[1].get('dry_run', False)
    else:
        assert fn is web.run_erbb2_batch
        assert args[1]['modality'] == 'ALL'
        assert args[1]['gene'] == 'ERBB2'
        assert args[1]['reuse_results'] is False
    with pytest.raises(HTTPException):
        asyncio.run(server.api_run({'scenario': 'unknown'}))
    assert len(jobs) == 1


def test_analysis_refreshes_contrasts_and_discovery_counts(tmp_path, monkeypatch):
    from ..discovery_agent import _build_evidence_summary
    out = tmp_path / 'analysis'; (out / '03_de').mkdir(parents=True)
    manifest = {'samples': 1018, 'group_counts': {'TNBC': 157, 'Non_TNBC': 861}}
    (out / 'input_manifest.json').write_text(json.dumps(manifest))
    (out / '03_de/de_summary.json').write_text(json.dumps({'direction_counts': {'Up_in_TNBC': 5288, 'Down_in_TNBC': 2952}}))
    (out / '03_de/DE_all_genes.tsv').write_text('gene_symbol\tlog2FoldChange\tpadj\nERBB2\t2\t0.01\n')
    result = r_tool._parse_r_results(str(out))
    monkeypatch.setattr(r_tool, 'run_r_analysis', lambda **kw: result)
    cfg = OmiCraftConfig(); cfg.data.base_dir = str(tmp_path / 'results')
    with configuration_scope(cfg):
        pending = orchestrator.resolve_contrasts_node({'research_question': RESEARCH_SCENARIOS['tnbc_discovery']})
        assert pending['contrasts'][0]['status'] == 'PARTIAL'
        update = orchestrator.r_analysis_node({'research_question': RESEARCH_SCENARIOS['tnbc_discovery'], **pending})
    contrast = update['contrasts'][0]
    assert (contrast['status'], contrast['actual_n_positive'], contrast['actual_n_reference']) == ('READY', 157, 861)
    assert update['research_plan']['contrasts'] == update['contrasts']
    assert update['sample_manifest'] == manifest
    assert _build_evidence_summary(result, {}, {}, {})['n_deg_significant'] == 8240


def test_provenance_uses_current_metadata_and_missing_checks_are_unknown(tmp_path, monkeypatch):
    from .. import falsification_attacks as attacks
    monkeypatch.setattr(attacks, 'PROCESSED_DIR', '/nonexistent/legacy')
    (tmp_path / 'sample_metadata.tsv').write_text('group\ter_status\tpr_status\ther2_status\nTNBC\tNegative\tNegative\tNegative\nNon_TNBC\tPositive\tNegative\tNegative\n')
    state = {'r_analysis_result': {'input_manifest': {'group_counts': {'TNBC': 1, 'Non_TNBC': 1}}, 'output_dir': str(tmp_path)}}
    result = attacks.run_all_attacks([], state)
    assert result['attacks']['provenance']['status'] == 'PASS'
    assert result['attacks']['null_calibration']['passed'] is None
    assert result['attacks']['confounder']['status'] == 'NOT_EVALUATED'
    assert result['verdict'] == 'HOLD'
    state['r_analysis_result']['input_manifest']['group_counts']['TNBC'] = 2
    result = attacks.run_all_attacks([], state)
    assert result['attacks']['provenance']['status'] == 'FAIL'
    assert result['verdict'] == 'REVISE'


def test_completed_execution_keeps_rejection_reasons_and_terminal_stage(tmp_path, monkeypatch):
    monkeypatch.setattr(modality_dispatch, 'readiness', lambda *a: {'tool_status': 'ready', 'blockers': []})
    monkeypatch.setattr(modality_dispatch, 'execute_binder', lambda *a: {
        'modality': 'DE_NOVO_BINDER', 'execution_status': 'COMPLETED', 'stage': 'binder',
        'is_mock': False, 'blockers': [], 'summary': {'candidates': [{'candidate_id': 'bad', 'is_mock': False,
            'final_verdict': 'FAIL', 'failure_reasons': ['IPTM_BELOW_HARD_LIMIT']}]}})
    result = modality_dispatch.execute({'modality': 'DE_NOVO_BINDER', 'input': {}}, {'output_dir': str(tmp_path), 'dry_run': False})
    execution = [r['execution_status'] for r in result['stages'] if r['stage'] == 'execution']
    assert execution == ['RUNNING', 'COMPLETED']
    assert result['execution_status'] == 'COMPLETED'
    assert result['critic']['decision'] == 'REJECT'
    assert result['critic']['reasons'] == ['bad: IPTM_BELOW_HARD_LIMIT']


def test_validation_warnings_do_not_change_an_advance_decision():
    result = critic_agent.review_modality({'modality': 'DE_NOVO_BINDER', 'execution_status': 'COMPLETED',
        'summary': {}, 'blockers': [], 'is_mock': False, 'validation_decision': 'ADVANCE',
        'validation': {'reasons': ['STRUCTURAL_WARNING']}})
    assert result['decision'] == 'ADVANCE'
    assert result['reasons'] == ['STRUCTURAL_WARNING']


def test_erbb2_batch_running_rows_are_not_labelled_not_run(tmp_path, monkeypatch):
    snapshots = []
    monkeypatch.setattr(server, 'run_state', {'running': True, 'history': []})
    monkeypatch.setattr(server, 'event_queue', deque())
    monkeypatch.setattr(server, 'write_competition_report', lambda data, directory: snapshots.append(json.loads(json.dumps(data))))
    def missing(body):
        raise ValueError('fixture missing input')
    monkeypatch.setattr(web, 'design_input', missing)
    cfg = OmiCraftConfig(); cfg.data.base_dir = str(tmp_path)
    web.run_erbb2_batch('design_20261001T000000Z_12345678', {'gene': 'ERBB2', 'parent_run_id': 'reference_erbb2'}, cfg)
    assert len(snapshots[0]['summary']['modalities']) == 4
    assert all(row['execution_status'] == 'RUNNING' for row in snapshots[0]['summary']['modalities'])
    assert all(row['execution_status'] == 'BLOCKED' for row in snapshots[-1]['summary']['modalities'])


def test_critic_receives_current_gene_context_and_analysis(monkeypatch):
    from types import SimpleNamespace
    from .. import falsification_attacks, llm_client
    seen = {}
    monkeypatch.setattr(falsification_attacks, 'run_all_attacks', lambda *args: {})
    def complete(**kwargs):
        seen.update(kwargs['evidence'])
        return SimpleNamespace(text='{"verdict":"HOLD"}', elapsed=0)
    monkeypatch.setattr(llm_client, 'complete', complete)
    state = {'advance_targets': [{'gene_name': 'MMP7'}],
        'cell_contexts': {'MMP7': {'context_status': 'SINGLE_SUPPORTED'}},
        'depmap_results': {'MMP7': {'dependency_status': 'CONCORDANT'}},
        'r_analysis_result': {'success': True, 'input_manifest': {'group_counts': {'TNBC': 157}},
                              'de_summary': {'n_significant': 8240}}}
    critic_agent._llm_falsification(state)
    assert seen['cell_contexts'] == state['cell_contexts']
    assert seen['depmap_results'] == state['depmap_results']
    assert seen['input_manifest']['group_counts']['TNBC'] == 157
    assert seen['de_summary']['n_significant'] == 8240
