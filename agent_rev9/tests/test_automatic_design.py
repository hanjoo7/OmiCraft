import copy
import json
from pathlib import Path

from .. import automatic_design as auto, orchestrator, web_execution as web, modality_dispatch, server, design_execution
from ..configuration import OmiCraftConfig, configuration_scope


def candidate(gene):
    return {'gene_name': gene, 'tier': '2B', 'modality': 'DE_NOVO_BINDER',
            'tier_results': [{'modality': 'DE_NOVO_BINDER', 'Tier': '2B'},
                             {'modality': 'ADC', 'Tier': '2B'},
                             {'modality': 'SMALL_MOLECULE', 'Tier': '3'},
                             {'modality': 'DEGRADER', 'Tier': 'NOT_APPLICABLE'}],
            'safety_veto': {'ADC': {'has_veto': True}}}


def test_bulk_selection_preserves_qualification_and_veto():
    state = {'advance_targets': [candidate('A')]}
    before = copy.deepcopy(state)
    events = auto.select_routes(state)
    assert state == before
    assert [(e['route_id'], e['selected_by'], e['selection_scope']) for e in events] == [
        ('A:DE_NOVO_BINDER', 'user', 'all_eligible')]
    assert auto.select_routes({**state, 'user_selection_events': events}) == events


def test_batch_continues_after_failure_and_records_blocked_routes(tmp_path, monkeypatch):
    parent = tmp_path/'upstream_20260923T000000Z_12345678'
    cfg = OmiCraftConfig();cfg.data.base_dir = str(parent)
    targets = [candidate('A'), candidate('B')]
    state = {'design_mode': 'all_eligible', 'advance_targets': targets,
             'qualification_results': {t['gene_name']:t for t in targets}, '_run_directory': str(parent)}
    calls = []
    import threading
    barrier = threading.Barrier(2)
    monkeypatch.setattr(auto, 'route_input', lambda *a: {'asset_id': 'fixture', 'target_chain': 'A'})
    monkeypatch.setattr(web, 'design_input', lambda body: ({}, {}, {}))
    monkeypatch.setattr(modality_dispatch, 'readiness', lambda *a: {'blockers': []})
    monkeypatch.setattr(server, 'write_competition_report', lambda data,directory: (Path(directory)/'report.html').write_text('report'))
    def run(identifier, body, data, model_cfg, q, app_cfg, finalize):
        calls.append(body['gene'])
        barrier.wait(timeout=3)
        child = Path(app_cfg.data.base_dir)/identifier;child.mkdir()
        (child/'run_summary.json').write_text(json.dumps({'critic': {'decision': 'HOLD'}}))
        return {'execution_status': 'FAILED' if body['gene']=='A' else 'COMPLETED',
                'scientific_validation_status': 'HOLD', 'errors': []}
    monkeypatch.setattr(design_execution, 'run_bounded_design', run)
    with configuration_scope(cfg):
        result = auto.automatic_design_node(state)
    assert sorted(calls) == ['A','B']
    assert len(result['automatic_design_results']) == 8
    assert result['execution_status'] == 'PARTIAL'
    assert result['design_outcome']['validated_successes']==0
    assert result['design_outcome']['completed_routes']==1
    assert result['design_outcome']['status']=='NO_VALIDATED_CANDIDATE'
    blocked = [r for r in result['automatic_design_results'] if r['execution_status']=='BLOCKED']
    assert len(blocked)==6
    assert all((tmp_path/r['run_id']/'report.html').is_file() for r in blocked)
    assert all(r['validation_status']=='NOT_EVALUATED' for r in blocked)
    assert (parent/'agent_results/automatic_design/results.json').is_file()


def test_graph_automatic_design_runs_once_and_dry_run_never_executes(monkeypatch):
    from .. import graph
    calls=[]
    for name in ['_planner_node','_apear_node','_cell_context_node','_depmap_node','_discovery_llm_node','_qualification_node','_output_tables_node']:
        monkeypatch.setattr(graph,name,lambda state: {})
    monkeypatch.setattr(graph,'_r_analysis_node',lambda state: {'r_analysis_result':{'success':True}})
    monkeypatch.setattr(graph,'_critic_node',lambda state: {'critic_verdict':'APPROVE'})
    monkeypatch.setattr(graph,'_user_gate_node',lambda state: {'gate_status':'ALLOWED'})
    original=auto.automatic_design_node
    def automatic(state):
        if state.get('dry_run'):return original(state)
        calls.append('automatic')
        return {'automatic_design_results': [], 'execution_status':'PARTIAL'}
    monkeypatch.setattr(auto,'automatic_design_node',automatic)
    for dry in [False,True]:
        g=graph.build_graph();cfg={'configurable':{'thread_id':str(dry)}}
        list(g.stream({'design_mode':'all_eligible','dry_run':dry},config=cfg))
        list(g.stream(None,config=cfg))
    assert calls==['automatic']


def test_no_ligand_or_antibody_is_invented(monkeypatch):
    import pytest
    monkeypatch.setattr(web,'assets',lambda: [])
    for modality in ('SMALL_MOLECULE','ADC','DEGRADER'):
        with pytest.raises(ValueError):auto.route_input('FOXA1',modality)


def test_erbb2_batch_attempts_all_four_even_when_inputs_are_missing(tmp_path, monkeypatch):
    from collections import deque
    calls=[]
    def missing(body):
        calls.append(body['modality'])
        raise ValueError('fixture input missing')
    monkeypatch.setattr(web,'design_input',missing)
    monkeypatch.setattr(server,'write_competition_report',lambda *a: None)
    monkeypatch.setattr(server,'run_state',{'running':True,'history':[]})
    monkeypatch.setattr(server,'event_queue',deque())
    cfg=OmiCraftConfig();cfg.data.base_dir=str(tmp_path)
    web.run_erbb2_batch('design_20260923T000000Z_87654321',
        {'parent_run_id':'reference_erbb2','gene':'ERBB2','modality':'ALL'},cfg)
    assert sorted(calls)==sorted(web.ERBB2_MODALITIES)
    assert not server.run_state['running']
    assert server.run_state['result']['execution_status']=='PARTIAL'


def test_run_api_forwards_explicit_design_mode(monkeypatch):
    import asyncio
    from collections import deque
    captured={}
    def worker(function,*args,**kwargs):captured['args']=args
    monkeypatch.setattr(server,'start_worker',worker)
    monkeypatch.setattr(server,'run_state',{'running':False,'history':[],'run_id':0})
    monkeypatch.setattr(server,'event_queue',deque())
    response=asyncio.run(server.api_run({'question':'TNBC','design_mode':'all_eligible'}))
    assert response.status_code==200
    assert captured['args'][1]['design_mode']=='all_eligible'
