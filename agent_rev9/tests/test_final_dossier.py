import json
from pathlib import Path

import pytest

from .. import graph
from ..dossier import dossier_node


def test_dossier_preserves_hold_and_counts_only_validated_real_designs(tmp_path):
    good = dict(gene='MMP7', execution_status='COMPLETED', validation_status='ADVANCE', critic_decision='ADVANCE', is_mock=False)
    state = dict(_run_directory=str(tmp_path), run_id='upstream_test', critic_verdict='HOLD',
                 review_issues=['MISSING_EVIDENCE'], evidence_cards=[{'id':'EV1'}], execution_status='PARTIAL',
                 automatic_design_results=[good, {**good,'critic_decision':'REJECT'}, {**good,'is_mock':True}, {**good,'execution_status':'BLOCKED'}])
    update = dossier_node(state)
    data = json.loads(Path(update['dossier_path']).read_text())
    assert data['critic_verdict'] == 'HOLD'
    assert data['validated_designs'] == [good]
    assert data['review_issues'] == ['MISSING_EVIDENCE']
    assert len(data['design_results']) == 4
    assert 'execution_status' not in update  # Assembly must not upgrade PARTIAL.


@pytest.mark.parametrize('allowed', [True, False])
def test_graph_reaches_dossier_after_design_or_blocked_gate(tmp_path, monkeypatch, allowed):
    if not graph._LANGGRAPH_AVAILABLE:
        pytest.skip('LangGraph not installed')
    for name in ['_planner_node','_apear_node','_cell_context_node','_depmap_node','_discovery_llm_node','_qualification_node','_output_tables_node']:
        monkeypatch.setattr(graph, name, lambda s: {})
    monkeypatch.setattr(graph, '_r_analysis_node', lambda s: {'r_analysis_result':{'success':True}})
    monkeypatch.setattr(graph, '_critic_node', lambda s: {'critic_verdict':'HOLD'})
    monkeypatch.setattr(graph, '_user_gate_node', lambda s: {'modality_assessments':[{'gate':{'allowed':allowed}}]})
    monkeypatch.setattr(graph, '_design_node', lambda s: {'execution_status':'PARTIAL','automatic_design_results':[]})
    runner=graph.build_graph();config={'configurable':{'thread_id':'test-dossier'}}
    list(runner.stream({'_run_directory':str(tmp_path),'research_question':'test'},config))
    events=list(runner.stream(None,config))
    assert list(events[-1]) == ['dossier_review']
    assert any('design' in e for e in events) == allowed
    assert runner.get_state(config).next == ()
    assert (tmp_path/'agent_results/final_dossier.json').is_file()
