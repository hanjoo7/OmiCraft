"""Real graph control with scientific tools mocked: no API or GPU work."""
from types import SimpleNamespace

import pytest

from .. import graph, planner_agent, critic_agent, contrast_spec, llm_client


def test_replans_stop_after_two_and_preserve_review_count(monkeypatch):
    if not graph._LANGGRAPH_AVAILABLE:
        pytest.skip('LangGraph not installed')
    def unavailable(*args, **kwargs):
        raise llm_client.LLMUnavailable('offline test')
    monkeypatch.setattr(llm_client, 'complete', unavailable)
    monkeypatch.setattr(planner_agent, '_validate_plan', lambda *a: ([], []))
    monkeypatch.setattr(contrast_spec, 'resolve_contrasts', lambda **kw: SimpleNamespace(contrasts=[], to_dict=lambda: {}))
    monkeypatch.setattr(critic_agent, '_apply_deterministic_rules', lambda s: {'verdict':'REPLAN','reason':'test','issues_found':[]})
    monkeypatch.setattr(critic_agent, '_llm_falsification', lambda s: None)
    for name in ['_apear_node','_cell_context_node','_depmap_node','_discovery_llm_node','_qualification_node','_output_tables_node']:
        monkeypatch.setattr(graph, name, lambda s: {})
    monkeypatch.setattr(graph, '_r_analysis_node', lambda s: {'r_analysis_result':{'success':True}})
    monkeypatch.setattr(graph, '_user_gate_node', lambda s: {'gate_status':'HOLD','modality_assessments':[]})
    g=graph.build_graph()
    cfg={'configurable':{'thread_id':'replan-regression'},'recursion_limit':64}
    events=list(g.stream({'research_question':'TNBC','design_mode':'none'},config=cfg))
    planners=[e['planner'] for e in events if 'planner' in e]
    critics=[e['critic'] for e in events if 'critic' in e]
    assert [p['replan_count'] for p in planners]==[0,1,2]
    assert [c['critic_count'] for c in critics]==[1,2,3]
    assert [c['critic_verdict'] for c in critics]==['REPLAN','REPLAN','HOLD']
    assert g.get_state(cfg).next==('user_gate',)
    resumed=list(g.stream(None,config=cfg))
    assert any('user_gate' in e for e in resumed)
    assert not any('design' in e for e in resumed)
    assert not g.get_state(cfg).next


def test_critic_review_issues_are_not_execution_errors(monkeypatch):
    monkeypatch.setattr(critic_agent, '_apply_deterministic_rules', lambda s: {'verdict':'HOLD','reason':'missing evidence','issues_found':['MISSING_EVIDENCE']})
    monkeypatch.setattr(critic_agent, '_llm_falsification', lambda s: None)
    result=critic_agent.critic_node({})
    assert result['critic_verdict']=='HOLD'
    assert result['review_issues']==['MISSING_EVIDENCE']
    assert result['errors']==[]


def test_critic_receives_actual_evidence_and_attack_signature(monkeypatch):
    from .. import falsification_attacks
    seen={}
    def attacks(targets, state):
        assert targets==state['advance_targets']
        return {'attacks':{'test':{'passed':False}}}
    def complete(**kwargs):
        seen.update(kwargs['evidence'])
        return SimpleNamespace(text='{"verdict":"HOLD"}',elapsed=0)
    monkeypatch.setattr(falsification_attacks,'run_all_attacks',attacks)
    monkeypatch.setattr(llm_client,'complete',complete)
    critic_agent._llm_falsification({'advance_targets':[{'gene_name':'MMP7'}], 'contrasts':[{'name':'test'}], 'evidence_cards':[{'id':'E1'}]})
    assert seen['falsification_attacks']['attacks']
    assert seen['contrasts']==[{'name':'test'}]
    assert seen['evidence_cards']==[{'id':'E1'}]
    assert seen['gate_status']=='NOT_YET_EVALUATED'
