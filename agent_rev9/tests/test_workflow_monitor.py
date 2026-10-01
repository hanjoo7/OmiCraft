import json
from typing import TypedDict

from langgraph.graph import StateGraph, END
from langgraph.checkpoint.memory import MemorySaver
from fastapi.testclient import TestClient

from ..workflow_monitor import WorkflowObserver, preview
from .. import server


class State(TypedDict, total=False):
    value: int


def graph(fail=False):
    def second(s):
        if fail:
            raise RuntimeError('deliberate node failure')
        return {'value': s['value'] + 1}
    builder=StateGraph(State)
    builder.add_node('first', lambda s: {'value': 1})
    builder.add_node('second',second)
    builder.add_edge('__start__','first');builder.add_edge('first','second');builder.add_edge('second',END)
    return builder.compile()


def test_callbacks_are_actual_starts_and_completions(tmp_path):
    g=graph();snapshots=[];o=WorkflowObserver(g,'run_a',tmp_path,snapshots.append)
    list(g.stream({},config={'callbacks':[o]}));o.finish({'execution_status':'COMPLETED'})
    events=o.data['events']
    assert [(e['type'],e['node']) for e in events if e['type'].startswith('node_')]==[
        ('node_start','first'),('node_complete','first'),('node_start','second'),('node_complete','second')]
    assert any(s['nodes']['first']['status']=='RUNNING' and s['nodes']['second']['status']=='WAITING' for s in snapshots)
    assert o.data['nodes']['second']['input']=={'value':1}
    assert o.data['nodes']['second']['output']=={'value':2}
    assert o.data['status']=='COMPLETED'
    assert json.loads((tmp_path/'workflow.json').read_text())['run_id']=='run_a'


def test_failure_marks_node_and_run_and_does_not_mix_runs(tmp_path):
    g=graph(True);first=WorkflowObserver(g,'run_fail',tmp_path/'fail',lambda s:None)
    try:list(g.stream({},config={'callbacks':[first]}))
    except RuntimeError:pass
    first.finish({'execution_status':'FAILED','message':'deliberate node failure'})
    assert first.data['nodes']['second']['status']=='FAILED'
    assert first.data['nodes']['second']['error']=='deliberate node failure'
    assert first.data['status']=='FAILED'
    other=WorkflowObserver(graph(),'run_next',tmp_path/'next',lambda s:None)
    assert all(e['run_id']=='run_next' for e in other.data['events'])
    assert all(n['status']=='WAITING' for n in other.data['nodes'].values())


def test_preview_is_bounded_and_redacted():
    data=preview({'api_key':'sensitive','matrix':[[1]*10000]*10000,'large':'x'*100000,'nan':float('nan')})
    assert data['api_key']=='[redacted]'
    assert len(json.dumps(data))<15000
    json.dumps(data,allow_nan=False)


def test_interrupt_is_paused_and_does_not_skip_future_nodes(tmp_path):
    b=StateGraph(State);b.add_node('first',lambda s:{'value':1});b.add_node('gate',lambda s:s)
    b.add_edge('__start__','first');b.add_edge('first','gate');b.add_edge('gate',END)
    g=b.compile(checkpointer=MemorySaver(),interrupt_before=['gate']);o=WorkflowObserver(g,'pause',tmp_path,lambda s:None)
    cfg={'callbacks':[o],'configurable':{'thread_id':'test'}}
    list(g.stream({},cfg));o.finish({'execution_status':'COMPLETED'},g.get_state(cfg).next)
    assert o.data['status']=='PAUSED' and o.data['nodes']['gate']['status']=='WAITING'


def test_api_persistent_runs_and_unknown_ids(tmp_path,monkeypatch):
    monkeypatch.setattr(server,'REPORT_RUNS_ROOT',tmp_path)
    monkeypatch.setattr(server,'run_state',{'workflow':None})
    o=WorkflowObserver(graph(),'workflow_test',tmp_path/'workflow_test',lambda s:None)
    o.finish({'execution_status':'COMPLETED'})
    with TestClient(server.app) as client:
        assert client.get('/api/workflow/state?run_id=workflow_test').json()['status']=='COMPLETED'
        assert client.get('/api/workflow/runs').json()[0]['run_id']=='workflow_test'
        assert client.get('/api/workflow/state?run_id=missing').status_code==404
        assert client.get('/api/workflow/state?run_id=../secret').status_code==422
        assert client.get('/workflow').status_code==200


def test_supervisor_cancellation_finishes_monitor(tmp_path,monkeypatch):
    from ..workflow_monitor import finish_interrupted
    monkeypatch.setattr(server,'REPORT_RUNS_ROOT',tmp_path)
    g=graph();o=WorkflowObserver(g,'workflow_cancel',tmp_path/'workflow_cancel',lambda s:None)
    o.on_chain_start(None,{'value':0},run_id='task',metadata={'langgraph_node':'first'},name='first')
    monkeypatch.setattr(server,'run_state',{'workflow':o.data})
    finish_interrupted(server,{'execution_status':'CANCELLED','message':'User stopped run'})
    assert server.run_state['workflow']['status']=='CANCELLED'
    assert server.run_state['workflow']['nodes']['first']['status']=='CANCELLED'
    assert server.run_state['workflow']['nodes']['second']['status']=='SKIPPED'
    assert json.loads((tmp_path/'workflow_cancel/workflow.json').read_text())['status']=='CANCELLED'


def test_observer_does_not_change_graph_result(tmp_path):
    g=graph();expected=list(g.stream({}))
    o=WorkflowObserver(g,'unchanged',tmp_path,lambda s:None)
    assert list(g.stream({},config={'callbacks':[o]}))==expected
