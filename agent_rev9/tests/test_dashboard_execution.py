import threading
import time
from collections import deque

import pytest
from fastapi.testclient import TestClient

from .. import server


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setattr(server, 'run_state', {
        'running': False, 'result': None, 'history': [], 'run_id': 0,
    })
    monkeypatch.setattr(server, 'event_queue', deque())
    # Graph behavior is tested in-process; process supervision has separate tests.
    def start_local(function, *args, **kwargs):
        threading.Thread(target=function, args=args, kwargs=kwargs, daemon=True).start()
    monkeypatch.setattr(server, 'start_worker', start_local)
    with TestClient(server.app) as client:
        yield client


def finished(client):
    deadline = time.monotonic() + 3
    while time.monotonic() < deadline:
        data = client.get('/api/events').json()
        if not data['running']:
            return data
        time.sleep(.01)
    pytest.fail('Worker did not finish')


@pytest.mark.parametrize('body', [{}, {'question': ''}, {'question': '  '},
                                  {'question': 123}, {'question': 'x' * 10001}])
def test_run_rejects_invalid_question(client, body):
    assert client.post('/api/run', json=body).status_code == 422
    assert client.get('/api/status').json()['run_id'] == 0


def test_execution_events_preserve_question_errors_and_cursor(client, monkeypatch):
    inputs = []

    class Graph:
        def stream(self, state, config=None):
            inputs.append(state)
            yield {'qualification_node': {'messages': ['checked'], 'errors': ['missing evidence']}}
            yield {'user_selection_gate_node': {'critic_verdict': 'HOLD'}}

    monkeypatch.setattr(server, 'build_graph', Graph)
    response = client.post('/api/run', json={'question': ' ERBB2 evidence ', 'dry_run': True}).json()
    data = finished(client)
    assert response['status'] == 'started'
    assert inputs[0]['research_question'] == 'ERBB2 evidence'
    assert inputs[0]['dry_run'] is True
    assert [event['type'] for event in data['events']] == ['node', 'node', 'done']
    assert data['result']['execution_status'] == 'COMPLETED_WITH_ERRORS'
    assert data['result']['final_verdict'] == 'HOLD'
    assert data['result']['errors'] == ['missing evidence']
    assert client.get('/api/events', params={'run_id': data['run_id'],
                                            'cursor': data['cursor']}).json()['events'] == []
    assert len(client.get('/api/events?run_id=0&cursor=999').json()['events']) == 3
    assert client.get('/api/events?cursor=-1').status_code == 422


def test_concurrent_start_is_reserved_before_worker_runs(client, monkeypatch):
    release = threading.Event()
    started = threading.Event()
    inputs = []

    class Graph:
        def stream(self, state, config=None):
            inputs.append(state)
            started.set()
            release.wait(3)
            yield {'qualification': {'messages': ['done']}}

    monkeypatch.setattr(server, 'build_graph', Graph)
    try:
        first = client.post('/api/run', json={'question': 'first TNBC study'}).json()
        assert started.wait(2)
        second = client.post('/api/run', json={'question': 'second TNBC study'}).json()
        assert second == {'status': 'already_running', 'run_id': first['run_id']}
        assert len(inputs) == 1
    finally:
        release.set()
        finished(client)


@pytest.mark.parametrize('failure', ['construction', 'stream'])
def test_worker_exception_is_terminal_and_retryable(client, monkeypatch, failure):
    class Graph:
        def __init__(self):
            if failure == 'construction':
                raise RuntimeError('model unavailable')

        def stream(self, state, config=None):
            yield {'planner': {'messages': ['started']}}
            raise RuntimeError('model unavailable')

    monkeypatch.setattr(server, 'build_graph', Graph)
    assert client.post('/api/run', json={'question': 'first TNBC study'}).json()['status'] == 'started'
    data = finished(client)
    assert data['events'][-1]['type'] == 'error'
    assert data['result']['execution_status'] == 'FAILED'
    assert data['result']['message'] == 'model unavailable'
    assert not data['running']
    assert client.post('/api/run', json={'question': 'retry TNBC study'}).json()['status'] == 'started'
    assert finished(client)['run_id'] == 2

@pytest.mark.parametrize('resume_failure', [False, True])
def test_automatic_design_resumes_from_request_mode(client, monkeypatch, resume_failure):
    from types import SimpleNamespace
    calls = []
    class Graph:
        def stream(self, state, config=None):
            calls.append(state)
            if state is not None:
                yield {'critic': {'messages': ['review finished']}}
            elif resume_failure:
                raise RuntimeError('design resume failed')
            else:
                yield {'design': {'execution_status': 'COMPLETED', 'messages': ['design finished']}}
        def get_state(self, config):
            return SimpleNamespace(next=('user_gate',))
    monkeypatch.setattr(server, 'build_graph', Graph)
    response = client.post('/api/run', json={'question': 'TNBC design', 'design_mode': 'all_eligible'})
    assert response.json()['status'] == 'started'
    data = finished(client)
    assert len(calls) == 2 and calls[1] is None
    assert data['result']['execution_status'] == ('FAILED' if resume_failure else 'COMPLETED')


def test_gate_event_preserves_eligibility_status(client, monkeypatch):
    class Graph:
        def stream(self, state, config=None):
            assert config['recursion_limit'] >= 30
            yield {'user_gate': {'gate_status': 'HOLD', 'modality_assessments': [], 'messages': ['No eligible target']}}
    monkeypatch.setattr(server, 'build_graph', Graph)
    response = client.post('/api/run', json={'question': 'TNBC design', 'dry_run': True})
    assert response.json()['status'] == 'started'
    events = finished(client)['events']
    gate = next(e for e in events if e.get('node') == 'user_gate')
    assert gate['gate_status'] == 'HOLD'
