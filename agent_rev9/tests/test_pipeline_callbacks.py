"""Progress must arrive during execution without serializing callback functions."""
import pytest
from .. import graph
from ..pipeline_callbacks import pipeline_callbacks


@pytest.mark.parametrize('fallback', [False, True])
def test_progress_during_r_analysis_and_serializable_checkpoint(monkeypatch, fallback):
    from .. import planner_agent, orchestrator
    events = []
    thread = 'progress-test'
    pipeline_callbacks[thread] = {'progress': events.append}
    def planner(state):
        assert events == ['planner']
        return {'messages': ['planned']}
    def r_analysis(state):
        assert events == ['planner', 'r_analysis']
        state['_progress_callback']('R:differential_expression')
        assert events[-1] == 'R:differential_expression'
        # Return a full legacy state to exercise stripping before checkpointing.
        return {**state, 'r_analysis_result': {'success': False}}
    monkeypatch.setattr(planner_agent, 'planner_node', planner)
    monkeypatch.setattr(orchestrator, 'r_analysis_node', r_analysis)
    try:
        runner = graph._FallbackPipelineGraph() if fallback else graph.build_graph()
        cfg = {'configurable': {'thread_id': thread}}
        initial = {'research_question': 'fixture', 'messages': [], 'errors': []}
        results = list(runner.stream(initial, config=cfg))
        assert events == ['planner', 'r_analysis', 'R:differential_expression']
        assert '_progress_callback' not in initial
        assert all('_progress_callback' not in output for row in results for output in row.values())
        if not fallback:
            assert graph._LANGGRAPH_AVAILABLE
            checkpoint = runner.get_state(cfg)
            assert '_progress_callback' not in checkpoint.values
            assert not checkpoint.next
    finally:
        pipeline_callbacks.pop(thread, None)


def test_callbacks_are_scoped_to_execution_and_design_is_forwarded():
    from ..pipeline_callbacks import invoke_node
    events = []
    pipeline_callbacks['active'] = {'progress': events.append, 'design_progress': lambda *args: events.append(args)}
    try:
        def design(state):
            state['_design_progress_callback']('route', 'RUNNING', 'designing', '/report')
            return state
        result = invoke_node('design', design, {}, {'configurable': {'thread_id': 'active'}})
        assert events == ['design', ('route', 'RUNNING', 'designing', '/report')]
        assert result == {}
        assert invoke_node('other', lambda state: state, {}, {'configurable': {'thread_id': 'other'}}) == {}
        assert len(events) == 2
    finally:
        pipeline_callbacks.pop('active', None)
