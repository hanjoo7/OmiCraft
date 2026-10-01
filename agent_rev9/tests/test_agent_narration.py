import copy
import io
import json
import urllib.error

import pytest

from .. import agent_narration as notes, llm_client, orchestrator
from ..configuration import OmiCraftConfig, configuration_scope


def settings(tmp_path, monkeypatch):
    key = tmp_path / 'key.txt'
    key.write_text('private-test-key')
    cfg = {'enabled': True, 'api_type': 'responses', 'auth_header': 'api-key',
           'api_key_file': str(key), 'base_url': 'https://example.invalid/openai/v1',
           'endpoint': '/responses', 'default_model': 'gpt-5.6-sol',
           'available_models': ['gpt-5.6-sol']}
    monkeypatch.setattr(llm_client, 'load_settings', lambda: cfg)
    return cfg


def test_responses_contract_and_secret_redaction(tmp_path, monkeypatch):
    settings(tmp_path, monkeypatch)
    class Opener:
        def open(self, request, timeout):
            body = json.loads(request.data)
            assert request.full_url.endswith('/v1/responses')
            assert request.get_header('Api-key') == 'private-test-key'
            assert 'Authorization' not in request.headers
            assert body['model'] == 'gpt-5.6-sol' and body['store'] is False
            assert 'private-test-key' not in request.data.decode()
            return io.BytesIO(json.dumps({'status': 'completed', 'output': [
                {'type': 'reasoning', 'summary': []},
                {'type': 'message', 'content': [{'type': 'output_text', 'text': 'private-test-key note'}]}]}).encode())
    monkeypatch.setattr(llm_client.urllib.request, 'build_opener', lambda *a: Opener())
    assert llm_client.complete('summary', {'data': 1}).text == '[REDACTED] note'


def test_openai_bearer_auth_and_legacy_key_file(tmp_path, monkeypatch):
    key = tmp_path / 'key.txt'
    key.write_text('api-key: sk-test-secret')
    cfg = {'enabled': True, 'api_type': 'responses',
           'auth_header': 'authorization_bearer', 'api_key_env': 'OPENAI_API_KEY',
           'api_key_file': str(key), 'base_url': 'https://api.openai.com/v1',
           'endpoint': '/responses', 'default_model': 'gpt-4o',
           'available_models': ['gpt-4o']}
    monkeypatch.delenv('OPENAI_API_KEY', raising=False)
    monkeypatch.setattr(llm_client, 'load_settings', lambda: cfg)

    class Opener:
        def open(self, request, timeout):
            assert request.get_header('Authorization') == 'Bearer sk-test-secret'
            assert request.get_header('Api-key') is None
            return io.BytesIO(json.dumps({'status': 'completed', 'output': [
                {'type': 'message', 'content': [
                    {'type': 'output_text', 'text': 'agent response'}]}]}).encode())

    monkeypatch.setattr(llm_client.urllib.request, 'build_opener', lambda *a: Opener())
    assert llm_client.complete('summary', {'data': 1}).text == 'agent response'


def test_environment_key_takes_priority(tmp_path, monkeypatch):
    key = tmp_path / 'key.txt'
    key.write_text('api-key: sk-file-key')
    monkeypatch.setenv('OPENAI_API_KEY', 'sk-environment-key')
    assert llm_client._load_api_key({
        'api_key_env': 'OPENAI_API_KEY', 'api_key_file': str(key)
    }) == 'sk-environment-key'


def test_api_failure_does_not_expose_body(tmp_path, monkeypatch):
    settings(tmp_path, monkeypatch)
    class Opener:
        def open(self, *a, **kw):
            raise urllib.error.HTTPError('https://example.invalid', 401, 'private-test-key', {}, io.BytesIO(b'private-test-key'))
    monkeypatch.setattr(llm_client.urllib.request, 'build_opener', lambda *a: Opener())
    with pytest.raises(llm_client.LLMUnavailable, match='HTTP 401') as exc:
        llm_client.complete('summary', {})
    assert 'private-test-key' not in str(exc.value)


def test_reviews_are_advisory_and_missing_evidence_is_not_pass(tmp_path, monkeypatch):
    cfg = OmiCraftConfig(agent_logging=True)
    cfg.data.base_dir = str(tmp_path)
    state = {'contrasts': [], 'advance_targets': [{'gene_name': 'ERBB2'}], 'critic_verdict': 'HOLD',
             'gate_status': 'PENDING', 'budget': {'max_candidates': 2}}
    before = copy.deepcopy(state)
    monkeypatch.setattr(notes, 'complete', lambda *a, **kw: llm_client.Completion('설명', 'gpt-5.6-sol', 'id', {}, .1))
    with configuration_scope(cfg):
        result = notes.planner_review(state)
    assert state == before
    assert set(result) == {'messages', 'agent_notes'}
    assert any('[Planner:V1]' in message for message in result['messages'])
    assert result['agent_notes']['planner']['checks'][0]['status'] == 'NOT_EVALUATED'
    assert not any('Self-validation 통과' in message for message in result['messages'])


def test_missing_llm_does_not_change_pipeline_status(tmp_path, monkeypatch):
    cfg = OmiCraftConfig(agent_logging=True)
    cfg.data.base_dir = str(tmp_path)
    def unavailable(*a, **kw):
        raise llm_client.LLMUnavailable('LLM disabled')
    monkeypatch.setattr(notes, 'complete', unavailable)
    with configuration_scope(cfg):
        result = notes.discovery_review({'execution_status': 'COMPLETED'})
    assert 'errors' not in result and 'execution_status' not in result and 'critic_verdict' not in result
    assert result['agent_notes']['discovery']['status'] == 'UNAVAILABLE'
    assert result['agent_notes']['discovery']['checks'][3]['status'] == 'NOT_EVALUATED'


@pytest.mark.parametrize('role', ['planner', 'discovery', 'qualification', 'critic', 'dossier'])
def test_annotations_do_not_change_calculation_results(role, monkeypatch):
    state = {'execution_status': 'COMPLETED', 'critic_verdict': 'HOLD',
             'r_analysis_result': {'success': True}, 'qualification_results': {}}
    before = copy.deepcopy(state)
    calls = []
    def narrate(current, selected_role, **kwargs):
        calls.append(selected_role)
        return {'messages': ['[' + selected_role + '] interpretation']}
    monkeypatch.setattr(notes, 'narrate', narrate)
    update = getattr(notes, role + '_review')(state)
    assert calls == [role]
    assert state == before
    assert set(update) == {'messages'}


def test_dry_run_never_calls_llm(monkeypatch):
    from ..graph import _FallbackPipelineGraph
    monkeypatch.setattr(notes, 'complete', lambda *a, **kw: pytest.fail('LLM called in dry run'))
    monkeypatch.setattr(llm_client, 'complete', lambda *a, **kw: pytest.fail('LLM called in dry run'))
    with configuration_scope(OmiCraftConfig(agent_logging=True)):
        assert notes.planner_review({'dry_run': True}) == {}
        events = list(_FallbackPipelineGraph().stream({'dry_run': True}))
    assert events == [{'qualification': {'execution_status': 'NOT_RUN',
        'messages': ['Dry-run: upstream analysis not executed']}}]


def test_gsea_summary_reads_clusterprofiler_columns_without_duplicates(tmp_path):
    directory = tmp_path / '05_gsea'
    directory.mkdir()
    table = 'ID\tDescription\tNES\tp.adjust\nE2F\tE2F Targets\t2.3\t0.001\nOTHER\tOther\t1.1\t0.4\n'
    (directory / 'Hallmark_GSEA_all.tsv').write_text(table)
    (directory / 'Hallmark_GSEA_significant.tsv').write_text(table)
    assert notes._top_pathways({'output_dir': str(tmp_path)}) == [
        {'pathway': 'E2F Targets', 'NES': 2.3, 'padj': 0.001}]
