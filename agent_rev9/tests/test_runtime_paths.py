import json
from pathlib import Path

from .. import runtime_paths


def test_external_config_overrides_local_and_example(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime_paths, 'PACKAGE', tmp_path)
    folder = tmp_path / 'configs'
    folder.mkdir()
    example = folder / 'models.example.json'
    example.write_text('{}')
    assert runtime_paths.config_path('models') == example
    local = folder / 'models.local.json'
    local.write_text('{}')
    assert runtime_paths.config_path('models') == local
    external = tmp_path / 'external.json'
    external.write_text('{}')
    monkeypatch.setenv('OMICRAFT_MODELS_CONFIG', str(external))
    assert runtime_paths.config_path('models', 'OMICRAFT_MODELS_CONFIG') == external


def test_nested_paths_expand_without_corrupting_json(tmp_path, monkeypatch):
    workspace = tmp_path / 'quoted " workspace'
    monkeypatch.setenv('OMICRAFT_WORKSPACE', str(workspace))
    config = tmp_path / 'example.json'
    config.write_text(json.dumps({'path': '${OMICRAFT_WORKSPACE}/data',
        'models': [{'path': '${OMICRAFT_PACKAGE}/configs/binder.json', 'enabled': False}],
        'literal': '${UNSET_OMICRAFT_TEST_VARIABLE}', 'count': 4}))
    value = runtime_paths.load_json(config)
    assert value['path'] == str(workspace / 'data')
    assert Path(value['models'][0]['path']).is_file()
    assert value['models'][0]['enabled'] is False and value['count'] == 4
    assert value['literal'] == '${UNSET_OMICRAFT_TEST_VARIABLE}'


def test_new_workspace_has_empty_history_and_no_stale_runs(tmp_path, monkeypatch):
    import asyncio
    from .. import server
    from ..report_history import report_history
    directory = tmp_path / 'not_created_yet'
    assert report_history(directory, '') == {'runs': []}
    assert not directory.exists()
    monkeypatch.setattr(server, 'run_state', {'run_id': 0, 'running': False,
                                            'result': None, 'history': []})
    events = asyncio.run(server.api_events())
    assert events['erbb2_report_run_id'] is None
    assert events['pipeline_report_run_id'] is None


def test_explicit_saved_report_can_have_undated_name(tmp_path):
    from ..report_history import report_history
    saved = tmp_path / 'reference'
    saved.mkdir()
    (saved / 'run_summary.json').write_text(json.dumps({'execution_status': 'COMPLETED'}))
    history = report_history(tmp_path, 'reference')['runs']
    assert [row['run_id'] for row in history] == ['reference']
