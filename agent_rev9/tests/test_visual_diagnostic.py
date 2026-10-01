import gzip
import hashlib
import json

import pytest

from ..run_records import (
    diagnostic_value,
    inspect_diagnostic_artifact,
    select_diagnostic_attempt,
)
from ..server import RUN_REPORT_SCRIPT, report_document


def test_nonfinite_values_are_missing_not_zero():
    result = diagnostic_value({'zero': 0.0, 'nan': float('nan'), 'inf': float('inf'), 'nested': [None, '']})
    assert result == {'zero': 0.0, 'nan': None, 'inf': None, 'nested': [None, '']}
    json.dumps(result, allow_nan=False)


def test_representative_uses_completeness_then_earliest_not_score():
    rows = [
        {'attempt_id': 'a', 'artifact_valid': True, 'required_metric_count': 4, 'started_at': '2026-01-01', 'score': 0.2},
        {'attempt_id': 'b', 'artifact_valid': True, 'required_metric_count': 4, 'started_at': '2026-01-02', 'score': 0.9},
        {'attempt_id': 'c', 'artifact_valid': False, 'required_metric_count': 7, 'started_at': '2025-01-01', 'score': 1.0},
    ]
    assert select_diagnostic_attempt(rows)['attempt_id'] == 'a'
    rows[1]['required_metric_count'] = 5
    assert select_diagnostic_attempt(rows)['attempt_id'] == 'b'
    assert select_diagnostic_attempt([]) is None


def test_source_checksum_and_corruption(tmp_path):
    path = tmp_path / 'source.json'
    path.write_text('{"metric": 0}')
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    result = inspect_diagnostic_artifact(path, digest)
    assert result['parser_success'] and result['checksum_status'] == 'matched_source_manifest'
    path.write_text('{"metric": 9}')
    assert inspect_diagnostic_artifact(path, digest)['error'] == 'source_artifact_checksum_mismatch'
    path.write_text('{broken')
    assert not inspect_diagnostic_artifact(path)['parser_success']


@pytest.mark.parametrize('content', ['', '{}', '[]'])
def test_empty_outputs_are_not_valid(tmp_path, content):
    path = tmp_path / 'empty.json'
    path.write_text(content)
    assert not inspect_diagnostic_artifact(path)['visualizable']


def test_gzip_json_and_table_preserve_zero(tmp_path):
    path = tmp_path / 'confidence.json.gz'
    with gzip.open(path, 'wt') as handle:
        json.dump({'iptm': 0, 'missing': None}, handle)
    assert inspect_diagnostic_artifact(path)['parser_success']
    table = tmp_path / 'metrics.csv'
    table.write_text('metric,value\nclash,0\nmissing,\n')
    assert inspect_diagnostic_artifact(table)['parsed_items'] == 2


def test_invalid_sdf_is_not_visualizable(tmp_path):
    path = tmp_path / 'bad.sdf'
    path.write_text('broken molecule\n$$$$\n')
    assert not inspect_diagnostic_artifact(path)['parser_success']


def test_report_runs_offline_and_retains_diagnostic_sections():
    html = report_document('<head><script src="https://cdn.example/a.js"></script></head>', True)
    assert 'https://cdn.example' not in html
    assert "d.execution_profile==='visual_diagnostic'" in RUN_REPORT_SCRIPT
    for label in ['Omics', 'Small molecule', 'Binder', 'ADC', 'Degrader', 'Issues', 'Provenance']:
        assert label in RUN_REPORT_SCRIPT
    assert 'FRESH ANALYSIS' in RUN_REPORT_SCRIPT
    assert 'Computationally ranked conjugation-site candidates' in RUN_REPORT_SCRIPT
    assert 'new model calls' in RUN_REPORT_SCRIPT
