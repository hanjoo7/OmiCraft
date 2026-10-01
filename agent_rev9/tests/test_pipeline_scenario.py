"""Recorded-data boundaries for the scenario-style pipeline report."""
import json
from pathlib import Path

from ..pipeline_report import assemble_pipeline_report
from ..run_records import RunRecord
from ..server import pipeline_report_document, write_competition_report


def child(root, name, parent, verdict='ADVANCE', reuse=None):
    record = RunRecord(root / name, 'DE_NOVO_BINDER')
    record.finish({'stage': 'complete', 'parent_run_id': parent, 'gene': 'MMP7',
        'running': False, 'execution_status': 'COMPLETED', 'validation_decision': verdict,
        'critic': {'decision': verdict}, 'reuse': reuse,
        'summary': {'candidates': [{'candidate_id': name, 'final_verdict': 'PASS' if verdict == 'ADVANCE' else 'FAIL'}]}})
    return record


def test_recovers_saved_routes_without_counting_followup_pass_as_original(tmp_path):
    original = child(tmp_path, 'design_original', 'upstream_fixture', 'REJECT')
    child(tmp_path, 'design_followup', 'upstream_fixture')
    child(tmp_path, 'design_unrelated', 'upstream_other')
    directory = tmp_path / 'upstream_fixture/agent_results/automatic_design'
    directory.mkdir(parents=True)
    (directory / 'results.json').write_text(json.dumps({'routes': [{'gene': 'MMP7', 'modality': 'DE_NOVO_BINDER', 'run_id': original.directory.name}]}))
    data = {'run_id': 'upstream_fixture', 'summary': {'modalities': []}}
    result = assemble_pipeline_report(data, tmp_path)
    assert data['summary']['modalities'] == []  # historical record not mutated
    assert result['route_source'].endswith('results.json')
    assert result['design_report']['completed_routes'] == 1
    assert result['design_report']['passing_routes'] == 0
    assert [p['run_id'] for p in result['design_report']['followups']] == ['design_followup']
    assert result['design_report']['followups'][0]['validation'] == 'ADVANCE'


def test_actual_evidence_is_deduplicated_and_missing_data_remains_empty(tmp_path):
    directory = tmp_path / 'upstream_fixture/agent_results/qualification'
    directory.mkdir(parents=True)
    row = {'id': 'real-evidence', 'gene': 'MMP7', 'status': 'CONTRADICTORY', 'claim': '<unsafe>'}
    (directory / 'evidence_cards.jsonl').write_text(json.dumps(row)+'\n'+json.dumps(row)+'\n{partial')
    result = assemble_pipeline_report({'run_id': 'upstream_fixture'}, tmp_path)
    assert result['pipeline_details']['evidence'] == [row]
    assert result['pipeline_details']['plan'] == {}
    assert result['pipeline_details']['qualification'] == {}
    missing = assemble_pipeline_report({'run_id': 'upstream_missing'}, tmp_path)
    assert missing['pipeline_details']['evidence'] == []
    assert missing['design_report']['panels'] == []


def test_cannot_read_evidence_outside_runs_root(tmp_path):
    root = tmp_path / 'runs'
    root.mkdir()
    outside = tmp_path / 'external/agent_results/qualification'
    outside.mkdir(parents=True)
    (outside/'evidence_cards.jsonl').write_text('{"id":"private"}')
    result = assemble_pipeline_report({'run_id': '../external'}, root)
    assert result['pipeline_details']['evidence'] == []


def test_standalone_includes_real_data_and_followup_links(tmp_path):
    source = child(tmp_path, 'design_source', 'upstream_fixture')
    artifact = source.directory / 'metric.csv'
    artifact.write_text('metric,value\niptm,0.74\n')
    source.data['artifacts'] = [{'path': str(artifact), 'label': 'Observed metrics'}]
    (source.directory/'run_summary.json').write_text(json.dumps(source.data))
    record = RunRecord(tmp_path/'upstream_fixture', 'UPSTREAM')
    record.finish({'stage':'complete','execution_status':'COMPLETED','execution_profile':'upstream_analysis',
        'summary': {'question':'Real </script><script>bad()</script> question','candidates':[]}})
    path = write_competition_report(record.data, record.directory)
    html = path.read_text()
    assert 'id="pipeline-scenario"' in html
    assert "await fetch('/api/report' + location.search)" not in html
    assert 'Real </script><script>bad()' not in html
    assert '../design_source/metric.csv' in html
    assert '../design_source/report.html' in html
    assert 'mulberry32' not in html
    assert '20,412' not in html


def test_template_has_all_sections_controls_and_no_synthetic_demo():
    html = pipeline_report_document()
    for key in ['question','omics','funnel','qualification','routing','dossiers','design','files']:
        assert 'id="s-'+key+'"' in html
    for removed in ['s-execution', 'pr-execution', 'pr-dag-template', 'pr-agent-']:
        assert removed not in html
    assert 'View 3D' not in html
    assert 'createViewer' not in html
    assert 'id="pr-print"' in html
    for key in ['play','reset','results','speed']:
        assert 'id="pr-'+key+'"' not in html
    assert 'pipelineReportPlayback' not in html
    assert 'id="s-matrix"' not in html
    assert '표적 발굴에서 구조 설계까지' not in html
    assert 'function renderPipelineScenario' in html
    assert 'function renderDesignPanel' in html
    assert 'function mulberry32' not in html
    assert 'TROP2' not in html


def test_critic_rejection_is_not_counted_as_validation_success(tmp_path):
    record = child(tmp_path, 'design_critic_reject', 'upstream_fixture')
    record.data['critic'] = {'decision': 'REJECT'}
    (record.directory/'run_summary.json').write_text(json.dumps(record.data))
    report = assemble_pipeline_report({'run_id': 'upstream_fixture', 'summary': {'modalities': [
        {'modality': 'DE_NOVO_BINDER', 'run_id': record.directory.name, 'gene': 'MMP7'}]}}, tmp_path)
    assert report['design_report']['completed_routes'] == 1
    assert report['design_report']['passing_routes'] == 0


def test_artifact_lookup_does_not_expand_the_whole_pipeline(tmp_path, monkeypatch):
    from .. import run_records
    record = RunRecord(tmp_path/'upstream_fixture', 'UPSTREAM')
    artifact = record.directory/'plot.png'
    artifact.write_bytes(b'fixture')
    record.finish({'stage':'complete','execution_status':'COMPLETED','execution_profile':'upstream_analysis',
        'artifacts':[{'path':str(artifact)}]})
    import agent_rev9.pipeline_report as pipeline
    monkeypatch.setattr(pipeline, 'assemble_pipeline_report', lambda *a: (_ for _ in ()).throw(AssertionError('unnecessary expansion')))
    assert run_records.artifact_path(tmp_path, record.directory.name, None, 0) == artifact


def test_public_default_uses_pipeline_results_and_keeps_saved_fallback(tmp_path, monkeypatch):
    import importlib.util
    saved = tmp_path/'saved'
    saved.mkdir()
    (saved/'run_summary.json').write_text('{"artifacts": []}')
    monkeypatch.setenv('OMICRAFT_PUBLIC_RUN', str(saved))
    path = Path(__file__).resolve().parents[1]/'scripts/serve_review.py'
    spec = importlib.util.spec_from_file_location('scenario_proxy_test', path)
    proxy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(proxy)
    assert proxy.default_pipeline_report() == 'saved'
    for name, candidates in [('upstream_20260101T000000Z_abcdef12', [{'gene':'MMP7'}]),
                             ('upstream_20260102T000000Z_abcdef12', [])]:
        folder = tmp_path/name
        folder.mkdir()
        (folder/'run_summary.json').write_text(json.dumps({'execution_profile':'upstream_analysis',
            'artifacts':[], 'summary':{'candidates':candidates}}))
    assert proxy.default_pipeline_report() == 'upstream_20260101T000000Z_abcdef12'


def test_volcano_preserves_observed_coordinates_and_skips_invalid_rows(tmp_path):
    directory = tmp_path/'upstream_fixture'
    directory.mkdir()
    source = directory/'volcano_data.tsv'
    source.write_text('gene_symbol\tlog2FoldChange\tminus_log10_padj\tdirection\n'
        'OUTLIER\t-28.5\t307.6\tDown_in_TNBC\n'
        'MMP7\t2.1\t12.5\tUp_in_TNBC\n'
        'MISSING\tNA\t5\tNot_DEG\n'
        'INFINITE\t1\tinf\tNot_DEG\n')
    data = {'run_id':'upstream_fixture', 'artifacts':[{'path':str(source)}]}
    result = assemble_pipeline_report(data,tmp_path)['pipeline_details']['volcano']
    assert result['points'] == [['OUTLIER',-28.5,307.6,'Down_in_TNBC'],['MMP7',2.1,12.5,'Up_in_TNBC']]
    assert result['skipped'] == 2
    assert result['source'] == 'upstream_fixture/volcano_data.tsv'


def test_volcano_missing_or_outside_root_is_not_read(tmp_path):
    root = tmp_path/'runs'
    root.mkdir()
    outside = tmp_path/'volcano_data.tsv'
    outside.write_text('gene_symbol\tlog2FoldChange\tminus_log10_padj\nPRIVATE\t1\t3\n')
    data = {'run_id':'upstream_fixture', 'artifacts':[{'path':str(outside)}]}
    assert assemble_pipeline_report(data,root)['pipeline_details']['volcano'] is None
    assert assemble_pipeline_report({'run_id':'missing'},root)['pipeline_details']['volcano'] is None


def test_integrated_design_uses_pipeline_shell_without_discovery_sections():
    html = pipeline_report_document(integrated=True)
    for key in ['question', 'overview', 'comparison', 'design', 'validation', 'files']:
        assert f'id="s-{key}"' in html
    for key in ['omics', 'funnel', 'qualification', 'routing', 'dossiers']:
        assert f'id="s-{key}"' not in html
    assert 'renderDesignScenario(d);' in html
    assert 'id="pr-print"' in html
    assert "getElementById('pr-notice')" not in html


def test_integrated_export_uses_same_design_renderer(tmp_path, monkeypatch):
    from .. import integrated_report
    monkeypatch.setattr(integrated_report, 'assemble_report', lambda data, root: {
        **data, 'integrated_report': {'panels': [], 'passing_candidates': []}})
    data = {'run_id': 'design_fixture', 'gene': 'ERBB2', 'modality': 'ALL',
            'execution_profile': 'therapeutic_design', 'artifacts': []}
    html = write_competition_report(data, tmp_path).read_text()
    assert 'id="s-comparison"' in html
    assert 'renderDesignScenario(d);' in html
    assert "await fetch('/api/report' + location.search)" not in html
    assert 'window.__COMPETITION_RUN__=' in html
