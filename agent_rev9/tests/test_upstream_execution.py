import json
import subprocess
from pathlib import Path

from .. import orchestrator, r_tool, server
from ..configuration import get_config

RDS_FIELDS = ('counts_rds', 'metadata_rds', 'annotation_rds', 'clinical_rds', 'gene_sets_rds')


def rds_inputs(tmp_path):
    inputs = {}
    for field in RDS_FIELDS:
        path = tmp_path / (field + '.rds')
        path.write_bytes(b'test fixture; never read by R')
        inputs[field] = str(path)
    return inputs


def test_r_node_matches_wrapper_and_isolates_output(tmp_path, monkeypatch):
    inputs = rds_inputs(tmp_path)
    calls = []
    monkeypatch.setattr(r_tool, 'find_rscript', lambda: '/fixture/Rscript')

    def execute(command, **kwargs):
        calls.append(kwargs['env'])
        output = Path(kwargs['env']['HR_OUTPUT_DIR']) / '03_de'
        output.mkdir(parents=True)
        (output / 'DE_all_genes.tsv').write_text('gene_symbol\tlog2FoldChange\nERBB2\t1.2\n')
        (output / 'DEG.tsv').write_text('gene_symbol\tlog2FoldChange\tpadj\nERBB2\t1.2\t0.001\n')
        return subprocess.CompletedProcess(command, 0, stdout='', stderr='')

    monkeypatch.setattr(subprocess, 'run', execute)
    first = orchestrator.r_analysis_node(inputs)
    second = orchestrator.r_analysis_node(inputs)
    assert first['r_analysis_result']['success']
    assert orchestrator._extract_gene_list(first) == ['ERBB2']
    assert calls[0]['HR_COUNTS_RDS'] == inputs['counts_rds']
    assert calls[0]['HR_METADATA_RDS'] == inputs['metadata_rds']
    assert first['r_analysis_result']['output_dir'] != second['r_analysis_result']['output_dir']


def test_missing_r_inputs_stops_graph_and_has_terminal_blockers(tmp_path, monkeypatch):
    monkeypatch.setattr(r_tool, 'find_rscript', lambda: None)
    for field in RDS_FIELDS:
        setattr(get_config().data, field, str(tmp_path / ('missing_' + field + '.rds')))
    from .. import graph, planner_agent
    monkeypatch.setattr(planner_agent, 'planner_node', lambda state: {})
    monkeypatch.setattr(graph, '_LANGGRAPH_AVAILABLE', False)
    state = {'research_question': 'TNBC candidate analysis'}
    events = list(orchestrator.build_graph().stream(state))
    assert [next(iter(event)) for event in events] == ['planner', 'r_analysis']
    result = events[-1]['r_analysis']
    assert result['execution_status'] == 'BLOCKED_INPUT'
    assert 'HR_COUNTS_RDS' in str(result['blockers'])
    assert 'Rscript' in str(result['blockers'])
    assert not (Path(get_config().data.agent_results) / 'output_tables').exists()

    monkeypatch.setattr(server, 'run_state', {'running': False, 'run_id': 0, 'result': None, 'history': []})
    server.run_pipeline_sync(state['research_question'])
    terminal = server.run_state['result']
    assert terminal['execution_status'] == 'BLOCKED_INPUT'
    assert terminal['blockers'] == result['blockers']
    assert not server.run_state['running']


def test_missing_r_runtime_with_existing_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(r_tool, 'find_rscript', lambda: None)
    result = r_tool.run_r_analysis(**rds_inputs(tmp_path), output_dir=str(tmp_path / 'out'))
    assert result['execution_status'] == 'BLOCKED_ENVIRONMENT'
    assert result['missing_inputs'] == []
    assert not (tmp_path / 'out').exists()


def test_r_failure_is_not_reported_as_success(tmp_path, monkeypatch):
    monkeypatch.setattr(r_tool, 'find_rscript', lambda: '/fixture/Rscript')
    monkeypatch.setattr(subprocess, 'run', lambda *a, **kw:
                        subprocess.CompletedProcess(a[0], 1, stdout='', stderr='DESeq2 unavailable'))
    result = orchestrator.r_analysis_node(rds_inputs(tmp_path))
    assert result['execution_status'] == 'FAILED'
    assert 'DESeq2 unavailable' in result['errors'][0]


def test_r_success_requires_gene_table(tmp_path, monkeypatch):
    monkeypatch.setattr(r_tool, 'find_rscript', lambda: '/fixture/Rscript')
    monkeypatch.setattr(subprocess, 'run', lambda *a, **kw:
                        subprocess.CompletedProcess(a[0], 0, stdout='', stderr=''))
    result = orchestrator.r_analysis_node(rds_inputs(tmp_path))
    assert result['execution_status'] == 'FAILED'
    assert result['r_analysis_result']['success'] is False


def test_r_dry_run_does_not_launch_process():
    assert orchestrator.r_analysis_node({'dry_run': True})['execution_status'] == 'NOT_RUN'


def test_apear_missing_gsea_reports_real_reason(tmp_path):
    state = {'r_analysis_result': {'success': True, 'output_dir': str(tmp_path)}}
    result = orchestrator.apear_node(state)
    assert result['execution_status'] == 'PARTIAL'
    assert result['apear_results']['_error'] in result['errors'][0]
    assert "has no attribute 'get'" not in str(result)


def test_apear_uses_r_gsea_subdirectory(tmp_path, monkeypatch):
    from .. import apear_adapter

    seen = []

    def network(gsea_output_dir, output_dir, contrast_ids):
        seen.append(gsea_output_dir)
        return {'test': {'status': 'COMPLETED', 'n_nodes': 4}}

    monkeypatch.setattr(apear_adapter, 'run_apear_for_contrasts', network)
    result = orchestrator.apear_node({'r_analysis_result': {'success': True, 'output_dir': str(tmp_path)}})
    assert Path(seen[0]) == tmp_path / '05_gsea'
    assert result['execution_status'] == 'COMPLETED'
    assert '4개 노드' in result['messages'][0]


def test_sample_manifest_path_is_loaded(tmp_path):
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps({'groups': {'TNBC': 20, 'Luminal_A': 30}}))
    result = orchestrator.resolve_contrasts_node({'research_question': 'TNBC', 'sample_manifest': str(manifest)})
    assert result['contrasts'][0]['status'] == 'READY'


def test_official_depmap_axes_and_gene_ids(tmp_path):
    from ..depmap_module import load_chronos, load_demeter2
    crispr = tmp_path / 'Chronos.csv'
    crispr.write_text(',ERBB2 (2064),ESR1 (2099)\nACH-001,-0.7,0.1\nACH-002,-0.2,nan\n')
    rnai = tmp_path / 'DEMETER2.csv'
    rnai.write_text(',LINE_BREAST\nERBB2 (2064),-0.6\nESR1 (2099),0.05\n')
    assert load_chronos(crispr, ['ERBB2']) == {'ERBB2': {'ACH-001': -0.7, 'ACH-002': -0.2}}
    assert load_demeter2(rnai, ['ERBB2'], {'LINE_BREAST':'ACH-001'}) == {'ERBB2': {'ACH-001': -0.6}}


def test_unknown_pr_is_not_tnbc():
    from ..depmap_module import evaluate_dependency, filter_breast_models
    metadata = {'ACH-1': {'OncotreeLineage':'Breast','ER':'negative','HER2':'negative'}}
    breast = filter_breast_models(metadata)
    assert breast['ACH-1'] != 'TNBC'
    result = evaluate_dependency('EGFR', {'EGFR':{'ACH-1':-0.8}}, breast)
    assert result['crispr_summary']['cohort_scope'] == 'BREAST_LINEAGE'
    assert result['crispr_summary']['n_tnbc_models'] == 0
    assert result['crispr_summary']['n_models'] == 1


def test_cell_context_uses_each_gene_and_needs_donor_evidence(tmp_path):
    from ..cellxgene_adapter import assign_context_status
    table=tmp_path/'scores.tsv'
    table.write_text('gene\tmalignant\tfibroblast\nA\t0.9\t0.1\nB\t0.1\t0.8\n')
    a=assign_context_status(table, {'malignant':5,'fibroblast':5}, gene='A')
    b=assign_context_status(table, {'malignant':5,'fibroblast':5}, gene='B')
    assert (a['top_class'],b['top_class']) == ('malignant','fibroblast')
    assert a['second_score'] == 0.1
    assert a['context_status'] == 'UNRESOLVED'
    assert assign_context_status(table, {}, gene='missing')['context_status'] == 'INSUFFICIENT'


def test_scoped_adc_veto_does_not_reject_all_modalities(monkeypatch):
    from .. import target_qualification as qualification
    monkeypatch.setattr(qualification,'query_uniprot',lambda gene: {
        'status':'FOUND','subcellular_location':['Nucleus'],'is_secreted':False,'has_transmembrane':False,'pdb_ids':[]})
    monkeypatch.setattr(qualification,'query_hpa_normal',lambda *a: {
        'status':'FOUND','normal_tissue':[],'high_expression_tissues':[],'critical_tissue_expression':False})
    monkeypatch.setattr(qualification,'query_chembl',lambda gene: {'status':'FOUND','max_phase':0,'n_compounds':0})
    result=qualification.qualify_gene('GENE',{'dual_contrast':'partial_evidence'}, {}, {})
    assert result['safety_veto']['ADC']['has_veto']
    assert result['bio_confidence'] != 'B0'
    assert result['shortlist']['judgment'] != 'REJECT'
    assert {row['modality']:row['Tier'] for row in result['tier_results']}['ADC'] == 'NOT_APPLICABLE'
    assert result['target_info']['has_structure'] is False


def test_receptor_definition_does_not_use_pam50():
    from ..scripts.prepare_tcga_brca_inputs import receptor_group
    row={'PAM50':'Basal','breast_carcinoma_estrogen_receptor_status':'Negative',
         'breast_carcinoma_progesterone_receptor_status':'Negative',
         'lab_proc_her2_neu_immunohistochemistry_receptor_status':'Negative'}
    assert receptor_group(row)[0]=='TNBC'
    row['breast_carcinoma_progesterone_receptor_status']='Unknown'
    assert receptor_group(row)[0] is None
    row['lab_proc_her2_neu_immunohistochemistry_receptor_status']='Positive'
    assert receptor_group(row)[0]=='Non_TNBC'


def test_r_cache_verifies_input_and_output_content(tmp_path, monkeypatch):
    inputs = rds_inputs(tmp_path)
    get_config().upstream.cache_dir = str(tmp_path/'cache')
    monkeypatch.setattr(r_tool, 'find_rscript', lambda:'/fixture/Rscript')
    calls=[]
    def execute(command, **kwargs):
        calls.append(command)
        out=Path(kwargs['env']['HR_OUTPUT_DIR'])
        (out/'03_de').mkdir()
        (out/'03_de/DE_all_genes.tsv').write_text('gene_symbol\nERBB2\n')
        (out/'report_status.json').write_text('{"status":"SUCCEEDED"}')
        return subprocess.CompletedProcess(command,0,stdout='',stderr='')
    monkeypatch.setattr(subprocess,'run',execute)
    first=r_tool.run_r_analysis(**inputs,output_dir=str(tmp_path/'a'))
    second=r_tool.run_r_analysis(**inputs,output_dir=str(tmp_path/'b'))
    assert len(calls)==1 and second['execution_mode']=='VERIFIED_CACHE_REUSE'
    cached=Path(get_config().upstream.cache_dir)/first['cache_key']/'03_de/DE_all_genes.tsv'
    cached.write_text('modified')
    r_tool.run_r_analysis(**inputs,output_dir=str(tmp_path/'c'))
    assert len(calls)==2
    Path(inputs['metadata_rds']).write_bytes(b'changed input')
    r_tool.run_r_analysis(**inputs,output_dir=str(tmp_path/'d'))
    assert len(calls)==3


def test_server_persists_failed_run_report(tmp_path, monkeypatch):
    monkeypatch.setattr(r_tool,'find_rscript',lambda:None)
    for field in RDS_FIELDS:
        setattr(get_config().data, field, str(tmp_path/(field+'.missing')))
    monkeypatch.setattr(server,'run_state',{'running':False,'run_id':0,'result':None,'history':[]})
    server.run_pipeline_sync('TNBC analysis')
    terminal=server.run_state['result']
    directory=Path(get_config().data.base_dir)/terminal['report_run_id']
    data=json.loads((directory/'run_summary.json').read_text())
    assert data['execution_status']=='BLOCKED'
    assert (directory/'report.html').is_file()
    assert any(event['type']=='progress' for event in server.run_state['history'])


def test_opentargets_current_schema_and_absent_disease(monkeypatch):
    from .. import target_qualification as qualification
    seen=[]
    def response(url, headers=None, payload=None):
        query=json.loads(payload)
        seen.append(query)
        return {'data':{'disease':{'associatedTargets':{'rows':[{'score':0.7,'datatypeScores':[{'id':'somatic_mutation','score':0.4}]}]}}}}
    monkeypatch.setattr(qualification,'_http_get_json',response)
    result=qualification.query_open_targets('ENSG1')
    assert result['status']=='FOUND' and result['somatic_mutation_score']==0.4
    assert seen[0]['variables']['efoId']=='MONDO_0004989'
    monkeypatch.setattr(qualification,'_http_get_json',lambda *a,**kw:{'data':{'disease':None}})
    assert qualification.query_open_targets('ENSG1')['status']=='NOT_FOUND'


def test_chembl_does_not_use_wrong_species_search_hit(monkeypatch):
    from .. import target_qualification as qualification
    hits={'targets':[{'target_chembl_id':'mouse_target','organism':'Mus musculus','target_type':'SINGLE PROTEIN',
                      'target_components':[{'target_component_synonyms':[{'component_synonym':'ERBB2'}]}]}]}
    monkeypatch.setattr(qualification,'_http_get_json',lambda *a,**kw:hits)
    result=qualification.query_chembl('ERBB2')
    assert result['status']=='NOT_FOUND' and result['chembl_target_id'] is None


def test_server_success_persists_candidate_table_and_four_routes(tmp_path, monkeypatch):
    from ..routes import MODALITY_ROUTES
    class FixtureGraph:
        def stream(self, initial, config=None):
            directory=Path(get_config().data.agent_results)
            directory.mkdir(parents=True)
            (directory/'fixture.tsv').write_text('gene\tscore\nG\t1\n')
            yield {'qualification_node':{'execution_status':'COMPLETED','qualification_results':{
                'G':{'bio_confidence':'B3','best_tier':'HOLD','shortlist':{'judgment':'HOLD','reason':'fixture'},
                     'tier_results':[{'modality':m,'Tier':'HOLD'} for m in MODALITY_ROUTES]}}}}
    monkeypatch.setattr(server,'build_graph',lambda:FixtureGraph())
    monkeypatch.setattr(server,'run_state',{'running':False,'run_id':0,'result':None,'history':[]})
    server.run_pipeline_sync('<fixture question>')
    terminal=server.run_state['result']
    assert terminal['execution_status']=='COMPLETED'
    directory=Path(get_config().data.base_dir)/terminal['report_run_id']
    summary=json.loads((directory/'run_summary.json').read_text())
    assert summary['running'] is False
    assert len(summary['summary']['candidates'][0]['routes'])==4
    import re
    text=(directory/'report.html').read_text()
    payload=json.loads(re.search(r'window.__COMPETITION_RUN__=(.*?);</script>',text,re.S).group(1))
    assert payload['question']=='<fixture question>'
    assert '<fixture question>' not in text
    assert len(summary['artifacts'])==2


def test_nonfinite_cell_scores_remain_insufficient(tmp_path):
    from ..cellxgene_adapter import assign_context_status
    path=tmp_path/'scores.tsv'
    path.write_text('gene\tmalignant\tfibroblast\nEMPTY\tnan\t0\n')
    result=assign_context_status(path, {'malignant':5,'fibroblast':5}, gene='EMPTY')
    assert result['context_status']=='INSUFFICIENT'
    json.dumps(result,allow_nan=False)


def test_chembl_phase_comes_from_direct_mechanism(monkeypatch):
    from .. import target_qualification as qualification
    def response(url):
        if 'target/search' in url:
            return {'targets':[{'target_chembl_id':'human_target','organism':'Homo sapiens','target_type':'SINGLE PROTEIN',
                'target_components':[{'target_component_synonyms':[{'component_synonym':'ERBB2'}]}]}]}
        if 'activity.json' in url:
            assert 'pchembl_value__gte=6' in url
            return {'activities':[{'molecule_chembl_id':'M1'}]}
        return {'mechanisms':[{'direct_interaction':1,'max_phase':4,'molecule_chembl_id':'M1'},
                              {'direct_interaction':0,'max_phase':4,'molecule_chembl_id':'M2'}]}
    monkeypatch.setattr(qualification,'_http_get_json',response)
    result=qualification.query_chembl('ERBB2')
    assert result['max_phase']==4 and result['approved_drugs']==['M1']


def test_live_r_logs_preserve_empty_output_contract(tmp_path, monkeypatch):
    inputs=rds_inputs(tmp_path)
    monkeypatch.setattr(r_tool,'find_rscript',lambda:'/fixture/Rscript')
    def execute(command, **kwargs):
        out=Path(kwargs['env']['HR_OUTPUT_DIR'])
        assert list(out.iterdir())==[]
        assert Path(kwargs['stderr'].name).parent!=out
        kwargs['stderr'].write('estimating dispersions\n')
        (out/'03_de').mkdir()
        (out/'03_de/DE_all_genes.tsv').write_text('gene_symbol\nG\n')
        return subprocess.CompletedProcess(command,0,stdout=None,stderr=None)
    monkeypatch.setattr(subprocess,'run',execute)
    result=r_tool.run_r_analysis(**inputs,output_dir=str(tmp_path/'out'),progress_callback=lambda stage:None)
    assert result['success']
    assert (tmp_path/'out/r_stderr.log').read_text()=='estimating dispersions\n'
    assert not list(tmp_path.glob('out-logs-*'))


def test_api_cache_recovers_from_incomplete_json(tmp_path, monkeypatch):
    import hashlib
    import io

    from .. import target_qualification as q
    get_config().upstream.cache_dir = str(tmp_path)
    url = 'https://example.invalid/cache-contract'
    path = tmp_path / 'api' / (hashlib.sha256(url.encode()).hexdigest() + '.json')
    path.parent.mkdir()
    path.write_text('{')
    calls = []
    def response(*args, **kwargs):
        calls.append(args)
        return io.BytesIO(b'{"value": 7}')
    monkeypatch.setattr(q.urllib.request, 'urlopen', response)
    assert q._http_get_json(url) == {'value': 7}
    assert q._http_get_json(url) == {'value': 7}
    assert len(calls) == 1
    assert list(path.parent.iterdir()) == [path]
    assert json.loads(path.read_text())['response'] == {'value': 7}


def test_chembl_unreachable_is_not_absent_target(monkeypatch):
    from .. import target_qualification as q
    monkeypatch.setattr(q, '_http_get_json', lambda *args, **kwargs: None)
    assert q.query_chembl('ESR1')['status'] == 'UNAVAILABLE'
    monkeypatch.setattr(q, '_http_get_json', lambda *args, **kwargs: {'targets': []})
    assert q.query_chembl('NO_MATCH')['status'] == 'NOT_FOUND'
