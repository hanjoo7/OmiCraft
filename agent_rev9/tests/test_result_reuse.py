import copy
import json
from pathlib import Path

from .. import result_reuse as reuse
from ..configuration import get_config, configuration_scope
from ..protein_design_route import AlphaFold3ValidationOutput
from ..rfdiffusion3_stage import RFdiffusion3Output


def historical(tmp_path):
    folder = tmp_path / 'design_20261001T000000Z_12345678'
    folder.mkdir()
    structure = tmp_path / 'target.pdb'
    structure.write_text('original target')
    output = folder / 'contacts.csv'
    output.write_text('residue,contact\n1,2\n')
    body = dict(gene='MMP7', modality='ADC', parent_run_id='upstream_x')
    data = dict(structure_path=str(structure), antigen_chain='A')
    config = dict(seed=11, output_dir=str(folder / 'execution'), gpu_device=4, dry_run=False)
    result = dict(run_id=folder.name, modality='ADC', execution_status='COMPLETED',
                  is_mock=False, running=False, blockers=[], summary={}, stage='adc',
                  artifacts=[{'path': str(output)}], validation_decision='HOLD')
    for name, value in [('request', body), ('input', data), ('config', config), ('run_summary', result)]:
        (folder / (name + '.json')).write_text(json.dumps(value))
    return body, data, config, result


def test_historical_hit_and_invalidated_files(tmp_path):
    body, data, config, result = historical(tmp_path)
    first = reuse.find_design(tmp_path, body, data, config)
    assert first['source_run_id'] == result['run_id'] and first['legacy_import']
    assert reuse.find_design(tmp_path, body, data, {**config, 'gpu_device': 7, 'output_dir': '/new'})
    assert reuse.find_design(tmp_path, {**body, 'reuse_results': False}, data, config) is None
    assert reuse.find_design(tmp_path, {**body, 'gene': 'ERBB2'}, data, config) is None
    assert reuse.find_design(tmp_path, body, data, {**config, 'seed': 12}) is None
    Path(result['artifacts'][0]['path']).unlink()
    assert reuse.find_design(tmp_path, body, data, config) is None


def test_changed_input_and_code_invalidate(tmp_path, monkeypatch):
    body, data, config, result = historical(tmp_path)
    assert reuse.find_design(tmp_path, body, data, config)
    Path(data['structure_path']).write_text('different target')
    assert reuse.find_design(tmp_path, body, data, config) is None
    Path(data['structure_path']).write_text('original target')
    monkeypatch.setattr(reuse, 'code_stamp', lambda names: {'revision': 2})
    assert reuse.find_design(tmp_path, body, data, config) is None


def test_mock_failed_and_running_results_never_imported(tmp_path):
    body, data, config, result = historical(tmp_path)
    path = tmp_path / result['run_id'] / 'run_summary.json'
    for change in [dict(is_mock=True), dict(execution_status='FAILED'), dict(running=True)]:
        path.write_text(json.dumps({**result, **change}))
        assert reuse.find_design(tmp_path, body, data, config) is None


def test_stage_cache_skips_calls_but_rechecks_artifacts_and_force(tmp_path):
    calls = []
    class Tool:
        @reuse.cached_tool('test', [])
        def run(self, *, output_dir, seed, dry_run=False):
            calls.append(seed)
            path = Path(output_dir) / 'model.cif'
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(str(seed))
            return {'status': 'success', 'path': str(path)}
    tool = Tool()
    with reuse.reuse_scope(tmp_path):
        first = tool.run(output_dir=str(tmp_path / 'one'), seed=11)
        second = tool.run(output_dir=str(tmp_path / 'two'), seed=11)
        assert second['provenance']['execution_mode'] == 'REUSED_RESULT'
        assert calls == [11]
        Path(first['path']).unlink()
        tool.run(output_dir=str(tmp_path / 'three'), seed=11)
        tool.run(output_dir=str(tmp_path / 'four'), seed=12)
    with reuse.reuse_scope(tmp_path, False):
        tool.run(output_dir=str(tmp_path / 'five'), seed=12)
    assert calls == [11, 11, 12, 12]


def test_node_cache_reuses_across_run_ids_and_changes_question(tmp_path):
    calls=[]
    def node(state):
        calls.append(state['research_question'])
        return {'advance_targets': [{'gene_name': 'MMP7'}], 'messages': ['finished']}
    cfg=get_config().model_copy(deep=True)
    cfg.data.base_dir=str(tmp_path / 'upstream_one')
    with configuration_scope(cfg):
        base={'research_question':'TNBC', '_run_directory':str(tmp_path / 'upstream_one'), 'run_id':'upstream_one'}
        reuse.cached_node('qualification',node,base)
        second=reuse.cached_node('qualification',node,{**base,'_run_directory':str(tmp_path / 'upstream_two'),'run_id':'upstream_two'})
        assert second['reused_stages'][0]['source_run_id']=='upstream_one'
        reuse.cached_node('qualification',node,{**base,'research_question':'other'})
        reuse.cached_node('qualification',node,{**base,'reuse_results':False})
    assert calls==['TNBC','other','TNBC']


def test_replay_retains_rejection_and_calls_current_validation(tmp_path, monkeypatch):
    from .. import validation_agent, critic_agent
    body,data,config,result=historical(tmp_path)
    entry=reuse.find_design(tmp_path,body,data,config)
    monkeypatch.setattr(validation_agent,'validate_modality',lambda r: {'decision':'REJECT','reasons':['current gate']})
    monkeypatch.setattr(critic_agent,'review_modality',lambda r: {'decision':r['validation_decision']})
    value=reuse.reused_design(entry)
    assert value['validation_decision']=='REJECT'
    assert value['critic']['decision']=='REJECT'
    assert value['reuse']['source_run_id']==result['run_id']


def test_real_graph_reuses_nodes_but_runs_gate_and_critic_again(tmp_path, monkeypatch):
    from .. import graph, planner_agent, orchestrator, discovery_agent, qualification_agent
    counts={}
    def node(name, result):
        def call(state):
            counts[name]=counts.get(name,0)+1
            return copy.deepcopy(result)
        return call
    monkeypatch.setattr(planner_agent,'planner_node',node('planner', {'research_plan':{}}))
    for name in ['apear','cell_context','depmap']:
        monkeypatch.setattr(orchestrator,name+'_node',node(name,{}))
    monkeypatch.setattr(discovery_agent,'discovery_llm_review',node('discovery',{}))
    monkeypatch.setattr(qualification_agent,'qualification_node',node('qualification',{}))
    monkeypatch.setattr(graph,'_r_analysis_node',node('r',{'r_analysis_result':{'success':True}}))
    monkeypatch.setattr(graph,'_output_tables_node',node('output',{}))
    monkeypatch.setattr(graph,'_critic_node',node('critic',{'critic_verdict':'HOLD'}))
    monkeypatch.setattr(graph,'_user_gate_node',node('gate',{'gate_status':'HOLD'}))
    cfg=get_config().model_copy(deep=True)
    for i in range(2):
        folder=tmp_path/('upstream_'+str(i));cfg.data.base_dir=str(folder)
        with configuration_scope(cfg):
            g=graph.build_graph();params={'configurable':{'thread_id':str(i)}}
            events=list(g.stream({'research_question':'TNBC','run_id':folder.name,
                                  '_run_directory':str(folder),'design_mode':'none'},config=params))
            list(g.stream(None,config=params))
            if i: assert any(e.get('planner',{}).get('reused_stages') for e in events)
    assert counts['planner']==counts['qualification']==counts['discovery']==1
    assert counts['critic']==counts['gate']==counts['r']==2


def test_tool_dataclasses_survive_cache_roundtrip(tmp_path):
    calls=[]
    class Tool:
        @reuse.cached_tool('dataclasses', [])
        def run(self, *, output_dir, kind, dry_run=False):
            calls.append(kind)
            path=Path(output_dir)/'model.cif';path.parent.mkdir(parents=True,exist_ok=True);path.write_text('structure')
            if kind=='rfd3':
                return RFdiffusion3Output(status='success',generated_structures=(str(path),),metadata={'is_mock':False})
            return AlphaFold3ValidationOutput(status='success',metrics={'model_path':str(path)})
    tool=Tool()
    with reuse.reuse_scope(tmp_path):
        for kind,cls in [('rfd3',RFdiffusion3Output),('af3',AlphaFold3ValidationOutput)]:
            tool.run(output_dir=str(tmp_path/kind),kind=kind)
            cached=tool.run(output_dir=str(tmp_path/(kind+'_new')),kind=kind)
            assert isinstance(cached,cls)
            metadata=cached.metadata if kind=='rfd3' else cached.runtime_metadata
            assert metadata['execution_mode']=='REUSED_RESULT'
    assert calls==['rfd3','af3']
