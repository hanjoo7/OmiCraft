import pytest
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('_launch_degrader_staged', root/'launch.py')
launch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launch)
launch.load_package()
from agent_rev9.af3_staged import StagedAF3Runner
from agent_rev9.af3_ternary_tool import AF3TernaryTool
from agent_rev9.af3_msa_cache import TargetMSACache
from agent_rev9.degrader_execution import prepare_reference
from agent_rev9.degrader_status import ternary_execution_status
from agent_rev9.design_execution import execution_budget
from agent_rev9.small_molecule_config import LigandAF3Config
from agent_rev9 import server, web_execution as web, modality_dispatch
from agent_rev9.run_records import read_report


def model(directory):
    script = directory/'af3.py';script.write_text('# Test fixture, not AF3')
    db = directory/'db';db.mkdir();(db/'db.fa').write_text('fixture')
    weights = directory/'weights';weights.mkdir();(weights/'af3.bin').write_text('fixture')
    return LigandAF3Config(python_path=sys.executable, script_path=str(script), model_dir=str(weights),
        database_dir=str(db), target_msa_cache_dir=str(directory/'cache'), seeds=(11,), num_samples=1,
        protein_chain_id='A', ligand_chain_id='C', data_pipeline_timeout_seconds=99,
        inference_timeout_seconds=17, timeout_seconds=1)


CANDIDATE = {'candidate_id':'fixture', 'canonical_smiles':'CCO'}
CONTEXT = {'target_sequence':'AAAA', 'partner_proteins':[
    {'chain':'B','sequence':'CCCC'}, {'chain':'D','sequence':'DDDD'}, {'chain':'E','sequence':'EEEE'}]}


def backend(calls, fail_chain=None, fail_inference=False):
    def execute(runner, command, log, timeout, *, inference):
        calls.append((command, timeout, inference))
        runner.audit['commands'].append(command)
        args = dict(arg[2:].split('=',1) for arg in command if arg.startswith('--'))
        payload = json.loads(Path(args['json_path']).read_text())
        if inference:
            if fail_inference:raise subprocess.TimeoutExpired(command, timeout)
            assert args['run_data_pipeline']=='false' and args['run_inference']=='true'
            assert all('unpairedMsa' in item['protein'] for item in payload['sequences'] if 'protein' in item)
            return
        if payload['name'] == 'msa_' + str(fail_chain).lower():
            raise subprocess.TimeoutExpired(command, timeout)
        assert args['run_data_pipeline']=='true' and args['run_inference']=='false'
        protein=payload['sequences'][0]['protein']
        protein.update(unpairedMsa='>query\n'+protein['sequence']+'\n',pairedMsa='',templates=[])
        folder=Path(args['output_dir']);folder.mkdir(parents=True)
        (folder/'fixture_data.json').write_text(json.dumps(payload))
    return execute


class DegraderStaged(unittest.TestCase):
    def test_each_protein_is_cached_and_second_run_only_infers(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);tool=AF3TernaryTool(model(folder));calls=[]
            with patch.object(StagedAF3Runner,'execute',backend(calls)),patch.object(tool,'parse_output',side_effect=lambda *args: {'status':'success','samples':[]}):
                first=tool.run(candidate=CANDIDATE,context=CONTEXT,output_dir=folder/'first')
                second=tool.run(candidate=CANDIDATE,context=CONTEXT,output_dir=folder/'second')
            self.assertEqual(first['stage_details']['data_pipeline_calls'],4)
            self.assertTrue(all(row['stored'] for row in first['protein_msa_cache']))
            self.assertEqual(second['stage_details']['data_pipeline_calls'],0)
            self.assertEqual([r['cache_status'] for r in second['protein_msa_cache']],['HIT']*4)
            self.assertEqual(len(calls),6)
            self.assertEqual([timeout for _,timeout,infer in calls if infer],[17,17])

    def test_timeout_keeps_completed_chain_caches_and_resumes_remaining(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);tool=AF3TernaryTool(model(folder));calls=[]
            with patch.object(StagedAF3Runner,'execute',backend(calls,fail_chain='D')):
                first=tool.run(candidate=CANDIDATE,context=CONTEXT,output_dir=folder/'first')
            self.assertEqual(first['phase_status'],'TIMED_OUT')
            self.assertFalse(first['inference_started'])
            self.assertEqual(first['protein_msa_cache'][2]['status'],'TIMED_OUT')
            cache=TargetMSACache(tool.config.target_msa_cache_dir,tool.config.script_path,tool.config.database_dir)
            self.assertIsNotNone(cache.load('AAAA'));self.assertIsNotNone(cache.load('CCCC'))
            self.assertIsNone(cache.load('DDDD'))
            with patch.object(StagedAF3Runner,'execute',backend(calls)),patch.object(tool,'parse_output',side_effect=lambda *args: {'status':'success','samples':[]}):
                resumed=tool.run(candidate=CANDIDATE,context=CONTEXT,output_dir=folder/'second')
            self.assertEqual(resumed['stage_details']['data_pipeline_calls'],2)
            self.assertEqual(resumed['status'],'success')

    def test_inference_timeout_is_distinct_and_cache_lock_is_bounded(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);tool=AF3TernaryTool(model(folder))
            with patch.object(StagedAF3Runner,'execute',backend([],fail_inference=True)):
                result=tool.run(candidate=CANDIDATE,context=CONTEXT,output_dir=folder/'first')
            self.assertEqual(result['error_code'],'AF3_INFERENCE_TIMEOUT')
            self.assertTrue(result['inference_started'])
            self.assertEqual(ternary_execution_status({'ternary_prediction':result}),'TIMED_OUT')
            cache=TargetMSACache(tool.config.target_msa_cache_dir,tool.config.script_path,tool.config.database_dir)
            with cache.lock('AAAA'):
                with self.assertRaises(TimeoutError):
                    with cache.lock('AAAA',timeout=.01):pass

    def test_chemical_results_are_published_before_prediction_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            data,config,_=web.design_input({'parent_run_id':'reference_erbb2','gene':'ERBB2','modality':'DEGRADER'})
            snapshots=[];config.update(output_dir=temp,_snapshot_callback=lambda output:snapshots.append(copy.deepcopy(output)))
            def fail(*args,**kwargs):
                self.assertTrue(snapshots)
                self.assertEqual(snapshots[0]['summary']['qc']['conformer_generation_status'],'COMPLETED')
                names={Path(a['path']).name for a in snapshots[0]['artifacts']}
                self.assertTrue({'conformer.sdf','degrader_2d.png','degrader_3d.png','molecular_properties.tsv'} <= names)
                raise RuntimeError('fixture backend failure')
            with patch.object(AF3TernaryTool,'run',side_effect=fail):output=prepare_reference(data,config)
            self.assertEqual(output['summary']['ternary_readiness']['status'],'READY')
            self.assertEqual(output['summary']['ternary_status'],'FAILED')
            self.assertEqual(output['execution_status'],'PARTIAL')
            self.assertGreater(len(output['artifacts']),4)
            self.assertTrue((Path(temp)/'live_result.json').is_file())

    def test_snapshot_is_readable_while_job_runs_and_artifact_urls_stay_stable(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg=web.get_config().model_copy(deep=True);cfg.data.base_dir=temp;cfg.agent_logging=False
            name='design_20260925T010000Z_12345678'
            def execute(data,config):
                folder=Path(config['output_dir']);folder.mkdir(parents=True)
                first=folder/'conformer.sdf';first.write_text('fixture')
                second=folder/'degrader_2d.png';second.write_bytes(b'fixture')
                output={'stage':'reference_evaluation','execution_status':'PARTIAL','validation_decision':'HOLD',
                    'summary':{'candidate_id':'fixture','ternary_prediction':{'status':'running'}},
                    'artifacts':[{'path':str(first),'label':'conformer'}]}
                config['_snapshot_callback'](output)
                live=read_report(temp,run_id=name)
                self.assertTrue(live['running']);self.assertTrue(live['degrader_report']['structures'])
                original=live['artifacts'][0]['href']
                self.assertTrue((Path(temp)/name/'report.html').is_file())
                output['artifacts']=[{'path':str(second),'label':'figure'},*output['artifacts']]
                output['summary']['ternary_prediction']={'status':'backend_failed','error_message':'fixture timeout'}
                return output
            body={'gene':'ERBB2','modality':'DEGRADER','parent_run_id':'reference_erbb2'}
            with patch.object(server,'run_state',{'history':[]}),patch.object(modality_dispatch,'execute',side_effect=execute):
                web.run_design(name,body,{}, {'run_ternary':False},None,cfg,finalize=False)
            final=read_report(temp,run_id=name)
            self.assertFalse(final['running'])
            self.assertTrue(final['artifacts'][0]['path'].endswith('conformer.sdf'))
            self.assertEqual(final['summary']['ternary_status'],'FAILED')

    @pytest.mark.local_io
    def test_timeout_stops_descendant_processes(self):
        import os
        import psutil
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);runner=StagedAF3Runner(AF3TernaryTool(model(folder)),folder)
            child_path=folder/'child.pid'
            script="import subprocess,sys,time; from pathlib import Path; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)']); Path(sys.argv[1]).write_text(str(p.pid)); time.sleep(60)"
            with patch('agent_rev9.af3_staged.gpu_environment',return_value=dict(os.environ)):
                with self.assertRaises(subprocess.TimeoutExpired):
                    runner.execute([sys.executable,'-c',script,str(child_path)],folder/'process.log',.5,inference=False)
            self.assertTrue(child_path.is_file())
            pid=int(child_path.read_text())
            self.assertTrue(not psutil.pid_exists(pid) or psutil.Process(pid).status()==psutil.STATUS_ZOMBIE)

    def test_wrong_sequence_from_preprocessor_is_never_cached(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);tool=AF3TernaryTool(model(folder))
            normal=backend([])
            def corrupt(runner,command,log,timeout,*,inference):
                normal(runner,command,log,timeout,inference=inference)
                for path in (folder/'first').rglob('*_data.json'):
                    data=json.loads(path.read_text());data['sequences'][0]['protein']['sequence']='WWWW';path.write_text(json.dumps(data))
            with patch.object(StagedAF3Runner,'execute',corrupt):
                result=tool.run(candidate=CANDIDATE,context=CONTEXT,output_dir=folder/'first')
            self.assertEqual(result['phase_status'],'FAILED')
            self.assertFalse(result['inference_started'])
            self.assertFalse(list((folder/'cache').glob('*.json')))

    def test_readiness_never_implies_prediction_completed(self):
        for status,expected in [('not_run','NOT_RUN'),('success','COMPLETED'),('backend_failed','FAILED'),('interrupted','INTERRUPTED')]:
            self.assertEqual(ternary_execution_status({'ternary_readiness':{'status':'READY'},'ternary_prediction':{'status':status}}),expected)
        self.assertEqual(ternary_execution_status({'ternary_readiness':{'status':'READY'},'ternary_prediction':{'status':'backend_failed','error_message':'timed out after 1800 seconds'}}),'TIMED_OUT')

    def test_worker_budget_covers_both_stages_and_explicit_limit(self):
        config={'run_ternary':True,'af3':{'data_pipeline_timeout_seconds':5400,'inference_timeout_seconds':1800}}
        self.assertEqual(execution_budget({'modality':'DEGRADER'},config),7800)
        self.assertEqual(execution_budget({'modality':'DEGRADER'},{**config,'execution_timeout_seconds':2}),2)
        with self.assertRaises(ValueError):execution_budget({}, {'execution_timeout_seconds':float('nan')})


if __name__=='__main__':unittest.main()
