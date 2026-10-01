import pytest
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('_launch_demo',root/'launch.py')
launch=importlib.util.module_from_spec(spec);spec.loader.exec_module(launch);launch.load_package()
from agent_rev9 import server, web_execution as web, design_execution as execution
from agent_rev9.configuration import get_config
from agent_rev9.run_records import RunRecord, read_report
from agent_rev9.integrated_report import assemble_report
from agent_rev9.small_molecule_io import write_json


class DemoCompletion(unittest.TestCase):
    def test_live_report_shows_pass_before_last_route_finishes(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory)/'design_20260924T120000Z_12345678'
            record=RunRecord(folder,'DE_NOVO_BINDER')
            record.data.update(execution_profile='therapeutic_design',running=False)
            record.finish(dict(stage='binder',execution_status='COMPLETED',validation_decision='ADVANCE',
                summary={'candidates':[{'candidate_id':'binder','final_verdict':'PASS','is_mock':False}]}))
            data=dict(running=True,summary={'modalities':[dict(modality='DE_NOVO_BINDER',run_id=folder.name),dict(modality='DEGRADER',running=True)]})
            report=assemble_report(data,directory)['integrated_report']
            self.assertEqual(len(report['panels']),4)
            self.assertEqual(report['workflow_status'],'RUNNING')
            self.assertEqual(report['completed_routes'],1)
            self.assertEqual(len(report['passing_candidates']),1)

    def test_ended_parent_does_not_keep_missing_child_running(self):
        with tempfile.TemporaryDirectory() as directory:
            name='design_20260924T120000Z_52345678'
            data=dict(running=False,summary={'modalities':[dict(modality='ADC',run_id=name,running=True)]})
            report=assemble_report(data,directory)['integrated_report']
            self.assertEqual(report['workflow_status'],'INTERRUPTED')
            self.assertEqual(report['report_status'],'READY')
            self.assertFalse(report['panels'][2]['running'])

    def test_invalid_worker_configuration_releases_run(self):
        cfg=get_config().model_copy(deep=True)
        with tempfile.TemporaryDirectory() as directory:
            cfg.data.base_dir=directory
            name='design_20260924T120000Z_62345678'
            body=dict(parent_run_id='reference_erbb2',gene='ERBB2',modality='ADC')
            with patch.object(server,'run_state',dict(history=[],running=True)):
                result=execution.run_bounded_design(name,body,{}, {'execution_timeout_seconds':'invalid'},None,cfg,finalize=True)
                self.assertFalse(server.run_state['running'])
            self.assertEqual(result['execution_status'],'FAILED')
            self.assertTrue((Path(directory)/name/'report.html').is_file())

    @pytest.mark.local_io
    def test_worker_handles_blocked_input_and_exports_report(self):
        cfg=get_config().model_copy(deep=True);cfg.agent_logging=False
        with tempfile.TemporaryDirectory() as directory:
            cfg.data.base_dir=directory
            name='design_20260924T120000Z_22345678'
            body=dict(parent_run_id='reference_erbb2',gene='ERBB2',modality='ADC')
            with patch.object(server,'run_state',dict(history=[])):
                result=execution.run_bounded_design(name,body,{}, {'execution_timeout_seconds':20},None,cfg)
            self.assertEqual(result['execution_status'],'BLOCKED')
            report=read_report(directory,run_id=name)
            self.assertFalse(report['running'])
            self.assertIn('antibody_complex_missing',report['blockers'])
            self.assertTrue((Path(directory)/name/'report.html').is_file())
            self.assertFalse((Path(directory)/name/'worker_input.json').exists())

    @pytest.mark.local_io
    def test_timeout_ends_worker_and_keeps_report(self):
        cfg=get_config().model_copy(deep=True);cfg.agent_logging=False
        real_popen=subprocess.Popen
        workers=[]
        def sleeper(*args,**kwargs):
            worker=real_popen([sys.executable,'-c','import time; time.sleep(30)'],**kwargs)
            workers.append(worker)
            return worker
        with tempfile.TemporaryDirectory() as directory:
            cfg.data.base_dir=directory
            name='design_20260924T120000Z_32345678'
            body=dict(parent_run_id='reference_erbb2',gene='ERBB2',modality='ADC')
            with patch.object(execution.subprocess,'Popen',side_effect=sleeper),patch.object(server,'run_state',dict(history=[],running=True)):
                result=execution.run_bounded_design(name,body,{}, {'execution_timeout_seconds':1},None,cfg,finalize=True)
                self.assertFalse(server.run_state['running'])
            self.assertIn('EXECUTION_TIMEOUT',result['errors'][0])
            self.assertIsNotNone(workers[0].poll())
            self.assertTrue((Path(directory)/name/'report.html').is_file())
            self.assertFalse(read_report(directory,run_id=name)['running'])

    def test_export_embeds_four_modalities_without_network_fetch(self):
        with tempfile.TemporaryDirectory() as directory:
            record=RunRecord(Path(directory)/'design_20260924T120000Z_42345678','ALL')
            record.data.update(execution_profile='therapeutic_design',running=False,summary={'modalities':[]})
            record.finish(dict(stage='batch',execution_status='PARTIAL'))
            path=server.write_competition_report(record.data,record.directory)
            text=path.read_text()
            self.assertIn('renderIntegratedRun',text)
            self.assertIn('window.__COMPETITION_RUN__=',text)
            self.assertNotIn("await fetch('/api/report' + location.search)",text)
            self.assertIn('DEGRADER',text)
            self.assertTrue((record.directory/'3Dmol-min.js').is_file())

    def test_atomic_json_rejects_nan_without_destroying_previous_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path=write_json(Path(directory)/'result.json',{'status':'COMPLETED'})
            with self.assertRaises(ValueError):write_json(path,{'score':float('nan')})
            self.assertEqual(json.loads(path.read_text()),{'status':'COMPLETED'})


if __name__=='__main__':unittest.main()
