import importlib.util
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('_launch_pipeline_report',root/'launch.py')
launch=importlib.util.module_from_spec(spec);spec.loader.exec_module(launch);launch.load_package()
from agent_rev9.configuration import OmiCraftConfig, configuration_scope
from agent_rev9.orchestrator import _extract_gene_list, qualification_node
from agent_rev9.pipeline_report import assemble_pipeline_report
from agent_rev9.run_records import RunRecord, read_report
from agent_rev9.upstream_reporting import finish_report
from agent_rev9 import target_qualification, upstream_reporting


class PipelineDesignReport(unittest.TestCase):
    def test_upregulated_only_input_respects_70_percent_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'deg.tsv'
            path.write_text('gene_symbol\tlog2FoldChange\tpadj\n'+''.join(f'G{i}\t2\t{i/100000:.6f}\n' for i in range(1,152))+'NOT_SIGNIFICANT\t2\t0.9\n')
            state={'r_analysis_result':{'output_files':{'deg_table':str(path)}}}
            with configuration_scope(OmiCraftConfig()):
                self.assertEqual(_extract_gene_list(state),[f'G{i}' for i in range(1,71)])
                self.assertEqual(_extract_gene_list({**state,'gene_list':['EXPLICIT','EXPLICIT']}),['EXPLICIT'])

    def test_qualification_forwards_all_100_and_records_actual_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg=OmiCraftConfig();cfg.data.base_dir=directory
            genes=[f'G{i}' for i in range(100)]
            result={'results':{g:{'best_tier':'HOLD','tier_results':[]} for g in genes},'advance':[],'hold':genes,'reject':[]}
            with configuration_scope(cfg),patch.object(target_qualification,'run_qualification_pipeline',return_value=result) as qualify:
                output=qualification_node({'gene_list':genes})
            self.assertEqual(qualify.call_args.kwargs['genes'],genes)
            self.assertEqual(output['candidate_search']['evaluated_targets'],100)
            self.assertEqual(len(output['hold_targets']),100)

    def test_previous_pass_is_separate_from_current_reject(self):
        with tempfile.TemporaryDirectory() as directory:
            def child(name, verdict):
                record=RunRecord(Path(directory)/name,'DE_NOVO_BINDER')
                record.finish({'stage':'binder','execution_status':'COMPLETED','validation_decision':verdict,'running':False,
                    'summary':{'candidates':[{'candidate_id':name,'final_verdict':'PASS' if verdict=='ADVANCE' else 'FAIL'}]}})
                return record
            current=child('design_current','REJECT');passed=child('design_previous','ADVANCE')
            reference=RunRecord(Path(directory)/'reference_batch','ALL')
            reference.finish({'stage':'batch','execution_status':'PARTIAL','execution_profile':'therapeutic_design','gene':'ERBB2','running':False,
                'summary':{'modalities':[{'gene':'ERBB2','modality':'DE_NOVO_BINDER','run_id':passed.directory.name}]}})
            data={'run_id':'upstream_current','summary':{'candidates':[{'gene':'FOXA1'}],
                'modalities':[{'gene':'FOXA1','modality':'DE_NOVO_BINDER','run_id':current.directory.name}],
                'design_reference_run':reference.directory.name}}
            report=assemble_pipeline_report(data,directory)['design_report']
            self.assertEqual(report['passing_routes'],0)
            self.assertEqual(report['completed_routes'],1)
            self.assertEqual(report['candidate_count'],1)
            self.assertEqual(report['panels'][0]['validation'],'REJECT')
            self.assertEqual(report['reference']['source_kind'],'PREVIOUSLY_RECORDED')
            self.assertEqual(report['reference']['panels'][0]['validation'],'ADVANCE')

    def test_final_export_survives_plot_failure_and_includes_design(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg=OmiCraftConfig();cfg.data.base_dir=directory
            child=RunRecord(Path(directory)/'design_fixture','DE_NOVO_BINDER')
            child.finish({'stage':'binder','execution_status':'COMPLETED','validation_decision':'REJECT','running':False,
                'summary':{'candidates':[{'candidate_id':'binder_fixture','final_verdict':'FAIL'}]}})
            record=RunRecord(Path(directory)/'upstream_fixture','UPSTREAM');record.data['question']='Fixture <question>'
            state={'design_mode':'all_eligible','cell_contexts':{'G':{}},
                'candidate_search':{'requested_limit':100,'evaluated_targets':1},
                'automatic_design_results':[{'gene':'G','modality':'DE_NOVO_BINDER','run_id':child.directory.name}]}
            with configuration_scope(cfg),patch.object(upstream_reporting,'plot_evidence',side_effect=ValueError('fixture plot error')):
                finish_report(record,state,{'execution_status':'PARTIAL'})
            page=(record.directory/'report.html').read_text()
            payload=json.loads(re.search(r'window.__COMPETITION_RUN__=(.*?);</script>',page,re.S).group(1))
            self.assertEqual(payload['design_report']['panels'][0]['tables'][0]['rows'][0]['candidate_id'],'binder_fixture')
            self.assertEqual(payload['design_report']['panels'][0]['report_url'],'../design_fixture/report.html')
            self.assertIn('renderPipelineDesign',page)
            self.assertNotIn("await fetch('/api/report' + location.search)",page)
            saved=read_report(directory,run_id=record.directory.name)
            self.assertEqual(saved['report_status'],'READY')
            self.assertIn('fixture plot error',saved['report_warnings'][0])
            self.assertEqual(saved['summary']['candidate_search']['requested_limit'],100)

    def test_unavailable_reference_keeps_current_report(self):
        with tempfile.TemporaryDirectory() as directory:
            report=assemble_pipeline_report({'summary':{'design_reference_run':'missing'}},directory)['design_report']
            self.assertIsNone(report['reference'])
            self.assertEqual(report['report_status'],'READY')
            self.assertIn('run_summary_missing',report['reference_warning'])


if __name__=='__main__':unittest.main()
