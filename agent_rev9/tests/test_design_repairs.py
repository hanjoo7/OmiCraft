import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('_launch_test', root / 'launch.py')
launch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launch)
launch.load_package()

from agent_rev9 import web_execution as web, modality_dispatch as dispatch
from agent_rev9 import server, screening_agent, antibody_mapping
from agent_rev9.af3_ternary_tool import AF3TernaryTool
from agent_rev9.degrader_execution import ternary_readiness
from agent_rev9.report_history import report_history
from agent_rev9.small_molecule_config import LigandAF3Config
from agent_rev9.structure_io import protein_residues


class DesignRepairs(unittest.TestCase):
    def test_reference_numbering(self):
        data, config, _ = web.design_input(dict(parent_run_id='reference_erbb2', gene='ERBB2', modality='DE_NOVO_BINDER'))
        design = data['protein_design_input']
        self.assertEqual(design['hotspot_residues'], ['A953'])
        residues = protein_residues(design['target_structure'])
        self.assertIn(('A', '953', 'M'), residues)
        self.assertEqual(''.join(aa for _, _, aa in residues), design['target_sequence'])
        self.assertFalse(dispatch.readiness('DE_NOVO_BINDER', data, config)['blockers'])
        design['hotspot_residues'] = ['A252']
        self.assertIn('Selection is not a target residue: A252', dispatch.readiness('DE_NOVO_BINDER', data, config)['blockers'])

    def test_binder_failure_reaches_blockers(self):
        output = dict(screening_result=dict(status='NOT_CONFIGURED', targets=[dict(failure_reasons=['invalid residue'])]), errors=[])
        with tempfile.TemporaryDirectory() as directory, patch.object(screening_agent, 'screening_node', return_value=output):
            result = dispatch.execute_binder({}, {'output_dir': directory})
        self.assertEqual(result['execution_status'], 'BLOCKED')
        self.assertEqual(result['blockers'], ['invalid residue'])

    def test_binder_scientific_rejection_is_completed(self):
        output = dict(screening_result=dict(status='COMPLETED', candidates=[dict(final_verdict='FAIL', is_mock=False)], targets=[dict(failure_reasons=['low iPTM'])]))
        with tempfile.TemporaryDirectory() as directory, patch.object(screening_agent, 'screening_node', return_value=output):
            result = dispatch.execute_binder({}, {'output_dir': directory})
        self.assertEqual(result['execution_status'], 'COMPLETED')
        self.assertFalse(result['blockers'])

    def test_adc_conjugation_is_connected_without_full_completion_claim(self):
        analysis = dict(structure_path='input.cif', cdr_mapping_status='success', mapping=[], contacts=[], blockers=[])
        conjugation = dict(status='COMPLETED', candidates=[dict(chain='H', residue_id=7, sasa_angstrom2=12)])
        with tempfile.TemporaryDirectory() as directory, patch.object(antibody_mapping, 'analyze', return_value=analysis), patch.object(antibody_mapping, 'conjugation_feasibility', return_value=conjugation) as analyze:
            result = dispatch.execute_adc({}, {'output_dir': directory})
            self.assertTrue((Path(directory) / 'conjugation_candidates.csv').exists())
        analyze.assert_called_once()
        self.assertEqual(result['execution_status'], 'COMPLETED')
        self.assertFalse(result['blockers'])
        self.assertIn('internalization', result['missing_steps'])

    def test_ternary_actual_input_and_dry_run(self):
        data, config, _ = web.design_input(dict(parent_run_id='reference_erbb2', gene='ERBB2', modality='DEGRADER'))
        self.assertEqual(ternary_readiness(data, config)['status'], 'READY')
        model = LigandAF3Config(**config['af3'])
        context = dict(target_sequence=config['ternary_sequences']['ERBB2']['sequence'],
                       partner_proteins=[config['ternary_sequences'][name] for name in ['VHL','ElonginB','ElonginC']])
        candidate = dict(candidate_id=data['candidate_id'], canonical_smiles=data['full_degrader_smiles'])
        with tempfile.TemporaryDirectory() as directory, patch('subprocess.run', side_effect=AssertionError('Inference must not run')):
            result = AF3TernaryTool(model).run(candidate=candidate, context=context, output_dir=directory, dry_run=True)
            payload = json.loads(Path(result['input_json']).read_text())
            self.assertEqual(result['status'], 'not_run')
            self.assertEqual([next(iter(row.values()))['id'] for row in payload['sequences']], ['A','B','D','E','C'])
            self.assertEqual(payload['sequences'][-1]['ligand']['smiles'], data['full_degrader_smiles'])
            context['partner_proteins'][0] = {**context['partner_proteins'][0], 'chain':'C'}
            result = AF3TernaryTool(model).run(candidate=candidate, context=context, output_dir=directory, dry_run=True)
            self.assertEqual(result['status'], 'invalid_input')

    def test_binary_prediction_cannot_pass_as_ternary(self):
        from agent_rev9.af3_ligand_tool import AF3LigandTool
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'summary.json'
            path.write_text(json.dumps({'chain_ids': ['A', 'C']}))
            tool = AF3TernaryTool(LigandAF3Config())
            tool.expected_chain_ids = {'A', 'B', 'C'}
            with patch.object(AF3LigandTool, '_sample', return_value={'af3_status': 'success', 'summary_path': str(path)}):
                result = tool._sample(None, 'test', 1, 0)
            self.assertEqual(result['af3_status'], 'backend_failed')
            self.assertEqual(result['error_message'], 'TERNARY_CHAIN_IDS_MISSING_OR_MISMATCH')

    def test_all_blocked_batch_has_four_child_reports(self):
        cfg = web.get_config().model_copy(deep=True)
        with tempfile.TemporaryDirectory() as directory:
            cfg.data.base_dir = directory
            name = 'design_20260924T130000Z_12345678'
            with patch.object(web, 'design_input', side_effect=ValueError('missing required input')), patch.object(server, 'write_competition_report'), patch.object(server, 'run_state', dict(history=[])):
                web.run_erbb2_batch(name, dict(parent_run_id='reference_erbb2', gene='ERBB2', modality='ALL'), cfg)
            report = json.loads((Path(directory)/name/'run_summary.json').read_text())
            self.assertEqual(len(report['summary']['modalities']),4)
            for row in report['summary']['modalities']:
                child = json.loads((Path(directory)/row['run_id']/'run_summary.json').read_text())
                self.assertEqual(child['execution_status'],'BLOCKED')
                self.assertEqual(child['blockers'],['missing required input'])

    def test_history_excludes_unapproved_names_and_outside_symlinks(self):
        with tempfile.TemporaryDirectory() as directory, tempfile.TemporaryDirectory() as outside:
            root = Path(directory)
            for name in ['design_20260924T130000Z_12345678', 'secrets']:
                (root/name).mkdir()
                (root/name/'run_summary.json').write_text(json.dumps(dict(execution_profile='therapeutic_design',modality='ALL')))
            (Path(outside)/'run_summary.json').write_text('{}')
            (root/'design_20260924T140000Z_12345678').symlink_to(outside)
            self.assertEqual(len(report_history(root,'reference')['runs']),1)


if __name__ == '__main__':
    unittest.main()
