import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('_launch_binder_analysis',root/'launch.py')
launch=importlib.util.module_from_spec(spec);spec.loader.exec_module(launch);launch.load_package()
from agent_rev9.binder_analysis import analyze_binders, candidate_metrics, fit_transform, rank_rows
from agent_rev9.binder_redesign import import_backbones
from agent_rev9.configuration import BinderAnalysisConfig, OmiCraftConfig, RFD3Config, configuration_scope
from agent_rev9.screening_agent import _screen_binders
from agent_rev9.rfdiffusion3_stage import RFdiffusion3Runner
from agent_rev9 import web_execution as web


def pdb(path, binder_sequence='GGGGG', shift=0, rotation=None, translation=None):
    names={'A':'ALA','C':'CYS','D':'ASP','G':'GLY','S':'SER'}
    rows=[];index=1
    for chain,seq,points in [('A','ACD',[[0,0,0],[4,0,0],[0,4,0]]),('B',binder_sequence,[[2,2,4],[6,2,4],[2,6,4],[6,6,4],[4,4,7]])]:
        for pos,(aa,xyz) in enumerate(zip(seq,points),1):
            ca=np.asarray(xyz,dtype=float)
            if chain=='B':ca[2]+=shift
            for atom,offset in [('N',[-1,0,0]),('CA',[0,0,0]),('C',[1,0,0]),('O',[1,1,0])]:
                coord=ca+offset
                if rotation is not None:coord=coord@rotation+translation
                x,y,z=coord
                rows.append(f'ATOM  {index:5d} {atom:^4s} {names[aa]} {chain}{pos:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 90.00           {atom[0]:>2s}\n');index+=1
    path.write_text(''.join(rows)+'END\n')
    return path


def fixture(directory, shift=2):
    ref=pdb(directory/'reference.pdb')
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    pred=pdb(directory/'prediction.pdb','SSSSS',shift,rotation,np.array([20,30,40]))
    pae=np.ones((8,8));pae[:3,3:]=4;pae[3:,:3]=8
    conf=directory/'confidences.json';conf.write_text(json.dumps({'token_chain_ids':['A']*3+['B']*5,'pae':pae.tolist()}))
    candidate={'candidate_id':'candidate_1','gene_name':'TEST','target_id':'target_1','rfd3_candidate_id':'backbone_1',
        'sequence':'SSSSS','binder_chain':'B','rfd3_structure_path':str(ref),'is_mock':False,'final_verdict':'PASS',
        'af3_metrics':{'model_path':str(pred),'confidences_path':str(conf),'iptm':.8,'ptm':.8,'plddt_summary':{'by_chain_mean':{'B':90}}}}
    context={'target_sequence':'ACD'}
    return candidate,context


class BinderAnalysis(unittest.TestCase):
    def test_target_aligned_rmsd_and_directional_pae(self):
        with tempfile.TemporaryDirectory() as temp:
            candidate,context=fixture(Path(temp))
            row=candidate_metrics(candidate,context,BinderAnalysisConfig())
            self.assertAlmostEqual(row['binder_rmsd'],2,places=5)
            self.assertAlmostEqual(row['binder_internal_rmsd'],0,places=5)
            self.assertAlmostEqual(row['target_rmsd'],0,places=5)
            self.assertEqual(row['i_pae_target_to_binder'],4)
            self.assertEqual(row['i_pae_binder_to_target'],8)
            self.assertEqual(row['i_pae'],6)
            self.assertEqual(row['quality_status'],'PASS')

    def test_sequence_mismatch_and_mock_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            candidate,context=fixture(Path(temp));candidate['sequence']='AAAAA'
            row=candidate_metrics(candidate,context,BinderAnalysisConfig())
            self.assertEqual(row['quality_status'],'NOT_EVALUATED')
            self.assertIn('AF3_sequence_mismatch',' '.join(row['errors']))
            candidate['is_mock']=True
            self.assertEqual(candidate_metrics(candidate,context,BinderAnalysisConfig())['quality_status'],'NOT_EVALUATED')
        with self.assertRaises(ValueError):fit_transform(np.zeros((3,3)),np.zeros((3,3)))

    def test_ranking_and_best_bundle_preserve_primary_verdict(self):
        with tempfile.TemporaryDirectory() as temp:
            candidate,context=fixture(Path(temp));duplicate=copy.deepcopy(candidate);duplicate['candidate_id']='candidate_2';duplicate['final_verdict']='FAIL'
            screening={'candidates':[candidate,duplicate],'targets':[{'target_id':'target_1','design_conditions':context}]}
            out=analyze_binders(screening,Path(temp)/'analysis',BinderAnalysisConfig())
            self.assertEqual(out['unique_sequence_count'],1)
            self.assertEqual(out['best_candidate_id'],'candidate_1')
            self.assertEqual(out['ranking'][1]['duplicate_of'],'candidate_1')
            self.assertEqual(duplicate['final_verdict'],'FAIL')
            self.assertTrue(Path(out['best_structure_path']).is_file())
            with zipfile.ZipFile(out['archive_path']) as archive:
                self.assertIn('ranking.csv',archive.namelist());self.assertIn('best.json',archive.namelist())

    def test_failed_additional_filter_does_not_change_primary_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            candidate,context=fixture(Path(temp))
            out=analyze_binders({'candidates':[candidate],'targets':[{'target_id':'target_1','design_conditions':context}]},Path(temp)/'analysis',BinderAnalysisConfig(ipae_max=5))
            self.assertEqual(candidate['final_verdict'],'PASS')
            self.assertEqual(out['ranking'][0]['quality_status'],'FAIL')
            self.assertIsNone(out['best_structure_path'])

    def test_reanalysis_removes_obsolete_best_and_exports_primary_pass_duplicate(self):
        with tempfile.TemporaryDirectory() as temp:
            directory=Path(temp);candidate,context=fixture(directory)
            duplicate=copy.deepcopy(candidate);duplicate['candidate_id']='candidate_0';duplicate['final_verdict']='FAIL'
            screening={'candidates':[candidate,duplicate],'targets':[{'target_id':'target_1','design_conditions':context}]}
            out=analyze_binders(screening,directory/'analysis',BinderAnalysisConfig())
            self.assertEqual(out['best_candidate_id'],'candidate_1')
            self.assertEqual(out['ranking'][1]['duplicate_of'],'candidate_0')
            out=analyze_binders(screening,directory/'analysis',BinderAnalysisConfig(ipae_max=5))
            self.assertFalse((directory/'analysis/best.json').exists())
            self.assertFalse((directory/'analysis/best.pdb').exists())
            with zipfile.ZipFile(out['archive_path']) as archive:
                self.assertNotIn('best.json',archive.namelist())
                self.assertNotIn('best.pdb',archive.namelist())

    def test_redesign_uses_supplied_backbone_and_skips_rfd3(self):
        with tempfile.TemporaryDirectory() as temp:
            directory=Path(temp);source=pdb(directory/'source.pdb')
            cfg=OmiCraftConfig();cfg.data.base_dir=temp
            settings=RFD3Config(enabled=True,backbone_paths=(str(source),),target_structure=str(source),target_sequence='ACD',target_chain='A',binder_chain='B',design_length=5,num_designs=1)
            with configuration_scope(cfg),patch.object(RFdiffusion3Runner,'run',side_effect=AssertionError('RFD3 must not run')):
                output,errors,_=_screen_binders({'use_rfd3':True,'dry_run':True,'protein_design_input':settings.model_dump()},[{'gene_name':'TEST'}],[],{})
            self.assertFalse(errors)
            self.assertEqual(output['design_protocol'],'sequence_redesign')
            self.assertEqual(output['provided_backbones'],1)
            self.assertEqual(output['generated_rfd3'],0)
            self.assertEqual(output['requested_rfd3'],0)
            self.assertTrue(all(c['is_mock'] for c in output['candidates']))

    def test_redesign_api_uses_registered_asset_and_sequence_budget(self):
        with tempfile.TemporaryDirectory() as temp:
            source=pdb(Path(temp)/'source.pdb')
            row={'id':'complex','gene':'ERBB2','path':str(source)}
            with patch.object(web,'asset',return_value=row),patch.object(web,'model_config',return_value={'rfd3':{},'af3_binder':{}}):
                data,config,_=web.design_input({'parent_run_id':'reference_erbb2','gene':'ERBB2','modality':'DE_NOVO_BINDER',
                    'input':{'asset_id':'complex','sequence_redesign':True,'target_chain':'A','binder_chain':'B','num_sequences':3}})
            self.assertEqual(data['protein_design_input']['backbone_paths'],[str(source)])
            self.assertEqual(config['protein_mpnn']['num_sequences_per_backbone'],3)
            self.assertEqual(config['screening_options']['max_af3_candidates'],3)


if __name__=='__main__':unittest.main()
