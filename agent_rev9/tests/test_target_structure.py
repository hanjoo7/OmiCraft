import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('_launch_target_structure',root/'launch.py')
launch=importlib.util.module_from_spec(spec);spec.loader.exec_module(launch);launch.load_package()
from agent_rev9 import target_structure as target, web_execution as web
from agent_rev9.configuration import OmiCraftConfig, configuration_scope

PDB=b'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \nATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 90.00           C  \nTER\nEND\n'

def metadata(**changes):
    return {'uniprotAccession':'Q12345','gene':'TEST','taxId':9606,'entryId':'AF-Q12345-F1',
            'uniprotSequence':'AG','uniprotStart':1,'uniprotEnd':2,'latestVersion':6,
            'pdbUrl':'https://alphafold.ebi.ac.uk/files/AF-Q12345-F1-model_v6.pdb',**changes}


class TargetStructure(unittest.TestCase):
    def test_download_validates_sequence_registers_and_reuses(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg=OmiCraftConfig();cfg.data.base_dir=directory
            with configuration_scope(cfg),patch.object(web,'ASSETS',Path(directory)/'missing.json'),patch.object(target,'_http_get_json',return_value=[metadata()]) as lookup,patch.object(target.urllib.request,'urlopen',return_value=io.BytesIO(PDB)) as download:
                first=target.ensure_target_structure('TEST',{'status':'FOUND','uniprot_id':'Q12345'})
                second=target.ensure_target_structure('TEST',{'status':'FOUND','uniprot_id':'Q12345'})
            self.assertEqual(first['id'],second['id'])
            self.assertEqual(first['sequence_end'],2)
            self.assertEqual(first['source_type'],'alphafold_db')
            self.assertEqual(download.call_count,1)
            self.assertEqual(lookup.call_count,1)
            self.assertTrue(Path(first['path']).is_file())

    def test_wrong_species_and_oversized_target_do_not_download(self):
        with patch.object(web,'assets',return_value=[]),patch.object(target.urllib.request,'urlopen') as download:
            for model in [metadata(taxId=10090),metadata(uniprotSequence='A'*1201,uniprotEnd=1201)]:
                with patch.object(target,'_http_get_json',return_value=[model]),self.assertRaises(ValueError):
                    target.ensure_target_structure('TEST',{'status':'FOUND','uniprot_id':'Q12345'})
            download.assert_not_called()

    def test_mismatched_download_is_not_registered(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg=OmiCraftConfig();cfg.data.base_dir=directory
            with configuration_scope(cfg),patch.object(web,'assets',return_value=[]),patch.object(target,'_http_get_json',return_value=[metadata(uniprotSequence='AA')]),patch.object(target.urllib.request,'urlopen',return_value=io.BytesIO(PDB)),self.assertRaisesRegex(ValueError,'sequence_mismatch'):
                target.ensure_target_structure('TEST',{'status':'FOUND','uniprot_id':'Q12345'})
            self.assertEqual(list(Path(directory).rglob('asset.json')),[])
            self.assertEqual(list(Path(directory).rglob('structure.pdb')),[])


if __name__=='__main__':unittest.main()
