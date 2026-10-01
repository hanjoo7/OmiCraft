"""Demo contracts use temporary fixtures; no heavy model execution."""

import json
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure
from matplotlib import pyplot as plt

from .. import modality_dispatch as dispatch
from .. import run_records
from ..antibody_mapping import conjugation_feasibility
from ..configuration import CompetitionDemoConfig
from ..critic_agent import review_demo_modality
from ..degrader_execution import execute as assemble
from ..orchestrator import run_screening
from ..run_records import RunRecord
from ..server import RUN_REPORT_SCRIPT
from ..user_selection import create_selection_event
from ..validation_agent import binder_developability, validate_demo_modality


@pytest.fixture
def source(tmp_path):
    root = tmp_path / 'real_source'
    root.mkdir()
    (root / 'run_summary.json').write_text(json.dumps({'is_mock': False, 'summary': {'verdict': 'VERIFIED_REAL_PARTIAL'}}))
    return root


@pytest.fixture
def replay_spy(monkeypatch):
    calls = []
    def replay(modality, source, record, config):
        calls.append(modality)
        row = dispatch.result(modality, 'BLOCKED' if modality == 'DEGRADER' else 'PARTIAL', 'fixture',
            blockers=['attachment_missing'] if modality == 'DEGRADER' else [])
        row.update(demo_metrics={}, real_tool_execution_status='REAL_ARTIFACT_REPLAY',
                   execution_mode='REAL_ARTIFACT_REPLAY', missing_steps=['experimental_evidence'],
                   source_provenance=[], lightweight_calculations=[])
        return row
    monkeypatch.setattr(dispatch, 'replay_demo_modality', replay)
    monkeypatch.setattr(run_records, 'demo_visualizations', lambda *a: [])
    return calls


def test_demo_complete_with_hold_and_blocked(source, replay_spy, tmp_path):
    state = {'execution_profile': 'competition_demo', 'output_dir': str(tmp_path / 'demo'),
             'competition_demo': {'source_run': str(source), 'selection_mode': 'showcase_all'}}
    result = run_screening(state)
    assert result['demo_pipeline_status'] == 'COMPLETE'
    assert result['scientific_validation_status'] == 'VERIFIED_REAL_PARTIAL'
    assert result['clinical_or_scientific_selection'] is False
    assert replay_spy == list(dispatch.DEMO_TOOL_CHAINS)
    assert result['summary']['heavy_model_inference_calls'] == 0
    for row in result['summary']['modalities']:
        expected = (['Curated reference selection', 'RDKit full molecule QC', 'ETKDG', 'MMFF/UFF', 'Ternary readiness']
                    if row['modality'] == 'DEGRADER' else dispatch.DEMO_TOOL_CHAINS[row['modality']])
        assert row['selected_tools'] == expected
        assert row['execution_mode'] == 'REAL_ARTIFACT_REPLAY'
        assert row['is_mock'] is False
        assert row['critic_decision'] in {'HOLD', 'NOT_EVALUATED'}
        stages = [s['stage'] for s in result['stages'] if s.get('modality') == row['modality']]
        assert stages == ['modality_evaluation', 'tool_planning', 'readiness_check', 'execution_or_replay',
                          'validation', 'critic', 'visualization', 'report_registration']
    assert result['summary']['modalities'][-1]['modality_execution_status'] == 'BLOCKED'
    html = Path(result['report_path']).read_text()
    assert 'window.__COMPETITION_RUN__=' in html and 'renderCompetitionDemo' in html
    assert "await fetch('/api/report' + location.search)" not in html
    assert json.loads((source / 'run_summary.json').read_text())['summary']['verdict'] == 'VERIFIED_REAL_PARTIAL'


def test_demo_explicit_selection_is_not_showcase(source, replay_spy, tmp_path):
    state = {'execution_profile': 'competition_demo', 'output_dir': str(tmp_path / 'unselected'),
        'competition_demo': {'source_run': str(source)}, 'advance_targets': [{'gene_name': 'ERBB2', 'tier': '1A'}]}
    result = run_screening(state)
    assert not replay_spy and all(r['execution_status'] == 'NOT_RUN' for r in result['summary']['modalities'])
    event = create_selection_event('ERBB2', 'ADC', 'ERBB2:ADC', '1A')
    run_screening({**state, 'output_dir': str(tmp_path / 'selected'), 'user_selection_events': [event]})
    assert replay_spy == ['ADC']
    assert dispatch.DEMO_TOOL_CHAINS['ADC'][0] == 'ANARCII'
    assert dispatch.DEMO_TOOL_CHAINS['DEGRADER'][0] == 'RDKit component QC'


def test_replay_exception_still_reaches_validation_critic_report(source, tmp_path, monkeypatch):
    def broken(*args):
        raise ValueError('checksum_mismatch')
    monkeypatch.setattr(dispatch, 'replay_demo_modality', broken)
    result = dispatch.run_competition_demo({'output_dir': tmp_path / 'blocked'},
        CompetitionDemoConfig(source_run=str(source), selection_mode='showcase_all'))
    assert result['demo_pipeline_status'] == 'COMPLETE'
    assert all(r['execution_status'] == 'BLOCKED' for r in result['summary']['modalities'])
    assert all(r['critic_decision'] in {'HOLD', 'NOT_EVALUATED'} for r in result['summary']['modalities'])
    assert not list((tmp_path / 'blocked').rglob('*.png'))


def test_mock_source_is_rejected(source, tmp_path):
    (source / 'run_summary.json').write_text(json.dumps({'is_mock': True, 'summary': {'verdict': 'VERIFIED_REAL_PARTIAL'}}))
    with pytest.raises(ValueError, match='verified_real_partial'):
        dispatch.run_competition_demo({'output_dir': tmp_path / 'bad'}, CompetitionDemoConfig(source_run=str(source)))


def test_binder_developability_rules_and_critic():
    development = binder_developability('AAAAAVVVVVNASC')
    assert development['status'] == 'PASS_WITH_WARNINGS'
    assert development['hydrophobic_stretches'] and development['cysteine_count'] == 1
    assert any(m['kind'] == 'N_linked_glycosylation' for m in development['liability_motifs'])
    assert development['estimated_pI'] is not None
    assert development['experimental_binding'] == 'NOT_EVALUATED'
    row = dispatch.result('DE_NOVO_BINDER', 'COMPLETED', 'test',
        summary={'candidates': [{'final_verdict': 'PASS', 'is_mock': False}]})
    row.update(demo_metrics={'developability': development, 'structural_screening': 'PASS'},
               missing_steps=['experimental_developability'])
    row['validation'] = validate_demo_modality(row)
    assert review_demo_modality(row)['decision'] == 'PASS_WITH_WARNINGS'
    assert binder_developability('XXXX')['status'] == 'HOLD'


def test_adc_candidate_mapping_and_plot_sources(tmp_path):
    structure = Structure.Structure('test')
    model = Model.Model(0)
    structure.add(model)
    mapping = []
    for chain_id, x in [('H', 0), ('L', 10), ('C', 20)]:
        chain = Chain.Chain(chain_id)
        model.add(chain)
        residue = Residue.Residue((' ', 1, ' '), 'LYS', '')
        chain.add(residue)
        for i, name in enumerate(['CA', 'NZ']):
            residue.add(Atom.Atom(name, np.array([x, i*2., 0.]), 80., 1., ' ', name, i, element='C' if name=='CA' else 'N'))
        if chain_id != 'C':
            mapping.append({'chain_id': chain_id, 'structure_residue_id': 1, 'insertion_code': '',
                            'residue_name': 'LYS', 'region': 'CDR1' if chain_id=='H' else 'framework',
                            'original_sequence_index': 1})
    path = tmp_path / 'adc.pdb'
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))
    analysis = {'structure_path': str(path), 'mapping': mapping, 'contacts': []}
    feasibility = conjugation_feasibility(analysis)
    assert len(feasibility['candidates']) == 2
    assert all(r['selected'] is False for r in feasibility['candidates'])
    assert next(r for r in feasibility['candidates'] if r['chain']=='H')['exclusion_reasons'] == ['inside_CDR']
    record = RunRecord(tmp_path / 'plots', 'ADC')
    plots = run_records.demo_visualizations(record, {'modality': 'ADC', 'demo_metrics': {
        'mapping': mapping, 'cdr_contacts': [{'region': 'CDR1', 'heavy_atom_contacts': 1}], 'conjugation': feasibility}})
    assert len(plots) == 3
    assert not plt.get_fignums()
    for plot in plots:
        assert Path(plot['path']).is_file() and Path(plot['source_table']).is_file()
        assert plot['is_mock'] is False


def test_degrader_explicit_assembly_and_missing_maps(tmp_path):
    data = {'warhead_smiles': 'CC[*:1]', 'e3_ligand_smiles': 'CC[*:2]',
            'linker_smiles': '[*:3]CC[*:4]', 'warhead_map': 1, 'e3_map': 2, 'linker_maps': [3, 4], 'expected_full_smiles': 'CCCCCC'}
    directory = tmp_path / 'assembly'
    directory.mkdir()
    out = assemble(data, {'output_dir': directory})
    assert out['summary']['sanitization_status'] == 'success'
    assert out['summary']['atom_count_after_hydrogens'] > out['summary']['atom_count_before_hydrogens']
    assert 'nonbonded_heavy_atom_clashes' in out['summary']
    missing = assemble({**data, 'warhead_map': None}, {'output_dir': tmp_path / 'missing'})
    assert missing['execution_status'] == 'BLOCKED'
    assert not (tmp_path / 'missing/degrader.sdf').exists()


def test_reference_demo_label_and_checksum(source, tmp_path):
    directory = source / 'degrader'
    directory.mkdir()
    summary = dispatch.result('DEGRADER', 'BLOCKED', 'inputs')
    (directory / 'run_summary.json').write_text(json.dumps(summary))
    (directory / 'input.json').write_text(json.dumps({'warhead_smiles': 'CCO'}))
    config = tmp_path / 'reference.json'
    config.write_text(json.dumps({'is_mock': False, 'target_context': 'BRD4–VHL reference demo',
        'warhead_smiles': 'CC[*:1]', 'e3_ligand_smiles': 'CC[*:2]', 'linker_smiles': '[*:3]CC[*:4]',
        'warhead_map': 1, 'e3_map': 2, 'linker_maps': [3, 4], 'expected_full_smiles': 'CCCCCC'}))
    record = RunRecord(tmp_path / 'reference', 'DEGRADER')
    out = dispatch.replay_demo_modality('DEGRADER', source, record, CompetitionDemoConfig(
        source_run=str(source), degrader_component_config=str(config)))
    assert out['execution_mode'] == 'REFERENCE_DEMO'
    assert out['demo_metrics']['ERBB2_specific_design'] == 'NOT_EVALUATED'
    assert out['demo_metrics']['target_context'] == 'BRD4–VHL reference demo'
    (source / 'events.jsonl').write_text(json.dumps({'output_sha256': {str((directory/'run_summary.json').resolve()): 'incorrect'}})+'\n')
    with pytest.raises(ValueError, match='checksum_mismatch'):
        dispatch.replay_demo_modality('DEGRADER', source, record, CompetitionDemoConfig(source_run=str(source)))


def test_report_has_visible_core_and_closed_details():
    assert 'Competition demo pipeline:' in RUN_REPORT_SCRIPT
    assert 'Scientific validation:' in RUN_REPORT_SCRIPT
    assert "add(parent,'details')" in RUN_REPORT_SCRIPT
    assert '.open=true' not in RUN_REPORT_SCRIPT
    assert 'All mapped conjugation candidates' in RUN_REPORT_SCRIPT
    assert 'is_mock=false' in RUN_REPORT_SCRIPT
    assert 'ProteinMPNN score:' in RUN_REPORT_SCRIPT
    assert 'not binding affinity' not in dispatch.DEMO_TOOL_CHAINS['DE_NOVO_BINDER']


def reference_config():
    from ..modality_dispatch import select_degrader_reference
    return select_degrader_reference({'target': 'ERBB2'})


def test_published_reference_qc_and_selection():
    from ..degrader_execution import full_molecule_qc, input_blockers
    data = reference_config()
    assert data['candidate_id'] == 'SJF1528'
    assert data['candidate_origin'] == 'published_reference' and data['novel_design'] is False
    assert data['target_context'] == 'ERBB2/EGFR'
    assert input_blockers(data) == []
    assert all(data[k] is None for k in ('warhead_map', 'e3_map', 'linker_maps'))
    mol, qc = full_molecule_qc(data)
    assert qc['original_smiles'] == data['full_degrader_smiles']
    assert qc['canonical_smiles'] != qc['isomeric_smiles']
    assert qc['stereocenter_count'] == 3
    assert qc['inchikey'] == 'ZSCOIFSUFMYZEZ-YSWDPXALSA-N'
    assert qc['property_profile'] == 'beyond_rule_of_five'
    assert mol.GetNumAtoms() == 73
    assert all(row['mode'] == 'RUNTIME_SELECTION' for row in data['selection_trace'])


@pytest.mark.parametrize('missing', ['full_degrader_smiles', 'primary_doi', 'source_sha256', 'primary_pmid'])
def test_reference_required_evidence(missing, tmp_path):
    data = reference_config()
    data.pop(missing)
    output = assemble(data, {'output_dir': tmp_path})
    assert output['execution_status'] == 'BLOCKED'
    assert 'missing:' + missing in output['blockers']
    assert not list(tmp_path.iterdir())


def test_reference_source_mismatch_and_stereo_rejected():
    from ..degrader_execution import full_molecule_qc
    data = reference_config()
    with pytest.raises(ValueError, match='MISMATCH'):
        full_molecule_qc({**data, 'full_degrader_smiles': data['full_degrader_smiles'].replace('[C@H]', '[C@@H]', 1)})
    with pytest.raises(ValueError, match='STEREOCHEMISTRY'):
        full_molecule_qc({**data, 'full_degrader_smiles': data['full_degrader_smiles'].replace('@', '')})


def test_component_expected_structure_is_required(tmp_path):
    data = {'warhead_smiles': 'CC[*:1]', 'e3_ligand_smiles': 'CC[*:2]',
            'linker_smiles': '[*:3]CC[*:4]', 'warhead_map': 1, 'e3_map': 2, 'linker_maps': [3, 4]}
    assert assemble(data, {'output_dir': tmp_path})['execution_status'] == 'BLOCKED'
    failed = assemble({**data, 'expected_full_smiles': 'CCCCC'}, {'output_dir': tmp_path})
    assert failed['execution_status'] == 'FAILED'
    assert failed['summary']['assembly_status'] == 'FAILED'
    assert not list(tmp_path.iterdir())
    for key in ('warhead_map', 'e3_map', 'linker_maps'):
        assert assemble({**data, 'expected_full_smiles': 'CCCCCC', key: None},
                        {'output_dir': tmp_path})['execution_status'] == 'BLOCKED'


def test_reference_conformer_without_assembly_or_af3(tmp_path, monkeypatch):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

    from ..critic_agent import review_modality
    from ..validation_agent import validate_modality
    data = reference_config()
    # A small chiral fixture exercises preparation without adding an expensive SJF conformer to unit tests.
    smiles = 'N[C@@H](C)C(=O)O'
    mol = Chem.MolFromSmiles(smiles)
    data.update(candidate_id='chiral_test_fixture', full_degrader_smiles=smiles,
        structure_crosscheck=dict(isomeric_smiles=smiles, inchikey=Chem.MolToInchiKey(mol),
            molecular_formula=rdMolDescriptors.CalcMolFormula(mol), molecular_weight=Descriptors.MolWt(mol), stereocenter_count=1))
    monkeypatch.setattr(AllChem, 'MMFFHasAllMoleculeParams', lambda mol: False)
    out = assemble(data, {'output_dir': tmp_path})
    assert out['execution_status'] == 'PARTIAL', out  # No ternary prediction was run.
    assert out['summary']['assembly_status'] == 'NOT_APPLICABLE_REFERENCE_INPUT'
    assert out['summary']['assembly_performed'] is False
    assert out['summary']['qc']['force_field'] == 'UFF'
    assert out['summary']['qc']['conformer_identity_preserved'] is True
    out['validation'] = validate_modality(out)
    out['validation_decision'] = out['validation']['decision']
    assert out['validation']['full_molecule_valid'] and out['validation']['conformer_valid']
    assert out['validation']['ternary_prediction_performed'] is False
    assert review_modality(out)['decision'] == 'HOLD'  # Chemistry alone cannot pass structural review.
    assert (tmp_path / 'degrader_2d.png').stat().st_size > 0
    assert (tmp_path / 'degrader_3d.png').stat().st_size > 0
    assert not list(tmp_path.glob('ternary*.png'))
    assert not (tmp_path / 'af3_input.json').exists()


def test_reference_report_labels_and_prior_attempt():
    for text in ('PUBLISHED DEGRADER', 'CURATED REFERENCE CANDIDATE',
                 'NOT A NOVEL AGENT-DESIGNED MOLECULE', 'Previous attempt: TAK-285-only input',
                 'Component decomposition: NOT VERIFIED', 'Ternary readiness'):
        assert text in RUN_REPORT_SCRIPT
    assert 'selective HER2 degrader' not in RUN_REPORT_SCRIPT


def test_reference_integration_preserves_other_results(tmp_path):
    from ..run_records import integrate_degrader_reference
    report = tmp_path / 'visual_diagnostic_test'
    reference = tmp_path / 'reference'
    reference.mkdir()
    record = RunRecord(report, 'FOUR_MODALITY_DEMO')
    other = dispatch.result('SMALL_MOLECULE', 'PARTIAL', 'pose_qc', summary={'candidate_id': 'TAK-285'})
    old = dispatch.result('DEGRADER', 'BLOCKED', 'input', blockers=['missing:e3_ligand_smiles'])
    original = {'summary': {'modalities': [other, old], 'issues': [], 'attempts': [], 'omics': ['preserved']},
                'execution_profile': 'visual_diagnostic', 'scientific_validation_status': 'VERIFIED_REAL_PARTIAL',
                'demo_pipeline_status': 'COMPLETE'}
    record.data.update(original)
    (report / 'run_summary.json').write_text(json.dumps(record.data))
    new = dispatch.result('DEGRADER', 'COMPLETED', 'reference_evaluation', summary={
        'candidate_origin': 'published_reference', 'full_molecule_validation': 'COMPLETED',
        'ternary_status': 'BLOCKED_ADAPTER_UNSUPPORTED', 'ternary_readiness': {'reason': 'unsupported'}})
    new.update(run_id=reference.name, validation={'decision': 'HOLD'}, critic={'decision': 'PASS_WITH_WARNINGS'},
               stages=[{'timestamp': '2026-01-01T00:00:00Z'}])
    (reference / 'run_summary.json').write_text(json.dumps(new))
    result = integrate_degrader_reference(report, reference)
    assert result['summary']['modalities'][0] == other
    assert result['summary']['omics'] == ['preserved']
    assert result['summary']['modalities'][1]['previous_attempt'] == old
    assert json.loads((reference / 'previous_degrader_attempt.json').read_text()) == old
    manifest = json.loads((report / 'report_manifest.json').read_text())
    assert manifest['degrader_reference_run_id'] == reference.name
    assert manifest['figures'] == []
    assert not list(report.glob('*ternary*.png'))
