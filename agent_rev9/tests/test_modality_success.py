import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..af3_ligand_tool import AF3LigandTool, confidence_chain_ids
from ..af3_ternary_tool import AF3TernaryTool, ternary_geometry
from ..modality_validation import structural_assessment
from ..validation_agent import validate_modality
from ..critic_agent import review_modality
from .test_small_molecule_route import ligand_case, run_case


def test_repeated_token_chain_ids_preserve_real_chain_metrics(tmp_path, ligand_case):
    candidate = run_case(tmp_path, ligand_case)['candidates'][0]
    sample = candidate['af3_ligand']['samples'][0]
    path = Path(sample['summary_path'])
    data = json.loads(path.read_text())
    data['chain_ids'] = ['A'] * 286 + ['B'] * 38
    path.write_text(json.dumps(data))
    parsed = AF3LigandTool(ligand_case[1].af3)._sample(path.parent, 'fixture', 11, 0)
    assert parsed['ligand_chain_iptm'] == .8
    assert parsed['protein_ligand_pae_min'] == 2
    assert parsed['af3_status'] == 'success'


@pytest.mark.parametrize('data', [
    {'chain_ids': ['A', 'B'], 'chain_iptm': [.9]},
    {'chain_ids': ['A', 'B'], 'chain_pair_iptm': [[.9, .8]]},
    {'chain_ids': ['A', None]},
])
def test_chain_dimension_errors_fail_closed(data):
    with pytest.raises(ValueError):
        confidence_chain_ids(data)


def test_chain_fallback_and_set_mismatch():
    assert confidence_chain_ids({'chain_iptm': [.5, .6]}, {'token_chain_ids': ['A', 'A', 'B']}) == ['A', 'B']
    with pytest.raises(ValueError):
        confidence_chain_ids({'chain_ids': ['A', 'B']}, {'token_chain_ids': ['A', 'C']})


def adc():
    return {'modality': 'ADC', 'execution_status': 'COMPLETED', 'is_mock': False, 'blockers': [],
            'missing_steps': ['internalization', 'experimental_efficacy'],
            'summary': {'cdr_mapping_status': 'success', 'cdr_heavy_atom_contacts': 12, 'severe_clashes': 0,
                        'conjugation': {'status': 'COMPLETED', 'candidates': [
                            {'inside_cdr': False, 'at_antigen_interface': False, 'possible_disulfide': False,
                             'exclusion_reasons': [], 'sasa_angstrom2': 25}]}}}


def degrader():
    sample = {'seed': 11, 'af3_status': 'success', 'protein_chain_id': 'A', 'ligand_chain_id': 'C',
              'has_clash': False, 'chirality_valid': True, 'iptm': .8,
              'ternary_geometry': {'status': 'COMPLETED', 'severe_clash_count': 0, 'ligand_contacts_by_chain': {'A': 12, 'B': 9}},
              'chain_pair_metrics': [{'chain_a': c, 'chain_b': 'C', 'chain_pair_iptm': .7, 'chain_pair_pae_min': 3} for c in ['A', 'B']]}
    return {'modality': 'DEGRADER', 'execution_status': 'COMPLETED', 'is_mock': False, 'blockers': [],
            'missing_steps': ['experimental_degradation'], 'summary': {'input_mode': 'published_full_molecule',
            'full_molecule_validation': 'COMPLETED', 'provenance': {'curation_status': 'VERIFIED'},
            'qc': {'conformer_generation_status': 'COMPLETED', 'conformer_identity_preserved': True},
            'ternary_prediction': {'status': 'success', 'samples': [sample, {**copy.deepcopy(sample), 'seed': 23}]}}}


@pytest.mark.parametrize('fixture', [adc, degrader])
def test_structural_pass_keeps_experimental_evidence_missing(fixture):
    output = fixture()
    validation = validate_modality(output)
    assert validation['decision'] == 'ADVANCE'
    assert validation['experimental_validation'] == 'NOT_EVALUATED'
    assert validation['missing_evidence'] == output['missing_steps']
    output.update(validation=validation, validation_decision=validation['decision'])
    critic = review_modality(output)
    assert critic['decision'] == 'ADVANCE'
    assert critic['experimental_efficacy_established'] is False
    output['is_mock'] = True
    assert validate_modality(output)['decision'] == 'HOLD'


def test_adc_no_sites_or_no_contact_cannot_advance():
    result = adc()
    result['summary']['conjugation']['candidates'][0]['inside_cdr'] = True
    assert structural_assessment(result)['decision'] == 'HOLD'
    result = adc();result['summary']['cdr_heavy_atom_contacts'] = 0
    assert structural_assessment(result)['decision'] == 'HOLD'
    result['summary']['severe_clashes'] = 1
    assert structural_assessment(result)['decision'] == 'REJECT'


@pytest.mark.parametrize('failure', ['single_seed', 'duplicate_seed', 'missing_metric', 'low_confidence', 'no_e3_contacts', 'chirality'])
def test_degrader_missing_or_bad_evidence_never_advances(failure):
    result = degrader();samples = result['summary']['ternary_prediction']['samples']
    if failure == 'single_seed': samples.pop()
    elif failure == 'duplicate_seed': samples[1]['seed'] = 11
    elif failure == 'missing_metric': samples[0]['chain_pair_metrics'][0]['chain_pair_iptm'] = None
    elif failure == 'low_confidence': samples[0]['iptm'] = .1
    elif failure == 'no_e3_contacts': samples[0]['ternary_geometry']['ligand_contacts_by_chain']['B'] = 0
    else: samples[0]['chirality_valid'] = False
    assert structural_assessment(result)['decision'] != 'ADVANCE'


def test_ternary_geometry_requires_both_partners(monkeypatch):
    from .. import screening_gates
    atoms = [SimpleNamespace(label_chain_id=chain, coord=xyz, is_hydrogen=False)
             for chain, xyz in [('A', (0,0,0)), ('C', (3,0,0)), ('B', (6,0,0))]]
    monkeypatch.setattr(screening_gates, 'cif_atoms', lambda path: atoms)
    result = ternary_geometry('fixture', 'C')
    assert result['ligand_contacts_by_chain'] == {'A': 1, 'B': 1}
    assert result['severe_clash_count'] == 0






@pytest.mark.parametrize('fallback', [False, True])
def test_ternary_repeated_ids_produce_unique_chain_pairs(tmp_path, monkeypatch, fallback):
    from ..small_molecule_config import LigandAF3Config
    path = tmp_path / 'summary.json'
    path.write_text(json.dumps({'chain_ids': ['A', 'A', 'B', 'B', 'C'],
                               'chain_pair_iptm': [[.9,.7,.8],[.7,.9,.6],[.8,.6,.9]],
                               'chain_pair_pae_min': [[1,2,3],[2,1,4],[3,4,1]]}))
    details = tmp_path / 'confidences.json'
    if fallback:
        data = json.loads(path.read_text())
        details.write_text(json.dumps({'token_chain_ids': data.pop('chain_ids')}))
        path.write_text(json.dumps(data))
    monkeypatch.setattr(AF3LigandTool, '_sample', lambda *args: {
        'summary_path': str(path), 'confidences_path': str(details) if fallback else None, 'af3_status': 'success'})
    tool = AF3TernaryTool(LigandAF3Config(ligand_chain_id='C'))
    tool.expected_chain_ids = {'A', 'B', 'C'}
    parsed = tool._sample(None, 'fixture', 11, 0)
    assert [(r['chain_a'], r['chain_b']) for r in parsed['chain_pair_metrics']] == [('A','B'),('A','C'),('B','C')]
    assert parsed['chain_pair_metrics'][-1]['chain_pair_iptm'] == .6






def test_binder_rejection_retains_candidate_failure_reasons():
    result = {'modality': 'DE_NOVO_BINDER', 'execution_status': 'COMPLETED',
              'is_mock': False, 'blockers': [], 'summary': {'candidates': [
                  {'candidate_id': 'binder1', 'final_verdict': 'FAIL', 'is_mock': False,
                   'failure_reasons': ['IPTM_BELOW_HARD_LIMIT', 'HOTSPOT_CONTACT_LOST']}]}}
    validation = validate_modality(result)
    assert validation['decision'] == 'REJECT'
    assert validation['reasons'] == [
        'binder1: IPTM_BELOW_HARD_LIMIT', 'binder1: HOTSPOT_CONTACT_LOST']
    result['summary']['candidates'].append(
        {'candidate_id': 'binder2', 'final_verdict': 'PASS', 'is_mock': False})
    validation = validate_modality(result)
    assert validation['decision'] == 'ADVANCE'
    assert validation['reasons'] == []
