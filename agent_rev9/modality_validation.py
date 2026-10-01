"""Explicit structural triage policies. Passing supports prioritization, not efficacy."""
import math


def number(value):
    return type(value) in (int, float) and math.isfinite(value)


def structural_assessment(result):
    mode, summary = result['modality'], result.get('summary', {})
    assessment = dict(decision='HOLD', passed_checks=[], reasons=[],
                      scope=result.get('execution_scope', 'computational_screening'),
                      policy_version='structural_triage_v1',
                      threshold_status='heuristic_not_experimentally_calibrated',
                      experimental_validation='NOT_EVALUATED')
    reasons, passed = assessment['reasons'], assessment['passed_checks']
    if result.get('is_mock'):
        reasons.append('mock_not_scientific_evidence')
    if result.get('execution_status') != 'COMPLETED':
        reasons.append('computational_steps_incomplete')
    reasons.extend(result.get('blockers', []))
    if mode == 'ADC':
        checks = {
            'cdr_mapping': summary.get('cdr_mapping_status') == 'success',
            'cdr_contacts': number(summary.get('cdr_heavy_atom_contacts')) and summary['cdr_heavy_atom_contacts'] > 0,
            'no_severe_clashes': summary.get('severe_clashes') == 0,
            'conjugation_analysis': summary.get('conjugation', {}).get('status') == 'COMPLETED',
            'accessible_noninterface_site': any(
                row.get('inside_cdr') is False and row.get('at_antigen_interface') is False
                and row.get('possible_disulfide') is False and not row.get('exclusion_reasons')
                and number(row.get('sasa_angstrom2')) and row['sasa_angstrom2'] > 0
                for row in summary.get('conjugation', {}).get('candidates', [])),
        }
        failed_geometry = number(summary.get('severe_clashes')) and summary['severe_clashes'] > 0
        assessment['scope'] = 'antibody_structure_and_conjugation_feasibility'
        assessment['thresholds'] = {'cdr_contacts_min': 1, 'severe_clashes_max': 0, 'site_sasa_min_exclusive': 0}
    elif mode == 'DEGRADER':
        qc = summary.get('qc', {})
        samples = summary.get('ternary_prediction', {}).get('samples', [])
        assessment['scope'] = 'reference_chemistry_and_ternary_structural_triage'
        assessment['thresholds'] = {'iptm_min': .6, 'ligand_pair_iptm_min': .5, 'ligand_pair_pae_max': 10,
                                    'independent_seeds_min': 2, 'contact_distance_angstrom': 4.5}
        checks = {
            'verified_molecule': summary.get('full_molecule_validation') == 'COMPLETED'
                and summary.get('provenance', {}).get('curation_status') == 'VERIFIED',
            'conformer_identity': qc.get('conformer_generation_status') == 'COMPLETED'
                and qc.get('conformer_identity_preserved') is True,
            'ternary_prediction': summary.get('ternary_prediction', {}).get('status') == 'success' and bool(samples),
            'independent_seeds': len({s.get('seed') for s in samples if s.get('seed') is not None}) >= 2,
        }
        supported = []
        for sample in samples:
            geometry = sample.get('ternary_geometry', {})
            target, ligand = sample.get('protein_chain_id'), sample.get('ligand_chain_id')
            contacts = geometry.get('ligand_contacts_by_chain', {})
            pairs = sample.get('chain_pair_metrics', [])
            def interface(other):
                return any({p.get('chain_a'), p.get('chain_b')} == {other, ligand}
                           and number(p.get('chain_pair_iptm')) and .5 <= p['chain_pair_iptm'] <= 1
                           and number(p.get('chain_pair_pae_min')) and 0 <= p['chain_pair_pae_min'] <= 10 for p in pairs)
            partners = [chain for chain in contacts if chain not in {target, ligand}]
            supported.append(sample.get('af3_status') == 'success' and sample.get('has_clash') is False
                and sample.get('chirality_valid') is True and geometry.get('status') == 'COMPLETED'
                and geometry.get('severe_clash_count') == 0
                and number(sample.get('iptm')) and .6 <= sample['iptm'] <= 1
                and contacts.get(target, 0) > 0 and interface(target)
                and any(contacts.get(chain, 0) > 0 and interface(chain) for chain in partners))
        checks['ternary_interfaces_and_geometry'] = bool(supported) and all(supported)
        failed_geometry = bool(samples) and all(s.get('chirality_valid') is False or s.get('has_clash') is True
            or s.get('ternary_geometry', {}).get('severe_clash_count', 0) > 0 for s in samples)
    else:
        raise ValueError('Unsupported structural triage modality')
    for name, valid in checks.items():
        (passed if valid else reasons).append(name if valid else name + '_not_supported')
    assessment['reasons'] = list(dict.fromkeys(reasons))
    if failed_geometry and not result.get('is_mock'):
        assessment['decision'] = 'REJECT'
    elif not reasons:
        assessment['decision'] = 'ADVANCE'
    return assessment
