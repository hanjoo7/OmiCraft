"""Validation Agent: interpret tool evidence through fixed policy, no metric generation."""


def validate_small_molecule(candidate, policy):
    consensus = candidate["consensus"]
    status = consensus["consensus_status"]
    validation = candidate["rdkit_validation"]
    af3 = candidate["af3_ligand"]
    missing = ["CONSENSUS:" + name for name in consensus.get("missing_metrics", [])]
    failed, passed = [], []
    warnings = list(
        dict.fromkeys(
            validation.get("warning", [])
            + consensus.get("warning", [])
            + consensus.get("reason_codes", [])
        )
    )
    for stage in (
        candidate.get("docking", {}),
        af3,
        *af3.get("samples", []),
        *candidate.get("pose_qc", []),
        *candidate.get("preparation", {}).values(),
    ):
        warnings.extend(stage.get("warning", []))
    warnings = list(dict.fromkeys(warnings))
    if not validation["smiles_valid"]:
        failed.append("INVALID_SMILES")
    else:
        passed.append("RDKIT_SMILES_VALID")
    samples = af3.get("samples", [])
    required = (
        "ligand_atom_plddt_mean",
        "ligand_chain_iptm",
        "protein_ligand_chain_pair_iptm",
        "protein_ligand_pae_min",
        "contact_probability_summary",
        "has_clash",
        "chirality_valid",
    )
    if not samples:
        missing.append("AF3_SAMPLES_MISSING")
    for sample in samples:
        for field in required:
            if sample.get(field) is None:
                missing.append(f"AF3_{sample.get('seed')}_{sample.get('sample_id')}:{field}")
    if any(
        s.get("has_clash") is True
        or s.get("chirality_valid") is False
        or s.get("geometry_status") != "success"
        or s.get("severe_clash_count", 0) > 0
        or s.get("pocket_occupancy") != 1
        for s in samples
    ):
        missing.append("AF3_SAMPLE_GEOMETRY_NOT_SUPPORTED")
    if policy.require_all_af3_seeds and (
        af3.get("status") not in {"success"}
        or set(af3.get("expected_seeds", [])) != set(af3.get("successful_seeds", []))
    ):
        missing.append("AF3_SEED_OR_SAMPLE_OUTPUT_INCOMPLETE")
    if consensus.get("af3_consistency_status") == "inconsistent":
        missing.append("AF3_SEED_DISAGREEMENT")
    if consensus.get("af3_cross_seed_pose_rmsd") is None:
        missing.append("AF3_CROSS_SEED_GEOMETRY_MISSING")
    if consensus.get("docking_af3_pose_rmsd") is None:
        missing.append("DOCKING_AF3_GEOMETRY_MISSING")
    if consensus.get("contact_residue_jaccard") is None:
        missing.append("CONTACT_OVERLAP_MISSING")
    if status == "STRUCTURALLY_UNSUPPORTED":
        failed.append("STRUCTURALLY_UNSUPPORTED")
    if failed:
        call = "FAIL"
    elif status == "DISCORDANT":
        call = policy.discordant_call
    elif missing:
        call = "INDETERMINATE"
    elif status == "CONSENSUS_SUPPORTED":
        call = "PASS_WITH_WARNING" if warnings else "PASS"
        passed.append("STRUCTURAL_CONSENSUS")
    elif status == "POCKET_SUPPORTED_POSE_UNCERTAIN":
        call = "PASS_WITH_WARNING"
        warnings.append("POSE_UNCERTAIN")
    elif status in {"AF3_SUPPORTED_ONLY", "DOCKING_SUPPORTED_ONLY"}:
        call = "REVISE"
    else:
        call = "INDETERMINATE"
    eligible = (
        candidate["route"] == "small_molecule"
        and validation["smiles_valid"]
        and not missing
        and call in policy.eligible_screening_calls
        and status in policy.eligible_consensus_statuses
    )
    return {
        "candidate_id": candidate["candidate_id"],
        "route": candidate["route"],
        "screening_call": call,
        "consensus_status": status,
        "passed_checks": passed,
        "failed_checks": failed,
        "warnings": warnings,
        "missing_evidence": missing,
        "reason_codes": list(dict.fromkeys(failed + missing + consensus.get("reason_codes", []))),
        "eligible_for_admet": eligible,
        "recommended_action": "run_admet" if eligible else "review_evidence",
        "evidence_paths": [sample["model_path"] for sample in samples if sample.get("model_path")],
    }


def validation_node(state):
    """Compatible node facade; Screening invokes the same interpretation per candidate."""
    from .configuration import get_config

    results = [
        validate_small_molecule(c, get_config().small_molecule.validation)
        for c in state.get("screening_result", {}).get("candidates", [])
        if c.get("route") == "small_molecule" and "consensus" in c
    ]
    return {
        "validation_results": results,
        "current_agent": "validation",
        "messages": [f"[Validation] {len(results)} ligand candidates reviewed"],
    }


def validate_modality(result):
    """Separate successful execution from support for a candidate."""
    status, summary = result['execution_status'], result['summary']
    decision = 'NOT_EVALUATED'
    reasons = list(result['blockers'])
    if status in {'COMPLETED', 'PARTIAL'}:
        decision = 'HOLD'
        if result['modality'] == 'ADC' and summary.get('severe_clashes', 0) > 0:
            decision = 'REJECT'
            reasons.append('severe_interface_clashes')
        if result['modality'] in {'SMALL_MOLECULE', 'DE_NOVO_BINDER'}:
            candidates = summary.get('candidates', [])
            if any(c.get('final_verdict') in {'PASS', 'PASS_WITH_WARNING', 'GO', 'GO_WITH_WARNING'}
                   and not c.get('is_mock', True) for c in candidates):
                decision = 'ADVANCE'
            elif candidates and all(c.get('final_verdict') in {'FAIL', 'NO_GO', 'REJECT'} for c in candidates):
                decision = 'REJECT'
            if decision != 'ADVANCE':
                for candidate in candidates:
                    candidate_id = candidate.get('candidate_id') or 'candidate'
                    for reason in candidate.get('failure_reasons', []):
                        message = f'{candidate_id}: {reason}'
                        if message not in reasons:
                            reasons.append(message)
    if result['is_mock'] and decision == 'ADVANCE':
        decision = 'HOLD'
        reasons.append('mock_not_scientific_evidence')
    evidence = {}
    if result['modality'] == 'DEGRADER' and summary.get('input_mode') == 'published_full_molecule':
        qc = summary.get('qc', {})
        evidence = dict(provenance_valid=summary.get('provenance', {}).get('curation_status') == 'VERIFIED',
            full_molecule_valid=summary.get('full_molecule_validation') == 'COMPLETED',
            component_mapping_available=summary.get('component_mapping_status') == 'VERIFIED',
            assembly_performed=summary.get('assembly_performed', False),
            conformer_valid=qc.get('conformer_generation_status') == 'COMPLETED' and qc.get('conformer_identity_preserved') is True,
            ternary_readiness=summary.get('ternary_status'),
            ternary_prediction_performed=summary.get('ternary_prediction_performed', False),
            degradation_efficacy='NOT_EVALUATED')
        evidence['reference_evaluation_supported'] = all(evidence[k] for k in ('provenance_valid', 'full_molecule_valid', 'conformer_valid')) and not result['is_mock']
    computational = {}
    if result['modality'] in {'ADC', 'DEGRADER'}:
        from .modality_validation import structural_assessment
        computational = structural_assessment(result)
        if status in {'COMPLETED', 'PARTIAL'}:
            decision = computational['decision']
        reasons = list(dict.fromkeys(reasons + computational['reasons']))
    return {'decision': decision, 'reasons': reasons, 'scope': 'computational_checks_only',
            'computational_validation': computational,
            'missing_evidence': list(result.get('missing_steps', [])),
            'experimental_validation': 'NOT_EVALUATED', **evidence}


def binder_developability(sequence, structure_path=None, binder_chain='B'):
    """Uncalibrated sequence/geometry screening; not an experimental prediction."""
    import re
    from collections import Counter

    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
        return {'status': 'HOLD', 'warnings': ['nonstandard_or_missing_sequence'],
                'experimental_developability': 'NOT_EVALUATED'}
    hydrophobic = set('AVILMFWY')
    analysis = ProteinAnalysis(sequence)
    stretches = [{'start': m.start()+1, 'end': m.end(), 'sequence': m.group()}
                 for m in re.finditer(r'[AVILMFWY]{5,}', sequence)]
    motifs = [{'start': m.start()+1, 'motif': m.group(), 'kind': label}
              for pattern, label in [(r'N[^P][ST]', 'N_linked_glycosylation'),
                                     (r'N[GS]', 'deamidation_liability'), (r'D[GP]', 'isomerization_liability')]
              for m in re.finditer(pattern, sequence)]
    fraction = sum(a in hydrophobic for a in sequence) / len(sequence)
    warnings = []
    if stretches:
        warnings.append('long_hydrophobic_stretch')
    if fraction > .45:
        warnings.append('high_hydrophobic_fraction')
    if sequence.count('C') % 2:
        warnings.append('odd_cysteine_count')
    warnings.extend(sorted({m['kind'] for m in motifs}))
    output = {'scope': 'rule_based_demo_screening_not_experimental', 'length': len(sequence),
              'composition': dict(sorted(Counter(sequence).items())), 'hydrophobic_fraction': fraction,
              'hydrophobic_residues': ''.join(sorted(hydrophobic)), 'net_charge_ph7': analysis.charge_at_pH(7),
              'estimated_pI': analysis.isoelectric_point(), 'charge_method': 'Bio.SeqUtils.ProtParam',
              'cysteine_count': sequence.count('C'), 'hydrophobic_stretches': stretches,
              'liability_motifs': motifs, 'aggregation_warning': bool(stretches or fraction > .45),
              'low_confidence_residue_fraction': None, 'exposed_hydrophobic_patch': 'not_available',
              'experimental_binding': 'NOT_EVALUATED', 'experimental_developability': 'NOT_EVALUATED',
              'heuristic_policy': {'hydrophobic_stretch_min_length': 5, 'hydrophobic_fraction_warning': .45,
                                   'low_confidence_plddt_below': 50, 'surface_sasa_above_angstrom2': 20,
                                   'patch_ca_distance_angstrom': 8, 'patch_warning_residues': 6,
                                   'threshold_status': 'uncalibrated_demo_heuristics'}}
    if structure_path:
        try:
            import numpy as np
            from Bio.PDB import MMCIFParser, PDBParser
            from Bio.PDB.SASA import ShrakeRupley
            from Bio.SeqUtils import seq1

            parser = MMCIFParser(QUIET=True) if str(structure_path).endswith('.cif') else PDBParser(QUIET=True)
            model = parser.get_structure('binder', str(structure_path))[0]
            residues = [r for r in model[binder_chain] if r.id[0] == ' ' and 'CA' in r]
            if ''.join(seq1(r.resname) for r in residues) != sequence:
                raise ValueError('binder_structure_sequence_mismatch')
            scores = [float(r['CA'].bfactor) for r in residues]
            output['low_confidence_residue_fraction'] = sum(v < 50 for v in scores) / len(scores)
            ShrakeRupley(n_points=100).compute(model, level='R')
            exposed = [r for r in residues if seq1(r.resname) in hydrophobic and r.sasa > 20]
            remaining = set(range(len(exposed)))
            sizes = []
            while remaining:
                frontier = [remaining.pop()]
                size = 0
                while frontier:
                    i = frontier.pop()
                    size += 1
                    near = {j for j in remaining if np.linalg.norm(exposed[i]['CA'].coord-exposed[j]['CA'].coord) <= 8}
                    remaining -= near
                    frontier.extend(near)
                sizes.append(size)
            output['exposed_hydrophobic_patch'] = {'method': 'ShrakeRupley SASA + CA-distance connected components',
                'largest_patch_residues': max(sizes, default=0), 'exposed_hydrophobic_residues': len(exposed),
                'scope': 'geometry_proxy_in_predicted_complex_not_aggregation_prediction'}
            if max(sizes, default=0) >= 6:
                warnings.append('exposed_hydrophobic_patch_proxy')
        except (ImportError, OSError, ValueError, KeyError) as exc:
            warnings.append('structure_developability_not_available:' + str(exc))
    else:
        warnings.append('structure_developability_not_available')
    output.update(warnings=warnings, status='PASS_WITH_WARNINGS' if warnings else 'PASS')
    return output


def validate_demo_modality(result):
    base = validate_modality(result)
    missing = list(result.get('missing_steps', []))
    if result['modality'] == 'DE_NOVO_BINDER':
        developability = result.get('demo_metrics', {}).get('developability', {})
        if developability.get('status') == 'HOLD':
            base['decision'] = 'HOLD'
        base['developability_screening'] = developability.get('status', 'NOT_EVALUATED')
        base['structural_screening'] = result.get('demo_metrics', {}).get('structural_screening', 'NOT_EVALUATED')
    base.update(missing_evidence=missing, experimental_binding='NOT_EVALUATED',
                experimental_developability='NOT_EVALUATED')
    return base
