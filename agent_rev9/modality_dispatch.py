"""Small executors sharing a result contract; no implicit model fallback."""

from pathlib import Path

from .routes import MODALITY_ROUTES, validate_route


def modality_name(value):
    route = validate_route(value)
    return next(name for name, target in MODALITY_ROUTES.items() if route == target)


def result(modality, status, stage, *, summary=None, artifacts=None, blockers=None, is_mock=False):
    return {'modality': modality_name(modality), 'execution_status': status,
            'validation_decision': 'NOT_EVALUATED', 'stage': stage, 'summary': summary or {},
            'artifacts': artifacts or [], 'blockers': blockers or [], 'is_mock': is_mock}


def binder_configuration(config):
    from .configuration import OmiCraftConfig, get_config
    values = get_config().model_dump()
    for name in ('rfd3', 'protein_mpnn', 'af3_binder', 'binder_analysis'):
        if name in config:
            values[name] = {**values[name], **config[name]}
    return OmiCraftConfig(**values)


def readiness(modality, data, config):
    from .docking_backend import tool_readiness

    name = modality_name(modality)
    blockers, tools = [], {}
    if name == 'ADC':
        if not data.get('structure_path') or not Path(data['structure_path']).is_file():
            blockers.append('antibody_complex_missing')
        if not data.get('antigen_chain'):
            blockers.append('antigen_chain_missing')
        from .antibody_numbering import AnarciiBackend
        backend = config.get('numbering_backend', 'anarcii')
        if backend == 'anarcii':
            tools['numbering'] = AnarciiBackend(config.get('numbering_python_path')).availability()
            if tools['numbering']['status'] != 'ready':
                blockers.append('numbering_backend_unavailable')
        elif backend == 'anarci':
            import importlib.util
            if importlib.util.find_spec('anarci') is None:
                blockers.append('numbering_backend_unavailable')
        else:
            blockers.append('unknown_numbering_backend')
    elif name == 'DEGRADER':
        from .degrader_execution import input_blockers
        blockers.extend(input_blockers(data))
    elif name == 'SMALL_MOLECULE':
        from .af3_ligand_tool import AF3LigandTool
        from .small_molecule_config import SmallMoleculeConfig
        model = SmallMoleculeConfig(**config.get('small_molecule', {}))
        if model.docking.backend == 'mock':
            if not config.get('is_mock'):
                blockers.append('mock_requires_explicit_flag')
        else:
            tools['docking'] = tool_readiness(model.docking.backend,
                configured=model.docking.executable, environment=model.docking.environment)
            if tools['docking']['blocker']:
                blockers.append(tools['docking']['blocker'])
        context = data.get('context', {})
        for key in ('target_sequence', 'prepared_receptor_path'):
            if not context.get(key):
                blockers.append('missing:' + key)
        if context.get('prepared_receptor_path') and not Path(context['prepared_receptor_path']).is_file():
            blockers.append('prepared_receptor_missing')
        if not data.get('candidates'):
            blockers.append('candidates_missing')
        if model.run_af3 and model.docking.backend != 'mock':
            tools['af3'] = AF3LigandTool(model.af3).dependency_status()
            if tools['af3'] != 'available':
                blockers.append('af3:' + str(tools['af3']))
    else:
        cfg = binder_configuration(config)
        redesign = bool((data.get('protein_design_input') or {}).get('backbone_paths'))
        tools['rfd3'] = {'path': None, 'required': not redesign}
        import shutil
        tools['rfd3']['path'] = shutil.which(cfg.rfd3.executable)
        if not redesign and not tools['rfd3']['path']:
            blockers.append('rfd3_executable_missing')
        for group, fields in ((cfg.protein_mpnn, ('python_path', 'script_path', 'model_weights_dir')),
                              (cfg.af3_binder, ('python_path', 'script_path', 'model_dir', 'database_dir'))):
            for key in fields:
                value = getattr(group, key)
                if not value or not Path(value).exists():
                    blockers.append(type(group).__name__ + ':' + key + '_missing')
        if data.get('design_results_path'):
            if not Path(data['design_results_path']).is_file():
                blockers.append('binder_design_input_missing')
        elif not data.get('protein_design_input'):
            blockers.append('binder_design_input_missing')
        design = data.get('protein_design_input') or {}
        if design:
            from .structure_io import protein_residues
            try:
                residues = protein_residues(design.get('target_structure', ''))
                allowed = {chain + number for chain, number, _ in residues
                           if chain == design.get('target_chain', 'A')}
                if not allowed:
                    blockers.append('binder_target_chain_missing')
                for field in ('hotspot_residues', 'motif_residues', 'binding_site'):
                    for token in design.get(field, []):
                        if token not in allowed:
                            blockers.append('Selection is not a target residue: ' + token)
            except (OSError, ValueError) as exc:
                blockers.append('binder_target_input: ' + str(exc))
        if redesign:
            from .binder_analysis import protein_chains
            for path in design['backbone_paths']:
                try:
                    chains = protein_chains(path)
                    if design.get('binder_chain') == design.get('target_chain'):
                        raise ValueError('distinct_target_and_binder_chains_required')
                    if chains.get(design.get('target_chain'), {}).get('sequence') != design.get('target_sequence'):
                        raise ValueError('target_sequence_mismatch')
                    if len(chains.get(design.get('binder_chain'), {}).get('sequence', '')) != design.get('design_length'):
                        raise ValueError('binder_sequence_length_mismatch')
                except (OSError, ValueError, KeyError) as exc:
                    blockers.append('redesign_backbone: ' + str(exc))
        checkpoint = Path(cfg.protein_mpnn.model_weights_dir) / (cfg.protein_mpnn.model_name + '.pt')
        if not checkpoint.is_file():
            blockers.append('protein_mpnn_weights_missing')
        if not redesign and cfg.rfd3.checkpoint_path and not Path(cfg.rfd3.checkpoint_path).is_file():
            blockers.append('rfd3_checkpoint_missing')
        if cfg.af3_binder.model_dir and not any(Path(cfg.af3_binder.model_dir).glob('af3.bin*')):
            blockers.append('af3_weights_missing')
    return {'tool_status': 'ready' if not blockers else 'not_available', 'tools': tools,
            'blockers': blockers}


def execute_adc(input_data, config):
    from .antibody_mapping import analyze, conjugation_feasibility
    from .small_molecule_io import write_csv

    import shutil
    source = Path(input_data.get('structure_path', ''))
    if source.is_file():
        directory = Path(config['output_dir'])
        directory.mkdir(parents=True, exist_ok=True)
        copied = directory / ('antibody_complex' + source.suffix)
        if source.resolve() != copied.resolve():
            shutil.copy2(source, copied)
        input_data = {**input_data, 'structure_path': str(copied.resolve())}
    analysis = analyze(input_data, config)
    artifacts = [{'path': analysis['structure_path'], 'label': 'Input antibody complex'}]
    if analysis['cdr_mapping_status'] == 'success':
        for key in ('mapping', 'contacts'):
            path = Path(config['output_dir']) / (key + '.csv')
            write_csv(path, analysis[key])
            artifacts.append({'path': str(path.resolve()), 'label': key})
    mapped = analysis['cdr_mapping_status'] == 'success'
    checklist = [dict(stage='cdr_mapping', execution_status='COMPLETED' if mapped else 'BLOCKED')]
    blockers = list(analysis['blockers'])
    if mapped:
        try:
            analysis['conjugation'] = conjugation_feasibility(analysis)
            path = Path(config['output_dir']) / 'conjugation_candidates.csv'
            write_csv(path, analysis['conjugation']['candidates'])
            artifacts.append({'path': str(path.resolve()), 'label': 'Conjugation candidates'})
            checklist.append(dict(stage='conjugation_feasibility', execution_status='COMPLETED'))
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            blockers.append('conjugation_analysis_failed: ' + str(exc))
            checklist.append(dict(stage='conjugation_feasibility', execution_status='FAILED'))
    missing = ['linker_payload_selection', 'DAR_assessment', 'internalization', 'experimental_efficacy']
    analysis.update(execution_scope='antibody_structure_and_conjugation_feasibility',
                    execution_checklist=checklist + [dict(stage=s, execution_status='NOT_EVALUATED') for s in missing])
    complete = mapped and not blockers and all(step['execution_status'] == 'COMPLETED' for step in checklist)
    output = result('ADC', 'COMPLETED' if complete else 'PARTIAL' if mapped else 'BLOCKED', 'conjugation_feasibility',
                    summary=analysis, artifacts=artifacts, blockers=blockers)
    output.update(execution_scope=analysis['execution_scope'], missing_steps=missing)
    return output


def execute_degrader(input_data, config):
    from .degrader_execution import execute as run
    return run(input_data, config)


def execute_small_molecule(input_data, config):
    from .screening_tools import ScreeningOptions
    from .small_molecule_config import SmallMoleculeConfig
    from .small_molecule_workflow import run_small_molecule_workflow

    directory = Path(config['output_dir']) / 'small_molecule'
    model = SmallMoleculeConfig(**config.get('small_molecule', {}))
    options = ScreeningOptions.from_state({'screening_options': config.get('screening_options', {})})
    output = run_small_molecule_workflow(input_data['candidates'], input_data['context'], model,
                                         directory, options=options, dry_run=False,
                                         qualification=input_data.get('qualification'),
                                         progress=config.get('_progress_callback'))
    candidates = output.get('candidates', [])
    checklist = []
    for candidate in candidates:
        checks = {
            'receptor_preparation': candidate.get('preparation', {}).get('receptor', {}).get('status') == 'success',
            'ligand_preparation': candidate.get('preparation', {}).get('ligand', {}).get('status') == 'success',
            'docking': not model.run_docking or candidate.get('docking', {}).get('status') == 'success',
            'pose_qc': not model.run_pose_qc or bool(candidate.get('pose_qc')),
            'af3': not model.run_af3 or candidate.get('af3_ligand', {}).get('status') == 'success',
            'structural_comparison': not model.run_consensus or bool(candidate.get('consensus', {}).get('comparisons')),
            'admet': not model.run_admet or candidate.get('admet', {}).get('status') == 'success',
            'validation_and_critic': bool(candidate.get('validation')) and bool(candidate.get('critic')),
            'report': Path(output['candidate_summary_path']).is_file(),
        }
        enabled = dict(docking=model.run_docking, pose_qc=model.run_pose_qc, af3=model.run_af3,
                       structural_comparison=model.run_consensus, admet=model.run_admet)
        checklist.extend(dict(candidate_id=candidate.get('candidate_id'), stage=key,
                              execution_status='NOT_REQUESTED' if not enabled.get(key, True)
                              else 'COMPLETED' if done else 'INCOMPLETE') for key, done in checks.items())
    output['execution_checklist'] = checklist
    complete = bool(checklist) and all(row['execution_status'] in {'COMPLETED', 'NOT_REQUESTED'} for row in checklist)
    blockers = [row['stage'] + '_incomplete:' + str(row['candidate_id'])
                for row in checklist if row['execution_status'] == 'INCOMPLETE']
    artifacts = [{'path': str(p.resolve()), 'label': p.name}
                 for p in directory.rglob('*') if p.is_file()]
    status = 'COMPLETED' if candidates and complete else 'PARTIAL'
    return result('SMALL_MOLECULE', status, 'small_molecule', summary=output,
                  artifacts=artifacts, blockers=blockers, is_mock=model.docking.backend == 'mock')


def execute_binder(input_data, config):
    from .configuration import configuration_scope
    from .screening_agent import screening_node

    state = {**input_data, 'route': 'protein_binder_rfd3', 'use_rfd3': True, 'dry_run': False}
    state['_progress_callback'] = config.get('_progress_callback', lambda *args: None)
    state['protein_design_input'] = {**state.get('protein_design_input', {}),
                                     'output_dir': str(Path(config['output_dir']) / 'binder')}
    state["screening_options"] = {**state.get("screening_options", {}), **config.get("screening_options", {})}
    with configuration_scope(binder_configuration(config)):
        output = screening_node(state)
    screening = output.get('screening_result', {})
    status = screening.get('status', 'FAILED')
    status = status if status in {'COMPLETED', 'PARTIAL', 'BLOCKED', 'FAILED', 'NOT_RUN'} else {
        'FAILED_EXECUTION': 'FAILED', 'NOT_CONFIGURED': 'BLOCKED', 'NOT_EXECUTED': 'NOT_RUN'
    }.get(status, 'PARTIAL')
    artifacts = [{'path': str(p.resolve()), 'label': p.name}
                 for p in (Path(config['output_dir']) / 'binder').rglob('*') if p.is_file()]
    blockers = list(output.get('errors') or [])
    if status != 'COMPLETED':
        def collect(value):
            if isinstance(value, dict):
                for key, item in value.items():
                    if key in {'error_message', 'failure_reasons', 'errors', 'blockers'}:
                        for reason in item if isinstance(item, list) else [item]:
                            if isinstance(reason, str) and reason:
                                blockers.append(reason)
                    elif isinstance(item, (dict, list)):
                        collect(item)
            elif isinstance(value, list):
                for item in value:
                    collect(item)
        collect(output)
        if not blockers:
            blockers.append('binder_execution_incomplete: ' + status)
    return result('DE_NOVO_BINDER', status, 'binder', summary=screening, artifacts=artifacts,
                  blockers=list(dict.fromkeys(blockers)),
                  is_mock=any(c.get('is_mock') for c in screening.get('candidates', [])))


def execute(input_data, config):
    from .critic_agent import review_modality
    from .run_records import RunRecord
    from .validation_agent import validate_modality

    modality = modality_name(input_data['modality'])
    record = RunRecord(config['output_dir'], modality)
    record.data['selection_event'] = input_data.get('selection_event')
    data = input_data.get('input', input_data)
    try:
        if modality == 'DEGRADER' and not any(k in data for k in ('input_mode', 'warhead_smiles', 'full_degrader_smiles')):
            data = select_degrader_reference(data)
        if modality == 'DEGRADER':
            for item in data.get('selection_trace', []):
                record.event(item['stage'], 'COMPLETED', execution_mode=item['mode'], message=item['message'])
        ready = readiness(modality, data, config)
        record.data['readiness'] = ready
        record.event('readiness', 'COMPLETED' if ready['tool_status'] == 'ready' else 'BLOCKED')
        if ready['blockers']:
            output = result(modality, 'BLOCKED', 'readiness', summary=ready, blockers=ready['blockers'])
            if modality == 'ADC':
                output['summary'].update(numbering_status='backend_unavailable'
                    if 'numbering_backend_unavailable' in ready['blockers'] else 'not_run', cdr_mapping_status='not_run')
        elif config.get('dry_run', True):
            output = result(modality, 'NOT_RUN', 'execution', summary={'readiness': ready})
        else:
            record.event('execution', 'RUNNING', event='started')
            executors = {'SMALL_MOLECULE': execute_small_molecule, 'DE_NOVO_BINDER': execute_binder,
                         'ADC': execute_adc, 'DEGRADER': execute_degrader}
            output = executors[modality](data, config)
            record.event('execution', output['execution_status'], event='finished')
        output['is_mock'] = bool(output['is_mock'] or config.get('is_mock') or data.get('is_mock'))
        validation = validate_modality(output)
        output.update(validation=validation, validation_decision=validation['decision'],
                      computational_validation=validation.get('computational_validation', {}),
                      experimental_validation='NOT_EVALUATED')
        record.event('validation', 'COMPLETED', validation_decision=validation['decision'])
        output['critic'] = review_modality(output)
        record.event('critic', 'COMPLETED', validation_decision=output['validation_decision'], critic_decision=output['critic']['decision'])
        return record.finish(output)
    except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
        output = result(modality, 'FAILED', 'execution', blockers=[str(exc)])
        output['validation'] = {'decision': 'NOT_EVALUATED', 'reasons': [str(exc)]}
        output['critic'] = review_modality(output)
        return record.finish(output)


DEMO_TOOL_CHAINS = {
    'SMALL_MOLECULE': ['RDKit', 'GNINA', 'pose QC', 'AF3', 'structural comparison', 'ADMET'],
    'DE_NOVO_BINDER': ['RFD3', 'ProteinMPNN', 'AF3', 'interface QC', 'developability screening'],
    'ADC': ['ANARCII', 'sequence/structure mapping', 'CDR/contact analysis', 'conjugation feasibility'],
    'DEGRADER': ['RDKit component QC', 'attachment validation', 'assembly', 'geometry/ternary readiness'],
}


def replay_demo_modality(modality, source_run, record, config):
    import csv
    import json

    from .antibody_mapping import conjugation_feasibility
    from .degrader_execution import component_qc
    from .run_records import register_demo_artifact
    from .small_molecule_io import sha256
    from .validation_agent import binder_developability

    folder = {'SMALL_MOLECULE': 'small_molecule', 'DE_NOVO_BINDER': 'binder',
              'ADC': 'adc', 'DEGRADER': 'degrader'}[modality]
    source_run = Path(source_run).resolve()
    hashes, provenance = {}, []
    events = source_run / 'events.jsonl'
    if events.is_file():
        for line in events.read_text().splitlines():
            hashes.update(json.loads(line).get('output_sha256', {}))

    def artifact(path, label=None):
        path = Path(path).resolve()
        if not path.is_file() or not path.is_relative_to(source_run.parent):
            raise ValueError('replay_artifact_missing_or_outside_runs:' + str(path))
        digest = sha256(path)
        if str(path) in hashes and hashes[str(path)] != digest:
            raise ValueError('source_artifact_checksum_mismatch:' + str(path))
        provenance.append({'source_run_id': source_run.name, 'source_artifact_path': str(path),
            'source_artifact_sha256': digest, 'checksum_verification': 'matched_source_manifest'
                if str(path) in hashes else 'current_snapshot_no_prior_checksum',
            'execution_mode': 'REAL_ARTIFACT_REPLAY', 'is_mock': False})
        register_demo_artifact(record, path, modality, 'artifact_replay', title=label,
            artifact_type='structure' if path.suffix in {'.cif', '.sdf', '.gz'} else 'table')
        return path

    original = json.loads(artifact(source_run / folder / 'run_summary.json', modality + ' original summary').read_text())
    if original.get('is_mock') is not False or original.get('modality') != modality:
        raise ValueError('real_source_modality_required')
    output = {**original, 'source_run_id': source_run.name, 'demo_metrics': {},
              'missing_steps': list(original.get('missing_steps', [])), 'blockers': list(original.get('blockers', [])),
              'source_provenance': provenance, 'real_tool_execution_status': 'REAL_ARTIFACT_REPLAY',
              'execution_mode': 'REAL_ARTIFACT_REPLAY', 'lightweight_calculations': [], 'is_mock': False}
    summary, metrics = original['summary'], output['demo_metrics']
    if modality in {'SMALL_MOLECULE', 'DE_NOVO_BINDER'}:
        candidate = summary['candidates'][0]
        if candidate.get('is_mock') is not False:
            raise ValueError('real_candidate_required')
    if modality == 'SMALL_MOLECULE':
        poses = candidate['docking']['poses']
        for pose in poses:
            artifact(pose['pose_path'], 'GNINA pose '+str(pose['pose_rank']))
        sample = candidate['af3_ligand']['samples'][0]
        artifact(sample['model_path'], 'Small molecule AF3 structure')
        artifact(sample['summary_path'], 'Small molecule AF3 confidence')
        path = artifact(source_run / 'small_molecule/admet_standalone/output/admet_predictions_long.csv', 'All ADMET endpoints')
        with path.open() as handle:
            endpoints = list(csv.DictReader(handle))
        if any(r.get('is_mock', '').lower() != 'false' for r in endpoints):
            raise ValueError('real_admet_required')
        core = []
        for label, name in {'solubility': 'Solubility_AqSolDB', 'BBB': 'BBB_Martins', 'hERG': 'hERG', 'AMES': 'AMES', 'DILI': 'DILI'}.items():
            row = next((r for r in endpoints if r['endpoint_name'] == name), {})
            core.append({'endpoint': label, 'backend_endpoint': name, 'value': row.get('endpoint_value') or 'not_available',
                         'unit': row.get('endpoint_unit') or 'not_available', 'status': row.get('prediction_status') or 'not_available'})
        metrics.update(poses=[{k: p.get(k) for k in ('pose_rank', 'docking_score', 'cnn_score', 'cnn_affinity')} for p in poses],
            pose_qc=candidate['pose_qc'], comparisons=candidate['consensus']['comparisons'],
            consensus=candidate['consensus']['consensus_status'], af3_plddt=sample.get('ligand_atom_plddt_mean'),
            admet_core=core, available_endpoint_count=sum(r['prediction_status'] == 'success' for r in endpoints),
            missing_endpoint_count=sum(r['prediction_status'] != 'success' for r in endpoints),
            admet_scope='standalone backend results; structural handoff remains blocked',
            preparation=candidate['preparation'], structural_decision=candidate['final_verdict'])
    elif modality == 'DE_NOVO_BINDER':
        artifact(candidate['rfd3_structure_path'], 'Actual RFD3 backbone')
        artifact(candidate['protein_mpnn']['sequence_path'], 'Actual ProteinMPNN sequence')
        af3 = candidate['af3_metrics']
        structure = artifact(af3['model_path'], 'Target–binder AF3 structure')
        artifact(af3['summary_path'], 'Binder AF3 confidence')
        development = binder_developability(candidate['sequence'], structure)
        metrics.update(sequence=candidate['sequence'], structure_path=str(structure),
            protein_mpnn_score=candidate['protein_mpnn']['score'], protein_mpnn_score_meaning='sequence design score; not binding affinity',
            structural_screening=candidate['primary_structural_screening']['verdict'], iptm=af3.get('iptm'),
            ptm=af3.get('ptm'), interface_pae_min=af3.get('chain_pair_pae_min'), confidence=af3.get('plddt_summary'),
            interface=candidate['structural_metrics'], developability=development, cross_validation_status='not_used',
            protenix_enabled=False, protenix_required=False, experimental_binding='NOT_EVALUATED')
        output['missing_steps'] = [v for v in output['missing_steps'] if v != 'developability_backend_not_implemented'] + ['experimental_developability']
        output['lightweight_calculations'].append({'tool': 'Bio.ProtParam + ShrakeRupley', 'execution_mode': 'FRESH_ANALYSIS', 'scope': 'rule_based_developability'})
    elif modality == 'ADC':
        artifact(summary['structure_path'], 'Observed antibody–ERBB2 structure')
        for name in ['antibody_numbering.tsv', 'sequence_structure_mapping.tsv', 'cdr_contacts.tsv']:
            artifact(source_run / folder / name, name)
        feasibility = conjugation_feasibility(summary)
        template = None
        if config.linker_payload_template:
            template = json.loads(Path(config.linker_payload_template).read_text())
            if not isinstance(template, dict) or not template.get('source') or template.get('is_mock') is not False:
                raise ValueError('curated_adc_template_requires_source_and_is_mock_false')
        metrics.update(mapping=summary['mapping'], mapping_coverage={k: v['coverage'] for k, v in summary['mapping_qc'].items()},
            cdr_contacts=[{'region': region, 'heavy_atom_contacts': sum(r['heavy_atom_contact_count'] for r in summary['contacts'] if r['region'] == region)} for region in ('CDR1', 'CDR2', 'CDR3')],
            contact_count=summary['cdr_heavy_atom_contacts'], clash_count=summary['severe_clashes'],
            conjugation=feasibility, linker_payload=template or 'not_configured', DAR='NOT_RUN',
            internalization='NOT_EVALUATED', cytotoxicity='NOT_EVALUATED', feasibility_execution_status='COMPLETED')
        output['lightweight_calculations'].append({'tool': 'Bio.PDB + ShrakeRupley', 'execution_mode': 'FRESH_ANALYSIS', 'scope': 'conjugation_feasibility'})
    else:
        data = json.loads(artifact(source_run / folder / 'input.json', 'Explicit degrader input').read_text())
        if config.degrader_component_config:
            data = json.loads(Path(config.degrader_component_config).read_text())
            if data.get('is_mock') is not False:
                raise ValueError('explicit_real_degrader_components_required')
        elif config.degrader_reference_config:
            data = select_degrader_reference(data, config.degrader_reference_config)
        if data.get('input_mode') == 'published_full_molecule':
            directory = record.directory / 'degrader'
            evaluated = execute_degrader(data, {'output_dir': str(directory)})
            output.update(evaluated, previous_attempt={'candidate_id': 'TAK-285', 'status': 'BLOCKED_INPUT',
                'reason': 'E3 ligand, linker and attachment evidence absent', 'source_run_id': source_run.name},
                demo_metrics=evaluated['summary'], execution_mode='FRESH_ANALYSIS',
                real_tool_execution_status='FRESH_ANALYSIS', missing_steps=['ternary_prediction', 'experimental_degradation'])
            for item in evaluated['artifacts']:
                register_demo_artifact(record, item['path'], modality, 'reference_evaluation',
                    artifact_type='plot' if Path(item['path']).suffix == '.png' else 'structure' if Path(item['path']).suffix == '.sdf' else 'table',
                    source_table=directory / 'rdkit_validation.json', execution_mode='FRESH_ANALYSIS')
            return output
        qc = component_qc('supplied_warhead', 'warhead', data.get('warhead_smiles'))
        metrics.update(warhead_smiles=data.get('warhead_smiles'), warhead_qc=qc,
            input_schema=['warhead_smiles', 'e3_ligand_smiles', 'linker_smiles', 'warhead_map', 'e3_map', 'linker_maps'],
            target_context=data.get('target_id', 'ERBB2'), ternary_status='NOT_RUN',
            next_action='Supply an E3 ligand, linker and explicit attachment maps; no exit vector will be inferred.')
        ready = readiness(modality, data, {})
        metrics['missing_inputs'] = ready['blockers']
        if ready['blockers']:
            output.update(execution_status='BLOCKED', blockers=ready['blockers'])
            metrics.update(assembly='NOT_RUN', geometry='NOT_RUN')
        else:
            context = data.get('target_context', data.get('target_id'))
            if context not in {'ERBB2', 'BRD4–VHL reference demo'}:
                raise ValueError('degrader_target_context_required')
            directory = record.directory / 'degrader'
            directory.mkdir(exist_ok=True)
            assembled = execute_degrader(data, {'output_dir': str(directory)})
            output.update(execution_status=assembled['execution_status'], blockers=assembled['blockers'])
            metrics.update(assembly=assembled['summary'], target_context=context,
                           attachment_maps={k: data[k] for k in ('warhead_map', 'e3_map', 'linker_maps')})
            output['real_tool_execution_status'] = 'REFERENCE_DEMO' if context == 'BRD4–VHL reference demo' else 'FRESH_INFERENCE'
            output['execution_mode'] = output['real_tool_execution_status']
            metrics['ERBB2_specific_design'] = 'NOT_EVALUATED' if context != 'ERBB2' else 'ASSEMBLY_ONLY'
            for item in assembled['artifacts']:
                register_demo_artifact(record, item['path'], modality, 'assembly', artifact_type='structure', execution_mode=output['execution_mode'])
        output['lightweight_calculations'].append({'tool': 'RDKit', 'execution_mode': 'FRESH_ANALYSIS', 'scope': 'component_QC'})
    return output


def run_competition_demo(state, config):
    import json

    from .critic_agent import review_demo_modality
    from .run_records import RunRecord, demo_visualizations, register_demo_artifact
    from .small_molecule_io import now, write_json
    from .user_selection import check_design_gate
    from .validation_agent import validate_demo_modality

    source = Path(config.source_run).resolve()
    root_source = json.loads((source / 'run_summary.json').read_text())
    scientific = root_source.get('scientific_validation_status') or root_source.get('summary', {}).get('verdict')
    if root_source.get('is_mock') is not False or scientific != 'VERIFIED_REAL_PARTIAL':
        raise ValueError('verified_real_partial_source_required')
    record = RunRecord(state['output_dir'], 'FOUR_MODALITY_DEMO')
    results = []
    qualification = state.get('advance_targets') or [{'gene_name': 'ERBB2', 'tier': None,
        'qualification_status': 'source_target_context_only; qualification package not present in source run'}]
    record.data.update(execution_profile='competition_demo', selection_mode=config.selection_mode,
                       clinical_or_scientific_selection=False, scientific_validation_status=scientific,
                       demo_pipeline_status='FAILED', source_run_id=source.name)

    def event(stage, modality, status, tool, message, mode='REAL_ARTIFACT_REPLAY', paths=None, started=None):
        record.event(stage, status.upper(), modality=modality, tool=tool,
            execution_mode=mode, is_mock=False, message=message, artifact_paths=paths or [],
            started_at=started or now(), finished_at=now())

    event('qualification', 'ALL', 'partial', 'Qualification context',
          'Existing target context reused; unavailable upstream tier is not inferred.', mode='NOT_RUN')
    for modality, chain in DEMO_TOOL_CHAINS.items():
        if modality == 'DEGRADER' and config.degrader_reference_config and not config.degrader_component_config:
            chain = ['Curated reference selection', 'RDKit full molecule QC', 'ETKDG', 'MMFF/UFF', 'Ternary readiness']
        event('modality_evaluation', modality, 'completed', 'Modality assessment', 'Evaluate independent route for comparative demonstration.')
        reason = {'SMALL_MOLECULE': 'Compare existing docking, AF3 and independent ADMET evidence.',
                  'DE_NOVO_BINDER': 'Review the designed protein interface and sequence liabilities.',
                  'ADC': 'Assess the observed antibody interface and mapped conjugation candidates.',
                  'DEGRADER': 'Select a curated published full molecule unless explicit assembly components are supplied; no attachment inference.'}[modality]
        event('tool_planning', modality, 'completed', ' → '.join(chain), reason, mode='NOT_RUN')
        start = now()
        allowed = config.selection_mode == 'showcase_all'
        if not allowed:
            allowed = any(check_design_gate(t.get('candidate_id') or t.get('gene_name', ''), modality,
                t.get('tier') or t.get('best_tier') or 'UNKNOWN', state.get('user_selection_events', []))['allowed']
                and t.get('safety_verdict') not in {'REJECT', 'VETO'} for t in qualification)
        try:
            if not allowed:
                item = result(modality, 'NOT_RUN', 'selection', blockers=['explicit_user_selection_required'])
                item.update(demo_metrics={}, real_tool_execution_status='NOT_RUN', execution_mode='NOT_RUN',
                            missing_steps=['user_selection'], source_provenance=[])
            else:
                item = replay_demo_modality(modality, source, record, config)
        except Exception as exc:
            item = result(modality, 'BLOCKED', 'artifact_replay', blockers=[str(exc)])
            item.update(demo_metrics={}, real_tool_execution_status='NOT_RUN', execution_mode='NOT_RUN',
                        missing_steps=['valid_real_artifacts_or_inputs'], source_provenance=[])
        if allowed:
            event('readiness_check', modality, 'blocked' if item['execution_status'] == 'BLOCKED' else 'completed',
                  'Artifact provenance and inputs', '; '.join(item['blockers']) or 'Real artifacts and checksums verified.',
                  item['execution_mode'], started=start)
        if not allowed:
            event('readiness_check', modality, 'blocked', 'Selection gate', 'No model or replay executed before explicit user selection.', mode='NOT_RUN')
        mode = item['execution_mode']
        event('execution_or_replay', modality, item['execution_status'], ' → '.join(chain),
              '; '.join(item['blockers']) or 'Actual artifacts reused; no heavy model inference.', mode, started=start)
        validation = validate_demo_modality(item)
        item.update(validation=validation, validation_decision=validation['decision'])
        event('validation', modality, 'completed', 'Validation', validation['decision'], mode)
        item['critic'] = review_demo_modality(item)
        event('critic', modality, 'completed', 'Critic', item['critic']['decision'], mode)
        item.update(demo_pipeline_status='COMPLETE', modality_execution_status=item['execution_status'],
                    scientific_validation_status=scientific, critic_decision=item['critic']['decision'],
                    agent_reason=reason, selected_tools=chain, tool_order=chain,
                    readiness_result='blocked' if item['execution_status'] in {'BLOCKED', 'NOT_RUN'} else 'available',
                    execution_result=item['execution_status'], validation_result=validation['decision'],
                    critic_result=item['critic']['decision'])
        try:
            figures = demo_visualizations(record, item)
            item['visualizations'] = [v['path'] for v in figures]
            event('visualization', modality, 'completed', 'Matplotlib / RDKit',
                  'Evidence-based plots registered.' if figures else 'No available measurements; show status cards.', mode,
                  paths=item['visualizations'])
        except Exception as exc:
            item.update(demo_pipeline_status='FAILED', visualization_error=str(exc), visualizations=[])
            event('visualization', modality, 'failed', 'Plotting', str(exc), mode)
        directory = record.directory / modality.lower()
        directory.mkdir(exist_ok=True)
        path = write_json(directory / 'run_summary.json', {**item, 'run_id': directory.name,
                          'stages': [s for s in record.data['stages'] if s.get('modality') == modality],
                          'artifacts': [a for a in record.data['artifacts'] if a.get('modality') == modality]})
        register_demo_artifact(record, path, modality, 'report_registration', title=modality + ' demo summary')
        event('report_registration', modality, 'completed', 'RunRecord', 'Summary and artifacts registered.', mode, paths=[str(path)])
        results.append(item)
    output = {'stage': 'report', 'modality': 'FOUR_MODALITY_DEMO', 'execution_status': 'PARTIAL',
              'validation_decision': 'HOLD', 'is_mock': False, 'execution_mode': 'REAL_ARTIFACT_REPLAY',
              'real_tool_execution_status': 'REAL_ARTIFACT_REPLAY',
              'demo_pipeline_status': 'COMPLETE' if all(r['demo_pipeline_status'] == 'COMPLETE' for r in results) else 'FAILED',
              'scientific_validation_status': scientific, 'critic': {'decision': 'HOLD', 'scope': 'scientific evidence remains partial'},
              'summary': {'modalities': results, 'qualification': qualification,
                          'selection_mode': config.selection_mode, 'clinical_or_scientific_selection': False,
                          'heavy_model_inference_calls': 0, 'source_run_id': source.name},
              'blockers': list(dict.fromkeys(b for r in results for b in r['blockers'])),
              'missing_steps': list(dict.fromkeys(s for r in results for s in r['missing_steps']))}
    record.finish(output)
    try:
        from .server import write_competition_report
        report_path = write_competition_report(record.data, record.directory)
        register_demo_artifact(record, report_path, 'ALL', 'report', artifact_type='report', title='Competition demo HTML')
        record.data['report_path'] = str(report_path)
    except Exception as exc:
        record.data.update(demo_pipeline_status='FAILED', report_error=str(exc))
    write_json(record.directory / 'run_summary.json', record.data)
    return record.data


def select_degrader_reference(context, config_path=None):
    import json
    path = Path(config_path) if config_path else Path(__file__).parent / 'configs/degrader_sjf1528_reference.json'
    target = context.get('target') or context.get('target_id') or 'ERBB2'
    if target != 'ERBB2':
        raise ValueError('SJF1528_reference_selection_requires_ERBB2_context')
    data = json.loads(path.read_text())
    if data.get('candidate_id') != 'SJF1528' or data.get('candidate_origin') != 'published_reference':
        raise ValueError('published_SJF1528_reference_required')
    return {**data, 'selection_trace': [
        {'stage': 'target_context', 'mode': 'RUNTIME_SELECTION', 'message': 'ERBB2 target confirmed.'},
        {'stage': 'input_evidence', 'mode': 'RUNTIME_SELECTION', 'message': 'TAK-285 has no curated degrader exit-vector evidence; no linker or attachment inferred.'},
        {'stage': 'candidate_selection', 'mode': 'RUNTIME_SELECTION', 'message': 'Selected SJF1528 from local curated reference config; literature research occurred at ingestion, not in this runtime.'},
        {'stage': 'input_mode', 'mode': 'RUNTIME_SELECTION', 'message': 'Use verified full molecule; component assembly is not applicable.'}]}


def run_degrader_reference(directory, config=None):
    from .critic_agent import review_modality
    from .run_records import RunRecord, register_demo_artifact
    from .small_molecule_io import write_json
    from .validation_agent import validate_modality

    config = dict(config or {})
    data = select_degrader_reference({'target': 'ERBB2'}, config.pop('reference_config', None))
    record = RunRecord(directory, 'DEGRADER')
    write_json(record.directory / 'input.json', data)
    record.event('curated_ingestion', 'COMPLETED', execution_mode='CURATED_INGESTION',
                 message=data['curation_message'], source_records=data.get('source_records', []))
    for item in data['selection_trace']:
        record.event(item['stage'], 'COMPLETED', execution_mode=item['mode'], message=item['message'])
    record.event('full_molecule_preparation', 'NOT_RUN', event='started', execution_mode='FRESH_ANALYSIS')
    output = execute_degrader(data, {**config, 'output_dir': str(record.directory)})
    output['validation'] = validate_modality(output)
    output['validation_decision'] = output['validation']['decision']
    output['critic'] = review_modality(output)
    output.update(execution_mode='FRESH_ANALYSIS', scientific_validation_status='VERIFIED_REAL_PARTIAL',
                  demo_pipeline_status='COMPLETE', decision_trace=data['selection_trace'])
    for name in ('validation', 'critic'):
        write_json(record.directory / (name + '.json'), output[name])
        record.event(name, 'COMPLETED', decision=output[name]['decision'])
    output['artifacts'] = []
    for path in sorted(record.directory.iterdir()):
        if path.is_file() and path.name not in {'run_summary.json', 'events.jsonl'}:
            register_demo_artifact(record, path, 'DEGRADER', 'reference_evaluation',
                artifact_type='plot' if path.suffix == '.png' else 'structure' if path.suffix == '.sdf' else 'table',
                execution_mode='FRESH_ANALYSIS', source_table=record.directory / 'rdkit_validation.json')
    output['artifacts'] = record.data['artifacts']
    return record.finish(output)
