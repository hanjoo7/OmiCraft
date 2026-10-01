"""Execute eligible routes selected by an explicit full-pipeline request."""
import copy
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .configuration import configuration_scope, get_config
from .routes import MODALITY_ROUTES
from .small_molecule_io import write_json
from .user_selection import check_design_gate, create_selection_event


def select_routes(state):
    events = list(state.get('user_selection_events', []))
    for target in state.get('advance_targets', []):
        gene = target.get('gene_name') or target.get('candidate_id')
        for route in target.get('tier_results', []):
            modality, tier = route.get('modality'), route.get('Tier')
            if modality not in MODALITY_ROUTES or tier not in {'1A', '1B', '2A', '2B'}:
                continue
            if target.get('safety_verdict') in {'VETO', 'REJECT'} or target.get('safety_veto', {}).get(modality, {}).get('has_veto'):
                continue
            key = gene + ':' + modality
            if any(e.get('route_id') == key and e.get('status') == 'SELECTED' and e.get('selected_by') == 'user' for e in events):
                continue
            events = [event for event in events if event.get('route_id') != key]
            event = create_selection_event(gene, modality, key, tier, selected_by='user',
                rationale='Run Pipeline request: execute all eligible routes',
                constraints={'design_mode': 'all_eligible', 'parallel_workers': 4, 'gpu_devices': [4, 5, 6, 7], 'binder_backbones': 1})
            event['selection_scope'] = 'all_eligible'
            events.append(event)
    return events


def route_input(gene, modality, supplied=None, uniprot=None):
    from .web_execution import assets, pdb_text
    if supplied:
        return copy.deepcopy(supplied)
    if modality == 'SMALL_MOLECULE':
        raise ValueError('ligand_smiles_and_docking_box_required')
    if modality == 'DEGRADER':
        raise ValueError('warhead_e3_linker_and_attachment_maps_required')
    rows = [row for row in assets() if row['gene'] == gene and row['kind'] == 'protein']
    if modality == 'ADC':
        rows = [row for row in rows if row.get('role') == 'antibody_complex']
        if not rows:
            raise ValueError('antibody_antigen_complex_and_antigen_chain_required')
        raise ValueError('explicit_antigen_chain_required')
    # Canonical single-chain predictions have an unambiguous target sequence.
    rows = [row for row in rows if row.get('source_type') == 'alphafold_db'
            and row.get('uniprot_id') and '-' not in row['uniprot_id']]
    if not rows:
        from .target_structure import ensure_target_structure
        rows = [ensure_target_structure(gene, uniprot)]
    row = sorted(rows, key=lambda item: item['id'])[0]
    from Bio.PDB import PDBParser
    from Bio.Data.PDBData import protein_letters_3to1
    import io
    model = PDBParser(QUIET=True).get_structure('target', io.StringIO(pdb_text(row)))[0]
    chains = [chain for chain in model if any(r.resname in protein_letters_3to1 for r in chain)]
    if len(chains) != 1:
        raise ValueError('ambiguous_target_chain: explicit selection required')
    return {'asset_id': row['id'], 'target_chain': chains[0].id,
            'design_length': 100, 'hotspot_residues': []}


def _blocked_report(directory, parent, gene, modality, status, reason, selection):
    from . import server
    from .run_records import RunRecord
    record = RunRecord(directory, modality)
    record.data.update(execution_profile='therapeutic_design', parent_run_id=parent,
                       gene=gene, execution_mode='NOT_RUN', running=False,
                       selection_event=selection, workflow_status='FINISHED', report_status='READY')
    record.finish({'stage': 'readiness', 'execution_status': status,
                   'validation_decision': 'NOT_EVALUATED', 'blockers': [reason],
                   'critic': {'decision': 'NOT_EVALUATED'}})
    server.write_competition_report(record.data, directory)


def automatic_design_node(state):
    if state.get('design_mode') != 'all_eligible' or state.get('dry_run'):
        return {}
    from .web_execution import design_input
    from .design_execution import run_bounded_design as run_design
    from .modality_dispatch import readiness
    from .orchestrator import assess_modalities

    cfg = get_config().model_copy(deep=True)
    directory = Path(state.get('_run_directory') or cfg.data.base_dir).resolve()
    if not re.fullmatch(r'upstream_\d{8}T\d{6}Z_[0-9a-f]{8}', directory.name):
        directory = directory / ('upstream_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_') + uuid.uuid4().hex[:8])
    directory.mkdir(parents=True, exist_ok=True)
    app_config = cfg.model_copy(deep=True)
    app_config.data.base_dir = str(directory.parent)
    selected = select_routes(state)
    working = {**state, 'user_selection_events': selected}
    snapshot = {k: v for k, v in working.items() if not k.startswith('_')}
    (directory / 'pipeline_state.json').write_text(json.dumps(snapshot, ensure_ascii=False, indent=2, default=str))
    from .runtime_paths import config_path, load_json
    input_file = config_path('automatic_design_inputs', 'OMICRAFT_DESIGN_INPUTS_CONFIG')
    configured = load_json(input_file).get('routes', {}) if input_file.is_file() else {}
    rows, messages = [], []
    def execute_route(target, modality):
        gene = target.get('gene_name') or target.get('candidate_id')
        key = gene + ':' + modality
        route = next((r for r in target.get('tier_results', []) if r.get('modality') == modality), {})
        tier = route.get('Tier', 'UNKNOWN')
        selection = next((e for e in selected if e.get('route_id') == key), None)
        gate = check_design_gate(gene, modality, tier, selected)
        if target.get('safety_verdict') in {'VETO', 'REJECT'} or target.get('safety_veto', {}).get(modality, {}).get('has_veto'):
            gate = {'allowed': False, 'gate_status': 'BLOCKED_REJECT', 'reason': 'target_safety_veto'}
        child_id = 'design_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_') + uuid.uuid4().hex[:8]
        child_dir = directory.parent / child_id
        row = {'gene': gene, 'modality': modality, 'tier': tier, 'run_id': child_id,
               'report_url': '/report?run_id=' + child_id, 'execution_status': 'BLOCKED',
               'execution_mode': 'NOT_RUN', 'validation_status': 'NOT_EVALUATED',
               'critic_decision': 'NOT_EVALUATED', 'is_mock': False, 'blocker': ''}
        started = False
        try:
            if not gate['allowed']:
                raise ValueError(gate['gate_status'] + ': ' + gate['reason'])
            with configuration_scope(app_config):
                supplied = configured.get(key)
                callback = state.get('_design_progress_callback')
                if callback:
                    callback(key, 'RUNNING', '[Design] ' + key + ': preparing inputs and target structure', row['report_url'])
                uniprot = (state.get('qualification_results') or {}).get(gene, {}).get('uniprot')
                body = {'parent_run_id': directory.name, 'gene': gene, 'modality': modality,
                        'input': route_input(gene, modality, supplied, uniprot), 'selection_scope': 'all_eligible',
                        'reuse_results': state.get('reuse_results', True)}
                data, model_cfg, qualification = design_input(body)
                blockers = readiness(modality, data, model_cfg)['blockers']
                if blockers:
                    raise ValueError('; '.join(blockers))
            started = True
            result = run_design(child_id, body, data, model_cfg, qualification, app_config, finalize=False)
            row.update(execution_status=result['execution_status'], execution_mode=result.get('execution_mode', 'FRESH_ANALYSIS'),
                       validation_status=result.get('scientific_validation_status', 'NOT_EVALUATED'),
                       blocker='; '.join(result.get('errors', [])))
            summary = child_dir / 'run_summary.json'
            if summary.is_file():
                child = json.loads(summary.read_text())
                row['critic_decision'] = child.get('critic', {}).get('decision', 'NOT_EVALUATED')
                row['gpu_device'] = child.get('gpu_device')
                row['is_mock'] = bool(child.get('is_mock', True))
        except (ValueError, FileNotFoundError) as exc:
            row['blocker'] = str(exc)
        except Exception as exc:
            row.update(execution_status='FAILED', blocker=str(exc))
        if not started or not (child_dir / 'run_summary.json').is_file():
            _blocked_report(child_dir, directory.name, gene, modality, row['execution_status'], row['blocker'], selection)
        message = f'[Design] {key}: {row["execution_status"]}' + (' · ' + row['blocker'] if row['blocker'] else '')
        return row, message

    jobs = [(target, modality) for target in state.get('advance_targets', []) for modality in MODALITY_ROUTES]
    with ThreadPoolExecutor(max_workers=4, thread_name_prefix='design') as pool:
        pending = {pool.submit(execute_route, target, modality): index
                   for index, (target, modality) in enumerate(jobs)}
        completed = {}
        for future in as_completed(pending):
            row, message = future.result()
            completed[pending[future]] = row
            rows = [completed[index] for index in sorted(completed)]
            messages.append(message)
            callback = state.get('_design_progress_callback')
            if callback:
                callback(row['gene'] + ':' + row['modality'], row['execution_status'], message, row['report_url'])
            write_json(directory / 'agent_results/automatic_design/results.json',
                       {'design_mode': 'all_eligible', 'parallel_workers': 4, 'routes': rows})
    if not rows:
        messages.append('[Design] No ADVANCE candidates; no model jobs started')
    allowed = [row for row in rows if any(e['route_id'] == row['gene'] + ':' + row['modality'] for e in selected)]
    status = 'COMPLETED' if allowed and all(r['execution_status'] == 'COMPLETED' for r in allowed) else 'PARTIAL'
    successful = [r for r in rows if r['execution_status'] == 'COMPLETED'
                  and r['validation_status'] == 'ADVANCE' and r['critic_decision'] == 'ADVANCE' and not r['is_mock']]
    outcome = {'completed_routes': sum(r['execution_status'] == 'COMPLETED' for r in rows),
               'validated_successes': len(successful), 'success_run_ids': [r['run_id'] for r in successful],
               'blocked_routes': sum(r['execution_status'] == 'BLOCKED' for r in rows),
               'failed_routes': sum(r['execution_status'] == 'FAILED' for r in rows),
               'status': 'VALIDATED_CANDIDATES_FOUND' if successful else 'NO_VALIDATED_CANDIDATE'}
    return {'design_outcome': outcome, 'automatic_design_results': rows, 'user_selection_events': selected,
            'modality_assessments': assess_modalities(working)['modality_assessments'],
            'gate_status': 'ALLOWED' if allowed else 'HOLD', 'execution_status': status,
            'design_mode': 'all_eligible', 'messages': ([f'[Design] {len(rows)}개 경로 처리 종료 · {status}'] if state.get('_design_progress_callback') else messages),
            'design_run_directory': str(directory), 'current_agent': 'design'}
