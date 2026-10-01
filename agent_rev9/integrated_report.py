"""Assemble an ERBB2 report from recorded child runs without changing their verdicts."""
import copy
from pathlib import Path

MODALITIES = ('SMALL_MOLECULE', 'DE_NOVO_BINDER', 'ADC', 'DEGRADER')
TITLES = dict(zip(MODALITIES, ('Small molecule', 'De novo binder', 'ADC', 'Degrader')))


def _table(title, rows, columns):
    return {'title': title, 'columns': columns, 'rows': rows}


def _step(name, status, message=''):
    return {'stage': name, 'execution_status': status, 'message': message}


def _status(value):
    return {'success': 'COMPLETED', 'MPNN_COMPLETED': 'COMPLETED', 'AF3_COMPLETED': 'COMPLETED',
            'prediction_failed': 'FAILED', 'backend_failed': 'FAILED', 'not_run': 'NOT_RUN',
            'invalid_input': 'BLOCKED', 'NOT_CONFIGURED': 'BLOCKED', 'EXECUTED': 'COMPLETED'}.get(value, value or 'NOT_RUN')


def panel(child, row):
    modality = row['modality']
    s = child.get('summary', {})
    running = bool(child.get('running', row.get('running', False)))
    p = {'modality': modality, 'title': TITLES[modality], 'run_id': row.get('run_id'),
         'running': running, 'execution_status': 'RUNNING' if running else child.get('execution_status', row.get('execution_status', 'NOT_RUN')),
         'validation': child.get('validation_decision', row.get('validation_status', 'NOT_EVALUATED')),
         'critic': child.get('critic', {}).get('decision', row.get('critic_decision', 'NOT_EVALUATED')),
         'execution_mode': child.get('execution_mode', row.get('execution_mode', 'NOT_RUN')),
         'is_mock': child.get('is_mock', row.get('is_mock', False)),
         'scope': child.get('execution_scope', s.get('execution_scope', 'computational_screening')),
         'blockers': child.get('blockers', [row['blocker']] if row.get('blocker') else []),
         'missing_steps': child.get('missing_steps', []), 'input_status': row.get('input_status'),
         'computational_validation': child.get('computational_validation', child.get('validation', {}).get('computational_validation', {})),
         'experimental_validation': child.get('experimental_validation', 'NOT_EVALUATED'), 'stages': [], 'tables': [], 'highlights': [],
         'artifacts': child.get('artifacts', []), 'events': child.get('stages', []),
         'agent_notes': s.get('agent_notes', {}), 'structures': [], 'figures': [],
         'report_url': '/report?run_id=' + row['run_id'] if row.get('run_id') else None}
    candidates = s.get('candidates', [])
    wanted = []
    if modality == 'SMALL_MOLECULE':
        p['scope'] = 'Ligand preparation → docking → AF3 → ADMET → review'
        p['stages'] = list(s.get('execution_checklist', []))
        for c in candidates:
            cid = c.get('candidate_id', '')
            p['highlights'].append({'label': cid, 'value': c.get('final_verdict', 'NOT_EVALUATED')})
            for title, values, columns in [
                ('Docking', c.get('docking', {}).get('poses', []), ['pose_rank','docking_score','cnn_score','cnn_affinity']),
                ('Pose QC', c.get('pose_qc', []), ['pose_rank','pose_qc_status','contact_residue_count','severe_clash_count']),
                ('AF3', c.get('af3_ligand', {}).get('samples', []), ['seed','sample_id','af3_status','iptm','ligand_atom_plddt_mean','protein_ligand_chain_pair_iptm']),
                ('Docking / AF3 comparison', c.get('consensus', {}).get('comparisons', []), ['docking_pose_rank','docking_af3_pose_rmsd','contact_residue_jaccard','consensus_status']),
                ('ADMET', c.get('admet', {}).get('endpoints', []), ['endpoint_name','endpoint_value','endpoint_unit','prediction_status'])]:
                p['tables'].append(_table(cid + ' · ' + title, values, columns))
            p['blockers'] += c.get('critic', {}).get('reason_codes', [])
            wanted += [v.get('model_path') for v in c.get('af3_ligand', {}).get('samples', [])]
            p['highlights'].append({'label':'ADMET endpoints', 'value':c.get('admet', {}).get('available_endpoint_count')})
    elif modality == 'DE_NOVO_BINDER':
        p['scope'] = 'RFD3 → ProteinMPNN → AF3 → interface screening'
        targets = s.get('targets', [])
        p['stages'] = [_step('backbone_import', 'COMPLETED') if target.get('design_protocol') == 'sequence_redesign' else _step('rfd3', _status(target.get('rfd3', {}).get('status'))) for target in targets]
        if s.get('design_protocol') == 'sequence_redesign':
            p['scope'] = 'Provided complex → ProteinMPNN sequence redesign → AF3 → interface screening'
        for c in candidates:
            metrics, geometry = c.get('af3_metrics', {}), c.get('structural_metrics', {})
            p['highlights'].extend([{'label':c.get('candidate_id','Binder'), 'value':c.get('final_verdict','NOT_EVALUATED')},
                                    {'label':'iPTM', 'value':metrics.get('iptm')}, {'label':'pTM', 'value':metrics.get('ptm')},
                                    {'label':'Interface contacts', 'value':geometry.get('interface_contact_count')}])
            p['tables'].append(_table('Binder candidates', [{'candidate_id':c.get('candidate_id'),
                'verdict':c.get('final_verdict'), 'iptm':metrics.get('iptm'),'ptm':metrics.get('ptm'),
                'binder_plddt':metrics.get('plddt_summary',{}).get('by_chain_mean',{}).get('B'),
                'contacts':geometry.get('interface_contact_count'),'hotspot_preserved':geometry.get('hotspot_preserved'),
                'sequence':c.get('sequence')}], ['candidate_id','verdict','iptm','ptm','binder_plddt','contacts','hotspot_preserved','sequence']))
            p['tables'].append(_table('AF3 samples', metrics.get('samples', []), ['seed','sample_id','iptm','ptm','has_clash','chain_pair_iptm','chain_pair_pae_min']))
            p['stages'] += [_step('protein_mpnn',_status(c.get('mpnn_status'))),_step('af3',_status(c.get('af3_status'))),
                            _step('interface_qc',_status(geometry.get('status')))]
            p['blockers'] += c.get('failure_reasons', [])
            wanted += [metrics.get('model_path'), c.get('rfd3_structure_path')]
        for target in targets:
            p['blockers'] += target.get('failure_reasons', [])
        analysis = s.get('binder_analysis', {})
        if analysis:
            p['stages'].append(_step('binder_analysis',analysis.get('status','NOT_RUN')))
            p['highlights'] += [{'label':'Additional quality-filter passes','value':analysis.get('quality_pass_count')},
                                {'label':'Unique sequences','value':analysis.get('unique_sequence_count')}]
            p['tables'].append(_table('Additional quality filters and ranking · primary verdict unchanged', analysis.get('ranking',[]),
                ['rank','candidate_id','original_verdict','quality_status','binder_plddt','i_pae','binder_rmsd','binder_internal_rmsd','complex_rmsd','duplicate_of','filter_reasons','errors']))
            p['tables'].append(_table('Analysis settings', [{'setting':k,'value':v} for k,v in analysis.get('settings',{}).items()],['setting','value']))
            p['analysis_notes'] = analysis.get('warnings', [])
            wanted += [analysis.get('best_structure_path')]
    elif modality == 'ADC':
        p['stages'] = list(s.get('execution_checklist', []))
        p['highlights'] = [{'label':key, 'value':s.get(key)} for key in ['cdr_mapping_status','cdr_heavy_atom_contacts','severe_clashes']]
        p['tables'] += [_table('Conjugation candidates', s.get('conjugation',{}).get('candidates',[]),
                            ['review_rank','chain','residue_id','residue','sasa_angstrom2','cdr_distance_angstrom','exclusion_reasons']),
                        _table('CDR contacts',s.get('contacts',[]),['antibody_chain','antibody_residue_id','antigen_chain','antigen_residue_id','minimum_distance'])]
        wanted = [s.get('structure_path')]
    else:
        from .degrader_status import ternary_execution_status
        ternary_status = ternary_execution_status(s)
        p['stages'] = [dict(step) for step in s.get('execution_checklist', [])]
        for step in p['stages']:
            if step['stage'] == 'ternary_prediction':
                step['execution_status'] = ternary_status
        p['highlights'] = [{'label':'Reference molecule','value':s.get('candidate_id')},
                           {'label':'Chemical preparation','value':s.get('full_molecule_validation')},
                           {'label':'Input readiness','value':s.get('ternary_readiness',{}).get('status','NOT_RUN')},
                           {'label':'Ternary prediction','value':ternary_status}]
        p['tables'] = [_table('Chemical properties',[s['qc']] if s.get('qc') else [],
                             ['molecular_formula','molecular_weight','logp','tpsa','conformer_generation_status','conformer_identity_preserved']),
                       _table('Ternary inputs',s.get('ternary_readiness',{}).get('sequences',[]),['entity','available','chain','sequence_length']),
                       _table('AF3 ternary samples',s.get('ternary_prediction',{}).get('samples',[]),['seed','sample_id','af3_status','iptm','ptm','chain_pair_metrics','error_message'])]
        p['tables'].append(_table('Protein MSA preparation', s.get('ternary_prediction',{}).get('protein_msa_cache',[]),
            ['chain','sequence_length','status','cache_status','stored','error']))
        wanted = [v.get('model_path') for v in s.get('ternary_prediction',{}).get('samples',[])]
    latest = {}
    for event in child.get('stages', []):
        if event.get('stage') not in {'readiness','therapeutic_design','execution'}:
            latest[event['stage']] = event
    existing = {step['stage'] for step in p['stages']}
    p['stages'] += [event for key,event in latest.items() if key not in existing]
    if not p['stages']:
        p['stages'] = list(latest.values()) or [_step('readiness',p['execution_status'],'; '.join(p['blockers']))]
    for stage in p['stages']:
        if stage.get('execution_status') in {'RUNNING','QUEUED'} and not running:
            stage['execution_status'] = 'INTERRUPTED'
    if not running and child:
        p['stages'].append(_step('report', 'COMPLETED'))
    by_path = {str(Path(a['path']).resolve()): a for a in p['artifacts'] if isinstance(a,dict) and a.get('path')}
    for path in wanted:
        if path and (a := by_path.get(str(Path(path).resolve()))) and a.get('available') and not str(path).endswith(('.gz','.zst')):
            p['structures'].append({**a,'format':'pdb' if str(path).endswith('.pdb') else 'cif'})
    for a in p['artifacts']:
        name = a.get('path','').lower()
        if a.get('available') and name.endswith(('.png','.jpg','.jpeg','.svg')):
            p['figures'].append(a)
        if modality == 'DEGRADER' and a.get('available') and name.endswith('conformer.sdf'):
            p['structures'].append({**a, 'format':'sdf'})
    p['blockers'] = list(dict.fromkeys(str(v) for v in p['blockers'] if v))
    return p


def assemble_report(data, root):
    from .run_records import read_report
    out = copy.deepcopy(data)
    rows = {row['modality']: row for row in data.get('summary',{}).get('modalities',[]) if row.get('modality') in MODALITIES}
    panels = []
    for modality in MODALITIES:
        row = rows.get(modality, {'modality':modality, 'execution_status':'NOT_RUN'})
        child = {}
        if row.get('run_id'):
            try:
                child = read_report(root, run_id=row['run_id'], expand_integrated=False)
                if child.get('modality') != modality:
                    child = {}
                    raise ValueError('child_modality_mismatch')
            except (OSError, ValueError, TypeError, KeyError) as exc:
                row = {**row, 'blocker':'Child report unavailable: '+str(exc)}
        if not child and not data.get('running') and row.get('running'):
            row = {**row, 'running':False, 'execution_status':'FAILED',
                   'blocker':row.get('blocker') or 'Execution ended before a child report was saved'}
        panels.append(panel(child, row))
    running = bool(data.get('running')) or any(p['running'] for p in panels)
    done = sum(bool(p['run_id']) and not p['running'] and p['execution_status'] not in {'NOT_RUN','RUNNING','QUEUED'} for p in panels)
    passes = [{'modality':p['modality'], 'run_id':p['run_id'], 'report_url':p['report_url'], 'metrics':p['highlights']}
              for p in panels if p['validation']=='ADVANCE' and not p['is_mock']]
    out['integrated_report'] = {'panels':panels, 'completed_routes':done, 'total_routes':4,
        'workflow_status':'RUNNING' if running else 'FINISHED' if done==4 else 'INTERRUPTED',
        'report_status':'LIVE' if running else 'READY', 'passing_candidates':passes,
        'scope':'Computational demonstration; execution completion and candidate assessment are separate.',
        'previous_report_url':'/report?run_id=visual_diagnostic_20260922T152713Z'}
    return out


def _single_modality_report(data, root, field, title):
    from urllib.parse import urlencode
    out = copy.deepcopy(data)
    directory = Path(root) / data['run_id']
    for index, artifact in enumerate(out.get('artifacts', [])):
        path = Path(artifact['path'])
        path = (path if path.is_absolute() else directory / path).resolve()
        available = path.is_file() and path.is_relative_to(Path(root).resolve())
        artifact.update(available=available, href='/api/report?' + urlencode({'run_id':data['run_id'],'artifact':index}) if available else None)
    out[field] = panel(out, {'modality':data['modality'],'run_id':data['run_id']})
    out[field]['title'] = (data.get('gene') or 'Target') + ' · ' + title
    return out


def assemble_binder_report(data, root):
    return _single_modality_report(data, root, 'binder_report', 'De novo binder / redesign')


def assemble_degrader_report(data, root):
    return _single_modality_report(data, root, 'degrader_report', 'Degrader')
