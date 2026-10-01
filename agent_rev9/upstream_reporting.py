"""Persist upstream analysis results and render a portable report."""

from .run_records import register_demo_artifact
from .configuration import get_config
from .small_molecule_io import write_json


def finish_report(record, state, terminal):
    r = state.get('r_analysis_result', {})
    candidates = []
    for gene, result in state.get('qualification_results', {}).items():
        candidates.append({'gene': gene, 'B': result['bio_confidence'], 'tier': result['best_tier'],
                           'decision': result['shortlist']['judgment'], 'reason': result['shortlist']['reason'],
                           'routes': result['tier_results']})
    summary = {'question': record.data['question'], 'de': r.get('de_summary', {}),
        'gsea': r.get('gsea_summary', {}), 'survival': r.get('survival_summary', {}),
        'input_manifest': r.get('input_manifest', {}), 'r_execution_mode': r.get('execution_mode', 'NOT_RUN'),
        'r_source_output_dir': r.get('source_output_dir', r.get('output_dir')), 'r_cache_key': r.get('cache_key'),
        'cell_context': state.get('cell_context_result', {}), 'depmap': state.get('depmap_result', {}),
        'candidates': candidates, 'contrasts': state.get('contrasts', []),
        'modality_assessments': state.get('modality_assessments', []),
        'selection': {'gate_status': state.get('gate_status'), 'events': state.get('user_selection_events', [])},
        'design_mode': state.get('design_mode', 'none'),
        'dossier': {'status': state.get('dossier_status', 'NOT_RUN'), 'path': state.get('dossier_path')},
        'candidate_search': state.get('candidate_search', {}),
        'design_outcome': state.get('design_outcome', {}),
        'reused_stages': state.get('reused_stages', []),
        'review_issues': terminal.get('review_issues', []),
        'design_reference_run': get_config().upstream.report_reference_run if state.get('design_mode') == 'all_eligible' else '',
        'modalities': state.get('automatic_design_results', []),
        'agent_notes': state.get('agent_notes', {}),
        'limitations': ['TNBC vs Non-TNBC association; no normal-tissue DE contrast or external replication',
                       'Census cell annotations alone do not establish malignant cell origin',
                       'This TNBC reference is one ICB-treated patient dataset; treatment and sampling limit generalization',
                       'Breast lineage dependency is not TNBC-specific without receptor metadata',
                       'Automatic Design executes eligible routes under the recorded bulk user selection; missing inputs remain blocked']}
    from .resource_usage import report_usage
    usage = report_usage(record.directory, [r['run_id'] for r in summary['modalities'] if r.get('run_id')])
    if usage is not None:
        summary['resource_usage'] = usage
    record.data.update(summary=summary, execution_profile='upstream_analysis',
                       execution_mode=('MIXED_VERIFIED_REUSE_AND_FRESH_ANALYSIS' if r.get('execution_mode') == 'VERIFIED_CACHE_REUSE' else r.get('execution_mode', 'NOT_RUN')))
    if state.get("cell_contexts") or state.get("depmap_results"):
        try:
            plot_evidence(record, state)
        except Exception as exc:
            record.data.setdefault('report_warnings', []).append('Evidence figures: ' + str(exc))
    sync_artifacts(record, r.get("execution_mode", "FRESH_ANALYSIS"))
    execution = terminal['execution_status']
    normalized = ('BLOCKED' if execution.startswith('BLOCKED') else 'PARTIAL' if execution in
                  {'COMPLETED_WITH_ERRORS', 'PARTIAL'} else execution)
    record.data.update(running=False, workflow_status='FINISHED', report_status='READY')
    record.finish({'stage': 'complete', 'execution_status': normalized,
                   'validation_decision': 'HOLD' if candidates else 'NOT_EVALUATED',
                   'blockers': terminal.get('blockers', []), 'errors': terminal.get('errors', [])})
    write_html(record)
    register_demo_artifact(record, record.directory / 'report.html', 'UPSTREAM', 'report',
                           artifact_type='report', title='Standalone HTML report', execution_mode='REPORT_RENDER')
    write_json(record.directory / 'run_summary.json', record.data)


def sync_artifacts(record, execution_mode="FRESH_ANALYSIS"):
    record.data["artifacts"] = [a for a in record.data["artifacts"] if a.get("artifact_type") == "report"]
    allowed = {'.png', '.svg', '.tsv', '.csv', '.json', '.jsonl', '.log', '.txt'}
    for path in sorted((record.directory / 'agent_results').rglob('*')):
        if path.is_file() and path.suffix in allowed:
            register_demo_artifact(record, path, 'UPSTREAM', path.parent.name,
                artifact_type='plot' if path.suffix in {'.png', '.svg'} else 'table',
                title=str(path.relative_to(record.directory / 'agent_results')),
                execution_mode=execution_mode if "r_analysis" in path.parts else "FRESH_ANALYSIS")


def write_html(record):
    from .server import write_competition_report
    return write_competition_report(record.data, record.directory)


def plot_evidence(record, state):
    import csv

    import matplotlib
    import numpy as np
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    directory=record.directory/'agent_results/visualizations'
    directory.mkdir(parents=True,exist_ok=True)
    contexts=state.get('cell_contexts', {})
    if contexts:
        genes=sorted(g for g,r in contexts.items() if r.get('all_classes'))
        classes=sorted({c for g in genes for c in contexts[g]['all_classes']})
        if genes and classes:
            matrix=np.array([[contexts[g]['all_classes'].get(c,np.nan) for c in classes] for g in genes])
            fig,ax=plt.subplots(figsize=(10,max(4,len(genes)*.22)))
            image=ax.imshow(matrix,aspect='auto',cmap='viridis')
            ax.set(yticks=range(len(genes)),yticklabels=genes,xticks=range(len(classes)),xticklabels=classes,
                   title='TNBC reference: genesorteR specificity (sampled donors)')
            plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
            fig.colorbar(image,ax=ax,label='Specificity score')
            fig.savefig(directory/'cell_specificity.png',dpi=140,bbox_inches='tight')
            plt.close(fig)
            with (directory/'cell_specificity.tsv').open('w',newline='') as handle:
                writer=csv.writer(handle,delimiter='\t')
                writer.writerow(['gene',*classes])
                writer.writerows([gene,*values] for gene,values in zip(genes,matrix))
    dependencies=state.get('depmap_results', {})
    if dependencies:
        genes=sorted(dependencies)
        fig,axes=plt.subplots(1,2,figsize=(10,max(4,len(genes)*.22)),sharey=True,sharex=True)
        rows=[]
        for ax,key,label in zip(axes,['crispr_summary','rnai_summary'],['CRISPR Chronos','RNAi DEMETER2']):
            for i,gene in enumerate(genes):
                values=dependencies[gene].get(key) or {}
                value=values.get('median_effect')
                if value is not None:
                    ax.barh(i,value,color='#3d98a7' if key=='crispr_summary' else '#d99748')
                    rows.append([gene,label,value,values.get('cohort_scope','UNKNOWN'),values.get('n_models',0)])
                else:
                    ax.text(.98,i,'not measured',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7,color='gray')
            ax.axvline(0,color='gray',linewidth=.6)
            ax.set(title=label,xlabel='Median gene effect')
            ax.set(yticks=range(len(genes)),yticklabels=genes)
        fig.suptitle('Breast-lineage dependency; TNBC subtype is unconfirmed')
        fig.savefig(directory/'dependency_effects.png',dpi=140,bbox_inches='tight')
        plt.close(fig)
        with (directory/'dependency_effects.tsv').open('w',newline='') as handle:
            writer=csv.writer(handle,delimiter='\t')
            writer.writerow(['gene','technology','median_effect','cohort_scope','n_models'])
            writer.writerows(rows)
