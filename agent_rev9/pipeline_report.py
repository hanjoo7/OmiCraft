"""Attach recorded Design outputs to a target-analysis report."""
import copy
import csv
import math
import json
from pathlib import Path

from .integrated_report import MODALITIES, panel


def _panels(rows, root):
    from .run_records import read_report
    result = []
    for row in rows:
        if row.get('modality') not in MODALITIES:
            continue
        child = {}
        if row.get('run_id'):
            try:
                child = read_report(root, run_id=row['run_id'], expand_integrated=False)
                if child.get('modality') != row['modality']:
                    child = {}
                    raise ValueError('child_modality_mismatch')
            except (OSError, ValueError, TypeError, KeyError) as exc:
                row = {**row, 'blocker': 'Child report unavailable: ' + str(exc)}
        item = panel(child, row)
        item['gene'] = row.get('gene', child.get('gene', 'Unknown target'))
        item['candidate_count'] = len(child.get('summary', {}).get('candidates', []))
        result.append(item)
    return result


def _volcano_data(data, root):
    """Read registered observations; changing the viewport never changes DE values."""
    root = Path(root).resolve()
    for artifact in data.get('artifacts', []):
        source = artifact.get('path', '')
        if not source.endswith('/volcano_data.tsv'):
            continue
        path = Path(source)
        if not path.is_absolute():
            path = root / str(data.get('run_id', '')) / path
        path = path.resolve()
        if not path.is_relative_to(root):
            continue
        points, skipped = [], 0
        try:
            with path.open(newline='') as stream:
                for row in csv.DictReader(stream, delimiter='\t'):
                    try:
                        x, y = float(row['log2FoldChange']), float(row['minus_log10_padj'])
                        if not math.isfinite(x) or not math.isfinite(y) or y < 0:
                            raise ValueError('invalid coordinate')
                    except (ValueError, KeyError, TypeError):
                        skipped += 1
                        continue
                    points.append([row.get('gene_symbol') or row.get('gene_label') or row.get('gene_id'),
                                   x, y, row.get('direction', 'Not_DEG')])
        except (OSError, ValueError, csv.Error):
            continue
        if points:
            return {'points': points, 'source': str(path.relative_to(root)), 'skipped': skipped}
    return None


def assemble_pipeline_report(data, root):
    from .run_records import read_report
    out = copy.deepcopy(data)
    summary = out.setdefault('summary', {})
    # Older runs saved route results on disk but omitted them from their summary.
    directory = Path(root).resolve() / str(out.get('run_id', ''))
    def saved(relative, default):
        path = (directory / relative).resolve()
        if not out.get('run_id') or not path.is_relative_to(Path(root).resolve()):
            return default
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError):
            return default
    if not summary.get('modalities'):
        routes = saved('agent_results/automatic_design/results.json', {}).get('routes', [])
        if routes:
            summary['modalities'] = routes
            out['route_source'] = 'agent_results/automatic_design/results.json'
    from .resource_usage import report_usage
    if directory.resolve().is_relative_to(Path(root).resolve()):
        usage = report_usage(directory, [row['run_id'] for row in summary.get('modalities', []) if row.get('run_id')])
        if usage is not None:
            out['resource_usage'] = summary['resource_usage'] = usage
    evidence = {}
    path = directory / 'agent_results/qualification/evidence_cards.jsonl'
    if out.get('run_id') and path.resolve().is_relative_to(Path(root).resolve()):
        try:
            for line in path.read_text().splitlines():
                try:
                    card = json.loads(line)
                    if isinstance(card, dict):
                        evidence[card.get('id') or card.get('content_hash') or line] = card
                except ValueError:
                    pass
        except OSError:
            pass
    out['pipeline_details'] = {
        'plan': saved('agent_results/contrast_spec/research_plan.json', {}),
        'qualification': saved('agent_results/qualification/qualification_results.json', {}),
        'evidence': list(evidence.values()),
        'volcano': _volcano_data(out, root),
    }
    from .candidate_display import hidden_genes, symbol
    hidden = hidden_genes()
    hidden.update(symbol(point[0]) for point in (out['pipeline_details']['volcano'] or {}).get('points', [])
                  if point[3] == 'Down_in_TNBC')
    out['candidate_display'] = {'hidden_genes': sorted(hidden)}
    panels = _panels(summary.get('modalities', []), root)
    completed = [p for p in panels if not p['running'] and p['execution_status'] == 'COMPLETED']
    passed = [p for p in completed if p['validation'] == 'ADVANCE'
              and p['critic'] == 'ADVANCE' and not p['is_mock']]
    out['design_report'] = {
        'panels': panels, 'completed_routes': len(completed), 'passing_routes': len(passed),
        'candidate_count': sum(p['candidate_count'] for p in panels),
        'evaluated_targets': len(summary.get('candidates', [])),
        'candidate_search': summary.get('candidate_search', {}),
        'workflow_status': 'RUNNING' if data.get('running') else 'FINISHED',
        'report_status': 'LIVE' if data.get('running') else 'READY',
        'reference': None, 'followups': [],
    }
    reference_id = summary.get('design_reference_run')
    if reference_id and reference_id != data.get('run_id'):
        try:
            reference = read_report(root, run_id=reference_id, expand_integrated=False)
            if (reference.get('modality') != 'ALL' or reference.get('gene') != 'ERBB2'
                    or reference.get('execution_profile') != 'therapeutic_design' or reference.get('running')):
                raise ValueError('reference_must_be_a_finished_ERBB2_batch')
            out['design_report']['reference'] = {
                'run_id': reference_id, 'source_kind': 'PREVIOUSLY_RECORDED',
                'report_url': '/report?run_id=' + reference_id,
                'panels': _panels(reference.get('summary', {}).get('modalities', []), root),
            }
        except (OSError, ValueError, TypeError, KeyError) as exc:
            out['design_report']['reference_warning'] = str(exc)
    # Follow-up experiments belong to this parent, but never alter its historical totals.
    current_ids = {row.get('run_id') for row in summary.get('modalities', [])}
    if out.get('run_id'):
        for path in sorted(Path(root).glob('design_*/run_summary.json')):
            if path.parent.name in current_ids or not path.resolve().is_relative_to(Path(root).resolve()):
                continue
            try:
                child = json.loads(path.read_text())
                if child.get('parent_run_id') != out['run_id'] or child.get('modality') not in MODALITIES:
                    continue
                rows = _panels([{'run_id': path.parent.name, 'gene': child.get('gene'),
                                 'modality': child['modality']}], root)
                for item in rows:
                    item['reuse'] = child.get('reuse')
                out['design_report']['followups'].extend(rows)
            except (OSError, ValueError, TypeError, KeyError):
                continue
    return out
