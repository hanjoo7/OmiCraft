"""Local report index. The public proxy does not expose this endpoint."""
import json
import re
from pathlib import Path

RUN_NAME = re.compile(r'(?:upstream|design|lab)_\d{8}T\d{6}Z_[0-9a-f]{8}')


def report_history(root, saved_run, limit=100):
    root = Path(root).resolve()
    rows = []
    if not root.is_dir():
        return {"runs": rows}
    for path in root.iterdir():
        if not (RUN_NAME.fullmatch(path.name) or path.name == saved_run):
            continue
        summary = path / 'run_summary.json'
        if not summary.is_file() or not summary.resolve().is_relative_to(root):
            continue
        try:
            data = json.loads(summary.read_text())
            if data.get('execution_profile') not in {'upstream_analysis', 'therapeutic_design', 'structure_lab'} and path.name != saved_run:
                continue
            row = {key: data.get(key) for key in ('gene', 'modality', 'execution_profile',
                    'execution_status', 'validation_decision', 'parent_run_id', 'running')}
            row.update(run_id=path.name, report_url='/report?run_id=' + path.name)
            rows.append(row)
        except (OSError, ValueError, TypeError, AttributeError):
            continue
    rows.sort(key=lambda row: ((re.search(r'\d{8}T\d{6}Z', row['run_id']) or [''])[0], row['run_id']), reverse=True)
    return {'runs': rows[:limit]}
