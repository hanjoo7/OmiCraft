"""Assemble the final evidence dossier without changing scientific decisions."""
from datetime import datetime, timezone
from pathlib import Path

from .configuration import get_config
from .small_molecule_io import write_json


def dossier_node(state):
    directory = Path(state['_run_directory']) / 'agent_results' if state.get('_run_directory') else Path(get_config().data.agent_results)
    path = directory / 'final_dossier.json'
    routes = state.get('automatic_design_results', [])
    passed = [row for row in routes if row.get('execution_status') == 'COMPLETED'
              and row.get('validation_status') == 'ADVANCE'
              and row.get('critic_decision') == 'ADVANCE' and not row.get('is_mock')]
    data = {
        'run_id': state.get('run_id'), 'generated_at': datetime.now(timezone.utc).isoformat(),
        'research_question': state.get('research_question'),
        'critic_verdict': state.get('critic_verdict', 'NOT_EVALUATED'),
        'critic_reason': state.get('critic_reason'), 'review_issues': state.get('review_issues', []),
        'qualification_results': state.get('qualification_results', {}),
        'advance_targets': state.get('advance_targets', []),
        'hold_targets': state.get('hold_targets', []), 'reject_targets': state.get('reject_targets', []),
        'design_results': routes, 'validated_designs': passed,
        'design_outcome': state.get('design_outcome', {}),
        'evidence_cards': state.get('evidence_cards', []),
        'errors': state.get('errors', []),
        'scope': 'Recorded evidence and computational results; experimental validation remains required',
    }
    write_json(path, data)
    return {'dossier_path': str(path), 'dossier_status': 'COMPLETED',
            'messages': [f'[Dossier] 근거 및 설계 {len(routes)}개 경로 정리 완료 · 계산 검증 통과 {len(passed)}개 · 전체 판정 {data["critic_verdict"]}']}
