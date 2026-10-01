"""Keep input readiness separate from the outcome of ternary prediction."""


def ternary_execution_status(summary):
    prediction = summary.get('ternary_prediction', {})
    status = prediction.get('status', 'not_run')
    if status == 'interrupted':
        return 'INTERRUPTED'
    if status == 'success':
        return 'COMPLETED'
    if status in {'running', 'RUNNING'}:
        return 'RUNNING'
    if status not in {'not_run', 'NOT_RUN', None}:
        if prediction.get('phase_status') == 'TIMED_OUT' or prediction.get('error_code', '').endswith('_TIMEOUT'):
            return 'TIMED_OUT'
        if 'timed out' in prediction.get('error_message', '').lower():
            return 'TIMED_OUT'
        return 'PARTIAL' if status == 'partial' else 'FAILED'
    if summary.get('ternary_readiness', {}).get('blockers'):
        return 'BLOCKED'
    return 'NOT_RUN'
