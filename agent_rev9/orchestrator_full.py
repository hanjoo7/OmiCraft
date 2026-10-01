"""Compatibility entry point for the v7 workflow."""
from .agent_narration import critic_review as omics_critic_node, dossier_review as dossier_node
from .user_selection import user_selection_gate_node as user_selection_node


def run_full_pipeline(question=None, verbose=True, **options):
    from .orchestrator import run_full_pipeline as run
    state = dict(question) if isinstance(question, dict) else {'research_question': question or 'TNBC 치료 표적 및 모달리티 평가'}
    output = run({**state, **options})
    if verbose:
        for message in output.get('messages', []):
            print(message)
    return output


if __name__ == '__main__':
    run_full_pipeline()
