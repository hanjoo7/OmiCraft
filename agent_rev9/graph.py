"""
OmiCraft LangGraph StateGraph — 에이전트 오케스트레이션

노드 구성:
  planner          LLM이 연구질문 파싱 → DAG/예산/contrast 결정
  r_analysis       DESeq2 + GSEA (R subprocess tool)
  apear            aPEAR 경로 네트워크 (tool)
  cell_context     CELLxGENE + genesorteR (tool)
  depmap           DepMap CRISPR/RNAi (tool)
  qualification    UniProt/HPA/OpenTargets/ChEMBL + Tier (tool + LLM)
  output_tables    TSV/JSONL 출력 테이블 (tool)
  critic           LLM falsification 6개 공격 → APPROVE/HOLD/REPLAN
  user_gate        사용자 선택 gate [interrupt_before]
  design           분자 설계 실행 (tool, parallel)

조건부 엣지:
  r_analysis  → apear        (success)
             → END           (failure)
  critic      → planner      (REPLAN, replan_count < 2)
             → user_gate     (APPROVE / HOLD)
  user_gate   → design       (ALLOWED 또는 all_eligible)
             → END           (PENDING / HOLD)

진입점: planner
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)

try:
    from langchain_core.runnables import RunnableConfig
    from langgraph.graph import StateGraph, END
    from langgraph.checkpoint.memory import MemorySaver
    _LANGGRAPH_AVAILABLE = True
except ImportError:
    RunnableConfig = dict
    _LANGGRAPH_AVAILABLE = False
    log.warning(
        "langgraph 패키지가 없습니다. "
        "pip install langgraph 후 재실행하세요. "
        "현재는 fallback PipelineGraph를 사용합니다."
    )

from .state import OmiCraftState
from .result_reuse import cached_node
from .pipeline_callbacks import invoke_node


def _runtime_node(name, function):
    def execute(state, config: RunnableConfig):
        return invoke_node(name, function, state, config)
    return execute


# ── 노드 함수 임포트 ───────────────────────────────────────

def _planner_node(state):
    from .planner_agent import planner_node
    return cached_node('planner', planner_node, state)


def _r_analysis_node(state):
    from .orchestrator import r_analysis_node
    return r_analysis_node(state)  # R has its own content-verified dataset cache.


def _apear_node(state):
    from .orchestrator import apear_node
    return cached_node('apear', apear_node, state)


def _cell_context_node(state):
    from .orchestrator import cell_context_node
    return cached_node('cell_context', cell_context_node, state)


def _depmap_node(state):
    from .orchestrator import depmap_node
    return cached_node('depmap', depmap_node, state)


def _discovery_llm_node(state):
    from .discovery_agent import discovery_llm_review
    return cached_node('discovery', discovery_llm_review, state)


def _qualification_node(state):
    from .qualification_agent import qualification_node
    return cached_node('qualification', qualification_node, state)


def _output_tables_node(state):
    from .orchestrator import output_tables_node
    return output_tables_node(state)


def _critic_node(state):
    from .critic_agent import critic_node
    return critic_node(state)


def _user_gate_node(state):
    from .orchestrator import user_selection_gate_node, assess_modalities
    assessed = assess_modalities(state)
    gated = user_selection_gate_node({**state, **assessed})
    return {**assessed, **gated}


def _dossier_node(state):
    from .dossier import dossier_node
    return dossier_node(state)


def _design_node(state):
    from .automatic_design import automatic_design_node
    return automatic_design_node(state)


# ── 조건부 라우팅 함수 ─────────────────────────────────────

def _route_after_r_analysis(state: dict) -> str:
    """R 분석 성공 여부에 따라 apear 진행 또는 파이프라인 종료."""
    if (state.get("r_analysis_result") or {}).get("success"):
        return "apear"
    return END


def _route_after_critic(state: dict) -> str:
    """Critic 판정에 따라 REPLAN(planner 회귀) 또는 user_gate 진행."""
    verdict = state.get("critic_verdict", "HOLD")
    replan_count = state.get("replan_count", 0)
    if verdict == "REPLAN" and replan_count < 2:
        log.info(
            "[Graph] Critic → REPLAN (횟수 %d/2) — planner로 회귀", replan_count
        )
        return "planner"
    return "user_gate"


def _route_after_gate(state: dict) -> str:
    """User selection gate 결과에 따라 design 실행 또는 종료."""
    if state.get("design_mode") == "all_eligible":
        return "design"
    assessments = state.get("modality_assessments") or []
    if any(a.get("gate", {}).get("allowed") for a in assessments):
        return "design"
    return "dossier_review"


# ── LangGraph 그래프 빌더 ──────────────────────────────────

def build_graph(checkpointer=None):
    """
    LangGraph StateGraph를 컴파일해 반환한다.

    langgraph가 설치되지 않은 경우 FallbackPipelineGraph를 반환한다.
    FallbackPipelineGraph는 동일한 .stream() 인터페이스를 제공하지만
    조건부 엣지와 human-in-the-loop는 지원하지 않는다.
    """
    if not _LANGGRAPH_AVAILABLE:
        log.warning("[Graph] langgraph 미설치 — FallbackPipelineGraph 사용")
        return _FallbackPipelineGraph()

    builder = StateGraph(OmiCraftState)

    # ── 노드 등록 ──────────────────────────────────────────
    builder.add_node("planner", _runtime_node("planner", _planner_node))
    builder.add_node("r_analysis", _runtime_node("r_analysis", _r_analysis_node))
    builder.add_node("apear", _runtime_node("apear", _apear_node))
    builder.add_node("cell_context", _runtime_node("cell_context", _cell_context_node))
    builder.add_node("depmap", _runtime_node("depmap", _depmap_node))
    builder.add_node("discovery", _runtime_node("discovery", _discovery_llm_node))
    builder.add_node("qualification", _runtime_node("qualification", _qualification_node))
    builder.add_node("output_tables", _runtime_node("output_tables", _output_tables_node))
    builder.add_node("critic", _runtime_node("critic", _critic_node))
    builder.add_node("user_gate", _runtime_node("user_gate", _user_gate_node))
    builder.add_node("design", _runtime_node("design", _design_node))
    builder.add_node("dossier_review", _runtime_node("dossier_review", _dossier_node))

    # ── 진입점 ────────────────────────────────────────────
    builder.set_entry_point("planner")

    # ── 순차 엣지 ─────────────────────────────────────────
    builder.add_edge("planner",       "r_analysis")
    builder.add_edge("apear",         "cell_context")
    builder.add_edge("cell_context",  "depmap")
    builder.add_edge("depmap",        "discovery")
    builder.add_edge("discovery",     "qualification")
    builder.add_edge("qualification", "output_tables")
    builder.add_edge("output_tables", "critic")
    builder.add_edge("design",        "dossier_review")
    builder.add_edge("dossier_review", END)

    # ── 조건부 엣지 ───────────────────────────────────────
    # R 분석 성공 → apear, 실패 → END
    builder.add_conditional_edges(
        "r_analysis",
        _route_after_r_analysis,
        {"apear": "apear", END: END},
    )

    # Critic 판정: REPLAN → planner, 나머지 → user_gate
    builder.add_conditional_edges(
        "critic",
        _route_after_critic,
        {"planner": "planner", "user_gate": "user_gate"},
    )

    # Gate 결과: ALLOWED → design, 아니면 → END
    builder.add_conditional_edges(
        "user_gate",
        _route_after_gate,
        {"design": "design", "dossier_review": "dossier_review"},
    )

    # ── 컴파일 (human-in-the-loop: user_gate 전 중단) ─────
    _checkpointer = checkpointer if checkpointer is not None else MemorySaver()
    graph = builder.compile(
        checkpointer=_checkpointer,
        interrupt_before=["user_gate"],
    )

    log.info("[Graph] LangGraph StateGraph 컴파일 완료")
    return graph


# ── Fallback: langgraph 미설치 환경용 ─────────────────────

class _FallbackPipelineGraph:
    """langgraph 미설치 환경에서 동일한 .stream() 인터페이스를 제공하는 대체 실행기.

    조건부 엣지와 REPLAN 루프는 지원하지 않는다.
    langgraph 설치 후 build_graph()를 재호출하면 완전 기능을 사용할 수 있다.
    """

    def stream(self, state: dict, config: dict = None):
        from .orchestrator import (
            r_analysis_node, apear_node, cell_context_node, depmap_node,
            output_tables_node, assess_modalities, user_selection_gate_node,
        )
        from .planner_agent import planner_node
        from .qualification_agent import qualification_node
        from .critic_agent import critic_node
        from .automatic_design import automatic_design_node
        from .discovery_agent import discovery_llm_review

        current = dict(state)
        nodes = [
            ("planner",       planner_node),
            ("r_analysis",    r_analysis_node),
            ("apear",         apear_node),
            ("cell_context",  cell_context_node),
            ("depmap",        depmap_node),
            ("discovery",     discovery_llm_review),
            ("qualification", qualification_node),
            ("output_tables", output_tables_node),
            ("critic",        critic_node),
            ("user_gate",     lambda s: {
                **assess_modalities(s),
                **user_selection_gate_node(s),
            }),
        ]

        if current.get("execution_profile") == "competition_demo" or current.get("advance_targets"):
            from .orchestrator import run_screening
            yield {"therapeutic_design": run_screening(current)}
            return

        if state.get("dry_run"):
            yield {"qualification": {"execution_status": "NOT_RUN",
                   "messages": ["Dry-run: upstream analysis not executed"]}}
            return

        if current.get("design_mode") == "all_eligible":
            nodes.append(("design", automatic_design_node))

        for node_name, fn in nodes:
            def execute(local_state):
                return cached_node(node_name, fn, local_state) if node_name in {
                    'planner', 'apear', 'cell_context', 'depmap', 'discovery', 'qualification'} else fn(local_state)
            update = invoke_node(node_name, execute, current, config)
            current.update(update)
            yield {node_name: update}

            # R 분석 실패 시 중단
            if node_name == "r_analysis" and not (update.get("r_analysis_result") or {}).get("success"):
                return

        # gate 통과 시 design 실행 (all_eligible 아닐 때)
        if current.get("design_mode") != "all_eligible":
            assessments = current.get("modality_assessments") or []
            if any(a.get("gate", {}).get("allowed") for a in assessments):
                from .orchestrator import run_screening
                update = run_screening(current)
                current.update(update)
                yield {"therapeutic_design": update}

        yield {"dossier_review": invoke_node("dossier_review", _dossier_node, current, config)}
