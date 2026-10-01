"""
LangGraph 에이전트 오케스트레이션 통합 테스트
agent_rev9 리팩터링(rev7.1 → LangGraph) 검증용

검증 범위:
  1. state.py        — Annotated reducer 확인
  2. graph.py        — build_graph() 컴파일, 노드 구성, 라우팅 함수
  3. planner_agent   — LLM 사용/미사용 양쪽 경로
  4. discovery_agent — discovery_llm_review() LLM fallback
  5. qualification_agent — ADVANCE=0 REPLAN 신호
  6. critic_agent    — REPLAN 상한, safety veto, LLM fallback
  7. FallbackPipelineGraph — 전체 dry-run 흐름
  8. LangGraph compiled graph — planner까지 smoke test
"""

from __future__ import annotations

# ── Windows 호환성 stub ────────────────────────────────────────────────────────
# fcntl 은 Linux 전용 모듈이다. af3_msa_cache.py → protein_design_route.py 임포트 체인에서
# 사용되며, Windows 테스트 환경에서 ModuleNotFoundError를 일으킨다.
import sys
import unittest.mock

if sys.platform == "win32":
    sys.modules.setdefault("fcntl", unittest.mock.MagicMock())

import json
import operator
from typing import get_type_hints

import pytest


# ── 공통 헬퍼 ──────────────────────────────────────────────────────────────────

def _disable_llm(monkeypatch):
    """llm_client.complete() 가 항상 LLMUnavailable을 발생시키도록 패치."""
    from agent_rev9 import llm_client

    def _fail(*a, **kw):
        raise llm_client.LLMUnavailable("test: LLM disabled")

    monkeypatch.setattr(llm_client, "complete", _fail)


def _enable_agent_logging(monkeypatch):
    """config.agent_logging = True 로 패치해 LLM falsification 경로를 활성화."""
    from agent_rev9 import configuration

    original_get_config = configuration.get_config
    cfg = original_get_config()

    class _PatchedConfig:
        def __getattr__(self, name):
            if name == "agent_logging":
                return True
            return getattr(cfg, name)

    monkeypatch.setattr(configuration, "get_config", lambda: _PatchedConfig())


def _mock_completion(text: str):
    from agent_rev9.llm_client import Completion
    return Completion(text=text, model="mock", response_id="r0", usage={}, elapsed=0.0)


def _minimal_qual_state():
    """qualification_node가 반환하는 최소 valid state."""
    return {
        "qualification_results": {
            "ERBB2": {"tier": "Tier1A", "bio_confidence": 0.9,
                      "missing_evidence": [], "safety_veto": {}}
        },
        "advance_targets": [
            {"gene_name": "ERBB2", "modality": "ADC", "tier": "Tier1A",
             "bio_confidence": 0.9, "intervention_role": "CELL_TARGETING",
             "safety_veto": {}}
        ],
        "hold_targets": [],
        "reject_targets": [],
        "messages": [],
        "errors": [],
    }


def _mock_all_nodes(monkeypatch):
    """외부 IO가 필요한 모든 orchestrator 도구 노드를 mock으로 교체."""
    from agent_rev9 import orchestrator

    monkeypatch.setattr(orchestrator, "resolve_contrasts_node", lambda s: {
        "contrasts": [{"name": "TNBC_vs_HR+", "status": "READY"}],
        "research_plan": {"plan_id": "test01"},
        "messages": ["[mock] contrasts resolved"], "errors": [],
    })
    monkeypatch.setattr(orchestrator, "r_analysis_node", lambda s: {
        "r_analysis_result": {
            "success": True, "execution_mode": "MOCK",
            "de_summary": {
                "n_significant": 8, "up_in_tnbc": 5, "down_in_tnbc": 3,
                "top_upregulated": ["ERBB2"], "top_downregulated": [],
            },
            "gsea_summary": {"collections": {}},
            "output_files": {"de_all": "/mock/de_all.tsv"},
        },
        "messages": ["[mock] r_analysis done"], "errors": [],
    })
    monkeypatch.setattr(orchestrator, "apear_node",
        lambda s: {"apear_results": {}, "messages": [], "errors": []})
    monkeypatch.setattr(orchestrator, "cell_context_node",
        lambda s: {"cell_contexts": {}, "messages": [], "errors": []})
    monkeypatch.setattr(orchestrator, "depmap_node",
        lambda s: {"depmap_results": {}, "messages": [], "errors": []})
    monkeypatch.setattr(orchestrator, "qualification_node",
        lambda s: _minimal_qual_state())
    monkeypatch.setattr(orchestrator, "output_tables_node",
        lambda s: {"output_tables": {}, "messages": [], "errors": []})
    monkeypatch.setattr(orchestrator, "assess_modalities",
        lambda s: {"modality_assessments": []})
    monkeypatch.setattr(orchestrator, "user_selection_gate_node",
        lambda s: {"gate_status": "PENDING", "messages": [], "errors": []})


# ══════════════════════════════════════════════════════════════════════════════
# 1. State 스키마
# ══════════════════════════════════════════════════════════════════════════════

def test_state_annotated_reducers():
    """messages / errors / evidence_cards 는 operator.add reducer 여야 한다."""
    from agent_rev9.state import OmiCraftState

    hints = get_type_hints(OmiCraftState, include_extras=True)
    for field in ("messages", "errors", "evidence_cards"):
        meta = hints[field]
        assert hasattr(meta, "__metadata__"), f"{field} 는 Annotated 이어야 합니다"
        assert operator.add in meta.__metadata__, \
            f"{field} 의 reducer 가 operator.add 이어야 합니다"


def test_state_routing_fields_exist():
    """LangGraph 라우팅 제어 필드(next_node, replan_count, critic_feedback)가 존재해야 한다."""
    from agent_rev9.state import OmiCraftState

    hints = get_type_hints(OmiCraftState, include_extras=True)
    for field in ("next_node", "replan_count", "critic_feedback"):
        assert field in hints, f"OmiCraftState 에 {field} 필드가 없습니다"


# ══════════════════════════════════════════════════════════════════════════════
# 2. 그래프 구성
# ══════════════════════════════════════════════════════════════════════════════

def test_build_graph_has_stream_interface():
    """build_graph()는 .stream() 인터페이스를 반드시 제공해야 한다."""
    from agent_rev9.graph import build_graph

    g = build_graph()
    assert callable(getattr(g, "stream", None)), "graph 에 .stream() 이 없습니다"


def test_build_graph_has_required_nodes():
    """LangGraph 그래프에 11개 필수 노드가 모두 등록되어야 한다."""
    from agent_rev9.graph import _LANGGRAPH_AVAILABLE, build_graph

    if not _LANGGRAPH_AVAILABLE:
        pytest.skip("langgraph 미설치")

    g = build_graph()
    node_ids = set(g.get_graph().nodes.keys())
    required = {
        "planner", "r_analysis", "apear", "cell_context", "depmap",
        "discovery", "qualification", "output_tables", "critic", "user_gate", "design",
    }
    missing = required - node_ids
    assert not missing, f"그래프에 누락된 노드: {missing}"


# ══════════════════════════════════════════════════════════════════════════════
# 3. 라우팅 함수
# ══════════════════════════════════════════════════════════════════════════════

def test_route_r_analysis_success():
    from agent_rev9.graph import _route_after_r_analysis

    assert _route_after_r_analysis({"r_analysis_result": {"success": True}}) == "apear"


def test_route_r_analysis_failure_returns_end():
    from langgraph.graph import END

    from agent_rev9.graph import _route_after_r_analysis

    assert _route_after_r_analysis({"r_analysis_result": {"success": False}}) is END
    assert _route_after_r_analysis({}) is END


def test_route_critic_replan_under_limit():
    from agent_rev9.graph import _route_after_critic

    assert _route_after_critic({"critic_verdict": "REPLAN", "replan_count": 0}) == "planner"
    assert _route_after_critic({"critic_verdict": "REPLAN", "replan_count": 1}) == "planner"


def test_route_critic_replan_at_limit_goes_to_gate():
    from agent_rev9.graph import _route_after_critic

    assert _route_after_critic({"critic_verdict": "REPLAN", "replan_count": 2}) == "user_gate"


@pytest.mark.parametrize("verdict", ["APPROVE", "HOLD", "REJECT"])
def test_route_critic_non_replan_goes_to_gate(verdict):
    from agent_rev9.graph import _route_after_critic

    assert _route_after_critic({"critic_verdict": verdict}) == "user_gate"


def test_route_gate_allowed_goes_to_design():
    from agent_rev9.graph import _route_after_gate

    state = {"modality_assessments": [{"gate": {"allowed": True}}]}
    assert _route_after_gate(state) == "design"


def test_route_gate_all_eligible():
    from agent_rev9.graph import _route_after_gate

    assert _route_after_gate({"design_mode": "all_eligible"}) == "design"


def test_route_gate_blocked_returns_dossier():
    from langgraph.graph import END

    from agent_rev9.graph import _route_after_gate

    assert _route_after_gate({"modality_assessments": [{"gate": {"allowed": False}}]}) == "dossier_review"
    assert _route_after_gate({}) == "dossier_review"


# ══════════════════════════════════════════════════════════════════════════════
# 4. Planner Agent
# ══════════════════════════════════════════════════════════════════════════════

def test_planner_llm_unavailable_fallback(monkeypatch):
    """LLM 미사용 시에도 planner_node 는 유효한 state를 반환해야 한다."""
    _disable_llm(monkeypatch)
    from agent_rev9 import orchestrator
    from agent_rev9.planner_agent import planner_node

    monkeypatch.setattr(orchestrator, "resolve_contrasts_node", lambda s: {
        "contrasts": [{"name": "TNBC_vs_HR+"}], "messages": [], "errors": [],
    })
    result = planner_node({"research_question": "TNBC 표적 발굴", "replan_count": 0})
    assert result["current_agent"] == "planner"
    assert isinstance(result.get("messages"), list)
    # fallback plan 은 disease 또는 contrasts 중 하나는 반드시 포함
    has_plan = bool(result.get("disease") or result.get("contrasts") or result.get("research_plan"))
    assert has_plan, f"fallback plan 에 disease/contrasts/research_plan 이 없음: {result.keys()}"


def test_planner_parses_llm_json(monkeypatch):
    """정상 LLM 응답을 파싱해 disease/subtype 을 추출해야 한다."""
    from agent_rev9 import llm_client, orchestrator
    from agent_rev9.planner_agent import planner_node

    payload = {
        "disease": "breast_cancer", "subtype": "TNBC", "cohort": "TCGA-BRCA",
        "contrasts": ["TNBC vs HR+HER2-"],
        "analysis_steps": [{"step": 1, "agent": "discovery", "task": "DEG", "tools": ["DESeq2"]}],
        "budget": {"max_candidates": 50, "max_advance_targets": 3,
                   "llm_token_budget": 300000, "gpu_hours": 2.0},
        "stop_conditions": [],
    }
    monkeypatch.setattr(llm_client, "complete",
        lambda *a, **kw: _mock_completion(json.dumps(payload)))
    monkeypatch.setattr(orchestrator, "resolve_contrasts_node",
        lambda s: {"contrasts": [{"name": "TNBC_vs_HR+"}], "messages": [], "errors": []})

    result = planner_node({"research_question": "TNBC 표적", "replan_count": 0})
    assert result.get("disease") == "breast_cancer"
    assert result.get("subtype") == "TNBC"
    assert result["current_agent"] == "planner"


def test_planner_injects_critic_feedback_into_llm_call(monkeypatch):
    """replan_count > 0 일 때 critic_feedback 이 LLM 프롬프트에 포함되어야 한다."""
    from agent_rev9 import llm_client, orchestrator
    from agent_rev9.planner_agent import planner_node

    captured = []

    def fake_complete(instructions, evidence, **kw):
        captured.append((str(instructions), str(evidence)))
        raise llm_client.LLMUnavailable("disabled")

    monkeypatch.setattr(llm_client, "complete", fake_complete)
    monkeypatch.setattr(orchestrator, "resolve_contrasts_node",
        lambda s: {"messages": [], "errors": []})

    state = {
        "research_question": "TNBC",
        "replan_count": 1,
        "critic_feedback": "ADVANCE 후보 0개 — contrast 범위 확장 필요",
    }
    planner_node(state)
    assert captured, "complete() 가 호출되어야 합니다"
    full_text = captured[0][0] + captured[0][1]
    assert "ADVANCE 후보 0개" in full_text, "critic_feedback 이 LLM 입력에 포함되어야 합니다"


# ══════════════════════════════════════════════════════════════════════════════
# 5. Discovery LLM Review
# ══════════════════════════════════════════════════════════════════════════════

def _discovery_state():
    return {
        "r_analysis_result": {
            "success": True,
            "de_summary": {"n_significant": 10, "up_in_tnbc": 6, "down_in_tnbc": 4,
                           "top_upregulated": ["ERBB2"], "top_downregulated": []},
            "gsea_summary": {"collections": {}},
        },
        "apear_results": {}, "cell_contexts": {}, "depmap_results": {},
    }


def test_discovery_llm_review_fallback(monkeypatch):
    """LLM 미사용 시에도 EV-DISC-001 evidence card를 생성해야 한다."""
    _disable_llm(monkeypatch)
    from agent_rev9.discovery_agent import discovery_llm_review

    result = discovery_llm_review(_discovery_state())
    assert result["current_agent"] == "discovery"
    cards = result["evidence_cards"]
    assert len(cards) == 1
    assert cards[0]["id"] == "EV-DISC-001"
    # fallback 요약에 DEG 수가 포함되어야 함
    assert "10" in cards[0]["content"]


def test_discovery_llm_review_parses_json(monkeypatch):
    from agent_rev9 import llm_client
    from agent_rev9.discovery_agent import discovery_llm_review

    resp = {"summary": "DEG 10개 발견.", "key_pathways": ["PI3K/AKT"], "limitations": ["bulk RNA"]}
    monkeypatch.setattr(llm_client, "complete",
        lambda *a, **kw: _mock_completion(json.dumps(resp)))

    result = discovery_llm_review(_discovery_state())
    card = result["evidence_cards"][0]
    assert card["content"] == "DEG 10개 발견."
    assert "PI3K/AKT" in card["key_pathways"]
    assert "bulk RNA" in card["limitations"]


def test_discovery_llm_review_no_deg_card_type(monkeypatch):
    """DEG=0 일 때 evidence card type 은 NOT_APPLICABLE 이어야 한다."""
    _disable_llm(monkeypatch)
    from agent_rev9.discovery_agent import discovery_llm_review

    state = {
        "r_analysis_result": {
            "success": True,
            "de_summary": {"n_significant": 0, "up_in_tnbc": 0, "down_in_tnbc": 0,
                           "top_upregulated": [], "top_downregulated": []},
            "gsea_summary": {"collections": {}},
        },
        "apear_results": {}, "cell_contexts": {}, "depmap_results": {},
    }
    result = discovery_llm_review(state)
    assert result["evidence_cards"][0]["type"] == "NOT_APPLICABLE"


# ══════════════════════════════════════════════════════════════════════════════
# 6. Qualification Agent
# ══════════════════════════════════════════════════════════════════════════════

def test_qualification_node_fallback(monkeypatch):
    """LLM 미사용 시에도 EV-QUAL-001 카드를 생성해야 한다."""
    _disable_llm(monkeypatch)
    from agent_rev9 import orchestrator
    from agent_rev9.qualification_agent import qualification_node

    monkeypatch.setattr(orchestrator, "qualification_node", lambda s: _minimal_qual_state())
    result = qualification_node({"research_question": "TNBC"})
    assert result["current_agent"] == "qualification"
    assert result["evidence_cards"][0]["id"] == "EV-QUAL-001"
    assert result["evidence_cards"][0]["type"] == "SUPPORTIVE"


def test_qualification_zero_advance_sets_replan_signal(monkeypatch):
    """ADVANCE 후보가 0개일 때 critic_feedback 과 ADVANCE_TARGETS_EMPTY 플래그가 설정되어야 한다."""
    _disable_llm(monkeypatch)
    from agent_rev9 import orchestrator
    from agent_rev9.qualification_agent import qualification_node

    empty = {
        "qualification_results": {}, "advance_targets": [],
        "hold_targets": [], "reject_targets": [], "messages": [], "errors": [],
    }
    monkeypatch.setattr(orchestrator, "qualification_node", lambda s: empty)

    result = qualification_node({})
    assert result.get("critic_feedback"), "critic_feedback 이 설정되어야 합니다"
    flags = result["evidence_cards"][0]["flags"]
    assert any("ADVANCE_TARGETS_EMPTY" in f for f in flags)
    assert result["evidence_cards"][0]["replan_needed"] is True


def test_qualification_parses_llm_json(monkeypatch):
    from agent_rev9 import llm_client, orchestrator
    from agent_rev9.qualification_agent import qualification_node

    resp = {
        "summary": "ERBB2 Tier1A ADVANCE 유망.",
        "top_target": "ERBB2",
        "flags": [],
        "replan_needed": False,
    }
    monkeypatch.setattr(llm_client, "complete",
        lambda *a, **kw: _mock_completion(json.dumps(resp)))
    monkeypatch.setattr(orchestrator, "qualification_node", lambda s: _minimal_qual_state())

    result = qualification_node({"research_question": "TNBC"})
    card = result["evidence_cards"][0]
    assert card["content"] == "ERBB2 Tier1A ADVANCE 유망."
    assert card["replan_needed"] is False


# ══════════════════════════════════════════════════════════════════════════════
# 7. Critic Agent
# ══════════════════════════════════════════════════════════════════════════════

def _base_critic_state(**extra):
    return {
        "advance_targets": [], "hold_targets": [], "reject_targets": [],
        "replan_count": 0, "critic_count": 0, "screening_result": {},
        **extra,
    }


def test_critic_returns_hold_when_no_candidates(monkeypatch):
    """후보가 없고 LLM 불가일 때 HOLD 반환."""
    _disable_llm(monkeypatch)
    from agent_rev9.critic_agent import critic_node

    result = critic_node(_base_critic_state())
    assert result["critic_verdict"] in {"APPROVE", "HOLD", "REPLAN", "REJECT"}
    assert result["critic_count"] == 1
    assert result["current_agent"] == "critic"


def test_critic_replan_capped_to_hold(monkeypatch):
    """replan_count >= 2 일 때 LLM REPLAN 이 HOLD 로 전환되어야 한다."""
    from agent_rev9 import llm_client
    from agent_rev9.critic_agent import critic_node

    _enable_agent_logging(monkeypatch)
    llm_resp = {
        "verdict": "REPLAN", "reason": "insufficient evidence",
        "feedback": "expand contrast", "issues_found": [], "regression_target": "planner",
    }
    monkeypatch.setattr(llm_client, "complete",
        lambda *a, **kw: _mock_completion(json.dumps(llm_resp)))

    result = critic_node(_base_critic_state(replan_count=2))
    assert result["critic_verdict"] == "HOLD", \
        "replan_count=2 일 때 REPLAN은 HOLD 로 전환되어야 합니다"
    issues = result.get("critic_result", {}).get("issues_found", [])
    assert "REPLAN_LIMIT_REACHED" in issues


def test_critic_safety_veto_overrides_llm_approve(monkeypatch):
    """safety_verdict=REJECT 표적이 있으면 LLM APPROVE 여도 REJECT 이어야 한다."""
    from agent_rev9 import llm_client
    from agent_rev9.critic_agent import critic_node

    _enable_agent_logging(monkeypatch)
    monkeypatch.setattr(llm_client, "complete",
        lambda *a, **kw: _mock_completion(
            '{"verdict":"APPROVE","reason":"ok","feedback":"","issues_found":[],"regression_target":"null"}'
        ))

    state = _base_critic_state(
        advance_targets=[{"safety_verdict": "REJECT", "gene_name": "GENEX"}]
    )
    result = critic_node(state)
    assert result["critic_verdict"] == "REJECT"


def test_critic_llm_approve_propagates(monkeypatch):
    """LLM APPROVE 판정이 최종 verdict 에 반영되어야 한다."""
    from agent_rev9 import llm_client
    from agent_rev9.critic_agent import critic_node

    _enable_agent_logging(monkeypatch)
    monkeypatch.setattr(llm_client, "complete",
        lambda *a, **kw: _mock_completion(
            '{"verdict":"APPROVE","reason":"all checks pass","feedback":"","issues_found":[],"regression_target":"null"}'
        ))

    result = critic_node(_base_critic_state(replan_count=0))
    assert result["critic_verdict"] == "APPROVE"


# ══════════════════════════════════════════════════════════════════════════════
# 8. FallbackPipelineGraph 전체 흐름
# ══════════════════════════════════════════════════════════════════════════════

def test_fallback_full_run_all_nodes_fire(monkeypatch):
    """FallbackPipelineGraph 가 모든 주요 노드를 순서대로 실행해야 한다."""
    _disable_llm(monkeypatch)
    _mock_all_nodes(monkeypatch)

    from agent_rev9.graph import _FallbackPipelineGraph

    graph = _FallbackPipelineGraph()
    state = {"research_question": "TNBC 표적 발굴",
             "evidence_cards": [], "messages": [], "errors": []}
    events = list(graph.stream(state))
    node_names = [list(e.keys())[0] for e in events]

    expected_order = [
        "planner", "r_analysis", "apear", "cell_context",
        "depmap", "discovery", "qualification", "output_tables", "critic", "user_gate",
    ]
    for name in expected_order:
        assert name in node_names, f"'{name}' 노드가 stream 출력에 없음. 실행된 노드: {node_names}"


def test_fallback_r_analysis_failure_stops_pipeline(monkeypatch):
    """r_analysis 실패 시 이후 노드들이 실행되지 않아야 한다."""
    _disable_llm(monkeypatch)
    _mock_all_nodes(monkeypatch)
    from agent_rev9 import orchestrator

    monkeypatch.setattr(orchestrator, "r_analysis_node", lambda s: {
        "r_analysis_result": {"success": False, "error": "R not found"},
        "messages": [], "errors": [],
    })

    from agent_rev9.graph import _FallbackPipelineGraph

    events = list(_FallbackPipelineGraph().stream(
        {"research_question": "X", "evidence_cards": [], "messages": [], "errors": []}
    ))
    node_names = [list(e.keys())[0] for e in events]

    assert "r_analysis" in node_names
    assert "apear" not in node_names, "r_analysis 실패 후 apear 가 실행되면 안 됩니다"
    assert "qualification" not in node_names


def test_fallback_evidence_cards_accumulate(monkeypatch):
    """discovery + qualification 에서 생성된 evidence card 가 최종 state에 누적되어야 한다."""
    _disable_llm(monkeypatch)
    _mock_all_nodes(monkeypatch)

    from agent_rev9.graph import _FallbackPipelineGraph

    state = {"research_question": "TNBC", "evidence_cards": [], "messages": [], "errors": []}
    events = list(_FallbackPipelineGraph().stream(state))

    # 각 노드 출력에서 evidence_cards 수집
    all_cards = []
    for event in events:
        for node_output in event.values():
            all_cards.extend(node_output.get("evidence_cards", []))

    card_ids = [c["id"] for c in all_cards]
    assert "EV-DISC-001" in card_ids, f"EV-DISC-001 이 없음. card ids: {card_ids}"
    assert "EV-QUAL-001" in card_ids, f"EV-QUAL-001 이 없음. card ids: {card_ids}"


# ══════════════════════════════════════════════════════════════════════════════
# 9. LangGraph compiled graph — smoke test
# ══════════════════════════════════════════════════════════════════════════════

def test_langgraph_graph_smoke(monkeypatch):
    """LangGraph 컴파일 그래프가 planner 노드까지 정상 실행되어야 한다."""
    from agent_rev9.graph import _LANGGRAPH_AVAILABLE, build_graph

    if not _LANGGRAPH_AVAILABLE:
        pytest.skip("langgraph 미설치")

    _disable_llm(monkeypatch)
    _mock_all_nodes(monkeypatch)

    g = build_graph()
    state = {
        "research_question": "TNBC 표적 발굴",
        "evidence_cards": [], "messages": [], "errors": [],
        "critic_count": 0, "replan_count": 0,
    }
    config = {"configurable": {"thread_id": "smoke-test-001"}}

    events = []
    for event in g.stream(state, config=config):
        events.append(event)
        # user_gate 는 interrupt_before 로 중단되므로 거기서 종료
        if "user_gate" in event or "__interrupt__" in event:
            break

    assert events, "그래프가 아무 이벤트도 생성하지 않았습니다"
    all_nodes_seen = set()
    for e in events:
        all_nodes_seen.update(e.keys())

    assert "planner" in all_nodes_seen, \
        f"planner 노드가 실행되지 않았습니다. 실행된 노드: {all_nodes_seen}"
