"""
Agent ② Omics Discovery
DEG / Pathway / Cell-context / DepMap 도구 실행 후 LLM이 결과를 해석한다.

tool 담당 (결정론적):
  r_analysis_node    DESeq2 + GSEA
  apear_node         aPEAR 경로 네트워크
  cell_context_node  CELLxGENE + genesorteR
  depmap_node        DepMap CRISPR/RNAi 의존성

LLM 담당:
  - DEG 방향 / pathway / 세포맥락 결과 해석 → evidence card 생성
  - R 분석 실패 시 critic_feedback 생성 → graph가 REPLAN 트리거 가능

출력 state 필드:
  r_analysis_result, apear_results, cell_context_results, depmap_results,
  evidence_cards, messages, errors, (next_node on failure)
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)

_DISCOVERY_SYSTEM = """당신은 OmiCraft의 Discovery 에이전트입니다.
제공된 오믹스 분석 결과만을 근거로 한국어 2~4문장으로 설명하세요.

규칙:
- 계산 결과와 기존 판정을 변경하지 마세요
- 누락된 근거를 통과로 표현하지 마세요
- 새 실험, 검증 통과, 결합 효능을 추정해서 주장하지 마세요
- DEG 방향, 경로 및 세포맥락의 의미와 한계를 객관적으로 요약하세요

반드시 아래 JSON으로만 응답하세요:
{
  "summary": "2~4문장 해석",
  "key_pathways": ["pathway 1", "pathway 2"],
  "limitations": ["한계 1", "한계 2"]
}"""


def _build_evidence_summary(r_result: dict, apear_result: dict, cell_result: dict, depmap_result: dict) -> dict:
    """LLM에 전달할 증거 요약 구성."""
    de_summary = r_result.get("de_summary", {})
    gsea_summary = r_result.get("gsea_summary", {})
    apear_nodes = sum(
        v.get("n_nodes", 0) for v in apear_result.values()
        if isinstance(v, dict)
    ) if apear_result else 0

    cell_counts = {}
    for gene, ctx in (cell_result or {}).items():
        cls = ctx.get("top_class") or ctx.get("context_status") or "UNKNOWN"
        cell_counts[cls] = cell_counts.get(cls, 0) + 1

    dep_concordant = sum(
        1 for v in (depmap_result or {}).values()
        if v.get("dependency_status") == "CONCORDANT"
    )

    return {
        "n_deg_significant": de_summary.get("n_significant", 0),
        "n_up_tnbc": de_summary.get("up_in_tnbc", 0),
        "n_down_tnbc": de_summary.get("down_in_tnbc", 0),
        "top_upregulated": de_summary.get("top_upregulated", [])[:5],
        "top_downregulated": de_summary.get("top_downregulated", [])[:5],
        "gsea_significant": sum(
            v.get("significant", 0) for v in (gsea_summary.get("collections") or {}).values()
        ),
        "apear_pathway_nodes": apear_nodes,
        "cell_context_classes": cell_counts,
        "depmap_concordant": dep_concordant,
        "execution_mode": r_result.get("execution_mode", "UNKNOWN"),
    }


def discovery_llm_review(state: dict) -> dict:
    """Discovery LLM 검토 노드: 4개 tool 노드의 결과를 읽어 LLM이 해석 + evidence card 생성.

    graph.py에서 depmap 노드 직후, qualification 노드 직전에 실행된다.
    tool 노드(r_analysis/apear/cell_context/depmap)는 이미 state에 결과를 기록했으므로
    이 함수는 tool을 다시 실행하지 않고 LLM 해석과 evidence card 생성만 담당한다.
    """
    from .llm_client import LLMUnavailable, complete
    import json, re

    r_result = state.get("r_analysis_result") or {}
    apear_results = state.get("apear_results") or {}
    cell_contexts = state.get("cell_contexts") or state.get("cell_context_results") or {}
    depmap_results = state.get("depmap_results") or {}

    evidence_summary = _build_evidence_summary(r_result, apear_results, cell_contexts, depmap_results)
    llm_summary = ""
    key_pathways: list[str] = []
    limitations: list[str] = []
    messages: list[str] = []

    try:
        response = complete(
            instructions=_DISCOVERY_SYSTEM,
            evidence=evidence_summary,
            role="discovery",
        )
        cleaned = re.sub(r"```json\s*|```\s*", "", response.text).strip()
        try:
            parsed = json.loads(cleaned)
            llm_summary = parsed.get("summary", response.text[:300])
            key_pathways = parsed.get("key_pathways", [])
            limitations = parsed.get("limitations", [])
        except (json.JSONDecodeError, AttributeError):
            llm_summary = response.text[:300]
        log.info("[Discovery] LLM 해석 완료 (%.1fs)", response.elapsed)
        messages.append(f"[Discovery] LLM 해석: {llm_summary[:80]}…")
    except LLMUnavailable as exc:
        log.warning("[Discovery] LLM 불가 — 해석 생략: %s", exc)
        llm_summary = (
            f"DEG {evidence_summary['n_deg_significant']}개 "
            f"(TNBC 상향 {evidence_summary['n_up_tnbc']}, "
            f"하향 {evidence_summary['n_down_tnbc']}). "
            f"GSEA 유의 {evidence_summary['gsea_significant']}개."
        )

    evidence_cards = [{
        "id": "EV-DISC-001",
        "source": "Discovery",
        "type": "SUPPORTIVE" if evidence_summary["n_deg_significant"] > 0 else "NOT_APPLICABLE",
        "content": llm_summary,
        "gene": "",
        "confidence": 0.8 if r_result.get("success") else 0.3,
        "key_pathways": key_pathways,
        "limitations": limitations,
    }]

    messages.append(
        f"[Discovery] 완료 — DEG={evidence_summary['n_deg_significant']}, "
        f"DepMap concordant={evidence_summary['depmap_concordant']}"
    )

    return {
        "evidence_cards": evidence_cards,
        "messages": messages,
        "errors": [],
        "current_agent": "discovery",
    }


def discovery_node(state: dict) -> dict:
    """Discovery 노드: tool 실행 → LLM 해석 → evidence card 생성."""
    from .orchestrator import r_analysis_node, apear_node, cell_context_node, depmap_node
    from .llm_client import LLMUnavailable, complete
    import json, re

    current = dict(state)
    output: dict = {}
    messages: list[str] = []
    errors: list[str] = []

    # ── Tool 실행 ──────────────────────────────────────────
    for node_fn in (r_analysis_node, apear_node, cell_context_node, depmap_node):
        result = node_fn(current)
        current.update(result)
        output.update(result)
        messages.extend(result.pop("messages", []))
        errors.extend(result.pop("errors", []))

        # R 분석 실패 시 중단 — critic에 feedback 제공
        if node_fn is r_analysis_node and not (result.get("r_analysis_result") or {}).get("success"):
            reason = (result.get("r_analysis_result") or {}).get("error", "R 분석 실패")
            log.warning("[Discovery] R 분석 실패: %s", reason)
            return {
                **output,
                "messages": messages + [f"[Discovery] R 분석 실패 → critic 피드백 생성"],
                "errors": errors + [f"[Discovery] {reason}"],
                "critic_feedback": f"R 분석 실패: {reason}. 데이터 경로와 입력 파일을 점검하세요.",
                "current_agent": "discovery",
            }

    # ── LLM 해석 ──────────────────────────────────────────
    r_result = output.get("r_analysis_result") or {}
    apear_results = output.get("apear_results") or {}
    cell_contexts = output.get("cell_contexts") or output.get("cell_context_results") or {}
    depmap_results = output.get("depmap_results") or {}

    evidence_summary = _build_evidence_summary(r_result, apear_results, cell_contexts, depmap_results)
    llm_summary = ""
    key_pathways: list[str] = []
    limitations: list[str] = []

    try:
        response = complete(
            instructions=_DISCOVERY_SYSTEM,
            evidence=evidence_summary,
            role="discovery",
        )
        # JSON 파싱 시도
        cleaned = re.sub(r"```json\s*|```\s*", "", response.text).strip()
        try:
            parsed = json.loads(cleaned)
            llm_summary = parsed.get("summary", response.text[:300])
            key_pathways = parsed.get("key_pathways", [])
            limitations = parsed.get("limitations", [])
        except (json.JSONDecodeError, AttributeError):
            llm_summary = response.text[:300]
        log.info("[Discovery] LLM 해석 완료 (%.1fs)", response.elapsed)
        messages.append(f"[Discovery] LLM 해석: {llm_summary[:80]}…")
    except LLMUnavailable as exc:
        log.warning("[Discovery] LLM 불가 — 해석 생략: %s", exc)
        llm_summary = (
            f"DEG {evidence_summary['n_deg_significant']}개 "
            f"(TNBC 상향 {evidence_summary['n_up_tnbc']}, "
            f"하향 {evidence_summary['n_down_tnbc']}). "
            f"GSEA 유의 {evidence_summary['gsea_significant']}개."
        )

    # ── Evidence card 생성 ────────────────────────────────
    evidence_cards = [{
        "id": "EV-DISC-001",
        "source": "Discovery",
        "type": "SUPPORTIVE" if evidence_summary["n_deg_significant"] > 0 else "NOT_APPLICABLE",
        "content": llm_summary,
        "gene": "",
        "confidence": 0.8 if r_result.get("success") else 0.3,
        "key_pathways": key_pathways,
        "limitations": limitations,
    }]

    messages.append(
        f"[Discovery] 완료 — DEG={evidence_summary['n_deg_significant']}, "
        f"DepMap concordant={evidence_summary['depmap_concordant']}"
    )

    return {
        **output,
        "evidence_cards": evidence_cards,
        "messages": messages,
        "errors": errors,
        "current_agent": "discovery",
    }
