"""
Agent ③ Qualification
UniProt/HPA/OpenTargets/ChEMBL + Tier 결정 후 LLM이 결과를 검토한다.

tool 담당 (결정론적):
  orchestrator.qualification_node  — run_qualification_pipeline() 호출
  orchestrator.output_tables_node  — TSV/JSONL 출력 테이블 생성 (qualification_node 내부에서 처리)

LLM 담당:
  - ADVANCE 후보 목록과 Tier 판정 검토 → 치료 전략적 의미 요약
  - 누락 근거 또는 의심스러운 판정에 대한 주의 플래그 → evidence card
  - critic 에게 전달할 qualification 요약 생성

출력 state 필드:
  qualification_results, advance_targets, hold_targets, reject_targets,
  evidence_cards, messages, errors
"""

from __future__ import annotations

import json
import logging
import re

log = logging.getLogger(__name__)

_QUALIFICATION_SYSTEM = """당신은 OmiCraft의 Qualification 에이전트입니다.
제공된 표적 적격성 평가 결과만을 근거로 한국어 2~4문장으로 설명하세요.

규칙:
- 기존 Tier 판정, 안전성 veto, ChEMBL phase를 변경하지 마세요
- 누락된 근거를 통과로 표현하지 마세요
- 표적이 0개인 경우 분석 설계 변경을 제안하세요

반드시 아래 JSON으로만 응답하세요:
{
  "summary": "2~4문장 핵심 인사이트",
  "top_target": "가장 유망한 표적 유전자명 또는 null",
  "flags": ["주의 사항 목록 — 의심스러운 판정, 누락 근거 등"],
  "replan_needed": true | false
}"""


def _build_qualification_summary(qr: dict, advance: list, hold: list, reject: list) -> dict:
    """LLM에 전달할 qualification 요약."""
    advance_brief = [
        {k: t.get(k) for k in ("gene_name", "modality", "tier", "bio_confidence",
                                 "intervention_role", "safety_veto")}
        for t in advance[:8]
    ]
    missing_evidence_genes = [
        gene for gene, rec in qr.items()
        if rec.get("missing_evidence")
    ]
    veto_genes = [
        gene for gene, rec in qr.items()
        if any((v or {}).get("has_veto") for v in (rec.get("safety_veto") or {}).values()
               if isinstance(v, dict))
    ]
    return {
        "n_advance": len(advance),
        "n_hold": len(hold),
        "n_reject": len(reject),
        "advance_targets": advance_brief,
        "missing_evidence_genes": missing_evidence_genes[:10],
        "safety_veto_genes": veto_genes[:10],
        "all_reject": len(advance) == 0,
    }


def qualification_node(state: dict) -> dict:
    """Qualification 노드: tool 실행 → LLM 검토 → evidence card 생성."""
    from .orchestrator import qualification_node as _run_qualification
    from .llm_client import LLMUnavailable, complete

    # ── Tool 실행 ──────────────────────────────────────────
    result = _run_qualification(state)
    messages: list[str] = list(result.pop("messages", []))
    errors: list[str] = list(result.pop("errors", []))

    qr = result.get("qualification_results") or {}
    advance = result.get("advance_targets") or []
    hold = result.get("hold_targets") or []
    reject = result.get("reject_targets") or []

    # ── LLM 검토 ──────────────────────────────────────────
    qual_summary = _build_qualification_summary(qr, advance, hold, reject)
    llm_summary = ""
    flags: list[str] = []
    replan_needed = False

    try:
        response = complete(
            instructions=_QUALIFICATION_SYSTEM,
            evidence=qual_summary,
            role="qualification",
        )
        cleaned = re.sub(r"```json\s*|```\s*", "", response.text).strip()
        try:
            parsed = json.loads(cleaned)
            llm_summary = parsed.get("summary", response.text[:300])
            flags = parsed.get("flags", [])
            replan_needed = bool(parsed.get("replan_needed", False))
        except (json.JSONDecodeError, AttributeError):
            llm_summary = response.text[:300]
        log.info("[Qualification] LLM 검토 완료 (%.1fs)", response.elapsed)
        messages.append(f"[Qualification] LLM 검토: {llm_summary[:80]}…")
    except LLMUnavailable as exc:
        log.warning("[Qualification] LLM 불가 — 기본 요약: %s", exc)
        llm_summary = (
            f"ADVANCE {len(advance)}개, HOLD {len(hold)}개, "
            f"REJECT {len(reject)}개. "
            + ("ADVANCE 후보 없음 — 파이프라인 재설계 검토 필요." if not advance else "")
        )
        if not advance:
            replan_needed = True

    # ADVANCE 0건 플래그
    if not advance:
        flags.append("ADVANCE_TARGETS_EMPTY: Planner 재계획 또는 contrast 범위 확장 필요")
        replan_needed = True

    # ── Evidence card 생성 ────────────────────────────────
    ev_type = "SUPPORTIVE" if advance else "CONTRADICTORY"
    evidence_cards = [{
        "id": "EV-QUAL-001",
        "source": "Qualification",
        "type": ev_type,
        "content": llm_summary,
        "gene": (advance[0].get("gene_name") if advance else ""),
        "confidence": 0.9 if len(advance) >= 1 else 0.4,
        "flags": flags,
        "replan_needed": replan_needed,
    }]

    # REPLAN 신호: LLM이 재계획 필요하다고 판단하면 critic_feedback 사전 설정
    extra: dict = {}
    if replan_needed:
        extra["critic_feedback"] = (
            f"Qualification 결과 ADVANCE 후보 {len(advance)}개. "
            + ("; ".join(flags) if flags else "")
            + " — contrast 범위를 확장하거나 분석 전략을 재검토하세요."
        )
        messages.append("[Qualification] REPLAN 신호 — critic으로 피드백 전달")

    messages.append(
        f"[Qualification] 완료 — ADVANCE={len(advance)}, "
        f"HOLD={len(hold)}, REJECT={len(reject)}"
    )

    return {
        **result,
        **extra,
        "evidence_cards": evidence_cards,
        "messages": messages,
        "errors": errors,
        "current_agent": "qualification",
    }
