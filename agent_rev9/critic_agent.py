from __future__ import annotations

import json
import logging
import re

from .small_molecule_critic import audit_candidate

log = logging.getLogger(__name__)


# ── LLM 기반 falsification 판정 ───────────────────────────

_CRITIC_SYSTEM = """당신은 신약개발 파이프라인의 독립 검증자(Critic)입니다.
제공된 분석 결과와 falsification 공격 결과를 비판적으로 검토하고 판정을 내립니다.

## 판정 종류
- APPROVE: 모든 근거가 일관되고 충분함
- REPLAN: 분석 설계에 근본적 문제 발견 → Planner로 회귀 필요
- HOLD: 불확실성이 높아 판단 보류 (재계획 불필요)
- REJECT: 치명적 오류 또는 safety concern

## 반드시 아래 JSON만 출력하세요
{
  "verdict": "APPROVE | REPLAN | HOLD | REJECT",
  "reason": "판정 근거 1~2문장",
  "feedback": "REPLAN 시 Planner에게 전달할 구체적 개선 방향 (REPLAN이 아니면 빈 문자열)",
  "issues_found": ["발견된 문제 목록"],
  "regression_target": "planner | discovery | qualification | null"
}

## 검증 체크리스트
1. Contrast 타당성: 아형 특이성을 검증하기에 충분한가?
2. Cell-of-origin: bulk DEG가 종양세포 유래인지 확인되었는가?
3. Safety veto: 모달리티별 veto와 HOLD 사유가 반영되었는가? 정상조직 고발현만으로 모든 모달리티를 REJECT하지 마세요.
4. Target-modality 정합성: 표적 위치와 모달리티가 일치하는가?
5. 무근거 주장: Evidence Card 없이 도출된 결론이 있는가?
6. ADVANCE 후보 0건: 파이프라인 재설계가 필요한가?

주의: APPROVE는 현재 단계의 계산 근거 검토 승인입니다. 실험 효능 확인으로 해석하지 마세요.
계산 결과와 기존 판정을 변경하지 마세요. 검증 미실행/자료 부족은 실제 반증과 구분하세요.
재계획은 실행 가능한 분석 수정이 있을 때만 요청하세요. 외부 검증 자료 부족만으로 동일 분석을 반복하지 마세요."""


def _parse_critic_json(text: str) -> dict | None:
    cleaned = re.sub(r"```json\s*", "", text)
    cleaned = re.sub(r"```\s*", "", cleaned).strip()
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        match = re.search(r"\{[\s\S]*\}", cleaned)
        if match:
            try:
                return json.loads(match.group())
            except json.JSONDecodeError:
                pass
    return None


def _llm_falsification(state: dict) -> dict:
    """LLM 기반 falsification 판정. LLM 불가 시 None 반환."""
    from .llm_client import LLMUnavailable, complete

    try:
        from .falsification_attacks import run_all_attacks
        attacks = run_all_attacks(state.get("advance_targets") or [], state)
    except Exception as exc:
        log.warning("[Critic] falsification_attacks 오류: %s", exc)
        attacks = {}

    advance = state.get("advance_targets") or []
    hold = state.get("hold_targets") or []
    reject = state.get("reject_targets") or []

    evidence = {
        "advance_targets": [
            {k: t.get(k) for k in ("gene_name", "modality", "tier", "bio_confidence",
                                    "intervention_role", "safety_veto")}
            for t in advance[:10]
        ],
        "counts": {"advance": len(advance), "hold": len(hold), "reject": len(reject)},
        "falsification_attacks": attacks,
        "gate_status": state.get("gate_status", "NOT_YET_EVALUATED"),
        "review_phase": "before_design_gate",
        "contrasts": state.get("contrasts", []),
        "research_plan": state.get("research_plan", {}),
        "evidence_cards": state.get("evidence_cards", [])[-30:],
        "qualification_evidence": {
            t.get("gene_name"): {k: (state.get("qualification_results", {}).get(t.get("gene_name"), {})).get(k)
                for k in ("intervention_role", "evidence_cards", "missing_evidence", "tier_results", "safety_veto")}
            for t in advance[:10]
        },
        "cell_contexts": {t.get("gene_name"): (state.get("cell_contexts") or state.get("cell_context_results") or {}).get(t.get("gene_name"), {}) for t in advance},
        "depmap_results": {t.get("gene_name"): (state.get("depmap_results") or {}).get(t.get("gene_name"), {}) for t in advance},
        "r_analysis_success": (state.get("r_analysis_result") or {}).get("success"),
        "input_manifest": (state.get("r_analysis_result") or {}).get("input_manifest", {}),
        "de_summary": (state.get("r_analysis_result") or {}).get("de_summary", {}),
        "analysis_parameters": (state.get("r_analysis_result") or {}).get("analysis_parameters", {}),
        "replan_count": state.get("replan_count", 0),
    }

    try:
        response = complete(
            instructions=_CRITIC_SYSTEM,
            evidence=evidence,
            role="critic",
        )
        result = _parse_critic_json(response.text)
        if result:
            log.info("[Critic] LLM 판정: %s (%.1fs)", result.get("verdict"), response.elapsed)
        return result
    except LLMUnavailable as exc:
        log.warning("[Critic] LLM 불가: %s", exc)
        return None


def _apply_deterministic_rules(state: dict) -> dict:
    if any(t.get("safety_verdict") == "REJECT" for t in state.get("advance_targets", [])):
        return {
            "verdict": "REJECT",
            "reason": "Target safety veto",
            "issues_found": ["TARGET_SAFETY_VETO"],
        }

    result = state.get("screening_result", {})
    candidates = result.get("candidates", [])
    accepted = []
    issues = []
    for candidate in candidates:
        if candidate.get("is_mock", True):
            issues.append("MOCK_CANDIDATE")
            continue
        if candidate.get("route") == "small_molecule" and "consensus" in candidate:
            audit = audit_candidate(candidate)
            passed = audit["verdict"] in {"GO", "GO_WITH_WARNING"}
            issues.extend(audit["reason_codes"])
        elif candidate.get("route") == "protein_binder_rfd3":
            passed = (
                candidate.get("af3_validation", {}).get("status") == "success"
                and candidate.get("primary_structural_screening", {}).get("verdict") == "PASS"
                and candidate.get("structural_validation", {}).get("passed") is True
            )
            issues.extend(candidate.get("failure_reasons", []))
        else:
            passed = (
                candidate.get("route") == "small_molecule"
                and candidate.get("pose_screening", {}).get("verdict") == "PASS"
                and candidate.get("cross_model_validation", {}).get("verdict") == "PASS"
                and candidate.get("admet_screening", {}).get("verdict") == "PASS"
            )
            issues.extend(candidate.get("failure_reasons", []))
        if passed:
            accepted.append(candidate["candidate_id"])

    verdict = "APPROVE" if accepted else "HOLD"
    return {
        "verdict": verdict,
        "reason": f"{len(accepted)} candidates passed the configured computational checks",
        "issues_found": list(dict.fromkeys(issues)),
        "accepted_candidate_ids": accepted,
    }


def critic_node(state: dict) -> dict:
    """Critic 노드: 결정론적 룰 + LLM falsification → APPROVE/HOLD/REPLAN/REJECT."""
    from .configuration import get_config

    # 1. 결정론적 룰 (safety veto, mock candidate 등)
    rule_result = _apply_deterministic_rules(state)

    # 설계 전 qualification-only 단계: 후보 목록만 있고 설계 결과 없음
    # → LLM이 ADVANCE 후보의 타당성을 검토
    verdict = rule_result["verdict"]
    reason = rule_result["reason"]
    feedback = ""
    issues = rule_result.get("issues_found", [])

    # 2. LLM falsification (agent_logging 설정 시, 또는 qualification 단계)
    config = get_config()
    run_llm = getattr(config, "agent_logging", True)
    if run_llm:
        llm_result = _llm_falsification(state)
        if llm_result:
            llm_verdict = llm_result.get("verdict", "HOLD")
            # 결정론적 REJECT은 LLM이 덮어쓸 수 없음
            if verdict != "REJECT":
                verdict = llm_verdict
            reason = llm_result.get("reason", reason)
            feedback = llm_result.get("feedback", "")
            issues = list(dict.fromkeys(issues + llm_result.get("issues_found", [])))

    # 3. REPLAN 상한: 동일 후보 2회 초과 회귀 방지
    replan_count = state.get("replan_count", 0)
    if verdict == "REPLAN" and replan_count >= 2:
        verdict = "HOLD"
        reason = f"REPLAN 상한 도달 (횟수={replan_count}) — HOLD로 전환"
        feedback = ""
        issues.append("REPLAN_LIMIT_REACHED")
        log.info("[Critic] REPLAN 상한 도달 → HOLD")

    final = {**rule_result, "verdict": verdict, "reason": reason,
             "issues_found": issues}

    log.info("[Critic] 최종 판정: %s", verdict)

    return {
        "critic_verdict":    verdict,
        "critic_reason":     reason,
        "critic_feedback":   feedback,
        "critic_count":      state.get("critic_count", 0) + 1,
        "critic_result":     final,
        "messages":          [f"[Critic] {verdict}: {reason[:120]}"],
        "review_issues":     issues,
        "errors":            [],
        "current_agent":     "critic",
    }


def review_modality(result):
    if result['modality'] == 'DEGRADER' and result['summary'].get('input_mode') == 'published_full_molecule':
        return review_degrader_reference(result)
    decision = result['validation_decision']
    reasons = list(dict.fromkeys(list(result.get('blockers', [])) + list(result.get('validation', {}).get('reasons', []))))
    if decision == 'ADVANCE' and (result['is_mock'] or result['execution_status'] != 'COMPLETED' or result.get('blockers')):
        decision = 'HOLD'
        reasons.append('incomplete_or_mock_evidence')
    return {'decision': decision, 'reasons': reasons,
            'experimental_efficacy_established': False}


def review_demo_modality(result):
    if result['modality'] == 'DEGRADER' and result['summary'].get('input_mode') == 'published_full_molecule':
        return review_degrader_reference(result)
    validation = result['validation']
    decision = {'ADVANCE': 'PASS_WITH_WARNINGS', 'HOLD': 'HOLD',
                'REJECT': 'NO_GO', 'NOT_EVALUATED': 'NOT_EVALUATED'}[validation['decision']]
    reasons = list(result.get('blockers', [])) + list(validation.get('missing_evidence', []))
    if result.get('is_mock'):
        decision = 'NOT_EVALUATED'
        reasons.append('mock_not_allowed_in_competition_demo')
    if result['modality'] == 'DE_NOVO_BINDER':
        development = result.get('demo_metrics', {}).get('developability', {})
        reasons.extend(development.get('warnings', []))
        if validation.get('structural_screening') != 'PASS' or development.get('status') == 'HOLD':
            decision = 'HOLD' if validation['decision'] != 'REJECT' else 'NO_GO'
    return {'decision': decision, 'reasons': list(dict.fromkeys(reasons)),
            'scope': 'computational_demo_screening_not_therapeutic_efficacy',
            'experimental_efficacy_established': False}


def review_degrader_reference(result):
    evidence = result.get('validation', {})
    assessment = evidence.get('computational_validation')
    if assessment:
        decision = assessment['decision']
        if result.get('is_mock') or result['execution_status'] != 'COMPLETED' or result.get('blockers'):
            decision = 'HOLD' if decision == 'ADVANCE' else decision
        return dict(decision=decision, reasons=assessment['reasons'], scope=assessment['scope'],
                    novel_design=False, degradation_efficacy='NOT_EVALUATED', experimental_efficacy_established=False)
    warnings = result['summary'].get('qc', {}).get('warning_message', [])
    if evidence.get('reference_evaluation_supported'):
        decision = 'PASS_WITH_WARNINGS' if warnings or result['blockers'] else 'PASS_FOR_REFERENCE_EVALUATION'
    else:
        decision = 'HOLD' if result['execution_status'] in {'COMPLETED', 'PARTIAL', 'FAILED'} else 'NOT_EVALUATED'
    return dict(decision=decision, reasons=list(dict.fromkeys(warnings + result['blockers'])),
                scope='published_reference_chemical_preparation_only', novel_design=False,
                degradation_efficacy='NOT_EVALUATED', experimental_efficacy_established=False)
