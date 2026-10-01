"""
Agent ① Coordinator-Planner
연구질문 → 질환/아형/비교군 추출 → 실행 DAG / 예산 / 중단조건 생성

LLM 담당:
  - 연구질문 파싱 → disease / subtype / contrasts / analysis_steps / budget
  - V1~V5 self-validation (질환명, contrast 수, DAG 구성, 예산, 데이터 가용성)
  - REPLAN 시 critic_feedback을 반영한 계획 수정

tool 담당 (결정론적):
  - contrast_spec.resolve_contrasts() → ContrastSpec / ResearchPlan 확정
  - config에서 데이터 경로 가용성 확인

출력 state 필드:
  disease, subtype, contrasts, research_plan, contrast_manifest,
  dag_plan, budget, stop_conditions, evidence_cards, messages, errors
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

from .state import OmiCraftState

log = logging.getLogger(__name__)


# ── 프롬프트 ──────────────────────────────────────────────

_SYSTEM = """당신은 신약개발 연구 계획을 수립하는 전문가입니다.
연구질문을 분석하여 오믹스 기반 표적 발굴 파이프라인의 실행 계획을 생성합니다.

반드시 아래 JSON 형식으로만 응답하세요:

{
  "disease": "질환명 (영문)",
  "subtype": "아형명 (영문, 없으면 null)",
  "cohort": "사용할 코호트 설명",
  "contrasts": [
    "contrast 1 설명",
    "contrast 2 설명"
  ],
  "analysis_steps": [
    {"step": 1, "agent": "discovery",      "task": "차등발현분석",        "tools": ["DESeq2"]},
    {"step": 2, "agent": "discovery",      "task": "경로분석",            "tools": ["fgsea"]},
    {"step": 3, "agent": "qualification",  "task": "세포맥락 판정",        "tools": ["CELLxGENE"]},
    {"step": 4, "agent": "qualification",  "task": "DepMap 의존성",        "tools": ["DepMap"]},
    {"step": 5, "agent": "qualification",  "task": "Tier 판정 + 모달리티", "tools": ["ChEMBL"]},
    {"step": 6, "agent": "critic",         "task": "최종 검증",            "tools": ["LLM"]}
  ],
  "budget": {
    "max_candidates": 100,
    "max_advance_targets": 5,
    "llm_token_budget": 600000,
    "gpu_hours": 6.0
  },
  "stop_conditions": [
    "safety veto 발생 시 해당 후보 즉시 REJECT",
    "동일 후보 Critic 회귀 2회 초과 시 HOLD",
    "ADVANCE 후보 0건이면 contrast 확장 후 재실행"
  ],
  "reasoning": "계획 수립 근거 설명"
}

주의사항:
- 비교군은 최소 2개 contrast를 설정하세요
- contrast에는 반드시 데이터에 존재하는 아형 이름을 사용하세요
  (사용 가능: TNBC, Luminal_A, Luminal_B, HER2_enriched)
- "normal tissue" / "adjacent normal"은 별도 코호트가 없으면 사용하지 마세요
- 비판적 피드백이 제공된 경우 반드시 해당 문제를 반영하세요
"""

_USER_TEMPLATE = """\
## 연구질문
{question}

## 사전 구축 데이터
- TCGA-BRCA RNA-Seq: {n_samples} 샘플 (TNBC {n_tnbc}명, non-TNBC {n_nontnbc}명)
- Clinical: ER/PR/HER2 수용체 상태, 생존 데이터
- Mutation (MAF), CNV
- 정상조직 참조: GTEx (54 tissues), HPA (20,162 genes)
- 의존성: DepMap CRISPR (1,208 cell lines)
- 약물: ChEMBL 35 (SQLite)

{feedback_section}위 데이터를 활용한 실행 계획을 JSON으로 생성하세요."""

_FEEDBACK_SECTION = """\
## Critic 피드백 (이전 계획의 문제점 — 반드시 반영하세요)
{feedback}

재계획 횟수: {replan_count}/2

"""


# ── JSON 파싱 헬퍼 ─────────────────────────────────────────

def _parse_json(text: str) -> dict | None:
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


# ── Self-validation ────────────────────────────────────────

def _validate_plan(plan: dict, config) -> tuple[list[str], list[str]]:
    """V1~V5 self-validation. 반환: (messages, errors)"""
    messages: list[str] = []
    errors: list[str] = []

    # V1: 질환명 추출 여부
    disease = plan.get("disease", "")
    if not disease:
        errors.append("[Planner:V1] 질환명 추출 실패 — 연구질문에서 질환을 식별하지 못했습니다")
    else:
        messages.append(f"[Planner:V1] 질환명 확인: {disease}")

    # V2: contrast 최소 2개 (자동 보완 — 실패가 아님)
    contrasts = plan.get("contrasts", [])
    if len(contrasts) < 2:
        if len(contrasts) == 1:
            contrasts.append("TNBC vs non-TNBC breast cancer (subtype specificity)")
        else:
            contrasts = [
                "TNBC tumor vs non-TNBC breast cancer",
                "TNBC vs Luminal_A breast cancer",
            ]
        plan["contrasts"] = contrasts
        messages.append(f"[Planner:V2] contrast {len(contrasts)}개로 보완")
    else:
        messages.append(f"[Planner:V2] Contrast {len(contrasts)}개 — 적절")

    # V3: 필수 에이전트 포함 여부
    steps = plan.get("analysis_steps", [])
    agent_types = {s.get("agent", "") for s in steps}
    missing = {"discovery", "qualification"} - agent_types
    if missing:
        errors.append(f"[Planner:V3] 필수 에이전트 누락: {missing}")
    else:
        messages.append(f"[Planner:V3] DAG {len(steps)} steps — 적절")

    # V4: 예산 범위
    budget = plan.get("budget", {})
    max_cand = budget.get("max_candidates", 0)
    if max_cand < 10:
        messages.append(f"[Planner:V4] 후보 예산 {max_cand}개 → 100으로 조정")
        budget["max_candidates"] = 100
        plan["budget"] = budget
    elif max_cand > 500:
        messages.append(f"[Planner:V4] 후보 예산 {max_cand}개 → 200으로 조정")
        budget["max_candidates"] = 200
        plan["budget"] = budget
    else:
        messages.append(f"[Planner:V4] 후보 예산 {max_cand}개 — 적절")

    # V5: 데이터 가용성 확인
    missing_data = []
    for attr in ("counts_rds", "metadata_rds"):
        val = getattr(config.data, attr, None)
        if val and not Path(val).is_file():
            missing_data.append(attr)
    if missing_data:
        errors.append(f"[Planner:V5] 필수 데이터 미확인: {missing_data} — 경로 점검 필요")
    else:
        messages.append("[Planner:V5] 필수 데이터 가용성 확인 완료")

    return messages, errors


# ── 기본 계획 (LLM 불가 시 fallback) ──────────────────────

def _default_plan(state: OmiCraftState) -> dict:
    return {
        "disease": state.get("disease", "breast_cancer"),
        "subtype": state.get("subtype", "TNBC"),
        "cohort": "TCGA-BRCA",
        "contrasts": [
            "TNBC vs non-TNBC breast cancer",
            "TNBC vs Luminal_A breast cancer",
        ],
        "analysis_steps": [
            {"step": 1, "agent": "discovery",     "task": "차등발현분석",        "tools": ["DESeq2"]},
            {"step": 2, "agent": "discovery",     "task": "경로분석",            "tools": ["fgsea"]},
            {"step": 3, "agent": "qualification", "task": "세포맥락 판정",        "tools": ["CELLxGENE"]},
            {"step": 4, "agent": "qualification", "task": "DepMap 의존성",        "tools": ["DepMap"]},
            {"step": 5, "agent": "qualification", "task": "Tier 판정 + 모달리티", "tools": ["ChEMBL"]},
            {"step": 6, "agent": "critic",        "task": "최종 검증",            "tools": ["LLM"]},
        ],
        "budget": {"max_candidates": 100, "max_advance_targets": 5},
        "stop_conditions": [
            "safety veto 발생 시 해당 후보 즉시 REJECT",
            "동일 후보 Critic 회귀 2회 초과 시 HOLD",
        ],
        "reasoning": "LLM 불가 — 기본 계획 적용",
    }


# ── 메인 노드 ──────────────────────────────────────────────

def planner_node(state: OmiCraftState) -> dict:
    """Coordinator-Planner node: 연구질문 → 실행 계획 확정."""
    from .configuration import get_config
    from .llm_client import LLMUnavailable, complete

    question = state.get("research_question", "")
    if not question:
        return {
            "errors": ["[Planner] 연구질문이 비어 있습니다"],
            "messages": ["[Planner] 연구질문 없음 — 파이프라인 중단"],
            "current_agent": "planner",
        }

    config = get_config()
    critic_feedback = state.get("critic_feedback", "")
    replan_count = state.get("replan_count", 0)
    if state.get("critic_verdict") == "REPLAN":
        replan_count += 1

    log.info("[Planner] 시작 — replan_count=%d", replan_count)

    # LLM 호출
    plan: dict | None = None
    try:
        feedback_section = (
            _FEEDBACK_SECTION.format(feedback=critic_feedback, replan_count=replan_count)
            if critic_feedback else ""
        )
        user_msg = _USER_TEMPLATE.format(
            question=question,
            n_samples=getattr(getattr(config, "cohort", None), "n_samples", "?"),
            n_tnbc=getattr(getattr(config, "cohort", None), "n_tnbc", "?"),
            n_nontnbc=getattr(getattr(config, "cohort", None), "n_nontnbc", "?"),
            feedback_section=feedback_section,
        )
        response = complete(
            instructions=_SYSTEM,
            evidence={"question": question, "feedback": critic_feedback,
                      "replan_count": replan_count, "prompt": user_msg},
            role="planner",
        )
        plan = _parse_json(response.text)
        log.info("[Planner] LLM 응답 수신 (%.1fs)", response.elapsed)
    except LLMUnavailable as exc:
        log.warning("[Planner] LLM 불가 — 기본 계획 사용: %s", exc)

    if not plan:
        plan = _default_plan(state)

    # Self-validation
    val_msgs, val_errors = _validate_plan(plan, config)

    # ContrastSpec 확정 (결정론적 tool)
    contrasts_dicts: list = []
    contrast_manifest = ""
    research_plan: dict = {}
    try:
        from .contrast_spec import resolve_contrasts
        output_dir = Path(config.data.agent_results) / "contrast_spec"
        resolved = resolve_contrasts(
            question=question,
            sample_manifest=(state.get("r_analysis_result") or {}).get("input_manifest") or state.get("sample_manifest"),
            data_flags=state.get("data_flags") or {"rnaseq": True},
            plan_id=state.get("run_id", "run001"),
            output_dir=str(output_dir),
        )
        contrasts_dicts = [c.to_dict() for c in resolved.contrasts]
        contrast_manifest = str(output_dir / "contrast_manifest.tsv")
        research_plan = resolved.to_dict()
        val_msgs.append(f"[Planner] ContrastSpec 확정 — {len(contrasts_dicts)}개 contrast")
    except Exception as exc:
        log.warning("[Planner] ContrastSpec 오류 (계속 진행): %s", exc)
        contrasts_dicts = plan.get("contrasts", [])
        val_errors.append(f"[Planner] ContrastSpec 오류: {exc}")

    # Evidence card
    validation_passed = not any(e for e in val_errors if "V1" in e or "V3" in e)
    evidence_cards = [{
        "id": f"EV-PLAN-{replan_count:02d}",
        "source": "Planner",
        "type": "SUPPORTIVE" if validation_passed else "CONTRADICTORY",
        "content": (
            f"연구계획: {plan.get('disease')}, contrast {len(contrasts_dicts)}개, "
            f"validation {'PASS' if validation_passed else 'FAIL'}, "
            f"replan_count={replan_count}"
        ),
        "gene": "",
        "confidence": 0.95 if validation_passed else 0.6,
    }]

    all_msgs = (
        [f"[Planner] 시작 (replan={replan_count})"]
        + val_msgs
        + (["[Planner] self-validation PASS"] if validation_passed
           else ["[Planner] self-validation 일부 실패 — 자동 보완 적용"])
    )

    log.info(
        "[Planner] 완료 — disease=%s, contrasts=%d, errors=%d",
        plan.get("disease"), len(contrasts_dicts), len(val_errors),
    )

    return {
        "disease":           plan.get("disease", ""),
        "subtype":           plan.get("subtype", ""),
        "contrasts":         contrasts_dicts,
        "research_plan":     research_plan,
        "contrast_manifest": contrast_manifest,
        "dag_plan":          {"steps": plan.get("analysis_steps", [])},
        "budget":            plan.get("budget", {}),
        "stop_conditions":   plan.get("stop_conditions", []),
        "evidence_cards":    evidence_cards,
        "messages":          all_msgs,
        "errors":            val_errors,
        "current_agent":     "planner",
        "replan_count":      replan_count,
        "critic_count":      state.get("critic_count", 0),
    }
