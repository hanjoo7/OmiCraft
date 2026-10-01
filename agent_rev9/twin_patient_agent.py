"""
시나리오 3: 쌍둥이 환자 대조 — 내성 우회 경로 + 병용 가설
동일 변이 · 다른 조직 → 서로 다른 표적 · 모달리티 · 병용 전략 도출

원칙:
  - 표적/모달리티/Tier 선택은 qualify_gene() + check_design_gate() 결과 기반 결정론적 선택
  - 병용/내성 근거는 qualification 결과의 intervention_role / safety_veto 필드에서 추론
  - LLM 호출 없음: 재현성 보장
"""

from __future__ import annotations

import json
import logging
import os

from .state import OmiCraftState
from .configuration import get_config
from .user_selection import check_design_gate

log = logging.getLogger(__name__)


# ── 환자 프로파일 구성 ────────────────────────────────────────

def _build_patients_from_state(state: OmiCraftState) -> dict:
    """state에서 환자 프로파일 구성.

    우선순위:
    1. state["twin_patients"] — 사용자가 직접 제공한 프로파일
    2. advance_targets 상위 2개의 서로 다른 expression/context로 대리 프로파일 생성
    3. 둘 다 없으면 빈 dict 반환
    """
    # 1. 사용자 제공 프로파일
    explicit = state.get("twin_patients")
    if explicit and isinstance(explicit, dict) and len(explicit) >= 2:
        return explicit

    # 2. advance_targets 기반 자동 생성
    advance = state.get("advance_targets", [])
    if not advance:
        return {}

    disease = state.get("disease", "unknown")
    subtype = state.get("subtype", "")

    up_genes   = [t for t in advance if "Up"   in t.get("direction", "")][:3]
    down_genes = [t for t in advance if "Down" in t.get("direction", "")][:3]

    up_expr   = ", ".join(f"{t['gene_name']} high" for t in up_genes)   or "발현 데이터 없음"
    down_expr = ", ".join(f"{t['gene_name']} low"  for t in down_genes) or "발현 데이터 없음"

    dep_genes = [
        t["gene_name"] for t in advance
        if (t.get("crispr_pct_dependent") or 0) > 50
    ][:2]
    dep_context = (
        f"{', '.join(dep_genes)} CRISPR 의존성 확인" if dep_genes else "CRISPR 의존성 미확인"
    )

    return {
        "A": {
            "id": "PT-A",
            "name": f"환자 A ({subtype} 전형)",
            "tissue": f"{disease} ({subtype})",
            "subtype": subtype,
            "mutation": "omics 기반 (state.advance_targets 참조)",
            "expression_profile": up_expr,
            "context": f"과발현 표적 중심 — {dep_context}",
            "preferred_targets": [t["gene_name"] for t in up_genes],
        },
        "B": {
            "id": "PT-B",
            "name": f"환자 B ({subtype} 대조)",
            "tissue": f"{disease} ({subtype} 저발현 아형)",
            "subtype": f"{subtype}_alt",
            "mutation": "omics 기반 (대조군 맥락)",
            "expression_profile": down_expr,
            "context": "저발현 표적 — 기능 회복 또는 대체 경로 의존",
            "preferred_targets": [t["gene_name"] for t in down_genes],
        },
    }


# ── Rule-based 표적/모달리티 선택 ────────────────────────────

_TIER_RANK = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4, "HOLD": 99, "REJECT": 100}
_B_RANK    = {"B1": 0, "B2": 1, "B3": 2, "B0": 99}


def _select_best_target(
    preferred_genes: list[str],
    advance_targets: list[dict],
    selection_events: list[dict],
) -> dict | None:
    """preferred_genes 순서 우선, 이후 Tier/B 기준으로 최선의 표적 선택.

    반환: CandidateTarget dict 또는 None
    """
    # preferred 유전자 우선
    preferred_set = set(preferred_genes)
    preferred = [t for t in advance_targets if t.get("gene_name") in preferred_set]
    rest       = [t for t in advance_targets if t.get("gene_name") not in preferred_set]

    def _score(t: dict) -> tuple:
        tier = t.get("tier") or t.get("best_tier") or t.get("Tier") or "HOLD"
        b    = t.get("bio_confidence") or t.get("B") or "B3"
        return (_TIER_RANK.get(tier, 50), _B_RANK.get(b, 50))

    candidates = sorted(preferred, key=_score) + sorted(rest, key=_score)

    # gate 통과하는 첫 번째 후보 선택
    for t in candidates:
        gene     = t.get("gene_name", "")
        modality = t.get("modality", "")
        tier     = t.get("tier") or t.get("best_tier") or t.get("Tier") or "HOLD"

        gate = check_design_gate(
            gene=gene,
            modality=modality,
            tier=tier,
            selection_events=selection_events,
            allow_tier3_auto=True,  # 환자 대조 분석은 탐색 목적 — Tier 3 허용
        )
        if gate["allowed"] or gate["gate_status"] == "BLOCKED_NO_SELECTION":
            # BLOCKED_NO_SELECTION: 유효한 후보이지만 선택 이벤트 없음 — 탐색용으로 허용
            return t

    return candidates[0] if candidates else None


def _build_combination_hint(target: dict, qr: dict) -> str:
    """qualification 결과에서 병용 파트너 힌트 추출."""
    gene = target.get("gene_name", "")
    gene_qr = qr.get(gene, {})

    intervention_role = (
        target.get("intervention_role")
        or gene_qr.get("intervention_role")
        or "UNKNOWN"
    )

    # pathway 정보에서 병용 힌트
    top_pathway = gene_qr.get("top_pathway") or ""

    hints: list[str] = []

    if intervention_role == "FUNCTION_MODULATION":
        hints.append("신호전달 경로 억제 — PI3K/AKT 또는 MAPK 경로 병용 고려")
    elif intervention_role == "CELL_TARGETING":
        hints.append("면역항암제 병용 (checkpoint inhibitor) 고려")
    elif intervention_role == "BIOMARKER_ONLY":
        hints.append("바이오마커 역할 — 동반 진단 활용 고려")

    if "cell cycle" in top_pathway.lower():
        hints.append("CDK4/6 억제제 병용 가능성")
    if "dna repair" in top_pathway.lower() or "brca" in top_pathway.lower():
        hints.append("PARP 억제제 병용 가능성")

    return "; ".join(hints) if hints else "추가 근거 확보 필요"


def _build_resistance_hint(target: dict, qr: dict) -> str:
    """qualification 결과에서 내성 우회 경로 힌트 추출."""
    gene = target.get("gene_name", "")
    gene_qr = qr.get(gene, {})

    modality = target.get("modality", "")
    tier = target.get("tier") or target.get("best_tier") or target.get("Tier") or ""

    hints: list[str] = []

    # 낮은 Tier → 내성 우회 경로 필요
    if tier in ("2B", "3"):
        hints.append(f"Tier {tier}: 내성 가능성 높음 — 대체 modality 탐색 권장")

    # depmap 정보
    dep_status = gene_qr.get("dependency_status") or gene_qr.get("depmap_status") or ""
    if dep_status == "DISCORDANT":
        hints.append("CRISPR/RNAi 불일치 — 의존성 맥락 의존적; 병용으로 내성 완화 가능")
    elif dep_status == "CRISPR_ONLY":
        hints.append("CRISPR 의존성만 확인 — 단일요법 내성 시 대안 경로 필요")

    # safety veto 정보
    safety = gene_qr.get("safety_veto") or {}
    if any(v.get("has_veto") for v in safety.values() if isinstance(v, dict)):
        hints.append("안전성 우려 — 선택적 modality 또는 국소 전달 고려")

    return "; ".join(hints) if hints else "현재 데이터 기반 내성 우회 경로 미확인"


# ── 메인 노드 ─────────────────────────────────────────────────

def twin_patient_node(state: OmiCraftState) -> dict:
    """쌍둥이 환자 대조 분석 — rule-based 결정론적 선택."""
    messages: list[str] = []

    messages.append("[TwinPatient] 쌍둥이 환자 대조 시작")

    config = get_config()

    # 환자 프로파일
    twin_patients = _build_patients_from_state(state)
    if not twin_patients:
        messages.append("[TwinPatient] state에 환자 프로파일 없음 — 건너뜀")
        messages.append("[TwinPatient] state['twin_patients'] 또는 'advance_targets'를 제공하세요")
        return {"messages": messages, "current_agent": "twin_patient"}

    advance_targets   = state.get("advance_targets") or []
    selection_events  = state.get("user_selection_events") or []
    qr                = state.get("qualification_results") or {}

    messages.append(f"[TwinPatient] 환자 {len(twin_patients)}명 분석 (rule-based)")

    patient_results: dict[str, dict] = {}
    evidence_cards: list[dict] = []

    for pid, patient in twin_patients.items():
        messages.append(f"\n[TwinPatient] === {patient['name']} ({patient['tissue']}) ===")
        messages.append(f"[TwinPatient] 발현: {patient['expression_profile']}")

        preferred_genes = patient.get("preferred_targets", [])
        target = _select_best_target(preferred_genes, advance_targets, selection_events)

        if not target:
            messages.append(f"[TwinPatient:{pid}] 적합한 ADVANCE 표적 없음")
            continue

        gene     = target.get("gene_name", "?")
        modality = target.get("modality", "?")
        tier     = target.get("tier") or target.get("best_tier") or target.get("Tier") or "?"
        b        = target.get("bio_confidence") or target.get("B") or "?"

        combination      = _build_combination_hint(target, qr)
        resistance_bypass = _build_resistance_hint(target, qr)
        rationale = (
            f"Tier={tier}, B={b}, Modality={modality} — "
            f"intervention_role={target.get('intervention_role', 'UNKNOWN')}"
        )

        result = {
            "primary_target": gene,
            "modality": modality,
            "tier": tier,
            "bio_confidence": b,
            "rationale": rationale,
            "combination": combination,
            "resistance_bypass": resistance_bypass,
            "selection_method": "rule_based",
        }
        patient_results[pid] = result

        messages.append(f"[TwinPatient:{pid}] 표적: {gene} (Tier={tier}, B={b})")
        messages.append(f"[TwinPatient:{pid}] 모달리티: {modality}")
        messages.append(f"[TwinPatient:{pid}] 근거: {rationale}")
        if combination:
            messages.append(f"[TwinPatient:{pid}] 병용: {combination[:100]}")
        if resistance_bypass:
            messages.append(f"[TwinPatient:{pid}] 내성 우회: {resistance_bypass[:100]}")

        evidence_cards.append({
            "id": f"EV-TWIN-{pid}",
            "source": f"TwinPatient ({patient['name']})",
            "type": "SUPPORTIVE",
            "content": f"{patient['name']}: {gene} → {modality} (Tier={tier})",
            "gene": gene,
            "confidence": {"B1": 0.9, "B2": 0.7, "B3": 0.5}.get(b, 0.5),
        })

    # 대조 분석
    if len(patient_results) == 2:
        messages.append("\n[TwinPatient] === 환자 A vs B 대조 ===")
        ra = patient_results.get("A", {})
        rb = patient_results.get("B", {})

        same_target   = ra.get("primary_target") == rb.get("primary_target")
        same_modality = ra.get("modality") == rb.get("modality")

        messages.append(
            f"[TwinPatient] 표적 동일: {'예' if same_target else '아니오'} "
            f"(A={ra.get('primary_target')}, B={rb.get('primary_target')})"
        )
        messages.append(
            f"[TwinPatient] 모달리티 동일: {'예' if same_modality else '아니오'} "
            f"(A={ra.get('modality')}, B={rb.get('modality')})"
        )
        if ra.get("combination") != rb.get("combination"):
            messages.append(
                f"[TwinPatient] 병용 전략 차이: "
                f"A={ra.get('combination', '')[:60]}, "
                f"B={rb.get('combination', '')[:60]}"
            )
        messages.append(
            "[TwinPatient] → 동일 암종이지만 환자 맥락에 따라 다른 전략 도출 — Precision Medicine 시연"
        )

    # 저장
    result_dir = config.data.agent_results
    os.makedirs(result_dir, exist_ok=True)
    with open(
        os.path.join(result_dir, "twin_patient_output.json"), "w", encoding="utf-8"
    ) as f:
        json.dump(
            {"patients": twin_patients, "results": patient_results},
            f, indent=2, ensure_ascii=False, default=str,
        )

    messages.append("[TwinPatient] 쌍둥이 환자 대조 완료")

    return {
        "twin_patient_results": patient_results,
        "evidence_cards": evidence_cards,
        "messages": messages,
        "current_agent": "twin_patient",
    }
