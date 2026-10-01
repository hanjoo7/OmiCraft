"""
User Selection Gate — 워크플로우 v5 Section 08

사용자가 표적·치료 방향·방식·자원을 확정한 후에만 고비용 설계 job이 시작된다.
선택 전 설계 job 시작 금지.

선택 이벤트:
  - route_id: modality 식별자 (예: "GENE1:SMALL_MOLECULE")
  - selected_by: "user" | "auto" (자동 실행 전용 테스트용)
  - timestamp: UTC ISO 형식
  - rationale: 선택 근거 (사용자 입력 또는 규칙 기반)
  - constraints: 자원 제약 (GPU 시간, 예산 등)

Tier에 따른 진입 허가:
  - 1A/1B/2A/2B: 사용자 선택 후 허가
  - 3:           Tier 3는 고비용 자동 설계 기본 비활성 (사용자가 명시적 승인 필요)
  - HOLD/REJECT: 진입 불가
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

SelectionStatus = Literal[
    "SELECTED",         # 사용자가 선택함
    "PENDING",          # 아직 선택 대기 중
    "HOLD",             # 필수 조건 미충족으로 선택 불가
    "REJECTED",         # REJECT Tier — 진입 불가
    "NOT_APPLICABLE",   # 해당 modality에 적용 불가
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ── 선택 이벤트 생성 ──────────────────────────────────────

def create_selection_event(
    gene: str,
    modality: str,
    route_id: str,
    tier: str,
    selected_by: str = "user",
    rationale: str = "",
    constraints: dict | None = None,
) -> dict:
    """선택 이벤트를 생성한다.

    Args:
        gene: 선택된 유전자.
        modality: 선택된 modality (SMALL_MOLECULE / DE_NOVO_BINDER / ADC 등).
        route_id: "{gene}:{modality}" 형식 route 식별자.
        tier: 현재 Tier (1A/1B/2A/2B/3/HOLD/REJECT).
        selected_by: "user" | "auto" (테스트 전용).
        rationale: 선택 근거.
        constraints: 자원 제약 dict.

    Returns:
        selection_event dict.
    """
    return {
        "gene": gene,
        "candidate": gene,
        "reason": rationale,
        "modality": modality,
        "route_id": route_id,
        "tier_at_selection": tier,
        "selected_by": selected_by,
        "timestamp": _utc_now(),
        "rationale": rationale,
        "constraints": constraints or {},
        "status": "SELECTED",
    }


# ── 진입 허가 판단 ────────────────────────────────────────

def check_design_gate(
    gene: str,
    modality: str,
    tier: str,
    selection_events: list[dict],
    allow_tier3_auto: bool = False,
) -> dict:
    """설계 job 시작 허가 여부를 판단한다.

    워크플로우 원칙:
      - 사용자 선택 전 고비용 설계 job 시작 금지.
      - Tier 3는 고비용 자동 설계 기본 비활성.
      - HOLD/REJECT는 설계 진입 불가.
      - 선택된 경로의 필수 요건 미충족 시 명확히 HOLD.

    Returns:
        {
          "allowed": bool,
          "reason": str,
          "gate_status": str,  # ALLOWED | BLOCKED_NO_SELECTION | BLOCKED_HOLD
                               # BLOCKED_REJECT | BLOCKED_TIER3 | BLOCKED_NOT_APPLICABLE
        }
    """
    # REJECT / NOT_APPLICABLE — 무조건 차단
    if tier == "REJECT":
        return {"allowed": False,
                "reason": f"{gene}:{modality} B0 또는 critical veto — REJECT Tier",
                "gate_status": "BLOCKED_REJECT"}
    if tier == "NOT_APPLICABLE":
        return {"allowed": False,
                "reason": f"{gene}:{modality} 해당 modality 적용 불가 (N/A)",
                "gate_status": "BLOCKED_NOT_APPLICABLE"}

    # HOLD — 필수 항목 미충족
    if tier == "HOLD":
        return {"allowed": False,
                "reason": f"{gene}:{modality} HOLD — 필수 항목 UNKNOWN. 추가 근거 확보 후 재평가",
                "gate_status": "BLOCKED_HOLD"}

    if tier not in {"1A", "1B", "2A", "2B", "3"}:
        return {"allowed": False, "reason": "Qualification Tier is unknown",
                "gate_status": "BLOCKED_UNKNOWN_TIER"}

    # Tier 3 — 기본 비활성
    if tier == "3" and not allow_tier3_auto:
        return {"allowed": False,
                "reason": (f"{gene}:{modality} Tier 3 탐색 경로 — "
                           "고비용 자동 설계 기본 비활성. 사용자 명시적 승인 필요."),
                "gate_status": "BLOCKED_TIER3"}

    # 사용자 선택 확인
    route_id = f"{gene}:{modality}"
    selected = any(
        ev.get("route_id") == route_id and
        ev.get("status") == "SELECTED" and
        ev.get("selected_by") == "user"  # auto는 테스트 전용
        for ev in selection_events
    )
    if not selected:
        return {"allowed": False,
                "reason": (f"{gene}:{modality} 사용자 선택 전 설계 실행 불가. "
                           f"route_id={route_id!r}를 선택 후 재실행하세요."),
                "gate_status": "BLOCKED_NO_SELECTION"}

    return {"allowed": True,
            "reason": f"{gene}:{modality} 사용자 선택 확인 — 설계 진입 허가",
            "gate_status": "ALLOWED"}


# ── 선택 이벤트 저장/로드 ─────────────────────────────────

def save_selection_events(
    events: list[dict],
    output_path: str | Path,
) -> str:
    """선택 이벤트를 JSON 파일에 저장한다."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"events": events, "updated_at": _utc_now()},
                   ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return str(path)


def load_selection_events(path: str | Path) -> list[dict]:
    """선택 이벤트를 파일에서 로드한다. 파일 없으면 빈 목록 반환."""
    p = Path(path)
    if not p.is_file():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return data.get("events", [])
    except (OSError, json.JSONDecodeError):
        return []


# ── State에서 선택 이벤트 조회 ────────────────────────────

def get_selection_events_from_state(state: dict) -> list[dict]:
    """OmiCraftState에서 선택 이벤트를 추출한다."""
    return state.get("user_selection_events", [])


def apply_user_selection(
    state: dict,
    gene: str,
    modality: str,
    rationale: str = "",
    constraints: dict | None = None,
    tier: str = "UNKNOWN",
) -> dict:
    """State에 사용자 선택 이벤트를 추가하고 업데이트된 State 조각을 반환.

    orchestrator에서 호출:
        state.update(apply_user_selection(state, gene, modality, ...))
    """
    route_id = f"{gene}:{modality}"
    event = create_selection_event(
        gene=gene,
        modality=modality,
        route_id=route_id,
        tier=tier,
        selected_by="user",
        rationale=rationale,
        constraints=constraints,
    )
    events = list(get_selection_events_from_state(state))
    # 동일 route_id 이전 이벤트 교체
    events = [e for e in events if e.get("route_id") != route_id]
    events.append(event)

    return {
        "user_selection_events": events,
        "selected_route": f"{gene}:{modality}",
        "messages": [
            f"[UserSelection] {gene}:{modality} 선택됨 (Tier={tier})"
        ],
    }


# ── Orchestrator용 gate 노드 ──────────────────────────────

def user_selection_gate_node(state: dict) -> dict:
    """Orchestrator 노드: 모든 ADVANCE 표적에 대해 gate를 확인한다.

    선택된 표적이 없으면 HOLD 반환.
    선택된 표적이 있으면 gate를 통과한 표적 목록을 반환.
    """
    events = get_selection_events_from_state(state)
    advance_targets = state.get("advance_targets", [])

    if not advance_targets:
        return {
            "gate_status": "HOLD",
            "gate_reason": "ADVANCE 표적 없음 — qualification 먼저 실행 필요",
            "design_allowed_targets": [],
            "messages": ["[UserSelectionGate] ADVANCE 표적 없음"],
        }

    if not events:
        # 선택 이벤트 없음 — 표적 목록 요약 후 HOLD
        target_summary = ", ".join(
            f"{t.get('gene_name', '?')}:{t.get('modality', '?')}"
            for t in advance_targets[:5]
        )
        return {
            "gate_status": "PENDING",
            "gate_reason": (
                "사용자 선택 대기 중. "
                f"선택 가능한 표적: [{target_summary}]"
            ),
            "design_allowed_targets": [],
            "messages": [
                "[UserSelectionGate] 사용자 선택 필요",
                f"  선택 가능: {target_summary}",
                "  선택 방법: apply_user_selection(state, gene, modality) 호출",
            ],
        }

    # gate 통과 여부 확인
    allowed: list[dict] = []
    blocked: list[dict] = []

    for target in advance_targets:
        gene     = target.get("gene_name", "")
        modality = target.get("modality", "")
        tier     = target.get("tier") or target.get("best_tier") or "UNKNOWN"

        gate = check_design_gate(
            gene=gene,
            modality=modality,
            tier=tier,
            selection_events=events,
        )
        if gate["allowed"]:
            allowed.append({**target, "gate": gate})
        else:
            blocked.append({**target, "gate": gate})

    messages = [
        f"[UserSelectionGate] 허가: {len(allowed)}개, 차단: {len(blocked)}개"
    ]
    for t in allowed:
        messages.append(
            f"  ALLOWED: {t.get('gene_name')}:{t.get('modality')} (Tier={t.get('tier')})"
        )
    for t in blocked[:3]:
        messages.append(
            f"  BLOCKED: {t.get('gene_name')}:{t.get('modality')} — "
            f"{t['gate']['reason'][:60]}"
        )

    return {
        "gate_status": "ALLOWED" if allowed else "PENDING",
        "gate_reason": f"{len(allowed)}개 표적 허가" if allowed
                       else "선택된 표적 없음 — 설계 실행 전 선택 필요",
        "design_allowed_targets": allowed,
        "blocked_targets": blocked,
        "messages": messages,
    }
