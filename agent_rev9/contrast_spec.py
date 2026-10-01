"""
ContrastSpec / ResearchPlan / resolve_contrasts()

워크플로우 v5 Section 01 구현.

ContrastSpec: 비교 단위(TNBC vs Non-TNBC 등)의 상태와 메타.
ResearchPlan: 전체 분석 계획을 담는 불변 객체.
resolve_contrasts(): 질문·데이터 가용성으로 ContrastSpec 목록을 결정.

상태 종류:
  READY            — 실행 가능한 비교
  UNAVAILABLE      — 필요한 샘플 그룹이 데이터에 없음
  UNSUPPORTED      — 플랫폼·배치 문제로 지원 불가
  INVALID_DESIGN   — 설계 행렬 rank 문제 등 통계 오류
  PARTIAL          — 두 비교 중 하나만 완료된 중간 상태

positive_direction: 비교 기준 집합(log2FC 양수 = 이 그룹이 높음)
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, Optional

ContrastStatus = Literal[
    "READY", "UNAVAILABLE", "UNSUPPORTED", "INVALID_DESIGN", "PARTIAL"
]

SUPPORTED_GROUPS = {"TNBC", "Luminal_A", "Luminal_B", "HER2_enriched", "Non_TNBC", "Normal"}


@dataclass(frozen=True)
class ContrastSpec:
    """단일 비교의 명세."""

    contrast_id: str
    name: str                           # 예: "TNBC_vs_NonTNBC"
    positive_group: str                 # log2FC 양수 방향 그룹
    reference_group: str                # 대조군
    status: ContrastStatus
    reason: str = ""                    # UNAVAILABLE / UNSUPPORTED 사유
    min_n_positive: int = 0
    min_n_reference: int = 0
    actual_n_positive: Optional[int] = None
    actual_n_reference: Optional[int] = None
    positive_direction: str = ""        # "TNBC > Non-TNBC 방향이면 양수"
    result_status: Optional[str] = None # 실행 후 채움: dual_contrast_supported / partial_evidence / discordant

    def is_ready(self) -> bool:
        return self.status == "READY"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ResearchPlan:
    """전체 분석 실행 계획."""

    plan_id: str
    disease: str
    subtype: Optional[str]
    cohort_description: str
    contrasts: list[ContrastSpec]
    analysis_steps: list[dict] = field(default_factory=list)
    budget: dict = field(default_factory=dict)
    stop_conditions: list[str] = field(default_factory=list)
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    ruleset_version: str = "omicraft-v5"

    # 상태 추적
    execution_status: dict = field(default_factory=dict)  # step -> status

    def ready_contrasts(self) -> list[ContrastSpec]:
        return [c for c in self.contrasts if c.is_ready()]

    def partial_contrasts(self) -> list[ContrastSpec]:
        return [c for c in self.contrasts if c.status == "PARTIAL"]

    def all_ready(self) -> bool:
        return all(c.is_ready() for c in self.contrasts)

    def contrast_by_id(self, contrast_id: str) -> Optional[ContrastSpec]:
        for c in self.contrasts:
            if c.contrast_id == contrast_id:
                return c
        return None

    def to_dict(self) -> dict:
        return {
            "plan_id": self.plan_id,
            "disease": self.disease,
            "subtype": self.subtype,
            "cohort_description": self.cohort_description,
            "contrasts": [c.to_dict() for c in self.contrasts],
            "analysis_steps": self.analysis_steps,
            "budget": self.budget,
            "stop_conditions": self.stop_conditions,
            "created_at": self.created_at,
            "ruleset_version": self.ruleset_version,
            "execution_status": self.execution_status,
        }

    def save(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )


# ── 기본 분석 단계 정의 ─────────────────────────────────────

DEFAULT_ANALYSIS_STEPS = [
    {"step": 1, "id": "qc_batch_de", "agent": "discovery",
     "task": "QC + 배치 평가 + DESeq2", "tools": ["DESeq2", "SVA", "vst"]},
    {"step": 2, "id": "gsea", "agent": "discovery",
     "task": "Hallmark + GO BP GSEA + NES plots", "tools": ["fgsea", "MSigDB"]},
    {"step": 3, "id": "apear_network", "agent": "discovery",
     "task": "GO BP aPEAR network", "tools": ["aPEAR"]},
    {"step": 4, "id": "cell_context", "agent": "qualification",
     "task": "CELLxGENE + genesorteR + cell-of-origin", "tools": ["cellxgene_census", "genesorteR"]},
    {"step": 5, "id": "depmap", "agent": "qualification",
     "task": "DepMap CRISPR + RNAi dependency", "tools": ["DepMap"]},
    {"step": 6, "id": "safety_direction", "agent": "qualification",
     "task": "안전성·치료방향·기능근거", "tools": ["UniProt", "HPA", "OpenTargets", "ChEMBL"]},
    {"step": 7, "id": "shortlist_tier", "agent": "qualification",
     "task": "Biological confidence + B/M/D Tier + shortlist", "tools": ["tier_system"]},
    {"step": 8, "id": "user_selection", "agent": "coordinator",
     "task": "사용자 modality 선택", "tools": []},
    {"step": 9, "id": "design", "agent": "design",
     "task": "전용 설계", "tools": ["RFdiffusion3", "ProteinMPNN", "AlphaFold3"]},
    {"step": 10, "id": "critic", "agent": "critic",
     "task": "독립 검증", "tools": ["critic"]},
]

DEFAULT_BUDGET = {
    "max_candidates": 100,
    "max_advance_targets": 5,
    "max_design_targets": 3,
    "gpu_hours": 6.0,
    "max_critic_regressions": 2,
}

DEFAULT_STOP_CONDITIONS = [
    "critical safety veto 시 해당 유전자×가설×경로 범위 REJECT",
    "동일 후보 Critic 회귀 2회 초과 시 HOLD + 종료",
    "ADVANCE 후보 0건이면 HOLD — contrast 확장 또는 데이터 보강 필요",
    "사용자 선택 전 고비용 설계 job 시작 금지",
]


# ── 샘플 수 확인 헬퍼 ──────────────────────────────────────

def _count_groups(sample_manifest: dict | None) -> dict[str, int]:
    """sample_manifest에서 그룹별 샘플 수를 반환한다.

    manifest 형식: {"groups": {"TNBC": 120, "Luminal_A": 200, ...}}
    또는 {"samples": [{"group": "TNBC"}, ...]}
    """
    if not sample_manifest:
        return {}
    # R input_manifest uses group_counts and an integer samples count.
    if "group_counts" in sample_manifest:
        return dict(sample_manifest["group_counts"])
    # 집계된 형식
    if "groups" in sample_manifest:
        return dict(sample_manifest["groups"])
    # 샘플 목록 형식
    counts: dict[str, int] = {}
    for s in sample_manifest.get("samples", []):
        g = s.get("group") or s.get("subtype") or ""
        if g:
            counts[g] = counts.get(g, 0) + 1
    return counts


def _non_tnbc_groups(group_counts: dict[str, int]) -> list[str]:
    if group_counts.get("Non_TNBC", 0) > 0:
        return ["Non_TNBC"]
    known = {"Luminal_A", "Luminal_B", "HER2_enriched"}
    return [g for g in known if group_counts.get(g, 0) > 0]


# ── 핵심 함수 ──────────────────────────────────────────────

def resolve_contrasts(
    question: str,
    sample_manifest: dict | None = None,
    data_flags: dict | None = None,
    plan_id: str | None = None,
    output_dir: str | Path | None = None,
) -> ResearchPlan:
    """질문과 데이터 가용성으로 ContrastSpec 목록을 결정하고 ResearchPlan을 반환.

    Args:
        question: 사용자 질문 문자열.
        sample_manifest: 그룹별 샘플 수 정보.
            {"groups": {"TNBC": 120, "Luminal_A": 200, ...}} 또는
            {"samples": [...]} 형식.
        data_flags: 추가 가용성 플래그.
            {"has_normal_tissue": bool, "paired": bool, "has_batch_info": bool, ...}
        plan_id: 계획 ID (없으면 타임스탬프 기반 자동 생성).
        output_dir: 결과 파일 저장 경로.

    Returns:
        ResearchPlan (contrast 목록 포함).
    """
    import re

    if plan_id is None:
        plan_id = f"plan_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"

    flags = data_flags or {}
    group_counts = _count_groups(sample_manifest)

    # ── TNBC 유무 확인 ───────────────────────────────────
    n_tnbc = group_counts.get("TNBC", 0)
    non_tnbc_groups = _non_tnbc_groups(group_counts)
    n_non_tnbc = sum(group_counts.get(g, 0) for g in non_tnbc_groups)
    has_normal = flags.get("has_normal_tissue", False)
    n_normal = group_counts.get("Normal", 0)

    contrasts: list[ContrastSpec] = []

    # ── Contrast 1: TNBC vs Non-TNBC ─────────────────────
    if sample_manifest is None:
        c1_status: ContrastStatus = "PARTIAL"
        c1_reason = "샘플 메타데이터 확인 대기 — R 분석 입력 확인 후 비교군 확정"
    elif n_tnbc == 0:
        c1_status: ContrastStatus = "UNAVAILABLE"
        c1_reason = f"TNBC 샘플 없음 (group_counts={group_counts})"
    elif n_non_tnbc == 0:
        c1_status = "UNAVAILABLE"
        c1_reason = f"Non-TNBC 비교군 없음 (확인된 그룹: {list(group_counts.keys())})"
    elif n_tnbc < 10 or n_non_tnbc < 10:
        c1_status = "UNSUPPORTED"
        c1_reason = f"최소 샘플 수 미달 (TNBC={n_tnbc}, Non-TNBC={n_non_tnbc}, 최소 10 필요)"
    else:
        c1_status = "READY"
        c1_reason = (
            f"TNBC={n_tnbc}, Non-TNBC={n_non_tnbc} "
            f"({', '.join(non_tnbc_groups) or '그룹 미확인'})"
        )

    contrasts.append(ContrastSpec(
        contrast_id="tnbc_vs_nontnbc",
        name="TNBC_vs_NonTNBC",
        positive_group="TNBC",
        reference_group="Non-TNBC",
        status=c1_status,
        reason=c1_reason,
        min_n_positive=10,
        min_n_reference=10,
        actual_n_positive=n_tnbc if n_tnbc else None,
        actual_n_reference=n_non_tnbc if n_non_tnbc else None,
        positive_direction="TNBC > Non-TNBC 방향이면 log2FC 양수",
    ))

    # ── Contrast 2: TNBC vs Normal ───────────────────────
    # 질문에 'normal' / '정상' 언급 여부
    wants_normal = bool(re.search(r"normal|정상|adjacent", question, re.IGNORECASE))

    if not wants_normal:
        # 요청하지 않았으면 스킵 — PARTIAL 처리
        c2_status: ContrastStatus = "PARTIAL"
        c2_reason = "질문에 정상조직 비교 요청 없음 — 향후 contrast 추가 시 재개"
    elif not has_normal and n_normal == 0:
        c2_status = "UNAVAILABLE"
        c2_reason = (
            "정상조직 샘플 없음. "
            "TCGA-BRCA adjacent normal은 별도 확인 필요; "
            "GTEx를 단순 결합해 질환 차이로 해석하지 않는다."
        )
    elif n_normal < 5:
        c2_status = "UNSUPPORTED"
        c2_reason = f"정상 샘플 수 부족 (n={n_normal}, 최소 5 필요)"
    else:
        c2_status = "READY"
        c2_reason = f"TNBC={n_tnbc}, Normal={n_normal}"

    contrasts.append(ContrastSpec(
        contrast_id="tnbc_vs_normal",
        name="TNBC_vs_Normal",
        positive_group="TNBC",
        reference_group="Normal",
        status=c2_status,
        reason=c2_reason,
        min_n_positive=5,
        min_n_reference=5,
        actual_n_positive=n_tnbc if n_tnbc else None,
        actual_n_reference=n_normal if n_normal else None,
        positive_direction="TNBC > Normal 방향이면 log2FC 양수",
    ))

    # ── 질환/아형 추출 (간단한 패턴 매칭) ────────────────
    disease = "breast_cancer"
    subtype = "TNBC"
    if re.search(r"lung|폐암", question, re.IGNORECASE):
        disease = "lung_cancer"
        subtype = None
    elif re.search(r"colorectal|대장", question, re.IGNORECASE):
        disease = "colorectal_cancer"
        subtype = None

    cohort_desc = "TCGA-BRCA RNA-Seq"
    if sample_manifest:
        total = sum(group_counts.values())
        cohort_desc = f"TCGA-BRCA RNA-Seq ({total} samples: " + ", ".join(
            f"{g}={n}" for g, n in sorted(group_counts.items())
        ) + ")"

    plan = ResearchPlan(
        plan_id=plan_id,
        disease=disease,
        subtype=subtype,
        cohort_description=cohort_desc,
        contrasts=contrasts,
        analysis_steps=DEFAULT_ANALYSIS_STEPS,
        budget=DEFAULT_BUDGET,
        stop_conditions=DEFAULT_STOP_CONDITIONS,
    )

    # 결과 저장
    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        plan.save(out / "research_plan.json")
        _write_contrast_manifest(contrasts, out / "contrast_manifest.tsv")

    return plan


def _write_contrast_manifest(contrasts: list[ContrastSpec], path: Path) -> None:
    """contrast_manifest.tsv 작성."""
    header = [
        "contrast_id", "name", "positive_group", "reference_group",
        "status", "actual_n_positive", "actual_n_reference", "reason",
    ]
    rows = ["\t".join(header)]
    for c in contrasts:
        rows.append("\t".join([
            c.contrast_id,
            c.name,
            c.positive_group,
            c.reference_group,
            c.status,
            str(c.actual_n_positive or ""),
            str(c.actual_n_reference or ""),
            c.reason.replace("\t", " "),
        ]))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def update_contrast_result(
    plan: ResearchPlan,
    contrast_id: str,
    result_status: str,
) -> ResearchPlan:
    """contrast 실행 결과 상태를 갱신한 새 ResearchPlan을 반환.

    result_status: dual_contrast_supported | partial_evidence | discordant
    """
    new_contrasts = []
    for c in plan.contrasts:
        if c.contrast_id == contrast_id:
            # frozen dataclass → 새 인스턴스
            c = ContrastSpec(
                contrast_id=c.contrast_id,
                name=c.name,
                positive_group=c.positive_group,
                reference_group=c.reference_group,
                status=c.status,
                reason=c.reason,
                min_n_positive=c.min_n_positive,
                min_n_reference=c.min_n_reference,
                actual_n_positive=c.actual_n_positive,
                actual_n_reference=c.actual_n_reference,
                positive_direction=c.positive_direction,
                result_status=result_status,
            )
        new_contrasts.append(c)

    return ResearchPlan(
        plan_id=plan.plan_id,
        disease=plan.disease,
        subtype=plan.subtype,
        cohort_description=plan.cohort_description,
        contrasts=new_contrasts,
        analysis_steps=plan.analysis_steps,
        budget=plan.budget,
        stop_conditions=plan.stop_conditions,
        created_at=plan.created_at,
        ruleset_version=plan.ruleset_version,
        execution_status=plan.execution_status,
    )


def dual_contrast_summary(plan: ResearchPlan, gene: str, results: dict) -> str:
    """두 비교의 결과를 요약.

    results: {contrast_id: {"log2FC": float, "padj": float, "direction": str}}
    Returns: "dual_contrast_supported" | "partial_evidence" | "discordant"
    """
    ready = {c.contrast_id: c for c in plan.ready_contrasts()}
    if not ready:
        return "partial_evidence"

    directions = {}
    for cid in ready:
        r = results.get(cid, {})
        if r.get("padj") is None or r.get("padj", 1.0) >= 0.05:
            continue
        lfc = r.get("log2FC") or r.get("log2FoldChange") or 0.0
        directions[cid] = "up" if lfc > 0 else "down"

    if len(directions) == len(ready) and len(set(directions.values())) == 1:
        return "dual_contrast_supported"
    if len(directions) > 0:
        if len(set(directions.values())) > 1:
            return "discordant"
        return "partial_evidence"
    return "partial_evidence"
