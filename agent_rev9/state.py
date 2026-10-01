"""
OmiCraft AgentState — 5개 에이전트가 공유하는 상태 스키마
워크플로우 v5 기반 전체 파이프라인 출력 필드 포함

LangGraph 호환:
  - messages / errors / evidence_cards 는 Annotated[list, operator.add] reducer 사용
    → 각 노드가 반환한 항목이 기존 리스트에 추가(append)됨
  - 나머지 필드는 노드 반환값으로 override
"""

import operator
from typing import Annotated, Literal, TypedDict


class EvidenceCard(TypedDict, total=False):
    """개별 근거 카드 (§4 evidence_cards.jsonl 행 단위)"""

    id: str
    source: str   # "DESeq2", "DepMap", "GTEx", "ChEMBL", "Critic", ...
    type: Literal["SUPPORTIVE", "CONTRADICTORY", "NOT_APPLICABLE", "UNKNOWN"]
    content: str
    gene: str
    confidence: float


class CandidateTarget(TypedDict, total=False):
    """후보 표적 — ADVANCE/HOLD/REJECT 공통 형식"""

    gene_id: str
    gene_name: str
    # str로 선언해 ADC, CAR_T, DEGRADER, NEUTRALIZING_AB 등 모든 modality 허용
    modality: str
    log2fc: float
    fdr: float
    direction: str          # "Up_in_TNBC" | "Down_in_TNBC"
    bio_tier: str           # "B1" | "B2" | "B3" | "B0"
    dev_tier: str           # "D1" | "D2" | "D3"
    M: str                  # "M1" | "M2" | "M3" | "N/A"
    D: str                  # "D1" | "D2" | "D3"
    Tier: str               # "Tier1A" | "Tier1B" | "Tier2A" | "Tier2B" | "Tier3"
    safety_verdict: str     # "PASS" | "CAUTION" | "HOLD" | "VETO"
    judgment: str           # "ADVANCE" | "HOLD" | "REJECT"
    modality_reason: str
    evidence_ids: list

    # cell-of-origin (WF 04-05)
    cell_origin_status: str     # "SINGLE_SUPPORTED" | "MULTI_SUPPORTED" | "UNRESOLVED" | "INSUFFICIENT"
    cell_origins: list          # ["TUMOR_EPITHELIAL"] or ["TUMOR_EPITHELIAL", "IMMUNE"]

    # DepMap (WF 05-B)
    depmap_concordance: str     # "CONCORDANT" | "CRISPR_ONLY" | "DISCORDANT" | ...
    crispr_pct_dependent: float

    # Intervention (WF 06)
    intervention_role: str      # "FUNCTION_MODULATION" | "CELL_TARGETING" | "BIOMARKER_ONLY"
    therapeutic_direction: str  # "Inhibit" | "Degrade" | "Silence" | "Neutralize" | ...


class OmiCraftState(TypedDict, total=False):
    """5-Agent DAG 공유 상태 — 워크플로우 01~09"""

    # ── 연구 계획 (WF 01) ──────────────────────────────────
    research_question: str
    sample_manifest: dict
    disease: str
    subtype: str
    contrasts: list         # list[ContrastSpec.to_dict()]
    contrast_status: dict   # {contrast_name: "READY"|"UNAVAILABLE"|"PARTIAL"|"DISCORDANT"}

    # ResearchPlan (contrast_spec.py)
    research_plan: dict     # ResearchPlan.to_dict() — plan_id, question, contrasts, status
    contrast_manifest: str  # 경로: contrast_manifest.tsv

    # ── Discovery 출력 (WF 02-03) ────────────────────────
    r_analysis_result: dict         # r_tool.run_r_analysis() 전체 출력
    apear_results: dict             # apear_adapter.run_apear_for_contrasts() 출력
    n_candidates: int
    candidate_pool_path: str

    # ── Qualification 출력 (WF 04-07) ───────────────────
    qualification_results: dict     # gene → qualify_gene() 전체 결과 (tier_results 포함)
    advance_targets: list           # list[CandidateTarget] — ADVANCE 판정
    hold_targets: list              # list[CandidateTarget] — HOLD 판정
    reject_targets: list            # list[CandidateTarget] — REJECT 판정
    qualification_results_path: str # 경로: qualification_results.json

    # Cell context (WF 04-05-A) — cellxgene_adapter.py
    cell_contexts: dict             # gene → assign_context_status() 결과
    cell_context_results: dict      # 하위호환: gene → {top_category, specificity, ...}
    cell_origin_ties: dict          # gene → {status, origins, gap}

    # DepMap (WF 05-B) — depmap_module.py
    depmap_results: dict            # gene → evaluate_dependency() 결과

    # Intervention + safety (WF 06)
    interventions: list             # [{gene, intervention_role, therapeutic_direction, ...}]
    safety_results: dict            # gene → {verdict, scope, reason}

    # ── 출력 테이블 경로 (WF §4) ────────────────────────
    gene_overview_path: str         # gene_overview.tsv
    de_evidence_path: str           # de_evidence.tsv
    modality_options_path: str      # modality_options.tsv
    evidence_cards_path: str        # evidence_cards.jsonl
    output_tables: dict             # {table_name: file_path}

    # ── 사용자 선택 gate (WF 08) — user_selection.py ───
    user_selection_events: list     # list[SelectionEvent dict] — create_selection_event() 출력
    selected_route: str             # 선택된 "{gene}:{modality}"
    gate_status: str                # "ALLOWED" | "PENDING" | "HOLD"
    gate_reason: str                # gate 상태 설명
    design_allowed_targets: list    # gate를 통과한 표적 목록
    blocked_targets: list           # gate 차단 표적 목록
    # (하위호환)
    selected_targets: list          # list[CandidateTarget] — 사용자 선택 후 확정
    user_selection_made: bool       # True가 되기 전 고비용 설계 차단

    # ── Evidence Store ───────────────────────────────────
    evidence_cards: Annotated[list, operator.add]   # list[EvidenceCard] — 노드마다 append

    # ── 설계 모드 (WF 09) ────────────────────────────────
    route: Literal["small_molecule", "protein_binder_rfd3", "adc", "degrader"]
    design_mode: Literal["small_molecule", "protein_binder", "adc", "degrader"]
    use_rfd3: bool
    dry_run: bool
    structural_backend: str         # "boltz2" | "alphafold3"
    protein_design_input: dict
    protein_design_result: dict
    screening_result: dict
    screening_options: dict
    design_results_path: str
    screening_results_path: str
    validation_results: list
    candidate_rankings: list
    candidate_summary_path: str
    budget: dict

    # ── Critic 판정 ──────────────────────────────────────
    critic_verdict: Literal["APPROVE", "HOLD", "REJECT", "REPLAN", "REQUEST_EVIDENCE"]
    _run_directory: str
    reuse_results: bool
    reused_stages: Annotated[list, operator.add]
    review_issues: Annotated[list, operator.add]
    critic_result: dict
    dossier_path: str
    dossier_status: str
    design_outcome: dict
    automatic_design_results: list
    critic_reason: str
    critic_feedback: str            # Critic → Planner 피드백 (REPLAN 시 개선 방향)
    critic_regression_target: str   # "planner" | "discovery" | "qualification"
    critic_count: int

    # ── LangGraph 라우팅 제어 ─────────────────────────────
    next_node: str                  # 조건부 엣지가 읽는 라우팅 힌트
    replan_count: int               # Planner 재호출 횟수 (무한루프 방지, max=2)

    # ── 실행 메타 ────────────────────────────────────────
    current_agent: str
    messages: Annotated[list, operator.add]   # 노드마다 append
    errors:   Annotated[list, operator.add]   # 노드마다 append
    run_id: str

    # ── 쌍둥이 환자 대조 (optional) ──────────────────────
    twin_patients: dict             # {pid: {id, name, tissue, subtype, mutation, ...}}
    twin_patient_results: dict      # {pid: {primary_target, modality, tier, rationale, ...}}

    modality_input: dict
    modality_config: dict
    modality_assessments: list
    modality_results: list
    output_dir: str
    execution_status: str
    validation_decision: str
