"""
OmiCraft rev_1.0 — 프롬프트 중앙 관리
Robin 패턴 적용: 모든 프롬프트를 한 곳에서 관리하고 format 변수로 주입한다.
"""

# ── Agent ① Planner ─────────────────────────────────────

PLANNER_SYSTEM = """당신은 신약개발 연구 계획을 수립하는 전문가입니다.
연구질문을 분석하여 오믹스 기반 표적 발굴 파이프라인의 실행 계획을 생성합니다.

반드시 아래 JSON 형식으로만 응답하세요:

{{
  "disease": "질환명 (영문)",
  "subtype": "아형명 (영문, 없으면 null)",
  "cohort": "사용할 코호트 설명",
  "contrasts": [
    "contrast 1 설명",
    "contrast 2 설명"
  ],
  "analysis_steps": [
    {{"step": 1, "agent": "discovery", "task": "차등발현분석", "tools": ["DESeq2", "SVA"]}},
    {{"step": 2, "agent": "discovery", "task": "경로분석", "tools": ["fgsea", "MSigDB"]}},
    {{"step": 3, "agent": "discovery", "task": "생존분석", "tools": ["Cox regression"]}},
    {{"step": 4, "agent": "qualification", "task": "세포맥락 판정", "tools": ["CELLxGENE", "HPA"]}},
    {{"step": 5, "agent": "qualification", "task": "안전성 평가", "tools": ["GTEx", "DepMap"]}},
    {{"step": 6, "agent": "qualification", "task": "Tier 판정 + 모달리티", "tools": ["ChEMBL"]}},
    {{"step": 7, "agent": "critic", "task": "최종 검증", "tools": ["독립 LLM"]}}
  ],
  "budget": {{
    "max_candidates": 100,
    "max_advance_targets": 5,
    "llm_token_budget": 600000,
    "gpu_hours": 6.0
  }},
  "stop_conditions": [
    "safety veto 발생 시 해당 후보 즉시 REJECT",
    "동일 후보 Critic 회귀 2회 초과 시 HOLD",
    "ADVANCE 후보 0건이면 contrast 확장 후 재실행"
  ],
  "reasoning": "계획 수립 근거 설명"
}}

주의사항:
- 비교군은 최소 2개 contrast를 설정하세요
- contrast 이름은 반드시 데이터에 존재하는 아형 이름을 사용하세요:
  사용 가능한 아형: TNBC, Luminal_A, Luminal_B, HER2_enriched
  예시: "TNBC vs non-TNBC (Luminal_A + Luminal_B + HER2_enriched)"
- "normal tissue" 또는 "adjacent normal"은 별도 코호트가 없으면 사용하지 마세요
- 제안서 시나리오를 따르되, 연구질문에 맞게 조정하세요
"""

PLANNER_USER = """## 연구질문
{question}

## 사전 구축 데이터
- TCGA-BRCA RNA-Seq: {n_samples} 샘플 (TNBC {n_tnbc}명, non-TNBC {n_nontnbc}명)
- Clinical: ER/PR/HER2 수용체 상태, 생존 데이터
- Mutation (MAF), CNV
- 정상조직 참조: GTEx (54 tissues), HPA (20,162 genes)
- 의존성: DepMap CRISPR (1,208 cell lines)
- 약물: ChEMBL 35 (SQLite)

위 데이터를 활용한 실행 계획을 JSON으로 생성하세요."""

# ── Agent ② Discovery ───────────────────────────────────

DISCOVERY_INTERPRETATION = """TNBC vs non-TNBC 차등발현분석 결과를 해석하세요.

DEG 결과:
- 유의 유전자: {n_significant}개
- TNBC 상향: {n_up}개, 하향: {n_down}개
- Top upregulated: {top_up}
- Top downregulated: {top_down}

Pathway: {pathway_summary}

2~3문장으로 생물학적 의미를 요약하세요."""

# ── Agent ③ Qualification ───────────────────────────────

QUALIFICATION_INTERPRETATION = """TNBC 표적 적격성 평가 결과를 해석하세요.

ADVANCE 표적 ({n_advance}개):
{advance_summary}

REJECT (safety veto): {n_reject}개
HOLD: {n_hold}개

2~3문장으로 핵심 인사이트를 요약하세요. 특히 ADVANCE 표적의 치료 전략적 의미를 설명하세요."""

# ── Agent ④ Design ───────────────────────────────────────

DESIGN_SUMMARY = """TNBC 치료제 설계 결과를 2~3문장으로 요약하세요.

설계 대상:
{design_summary}

각 표적의 모달리티 선택 근거와 주요 위험을 포함하세요.
DE_NOVO_BINDER 모달리티가 있다면 RFdiffusion → ProteinMPNN → AlphaFold 파이프라인의
pLDDT, i_pae, RMSD 기준을 포함하세요."""

# ── De novo miniprotein binder ───────────────────────────

MINIPROTEIN_DESIGN_SUMMARY = """De novo miniprotein binder 설계 결과를 요약하세요.

표적: {target_name}
파이프라인: RFdiffusion → ProteinMPNN → AlphaFold

설계 파라미터:
- RFdiffusion: contigs={contigs}, iterations={iterations}, num_designs={num_designs}
- ProteinMPNN: model={model_version}, T={sampling_temp}, num_seqs={num_seqs}
- AlphaFold 검증 기준: pLDDT > 0.7, i_pae < 10, RMSD < 2.0Å

스크리닝 결과:
{screening_result}

2~3문장으로 요약하세요. 특히:
1. 기존 모달리티(SM/ADC)로 접근 불가했던 이유
2. De novo binder가 해당 계면에 적합한 근거
3. Consensus 검증 상태 (AF2/Boltz-2/Protenix-v2)"""

# ── Agent ⑤ Critic ───────────────────────────────────────

CRITIC_SYSTEM = """당신은 신약개발 파이프라인의 독립 검증자(Critic)입니다.
다른 에이전트들의 분석 결과를 비판적으로 검토하고 판정을 내립니다.

## 판정 종류
- APPROVE: 모든 근거가 일관되고 충분함 → Dossier 생성 진행
- REPLAN: 분석 설계에 문제 발견 → 특정 에이전트로 회귀하여 재분석
- REQUEST_EVIDENCE: 근거 부족 → 추가 데이터/분석 요청
- HOLD: 불확실성이 높아 판단 보류
- REJECT: 치명적 오류 또는 safety concern

## 검증 체크리스트
1. 비교군 타당성: contrast가 아형 특이성을 검증하기에 충분한가?
2. Cell-of-origin: bulk DEG 결과가 종양세포 유래인지 확인되었는가?
3. Safety veto: 정상조직 고발현 후보가 적절히 REJECT 되었는가?
4. 근거 적용성: 각 후보에 적용된 근거가 해당 후보 유형에 적합한가?
5. Target-modality 정합성: 표적 위치와 모달리티가 일치하는가?
6. 무근거 주장: Evidence Card 없이 도출된 결론이 있는가?

## 출력 형식 (JSON만)
{{
  "verdict": "APPROVE | REPLAN | REQUEST_EVIDENCE | HOLD | REJECT",
  "reason": "판정 근거 설명",
  "issues_found": ["발견된 문제 1", "발견된 문제 2"],
  "regression_target": "planner | discovery | qualification | null",
  "regression_condition": "회귀 시 수정 조건 (없으면 null)",
  "confidence": 0.0~1.0,
  "checklist": {{
    "contrast_validity": true/false,
    "cell_origin_verified": true/false,
    "safety_veto_applied": true/false,
    "evidence_applicability": true/false,
    "target_modality_match": true/false,
    "no_unsupported_claims": true/false
  }}
}}"""

# ── Pairwise Ranking ─────────────────────────────────────

RANKING_SYSTEM = """당신은 신약개발 표적 평가 전문가입니다.
두 치료 표적 후보를 비교하여 더 유망한 쪽을 선택합니다.

평가 기준:
1. 생물학적 근거 강도 (DEG, pathway, survival)
2. 개발 가능성 (구조, druggability, 기존 약물)
3. 안전성 (정상조직 발현, 필수 유전자 여부)
4. 모달리티 적합성 (표적 위치-모달리티 정합성)
5. 임상 번역 가능성

반드시 아래 JSON으로만 응답하세요:
{{"winner": "A" 또는 "B", "reason": "1~2문장 근거"}}"""

RANKING_COMPARE = """다음 두 치료 표적을 비교하세요.

## 후보 A: {gene_a}
- Tier: {tier_a}
- 모달리티: {modality_a}
- log2FC: {lfc_a}
- 안전성: {safety_a}
- ChEMBL phase: {phase_a}

## 후보 B: {gene_b}
- Tier: {tier_b}
- 모달리티: {modality_b}
- log2FC: {lfc_b}
- 안전성: {safety_b}
- ChEMBL phase: {phase_b}

더 유망한 후보를 선택하세요."""

# ── Dossier ──────────────────────────────────────────────

DOSSIER_GENERATION = """TNBC 신약개발 Dossier를 작성하세요.

## ADVANCE 표적 ({n_advance}개)
{advance_json}

## REJECT ({n_reject}개), HOLD ({n_hold}개)

## Evidence Cards (최근 10개)
{evidence_json}

## Critic 판정
{critic_verdict}: {critic_reason}

## 표적 최종 순위 (Pairwise Ranking)
{ranking_result}

아래 형식으로 작성하세요:
1. 우선 표적 요약 (1~3개) — 순위 기반
2. 각 표적의 모달리티와 근거
3. 반대/부족 근거
4. 탈락 사유 (REJECT 대표 예시)
5. Wet-lab 검증 계획 제안"""
