"""
B/M/D/Tier 시스템 — 운영 워크플로우 v5 Section 07~08

B: Biological confidence (B0~B3)
M: Modality suitability (M1~M3, N/A)
D: Development readiness (D1~D3)
Tier: 최종 등급 (1A, 1B, 2A, 2B, 3, HOLD, REJECT)

gene × context × intervention hypothesis → B
hypothesis × modality → M, D, Tier
"""

from typing import Literal


BioConfidence = Literal["B0", "B1", "B2", "B3"]
ModalitySuitability = Literal["M1", "M2", "M3", "N/A", "UNKNOWN"]
DevReadiness = Literal["D1", "D2", "D3"]
TierLevel = Literal["1A", "1B", "2A", "2B", "3", "HOLD", "NOT_APPLICABLE", "REJECT"]


def assign_bio_confidence(evidence: dict) -> tuple:
    """Biological confidence 등급 산출.

    evidence keys:
        dual_contrast: str (dual_contrast_supported/partial_evidence/discordant)
        cell_origin: str (SINGLE_SUPPORTED/MULTI_SUPPORTED/UNRESOLVED/INSUFFICIENT)
        depmap_status: str (CONCORDANT/CRISPR_ONLY/DISCORDANT/UNAVAILABLE)
        has_functional_evidence: bool
        has_safety_veto: bool
        has_external_replication: bool
        has_causal_evidence: bool

    Returns: (tier, reasons)
    """
    reasons = []

    # B0: 견고한 반증
    if evidence.get("has_safety_veto"):
        return "B0", ["Critical safety veto — 기각"]
    if evidence.get("dual_contrast") == "discordant":
        reasons.append("두 비교 방향 상충")
        return "B0", reasons

    # B1: 강한 재현 + 세포 맥락 + 기능/인과 + 안전성 지지
    score = 0
    if evidence.get("dual_contrast") == "dual_contrast_supported":
        score += 2
        reasons.append("두 비교 재현")
    elif evidence.get("dual_contrast") == "partial_evidence":
        score += 1
        reasons.append("한 비교만 지지")

    if evidence.get("cell_origin") == "SINGLE_SUPPORTED":
        score += 2
        reasons.append("종양세포 맥락 확인")
    elif evidence.get("cell_origin") == "MULTI_SUPPORTED":
        score += 1
        reasons.append("복수 세포 발현")
    elif evidence.get("cell_origin") == "UNRESOLVED":
        reasons.append("세포 출처 미해결")

    if evidence.get("has_functional_evidence"):
        score += 2
        reasons.append("기능/인과 근거")
    if evidence.get("has_external_replication"):
        score += 1
        reasons.append("외부 재현")
    if evidence.get("depmap_status") == "CONCORDANT":
        score += 1
        reasons.append("DepMap CRISPR+RNAi 일치")
    elif evidence.get("depmap_status") == "CRISPR_ONLY":
        crispr_effect = evidence.get("crispr_effect", 0) or 0
        if crispr_effect < -0.5:
            score += 1
            reasons.append(f"CRISPR 강한 의존성 (median={crispr_effect:.2f})")

    # Association supports confidence without implying experimental causality.
    for key, label in (("ot_genetic_association", "genetic association"),
                       ("ot_somatic_mutation", "somatic mutation")):
        value = evidence.get(key)
        if value is not None and value > 0:
            score += 1
            reasons.append(f"OT {label} ({value:.3f})")

    # cell_origin=INSUFFICIENT는 항상 +0 (scRNA-seq 미보유 시 패널티 없음)
    # B1 임계값: 5 (scRNA-seq 없어도 dual+functional+CRISPR_strong으로 달성 가능)
    if score >= 5:
        return "B1", reasons
    elif score >= 3:
        return "B2", reasons
    else:
        return "B3", reasons


def assign_modality_suitability(target_info: dict, modality: str) -> tuple:
    """Modality suitability (M) 산출.

    target_info: location, has_pocket, has_structure, surface_accessible, etc.
    """
    location = target_info.get("location", "UNKNOWN")
    has_pocket = target_info.get("has_pocket", False)
    has_structure = target_info.get("has_structure", False)
    surface = target_info.get("surface_accessible", False)

    modality_rules = {
        "SMALL_MOLECULE": lambda: ("M1" if has_pocket else "M2" if location == "INTRACELLULAR" else "M3",
                                   "pocket 가용" if has_pocket else "pocket 미확인"),
        "ADC": lambda: ("M1" if surface else "N/A",
                        "표면 접근 가능" if surface else "표면 비접근"),
        "DE_NOVO_BINDER": lambda: ("M1" if has_structure else "M2",
                                    "구조 가용" if has_structure else "구조 필요"),
        "DEGRADER": lambda: ("M1" if has_pocket and location == "INTRACELLULAR" else "M2",
                             "세포 내 + binder"),
        "RNA_THERAPEUTIC": lambda: ("M2", "transcript 서열 확인 필요"),
        "BIOMARKER": lambda: ("M1", "진단 assay"),
    }

    func = modality_rules.get(modality, lambda: ("UNKNOWN", "modality 규칙 없음"))
    m, reason = func()
    return m, reason


def assign_dev_readiness(asset_info: dict) -> tuple:
    """Development readiness (D) 산출."""
    max_phase = asset_info.get("chembl_max_phase", 0)
    n_compounds = asset_info.get("chembl_n_compounds", 0)
    has_structure = asset_info.get("has_structure", False)

    if max_phase >= 4 or (max_phase >= 2 and n_compounds > 50):
        return "D1", f"임상 자산 존재 (phase {max_phase})"
    elif has_structure and n_compounds > 10:
        return "D2", "구조+화합물 있으나 보완 필요"
    elif has_structure:
        return "D2", "구조 가용, ligand 부족"
    else:
        return "D3", "경로 불명확"


def assign_tier(b: BioConfidence, m: ModalitySuitability, d: DevReadiness) -> TierLevel:
    """B + M + D → 최종 Tier."""
    if b == "B0":
        return "REJECT"
    if m == "N/A":
        return "NOT_APPLICABLE"
    if m == "UNKNOWN" or b == "B3":
        return "3"

    tier_map = {
        ("B1", "M1", "D1"): "1A",
        ("B1", "M1", "D2"): "1B",
        ("B1", "M2", "D1"): "2A",
        ("B1", "M2", "D2"): "2B",
        ("B2", "M1", "D1"): "2A",
        ("B2", "M1", "D2"): "2B",
        ("B2", "M2", "D1"): "2B",
        ("B2", "M2", "D2"): "3",
    }
    return tier_map.get((b, m, d), "3")


def evaluate_target_modality(gene: str, bio_evidence: dict,
                              target_info: dict, asset_info: dict,
                              modalities: list) -> list:
    """표적 × 모달리티별 Tier 산출.

    Returns: [{gene, modality, B, M, D, Tier, reasons}]
    """
    b, b_reasons = assign_bio_confidence(bio_evidence)
    results = []

    for mod in modalities:
        m, m_reason = assign_modality_suitability(target_info, mod)
        d, d_reason = assign_dev_readiness(asset_info)
        tier = assign_tier(b, m, d)

        results.append({
            "gene": gene,
            "modality": mod,
            "B": b,
            "M": m,
            "D": d,
            "Tier": tier,
            "bio_reasons": b_reasons,
            "modality_reason": m_reason,
            "readiness_reason": d_reason,
        })

    return results
