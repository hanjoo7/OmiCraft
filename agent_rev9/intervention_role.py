"""
⑤ Intervention Role + Therapeutic Direction
워크플로우 06: 발현 부호에서 직접 결정하지 않는다

intervention_role: FUNCTION_MODULATION / CELL_TARGETING / BIOMARKER_ONLY
therapeutic_direction: Inhibit / Degrade / Neutralize / Silence / Activate / NOT_APPLICABLE / UNRESOLVED
"""

from typing import Literal

InterventionRole = Literal["FUNCTION_MODULATION", "CELL_TARGETING", "BIOMARKER_ONLY"]
TherapeuticDirection = Literal["Inhibit", "Degrade", "Neutralize", "Silence", "Activate", "NOT_APPLICABLE", "UNRESOLVED"]


def assign_intervention(modality: str, location: str, deg_direction: str,
                        has_functional_evidence: bool = False) -> dict:
    """Intervention role과 therapeutic direction 결정.

    Args:
        modality: ADC, SMALL_MOLECULE, DEGRADER, RNA_THERAPEUTIC, etc.
        location: SURFACE, INTRACELLULAR, NUCLEAR, SECRETED
        deg_direction: Up_in_TNBC, Down_in_TNBC
        has_functional_evidence: 기능 실험 근거 유무

    Returns:
        {intervention_role, therapeutic_direction, intended_effect, reason}
    """
    # ADC/CAR-T → CELL_TARGETING (기능 억제가 아닌 세포 제거)
    if modality in ("ADC", "CAR_T", "BISPECIFIC"):
        return {
            "intervention_role": "CELL_TARGETING",
            "therapeutic_direction": "NOT_APPLICABLE",
            "intended_effect": "tumor_cell_elimination",
            "reason": f"{modality}는 항원 기능 억제가 아닌 세포 표적화 · 제거",
        }

    # BIOMARKER → BIOMARKER_ONLY
    if modality == "BIOMARKER":
        return {
            "intervention_role": "BIOMARKER_ONLY",
            "therapeutic_direction": "NOT_APPLICABLE",
            "intended_effect": "patient_stratification",
            "reason": "진단/선별 용도, 치료 modality 아님",
        }

    # FUNCTION_MODULATION — 방향은 근거 기반
    direction = "UNRESOLVED"
    reason = ""

    if modality == "SMALL_MOLECULE" and location in ("INTRACELLULAR", "NUCLEAR"):
        if has_functional_evidence:
            direction = "Inhibit"
            reason = "세포 내 기능 억제 (perturbation 근거 확인)"
        else:
            direction = "UNRESOLVED"
            reason = "기능 실험 근거 미확인 — 방향 보류"

    elif modality == "DEGRADER":
        direction = "Degrade"
        reason = "단백질 전체 제거 (scaffold 포함)"

    elif modality == "NEUTRALIZING_AB":
        direction = "Neutralize"
        reason = "분비 ligand/외부 신호 중화"

    elif modality in ("RNA_THERAPEUTIC", "SIRNA", "ASO"):
        if deg_direction == "Up_in_TNBC":
            direction = "Silence"
            reason = "TNBC 과발현 transcript 감소"
        else:
            direction = "UNRESOLVED"
            reason = "하향 발현 유전자에 silencing은 부적절 — 방향 재검토 필요"

    elif modality == "DE_NOVO_BINDER":
        if location == "SURFACE":
            direction = "Inhibit"
            reason = "표면 결합 → 기능 차단"
        else:
            direction = "UNRESOLVED"
            reason = "binder 결합 후 기대 효과 미확인"

    return {
        "intervention_role": "FUNCTION_MODULATION",
        "therapeutic_direction": direction,
        "intended_effect": f"{direction.lower()}_target" if direction != "UNRESOLVED" else "unknown",
        "reason": reason,
    }


def batch_assign(targets: list) -> list:
    """여러 표적에 대해 intervention 배정."""
    results = []
    for t in targets:
        result = assign_intervention(
            modality=t.get("modality", ""),
            location=t.get("location", "UNKNOWN"),
            deg_direction=t.get("direction", ""),
            has_functional_evidence=t.get("has_functional_evidence", False),
        )
        result["gene"] = t.get("gene_name", "")
        result["modality"] = t.get("modality", "")
        results.append(result)

    # 요약
    roles = {}
    for r in results:
        role = r["intervention_role"]
        roles[role] = roles.get(role, 0) + 1
    print(f"[Intervention] {len(results)} targets: {roles}")

    return results


if __name__ == "__main__":
    import sys
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")

    test = [
        {"gene_name": "ESR1", "modality": "ADC", "location": "SURFACE", "direction": "Down_in_TNBC"},
        {"gene_name": "PARP1", "modality": "SMALL_MOLECULE", "location": "NUCLEAR", "direction": "Up_in_TNBC", "has_functional_evidence": True},
        {"gene_name": "RRM2", "modality": "SMALL_MOLECULE", "location": "NUCLEAR", "direction": "Up_in_TNBC"},
        {"gene_name": "GATA3", "modality": "BIOMARKER", "location": "NUCLEAR", "direction": "Down_in_TNBC"},
        {"gene_name": "TRPM8", "modality": "RNA_THERAPEUTIC", "location": "SURFACE", "direction": "Up_in_TNBC"},
        {"gene_name": "VEGFA", "modality": "NEUTRALIZING_AB", "location": "SECRETED", "direction": "Up_in_TNBC"},
    ]

    results = batch_assign(test)
    print("\nIntervention Results:")
    for r in results:
        print(f"  {r['gene']:12s} {r['modality']:15s} → role={r['intervention_role']:20s} "
              f"dir={r['therapeutic_direction']:15s} | {r['reason'][:50]}")
