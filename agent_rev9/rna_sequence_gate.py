"""
⑧ RNA Sequence Gate — 공개 데이터 전용
워크플로우 08: "TCGA STAR count는 gene-level 발현 근거. RNA modality 진입에는
GENCODE/RefSeq의 reference transcript와 isoform/splice-junction 근거를 별도 연결"

공개 데이터에서 확인 가능한 수준:
  gene-level 발현 + GENCODE transcript → 조건부 D2
  isoform-specific → junction 근거 필요, 없으면 D3/HOLD
"""

from typing import Literal

SequenceConfidence = Literal["CONFIRMED", "CONDITIONAL", "UNRESOLVED", "INSUFFICIENT"]


def evaluate_rna_gate(
    gene_symbol: str,
    target_type: str = "gene_whole",  # gene_whole, isoform, fusion, mutant_allele
    has_gencode_transcript: bool = True,
    has_junction_evidence: bool = False,
    has_independent_cohort: bool = False,
) -> dict:
    """RNA sequence gate 평가.

    Returns:
        {eligible, sequence_confidence, readiness, transcript_id, missing_evidence, reason}
    """
    if target_type == "gene_whole":
        if has_gencode_transcript:
            return {
                "eligible": True,
                "sequence_confidence": "CONDITIONAL",
                "readiness": "D2",
                "reason": f"Gene-level 발현 확인, GENCODE transcript 존재 → 공유 exon 대상 조건부 설계",
                "missing_evidence": ["isoform-specific 확인 필요"] if not has_junction_evidence else [],
            }
        return {
            "eligible": False,
            "sequence_confidence": "INSUFFICIENT",
            "readiness": "HOLD",
            "reason": "GENCODE transcript 미확인",
            "missing_evidence": ["GENCODE/RefSeq transcript ID 필요"],
        }

    elif target_type == "isoform":
        if has_junction_evidence and has_independent_cohort:
            return {
                "eligible": True,
                "sequence_confidence": "CONFIRMED",
                "readiness": "D2",
                "reason": "질환 관련 exon/junction + 독립 코호트 재현",
                "missing_evidence": [],
            }
        elif has_junction_evidence:
            return {
                "eligible": True,
                "sequence_confidence": "CONDITIONAL",
                "readiness": "D2",
                "reason": "Junction 근거 있으나 독립 코호트 미확인",
                "missing_evidence": ["독립 코호트 재현 필요"],
            }
        return {
            "eligible": False,
            "sequence_confidence": "UNRESOLVED",
            "readiness": "HOLD",
            "reason": "Gene-level count만으로 isoform-specific 설계 불가",
            "missing_evidence": ["질환 관련 exon/junction 근거", "long-read RNA-seq 또는 RT-PCR"],
        }

    elif target_type == "fusion":
        if has_junction_evidence:
            return {
                "eligible": True,
                "sequence_confidence": "CONDITIONAL",
                "readiness": "D2",
                "reason": "공개 breakpoint/junction 서열 확인",
                "missing_evidence": ["환자 prevalence 확인 필요"] if not has_independent_cohort else [],
            }
        return {
            "eligible": False,
            "sequence_confidence": "INSUFFICIENT",
            "readiness": "HOLD",
            "reason": "정확한 breakpoint/junction 서열 미확인",
            "missing_evidence": ["공개 breakpoint 서열", "환자 prevalence"],
        }

    elif target_type == "mutant_allele":
        return {
            "eligible": False,
            "sequence_confidence": "UNRESOLVED",
            "readiness": "HOLD",
            "reason": "Allele-specific 설계는 variant + phasing 근거 필요 → 현재 HOLD",
            "missing_evidence": ["공개 variant RNA 발현 근거", "allele phasing"],
        }

    return {
        "eligible": False,
        "sequence_confidence": "INSUFFICIENT",
        "readiness": "HOLD",
        "reason": f"지원하지 않는 target_type: {target_type}",
        "missing_evidence": [],
    }


if __name__ == "__main__":
    import sys
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")

    cases = [
        ("TRPM8", "gene_whole", True, False, False),
        ("ESR1", "isoform", False, False, False),
        ("BCR-ABL1", "fusion", False, True, True),
        ("KRAS_G12C", "mutant_allele", True, False, False),
    ]

    print("RNA Sequence Gate:")
    for gene, tt, gc, junc, ind in cases:
        r = evaluate_rna_gate(gene, tt, gc, junc, ind)
        print(f"  {gene:15s} type={tt:15s} → eligible={r['eligible']} "
              f"conf={r['sequence_confidence']:12s} D={r['readiness']} | {r['reason'][:50]}")
