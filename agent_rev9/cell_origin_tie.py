"""
③ Cell-of-origin tie 처리
워크플로우 05-A: SINGLE_SUPPORTED / MULTI_SUPPORTED / UNRESOLVED / INSUFFICIENT

gap = (s1 - s2) / max(s1, epsilon)
gap > 0.10 + bootstrap 안정 → SINGLE_SUPPORTED
두 카테고리 모두 안정 발현 → MULTI_SUPPORTED
gap <= 0.10 또는 불안정 → UNRESOLVED
데이터 부족 → INSUFFICIENT
"""

from typing import Literal

CellOriginStatus = Literal["SINGLE_SUPPORTED", "MULTI_SUPPORTED", "UNRESOLVED", "INSUFFICIENT"]

TIE_GAP_THRESHOLD = 0.10  # gap 이하면 tie review 대상


def classify_cell_origin(cell_context: dict) -> dict:
    """Cell-of-origin tie 분류.

    Args:
        cell_context: cell_context_enhanced.run_cell_context() 결과의 단일 유전자

    Returns:
        {status, top_category, origins: list, gap, reason}
    """
    if not cell_context or cell_context.get("status") == "INSUFFICIENT":
        return {
            "status": "INSUFFICIENT",
            "top_category": None,
            "origins": [],
            "gap": None,
            "reason": "HPA single-cell 데이터 부족 또는 미측정",
        }

    categories = cell_context.get("categories", {})
    if not categories:
        return {"status": "INSUFFICIENT", "top_category": None, "origins": [], "gap": None,
                "reason": "카테고리별 발현 정보 없음"}

    # specificity 순 정렬
    sorted_cats = sorted(categories.items(), key=lambda x: -x[1].get("specificity", 0))
    top_cat, top_scores = sorted_cats[0]
    top_spec = top_scores.get("specificity", 0)

    if len(sorted_cats) < 2:
        return {
            "status": "SINGLE_SUPPORTED",
            "top_category": top_cat,
            "origins": [top_cat],
            "gap": 1.0,
            "reason": f"단일 카테고리 ({top_cat})",
        }

    second_cat, second_scores = sorted_cats[1]
    second_spec = second_scores.get("specificity", 0)

    epsilon = 0.01
    gap = (top_spec - second_spec) / max(top_spec, epsilon)

    # 분류 규칙
    if gap > TIE_GAP_THRESHOLD and top_spec > 0:
        # 명확한 우세
        return {
            "status": "SINGLE_SUPPORTED",
            "top_category": top_cat,
            "origins": [top_cat],
            "gap": round(gap, 4),
            "reason": f"{top_cat} 우세 (gap={gap:.3f} > {TIE_GAP_THRESHOLD})",
        }

    # 두 카테고리 모두 양의 specificity
    if top_spec > 0 and second_spec > 0:
        return {
            "status": "MULTI_SUPPORTED",
            "top_category": top_cat,
            "origins": [top_cat, second_cat],
            "gap": round(gap, 4),
            "reason": f"{top_cat}+{second_cat} 복수 발현 (gap={gap:.3f} <= {TIE_GAP_THRESHOLD})",
        }

    # gap이 작고 불안정
    return {
        "status": "UNRESOLVED",
        "top_category": top_cat,
        "origins": [top_cat, second_cat],
        "gap": round(gap, 4),
        "reason": f"tie 미해결: {top_cat} vs {second_cat} (gap={gap:.3f})",
    }


def batch_classify(cell_contexts: dict) -> dict:
    """여러 유전자에 대해 tie 분류."""
    results = {}
    status_counts = {"SINGLE_SUPPORTED": 0, "MULTI_SUPPORTED": 0, "UNRESOLVED": 0, "INSUFFICIENT": 0}

    for gene, ctx in cell_contexts.items():
        result = classify_cell_origin(ctx)
        results[gene] = result
        status_counts[result["status"]] += 1

    print(f"[CellOriginTie] {len(results)} genes classified: {status_counts}")
    return results


if __name__ == "__main__":
    import os
    import sys
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

    import pandas as pd

    from .cell_context_enhanced import run_cell_context

    dual_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed",
        "deg_results", "dual_contrast_classification.csv"
    )
    df = pd.read_csv(dual_path, index_col=0)
    dual_genes = df[df["dual_status"] == "dual_contrast_supported"]
    test_genes = dual_genes["gene_name"].dropna().head(30).tolist()

    contexts = run_cell_context(test_genes)
    ties = batch_classify(contexts)

    print("\nSample tie results:")
    for gene in list(ties.keys())[:10]:
        t = ties[gene]
        print(f"  {gene:12s}: {t['status']:20s} origins={t['origins']} gap={t.get('gap','N/A')}")
