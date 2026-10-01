"""
② Cell Context Enhanced — CELLxGENE REST API + HPA single-cell 강화
Census SDK가 Windows에서 빌드 불가하므로 REST API로 대체하고
HPA single-cell 데이터로 genesorteR와 유사한 specificity 계산을 수행한다.

출력: cell_context.tsv (gene × cell_class × specificity)
"""

import os
import sys
import csv
import zipfile
import io
from collections import defaultdict

try:
    import requests
except ImportError:
    requests = None

DATASET_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "dataset")
HPA_SC_PATH = os.path.join(DATASET_DIR, "hpa", "rna_single_cell_type.tsv.zip")
PROCESSED_DIR = os.path.join(DATASET_DIR, "tcga_brca", "processed")

# 세포 분류
CELL_CATEGORIES = {
    "TUMOR_EPITHELIAL": {
        "breast hormone-responsive cells", "breast secretory cells",
        "breast myoepithelial cells", "breast lactating cells",
    },
    "IMMUNE": {
        "t-cells", "b-cells", "nk-cells", "macrophages", "monocytes",
        "neutrophils", "mast cells", "plasma cells", "pdcs", "cdc",
    },
    "STROMA": {
        "fibroblasts", "pericytes", "smooth muscle cells",
        "vascular smooth muscle cells", "adipocytes",
    },
    "ENDOTHELIAL": {
        "vascular endothelial cells", "lymphatic endothelial cells",
    },
}


def load_hpa_singlecell(genes: set) -> dict:
    """HPA single-cell 데이터 로드 → {gene: {cell_type: nCPM}}"""
    if not os.path.exists(HPA_SC_PATH):
        return {}
    result = defaultdict(dict)
    with zipfile.ZipFile(HPA_SC_PATH, "r") as zf:
        with zf.open(zf.namelist()[0]) as f:
            reader = csv.DictReader(io.TextIOWrapper(f, encoding="utf-8"), delimiter="\t")
            for row in reader:
                gene = row.get("Gene name", "")
                if gene not in genes:
                    continue
                ct = row.get("Cell type", "")
                try:
                    ncpm = float(row.get("nCPM", 0))
                except ValueError:
                    ncpm = 0
                result[gene][ct] = ncpm
    return dict(result)


def classify_cell_type(ct: str) -> str:
    ct_lower = ct.lower().strip()
    for category, members in CELL_CATEGORIES.items():
        if ct_lower in members:
            return category
    return "OTHER"


def compute_specificity(gene_expr: dict) -> dict:
    """genesorteR-like specificity 계산.
    specificity = (해당 세포 발현 - 다른 세포 평균) / (해당 세포 발현 + 다른 세포 평균 + epsilon)

    Returns: {category: {specificity, ncpm, expressing_fraction_proxy, n_cell_types}}
    """
    # 카테고리별 총 발현
    cat_sums = defaultdict(float)
    cat_counts = defaultdict(int)
    for ct, ncpm in gene_expr.items():
        cat = classify_cell_type(ct)
        cat_sums[cat] += ncpm
        cat_counts[cat] += 1

    total_expr = sum(cat_sums.values())
    if total_expr == 0:
        return {}

    result = {}
    for cat in cat_sums:
        this_expr = cat_sums[cat]
        other_expr = total_expr - this_expr
        n_other = sum(cat_counts[c] for c in cat_counts if c != cat)
        other_mean = other_expr / max(n_other, 1)
        epsilon = 0.01

        specificity = (this_expr - other_mean) / (this_expr + other_mean + epsilon)
        fraction_proxy = cat_counts[cat] / max(sum(cat_counts.values()), 1)

        result[cat] = {
            "specificity": round(specificity, 4),
            "total_ncpm": round(this_expr, 2),
            "n_cell_types": cat_counts[cat],
            "fraction_of_categories": round(fraction_proxy, 3),
        }

    return result


def cellxgene_breast_datasets() -> dict:
    """CELLxGENE REST API로 breast tissue 데이터셋 정보 조회."""
    if not requests:
        return {"error": "requests not installed"}
    try:
        resp = requests.get(
            "https://api.cellxgene.cziscience.com/dp/v1/datasets",
            params={"tissue": "breast", "organism": "Homo sapiens"},
            timeout=15,
        )
        if resp.status_code == 200:
            data = resp.json()
            n = len(data) if isinstance(data, list) else 1
            return {"n_datasets": n, "source": "CELLxGENE REST API"}
        return {"n_datasets": 0, "error": f"HTTP {resp.status_code}"}
    except Exception as e:
        return {"n_datasets": 0, "error": str(e)}


def run_cell_context(gene_names: list, output_dir: str = None) -> dict:
    """전체 cell context 분석 실행.

    Returns:
        {gene: {top_category, specificity, status, categories: {cat: scores}}}
    """
    if output_dir is None:
        output_dir = os.path.join(PROCESSED_DIR, "cell_context")
    os.makedirs(output_dir, exist_ok=True)

    print(f"[CellContext] {len(gene_names)} genes 분석 시작")

    # CELLxGENE dataset 확인
    census_info = cellxgene_breast_datasets()
    print(f"[CellContext] CELLxGENE breast: {census_info}")

    # HPA single-cell 로드
    print(f"[CellContext] HPA single-cell 로드 중...")
    sc_data = load_hpa_singlecell(set(gene_names))
    print(f"[CellContext] {len(sc_data)}/{len(gene_names)} genes에 데이터")

    # 유전자별 specificity 계산
    results = {}
    rows = []
    for gene in gene_names:
        expr = sc_data.get(gene, {})
        if not expr:
            results[gene] = {"status": "INSUFFICIENT", "top_category": None, "categories": {}}
            continue

        spec = compute_specificity(expr)
        if not spec:
            results[gene] = {"status": "INSUFFICIENT", "top_category": None, "categories": {}}
            continue

        # 최고 specificity 카테고리
        top_cat = max(spec, key=lambda c: spec[c]["specificity"])
        top_spec = spec[top_cat]["specificity"]

        # 두 번째 카테고리
        sorted_cats = sorted(spec.items(), key=lambda x: -x[1]["specificity"])
        second_spec = sorted_cats[1][1]["specificity"] if len(sorted_cats) > 1 else -1
        gap = top_spec - second_spec

        results[gene] = {
            "top_category": top_cat,
            "top_specificity": top_spec,
            "gap": round(gap, 4),
            "categories": spec,
        }

        # TSV 행 추가
        for cat, scores in spec.items():
            rows.append({
                "gene": gene,
                "cell_class": cat,
                "specificity": scores["specificity"],
                "total_ncpm": scores["total_ncpm"],
                "n_cell_types": scores["n_cell_types"],
            })

    # cell_context.tsv 저장
    tsv_path = os.path.join(output_dir, "cell_context.tsv")
    if rows:
        with open(tsv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["gene", "cell_class", "specificity",
                                                    "total_ncpm", "n_cell_types"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    print(f"[CellContext] Saved: {tsv_path} ({len(rows)} rows)")

    # 요약
    cat_counts = defaultdict(int)
    for r in results.values():
        tc = r.get("top_category")
        if tc:
            cat_counts[tc] += 1
    print(f"[CellContext] Top categories: {dict(cat_counts)}")

    return results


if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")

    # dual_contrast_supported 유전자 상위 50개로 테스트
    import pandas as pd
    dual_path = os.path.join(PROCESSED_DIR, "deg_results", "dual_contrast_classification.csv")
    if os.path.exists(dual_path):
        df = pd.read_csv(dual_path, index_col=0)
        dual_genes = df[df["dual_status"] == "dual_contrast_supported"]
        if "gene_name" in dual_genes.columns:
            test_genes = dual_genes["gene_name"].dropna().head(50).tolist()
        else:
            test_genes = dual_genes.index[:50].tolist()
    else:
        test_genes = ["ESR1", "PARP1", "GATA3", "TACSTD2", "RRM2"]

    print(f"Testing with {len(test_genes)} genes")
    results = run_cell_context(test_genes)
    print(f"\nSample results:")
    for gene in list(results.keys())[:5]:
        r = results[gene]
        print(f"  {gene}: top={r.get('top_category')} spec={r.get('top_specificity', 'N/A')} gap={r.get('gap', 'N/A')}")
