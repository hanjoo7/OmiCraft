"""
Dual Contrast DEG Analysis
TNBC vs Non-TNBC + TNBC vs Normal → 통합 분류

워크플로우 02: "두 비교는 별도 표로 보존한다. 일관된 방향과 사전 설정된
유의성·효과크기 기준을 모두 만족하면 dual_contrast_supported"
"""

import os
import sys
import pandas as pd
import numpy as np

PROCESSED_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed"
)
NORMAL_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "rnaseq_normal"
)


def build_normal_count_matrix():
    """Normal RNA-Seq 파일들을 통합 count matrix로."""
    import glob

    tsv_files = sorted(glob.glob(os.path.join(NORMAL_DIR, "TCGA-*.tsv")))
    if not tsv_files:
        return None, []

    # 첫 파일에서 유전자 목록
    first = pd.read_csv(tsv_files[0], sep="\t", comment="#", skiprows=1)
    first = first[first["gene_type"] == "protein_coding"]
    gene_ids = first["gene_id"].apply(lambda x: x.split(".")[0]).values
    gene_names = first["gene_name"].values

    # 전체 matrix
    sample_ids = []
    all_counts = []

    for fpath in tsv_files:
        basename = os.path.basename(fpath)
        sid = basename.split("_normal_")[0] if "_normal_" in basename else basename.split("_")[0]

        df = pd.read_csv(fpath, sep="\t", comment="#", skiprows=1)
        df = df[df["gene_type"] == "protein_coding"]
        counts = df["unstranded"].values

        if len(counts) == len(gene_ids):
            sample_ids.append(sid + "_N")
            all_counts.append(counts)

    if not all_counts:
        return None, []

    matrix = pd.DataFrame(
        np.array(all_counts).T,
        index=gene_ids,
        columns=sample_ids,
    )
    return matrix, gene_names


def run_dual_contrast_deg(fdr_threshold=0.01, lfc_threshold=1.0):
    """Dual contrast DEG: TNBC vs Non-TNBC + TNBC vs Normal."""

    print("=" * 60)
    print("Dual Contrast DEG Analysis")
    print("=" * 60)

    # 1. 기존 TNBC vs Non-TNBC 결과 로드
    deg1_path = os.path.join(PROCESSED_DIR, "deg_results", "deg_tnbc_vs_nontnbc_full.csv")
    if not os.path.exists(deg1_path):
        print("[ERROR] TNBC vs Non-TNBC DEG 결과 없음")
        return None

    deg1 = pd.read_csv(deg1_path, index_col=0)
    print(f"[Contrast 1] TNBC vs Non-TNBC: {len(deg1)} genes tested")
    sig1 = deg1[(deg1["padj"] < fdr_threshold) & (deg1["log2FoldChange"].abs() > lfc_threshold)]
    print(f"  Significant: {len(sig1)} (Up={sum(sig1['log2FoldChange']>0)}, Down={sum(sig1['log2FoldChange']<0)})")

    # 2. Normal 데이터 확인
    normal_files = len([f for f in os.listdir(NORMAL_DIR) if f.startswith("TCGA-")] if os.path.exists(NORMAL_DIR) else [])
    if normal_files == 0:
        print("\n[Contrast 2] TNBC vs Normal: UNAVAILABLE (Normal 샘플 없음)")
        print("→ 'TNBC 특이성' 주장은 보류 (partial_evidence)")

        # 단일 contrast 결과 저장
        result = _classify_single_contrast(deg1, fdr_threshold, lfc_threshold)
        return result

    # 3. TNBC vs Normal DEG 실행
    print(f"\n[Contrast 2] TNBC vs Normal: {normal_files} Normal 샘플 발견")
    print("  Normal count matrix 구축 중...")

    normal_matrix, gene_names_n = build_normal_count_matrix()
    if normal_matrix is None:
        print("  [WARN] Normal matrix 구축 실패 → partial_evidence")
        result = _classify_single_contrast(deg1, fdr_threshold, lfc_threshold)
        return result

    print(f"  Normal matrix: {normal_matrix.shape[0]} genes × {normal_matrix.shape[1]} samples")

    # TNBC count matrix 로드
    counts_path = os.path.join(PROCESSED_DIR, "tcga_brca_counts.csv")
    meta_path = os.path.join(PROCESSED_DIR, "sample_metadata.csv")

    counts = pd.read_csv(counts_path, index_col=0)
    counts = counts.drop(columns=["gene_name"], errors="ignore")
    meta = pd.read_csv(meta_path)
    tnbc_samples = meta[meta["subtype"] == "TNBC"]["sample_id"].tolist()

    # TNBC 샘플만 추출
    tnbc_counts = counts[[c for c in tnbc_samples if c in counts.columns]]
    print(f"  TNBC samples: {tnbc_counts.shape[1]}")

    # 중복 샘플 제거 (같은 환자의 복수 Normal)
    normal_matrix = normal_matrix.loc[:, ~normal_matrix.columns.duplicated()]
    tnbc_counts = tnbc_counts.loc[:, ~tnbc_counts.columns.duplicated()]
    # 컬럼 겹침 방지
    overlap = set(tnbc_counts.columns) & set(normal_matrix.columns)
    if overlap:
        normal_matrix = normal_matrix.drop(columns=list(overlap), errors="ignore")

    # 공통 유전자
    common_genes = tnbc_counts.index.intersection(normal_matrix.index)
    print(f"  Common genes: {len(common_genes)}")

    if len(common_genes) < 5000:
        print("  [WARN] 공통 유전자 부족 → partial_evidence")
        result = _classify_single_contrast(deg1, fdr_threshold, lfc_threshold)
        return result

    # DESeq2 실행 (TNBC vs Normal)
    print("  PyDESeq2 실행 (TNBC vs Normal)...")
    try:
        os.environ["LOKY_MAX_CPU_COUNT"] = "8"
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        # 유전자 중복 제거
        tnbc_sub = tnbc_counts.loc[common_genes]
        normal_sub = normal_matrix.loc[common_genes]
        tnbc_sub = tnbc_sub[~tnbc_sub.index.duplicated(keep="first")]
        normal_sub = normal_sub[~normal_sub.index.duplicated(keep="first")]
        common_genes2 = tnbc_sub.index.intersection(normal_sub.index)

        # Merge
        merged = pd.concat([tnbc_sub.loc[common_genes2], normal_sub.loc[common_genes2]], axis=1)
        merged = merged.fillna(0).astype(int)

        # Metadata
        groups = (["TNBC"] * tnbc_counts.shape[1] + ["Normal"] * normal_matrix.shape[1])
        clinical = pd.DataFrame({"group": groups}, index=merged.columns)

        # Filter low counts
        keep = (merged > 0).sum(axis=1) >= 5
        merged = merged[keep]
        print(f"  After filter: {merged.shape[0]} genes × {merged.shape[1]} samples")

        # DESeq2
        dds = DeseqDataSet(counts=merged.T, metadata=clinical, design_factors="group", n_cpus=8)
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["group", "TNBC", "Normal"])
        stat.summary()

        deg2 = stat.results_df.copy()
        deg2 = deg2.sort_values("padj")

        # gene_name 추가
        if "gene_name" in deg1.columns:
            gene_name_map = dict(zip(deg1.index, deg1["gene_name"]))
            deg2["gene_name"] = deg2.index.map(lambda x: gene_name_map.get(x, ""))

        sig2 = deg2[(deg2["padj"] < fdr_threshold) & (deg2["log2FoldChange"].abs() > lfc_threshold)]
        print(f"  Significant: {len(sig2)} (Up={sum(sig2['log2FoldChange']>0)}, Down={sum(sig2['log2FoldChange']<0)})")

        # 저장
        out_dir = os.path.join(PROCESSED_DIR, "deg_results")
        os.makedirs(out_dir, exist_ok=True)
        deg2.to_csv(os.path.join(out_dir, "deg_tnbc_vs_normal_full.csv"))
        sig2.to_csv(os.path.join(out_dir, "deg_tnbc_vs_normal_significant.csv"))
        print(f"  Saved: deg_tnbc_vs_normal_*.csv")

    except Exception as e:
        print(f"  [ERROR] DESeq2 실패: {e}")
        result = _classify_single_contrast(deg1, fdr_threshold, lfc_threshold)
        return result

    # 4. Dual contrast 통합 분류
    result = _classify_dual(deg1, deg2, fdr_threshold, lfc_threshold)
    return result


def _classify_single_contrast(deg1, fdr_thr, lfc_thr):
    """단일 contrast 분류 (Normal 없음)."""
    sig = deg1[(deg1["padj"] < fdr_thr) & (deg1["log2FoldChange"].abs() > lfc_thr)]
    classification = pd.DataFrame(index=sig.index)
    classification["contrast1_direction"] = sig["log2FoldChange"].apply(lambda x: "Up" if x > 0 else "Down")
    classification["contrast1_fdr"] = sig["padj"]
    classification["contrast2_direction"] = "UNAVAILABLE"
    classification["contrast2_fdr"] = np.nan
    classification["dual_status"] = "single_contrast_only"
    if "gene_name" in sig.columns:
        classification["gene_name"] = sig["gene_name"]

    out_path = os.path.join(PROCESSED_DIR, "deg_results", "dual_contrast_classification.csv")
    classification.to_csv(out_path)
    print(f"\nDual contrast: single_contrast_only ({len(classification)} genes)")
    print(f"Saved: {out_path}")

    return {"status": "single_contrast_only", "n_genes": len(classification),
            "contrast1_sig": len(sig), "contrast2_sig": 0}


def _classify_dual(deg1, deg2, fdr_thr, lfc_thr):
    """두 비교 결과 통합 분류."""
    # 공통 유전자
    common = deg1.index.intersection(deg2.index)
    print(f"\n[Dual Classification] {len(common)} common genes")

    results = []
    for gene in common:
        r1 = deg1.loc[gene] if gene in deg1.index else {}
        r2 = deg2.loc[gene] if gene in deg2.index else {}

        fdr1 = r1.get("padj", 1)
        fdr2 = r2.get("padj", 1)
        lfc1 = r1.get("log2FoldChange", 0)
        lfc2 = r2.get("log2FoldChange", 0)

        sig1 = (pd.notna(fdr1) and fdr1 < fdr_thr and abs(lfc1) > lfc_thr)
        sig2 = (pd.notna(fdr2) and fdr2 < fdr_thr and abs(lfc2) > lfc_thr)

        dir1 = "Up" if lfc1 > 0 else "Down" if lfc1 < 0 else "NS"
        dir2 = "Up" if lfc2 > 0 else "Down" if lfc2 < 0 else "NS"

        if sig1 and sig2 and dir1 == dir2:
            status = "dual_contrast_supported"
        elif sig1 and sig2 and dir1 != dir2:
            status = "discordant"
        elif sig1 or sig2:
            status = "partial_evidence"
        else:
            status = "not_significant"

        results.append({
            "gene_id": gene,
            "gene_name": r1.get("gene_name", ""),
            "c1_lfc": round(lfc1, 3) if pd.notna(lfc1) else None,
            "c1_fdr": fdr1,
            "c1_direction": dir1 if sig1 else "NS",
            "c2_lfc": round(lfc2, 3) if pd.notna(lfc2) else None,
            "c2_fdr": fdr2,
            "c2_direction": dir2 if sig2 else "NS",
            "dual_status": status,
        })

    df = pd.DataFrame(results).set_index("gene_id")

    # 요약
    status_counts = df["dual_status"].value_counts()
    print(f"\nDual Contrast Classification:")
    for s, n in status_counts.items():
        print(f"  {s}: {n}")

    # 저장
    out_path = os.path.join(PROCESSED_DIR, "deg_results", "dual_contrast_classification.csv")
    df.to_csv(out_path)
    print(f"\nSaved: {out_path}")

    return {"status": "dual_contrast_completed",
            "classification": dict(status_counts),
            "n_dual_supported": int(status_counts.get("dual_contrast_supported", 0)),
            "n_discordant": int(status_counts.get("discordant", 0)),
            "n_partial": int(status_counts.get("partial_evidence", 0))}


if __name__ == "__main__":
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")
    run_dual_contrast_deg()
