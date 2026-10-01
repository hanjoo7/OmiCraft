"""
CELLxGENE + genesorteR + cell-of-origin 판정 — 워크플로우 v5 Section 04–05A

단계:
  1. CELLxGENE Census API로 참조 slice 조회 (Python)
  2. genesorteR::sortGenes()를 R subprocess로 실행
  3. donor bootstrap 안정성 평가
  4. SINGLE_SUPPORTED / MULTI_SUPPORTED / UNRESOLVED / INSUFFICIENT 판정
  5. cell_context.tsv 출력

주의사항 (워크플로우 §04):
  - epithelial을 malignant로 자동 변환하지 않는다.
  - 음수 값 포함 scaled/integrated matrix 대신 raw/비음수 layer를 사용한다.
  - cell 수를 biological replicate로 계산하지 않는다.
  - specificity score만으로 cell-of-origin을 단정하지 않는다.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

# ── cell-of-origin 상태 ────────────────────────────────────

CONTEXT_STATUS = {
    "SINGLE_SUPPORTED": "한 broad class가 donor 반복에서 우세하고 검출률·annotation 지지",
    "MULTI_SUPPORTED":  "둘 이상의 class에서 반복 가능한 발현 확인",
    "UNRESOLVED":       "상위 점수 비슷하거나 donor별 순위 교차 또는 annotation 충돌",
    "INSUFFICIENT":     "유전자 미측정, donor/세포 부족, 정상 reference만 존재",
}

# broad cell type label 표준화 맵 (census 다양한 표기 → broad class)
BROAD_CLASS_MAP = {
    "malignant": "malignant",
    "tumor": "malignant",
    "cancer": "malignant",
    "epithelial": "epithelial",
    "luminal": "epithelial",
    "basal": "epithelial",
    "fibroblast": "fibroblast",
    "stromal": "fibroblast",
    "macrophage": "macrophage",
    "monocyte": "macrophage",
    "dendritic": "macrophage",
    "endothelial": "endothelial",
    "t cell": "lymphocyte",
    "b cell": "lymphocyte",
    "nk cell": "lymphocyte",
    "lymphocyte": "lymphocyte",
    "plasma": "lymphocyte",
    "mast": "mast_cell",
    "pericyte": "pericyte",
    "smooth muscle": "smooth_muscle",
}

# ── genesorteR R 스크립트 ─────────────────────────────────

_GENESORTER_R_SCRIPT = r"""
suppressPackageStartupMessages({
  library(genesorteR)
  library(Matrix)
  library(jsonlite)
})

args       <- commandArgs(trailingOnly = TRUE)
expr_rds   <- args[1]   # sparse matrix RDS (genes x cells, non-negative)
meta_rds   <- args[2]   # cell metadata RDS: data.frame with cell_class column
output_dir <- args[3]
genes_json <- if (length(args) >= 4) args[4] else NULL  # gene 목록 JSON

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

write_result <- function(status, reason, data=NULL) {
  out <- list(status=status, reason=reason)
  if (!is.null(data)) out <- c(out, data)
  write_json(out, file.path(output_dir, "genesorter_status.json"), auto_unbox=TRUE)
  quit(status=0)
}

# 입력 로드
expr <- tryCatch(readRDS(expr_rds), error=function(e) NULL)
meta <- tryCatch(readRDS(meta_rds), error=function(e) NULL)
if (is.null(expr)) write_result("FAILED", paste("발현 행렬 로드 실패:", expr_rds))
if (is.null(meta)) write_result("FAILED", paste("메타데이터 로드 실패:", meta_rds))

if (!"cell_class" %in% colnames(meta))
  write_result("FAILED", "메타데이터에 cell_class 컬럼 없음")

# 세포 순서 일치
common_cells <- intersect(colnames(expr), rownames(meta))
if (length(common_cells) < 20)
  write_result("INSUFFICIENT", sprintf("공통 세포 %d개 (최소 20 필요)", length(common_cells)))

expr <- expr[, common_cells]
meta <- meta[common_cells, , drop=FALSE]

# 특정 유전자만 사용
if (!is.null(genes_json) && file.exists(genes_json)) {
  genes <- fromJSON(genes_json)
  keep  <- intersect(genes, rownames(expr))
  if (!length(keep)) write_result("INSUFFICIENT", "No requested genes measured")
  expr <- expr[keep, , drop=FALSE]
}

# 음수 체크
if (any(expr < 0, na.rm=TRUE))
  write_result("FAILED", "음수 발현 값 감지 — raw/비음수 layer를 사용해야 합니다")

cell_class <- as.character(meta$cell_class)
n_classes  <- length(unique(cell_class))
if (n_classes < 2)
  write_result("INSUFFICIENT", sprintf("cell class %d개 (최소 2 필요)", n_classes))

# genesorteR 실행
sg <- tryCatch(
  sortGenes(expr, cell_class, binarizeMethod="naive", cores=1),
  error=function(e) write_result("FAILED", paste("sortGenes 오류:", conditionMessage(e)))
)

spec_df <- as.data.frame(as.matrix(sg$specScore))
spec_df$gene <- rownames(spec_df)

# TSV 저장
spec_path <- file.path(output_dir, "spec_score.tsv")
write.table(spec_df, spec_path, sep="\t", quote=FALSE, row.names=FALSE)

# Donor-level detection preserves the biological replicate as the donor.
detection_rows <- list()
set.seed(11)
if ("donor_id" %in% names(meta)) {
  for (g in rownames(expr)) {
    detected <- as.numeric(expr[g, ] > 0)
    by_donor <- aggregate(detected, list(donor=meta$donor_id, class=meta$cell_class), mean)
    values <- split(by_donor$x, by_donor$class)
    means <- vapply(values, mean, numeric(1))
    boot <- replicate(100, {
      rates <- vapply(values, function(x) mean(sample(x, length(x), replace=TRUE)), numeric(1))
      if (sum(rates == max(rates)) != 1L) NA_character_ else names(which.max(rates))
    })
    for (cl in names(values)) detection_rows[[length(detection_rows)+1L]] <- data.frame(
      gene=g, cell_class=cl, detection_fraction=means[[cl]],
      expressing_donors=sum(values[[cl]] >= 0.05), n_donors=length(values[[cl]]),
      donor_bootstrap_top_fraction=mean(!is.na(boot) & boot == cl))
  }
}
detection_path <- file.path(output_dir, "donor_detection.tsv")
if (length(detection_rows)) write.table(do.call(rbind,detection_rows), detection_path,
                                       sep="\t",quote=FALSE,row.names=FALSE)

# donor 분포 요약
donor_col <- intersect(c("donor_id","donor","sample","patient_id"), colnames(meta))[1]
donor_summary <- list()
if (!is.na(donor_col)) {
  donor_summary <- tapply(meta[[donor_col]], meta$cell_class,
                          function(x) length(unique(x)))
}

write_result("COMPLETED", "genesorteR 완료",
  list(
    spec_path      = spec_path,
    detection_path = detection_path,
    n_genes        = nrow(spec_df),
    n_classes      = n_classes,
    class_names    = unique(cell_class),
    n_cells        = ncol(expr),
    donor_per_class = as.list(donor_summary)
  )
)
"""


# ── 유틸 ──────────────────────────────────────────────────

def _find_rscript() -> Optional[str]:
    from .r_tool import find_rscript
    return find_rscript()


def _broad_class(label: str) -> str:
    """세포 레이블을 broad class로 표준화한다."""
    lbl = label.lower()
    for key, broad in BROAD_CLASS_MAP.items():
        if key in lbl:
            return broad
    return "other"


def _save_json(data: dict, path: Path) -> None:
    try:
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError:
        pass


# ── CELLxGENE Census 조회 ─────────────────────────────────

def query_cellxgene_census(
    genes: list[str],
    tissue: str = "breast",
    disease: str = "breast carcinoma",
    census_version: str = "stable",
    max_cells: int = 50000,
) -> dict:
    """CELLxGENE Census에서 관심 유전자 slice를 조회한다.

    Returns:
        {
          "status": "COMPLETED" | "UNAVAILABLE" | "FAILED",
          "anndata_path": str | None,  # 임시 h5ad 경로
          "n_cells": int,
          "n_donors": int,
          "datasets": list[str],
          "measured_genes": list[str],
          "unmeasured_genes": list[str],
          "reason": str,
        }
    """
    try:
        import cellxgene_census
    except ImportError:
        return {
            "status": "UNAVAILABLE",
            "reason": "cellxgene_census 패키지 없음 — pip install cellxgene-census",
            "anndata_path": None, "n_cells": None, "n_donors": None,
            "datasets": [], "measured_genes": [], "unmeasured_genes": genes,
        }

    try:
        census = cellxgene_census.open_soma(census_version=census_version)
    except Exception as exc:
        return {
            "status": "UNAVAILABLE",
            "reason": f"Census 연결 실패: {exc}",
            "anndata_path": None, "n_cells": None, "n_donors": None,
            "datasets": [], "measured_genes": [], "unmeasured_genes": genes,
        }

    try:
        import hashlib
        exp = census["census_data"]["homo_sapiens"]
        var_df = exp.ms["RNA"].var.read(column_names=["feature_id", "feature_name"]).concat().to_pandas()
        measured_genes = [g for g in genes if g in set(var_df["feature_name"])]
        unmeasured_genes = [g for g in genes if g not in measured_genes]
        if not measured_genes:
            raise ValueError("No requested genes in Census reference")
        query = f"tissue_general == {tissue!r} and disease == {disease!r} and is_primary_data == True"
        obs = exp.obs.read(value_filter=query, column_names=["soma_joinid"]).concat().to_pandas()
        obs["priority"] = obs.soma_joinid.map(lambda x: hashlib.sha256(f"11:{x}".encode()).hexdigest())
        ids = obs.sort_values("priority").head(max_cells).soma_joinid.to_numpy()
        if not len(ids):
            raise ValueError("No cells matching disease and tissue")
        adata = cellxgene_census.get_anndata(census, organism="Homo sapiens", obs_coords=ids,
            var_value_filter=f"feature_name in {measured_genes!r}",
            obs_column_names=["cell_type", "tissue", "disease", "donor_id", "dataset_id", "assay"])
        # 임시 h5ad 저장
        with tempfile.NamedTemporaryFile(
            suffix=".h5ad", delete=False, prefix="cellxgene_"
        ) as tmp:
            tmp_path = tmp.name
        adata.write_h5ad(tmp_path)

        n_donors = adata.obs["donor_id"].nunique()
        datasets = list(adata.obs["dataset_id"].unique())

        return {
            "status": "COMPLETED",
            "reason": f"{len(adata)} cells, {n_donors} donors",
            "anndata_path": tmp_path,
            "n_cells": len(adata),
            "n_donors": n_donors,
            "datasets": datasets,
            "measured_genes": measured_genes,
            "query_provenance": {"census_version": census_version, "requested_disease": disease,
                "disease_filter_applied": True, "tissue": tissue, "max_cells": max_cells},
            "unmeasured_genes": unmeasured_genes,
        }
    except Exception as exc:
        return {
            "status": "FAILED",
            "reason": str(exc),
            "anndata_path": None, "n_cells": None, "n_donors": None,
            "datasets": [], "measured_genes": [], "unmeasured_genes": genes,
        }
    finally:
        try:
            census.close()
        except Exception:
            pass


# ── genesorteR 실행 ───────────────────────────────────────

def run_genesorter(
    expr_rds: str | Path,
    meta_rds: str | Path,
    output_dir: str | Path,
    genes: list[str] | None = None,
    timeout: int = 600,
) -> dict:
    """genesorteR::sortGenes()를 R subprocess로 실행한다.

    Args:
        expr_rds: gene × cell sparse matrix RDS (비음수 raw counts).
        meta_rds: cell metadata RDS (cell_class, donor_id 등 포함).
        output_dir: 출력 디렉토리.
        genes: 대상 유전자 목록. None이면 전체.
        timeout: R 실행 제한시간 (초).

    Returns:
        genesorter_status.json 내용.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rscript = _find_rscript()
    if not rscript:
        result = {"status": "UNAVAILABLE",
                  "reason": "Rscript 없음 — R을 설치하거나 conda 환경 활성화 필요"}
        _save_json(result, output_dir / "genesorter_status.json")
        return result

    genes_json_path = None
    if genes:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False, encoding="utf-8"
        ) as tmp:
            json.dump(genes, tmp)
            genes_json_path = tmp.name

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".R", delete=False, encoding="utf-8"
    ) as tmp:
        tmp.write(_GENESORTER_R_SCRIPT)
        r_tmp = tmp.name

    args = [rscript, r_tmp, str(expr_rds), str(meta_rds), str(output_dir)]
    if genes_json_path:
        args.append(genes_json_path)

    try:
        proc = subprocess.run(
            args, capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        result = {"status": "FAILED",
                  "reason": f"genesorteR 실행 제한시간 초과 ({timeout}s)"}
        _save_json(result, output_dir / "genesorter_status.json")
        return result
    finally:
        for p in [r_tmp, genes_json_path]:
            if p and os.path.exists(p):
                os.unlink(p)

    if proc.returncode != 0:
        err = (proc.stderr or "")[-400:]
        result = {"status": "FAILED",
                  "reason": f"R returncode={proc.returncode}: {err}"}
        _save_json(result, output_dir / "genesorter_status.json")
        return result

    status_path = output_dir / "genesorter_status.json"
    if status_path.is_file():
        try:
            return json.loads(status_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            pass

    return {"status": "FAILED", "reason": "genesorter_status.json 없음"}


# ── cell-of-origin 판정 ───────────────────────────────────

def assign_context_status(
    spec_tsv: str | Path,
    donor_per_class: dict[str, int],
    gap_threshold: float = 0.10,
    min_donors: int = 3,
    min_detection: float = 0.05,
    gene: str | None = None,
    detection_tsv: str | Path | None = None,
) -> dict:
    """specScore TSV를 읽고 cell-of-origin 상태를 판정한다.

    판정 로직 (워크플로우 §05A):
      1. donor_per_class로 각 class의 donor 수 확인
      2. 상위 class의 specScore gap 계산
      3. gap > gap_threshold & donor 충분 → SINGLE_SUPPORTED
      4. 복수 class 안정 발현 → MULTI_SUPPORTED
      5. gap 좁거나 donor 부족 → UNRESOLVED
      6. 미측정 / 데이터 부족 → INSUFFICIENT

    Returns:
        {
          "context_status": str,
          "top_class": str | None,
          "top_score": float,
          "second_class": str | None,
          "second_score": float,
          "gap": float,
          "all_classes": dict[str, float],
          "donor_per_class": dict[str, int],
          "annotation_note": str,
        }
    """
    spec_path = Path(spec_tsv)
    if not spec_path.is_file():
        return {
            "context_status": "INSUFFICIENT",
            "reason": f"spec_score.tsv 없음: {spec_path}",
            "top_class": None, "top_score": 0.0,
            "second_class": None, "second_score": 0.0,
            "gap": 0.0, "all_classes": {}, "donor_per_class": donor_per_class,
        }

    # specScore 읽기 (행: gene, 열: cell class)
    try:
        import csv
        with spec_path.open(encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = [row for row in reader if gene is None or row.get("gene") == gene]
    except OSError as exc:
        return {"context_status": "INSUFFICIENT", "reason": str(exc),
                "top_class": None, "top_score": 0.0,
                "second_class": None, "second_score": 0.0,
                "gap": 0.0, "all_classes": {}, "donor_per_class": donor_per_class}

    if not rows:
        return {"context_status": "INSUFFICIENT", "reason": "spec_score.tsv 행 없음",
                "top_class": None, "top_score": 0.0,
                "second_class": None, "second_score": 0.0,
                "gap": 0.0, "all_classes": {}, "donor_per_class": donor_per_class}

    class_cols = [c for c in rows[0].keys() if c != "gene"]
    if not class_cols:
        return {"context_status": "INSUFFICIENT", "reason": "specScore 클래스 컬럼 없음",
                "top_class": None, "top_score": 0.0,
                "second_class": None, "second_score": 0.0,
                "gap": 0.0, "all_classes": {}, "donor_per_class": donor_per_class}

    # A per-gene call uses that gene's complete score vector, including the runner-up.
    class_wins = {c: 0.0 for c in class_cols}
    for row in rows:
        for c in class_cols:
            try:
                value = float(row[c])
                if math.isfinite(value):
                    class_wins[c] += value
            except (ValueError, TypeError, KeyError):
                pass

    if not class_wins:
        return {"context_status": "INSUFFICIENT", "reason": "specScore 데이터 없음",
                "top_class": None, "top_score": 0.0,
                "second_class": None, "second_score": 0.0,
                "gap": 0.0, "all_classes": {}, "donor_per_class": donor_per_class}

    sorted_classes = sorted(class_wins.items(), key=lambda x: x[1], reverse=True)
    top_class, top_score = sorted_classes[0]
    second_class = sorted_classes[1][0] if len(sorted_classes) > 1 else None
    second_score = sorted_classes[1][1] if len(sorted_classes) > 1 else 0.0

    eps = 1e-9
    gap = (top_score - second_score) / max(top_score, eps)
    all_classes = {c: round(s, 4) for c, s in sorted_classes}

    # donor 수 확인
    top_donors = donor_per_class.get(top_class, 0)

    # 판정
    note_parts = []

    # malignant/epithelial 주의 (워크플로우 §04)
    if _broad_class(top_class) == "epithelial":
        note_parts.append(
            "epithelial class — malignant로 자동 변환하지 않음; "
            "악성 주석·copy-number 근거 별도 확인 필요"
        )

    if top_donors < min_donors:
        context_status = "INSUFFICIENT"
        note_parts.append(f"상위 class ({top_class}) donor {top_donors}개 < 최소 {min_donors}")
    elif gap > gap_threshold and top_donors >= min_donors:
        context_status = "SINGLE_SUPPORTED"
    elif second_class and second_score > 0:
        # 둘 다 안정적 발현인지 확인 (donor 기준)
        second_donors = donor_per_class.get(second_class, 0)
        if second_donors >= min_donors:
            context_status = "MULTI_SUPPORTED"
            note_parts.append(
                f"{top_class}(donors={top_donors}) + {second_class}(donors={second_donors}) "
                "둘 다 충분한 donor에서 확인"
            )
        else:
            context_status = "UNRESOLVED"
            note_parts.append("gap 좁거나 2위 class donor 부족")
    else:
        context_status = "UNRESOLVED"
        note_parts.append(f"gap={gap:.3f} — 추가 dataset/protein/spatial 검증 필요")

    detection = {}
    if gene is not None:
        if detection_tsv and Path(detection_tsv).is_file():
            with Path(detection_tsv).open() as handle:
                detection = next((r for r in csv.DictReader(handle, delimiter="\t")
                                  if r["gene"] == gene and r["cell_class"] == top_class), {})
        if not detection:
            context_status = "UNRESOLVED" if top_score > 0 else "INSUFFICIENT"
            note_parts.append("Donor detection/bootstrap evidence unavailable")
        elif float(detection["detection_fraction"]) < min_detection or int(detection["expressing_donors"]) < min_donors:
            context_status = "INSUFFICIENT"
            note_parts.append("Insufficient detection across biological donors")
        elif float(detection["donor_bootstrap_top_fraction"]) < 0.8:
            context_status = "UNRESOLVED"
            note_parts.append("Top class is not stable in donor detection bootstrap")

    return {
        "detection": detection,
        "context_status": context_status,
        "top_class": top_class,
        "broad_top_class": _broad_class(top_class),
        "top_score": round(top_score, 4),
        "second_class": second_class,
        "second_score": round(second_score, 4),
        "gap": round(gap, 4),
        "gap_threshold_used": gap_threshold,
        "all_classes": all_classes,
        "donor_per_class": donor_per_class,
        "top_class_donors": top_donors,
        "annotation_note": "; ".join(note_parts) if note_parts else "",
        "status_description": CONTEXT_STATUS.get(context_status, ""),
    }


# ── cell_context.tsv 작성 ─────────────────────────────────

def write_cell_context_tsv(
    genes: list[str],
    spec_tsv: str | Path,
    donor_per_class: dict[str, int],
    output_path: str | Path,
    dataset_info: dict | None = None,
) -> str:
    """cell_context.tsv를 작성한다.

    컬럼: gene, cell_class, specificity, detection_fraction,
           measured_status, n_donors, stability, annotation_provenance

    워크플로우 §04 출력 형식.
    """
    spec_path = Path(spec_tsv)
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "gene", "cell_class", "specificity", "detection_fraction",
        "measured_status", "n_donors", "annotation_provenance",
        "context_status", "note",
    ]

    rows = []
    if not spec_path.is_file():
        # 미측정 기록
        for gene in genes:
            rows.append({
                "gene": gene,
                "cell_class": "UNKNOWN",
                "specificity": "",
                "detection_fraction": "",
                "measured_status": "NOT_MEASURED",
                "n_donors": "",
                "annotation_provenance": (dataset_info.get("datasets") or ["unknown"])[0]
                    if dataset_info else "unknown",
                "context_status": "INSUFFICIENT",
                "note": "spec_score.tsv 없음",
            })
    else:
        try:
            import csv
            with spec_path.open(encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                spec_rows = {r["gene"]: r for r in reader if "gene" in r}
        except OSError:
            spec_rows = {}

        class_cols = []
        if spec_rows:
            first = next(iter(spec_rows.values()))
            class_cols = [c for c in first.keys() if c != "gene"]

        for gene in genes:
            if gene not in spec_rows:
                rows.append({
                    "gene": gene, "cell_class": "UNKNOWN",
                    "specificity": "", "detection_fraction": "",
                    "measured_status": "NOT_IN_SPEC_SCORE",
                    "n_donors": "", "annotation_provenance": "unknown",
                    "context_status": "INSUFFICIENT", "note": "specScore 없음",
                })
                continue

            row = spec_rows[gene]
            scores = {}
            for c in class_cols:
                try:
                    value = float(row[c])
                    if math.isfinite(value):
                        scores[c] = value
                except (ValueError, KeyError):
                    pass

            for cell_class, spec in sorted(scores.items(),
                                           key=lambda x: x[1], reverse=True):
                n_donors = donor_per_class.get(cell_class, 0)
                rows.append({
                    "gene": gene,
                    "cell_class": cell_class,
                    "specificity": round(spec, 4),
                    "detection_fraction": "",  # genesorteR는 직접 제공 안 함
                    "measured_status": "MEASURED",
                    "n_donors": n_donors,
                    "annotation_provenance": (dataset_info.get("datasets") or ["unknown"])[0]
                        if dataset_info else "unknown",
                    "context_status": "",
                    "note": "",
                })

    # TSV 쓰기
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(r.get(k, "")) for k in header))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(out_path)


# ── 통합 실행 함수 ────────────────────────────────────────

from .resource_usage import measured_tool


@measured_tool
def run_cell_context_analysis(
    genes: list[str],
    output_dir: str | Path,
    expr_rds: str | Path | None = None,
    meta_rds: str | Path | None = None,
    census_tissue: str = "breast",
    census_disease: str = "breast carcinoma",
    gap_threshold: float = 0.10,
    min_donors: int = 3,
) -> dict:
    """CELLxGENE + genesorteR 통합 실행.

    expr_rds / meta_rds가 제공되면 직접 사용.
    없으면 CELLxGENE Census에서 slice를 조회한 뒤 genesorteR를 실행한다.

    Returns:
        {
          "status": str,
          "context_results": {gene: context_status_dict},
          "cell_context_tsv": str | None,
          "genesorter_result": dict,
          "census_result": dict,
        }
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    from .configuration import get_config
    provenance = get_config().upstream.census_provenance_path
    census_result: dict = json.loads(Path(provenance).read_text()) if provenance and Path(provenance).is_file() else {"status": "SKIPPED", "datasets": ["local_reference"]}
    genesorter_result: dict = {"status": "SKIPPED"}

    # ── 1. 발현 데이터 확보 ───────────────────────────────
    if expr_rds is None or not Path(expr_rds).is_file():
        census_result = query_cellxgene_census(
            genes=genes,
            tissue=census_tissue,
            disease=census_disease,
        )
        if census_result["status"] not in ("COMPLETED",):
            result = {
                "status": census_result["status"],
                "reason": census_result["reason"],
                "context_results": {g: {"context_status": "INSUFFICIENT",
                                        "reason": census_result["reason"]}
                                    for g in genes},
                "cell_context_tsv": None,
                "genesorter_result": genesorter_result,
                "census_result": census_result,
            }
            _save_json(result, output_dir / "cell_context_status.json")
            return result

        expr_rds, meta_rds = convert_h5ad_to_rds(census_result["anndata_path"], output_dir / "reference")

    reference_measured = census_result.pop("measured_genes", None)
    analysis_genes = [gene for gene in genes if reference_measured is None or gene in set(reference_measured)]
    census_result["measured_genes"] = analysis_genes
    census_result["unmeasured_genes"] = [gene for gene in genes if gene not in analysis_genes]

    # ── 2. genesorteR 실행 ────────────────────────────────
    gs_dir = output_dir / "genesorter"
    if not analysis_genes:
        result = {"status": "INSUFFICIENT", "reason": "Requested genes unmeasured in reference datasets",
                  "context_results": {g: {"context_status": "INSUFFICIENT", "measured_status": "NOT_MEASURED"} for g in genes}}
        _save_json(result, output_dir / "cell_context_status.json")
        return result
    genesorter_result = run_genesorter(
        expr_rds=expr_rds,
        meta_rds=meta_rds,
        output_dir=gs_dir,
        genes=analysis_genes,
    )

    if genesorter_result.get("status") != "COMPLETED":
        result = {
            "status": genesorter_result.get("status", "FAILED"),
            "reason": genesorter_result.get("reason", "genesorteR 실패"),
            "context_results": {g: {"context_status": "INSUFFICIENT",
                                    "reason": genesorter_result.get("reason", "")}
                                 for g in genes},
            "cell_context_tsv": None,
            "genesorter_result": genesorter_result,
            "census_result": census_result,
        }
        _save_json(result, output_dir / "cell_context_status.json")
        return result

    # ── 3. context 판정 ───────────────────────────────────
    spec_tsv = Path(genesorter_result.get("spec_path", ""))
    donor_per_class = genesorter_result.get("donor_per_class", {})

    context_results: dict[str, dict] = {}
    for gene in genes:
        # 유전자 단위로 특이도가 있는 클래스 집계
        ctx = assign_context_status(
            spec_tsv=spec_tsv,
            donor_per_class=donor_per_class,
            gap_threshold=gap_threshold,
            min_donors=min_donors,
            gene=gene, detection_tsv=genesorter_result.get("detection_path"),
        )
        ctx["gene"] = gene
        context_results[gene] = ctx

    # ── 4. cell_context.tsv 저장 ──────────────────────────
    tsv_path = write_cell_context_tsv(
        genes=genes,
        spec_tsv=spec_tsv,
        donor_per_class=donor_per_class,
        output_path=output_dir / "cell_context.tsv",
        dataset_info=census_result,
    )

    result = {
        "status": "COMPLETED",
        "reason": f"{len(genes)} 유전자 cell context 분석 완료",
        "context_results": context_results,
        "cell_context_tsv": tsv_path,
        "genesorter_result": genesorter_result,
        "census_result": census_result,
    }
    _save_json(result, output_dir / "cell_context_status.json")
    return result


def convert_h5ad_to_rds(path, output_dir):
    import anndata
    from scipy.io import mmwrite
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = anndata.read_h5ad(path)
    names = data.var["feature_name"]
    keep = ~names.duplicated() & names.notna()
    mmwrite(output / "expression.mtx", data.X[:, keep.to_numpy()].T)
    data.var.loc[keep].to_csv(output / "genes.tsv", sep="\t", index=False)
    cells = data.obs.copy()
    cells["cell_id"] = cells.index.astype(str)
    cells["cell_class"] = cells["cell_type"].astype(str)
    cells["donor_id"] = cells["dataset_id"].astype(str) + ":" + cells["donor_id"].astype(str)
    cells.to_csv(output / "cells.tsv", sep="\t", index=False)
    script = Path(__file__).parent / "scripts/census_to_rds.R"
    subprocess.run([_find_rscript(), "--vanilla", str(script), str(output)], check=True,
                   capture_output=True, text=True, timeout=600)
    return str(output / "expression.rds"), str(output / "metadata.rds")
