"""
aPEAR GO network 어댑터 — 워크플로우 v5 Section 03

R aPEAR::enrichmentNetwork()를 subprocess로 호출한다.
R 또는 aPEAR 패키지가 없으면 Python leading-edge Jaccard fallback으로 대체한다.

필수 입력 컬럼: ID, Description, NES, p.adjust, core_enrichment
  (clusterProfiler / fgsea wrapper 표준 출력 형식)

출력:
  apear_status.json  — 실행 상태 및 파일 경로
  GOBP_network_{contrast}.png — 네트워크 플롯 (R 경로)
  apear_nodes.tsv / apear_edges.tsv — 노드·엣지 데이터
  pathway_network.json — Python fallback 결과
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional

# ── R 스크립트 (임시 파일로 저장 후 실행) ─────────────────

_APEAR_R_SCRIPT = r"""
suppressPackageStartupMessages({
  library(aPEAR)
  library(jsonlite)
  library(ggplot2)
})
args            <- commandArgs(trailingOnly = TRUE)
input_tsv       <- args[1]
output_dir      <- args[2]
contrast        <- if (length(args) >= 3) args[3] else "contrast"
min_cluster     <- if (length(args) >= 4) as.integer(args[4]) else 2
p_cutoff        <- if (length(args) >= 5) as.numeric(args[5])  else 0.05

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_empty <- function(reason, status="EMPTY") {
  write_json(list(status=status, reason=reason, n_sig_paths=0L,
                  n_nodes=0L, n_edges=0L),
             file.path(output_dir, "apear_status.json"), auto_unbox=TRUE)
  quit(status=0)
}

df <- tryCatch(read.delim(input_tsv, stringsAsFactors=FALSE, check.names=FALSE),
               error=function(e) NULL)
if (is.null(df) || nrow(df) == 0) write_empty("입력 GSEA 결과 없음")

required <- c("ID","Description","NES","p.adjust","core_enrichment")
missing  <- setdiff(required, colnames(df))
if (length(missing) > 0) write_empty(paste("필수 컬럼 없음:", paste(missing, collapse=", ")), "INVALID_INPUT")

sig <- df[!is.na(df$p.adjust) & df$p.adjust < p_cutoff, ]
if (nrow(sig) < 2) write_empty(sprintf("유의한 GO BP %d건 (최소 2, p.adjust<%.2f)", nrow(sig), p_cutoff))

n_total_sig <- nrow(sig)
sig <- head(sig[order(sig$p.adjust, -abs(sig$NES), sig$ID), ], 200)
write.table(sig, file.path(output_dir, "network_input.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
set.seed(11)
net <- tryCatch(
  aPEAR::enrichmentNetwork(sig, colorBy="NES", drawEllipses=FALSE,
                           minClusterSize=min_cluster, fontSize=3.2,
                           repelLabels=TRUE, plotOnly=FALSE),
  error=function(e) write_empty(paste("aPEAR 오류:", conditionMessage(e)), "FAILED")
)

for (index in seq_along(net$plot$layers)) {
  layer <- net$plot$layers[[index]]
  if (is.data.frame(layer$data) && "label" %in% names(layer$data)) {
    layer$data$label <- vapply(layer$data$label, function(label) {
      paste(strwrap(gsub("_", " ", sub("^GOBP_", "", label)), width=32), collapse="\n")
    }, character(1))
    layer$geom_params$max.overlaps <- Inf
    layer$geom_params$seed <- 11
    layer$geom_params$box.padding <- 0.6
    net$plot$layers[[index]] <- layer
  }
}
net$plot <- net$plot +
  scale_x_continuous(expand=expansion(mult=0.12)) +
  scale_y_continuous(expand=expansion(mult=0.12)) +
  labs(title="GO biological process network: TNBC vs Non-TNBC",
       caption="Positive NES: TNBC | Top 200 significant terms; full pathway names retained in tables") +
  theme(plot.margin=margin(20, 20, 20, 20), plot.title=element_text(size=14),
        plot.caption=element_text(size=9))
plot_path <- file.path(output_dir, paste0("GOBP_network_", contrast, ".png"))
tryCatch(ggsave(plot_path, plot=net$plot, width=14, height=11, dpi=150),
         error=function(e) write_empty(paste("Plot save failed:", conditionMessage(e)), "FAILED"))

nodes_path <- file.path(output_dir, "apear_nodes.tsv")
edges_path <- file.path(output_dir, "apear_edges.tsv")
nodes <- merge(sig, as.data.frame(net$clusters), by.x="Description", by.y="ID")
if (!nrow(nodes)) write_empty("No pathways retained by clustering")
sim <- aPEAR::pathwaySimilarity(sig, geneCol="core_enrichment", method="jaccard")
sim <- sim[nodes$Description, nodes$Description, drop=FALSE]
cluster_map <- setNames(nodes$Cluster, nodes$Description)
same <- outer(cluster_map[rownames(sim)], cluster_map[colnames(sim)], "==")
idx <- which(upper.tri(sim) & sim >= ifelse(same, 0.1, 0.5), arr.ind=TRUE)
edges <- data.frame(source=rownames(sim)[idx[,1]], target=colnames(sim)[idx[,2]], similarity=sim[idx])
write.table(nodes, nodes_path, sep="\t",quote=FALSE,row.names=FALSE)
write.table(edges, edges_path, sep="\t",quote=FALSE,row.names=FALSE)
n_nodes <- nrow(nodes)
n_edges <- nrow(edges)
write_json(
  list(status="COMPLETED", contrast=contrast, n_sig_paths=nrow(sig),
       n_nodes=n_nodes, n_edges=n_edges, p_cutoff=p_cutoff,
       total_significant=n_total_sig, network_limit=200L,
       selection="Top 200 significant terms by adjusted p, then absolute NES; full GSEA table retained",
       plot_path=plot_path,
       nodes_path=if (file.exists(nodes_path)) nodes_path else NULL,
       edges_path=if (file.exists(edges_path)) edges_path else NULL),
  file.path(output_dir, "apear_status.json"), auto_unbox=TRUE
)
"""


# ── Python fallback ───────────────────────────────────────

def _jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a and b) else 0.0


def _python_network(gobp_tsv: Path, output_dir: Path,
                    p_cutoff: float = 0.05,
                    sim_threshold: float = 0.2) -> dict:
    """clusterProfiler 형식 TSV에서 Python Jaccard 네트워크 구축."""
    try:
        import csv
        rows = []
        with gobp_tsv.open(encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rows.append(row)
    except OSError as exc:
        return {"status": "FAILED", "reason": str(exc), "n_nodes": 0, "n_edges": 0}

    sig = []
    for r in rows:
        try:
            padj = float(r.get("p.adjust", 1.0) or 1.0)
        except ValueError:
            continue
        if padj < p_cutoff:
            sig.append(r)

    if len(sig) < 2:
        result = {"status": "EMPTY",
                  "reason": f"유의한 GO BP {len(sig)}건 (Python fallback)",
                  "n_sig_paths": len(sig), "n_nodes": 0, "n_edges": 0,
                  "method": "python_jaccard"}
        _save_json(result, output_dir / "apear_status.json")
        return result

    # leading-edge 유전자 파싱 (슬래시 또는 슬래시/공백 구분)
    pw_genes: dict[str, set] = {}
    for r in sig:
        pid = r.get("ID") or r.get("Description", "")
        le  = r.get("core_enrichment") or r.get("leading_edge") or ""
        pw_genes[pid] = set(le.replace("/", ";").split(";")) if le else set()

    terms = list(pw_genes.keys())
    edges = []
    for i in range(len(terms)):
        for j in range(i + 1, len(terms)):
            sim = _jaccard(pw_genes[terms[i]], pw_genes[terms[j]])
            if sim >= sim_threshold:
                edges.append({"source": terms[i], "target": terms[j],
                               "similarity": round(sim, 3)})

    # connected components clustering
    adj: dict[str, set] = defaultdict(set)
    for e in edges:
        adj[e["source"]].add(e["target"])
        adj[e["target"]].add(e["source"])
    visited: set = set()
    clusters = []
    for node in adj:
        if node in visited:
            continue
        comp, stack = [], [node]
        while stack:
            n = stack.pop()
            if n in visited:
                continue
            visited.add(n)
            comp.append(n)
            stack.extend(adj[n] - visited)
        clusters.append(sorted(comp))

    nodes = []
    for r in sig:
        pid = r.get("ID") or r.get("Description", "")
        try:
            nes = round(float(r.get("NES", 0) or 0), 3)
            padj = round(float(r.get("p.adjust", 1) or 1), 5)
        except ValueError:
            nes, padj = 0.0, 1.0
        nodes.append({"ID": pid,
                      "Description": r.get("Description", pid),
                      "NES": nes, "p.adjust": padj,
                      "n_leading_edge": len(pw_genes[pid])})

    result = {
        "status": "COMPLETED",
        "method": "python_jaccard",
        "n_sig_paths": len(sig),
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "n_clusters": len(clusters),
        "node_list": nodes,
        "edge_list": edges[:100],
        "clusters": clusters[:20],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _save_json(result, output_dir / "pathway_network.json")
    _save_json(result, output_dir / "apear_status.json")

    # TSV 저장
    _write_tsv([{"ID": n["ID"], "Description": n["Description"],
                 "NES": n["NES"], "p.adjust": n["p.adjust"]} for n in nodes],
               output_dir / "apear_nodes.tsv")
    _write_tsv(edges, output_dir / "apear_edges.tsv")

    return result


# ── R 경로 ────────────────────────────────────────────────

def _find_rscript() -> Optional[str]:
    from .r_tool import find_rscript
    return find_rscript()


def _validate_tsv(path: Path) -> tuple[bool, str, int]:
    required = {"ID", "Description", "NES", "p.adjust", "core_enrichment"}
    try:
        with path.open(encoding="utf-8") as f:
            header = set(f.readline().rstrip("\n").split("\t"))
        missing = required - header
        if missing:
            return False, f"필수 컬럼 없음: {missing}", 0
        n = sum(1 for _ in path.open(encoding="utf-8")) - 1
        return True, "OK", max(0, n)
    except OSError as exc:
        return False, str(exc), 0


# ── 공개 인터페이스 ────────────────────────────────────────

from .resource_usage import measured_tool


@measured_tool
def run_apear(
    gobp_tsv_path: str | Path,
    output_dir: str | Path,
    contrast: str = "contrast",
    p_cutoff: float = 0.05,
    min_cluster_size: int = 2,
    timeout: int = 300,
) -> dict:
    """aPEAR GO network 분석.

    R + aPEAR 패키지가 있으면 R 경로로, 없으면 Python Jaccard fallback 사용.

    Returns:
        dict with keys: status, reason, n_sig_paths, n_nodes, n_edges,
                        plot_path, nodes_path, edges_path, method
    """
    gobp_tsv_path = Path(gobp_tsv_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not gobp_tsv_path.is_file():
        result = _empty("UNAVAILABLE", f"GOBP GSEA TSV 없음: {gobp_tsv_path}")
        _save_json(result, output_dir / "apear_status.json")
        return result

    ok, reason, n_rows = _validate_tsv(gobp_tsv_path)
    if not ok:
        result = _empty("UNAVAILABLE", reason)
        _save_json(result, output_dir / "apear_status.json")
        return result

    if n_rows == 0:
        result = _empty("EMPTY", "GOBP GSEA 결과 0행")
        _save_json(result, output_dir / "apear_status.json")
        return result

    # R 경로 시도
    rscript = _find_rscript()
    if rscript:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".R", delete=False, encoding="utf-8"
        ) as tmp:
            tmp.write(_APEAR_R_SCRIPT)
            r_tmp = tmp.name
        try:
            proc = subprocess.run(
                [rscript, r_tmp,
                 str(gobp_tsv_path), str(output_dir),
                 contrast, str(min_cluster_size), str(p_cutoff)],
                capture_output=True, text=True, timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            os.unlink(r_tmp)
            # timeout → Python fallback
            result = _python_network(gobp_tsv_path, output_dir, p_cutoff)
            result["warning"] = f"R 제한시간 초과({timeout}s) — Python fallback"
            return result
        finally:
            if os.path.exists(r_tmp):
                os.unlink(r_tmp)

        if proc.returncode == 0:
            status_path = output_dir / "apear_status.json"
            if status_path.is_file():
                try:
                    return json.loads(status_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    pass
        # R 실패 → Python fallback
        err = (proc.stderr or "")[-300:]
        result = _python_network(gobp_tsv_path, output_dir, p_cutoff)
        result["warning"] = f"R 실행 실패 — Python fallback. R stderr: {err}"
        return result

    # Rscript 없음 → Python fallback
    result = _python_network(gobp_tsv_path, output_dir, p_cutoff)
    result["warning"] = "Rscript 없음 — Python Jaccard fallback 사용"
    return result


def run_apear_for_contrasts(
    gsea_output_dir: str | Path,
    output_dir: str | Path,
    contrast_ids: list[str] | None = None,
    p_cutoff: float = 0.05,
) -> dict:
    """복수 contrast의 GOBP GSEA 출력에 대해 aPEAR를 실행.

    Returns: {contrast_id: apear_result}
    """
    gsea_dir = Path(gsea_output_dir)
    out_dir  = Path(output_dir)

    candidates: dict[str, Path] = {}
    if contrast_ids:
        for cid in contrast_ids:
            for fname in [f"GOBP_GSEA_{cid}.tsv", "GOBP_GSEA_all.tsv"]:
                p = gsea_dir / fname
                if p.is_file():
                    candidates[cid] = p
                    break
    else:
        for p in gsea_dir.glob("GOBP_GSEA*.tsv"):
            cid = p.stem.replace("GOBP_GSEA_", "").replace("GOBP_GSEA", "all")
            candidates[cid] = p

    if not candidates:
        default = gsea_dir / "GOBP_GSEA_all.tsv"
        if default.is_file():
            candidates["all"] = default
        else:
            return {"_error": "GOBP GSEA TSV 없음"}

    return {
        cid: run_apear(tsv, out_dir / cid, contrast=cid, p_cutoff=p_cutoff)
        for cid, tsv in candidates.items()
    }


# ── 구버전 호환 ───────────────────────────────────────────

def build_pathway_network(gsea_path: str = None, sim_threshold: float = 0.2) -> dict:
    """이전 버전 호환용 wrapper. run_apear() 사용을 권장."""
    if gsea_path is None:
        return {"status": "UNAVAILABLE", "reason": "gsea_path 미지정", "nodes": 0, "edges": 0}
    result = _python_network(Path(gsea_path),
                             Path(gsea_path).parent / "apear_output",
                             sim_threshold=sim_threshold)
    return result


# ── 유틸 ──────────────────────────────────────────────────

def _empty(status: str, reason: str) -> dict:
    return {"status": status, "reason": reason,
            "n_sig_paths": 0, "n_nodes": 0, "n_edges": 0,
            "plot_path": None, "nodes_path": None, "edges_path": None}


def _save_json(data: dict, path: Path) -> None:
    try:
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError:
        pass


def _write_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    try:
        keys = list(rows[0].keys())
        lines = ["\t".join(keys)]
        for r in rows:
            lines.append("\t".join(str(r.get(k, "")) for k in keys))
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except OSError:
        pass
