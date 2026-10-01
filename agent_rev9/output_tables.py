"""
Output Tables — 워크플로우 v5 Section 4

워크플로우 §4에서 정의한 모든 출력 TSV / JSONL 파일을 생성한다.

파일 목록:
  gene_overview.tsv       — gene 행 단위 개요
  de_evidence.tsv         — gene × contrast × cohort 차등발현
  cell_context.tsv        — gene × cell class × dataset (cellxgene_adapter에서도 생성)
  dependency.tsv          — gene × model × technology (depmap_module에서도 생성)
  hypotheses.tsv          — gene × context × intervention (target_qualification에서도 생성)
  modality_options.tsv    — hypothesis × modality (M, D, Tier, 다음실험 등)
  evidence_cards.jsonl    — 독립 claim/evidence JSONL

빈 행이 없도록: 분석이 수행되지 않은 필드는 "" 대신 "UNKNOWN" / "NOT_EXECUTED" 사용.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

# ── 공통 헬퍼 ─────────────────────────────────────────────

def _tsv(rows: list[dict], keys: list[str], path: Path) -> str:
    """rows를 TSV로 저장. 경로 문자열 반환."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(keys)]
    for r in rows:
        lines.append("\t".join(
            str(r.get(k, "")).replace("\t", " ").replace("\n", " ")
            for k in keys
        ))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def _val(v: Any, default: str = "") -> str:
    if v is None:
        return default
    if isinstance(v, list):
        return "; ".join(str(x) for x in v)
    if isinstance(v, bool):
        return str(v).upper()
    return str(v)


# ── gene_overview.tsv ─────────────────────────────────────

def write_gene_overview(
    genes: list[str],
    qualification_results: dict[str, dict],
    de_evidences: dict[str, dict],
    cell_contexts: dict[str, dict],
    depmap_results: dict[str, dict],
    output_path: str | Path,
    user_selection: dict | None = None,
) -> str:
    """gene_overview.tsv 작성 (gene 단위 1행).

    워크플로우 §4: gene_id/symbol, 두 비교 요약, pathway 요약,
    origins, 대표 B, 후보 상태, route IDs, 최종 선택.
    """
    keys = [
        "gene",
        "dual_contrast_status",     # dual_contrast_supported | partial_evidence | discordant
        "log2FC_tnbc_nontnbc",
        "padj_tnbc_nontnbc",
        "log2FC_tnbc_normal",
        "padj_tnbc_normal",
        "top_pathway",
        "cell_origin_status",       # SINGLE_SUPPORTED | MULTI_SUPPORTED | UNRESOLVED | INSUFFICIENT
        "top_cell_class",
        "depmap_status",            # CONCORDANT | CRISPR_ONLY | DISCORDANT | UNAVAILABLE
        "crispr_median_effect",
        "bio_confidence",           # B0 | B1 | B2 | B3
        "judgment",                 # ADVANCE | HOLD | REJECT
        "best_tier",                # 1A | 1B | 2A | 2B | 3 | HOLD | REJECT
        "intervention_role",
        "safety_veto",
        "route_ids",
        "selected_route",           # 사용자 선택 후 채움
        "missing_evidence",
    ]
    rows = []
    for gene in genes:
        qr  = qualification_results.get(gene, {})
        de  = de_evidences.get(gene, {})
        ctx = cell_contexts.get(gene, {})
        dep = depmap_results.get(gene, {})
        crispr = dep.get("crispr_summary", {}) or {}
        sel = (user_selection or {}).get(gene, {})
        tier_results = qr.get("tier_results", [])
        route_ids = [f"{t['modality']}:{t['Tier']}" for t in tier_results]
        has_veto = any(
            sv.get("has_veto")
            for sv in (qr.get("safety_veto") or {}).values()
        )
        rows.append({
            "gene": gene,
            "dual_contrast_status": _val(de.get("dual_contrast"), "NOT_EVALUATED"),
            "log2FC_tnbc_nontnbc": _val(de.get("log2FC_tnbc_nontnbc")),
            "padj_tnbc_nontnbc": _val(de.get("padj_tnbc_nontnbc")),
            "log2FC_tnbc_normal": _val(de.get("log2FC_tnbc_normal")),
            "padj_tnbc_normal": _val(de.get("padj_tnbc_normal")),
            "top_pathway": _val(de.get("top_pathway")),
            "cell_origin_status": _val(ctx.get("context_status"), "INSUFFICIENT"),
            "top_cell_class": _val(ctx.get("top_class")),
            "depmap_status": _val(dep.get("dependency_status"), "UNAVAILABLE"),
            "crispr_median_effect": _val(crispr.get("median_effect")),
            "bio_confidence": _val(qr.get("bio_confidence"), "B3"),
            "judgment": _val(qr.get("shortlist", {}).get("judgment"), "HOLD"),
            "best_tier": _val(qr.get("best_tier"), "HOLD"),
            "intervention_role": _val(qr.get("intervention_role"), "UNKNOWN"),
            "safety_veto": "YES" if has_veto else "NO",
            "route_ids": _val(route_ids),
            "selected_route": _val(sel.get("selected_route")),
            "missing_evidence": _val(qr.get("missing_evidence")),
        })
    return _tsv(rows, keys, Path(output_path))


# ── de_evidence.tsv ───────────────────────────────────────

def write_de_evidence(
    de_records: list[dict],
    output_path: str | Path,
) -> str:
    """de_evidence.tsv 작성 (gene × contrast × cohort 행).

    de_records 항목 형식:
        gene, contrast_id, cohort, log2FC, SE, padj, wald_stat,
        prevalence, tested_status, external_replication, source_id
    """
    keys = [
        "gene", "contrast_id", "cohort",
        "log2FC", "SE", "padj", "wald_stat",
        "prevalence_positive", "tested_status",
        "external_replication", "source_id",
    ]
    rows = []
    for r in de_records:
        rows.append({k: _val(r.get(k), "UNKNOWN") for k in keys})
    return _tsv(rows, keys, Path(output_path))


def de_records_from_r_output(
    r_output: dict,
    contrast_id: str = "tnbc_vs_nontnbc",
    cohort: str = "TCGA-BRCA",
) -> list[dict]:
    """run_r_analysis() 결과에서 de_records를 생성한다."""
    records: list[dict] = []
    output_files = r_output.get("output_files", {})

    # DEG table에서 읽기
    deg_path = output_files.get("de_all") or output_files.get("deg_table")
    if not deg_path:
        return records

    import csv
    try:
        with open(deg_path, encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                records.append({
                    "gene": row.get("gene_symbol") or row.get("gene_name") or row.get("gene") or row.get(""),
                    "contrast_id": contrast_id,
                    "cohort": cohort,
                    "log2FC": row.get("log2FoldChange") or row.get("log2FC"),
                    "SE": row.get("lfcSE"),
                    "padj": row.get("padj"),
                    "wald_stat": row.get("stat"),
                    "prevalence_positive": row.get("prevalence"),
                    "tested_status": "TESTED",
                    "external_replication": "UNKNOWN",
                    "source_id": f"DESeq2_{contrast_id}_{cohort}",
                })
    except (OSError, KeyError):
        pass

    return records


# ── modality_options.tsv ──────────────────────────────────

def write_modality_options(
    genes: list[str],
    qualification_results: dict[str, dict],
    output_path: str | Path,
    user_selection: dict | None = None,
) -> str:
    """modality_options.tsv 작성 (hypothesis × modality 행).

    워크플로우 §08: M, D, Tier, entry criteria, missing evidence,
    asset, cost, next experiment, status.
    """
    keys = [
        "gene", "modality",
        "bio_confidence", "M", "D", "Tier",
        "entry_criteria", "missing_evidence",
        "asset", "estimated_cost",
        "next_experiment", "execution_status",
        "selected_by", "selected_at",
    ]

    # 다음 실험 제안 맵
    next_experiment_map = {
        "SMALL_MOLECULE": "pocket 검증 → AutoDock Vina/GNINA docking → ADMET",
        "DE_NOVO_BINDER": "RFdiffusion3 골격 생성 → ProteinMPNN 서열 설계 → AF3 복합체 검증",
        "ADC": "내재화 assay → linker/payload 선택 → 접합체 평가",
        "CAR_T": "scFv 설계 → 살상 assay → on-target/off-tumor 평가",
        "SIRNA": "sequence gate → Bowtie off-target → 화학변형·전달 평가",
        "PROTAC": "E3 ligase 공존 확인 → linker 최적화 → 삼자 복합체 모델",
    }
    # 비용 수준 맵
    cost_map = {
        "SMALL_MOLECULE": "LOW-MEDIUM",
        "DE_NOVO_BINDER": "MEDIUM-HIGH",
        "ADC":     "HIGH",
        "CAR_T":   "HIGH",
        "SIRNA":   "MEDIUM",
        "PROTAC":  "MEDIUM-HIGH",
    }

    rows = []
    sel = user_selection or {}

    for gene in genes:
        qr = qualification_results.get(gene, {})
        tier_results = qr.get("tier_results", [])
        for tr in tier_results:
            modality = tr["modality"]
            gene_sel = sel.get(gene, {})
            is_selected = gene_sel.get("selected_route") == modality

            # entry criteria
            tier = tr["Tier"]
            entry_criteria = ""
            if tier == "HOLD":
                entry_criteria = "필수 항목 UNKNOWN — 추가 근거 확보 후 재평가"
            elif tier == "REJECT":
                entry_criteria = "B0 또는 critical veto — 진입 불가"
            elif tier in ("1A", "2A"):
                entry_criteria = "현 단계 진입 가능 (D1 자산 확인됨)"
            elif tier in ("1B", "2B"):
                entry_criteria = "feasibility 보완 과제 선행 필요"
            else:
                entry_criteria = "탐색 단계 — 고비용 자동 설계 기본 비활성"

            rows.append({
                "gene": gene,
                "modality": modality,
                "bio_confidence": tr.get("B", "B3"),
                "M": tr.get("M", "UNKNOWN"),
                "D": tr.get("D", "D3"),
                "Tier": tier,
                "entry_criteria": entry_criteria,
                "missing_evidence": _val(tr.get("modality_reason")),
                "asset": _val(qr.get("chembl", {}).get("approved_drugs")),
                "estimated_cost": cost_map.get(modality, "UNKNOWN"),
                "next_experiment": next_experiment_map.get(modality, "추가 설계 필요"),
                "execution_status": "SELECTED" if is_selected
                                    else "HOLD" if tier in ("HOLD", "REJECT")
                                    else "PENDING",
                "selected_by": _val(gene_sel.get("selected_by")),
                "selected_at": _val(gene_sel.get("selected_at")),
            })

    return _tsv(rows, keys, Path(output_path))


# ── evidence_cards.jsonl ──────────────────────────────────

def write_evidence_cards_jsonl(
    cards: list[dict],
    output_path: str | Path,
) -> str:
    """evidence_cards.jsonl 작성 (1행 1 Evidence Card)."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for card in cards:
            f.write(json.dumps(card, ensure_ascii=False) + "\n")
    return str(path)


# ── 일괄 생성 ─────────────────────────────────────────────

def write_all_tables(
    genes: list[str],
    qualification_results: dict[str, dict],
    de_evidences: dict[str, dict],
    cell_contexts: dict[str, dict],
    depmap_results: dict[str, dict],
    output_dir: str | Path,
    de_records: list[dict] | None = None,
    all_evidence_cards: list[dict] | None = None,
    user_selection: dict | None = None,
) -> dict:
    """모든 출력 테이블을 output_dir에 저장.

    Returns: {table_name: file_path}
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    paths: dict[str, str] = {}

    paths["gene_overview"] = write_gene_overview(
        genes=genes,
        qualification_results=qualification_results,
        de_evidences=de_evidences,
        cell_contexts=cell_contexts,
        depmap_results=depmap_results,
        output_path=out / "gene_overview.tsv",
        user_selection=user_selection,
    )

    if de_records is not None:
        paths["de_evidence"] = write_de_evidence(
            de_records=de_records,
            output_path=out / "de_evidence.tsv",
        )

    paths["modality_options"] = write_modality_options(
        genes=genes,
        qualification_results=qualification_results,
        output_path=out / "modality_options.tsv",
        user_selection=user_selection,
    )

    if all_evidence_cards is not None:
        paths["evidence_cards"] = write_evidence_cards_jsonl(
            cards=all_evidence_cards,
            output_path=out / "evidence_cards.jsonl",
        )

    # 상태 기록
    status = {
        "status": "COMPLETED",
        "n_genes": len(genes),
        "output_dir": str(out),
        "tables": paths,
    }
    (out / "output_tables_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return paths
