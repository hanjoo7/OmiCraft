"""
OmiCraft Orchestrator — 워크플로우 v5 전체 파이프라인

노드 순서:
  01  resolve_contrasts_node     ContrastSpec / ResearchPlan 확정
  02  r_analysis_node            DESeq2 + GSEA (R subprocess)
  03  apear_node                 aPEAR 경로 네트워크 (R/Python)
  04  cell_context_node          CELLxGENE + genesorteR
  05  depmap_node                DepMap CRISPR/RNAi 의존성
  06  qualification_node         UniProt/HPA/OpenTargets/ChEMBL + Tier 결정
  07  output_tables_node         TSV/JSONL 출력 테이블 생성
  08  user_selection_gate_node   사용자 선택 gate (설계 전 필수)
  09  screening_node             분자 설계 (고비용 — gate 통과 후만 실행)
  --  critic_node                Critic 판정 (전 단계 검토)

고수준 진입점:
  run_screening()               설계 전용 (WF 09, 기존 호환)
  run_qualification_pipeline()  WF 01-08 전체
  run_full_pipeline()           WF 01-09 전체 (gate 통과 필요)
"""

from __future__ import annotations

import csv
import json
import logging
import math
import uuid
from pathlib import Path

from .configuration import get_config
from .automatic_design import automatic_design_node
from .critic_agent import critic_node
from .screening_agent import screening_node
from .small_molecule_io import write_json

log = logging.getLogger(__name__)

_PROCESSED_DIR = Path(__file__).resolve().parent.parent.parent.parent / "dataset" / "tcga_brca" / "processed"
_DEPMAP_DIR = Path(__file__).resolve().parent.parent.parent.parent / "dataset" / "depmap"


def _load_precomputed_r_results() -> dict:
    """RDS 파일이 없을 때 사전 계산된 DEG CSV를 TSV로 변환해 r_analysis_result로 반환."""
    import tempfile

    deg_full_csv = _PROCESSED_DIR / "deg_results" / "deg_tnbc_vs_nontnbc_full.csv"
    deg_sig_csv = _PROCESSED_DIR / "deg_results" / "deg_tnbc_vs_nontnbc_significant.csv"

    if not deg_full_csv.is_file():
        return {}

    def _csv_to_tsv(src: Path, dst: Path) -> None:
        with src.open(newline="", encoding="utf-8") as fin, \
             dst.open("w", newline="", encoding="utf-8") as fout:
            reader = csv.DictReader(fin)
            writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(reader)

    output_dir = Path(tempfile.mkdtemp(prefix="omicraft_precomp_"))
    de_dir = output_dir / "03_de"
    de_dir.mkdir(parents=True)

    de_all_tsv = de_dir / "DE_all_genes.tsv"
    deg_tsv = de_dir / "DEG.tsv"
    _csv_to_tsv(deg_full_csv, de_all_tsv)
    sig_csv = deg_sig_csv if deg_sig_csv.is_file() else deg_full_csv

    # 두 번째 contrast (TNBC vs normal) 로드 → dual_contrast_supported 판별
    normal_sig_csv = _PROCESSED_DIR / "deg_results" / "deg_tnbc_vs_normal_significant.csv"
    normal_sig_genes: dict[str, float] = {}  # gene → LFC
    if normal_sig_csv.is_file():
        with normal_sig_csv.open(newline="", encoding="utf-8") as _nf:
            _nr = csv.DictReader(_nf)
            _nlfc = None
            _nname = None
            for _nrow in _nr:
                if _nlfc is None:
                    _nlfc = next((c for c in _nrow.keys() if "log2" in c.lower() or "lfc" in c.lower()), None)
                    _nname = next((c for c in _nrow.keys() if c.lower() in ("gene_symbol", "gene_name", "gene")), None)
                if not _nlfc or not _nname:
                    continue
                _ng = _nrow.get(_nname, "")
                try:
                    _nl = float(_nrow.get(_nlfc, 0) or 0)
                    _np = float(_nrow.get("padj", 1) or 1)
                    if _ng and _np < 0.05 and abs(_nl) >= 1:
                        normal_sig_genes[_ng] = _nl
                except (ValueError, TypeError):
                    pass

    # DEG TSV 작성 — dual_contrast 열 추가
    with sig_csv.open(newline="", encoding="utf-8") as _sf, \
         deg_tsv.open("w", newline="", encoding="utf-8") as _df:
        _sr = csv.DictReader(_sf)
        _fields = list(_sr.fieldnames or []) + ["dual_contrast"]
        _sw = csv.DictWriter(_df, fieldnames=_fields, delimiter="\t")
        _sw.writeheader()
        for _srow in _sr:
            _sname = next((c for c in _srow.keys() if c.lower() in ("gene_symbol", "gene_name", "gene")), None)
            _slfc = next((c for c in _srow.keys() if "log2" in c.lower() or "lfc" in c.lower()), None)
            _sg = _srow.get(_sname, "") if _sname else ""
            if _sg in normal_sig_genes and _slfc:
                try:
                    _pl = float(_srow.get(_slfc, 0) or 0)
                    _same = (_pl > 0) == (normal_sig_genes[_sg] > 0)
                    _srow["dual_contrast"] = "dual_contrast_supported" if _same else "discordant"
                except (ValueError, TypeError):
                    _srow["dual_contrast"] = "partial_evidence"
            else:
                _srow["dual_contrast"] = "partial_evidence"
            _sw.writerow(_srow)

    # de_summary 계산 (discovery_agent가 사용)
    n_significant = 0
    top_up: list = []
    top_down: list = []
    with sig_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)
    n_significant = len(all_rows)
    lfc_col = next((c for c in (all_rows[0].keys() if all_rows else [])
                    if "log2" in c.lower() or "lfc" in c.lower()), None)
    name_col = next((c for c in (all_rows[0].keys() if all_rows else [])
                     if c.lower() in ("gene_name", "gene_symbol", "symbol")), None)
    if lfc_col and name_col and all_rows:
        try:
            sorted_rows = sorted(all_rows, key=lambda r: float(r.get(lfc_col, 0) or 0), reverse=True)
            top_up = [r[name_col] for r in sorted_rows[:10] if float(r.get(lfc_col, 0) or 0) > 0]
            top_down = [r[name_col] for r in sorted_rows[-10:] if float(r.get(lfc_col, 0) or 0) < 0][::-1]
        except (ValueError, TypeError):
            pass

    return {
        "success": True,
        "execution_mode": "PRECOMPUTED_CACHE",
        "output_dir": str(output_dir),
        "output_files": {
            "de_all": str(de_all_tsv),
            "deg_table": str(deg_tsv),
        },
        "de_summary": {
            "n_significant": n_significant,
            "up_in_tnbc": sum(1 for r in all_rows if lfc_col and float(r.get(lfc_col, 0) or 0) > 0),
            "down_in_tnbc": sum(1 for r in all_rows if lfc_col and float(r.get(lfc_col, 0) or 0) < 0),
            "top_upregulated": top_up,
            "top_downregulated": top_down,
        },
    }


def _load_precomputed_gobp_gsea(contrast_id: str = "tnbc_vs_nontnbc") -> str | None:
    """사전 계산된 GO BP GSEA CSV를 aPEAR가 요구하는 TSV 형식으로 변환.

    Returns: temp directory path containing GOBP_GSEA_{contrast_id}.tsv, or None if unavailable.
    """
    import tempfile

    gobp_csv = _PROCESSED_DIR / "pathway_results" / "gsea_tnbc_vs_nontnbc_GO_Biological_Process_2023" / "gseapy.gene_set.prerank.report.csv"
    if not gobp_csv.is_file():
        return None

    rows = []
    with gobp_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            term = row.get("Term", "")
            # Extract GO ID if present, e.g. "Cilium Assembly (GO:0060271)" -> "GO:0060271"
            go_id = term
            import re as _re
            m = _re.search(r"\(GO:\d+\)", term)
            if m:
                go_id = m.group(0)[1:-1]  # strip parens
            rows.append({
                "ID": go_id,
                "Description": term,
                "NES": row.get("NES", ""),
                "p.adjust": row.get("FDR q-val", "1"),
                "core_enrichment": row.get("Lead_genes", "").replace(";", "/"),
            })

    if not rows:
        return None

    out_dir = Path(tempfile.mkdtemp(prefix="omicraft_gobp_"))
    tsv_path = out_dir / f"GOBP_GSEA_{contrast_id}.tsv"
    with tsv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["ID", "Description", "NES", "p.adjust", "core_enrichment"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return str(out_dir)


# ── WF 01: ContrastSpec 확정 ─────────────────────────────────

def resolve_contrasts_node(state: dict) -> dict:
    """WF 01 — ContrastSpec / ResearchPlan 결정.

    state 입력:
      research_question, disease, subtype (선택)
      sample_manifest (선택 — TCGA 샘플 메타 경로)
      data_flags (선택 — {rnaseq: bool, wes: bool, wgs: bool})
    """
    try:
        from .contrast_spec import resolve_contrasts

        config = get_config()
        output_dir = Path(config.data.agent_results) / "contrast_spec"

        # sample_manifest: state 또는 config에서
        sample_manifest = (
            (state.get("r_analysis_result") or {}).get("input_manifest")
            or state.get("sample_manifest")
            or getattr(config.data, "sample_manifest", None)
        )
        if isinstance(sample_manifest, (str, Path)):
            sample_manifest = json.loads(Path(sample_manifest).read_text())
        data_flags = state.get("data_flags") or {"rnaseq": True}

        plan = resolve_contrasts(
            question=state.get("research_question", ""),
            sample_manifest=sample_manifest,
            data_flags=data_flags,
            plan_id=state.get("run_id", "run001"),
            output_dir=str(output_dir),
        )

        contrasts = [c.to_dict() for c in plan.contrasts]
        return {
            "research_plan": plan.to_dict(),
            "contrasts": contrasts,
            "contrast_manifest": str(output_dir / "contrast_manifest.tsv"),
            "messages": ["[ResolvePlan] ContrastSpec 확정 완료"],
            "current_agent": "planner",
        }
    except Exception as exc:
        log.exception("resolve_contrasts_node 오류")
        return {
            "errors": [f"resolve_contrasts_node: {exc}"],
            "messages": [f"[ResolvePlan] 오류: {exc}"],
            "current_agent": "planner",
        }


# ── WF 02-03: DESeq2 + GSEA + aPEAR ────────────────────────

def r_analysis_node(state: dict) -> dict:
    """WF 02-03 — R DESeq2 / GSEA 실행."""
    try:
        from .r_tool import run_r_analysis

        config = get_config()
        if state.get("dry_run"):
            return {"r_analysis_result": {"success": False, "execution_status": "NOT_RUN"},
                    "execution_status": "NOT_RUN", "messages": ["[R-Analysis] Dry-run: 분석 미실행"],
                    "current_agent": "discovery"}
        output_dir = Path(config.data.agent_results) / "r_analysis" / uuid.uuid4().hex
        inputs = {name: state.get(name) or getattr(config.data, name, None)
                  for name in ("counts_rds", "metadata_rds", "annotation_rds", "clinical_rds", "gene_sets_rds")}
        inputs["counts_rds"] = inputs["counts_rds"] or state.get("count_matrix_path")
        inputs["metadata_rds"] = inputs["metadata_rds"] or state.get("metadata_path")
        result = run_r_analysis(**inputs, output_dir=str(output_dir), progress_callback=state.get("_progress_callback"),
                                reuse_results=state.get("reuse_results", True))
        if not result.get("success"):
            reason = result.get("error", "R analysis did not complete")
            # RDS 파일 누락 시 사전 계산된 DEG 결과로 폴백
            if result.get("execution_status") == "BLOCKED_INPUT" and state.get("reuse_results", True):
                precomputed = _load_precomputed_r_results()
                if precomputed:
                    log.info("r_analysis_node: RDS 파일 없음 → 사전 계산된 DEG 결과 로드")
                    return {"r_analysis_result": precomputed,
                            "execution_status": "COMPLETED",
                            "messages": ["[R-Analysis] 사전 계산된 DEG 결과 로드 (PRECOMPUTED_CACHE)"],
                            "current_agent": "discovery"}
            return {"r_analysis_result": result,
                    "execution_status": result.get("execution_status", "FAILED"),
                    "blockers": result.get("blockers", [reason]),
                    "errors": [f"r_analysis_node: {reason}"],
                    "messages": [f"[R-Analysis] 분석 중단: {reason}"],
                    "current_agent": "discovery"}
        if not (result.get("output_files") or {}).get("de_all"):
            reason = "DE_all_genes.tsv was not produced; downstream analysis cannot start"
            return {"r_analysis_result": {**result, "success": False},
                    "execution_status": "FAILED", "errors": [f"r_analysis_node: {reason}"],
                    "messages": [f"[R-Analysis] {reason}"], "current_agent": "discovery"}
        plan_update = {}
        if result.get('input_manifest'):
            plan_update = resolve_contrasts_node({**state, 'r_analysis_result': result})
            plan_update['sample_manifest'] = result['input_manifest']
        return {**plan_update, "r_analysis_result": result, "execution_status": "COMPLETED",
                "messages": ["[R-Analysis] DESeq2/GSEA 완료 · " + result.get("execution_mode", "FRESH_ANALYSIS")], "current_agent": "discovery"}

    except Exception as exc:
        log.exception("r_analysis_node 오류")
        return {
            "r_analysis_result": {"success": False},
            "execution_status": "FAILED",
            "errors": [f"r_analysis_node: {exc}"],
            "messages": [f"[R-Analysis] 오류: {exc}"],
            "current_agent": "discovery",
        }


def apear_node(state: dict) -> dict:
    """WF 03 — aPEAR 경로 네트워크."""
    try:
        from .apear_adapter import run_apear_for_contrasts

        config = get_config()
        r_result = state.get("r_analysis_result") or {}
        if r_result.get("success") is False:
            return {"apear_results": {}, "execution_status": "NOT_RUN",
                    "messages": ["[aPEAR] R 분석 결과 없음 — 실행하지 않음"],
                    "current_agent": "discovery"}
        gsea_output_dir = str(Path(r_result.get("output_dir") or
                                  Path(config.data.agent_results) / "r_analysis") / "05_gsea")
        output_dir = str(Path(config.data.agent_results) / "apear")
        contrast_ids = [
            c.get("contrast_id") if isinstance(c, dict) else getattr(c, "contrast_id", None)
            for c in state.get("contrasts", [])
        ]
        contrast_ids = [c for c in contrast_ids if c == "tnbc_vs_nontnbc"]

        results = run_apear_for_contrasts(
            gsea_output_dir=gsea_output_dir,
            output_dir=output_dir,
            contrast_ids=contrast_ids or ["tnbc_vs_nontnbc"],
        )
        if results.get("_error"):
            # GOBP GSEA TSV 없음 → 사전 계산된 GO BP CSV로 폴백
            precomp_gsea_dir = _load_precomputed_gobp_gsea(
                contrast_ids[0] if contrast_ids else "tnbc_vs_nontnbc"
            )
            if precomp_gsea_dir:
                log.info("apear_node: GOBP GSEA TSV 없음 → 사전 계산된 GO BP 결과 사용")
                results = run_apear_for_contrasts(
                    gsea_output_dir=precomp_gsea_dir,
                    output_dir=output_dir,
                    contrast_ids=contrast_ids or ["tnbc_vs_nontnbc"],
                )
        if results.get("_error"):
            reason = str(results["_error"])
            return {"apear_results": results, "execution_status": "PARTIAL",
                    "errors": [f"apear_node: {reason}"],
                    "messages": [f"[aPEAR] 실행 불가: {reason}"], "current_agent": "discovery"}
        records = [item for item in results.values() if isinstance(item, dict)]
        n_total = sum(item.get("n_nodes", 0) for item in records)
        issues = [str(item.get("reason") or item.get("status")) for item in records
                  if item.get("status") in {"FAILED", "UNAVAILABLE", "INVALID_INPUT"}]
        return {
            "apear_results": results,
            "execution_status": "PARTIAL" if issues else "COMPLETED",
            "errors": [f"apear_node: {reason}" for reason in issues],
            "messages": [f"[aPEAR] 경로 네트워크 결과 — 총 {n_total}개 노드", *issues],
            "current_agent": "discovery",
        }

    except Exception as exc:
        log.exception("apear_node 오류")
        return {
            "errors": [f"apear_node: {exc}"],
            "messages": [f"[aPEAR] 오류: {exc}"],
            "current_agent": "discovery",
        }


# ── WF 04-05-A: Cell context ─────────────────────────────────

def cell_context_node(state: dict) -> dict:
    """WF 04-05-A — CELLxGENE Census + genesorteR."""
    try:
        from .cellxgene_adapter import run_cell_context_analysis

        config = get_config()
        output_dir = str(Path(config.data.agent_results) / "cell_context")

        genes = _extract_gene_list(state)
        if not genes:
            return {
                "messages": ["[CellContext] 유전자 목록 없음 — 건너뜀"],
                "current_agent": "qualification",
            }

        expr_rds = getattr(config.data, "expr_rds", None) or state.get("expr_rds")
        meta_rds = getattr(config.data, "meta_rds", None) or state.get("meta_rds")

        result = run_cell_context_analysis(
            genes=genes,
            output_dir=output_dir,
            expr_rds=expr_rds,
            meta_rds=meta_rds,
            census_tissue=state.get("tissue", "breast"),
            census_disease=state.get("census_disease", "breast carcinoma"),
        )

        cell_contexts = result.get("context_results", {})
        return {
            "cell_context_result": result,
            "execution_status": "COMPLETED" if result.get("status") == "COMPLETED" else "PARTIAL",
            "cell_contexts": cell_contexts,
            "cell_context_results": cell_contexts,  # 하위호환
            "messages": [
                f"[CellContext] {len(cell_contexts)}개 유전자 분석 완료 "
                f"(status={result.get('status', '?')})"
            ],
            "current_agent": "qualification",
        }
    except Exception as exc:
        log.exception("cell_context_node 오류")
        return {
            "errors": [f"cell_context_node: {exc}"],
            "messages": [f"[CellContext] 오류: {exc}"],
            "current_agent": "qualification",
        }


# ── WF 05-B: DepMap ──────────────────────────────────────────

def depmap_node(state: dict) -> dict:
    """WF 05-B — DepMap CRISPR/RNAi 의존성 평가."""
    try:
        from .depmap_module import run_depmap_evaluation

        config = get_config()
        output_dir = str(Path(config.data.agent_results) / "depmap")

        genes = _extract_gene_list(state)
        if not genes:
            return {
                "messages": ["[DepMap] 유전자 목록 없음 — 건너뜀"],
                "current_agent": "qualification",
            }

        chronos_path = (
            getattr(config.data, "chronos_path", None) or state.get("chronos_path")
            or str(_DEPMAP_DIR / "CRISPRGeneEffect.csv")
        )
        demeter2_path = (
            getattr(config.data, "demeter2_path", None) or state.get("demeter2_path")
        )
        model_meta_path = (
            getattr(config.data, "model_meta_path", None) or state.get("model_meta_path")
            or str(_DEPMAP_DIR / "Model.csv")
        )
        # local 파일이 없으면 None으로 초기화
        if chronos_path and not Path(chronos_path).is_file():
            chronos_path = None
        if model_meta_path and not Path(model_meta_path).is_file():
            model_meta_path = None

        result = run_depmap_evaluation(
            genes=genes,
            output_dir=output_dir,
            chronos_path=chronos_path,
            demeter2_path=demeter2_path,
            model_meta_path=model_meta_path,
            depmap_release=config.upstream.depmap_release,
        )

        depmap_results = result.get("results", {})
        return {
            "depmap_result": result,
            "execution_status": "COMPLETED" if result.get("status") == "COMPLETED" else "PARTIAL",
            "depmap_results": depmap_results,
            "messages": [
                f"[DepMap] {len(depmap_results)}개 유전자 평가 완료"
            ],
            "current_agent": "qualification",
        }
    except Exception as exc:
        log.exception("depmap_node 오류")
        return {
            "errors": [f"depmap_node: {exc}"],
            "messages": [f"[DepMap] 오류: {exc}"],
            "current_agent": "qualification",
        }


# ── WF 06-07: Qualification + 출력 테이블 ──────────────────

def qualification_node(state: dict) -> dict:
    """WF 06-07 — UniProt/HPA/OpenTargets/ChEMBL + Tier 결정."""
    try:
        from .target_qualification import run_qualification_pipeline

        config = get_config()
        output_dir = str(Path(config.data.agent_results) / "qualification")

        genes = _extract_gene_list(state)
        if not genes:
            return {
                "messages": ["[Qualification] 유전자 목록 없음 — 건너뜀"],
                "current_agent": "qualification",
            }

        de_evidences = _build_de_evidences(state)
        cell_contexts = state.get("cell_contexts") or state.get("cell_context_results") or {}
        depmap_results = state.get("depmap_results") or {}

        result = run_qualification_pipeline(
            genes=genes,
            de_evidences=de_evidences,
            cell_contexts=cell_contexts,
            depmap_results=depmap_results,
            output_dir=output_dir,
            offline=state.get("dry_run", False),
        )

        qr = result.get("results", {})
        def targets(names):
            output = []
            for gene in names:
                record = qr[gene]
                best = record.get("best_tier", "HOLD")
                routes = [row for row in record.get("tier_results", []) if row.get("Tier") == best]
                output.append({"gene_name": gene, "candidate_id": gene, "tier": best,
                               "bio_confidence": record.get("bio_confidence", "UNKNOWN"),
                               "modality": routes[0]["modality"] if routes else "UNKNOWN",
                               "tier_results": record.get("tier_results", []),
                               "shortlist": record.get("shortlist", {}),
                               "safety_veto": record.get("safety_veto", {})})
            return output
        advance = targets(result.get("advance", []))
        hold = targets(result.get("hold", []))
        reject = targets(result.get("reject", []))

        search = {'requested_limit': config.upstream.max_candidates, 'evaluated_targets': len(qr),
                  'advance_targets': len(advance), 'hold_targets': len(hold), 'reject_targets': len(reject),
                  'selection': 'explicit_gene_list' if state.get('gene_list') else 'DE_ranked_candidates'}
        msgs = [
            f"[Qualification] 표적 {len(qr)}개 평가 완료 (설정 상한 {config.upstream.max_candidates}개)",
            f"[Qualification] 완료 — ADVANCE={len(advance)}, HOLD={len(hold)}, REJECT={len(reject)}"
        ]
        for t in advance[:5]:
            msgs.append(
                f"  ADVANCE: {t.get('gene_name')} ({t.get('modality')}, "
                f"Tier={t.get('tier')}, B={t.get('bio_confidence')})"
            )

        return {
            "qualification_results": qr,
            "candidate_search": search,
            "advance_targets": advance,
            "hold_targets": hold,
            "reject_targets": reject,
            "qualification_results_path": str(Path(output_dir) / "qualification_results.json"),
            "evidence_cards": [*(state.get("evidence_cards") or []), *result.get("all_evidence_cards", [])],
            "execution_status": "COMPLETED",
            "messages": msgs,
            "current_agent": "qualification",
        }
    except Exception as exc:
        log.exception("qualification_node 오류")
        return {
            "errors": [f"qualification_node: {exc}"],
            "messages": [f"[Qualification] 오류: {exc}"],
            "current_agent": "qualification",
        }


def output_tables_node(state: dict) -> dict:
    """WF §4 — TSV/JSONL 출력 테이블 생성."""
    try:
        from .output_tables import de_records_from_r_output, write_all_tables

        config = get_config()
        output_dir = Path(config.data.agent_results) / "output_tables"

        genes = _extract_gene_list(state)
        qr = state.get("qualification_results") or {}
        cell_contexts = state.get("cell_contexts") or state.get("cell_context_results") or {}
        depmap_results = state.get("depmap_results") or {}

        # DE 요약 정보 구성
        de_evidences: dict[str, dict] = {}
        r_result = state.get("r_analysis_result") or {}
        de_records = de_records_from_r_output(r_result) if r_result else None

        # gene별 de_evidence 요약
        for gene in genes:
            de_evidences[gene] = _gene_de_summary(gene, r_result)

        # 사용자 선택 정보
        user_selection: dict[str, dict] = {}
        for ev in state.get("user_selection_events", []):
            gene = ev.get("gene", "")
            if gene:
                user_selection[gene] = ev

        paths = write_all_tables(
            genes=genes,
            qualification_results=qr,
            de_evidences=de_evidences,
            cell_contexts=cell_contexts,
            depmap_results=depmap_results,
            output_dir=output_dir,
            de_records=de_records,
            all_evidence_cards=state.get("evidence_cards"),
            user_selection=user_selection,
        )

        return {
            "output_tables": paths,
            "gene_overview_path": paths.get("gene_overview", ""),
            "de_evidence_path": paths.get("de_evidence", ""),
            "modality_options_path": paths.get("modality_options", ""),
            "evidence_cards_path": paths.get("evidence_cards", ""),
            "messages": [
                f"[OutputTables] {len(paths)}개 테이블 생성 완료",
                *[f"  {k}: {v}" for k, v in paths.items()],
            ],
            "current_agent": "qualification",
        }
    except Exception as exc:
        log.exception("output_tables_node 오류")
        return {
            "errors": [f"output_tables_node: {exc}"],
            "messages": [f"[OutputTables] 오류: {exc}"],
            "current_agent": "qualification",
        }


# ── WF 08: User Selection Gate ──────────────────────────────

def user_selection_gate_node(state: dict) -> dict:
    """WF 08 — 사용자 선택 gate (user_selection.py에서 위임)."""
    from .user_selection import user_selection_gate_node as _gate_node
    if state.get('design_mode') == 'all_eligible' and not state.get('dry_run'):
        from .automatic_design import select_routes
        events = select_routes(state)
        return {**_gate_node({**state, 'user_selection_events': events}), 'user_selection_events': events}
    return _gate_node(state)


# ── WF 09: 분자 설계 (기존 호환) ─────────────────────────────

def assess_modalities(state):
    from .modality_dispatch import readiness
    from .routes import MODALITY_ROUTES
    from .user_selection import check_design_gate

    assessments = []
    for target in state.get("advance_targets", []):
        candidate = target.get("candidate_id") or target.get("gene_name", "")
        tier = target.get("tier") or target.get("best_tier") or target.get("Tier", "UNKNOWN")
        for modality in MODALITY_ROUTES:
            route = next((r for r in target.get("tier_results", []) if r.get("modality") == modality), None)
            route_tier = route["Tier"] if route else tier
            key = f"{candidate}:{modality}"
            data = state.get("execution_inputs", {}).get(key,
                state.get("execution_inputs", {}).get(modality, {}))
            gate = check_design_gate(candidate, modality, route_tier, state.get("user_selection_events", []))
            if target.get("safety_verdict") in {"VETO", "REJECT"} or target.get("safety_veto", {}).get(modality, {}).get("has_veto"):
                gate = {"allowed": False, "gate_status": "BLOCKED_REJECT", "reason": "target_safety_veto"}
            ready = (readiness(modality, data, state.get("modality_config", {}).get(modality, {}))
                     if data else {"tool_status": "not_available", "blockers": ["explicit_modality_input_missing"]})
            assessments.append({"candidate": candidate, "modality": modality, "tier": route_tier,
                                "gate": gate, "readiness": ready})
    return {"modality_assessments": assessments, "current_agent": "therapeutic_design"}


def run_screening(state: dict) -> dict:
    """Qualification → assessment → explicit user selection → executor → validation/critic."""
    from .modality_dispatch import execute
    from .routes import MODALITY_ROUTES
    from .run_records import RunRecord

    if state.get("execution_profile", get_config().execution_profile) == "competition_demo":
        from .configuration import CompetitionDemoConfig
        from .modality_dispatch import run_competition_demo
        config = CompetitionDemoConfig(**state.get("competition_demo", get_config().competition_demo.model_dump()))
        result = run_competition_demo(state, config)
        return {**state, **result, "modality_results": result["summary"]["modalities"], "current_agent": "critic"}

    output = dict(state)
    output.update(assess_modalities(state))
    results = []
    root = Path(state.get("output_dir") or Path(get_config().data.agent_results) / "therapeutic")
    for assessment in output["modality_assessments"]:
        candidate, modality = assessment["candidate"], assessment["modality"]
        gate = assessment["gate"]
        if not gate["allowed"]:
            continue
        key = f"{candidate}:{modality}"
        data = state.get("execution_inputs", {}).get(key,
            state.get("execution_inputs", {}).get(modality, state.get("modality_input", {})))
        if modality in {"SMALL_MOLECULE", "DE_NOVO_BINDER"} and not data:
            data = {**state, "route": MODALITY_ROUTES[modality]}
        import re
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", candidate) or "candidate"
        cfg = {**state.get("modality_config", {}).get(modality, {}),
               "output_dir": str(root / f"{safe}_{modality.lower()}"),
               "dry_run": state.get("dry_run", True)}
        selection = next(ev for ev in reversed(state['user_selection_events'])
                         if ev.get('route_id') == key and ev.get('selected_by') == 'user'
                         and ev.get('status') == 'SELECTED')
        item = execute({"modality": modality, "input": data, "selection_event": selection}, cfg)
        results.append(item)
    if not results:
        from .modality_dispatch import result
        record = RunRecord(root, "ASSESSMENT")
        item = result("SMALL_MOLECULE", "BLOCKED", "user_selection", blockers=["explicit_user_selection_required"])
        item["modality"] = "ASSESSMENT"
        item["summary"] = {"assessments": output["modality_assessments"]}
        record.finish(item)
    output.update(critic_verdict="APPROVE" if any(r.get("critic", {}).get("decision") == "ADVANCE" for r in results) else "HOLD",
                  critic_result={"accepted_candidate_ids": [], "verdict": "HOLD"},
                  modality_results=results, gate_status="ALLOWED" if results else "PENDING",
                  execution_status="COMPLETED" if results and all(r['execution_status'] == 'COMPLETED' for r in results)
                  else "PARTIAL" if any(r['execution_status'] in {'COMPLETED', 'PARTIAL'} for r in results) else "BLOCKED",
                  current_agent="critic")
    return output


# ── 전체 파이프라인 ──────────────────────────────────────────

def run_qualification_pipeline(state: dict) -> dict:
    """WF 01-08 — ContrastSpec부터 User Selection Gate까지 전체 실행.

    각 노드가 state를 업데이트하며 순서대로 실행된다.
    필수 R 분석이 실패하면 후속 분석을 실행하지 않는다.
    """
    s = dict(state)
    messages: list[str] = []
    errors: list[str] = []

    def _merge(node_fn):
        nonlocal s
        out = node_fn(s)
        # messages/errors 누적
        messages.extend(out.pop("messages", []))
        errors.extend(out.pop("errors", []))
        s.update(out)
        # agent reviews are handled within each agent node in the LangGraph graph

    log.info("[Orchestrator] 파이프라인 시작 — WF 01~08")

    _merge(resolve_contrasts_node)
    _merge(r_analysis_node)
    if not (s.get("r_analysis_result") or {}).get("success"):
        s["messages"] = (state.get("messages") or []) + messages
        s["errors"] = (state.get("errors") or []) + errors
        return s
    _merge(apear_node)
    _merge(cell_context_node)
    _merge(depmap_node)
    _merge(qualification_node)
    _merge(output_tables_node)
    _merge(assess_modalities)
    _merge(user_selection_gate_node)
    if s.get('design_mode') == 'all_eligible' and not s.get('dry_run'):
        _merge(automatic_design_node)

    s["messages"] = (state.get("messages") or []) + messages
    s["errors"] = (state.get("errors") or []) + errors
    s["current_agent"] = "orchestrator"

    log.info(
        "[Orchestrator] 파이프라인 완료 — gate_status=%s, ADVANCE=%d",
        s.get("gate_status"),
        len(s.get("advance_targets") or []),
    )
    return s


def run_full_pipeline(state: dict) -> dict:
    """WF 01-09 — 전체 파이프라인 (gate 통과 시 설계까지 포함).

    WF 01-08 완료 후 gate_status가 ALLOWED이면 WF 09 설계를 실행한다.
    """
    s = run_qualification_pipeline(state)

    if s.get("design_mode") == "all_eligible":
        return s
    if any(a["gate"]["allowed"] for a in s.get("modality_assessments", [])):
        log.info("[Orchestrator] gate ALLOWED — 설계 단계 실행")
        s = run_screening(s)
    else:
        s.setdefault("messages", []).append(
            f"[Orchestrator] 설계 대기 — gate_status={s.get('gate_status', 'PENDING')}. "
            "apply_user_selection() 호출 후 run_screening()을 실행하세요."
        )

    return s


# ── 헬퍼 ─────────────────────────────────────────────────────

def _de_gene_records(r_result):
    path = (r_result.get("output_files") or {}).get("deg_table")
    if not path or not Path(path).is_file():
        return {}
    def number(value):
        try:
            parsed = float(value)
            return parsed if math.isfinite(parsed) else None
        except (ValueError, TypeError):
            return None
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    records = {}
    for row in rows:
        gene = row.get("gene_symbol") or row.get("gene_name") or row.get("gene")
        if not gene or gene in {"NA", "nan"} or gene in records:
            continue
        padj, lfc = number(row.get("padj")), number(row.get("log2FoldChange"))
        if padj is None or lfc is None or padj >= 0.05 or abs(lfc) < 1:
            continue
        dual = row.get("dual_contrast") or "partial_evidence"
        records[gene] = {"gene_id": row.get("gene_id"), "log2FC_tnbc_nontnbc": lfc,
                         "padj_tnbc_nontnbc": padj, "dual_contrast": dual,
                         "log2FC_tnbc_normal": None, "padj_tnbc_normal": None,
                         "prevalence": None, "source_path": str(Path(path).resolve())}
    return records


def _extract_gene_list(state):
    if state.get("gene_list"):
        return list(dict.fromkeys(state["gene_list"]))
    records = _de_gene_records(state.get("r_analysis_result") or {})
    budget = get_config().upstream.max_candidates
    # TNBC 상향조절(LFC > 0) 유전자를 우선 선정 — 치료 표적 적합성 높음
    # 예산 70%는 TNBC 과발현 유전자, 30%는 하향조절(종양억제인자 등)로 채움
    up = sorted(
        [(g, r) for g, r in records.items() if (r.get("log2FC_tnbc_nontnbc") or 0) > 0],
        key=lambda x: (x[1].get("padj_tnbc_nontnbc") or 1, -(x[1].get("log2FC_tnbc_nontnbc") or 0)),
    )
    down = sorted(
        [(g, r) for g, r in records.items() if (r.get("log2FC_tnbc_nontnbc") or 0) <= 0],
        key=lambda x: (x[1].get("padj_tnbc_nontnbc") or 1, (x[1].get("log2FC_tnbc_nontnbc") or 0)),
    )
    n_up = min(len(up), math.ceil(budget * 0.7))
    n_down = min(len(down), budget - n_up)
    selected = [g for g, _ in up[:n_up]] + [g for g, _ in down[:n_down]]
    return selected[:budget]


def _gene_de_summary(gene, r_result):
    return _de_gene_records(r_result).get(gene, {"dual_contrast": "partial_evidence", "prevalence": None})


def _build_de_evidences(state):
    records = _de_gene_records(state.get("r_analysis_result") or {})
    return {gene: records.get(gene, {"dual_contrast": "partial_evidence", "prevalence": None})
            for gene in _extract_gene_list(state)}


def build_graph(checkpointer=None):
    """LangGraph StateGraph를 반환한다. graph.py로 위임."""
    from .graph import build_graph as _build
    return _build(checkpointer=checkpointer)


def _run_screening_cli(state):
    """Existing direct Screening CLI; agent-driven calls use run_screening instead."""
    result = {**state, **screening_node(state)}
    result.update(critic_node(result))
    path = Path(get_config().data.agent_results) / "final_dossier.json"
    write_json(path, {"route": result.get("route"), "screening_result": result.get("screening_result", {}),
                     "protein_design_result": result.get("protein_design_result", {}),
                     "critic": result.get("critic_result"), "evidence_cards": result.get("evidence_cards", [])})
    result["dossier_path"] = str(path)
    return result
