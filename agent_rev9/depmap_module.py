"""
DepMap 의존성 평가 — 워크플로우 v5 Section 05-B

CRISPR Chronos gene effect를 1차 근거, RNAi DEMETER2를 보강 근거로 사용.

결과 상태:
  CRISPR_ONLY    — CRISPR만 있음
  CONCORDANT     — CRISPR + RNAi 일치
  DISCORDANT     — CRISPR + RNAi 불일치 (원인 검토 필요, 자동 탈락 사유 아님)
  RNAI_ONLY      — RNAi만 있음
  UNAVAILABLE    — 데이터 없음 / 해당 없음

주의사항 (워크플로우 §05-B):
  - Chronos −0.5 등의 cutoff는 탐색 설정값; 절대 target cutoff로 사용 불가.
  - RNAi 데이터 없음만으로 감점하지 않는다.
  - Common essential과 종양 선택적 dependency를 구분한다.
  - 정상조직 안전성을 DepMap만으로 입증하지 않는다.
  - ADC/CAR-T 세포 표적화 modality에서 비필수성을 탈락 근거로 전이하지 않는다.
"""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Optional

# ── 의존성 임계값 (탐색 설정값 — 절대 cutoff 아님) ─────────
_CHRONOS_DEPENDENCY_THRESHOLD = -0.5   # 탐색 참고용
_COMMON_ESSENTIAL_THRESHOLD   = -0.9   # common essential 탐색 기준

DependencyStatus = str  # CRISPR_ONLY | CONCORDANT | DISCORDANT | RNAI_ONLY | UNAVAILABLE


def _parse_float(val: str | None) -> Optional[float]:
    if val is None or str(val).strip() in ("", "NA", "NaN", "nan"):
        return None
    try:
        return float(val)
    except ValueError:
        return None


# ── 파일 로더 ─────────────────────────────────────────────

def load_chronos(path: str | Path, genes=None, model_aliases=None) -> dict[str, dict[str, float]]:
    """Read official model-major or gene-major effect matrices without losing IDs."""
    p = Path(path)
    if not p.is_file():
        return {}
    wanted = set(genes) if genes is not None else None
    aliases = model_aliases or {}
    def symbol(label):
        return re.sub(r" \(\d+\)$", "", label.strip())
    def model(label):
        return aliases.get(label, label)
    result = {}
    with p.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        first = next(reader, None)
        if not first or len(header) < 2:
            return result
        import itertools
        model_major = header[0].lower() in {"modelid", "depmap_id", "cell_line_name"} or first[0].startswith("ACH-")
        columns = [(i, symbol(g)) for i, g in enumerate(header[1:], 1)
                   if wanted is None or symbol(g) in wanted] if model_major else []
        for row in itertools.chain([first], reader):
            if len(row) != len(header):
                raise ValueError("Malformed dependency matrix row: " + row[0])
            values = ((gene, model(row[0]), row[i]) for i, gene in columns) if model_major else (
                (symbol(row[0]), model(mid), value) for mid, value in zip(header[1:], row[1:]))
            if not model_major and wanted is not None and symbol(row[0]) not in wanted:
                continue
            for gene, mid, value in values:
                effect = _parse_float(value)
                if effect is not None and math.isfinite(effect):
                    if mid in result.get(gene, {}):
                        raise ValueError("Duplicate gene/model after ID mapping: " + gene + "/" + mid)
                    result.setdefault(gene, {})[mid] = effect
    return result


def load_demeter2(path: str | Path, genes=None, model_aliases=None) -> dict[str, dict[str, float]]:
    return load_chronos(path, genes, model_aliases)


def load_model_metadata(path: str | Path) -> dict[str, dict]:
    """Model metadata CSV를 로드.

    Returns: {model_id: {lineage, cancer_type, sample_collection_site, ...}}
    """
    p = Path(path)
    if not p.is_file():
        return {}
    result: dict[str, dict] = {}
    with p.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mid = (row.get("ModelID") or row.get("DepMap_ID") or
                   row.get("CCLE_Name") or "")
            if mid:
                result[mid] = dict(row)
    return result


# ── 유방암/TNBC 모델 필터 ─────────────────────────────────

def filter_breast_models(
    model_meta: dict[str, dict],
    cancer_type_col: str = "OncotreeLineage",
    tnbc_label: str | None = None,
) -> dict[str, str]:
    """유방암 모델을 필터링.

    Returns: {model_id: subtype_label}
    subtype_label: "TNBC" | "ER+" | "HER2+" | "breast_other" | "breast_unknown"
    """
    breast_models: dict[str, str] = {}
    for mid, meta in model_meta.items():
        lineage = (meta.get(cancer_type_col) or
                   meta.get("OncotreePrimaryDisease") or
                   meta.get("cancer_type") or "").lower()
        if "breast" not in lineage:
            continue

        er  = (meta.get("ER_status") or meta.get("ER") or "").lower()
        pr  = (meta.get("PR_status") or meta.get("PR") or "").lower()
        her = (meta.get("HER2_status") or meta.get("HER2") or "").lower()

        if all(value in ("negative", "neg", "0") for value in (er, pr, her)):
            subtype = "TNBC"
        elif her in ("positive", "pos", "3+", "amplified"):
            subtype = "HER2+"
        elif er in ("positive", "pos"):
            subtype = "ER+"
        else:
            subtype = "breast_other"

        breast_models[mid] = subtype
    return breast_models


# ── 단일 유전자 의존성 평가 ───────────────────────────────

def evaluate_dependency(
    gene: str,
    chronos_data: dict[str, dict[str, float]],
    breast_models: dict[str, str],
    demeter2_data: dict[str, dict[str, float]] | None = None,
    chronos_threshold: float = _CHRONOS_DEPENDENCY_THRESHOLD,
) -> dict:
    """단일 유전자의 CRISPR + RNAi 의존성을 평가.

    Args:
        gene: 평가할 유전자명.
        chronos_data: load_chronos() 결과.
        breast_models: filter_breast_models() 결과.
        demeter2_data: load_demeter2() 결과 (없으면 None).
        chronos_threshold: 탐색 기준 threshold (참고용).

    Returns:
        {
          "gene": str,
          "dependency_status": DependencyStatus,
          "crispr_summary": dict,
          "rnai_summary": dict | None,
          "concordance_note": str,
          "applicable": bool,
          "applicability_reason": str,
          "evidence_ids": list[str],
        }
    """
    gene_effects_crispr = chronos_data.get(gene, {})
    gene_effects_rnai   = (demeter2_data or {}).get(gene, {})

    # 유방암 모델만 필터
    breast_crispr = {mid: eff for mid, eff in gene_effects_crispr.items()
                     if mid in breast_models}
    breast_rnai   = {mid: eff for mid, eff in gene_effects_rnai.items()
                     if mid in breast_models}
    tnbc_crispr   = {mid: eff for mid, eff in breast_crispr.items()
                     if breast_models.get(mid) == "TNBC"}
    tnbc_rnai     = {mid: eff for mid, eff in breast_rnai.items()
                     if breast_models.get(mid) == "TNBC"}

    # CRISPR 요약
    crispr_summary = _summarize_effects(gene, tnbc_crispr, breast_crispr,
                                        threshold=chronos_threshold,
                                        technology="CRISPR_Chronos")

    # RNAi 요약
    rnai_summary: dict | None = None
    if gene_effects_rnai:
        rnai_summary = _summarize_effects(gene, tnbc_rnai, breast_rnai,
                                          threshold=chronos_threshold,
                                          technology="RNAi_DEMETER2")

    # 상태 결정
    has_crispr = crispr_summary["n_models"] > 0
    has_rnai   = rnai_summary is not None and rnai_summary["n_models"] > 0

    if has_crispr and has_rnai:
        # 방향 일치 여부 판단 (탐색 기준으로만)
        crispr_dep = crispr_summary["median_effect"] is not None and \
                     crispr_summary["median_effect"] < chronos_threshold
        rnai_dep   = rnai_summary["median_effect"] is not None and \
                     rnai_summary["median_effect"] < chronos_threshold
        if crispr_dep == rnai_dep:
            dep_status: DependencyStatus = "CONCORDANT"
            note = ("CRISPR + RNAi 모두 의존성 탐색 기준 충족" if crispr_dep else
                    "CRISPR + RNAi 모두 의존성 탐색 기준 미충족 — 기능 지지 아님")
        else:
            dep_status = "DISCORDANT"
            note = (
                "CRISPR ↔ RNAi 불일치 — "
                f"CRISPR median={crispr_summary['median_effect']:.3f}, "
                f"RNAi median={rnai_summary['median_effect']:.3f}. "
                "Knockout vs knockdown 차이, screen 품질, legacy release 검토 필요. "
                "불일치는 자동 탈락 사유가 아님."
            )
    elif has_crispr:
        dep_status = "CRISPR_ONLY"
        note = "RNAi 데이터 없음 — CRISPR 단독. RNAi 없음이 감점 사유 아님."
    elif has_rnai:
        dep_status = "RNAI_ONLY"
        note = "CRISPR 없음 — RNAi만 있음. 보강 근거로만 사용."
    else:
        dep_status = "UNAVAILABLE"
        note = f"유방암/TNBC 모델에서 {gene} 의존성 데이터 없음"

    # common essential 여부 (종양 선택적 dependency와 구분)
    common_essential_note = ""
    if crispr_summary.get("median_effect") is not None:
        if crispr_summary["median_effect"] < _COMMON_ESSENTIAL_THRESHOLD:
            common_essential_note = (
                f"median effect {crispr_summary['median_effect']:.3f} < "
                f"{_COMMON_ESSENTIAL_THRESHOLD} — common essential 가능성 검토 필요"
            )

    return {
        "gene": gene,
        "dependency_status": dep_status,
        "crispr_summary": crispr_summary,
        "rnai_summary": rnai_summary,
        "concordance_note": note,
        "common_essential_note": common_essential_note,
        "applicable": has_crispr or has_rnai,
        "applicability_reason": note,
        "evidence_ids": _make_evidence_ids(gene, crispr_summary, rnai_summary),
    }


def _summarize_effects(
    gene: str,
    tnbc_effects: dict[str, float],
    all_breast_effects: dict[str, float],
    threshold: float,
    technology: str,
) -> dict:
    all_vals = list(all_breast_effects.values())
    vals = list(tnbc_effects.values()) or all_vals
    cohort_scope = "TNBC" if tnbc_effects else "BREAST_LINEAGE"

    if not vals:
        return {
            "technology": technology,
            "n_models": 0,
            "n_tnbc_models": 0,
            "median_effect": None,
            "mean_effect": None,
            "dependency_fraction": None,
            "min_effect": None,
            "max_effect": None,
        }

    vals_sorted = sorted(vals)
    n = len(vals)
    median = vals_sorted[n // 2] if n % 2 == 1 else \
             (vals_sorted[n // 2 - 1] + vals_sorted[n // 2]) / 2
    mean = sum(vals) / n
    dep_fraction = sum(1 for v in vals if v < threshold) / n

    return {
        "technology": technology,
        "n_models": n,
        "n_tnbc_models": len(tnbc_effects),
        "cohort_scope": cohort_scope,
        "n_breast_models": len(all_vals),
        "median_effect": round(median, 4),
        "mean_effect": round(mean, 4),
        "dependency_fraction": round(dep_fraction, 4),
        "dependency_threshold_used": threshold,
        "dependency_threshold_note": "탐색 설정값 — 절대 target cutoff 아님",
        "min_effect": round(min(vals), 4),
        "max_effect": round(max(vals), 4),
    }


def _make_evidence_ids(
    gene: str,
    crispr: dict,
    rnai: dict | None,
) -> list[str]:
    ids = []
    if crispr.get("n_models", 0) > 0:
        ids.append(f"EV-DEPMAP-CRISPR-{gene}")
    if rnai and rnai.get("n_models", 0) > 0:
        ids.append(f"EV-DEPMAP-RNAI-{gene}")
    return ids


# ── 복수 유전자 평가 + dependency.tsv 출력 ────────────────

from .resource_usage import measured_tool


@measured_tool
def run_depmap_evaluation(
    genes: list[str],
    output_dir: str | Path,
    chronos_path: str | Path | None = None,
    demeter2_path: str | Path | None = None,
    model_meta_path: str | Path | None = None,
    depmap_release: str = "unknown",
) -> dict:
    """복수 유전자에 대해 DepMap 의존성을 평가하고 dependency.tsv를 저장.

    Args:
        genes: 평가할 유전자 목록.
        output_dir: 출력 디렉토리.
        chronos_path: CRISPR Chronos gene effect CSV.
        demeter2_path: RNAi DEMETER2 gene dependency CSV.
        model_meta_path: DepMap model metadata CSV.
        depmap_release: DepMap 릴리스 버전 (기록용).

    Returns:
        {
          "status": str,
          "depmap_release": str,
          "n_genes_evaluated": int,
          "results": {gene: evaluation_dict},
          "dependency_tsv": str | None,
          "summary": {status: count},
          "missing_evidence": list[str],
        }
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 파일 가용성 확인
    has_chronos  = chronos_path and Path(chronos_path).is_file()
    has_demeter2 = demeter2_path and Path(demeter2_path).is_file()
    has_meta     = model_meta_path and Path(model_meta_path).is_file()

    if not has_chronos:
        result = {
            "status": "UNAVAILABLE",
            "reason": f"CRISPR Chronos 파일 없음: {chronos_path}",
            "depmap_release": depmap_release,
            "n_genes_evaluated": 0,
            "results": {g: {"gene": g, "dependency_status": "UNAVAILABLE",
                            "reason": "Chronos 파일 없음"} for g in genes},
            "dependency_tsv": None,
            "summary": {"UNAVAILABLE": len(genes)},
            "missing_evidence": [f"DepMap_CRISPR_{g}" for g in genes],
        }
        _save_json(result, output_dir / "depmap_status.json")
        return result

    # 데이터 로드
    model_meta = load_model_metadata(model_meta_path) if has_meta else {}
    aliases = {}
    for mid, metadata in model_meta.items():
        alias = metadata.get("CCLEName") or metadata.get("CCLE_Name")
        if alias:
            if alias in aliases and aliases[alias] != mid:
                raise ValueError("Ambiguous CCLE name: " + alias)
            aliases[alias] = mid
    chronos_data = load_chronos(chronos_path, genes)
    demeter2_data = load_demeter2(demeter2_path, genes, aliases) if has_demeter2 else None
    breast_models = filter_breast_models(model_meta) if model_meta else {}

    # 모델 없으면 경고
    missing_evidence = []
    if not breast_models:
        result = {"status": "UNAVAILABLE", "reason": "Breast lineage metadata missing or no breast models",
                  "results": {}, "missing_evidence": ["breast_model_metadata"], "dependency_tsv": None}
        _save_json(result, output_dir / "depmap_status.json")
        return result

    # 유전자별 평가
    eval_results: dict[str, dict] = {}
    status_counts: dict[str, int] = {}

    for gene in genes:
        ev = evaluate_dependency(
            gene=gene,
            chronos_data=chronos_data,
            breast_models=breast_models,
            demeter2_data=demeter2_data,
        )
        eval_results[gene] = ev
        ds = ev["dependency_status"]
        status_counts[ds] = status_counts.get(ds, 0) + 1

        if not ev["applicable"]:
            missing_evidence.append(f"DepMap_{gene}")

    import csv
    with (output_dir / "depmap_measurements.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "model_id", "technology", "value", "release", "source"])
        for technology, measurements, source in (("CRISPR_Chronos", chronos_data, chronos_path),
                                                 ("RNAi_DEMETER2", demeter2_data or {}, demeter2_path)):
            for gene in genes:
                for model, value in measurements.get(gene, {}).items():
                    if model in breast_models:
                        writer.writerow([gene, model, technology, value, depmap_release, source])

    # dependency.tsv 저장
    tsv_path = _write_dependency_tsv(
        genes=genes,
        eval_results=eval_results,
        breast_models=breast_models,
        output_path=output_dir / "dependency.tsv",
        depmap_release=depmap_release,
    )

    result = {
        "status": "COMPLETED",
        "depmap_release": depmap_release,
        "n_genes_evaluated": len(genes),
        "n_breast_models": len(breast_models),
        "n_tnbc_models": sum(1 for s in breast_models.values() if s == "TNBC"),
        "results": eval_results,
        "dependency_tsv": tsv_path,
        "summary": status_counts,
        "missing_evidence": missing_evidence,
    }
    _save_json(result, output_dir / "depmap_status.json")
    return result


def _write_dependency_tsv(
    genes: list[str],
    eval_results: dict[str, dict],
    breast_models: dict[str, str],
    output_path: Path,
    depmap_release: str,
) -> str:
    """dependency.tsv 작성 (워크플로우 §4 출력 형식)."""
    header = [
        "gene", "dependency_status", "technology",
        "n_tnbc_models", "median_effect", "dependency_fraction",
        "common_essential_note", "concordance_note",
        "depmap_release", "applicability",
    ]
    rows = []
    for gene in genes:
        ev = eval_results.get(gene, {})
        crispr = ev.get("crispr_summary", {})
        rnai   = ev.get("rnai_summary")

        # CRISPR 행
        rows.append({
            "gene": gene,
            "dependency_status": ev.get("dependency_status", "UNAVAILABLE"),
            "technology": "CRISPR_Chronos",
            "n_tnbc_models": crispr.get("n_tnbc_models", 0),
            "median_effect": crispr.get("median_effect", ""),
            "dependency_fraction": crispr.get("dependency_fraction", ""),
            "common_essential_note": ev.get("common_essential_note", ""),
            "concordance_note": ev.get("concordance_note", ""),
            "depmap_release": depmap_release,
            "applicability": "APPLICABLE" if crispr.get("n_models", 0) > 0 else "UNKNOWN",
        })
        # RNAi 행 (있을 때만)
        if rnai and rnai.get("n_models", 0) > 0:
            rows.append({
                "gene": gene,
                "dependency_status": ev.get("dependency_status", "UNAVAILABLE"),
                "technology": "RNAi_DEMETER2",
                "n_tnbc_models": rnai.get("n_tnbc_models", 0),
                "median_effect": rnai.get("median_effect", ""),
                "dependency_fraction": rnai.get("dependency_fraction", ""),
                "common_essential_note": "",
                "concordance_note": "(보강 근거)",
                "depmap_release": depmap_release,
                "applicability": "SUPPORTIVE",
            })

    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(r.get(k, "")) for k in header))
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(output_path)


def _save_json(data: dict, path: Path) -> None:
    try:
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError:
        pass
