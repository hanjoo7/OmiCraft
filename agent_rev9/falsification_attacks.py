"""
FalsiTarget — 6가지 결정론적 공격 도구
논문: "Attack Your Own Findings: A Falsification Agent Framework"

모든 공격은 결정론적 함수 — LLM을 사용하지 않으므로 재현 가능하다.
각 공격은 (passed: bool, evidence: str) 튜플을 반환한다.

공격 목록:
  1. Scale invariance: 정규화 방식을 바꿔도 순위가 유지되는가
  2. Null calibration: 각 근거 구성요소가 chance 이상인가
  3. Negative control: 알려진 non-target이 후보에 포함되는가
  4. Confounder: 후보 점수가 교란 변수와 과도하게 상관하는가
  5. Robustness grid: 집계 방식을 바꿔도 결과가 안정한가
  6. Provenance: 메타데이터 변수의 의미가 문서와 일치하는가
"""

import os
import sys
import numpy as np
import pandas as pd
from typing import Tuple
from scipy import stats

PROCESSED_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed"
)

# ── 알려진 non-target (negative control) ──────────────────
# 유방암/TNBC에서 표적이 될 수 없는 유전자
NEGATIVE_CONTROLS = {
    # Housekeeping — 모든 세포에서 필수, 표적 불가
    "GAPDH": "glycolysis housekeeping",
    "ACTB": "cytoskeleton housekeeping",
    "RPL11": "ribosomal subunit",
    "EEF1A1": "translation factor",
    "TUBB": "tubulin essential",
    # 혈액계 — 고형종양 표적 부적합
    "HBB": "hemoglobin — erythrocyte specific",
    "CD3E": "T-cell marker — immune, not tumor",
    "MS4A1": "CD20 — B-cell marker",
    # 조직 부적합
    "INS": "insulin — pancreatic beta cell",
    "ALB": "albumin — liver specific",
}


def attack_scale_invariance(
    advance_genes: list,
    deg_path: str = None,
    rank_corr_threshold: float = 0.8,
) -> Tuple[bool, str]:
    """공격 1: Scale Invariance
    DEG 결과를 다른 정규화(log2FC vs stat)로 재순위화했을 때
    ADVANCE 표적의 순위가 유지되는지 확인.

    통과: 순위 상관 >= threshold
    실패: 정규화 변경 시 순위가 크게 바뀜 → 결과가 정규화에 의존
    """
    if deg_path is None:
        deg_path = os.path.join(PROCESSED_DIR, "deg_results", "deg_tnbc_vs_nontnbc_full.csv")

    if not os.path.exists(deg_path):
        return False, "DEG 결과 파일 없음 — 검증 불가"

    df = pd.read_csv(deg_path, sep="\t" if str(deg_path).endswith(".tsv") else ",")
    gene_column = next((key for key in ("gene_symbol", "gene_name", "gene") if key in df), df.columns[0])
    df = df.set_index(gene_column, drop=False)
    df = df.dropna(subset=["log2FoldChange", "padj"])

    # 방법 A: |log2FC| * -log10(padj) (기존 raw score)
    df["score_raw"] = df["log2FoldChange"].abs() * (-np.log10(df["padj"].clip(1e-300)))
    df["rank_raw"] = df["score_raw"].rank(ascending=False)

    # 방법 B: percentile rank 기반 (편향 제거 점수)
    n = len(df)
    df["score_pct"] = (df["log2FoldChange"].abs().rank() / n) * ((-np.log10(df["padj"].clip(1e-300))).rank() / n)
    df["rank_pct"] = df["score_pct"].rank(ascending=False)

    # ADVANCE 표적의 순위 비교
    advance_set = set(advance_genes)
    if "gene_name" in df.columns:
        mask = df["gene_name"].isin(advance_set)
    else:
        mask = df.index.isin(advance_set)

    if mask.sum() < 3:
        return True, f"ADVANCE 표적 {mask.sum()}개만 매칭 — 검증 생략"

    sub = df[mask]
    corr, pval = stats.spearmanr(sub["rank_raw"], sub["rank_pct"])

    passed = corr >= rank_corr_threshold
    evidence = (f"순위 상관 ρ={corr:.3f} (p={pval:.3f}), threshold={rank_corr_threshold}. "
                f"{'PASS: 정규화에 강건' if passed else 'FAIL: 정규화 변경 시 순위 불안정'}")
    return passed, evidence


def attack_null_calibration(
    advance_targets: list,
    qualification_path: str = None,
) -> Tuple[bool, str]:
    """공격 2: Null Calibration
    각 ADVANCE 표적의 근거 구성요소(bio_score, dev_score)가
    chance level 이상인지 확인.

    통과: 모든 ADVANCE 표적의 모든 구성요소 > 0
    실패: 어떤 구성요소가 0 또는 음수 → anti-predictive 구성요소가 희석되어 숨어 있음
    """
    if qualification_path is None:
        qualification_path = os.path.join(PROCESSED_DIR, "qualification_results", "qualification_full.csv")

    if not os.path.exists(qualification_path):
        return False, "Qualification 결과 없음"

    df = pd.read_csv(qualification_path)
    advance_df = df[df["judgment"] == "ADVANCE"]

    if len(advance_df) == 0:
        return False, "ADVANCE 표적 0건"

    # 각 구성요소 체크
    anti_predictive = []
    for _, row in advance_df.iterrows():
        gene = row.get("gene_name", "?")
        bio = row.get("bio_score", 0)
        dev = row.get("dev_score", 0)

        if bio <= 0:
            anti_predictive.append(f"{gene}: bio_score={bio} (chance 이하)")
        if dev <= 0:
            anti_predictive.append(f"{gene}: dev_score={dev} (chance 이하)")

    passed = len(anti_predictive) == 0
    if passed:
        evidence = f"모든 ADVANCE ({len(advance_df)}건) 구성요소 > 0 — 희석된 anti-predictive 없음"
    else:
        evidence = f"Anti-predictive 구성요소 {len(anti_predictive)}건: {'; '.join(anti_predictive[:3])}"
    return passed, evidence


def attack_negative_control(
    advance_genes: list,
    negative_controls: dict = None,
) -> Tuple[bool, str]:
    """공격 3: Negative Control
    알려진 non-target 유전자가 ADVANCE 목록에 포함되어 있으면 FAIL.

    논문: "synthesised tissue-inappropriate antigens reach the real candidates"
    → 가짜 표적이 진짜 후보에 섞이면 파이프라인 자체가 편향
    """
    if negative_controls is None:
        negative_controls = NEGATIVE_CONTROLS

    advance_set = set(advance_genes)
    contaminated = {g: reason for g, reason in negative_controls.items() if g in advance_set}

    passed = len(contaminated) == 0
    if passed:
        evidence = f"Negative control {len(negative_controls)}개 중 ADVANCE 침투 0건"
    else:
        details = "; ".join(f"{g} ({r})" for g, r in contaminated.items())
        evidence = f"FAIL: Non-target {len(contaminated)}개가 ADVANCE에 침투 — {details}"
    return passed, evidence


def attack_confounder(
    advance_genes: list,
    deg_path: str = None,
    corr_threshold: float = 0.3,
) -> Tuple[bool, str]:
    """공격 4: Confounder
    후보 유전자의 composite score가 평균 발현량과 과도하게 상관하면 FAIL.

    논문: "fails if a component correlates with expression level above |ρ|=0.3 across all genes"
    → 발현량이 높으면 자동으로 높은 점수 → 발현량 편향
    """
    pool_path = os.path.join(PROCESSED_DIR, "candidate_pool", "omics_candidate_pool.csv")
    if not os.path.exists(pool_path):
        return True, "Candidate Pool 없음 — 검증 생략"

    df = pd.read_csv(pool_path, index_col=0)
    if "composite_score" not in df.columns or "baseMean" not in df.columns:
        return True, "composite_score/baseMean 컬럼 없음 — 검증 생략"

    df = df.dropna(subset=["composite_score", "baseMean"])
    if len(df) < 10:
        return True, "데이터 부족 — 검증 생략"

    corr, pval = stats.spearmanr(df["composite_score"], df["baseMean"])

    passed = abs(corr) <= corr_threshold
    evidence = (f"Score-expression 상관 ρ={corr:.3f} (p={pval:.1e}), "
                f"threshold=|{corr_threshold}|. "
                f"{'PASS: 발현량 편향 없음' if passed else 'FAIL: score가 발현량에 편향'}")
    return passed, evidence


def attack_robustness_grid(
    advance_genes: list,
    deg_path: str = None,
    spread_threshold: float = 0.2,
) -> Tuple[bool, str]:
    """공격 5: Robustness Grid
    DEG 분석에서 다른 통계 방법(log2FC vs stat vs -log10(padj))으로
    ADVANCE 표적의 정규화된 순위가 안정한지 확인.

    통과: 방법 간 순위 spread < threshold
    실패: 방법에 따라 결과가 크게 달라짐
    """
    if deg_path is None:
        deg_path = os.path.join(PROCESSED_DIR, "deg_results", "deg_tnbc_vs_nontnbc_full.csv")

    if not os.path.exists(deg_path):
        return True, "DEG 결과 없음 — 검증 생략"

    df = pd.read_csv(deg_path, sep="\t" if str(deg_path).endswith(".tsv") else ",")
    gene_column = next((key for key in ("gene_symbol", "gene_name", "gene") if key in df), df.columns[0])
    df = df.set_index(gene_column, drop=False)
    df = df.dropna(subset=["log2FoldChange", "padj"])

    n_total = len(df)
    if n_total == 0:
        return True, "DEG 0건"

    # 3가지 방법으로 정규화 순위 (0~1)
    df["pct_lfc"] = df["log2FoldChange"].abs().rank() / n_total
    df["pct_nlogp"] = (-np.log10(df["padj"].clip(1e-300))).rank() / n_total
    if "stat" in df.columns:
        df["pct_stat"] = df["stat"].abs().rank() / n_total
    else:
        df["pct_stat"] = df["pct_nlogp"]

    # ADVANCE 표적의 순위 spread
    if "gene_name" in df.columns:
        mask = df["gene_name"].isin(set(advance_genes))
    else:
        mask = df.index.isin(set(advance_genes))

    sub = df[mask]
    if len(sub) < 2:
        return True, f"ADVANCE 매칭 {len(sub)}건 — 검증 생략"

    spreads = []
    for _, row in sub.iterrows():
        vals = [row["pct_lfc"], row["pct_nlogp"], row["pct_stat"]]
        spreads.append(max(vals) - min(vals))

    max_spread = max(spreads)
    mean_spread = np.mean(spreads)

    passed = max_spread <= spread_threshold
    evidence = (f"순위 spread: mean={mean_spread:.3f}, max={max_spread:.3f}, "
                f"threshold={spread_threshold}. "
                f"{'PASS: 방법 간 안정' if passed else 'FAIL: 방법에 따라 순위 불안정'}")
    return passed, evidence


def attack_provenance(
    state_contrasts: list,
    state_cohort: str = "",
) -> Tuple[bool, str]:
    """공격 6: Provenance
    메타데이터 변수의 의미가 실제 데이터 문서와 일치하는지 확인.

    논문 E4: "cohort 변수를 primary-vs-metastatic으로 해석했지만
    실제로는 두 개의 다른 임상시험이었음"

    여기서는 contrast 설정이 데이터와 일치하는지 확인:
    - TNBC 아형 분류 기준이 ER-/PR-/HER2-와 일치하는가
    - 비교군이 실제 존재하는 샘플 그룹인가
    """
    meta_path = os.path.join(PROCESSED_DIR, "sample_metadata.csv")
    if not os.path.exists(meta_path):
        return False, "sample_metadata.csv 없음 — provenance 검증 불가"

    meta = pd.read_csv(meta_path)

    issues = []

    # P1: TNBC 정의 확인 (ER-/PR-/HER2-)
    if "is_tnbc" in meta.columns and "er_status" in meta.columns:
        tnbc_samples = meta[meta["is_tnbc"] == 1]
        if len(tnbc_samples) > 0:
            er_pos = (tnbc_samples["er_status"].str.lower() == "positive").sum()
            pr_pos = (tnbc_samples["pr_status"].str.lower() == "positive").sum()
            her2_pos = (tnbc_samples["her2_status"].str.lower() == "positive").sum()

            if er_pos > 0:
                issues.append(f"TNBC 중 ER+ {er_pos}건 — 아형 분류 오류 가능")
            if pr_pos > 0:
                issues.append(f"TNBC 중 PR+ {pr_pos}건 — 아형 분류 오류 가능")
            if her2_pos > 0:
                issues.append(f"TNBC 중 HER2+ {her2_pos}건 — 아형 분류 오류 가능")

    # P2: 비교군 존재 확인
    subtypes = set(meta.get("subtype", pd.Series()).unique())
    if "TNBC" not in subtypes:
        issues.append("TNBC 아형이 metadata에 없음")

    # P3: contrast에 언급된 그룹이 실제 존재하는지
    for contrast in state_contrasts:
        if "normal" in contrast.lower() and "Normal" not in subtypes and "normal" not in str(meta.columns):
            issues.append(f"Contrast '{contrast}' — 정상조직 샘플이 metadata에 없을 수 있음 (확인 필요)")

    passed = len(issues) == 0
    if passed:
        evidence = "Provenance 검증 통과: TNBC 정의(ER-/PR-/HER2-) 일치, 비교군 존재 확인"
    else:
        evidence = f"Provenance 이슈 {len(issues)}건: {'; '.join(issues)}"
    return passed, evidence


def _current_run_attacks(advance_targets, state):
    """Audit current artifacts; unavailable tests are not biological contradictions."""
    from pathlib import Path
    result = state['r_analysis_result']
    genes = [t['gene_name'] for t in advance_targets if t.get('gene_name')]
    attacks = {}
    def checked(name, fn, *args):
        try:
            passed, evidence = fn(*args)
            skipped = '생략' in evidence or '검증 불가' in evidence
            attacks[name] = {'passed': None if skipped else bool(passed),
                             'status': 'NOT_EVALUATED' if skipped else 'PASS' if passed else 'FAIL',
                             'evidence': evidence}
        except (OSError, ValueError, KeyError, TypeError) as exc:
            attacks[name] = {'passed': None, 'status': 'NOT_EVALUATED', 'evidence': str(exc)}
    de_path = (result.get('output_files') or {}).get('de_all')
    for name, fn in [('scale_invariance', attack_scale_invariance), ('robustness_grid', attack_robustness_grid)]:
        if de_path:
            checked(name, fn, genes, de_path)
        else:
            attacks[name] = {'passed': None, 'status': 'NOT_EVALUATED', 'evidence': '현재 실행의 DE 결과 미확보'}
    checked('negative_control', attack_negative_control, genes)
    # The current qualification contract has categorical B/M/D grades, not numeric
    # bio_score/dev_score or a calibrated null distribution. Do not read a legacy CSV.
    attacks['null_calibration'] = {'passed': None, 'status': 'NOT_EVALUATED',
        'evidence': '현재 B/M/D 등급에 대한 null calibration 자료 미확보; 독립 검증 필요'}
    attacks['confounder'] = {'passed': None, 'status': 'NOT_EVALUATED',
        'evidence': '현재 후보의 composite score–발현량 교란 검증 자료 미확보; DE 공변량 보정과 별도 검증 필요'}
    def provenance():
        manifest = result.get('input_manifest') or {}
        counts = manifest.get('group_counts', {})
        metadata = Path(result.get('output_dir') or '.') / 'sample_metadata.tsv'
        if not counts or not metadata.is_file():
            return False, '현재 실행의 샘플 메타데이터 미확보 — 검증 불가'
        meta = pd.read_csv(metadata, sep='\t')
        required = {'group', 'er_status', 'pr_status', 'her2_status'}
        if not required.issubset(meta.columns):
            return False, '수용체 상태/그룹 열 미확보 — 검증 불가'
        observed = {str(k): int(v) for k, v in meta['group'].value_counts().items()}
        if observed != counts:
            return False, '현재 샘플 메타데이터와 input_manifest 그룹 수 불일치'
        if not observed.get('TNBC') or not observed.get('Non_TNBC'):
            return False, '현재 분석에 필요한 TNBC/Non_TNBC 비교군 누락'
        tnbc = meta.loc[meta['group'] == 'TNBC', ['er_status', 'pr_status', 'her2_status']]
        normalized = tnbc.apply(lambda col: col.astype(str).str.strip().str.lower())
        if normalized.isin(['positive', 'pos', '+']).any().any():
            return False, 'TNBC 그룹에 수용체 양성 샘플 존재'
        if not normalized.isin(['negative', 'neg', '-']).all().all():
            return False, 'TNBC 그룹의 수용체 음성 상태 일부 미확인 — 검증 불가'
        for contrast in state.get('contrasts', []):
            if isinstance(contrast, dict) and contrast.get('status') == 'READY':
                for key in ('positive_group', 'reference_group'):
                    group = str(contrast.get(key, '')).replace('Non-TNBC', 'Non_TNBC')
                    if not observed.get(group):
                        return False, 'READY contrast의 비교군이 현재 메타데이터에 없음: ' + group
        return True, f"현재 메타데이터 검증: TNBC={observed['TNBC']}, Non_TNBC={observed['Non_TNBC']}; TNBC 수용체 음성 확인"
    checked('provenance', provenance)
    failed = [name for name, row in attacks.items() if row['status'] == 'FAIL']
    unknown = [name for name, row in attacks.items() if row['status'] == 'NOT_EVALUATED']
    passed = sum(row['status'] == 'PASS' for row in attacks.values())
    return {'attacks': attacks, 'n_passed': passed, 'n_failed': len(failed),
            'n_not_evaluated': len(unknown), 'not_evaluated': unknown, 'failed': failed,
            'n_total': len(attacks), 'all_passed': passed == len(attacks),
            'single_points_of_failure': failed if len(failed) == 1 else [],
            'verdict': 'REVISE' if failed else 'HOLD' if unknown else 'ACCEPT',
            'coverage': f'{passed}/{len(attacks)} attacks passed'}


# ── 전체 공격 실행 ────────────────────────────────────────

def run_all_attacks(
    advance_targets: list,
    state: dict = None,
) -> dict:
    """6가지 공격 전부 실행.

    Args:
        advance_targets: ADVANCE 판정된 표적 리스트
        state: OmiCraftState (contrasts, cohort 등)

    Returns:
        {
            "attacks": {name: {"passed": bool, "evidence": str}},
            "n_passed": int,
            "n_failed": int,
            "all_passed": bool,
            "single_points_of_failure": list,
            "verdict": "ACCEPT" | "REVISE" | "REJECT",
        }
    """
    if state is None:
        state = {}

    if state.get('r_analysis_result'):
        return _current_run_attacks(advance_targets, state)

    advance_genes = [t.get("gene_name", "") for t in advance_targets if t.get("gene_name")]

    attacks = {}

    # 1. Scale invariance
    p, e = attack_scale_invariance(advance_genes)
    attacks["scale_invariance"] = {"passed": p, "evidence": e}

    # 2. Null calibration
    p, e = attack_null_calibration(advance_targets)
    attacks["null_calibration"] = {"passed": p, "evidence": e}

    # 3. Negative control
    p, e = attack_negative_control(advance_genes)
    attacks["negative_control"] = {"passed": p, "evidence": e}

    # 4. Confounder
    p, e = attack_confounder(advance_genes)
    attacks["confounder"] = {"passed": p, "evidence": e}

    # 5. Robustness grid
    p, e = attack_robustness_grid(advance_genes)
    attacks["robustness_grid"] = {"passed": p, "evidence": e}

    # 6. Provenance
    p, e = attack_provenance(state.get("contrasts", []), state.get("cohort", ""))
    attacks["provenance"] = {"passed": p, "evidence": e}

    # 집계
    n_passed = sum(1 for a in attacks.values() if a["passed"])
    n_failed = sum(1 for a in attacks.values() if not a["passed"])
    failed_names = [name for name, a in attacks.items() if not a["passed"]]

    # Arbiter 규칙 (논문 Section 3.4)
    # - confounder 실패 → REJECT (즉시)
    # - negative_control, null_calibration, provenance 실패 → REVISE (회귀 필요)
    # - scale_invariance, robustness_grid 실패 → WARNING (flagged but not blocking)
    blocking_failures = [n for n in failed_names if n not in ("scale_invariance", "robustness_grid")]
    warning_failures = [n for n in failed_names if n in ("scale_invariance", "robustness_grid")]

    if "confounder" in failed_names:
        verdict = "REJECT"
    elif len(blocking_failures) > 0:
        verdict = "REVISE"
    else:
        verdict = "ACCEPT"  # warning만 있으면 ACCEPT (flagged)

    # 단일 실패점 분석: 이 공격 하나만 빼면 에러가 통과하는가
    single_points = []
    for name in failed_names:
        other_failures = [n for n in failed_names if n != name]
        if len(other_failures) == 0:
            single_points.append(name)

    return {
        "attacks": attacks,
        "n_passed": n_passed,
        "n_failed": n_failed,
        "n_total": len(attacks),
        "all_passed": n_failed == 0,
        "failed": failed_names,
        "single_points_of_failure": single_points,
        "verdict": verdict,
        "coverage": f"{n_passed}/{len(attacks)} attacks passed",
    }


if __name__ == "__main__":
    import sys
    if sys.stdout.encoding != "utf-8":
        sys.stdout.reconfigure(encoding="utf-8")

    print("=" * 60)
    print("FalsiTarget — 6공격 독립 테스트")
    print("=" * 60)

    # 더미 ADVANCE 표적
    test_targets = [
        {"gene_name": "ESR1", "bio_score": 3, "dev_score": 4},
        {"gene_name": "AGTR1", "bio_score": 3, "dev_score": 5},
        {"gene_name": "GATA3", "bio_score": 5, "dev_score": 2},
        {"gene_name": "GAPDH", "bio_score": 2, "dev_score": 3},  # negative control!
    ]

    result = run_all_attacks(test_targets, state={"contrasts": ["TNBC vs normal", "TNBC vs non-TNBC"]})

    for name, attack in result["attacks"].items():
        icon = "✓" if attack["passed"] else "✗"
        print(f"  {icon} {name:25s}: {attack['evidence'][:80]}")

    print(f"\n  Verdict: {result['verdict']}")
    print(f"  Coverage: {result['coverage']}")
    if result["single_points_of_failure"]:
        print(f"  Single points of failure: {result['single_points_of_failure']}")
    if result["failed"]:
        print(f"  Failed attacks: {result['failed']}")
