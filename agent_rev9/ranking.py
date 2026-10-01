from __future__ import annotations
"""
OmiCraft rev_1.0 — Pairwise Ranking (Choix Bradley-Terry)
Robin 패턴 적용: LLM이 후보 쌍을 비교 → 통계적 랭킹 (단순 점수합 대체)

흐름:
1. ADVANCE 표적들로 pairwise 조합 생성
2. LLM이 각 쌍을 비교하여 winner 선택
3. Choix (Bradley-Terry) 모델로 통계적 strength score 산출
4. 최종 순위 결정
"""

import json
import re
import itertools
import random
import logging
from typing import Optional

logger = logging.getLogger(__name__)


def generate_pairs(n_candidates: int, pairs_per_candidate: int = 3) -> list[tuple[int, int]]:
    """후보 수에 맞는 비교 쌍 생성.
    각 후보가 최소 pairs_per_candidate번 비교에 참여하도록 한다.
    """
    if n_candidates < 2:
        return []

    all_pairs = list(itertools.combinations(range(n_candidates), 2))
    random.shuffle(all_pairs)

    # 각 후보가 최소 N번 참여하도록 필터
    count = {i: 0 for i in range(n_candidates)}
    selected = []
    for a, b in all_pairs:
        if count[a] < pairs_per_candidate or count[b] < pairs_per_candidate:
            selected.append((a, b))
            count[a] += 1
            count[b] += 1

    # 최소한 모든 후보가 1번은 참여
    for i in range(n_candidates):
        if count[i] == 0 and selected:
            partner = (i + 1) % n_candidates
            selected.append((min(i, partner), max(i, partner)))

    return selected


def run_pairwise_comparison(targets: list[dict], llm, pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """LLM으로 pairwise 비교 실행.

    Args:
        targets: ADVANCE 표적 리스트
        llm: LangChain ChatOpenAI instance
        pairs: 비교할 쌍 인덱스

    Returns: [(winner_idx, loser_idx), ...] 게임 결과
    """
    from langchain_core.messages import HumanMessage
    from .prompts import RANKING_SYSTEM, RANKING_COMPARE

    games = []

    for a_idx, b_idx in pairs:
        ta = targets[a_idx]
        tb = targets[b_idx]

        prompt = RANKING_COMPARE.format(
            gene_a=ta.get("gene_name", ""),
            tier_a=f"{ta.get('bio_tier', '')}+{ta.get('dev_tier', '')}",
            modality_a=ta.get("modality", ""),
            lfc_a=f"{ta.get('log2fc', 0):.2f}",
            safety_a=ta.get("safety_verdict", "PASS"),
            phase_a=ta.get("chembl_max_phase", 0),
            gene_b=tb.get("gene_name", ""),
            tier_b=f"{tb.get('bio_tier', '')}+{tb.get('dev_tier', '')}",
            modality_b=tb.get("modality", ""),
            lfc_b=f"{tb.get('log2fc', 0):.2f}",
            safety_b=tb.get("safety_verdict", "PASS"),
            phase_b=tb.get("chembl_max_phase", 0),
        )

        try:
            response = llm.invoke([
                HumanMessage(content=RANKING_SYSTEM),
                HumanMessage(content=prompt),
            ])

            # JSON 파싱
            cleaned = re.sub(r'```json\s*', '', response.content)
            cleaned = re.sub(r'```\s*', '', cleaned).strip()
            match = re.search(r'\{[\s\S]*?\}', cleaned)

            if match:
                result = json.loads(match.group())
                winner = result.get("winner", "A")
                if winner == "A":
                    games.append((a_idx, b_idx))
                else:
                    games.append((b_idx, a_idx))

                logger.info(f"  {ta.get('gene_name','')} vs {tb.get('gene_name','')}"
                           f" → winner: {'A' if winner=='A' else 'B'}"
                           f" ({targets[games[-1][0]].get('gene_name','')})")
            else:
                # 파싱 실패 → 스킵
                logger.warning(f"  Ranking parse failed for {ta.get('gene_name','')} vs {tb.get('gene_name','')}")

        except Exception as e:
            logger.warning(f"  Ranking LLM error: {e}")

    return games


def compute_ranking_scores(n_candidates: int, games: list[tuple[int, int]]) -> Optional[list[float]]:
    """Choix Bradley-Terry 모델로 strength score 산출.

    Returns: 각 후보의 strength score 리스트 (높을수록 강함), 실패 시 None
    """
    if not games or n_candidates < 2:
        return None

    try:
        import choix
        params = choix.ilsr_pairwise(n_candidates, games, alpha=0.1)
        return params.tolist()
    except ImportError:
        logger.warning("choix 패키지 미설치 — 단순 승률로 대체")
        # fallback: 단순 승률
        wins = {i: 0 for i in range(n_candidates)}
        total = {i: 0 for i in range(n_candidates)}
        for w, l in games:
            wins[w] += 1
            total[w] += 1
            total[l] += 1
        return [wins[i] / max(total[i], 1) for i in range(n_candidates)]
    except Exception as e:
        logger.error(f"Choix ranking failed: {e}")
        return None


def rank_advance_targets(targets: list[dict], llm, pairs_per_candidate: int = 3) -> list[dict]:
    """ADVANCE 표적을 pairwise ranking으로 순위화.

    Returns: strength_score가 추가된 표적 리스트 (내림차순 정렬)
    """
    n = len(targets)
    if n < 2:
        for t in targets:
            t["strength_score"] = 1.0
            t["rank"] = 1
        return targets

    logger.info(f"[Ranking] {n}개 표적 pairwise ranking 시작")

    pairs = generate_pairs(n, pairs_per_candidate)
    logger.info(f"[Ranking] {len(pairs)}개 비교 쌍 생성")

    games = run_pairwise_comparison(targets, llm, pairs)
    logger.info(f"[Ranking] {len(games)}개 게임 완료")

    scores = compute_ranking_scores(n, games)

    if scores:
        for i, t in enumerate(targets):
            t["strength_score"] = round(scores[i], 4)

        targets.sort(key=lambda x: x.get("strength_score", 0), reverse=True)

        for i, t in enumerate(targets):
            t["rank"] = i + 1

        logger.info("[Ranking] 최종 순위:")
        for t in targets:
            logger.info(f"  #{t['rank']} {t.get('gene_name',''):12s} "
                       f"score={t['strength_score']:+.3f} → {t.get('modality','')}")
    else:
        logger.warning("[Ranking] 점수 산출 실패 — 기존 순서 유지")
        for i, t in enumerate(targets):
            t["strength_score"] = 0.0
            t["rank"] = i + 1

    return targets
