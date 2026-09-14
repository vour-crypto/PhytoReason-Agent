"""
evidence_ranker.py — 证据重要性排序器。

建立证据层次结构、计算收敛度和多样性。
"""

from __future__ import annotations

from phyto_reason.fusion.fusion_result import (
    EvidenceHierarchy, FusionEvidence,
)


# ── 证据重要性层次 ──────────────────────────────────────────
# 从生物学意义出发的固有重要性
# 直接调控证据 > 间接关联 > 先验知识 > 表达相关性

EVIDENCE_TIER: dict[str, int] = {
    "motif": 5,            # 直接结合证据
    "pathway": 4,          # 通路位置证据
    "tf_prior": 3,         # 家族先验
    "ortholog": 3,         # 进化保守性
    "tissue": 2,           # 组织共定位
    "correlation": 2,      # 表达相关
    "module_membership": 1,  # 模块归属
}


class EvidenceRanker:
    """证据排序器。

    对融合后的证据列表进行:
      1. 重要性分级
      2. 收敛度计算 (证据间一致性)
      3. 多样性计算 (独立证据源数量)
      4. 构建层次结构
    """

    @staticmethod
    def get_tier(source: str) -> int:
        return EVIDENCE_TIER.get(source, 1)

    @staticmethod
    def rank_by_contribution(
        evidence_list: list[FusionEvidence],
    ) -> list[FusionEvidence]:
        """按贡献度降序排列。"""
        return sorted(evidence_list, key=lambda e: e.contribution, reverse=True)

    @staticmethod
    def rank_by_tier(
        evidence_list: list[FusionEvidence],
    ) -> list[FusionEvidence]:
        """按证据层次降序排列 (同层内按贡献)。"""
        return sorted(
            evidence_list,
            key=lambda e: (EVIDENCE_TIER.get(e.source, 1), e.contribution),
            reverse=True,
        )

    @staticmethod
    def compute_convergence(
        evidence_list: list[FusionEvidence],
    ) -> float:
        """计算证据收敛度。

        多个独立证据指向相同结论 → 高收敛。
        使用 pairwise score consistency。
        """
        if len(evidence_list) < 2:
            return 0.0

        scores = [e.normalized_score for e in evidence_list]
        mean = sum(scores) / len(scores)

        # 方差越低 → 收敛度越高
        variance = sum((s - mean) ** 2 for s in scores) / len(scores)
        # 将方差映射到 0-1: 低方差=高收敛
        convergence = max(0.0, 1.0 - variance * 3)
        return round(convergence, 3)

    @staticmethod
    def compute_diversity(
        evidence_list: list[FusionEvidence],
    ) -> int:
        """计算独立证据源数量。"""
        sources = set()
        for e in evidence_list:
            if e.normalized_score >= 0.20:  # 仅计算有意义的证据
                sources.add(e.source)
        return len(sources)

    @staticmethod
    def build_hierarchy(
        evidence_list: list[FusionEvidence],
    ) -> EvidenceHierarchy:
        """构建完整证据层次。"""
        if not evidence_list:
            return EvidenceHierarchy()

        ranked = EvidenceRanker.rank_by_tier(evidence_list)
        top_sources = [e.source for e in ranked[:5]]

        convergence = EvidenceRanker.compute_convergence(evidence_list)
        diversity = EvidenceRanker.compute_diversity(evidence_list)

        return EvidenceHierarchy(
            top_sources=top_sources,
            convergence_score=convergence,
            diversity_score=min(diversity / 7.0, 1.0),
            n_independent_sources=diversity,
        )
