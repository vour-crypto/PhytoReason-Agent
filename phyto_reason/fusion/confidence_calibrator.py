"""
confidence_calibrator.py — 置信度校准器。

对融合后的原始分数进行校准:
  - 矛盾检测惩罚
  - 收敛度奖励
  - 多样性最低要求
  - Gold/Silver/Weak 分级
"""

from __future__ import annotations

from phyto_reason.models.enums import ConfidenceLevel
from phyto_reason.fusion.fusion_result import (
    ContradictionFlag, ContradictionLevel, EvidenceHierarchy, FusionEvidence,
)
from phyto_reason.fusion.scoring_policy import ScoringPolicy, DEFAULT_POLICY


class ConfidenceCalibrator:
    """置信度校准器。

    策略:
      1. 从 raw_fused_score 开始
      2. 应用矛盾惩罚 (多次矛盾累加, capped)
      3. 收敛度奖励 (高收敛 + 加分)
      4. 多样性惩罚 (证据源不足 → 降级)
      5. 确定 Gold/Silver/Weak 分级
    """

    def __init__(self, policy: ScoringPolicy | None = None) -> None:
        self.policy = policy or DEFAULT_POLICY

    def calibrate(
        self,
        raw_fused_score: float,
        evidence_list: list[FusionEvidence],
        contradictions: list[ContradictionFlag],
        hierarchy: EvidenceHierarchy,
    ) -> tuple[float, ConfidenceLevel, str]:
        """执行完整校准。

        Returns:
            (calibrated_score, confidence_level, explanation)
        """
        score = raw_fused_score

        # 1. 基础检查: 证据不足
        if not evidence_list:
            return (0.0, ConfidenceLevel.NONE, "无可用证据。")

        # 2. 矛盾惩罚
        total_penalty = 0.0
        penalty_parts: list[str] = []
        for c in contradictions:
            penalty = c.penalty
            total_penalty += penalty
            penalty_parts.append(
                f"  - [{c.level.value}] {c.pattern}: -{penalty:.2f}"
            )

        # 对 penaly 设上限: 不超过 40% 原始分
        total_penalty = min(total_penalty, raw_fused_score * 0.4)
        score -= total_penalty

        # 3. 多样性检查
        n_active = hierarchy.n_independent_sources
        if n_active < self.policy.min_evidence_sources:
            diversity_ratio = n_active / self.policy.min_evidence_sources
            score *= diversity_ratio

        # 4. 收敛奖励
        if hierarchy.convergence_score >= 0.70 and n_active >= 2:
            bonus = self.policy.convergence_bonus_max * hierarchy.convergence_score
            score += bonus

        # 5. 确保非负且 ≤ 1.0
        score = max(0.0, min(score, 1.0))

        # 6. 分级
        level = ConfidenceLevel.NONE
        if score >= self.policy.gold_threshold and n_active >= self.policy.min_evidence_sources:
            level = ConfidenceLevel.GOLD
        elif score >= self.policy.silver_threshold:
            level = ConfidenceLevel.SILVER
        elif score >= 0.10:
            level = ConfidenceLevel.WEAK

        # 构建解释
        parts = [f"原始融合分: {raw_fused_score:.3f}"]
        if penalty_parts:
            parts.append("矛盾惩罚:")
            parts.extend(penalty_parts)
        if hierarchy.convergence_score >= 0.70:
            parts.append(f"收敛奖励: +{self.policy.convergence_bonus_max * hierarchy.convergence_score:.3f}")
        if n_active < self.policy.min_evidence_sources:
            parts.append(f"多样性不足 ({n_active}/{self.policy.min_evidence_sources}), 分数调整 ×{diversity_ratio:.2f}")
        parts.append(f"校准后分数: {score:.3f}")
        parts.append(f"置信级别: {level.value}")

        return (score, level, "\n".join(parts))
