"""
scoring_policy.py — 统一权重配置系统。

定义每个证据维度的权重、校准参数和级别阈值。
支持 configurable policies。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class ScoringPolicy(BaseModel):
    """证据融合权重配置。

    各维度权重总和应 ≤ 1.0 (实际融合时会在有效维度间归一化)。
    """
    correlation: float = Field(default=0.15, ge=0.0, le=1.0)
    module_membership: float = Field(default=0.10, ge=0.0, le=1.0)
    motif: float = Field(default=0.08, ge=0.0, le=1.0)
    pathway: float = Field(default=0.15, ge=0.0, le=1.0)
    tf_prior: float = Field(default=0.15, ge=0.0, le=1.0)
    tissue: float = Field(default=0.10, ge=0.0, le=1.0)
    ortholog: float = Field(default=0.10, ge=0.0, le=1.0)
    literature: float = Field(default=0.10, ge=0.0, le=1.0)

    gold_threshold: float = Field(default=0.70, ge=0.0, le=1.0)
    silver_threshold: float = Field(default=0.40, ge=0.0, le=1.0)

    min_evidence_sources: int = Field(default=2, ge=1, le=7)

    contradiction_high_penalty: float = Field(default=0.30, ge=0.0, le=1.0)
    contradiction_medium_penalty: float = Field(default=0.15, ge=0.0, le=1.0)
    contradiction_low_penalty: float = Field(default=0.05, ge=0.0, le=1.0)

    convergence_bonus_max: float = Field(default=0.10, ge=0.0, le=1.0)
    diversity_minimum: int = Field(default=2, ge=1, le=7)

    def list_active_weights(self) -> dict[str, float]:
        """返回所有非零权重。"""
        return {
            k: v for k, v in self.model_dump().items()
            if k in ("correlation", "module_membership", "motif", "pathway",
                     "tf_prior", "tissue", "ortholog", "literature")
            and v > 0
        }

    def get_weight(self, dimension: str) -> float:
        return getattr(self, dimension, 0.0)


DEFAULT_POLICY = ScoringPolicy()


class PolicyLibrary:
    """预定义策略库，支持按需切换。"""

    @staticmethod
    def balanced() -> ScoringPolicy:
        return ScoringPolicy(
            correlation=0.15, module_membership=0.10,
            motif=0.08, pathway=0.15,
            tf_prior=0.15, tissue=0.10, ortholog=0.10, literature=0.15,
        )

    @staticmethod
    def motif_focused() -> ScoringPolicy:
        return ScoringPolicy(
            correlation=0.10, module_membership=0.05,
            motif=0.35, pathway=0.20,
            tf_prior=0.10, tissue=0.10, ortholog=0.10,
        )

    @staticmethod
    def prior_knowledge_focused() -> ScoringPolicy:
        return ScoringPolicy(
            correlation=0.10, module_membership=0.05,
            motif=0.10, pathway=0.20,
            tf_prior=0.30, tissue=0.10, ortholog=0.15,
        )

    @staticmethod
    def expression_only() -> ScoringPolicy:
        return ScoringPolicy(
            correlation=0.40, module_membership=0.30,
            motif=0.10, pathway=0.10,
            tf_prior=0.05, tissue=0.05, ortholog=0.0,
        )

    @staticmethod
    def strict() -> ScoringPolicy:
        return ScoringPolicy(
            correlation=0.10, module_membership=0.05,
            motif=0.25, pathway=0.25,
            tf_prior=0.10, tissue=0.10, ortholog=0.15,
            gold_threshold=0.80, silver_threshold=0.50,
            min_evidence_sources=3,
        )
