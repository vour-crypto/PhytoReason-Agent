"""
plausibility_scorer.py — 生物学可信度评分器。

融合多个推理维度的分数，输出最终的生物学可信度评估。
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from phyto_reason.models.enums import BiologicalPlausibility
from phyto_reason.reasoning.reasoning_result import ReasoningResult


class PlausibilityWeights(BaseModel):
    """每个推理维度的权重配置。"""
    tf_prior: float = 0.25
    pathway: float = 0.25
    tissue: float = 0.20
    ortholog: float = 0.15
    motif: float = 0.15

    gold_threshold: float = 0.70
    silver_threshold: float = 0.40


DEFAULT_WEIGHTS = PlausibilityWeights()


class PlausibilityScore(BaseModel):
    """最终可信度评分。"""
    fused_score: float = Field(default=0.0, ge=0.0, le=1.0)
    biological_plausibility: BiologicalPlausibility = BiologicalPlausibility.UNKNOWN
    dimension_scores: dict[str, float] = Field(default_factory=dict)
    explanation: str = ""
    convergence_bonus: float = 0.0


class PlausibilityScorer:
    """生物学可信度评分器。

    融合多个推理器输出的 ReasoningResult，生成统一的
    PlausibilityScore，包含 fused_score 和 biological_plausibility。
    """

    def __init__(self, weights: PlausibilityWeights | None = None) -> None:
        self.weights = weights or DEFAULT_WEIGHTS

    def score(
        self,
        tf_prior_result: ReasoningResult | None = None,
        pathway_result: ReasoningResult | None = None,
        tissue_result: ReasoningResult | None = None,
        ortholog_result: ReasoningResult | None = None,
        motif_score: float | None = None,
    ) -> PlausibilityScore:
        """融合多维度推理结果。

        Args:
            tf_prior_result: TF家族先验推理结果
            pathway_result: 通路一致性推理结果
            tissue_result: 组织特异性推理结果
            ortholog_result: 同源基因推理结果
            motif_score: Motif结合分数 (来自外部)

        Returns:
            PlausibilityScore
        """
        dim_scores: dict[str, float] = {}
        dim_sum = 0.0
        weight_sum = 0.0

        # 1. TF prior
        if tf_prior_result is not None:
            dim_scores["tf_prior"] = tf_prior_result.score
            dim_sum += tf_prior_result.score * self.weights.tf_prior
            weight_sum += self.weights.tf_prior

        # 2. Pathway
        if pathway_result is not None:
            dim_scores["pathway"] = pathway_result.score
            dim_sum += pathway_result.score * self.weights.pathway
            weight_sum += self.weights.pathway

        # 3. Tissue
        if tissue_result is not None:
            dim_scores["tissue"] = tissue_result.score
            dim_sum += tissue_result.score * self.weights.tissue
            weight_sum += self.weights.tissue

        # 4. Ortholog
        if ortholog_result is not None:
            dim_scores["ortholog"] = ortholog_result.score
            dim_sum += ortholog_result.score * self.weights.ortholog
            weight_sum += self.weights.ortholog

        # 5. Motif (from external)
        if motif_score is not None:
            dim_scores["motif"] = motif_score
            dim_sum += motif_score * self.weights.motif
            weight_sum += self.weights.motif

        if weight_sum == 0:
            return PlausibilityScore(
                fused_score=0.0,
                biological_plausibility=BiologicalPlausibility.UNKNOWN,
                explanation="无可用推理数据。",
            )

        fused = dim_sum / weight_sum

        # 收敛加分: 多个独立证据得分都较高 → 额外加分
        high_scores = [s for s in dim_scores.values() if s >= 0.5]
        convergence_bonus = 0.0
        if len(high_scores) >= 2:
            convergence_bonus = min(len(high_scores) * 0.03, 0.12)
        elif len(high_scores) >= 1 and fused >= 0.5:
            convergence_bonus = 0.02

        fused = min(fused + convergence_bonus, 1.0)

        # 判定 final 可信度
        if fused >= self.weights.gold_threshold:
            plausibility = BiologicalPlausibility.PLAUSIBLE
            label = "生物学合理 (Biologically Plausible)"
        elif fused >= self.weights.silver_threshold:
            plausibility = BiologicalPlausibility.WEAK
            label = "弱支持 (Weakly Supported)"
        else:
            plausibility = BiologicalPlausibility.UNLIKELY
            label = "可能性低 (Unlikely)"

        dim_str = ", ".join(
            f"{k}={v:.2f}" for k, v in sorted(dim_scores.items())
        )

        return PlausibilityScore(
            fused_score=round(fused, 3),
            biological_plausibility=plausibility,
            dimension_scores=dim_scores,
            convergence_bonus=round(convergence_bonus, 3),
            explanation=(
                f"{label} (融合分: {fused:.2f}, 维度: [{dim_str}], "
                f"收敛加分: {convergence_bonus:.3f})"
            ),
        )
