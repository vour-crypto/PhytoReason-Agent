"""
confidence_evaluator.py — 结果可信度评估器。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class EvaluationResult(BaseModel):
    overall_score: float = Field(default=0.0, ge=0.0, le=1.0)
    label: str = "not_assessed"
    n_candidates: int = 0
    n_evidence_sources: int = 0
    n_contradictions: int = 0
    convergence_score: float = 0.0
    data_quality_score: float = 0.0
    explanation: str = ""


class ConfidenceEvaluator:
    """结果可信度评估器。"""

    @staticmethod
    def evaluate(
        n_candidates: int = 0,
        n_high_confidence: int = 0,
        n_evidence_sources: int = 0,
        n_contradictions: int = 0,
        convergence_score: float = 0.0,
        data_quality_score: float = 1.0,
        sample_count: int = 0,
    ) -> EvaluationResult:
        score = 0.0
        parts: list[str] = []

        cand_score = 0.0
        if n_candidates >= 10: cand_score = 0.20
        elif n_candidates >= 3: cand_score = 0.15
        elif n_candidates >= 1: cand_score = 0.08
        if n_high_confidence >= 1: cand_score += 0.10
        elif n_high_confidence >= 3: cand_score += 0.15
        cand_score = min(cand_score, 0.30)
        score += cand_score
        parts.append(f"候选TF: {cand_score:.2f}")

        ev_score = min(n_evidence_sources * 0.06, 0.30)
        score += ev_score
        parts.append(f"证据多样性: {ev_score:.2f}")

        conv_score = convergence_score * 0.15
        score += conv_score
        parts.append(f"收敛度: {conv_score:.2f}")

        contradiction_penalty = min(n_contradictions * 0.05, 0.15)
        score -= contradiction_penalty
        if contradiction_penalty > 0:
            parts.append(f"矛盾惩罚: -{contradiction_penalty:.2f}")

        dq_score = data_quality_score * 0.10
        score += dq_score
        parts.append(f"数据质量: {dq_score:.2f}")

        if sample_count < 4:
            score *= 0.5
            parts.append(f"样本不足惩罚: ×0.5")
        elif sample_count < 10:
            score *= 0.8
            parts.append(f"样本有限调整: ×0.8")

        score = max(0.0, min(score, 1.0))

        if score >= 0.70: label = "high"
        elif score >= 0.40: label = "moderate"
        elif score >= 0.10: label = "low"
        else: label = "unreliable"

        return EvaluationResult(
            overall_score=round(score, 3), label=label,
            n_candidates=n_candidates, n_evidence_sources=n_evidence_sources,
            n_contradictions=n_contradictions, convergence_score=convergence_score,
            data_quality_score=data_quality_score, explanation="; ".join(parts),
        )
