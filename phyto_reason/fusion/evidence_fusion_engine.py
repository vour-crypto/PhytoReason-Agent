"""
evidence_fusion_engine.py — 证据融合引擎主入口。

整合所有子模块，对 CandidateGene + 推理结果执行完整的
多证据融合、矛盾检测、校准和排序。
"""

from __future__ import annotations

from phyto_reason.fusion.fusion_result import (
    ContradictionFlag, ContradictionLevel, EvidenceHierarchy,
    FusionEvidence, UnifiedFusionResult,
)
from phyto_reason.fusion.scoring_policy import ScoringPolicy, DEFAULT_POLICY
from phyto_reason.fusion.contradiction_detector import ContradictionDetector
from phyto_reason.fusion.evidence_ranker import EvidenceRanker
from phyto_reason.fusion.confidence_calibrator import ConfidenceCalibrator


class EvidenceFusionEngine:
    """证据融合引擎。

    对单个基因执行完整证据融合流水线:
      1. 收集各维度证据分数
      2. 统一归一化
      3. 加权融合
      4. 矛盾检测
      5. 校准
      6. 构建 UnifiedFusionResult
    """

    def __init__(self, policy: ScoringPolicy | None = None) -> None:
        self.policy = policy or DEFAULT_POLICY
        self.detector = ContradictionDetector()
        self.ranker = EvidenceRanker()
        self.calibrator = ConfidenceCalibrator(policy)

    def fuse(
        self,
        gene_id: str,
        correlation_score: float = 0.0,
        module_membership: float = 0.0,
        motif_score: float = 0.0,
        pathway_score: float = 0.0,
        tf_prior_score: float = 0.0,
        tissue_score: float = 0.0,
        ortholog_score: float = 0.0,
        literature_score: float = 0.0,
        biological_plausibility: str | None = None,
    ) -> UnifiedFusionResult:
        """执行单基因完整证据融合。

        Args:
            gene_id: 基因ID
            correlation_score: 相关性分数 (0-1)
            module_membership: 模块归属分数 (0-1)
            motif_score: Motif结合分数 (0-1)
            pathway_score: 通路一致性分数 (0-1)
            tf_prior_score: TF家族先验分数 (0-1)
            tissue_score: 组织特异性分数 (0-1)
            ortholog_score: 同源保守性分数 (0-1)
            literature_score: 文献支持分数 (0-1)
            biological_plausibility: 生物学可信度字符串

        Returns:
            UnifiedFusionResult
        """
        # 1. 收集证据维度
        raw_scores: dict[str, float] = {
            "correlation": correlation_score,
            "module_membership": module_membership,
            "motif": motif_score,
            "pathway": pathway_score,
            "tf_prior": tf_prior_score,
            "tissue": tissue_score,
            "ortholog": ortholog_score,
            "literature": literature_score,
        }

        # 2. 归一化 & 构建 FusionEvidence
        evidence_list = self._build_evidence_list(raw_scores)

        # 3. 加权融合
        raw_fused = self._compute_raw_fused(evidence_list)

        # 4. 矛盾检测
        contradictions = self.detector.detect(
            correlation_score=correlation_score,
            motif_score=motif_score,
            pathway_score=pathway_score,
            tf_prior_score=tf_prior_score,
            tissue_score=tissue_score,
            ortholog_score=ortholog_score,
            module_membership=module_membership,
            biological_plausibility=biological_plausibility,
            evidence_list=evidence_list,
        )

        # 5. 证据排序 + 层次
        hierarchy = self.ranker.build_hierarchy(evidence_list)

        # 6. 置信度校准
        calibrated_score, confidence_level, explanation = self.calibrator.calibrate(
            raw_fused_score=raw_fused,
            evidence_list=evidence_list,
            contradictions=contradictions,
            hierarchy=hierarchy,
        )

        # 7. 构建结果
        return UnifiedFusionResult(
            target_gene_id=gene_id,
            raw_fused_score=round(raw_fused, 4),
            calibrated_score=round(calibrated_score, 4),
            confidence_level=confidence_level,
            evidence_list=evidence_list,
            contradictions=contradictions,
            hierarchy=hierarchy,
            explanation=explanation,
        )

    def _build_evidence_list(
        self, raw_scores: dict[str, float],
    ) -> list[FusionEvidence]:
        """构建FusionEvidence列表（归一化 + 权重赋值）。"""
        active = {k: v for k, v in raw_scores.items() if v > 0}
        if not active:
            return []

        evidence_list: list[FusionEvidence] = []
        for source, raw_score in active.items():
            weight = self.policy.get_weight(source)
            norm = min(raw_score, 1.0)  # already 0-1

            evidence_list.append(FusionEvidence(
                source=source,
                raw_score=raw_score,
                normalized_score=norm,
                weight=weight,
                contribution=norm * weight,
            ))

        return evidence_list

    def _compute_raw_fused(
        self, evidence_list: list[FusionEvidence],
    ) -> float:
        """加权融合: Σ(score × weight) / Σ(weight)。"""
        if not evidence_list:
            return 0.0

        total_weight = sum(e.weight for e in evidence_list)
        if total_weight == 0:
            return 0.0

        weighted = sum(e.normalized_score * e.weight for e in evidence_list)
        return weighted / total_weight
