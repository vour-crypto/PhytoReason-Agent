"""
biological_reasoner.py — 生物学推理引擎入口。

统筹所有子推理器，对候选基因进行多维度生物学合理性评估。
"""

from __future__ import annotations

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.reasoning.tf_prior_reasoner import TfPriorReasoner
from phyto_reason.reasoning.pathway_reasoner import PathwayReasoner
from phyto_reason.reasoning.tissue_reasoner import TissueReasoner
from phyto_reason.reasoning.ortholog_reasoner import OrthologReasoner
from phyto_reason.reasoning.plausibility_scorer import PlausibilityScorer
from phyto_reason.reasoning.reasoning_result import ReasoningResult


class BiologicalReasoner:
    """生物学推理引擎主入口。

    对一个或多个 CandidateGene 执行完整的生物学推理流水线:
      1. TF家族先验 (TfPriorReasoner)
      2. 通路一致性 (PathwayReasoner)
      3. 组织特异性 (TissueReasoner)
      4. 同源保守性 (OrthologReasoner)
      5. 可信度融合 (PlausibilityScorer)
    """

    def __init__(
        self,
        weights: PlausibilityScore | None = None,
    ) -> None:
        self.scorer = PlausibilityScorer()

    def reason(
        self,
        gene: CandidateGene,
        target_metabolite: str,
        upstream_genes_found: list[str] | None = None,
        bottleneck_hit: bool | None = None,
        tf_expression: dict[str, float] | None = None,
        has_ortholog: bool | None = None,
        ortholog_identity: float | None = None,
    ) -> CandidateGene:
        """对单个候选基因执行完整生物学推理。

        Args:
            gene: 待评估的候选基因
            target_metabolite: 目标代谢物名称或类别
            upstream_genes_found: 在表达数据中发现的通路上游酶
            bottleneck_hit: 是否命中限速酶
            tf_expression: 该基因的组织表达谱 {tissue: value}
            has_ortholog: 是否存在拟南芥同源基因
            ortholog_identity: 同源基因序列一致性

        Returns:
            更新的 CandidateGene (包含推理分数和evidence)
        """
        tf_family = gene.tf_family

        # 1. TF家族先验
        prior_result = TfPriorReasoner.reason(tf_family, target_metabolite)

        # 2. 通路一致性
        pathway_result = PathwayReasoner.reason(
            tf_family=tf_family,
            target_metabolite=target_metabolite,
            upstream_genes_found=upstream_genes_found,
            bottleneck_hit=bottleneck_hit,
        )

        # 3. 组织特异性
        if tf_expression:
            tissue_result = TissueReasoner.reason_by_expression(
                target_metabolite, tf_expression
            )
        else:
            tissue_result = TissueReasoner.reason_by_prior(target_metabolite)

        # 4. 同源保守性
        ortholog_result = OrthologReasoner.reason(
            tf_family=tf_family,
            target_metabolite=target_metabolite,
            has_ortholog=has_ortholog,
            ortholog_identity=ortholog_identity,
        )

        # 5. 可信度融合
        plausibility = self.scorer.score(
            tf_prior_result=prior_result,
            pathway_result=pathway_result,
            tissue_result=tissue_result,
            ortholog_result=ortholog_result,
            motif_score=gene.motif_score if gene.motif_score > 0 else None,
        )

        # 6. 更新 CandidateGene
        self._apply_to_gene(gene, prior_result, pathway_result, tissue_result,
                            ortholog_result, plausibility)

        return gene

    def reason_batch(
        self,
        genes: list[CandidateGene],
        target_metabolite: str,
        upstream_genes_found: list[str] | None = None,
        tf_expression_map: dict[str, dict[str, float]] | None = None,
    ) -> list[CandidateGene]:
        """批量推理多个候选基因。"""
        results: list[CandidateGene] = []
        for gene in genes:
            expression = None
            if tf_expression_map:
                expression = tf_expression_map.get(gene.gene_id)

            result = self.reason(
                gene=gene,
                target_metabolite=target_metabolite,
                upstream_genes_found=upstream_genes_found,
                tf_expression=expression,
            )
            results.append(result)

        # 按融合分数降序排列
        results.sort(key=lambda g: g.confidence_score, reverse=True)
        return results

    @staticmethod
    def _apply_to_gene(
        gene: CandidateGene,
        prior_result: ReasoningResult,
        pathway_result: ReasoningResult,
        tissue_result: ReasoningResult,
        ortholog_result: ReasoningResult,
        plausibility: PlausibilityScore,
    ) -> None:
        """将推理结果写入 CandidateGene 对象。"""
        gene.pathway_score = pathway_result.score
        gene.tissue_specificity_score = tissue_result.score
        gene.ortholog_score = ortholog_result.score

        gene.literature_score = max(
            gene.literature_score, prior_result.score
        )

        gene.biological_plausibility = plausibility.biological_plausibility
        gene.confidence_score = plausibility.fused_score
        gene._assign_level()

        etype_map = {
            "tf_prior": EvidenceType.TF_FAMILY_PRIOR,
            "pathway": EvidenceType.PATHWAY_CONSISTENCY,
            "tissue": EvidenceType.TISSUE_SPECIFICITY,
            "ortholog": EvidenceType.ORTHOLOG,
        }

        for result, label in [
            (prior_result, "tf_prior"),
            (pathway_result, "pathway"),
            (tissue_result, "tissue"),
            (ortholog_result, "ortholog"),
        ]:
            etype = etype_map[label]
            for ev in result.evidence_list:
                gene.add_evidence(Evidence(
                    evidence_type=etype,
                    source=f"{label}_reasoner",
                    score=ev.score,
                    description=ev.description,
                ))
