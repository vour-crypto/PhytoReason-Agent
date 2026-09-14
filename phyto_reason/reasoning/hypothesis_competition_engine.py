"""
hypothesis_competition_engine.py — 多假设竞争引擎。

系统必须生成多个 competing explanations 而非单一结论。
让 evidence 在假设之间竞争。
"""

from __future__ import annotations

import logging

from phyto_reason.models.mechanistic_hypothesis import (
    MechanisticHypothesis, MechanismType, EvidenceRecord, EvidenceRelationType,
)
from phyto_reason.ontology.metabolite_metadata import resolve_ontology

logger = logging.getLogger("hypothesis_competition_engine")


class HypothesisCompetitionEngine:
    """多假设竞争引擎。

    根据当前证据生成一系列 competing hypotheses，
    并评估每个假设的相对证据支持强度。
    """

    def generate_competing_hypotheses(
        self,
        target_metabolite: str = "",
        cv: float = 0.0,
        n_samples: int = 0,
        has_transcriptional_evidence: bool = False,
        top_tf_candidates: list[str] | None = None,
    ) -> list[MechanisticHypothesis]:
        """生成多个竞争假设。"""
        hypotheses: list[MechanisticHypothesis] = []
        onto = resolve_ontology(target_metabolite)

        # H1: 转录调控（默认假设）
        h1 = MechanisticHypothesis(
            id="H1",
            title="Transcriptional activation hypothesis",
            mechanism_type=MechanismType.TRANSCRIPTIONAL,
            target_metabolites=[target_metabolite] if target_metabolite else [],
            proposed_regulators=top_tf_candidates or [],
        )
        if has_transcriptional_evidence:
            h1.add_support(EvidenceRecord(
                source="competition_engine",
                evidence_type="transcriptional_evidence",
                score=0.4,
                relation_type=EvidenceRelationType.SUPPORTS,
                description="检测到共表达和 motif 证据，转录调控可能",
            ))
        else:
            h1.add_contradiction(EvidenceRecord(
                source="competition_engine",
                evidence_type="transcriptional_evidence",
                score=0.0,
                relation_type=EvidenceRelationType.CONTRADICTS,
                description="缺少转录调控的直接证据",
            ))
            h1.uncertainty_level = "high"
        hypotheses.append(h1)

        # H2: 胁迫响应（如果代谢物有胁迫关联）
        if onto and onto.stress_associations:
            h2 = MechanisticHypothesis(
                id="H2",
                title="Stress-induced metabolic redistribution",
                mechanism_type=MechanismType.STRESS_RESPONSE,
                target_metabolites=[target_metabolite] if target_metabolite else [],
                upstream_signal=f"stress response associated with {target_metabolite}",
            )
            h2.add_support(EvidenceRecord(
                source="metabolite_ontology",
                evidence_type="stress_association",
                score=0.5,
                relation_type=EvidenceRelationType.SUPPORTS,
                description=f"代谢物与胁迫响应相关: {', '.join(onto.stress_associations)}",
            ))
            h2.missing_evidence = ["胁迫响应 marker 基因表达数据"]
            hypotheses.append(h2)

        # H3: 转运重分布（如果代谢物有转运机制）
        if onto and onto.transport_modes:
            h3 = MechanisticHypothesis(
                id="H3",
                title="Transport-mediated redistribution",
                mechanism_type=MechanismType.TRANSPORT,
                target_metabolites=[target_metabolite] if target_metabolite else [],
            )
            h3.add_support(EvidenceRecord(
                source="metabolite_ontology",
                evidence_type="transport_evidence",
                score=0.3,
                relation_type=EvidenceRelationType.CONDITIONALLY_SUPPORTS,
                description=f"代谢物已知有转运机制，但需要转运蛋白表达数据",
            ))
            h3.missing_evidence = ["转运蛋白表达数据", "组织特异性积累数据"]
            hypotheses.append(h3)

        # H4: 前体限制（如果有前体途径信息）
        if onto and onto.precursor_pathways:
            h4 = MechanisticHypothesis(
                id="H4",
                title="Precursor limitation hypothesis",
                mechanism_type=MechanismType.PRECURSOR_LIMITATION,
                target_metabolites=[target_metabolite] if target_metabolite else [],
                upstream_signal=f"precursor ({', '.join(onto.precursor_pathways)}) availability",
            )
            h4.add_support(EvidenceRecord(
                source="metabolite_ontology",
                evidence_type="precursor_pathway",
                score=0.2,
                relation_type=EvidenceRelationType.CONDITIONALLY_SUPPORTS,
                description=f"前体途径存在，但需要前体丰度数据",
            ))
            h4.missing_evidence = ["前体代谢物丰度数据"]
            hypotheses.append(h4)

        return hypotheses

    def rank_by_evidence(
        self,
        hypotheses: list[MechanisticHypothesis],
    ) -> list[MechanisticHypothesis]:
        """按支持/反驳证据比率排序假设。

        v3.0: 使用 support ÷ (support + contra + missing) 比率，
        而非简单线性组合。比率 < 0.5 说明反驳证据过多。
        """
        scored = []
        for h in hypotheses:
            n_support = len(h.supporting_evidence)
            n_contra = len(h.contradictory_evidence)
            n_missing = len(h.missing_evidence)

            total = n_support + n_contra + n_missing
            if total == 0:
                ratio = 0.0
            else:
                ratio = n_support / max(total, 1)

            high_contra = sum(
                1 for c in h.contradictory_evidence
                if c.contradiction_severity == "high"
            )
            contra_penalty = min(high_contra * 0.2, 0.6)

            score = max(0.0, ratio - contra_penalty)
            scored.append((score, h))

        scored.sort(key=lambda x: x[0], reverse=True)
        return [h for _, h in scored]
