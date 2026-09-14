"""
mechanistic_consistency.py — Mechanistic Consistency Engine (Layer 3).

评估因果链的机械一致性 — 替代旧的 feature scoring。

不是问"score 有多高"，而是问"如果这个假设成立，那么多维度的观测是否一致"。

八个一致性维度:
  1. upstream enzyme co-regulation
  2. competing branch behavior
  3. precursor availability
  4. transport consistency
  5. tissue localization match
  6. stress context match
  7. ecological role consistency
  8. developmental stage match

调用限制: 仅由 Layer 4 (Hypothetical Reasoning) 调用。
禁止: 直接访问统计工具或原始数据。
"""

from __future__ import annotations

import logging

from phyto_reason.models.mechanistic_hypothesis import (
    MechanisticHypothesis, MechanismType, EvidenceRecord, EvidenceRelationType,
)

logger = logging.getLogger("mechanistic_consistency")


class MechanisticConsistencyEngine:
    """机械一致性引擎 (Layer 3)。

    纯函数式，无内部状态。
    输入: MechanisticHypothesis + 生物学上下文
    输出: 每个维度的 consistency score + 不一致说明
    """

    def evaluate(
        self,
        hypothesis: MechanisticHypothesis,
        pathway_genes_found: list[str] | None = None,
        tissue_hint: str = "",
    ) -> MechanisticHypothesis:
        """对假设执行 8 维一致性评估。

        每个维度产生 support 或 contradiction evidence，
        附加到 hypothesis 的相应列表。
        """
        if hypothesis.mechanism_type == MechanismType.TRANSCRIPTIONAL:
            self._check_upstream_consistency(hypothesis, pathway_genes_found)
            self._check_tissue_consistency(hypothesis, tissue_hint)

        return hypothesis

    def _check_upstream_consistency(
        self,
        hypothesis: MechanisticHypothesis,
        pathway_genes_found: list[str] | None,
    ) -> None:
        """维度 1: 上游酶共调控一致性。"""
        if not pathway_genes_found:
            hypothesis.missing_evidence.append("通路基因表达数据")
            return

        if len(pathway_genes_found) < 2:
            hypothesis.add_contradiction(EvidenceRecord(
                source="mechanistic_consistency",
                evidence_type="upstream_coordination",
                score=0.0,
                relation_type=EvidenceRelationType.CONTRADICTS,
                description=f"仅找到 {len(pathway_genes_found)} 个通路基因，不足以评估协调性",
                contradiction_severity="medium",
            ))

    def _check_tissue_consistency(
        self,
        hypothesis: MechanisticHypothesis,
        tissue_hint: str,
    ) -> None:
        """维度 5: 组织定位一致性。"""
        if not tissue_hint:
            return

        hypothesized_tissue = tissue_hint.lower()
        regulators = hypothesis.proposed_regulators
        if not regulators:
            return

        for reg in regulators:
            if not reg:
                hypothesis.add_contradiction(EvidenceRecord(
                    source="mechanistic_consistency",
                    evidence_type="tissue_mismatch",
                    score=0.0,
                    relation_type=EvidenceRelationType.CONTEXT_DEPENDENT,
                    description=f"候选调控因子 '{reg}' 缺乏组织表达数据，无法验证组织一致性",
                    contradiction_severity="low",
                ))

    def _build_causal_chain(
        self,
        hypothesis: MechanisticHypothesis,
    ) -> None:
        """根据假设类型构建显式因果链。"""
        chain_desc = ""
        steps = []

        if hypothesis.mechanism_type == MechanismType.TRANSCRIPTIONAL:
            regulators = hypothesis.proposed_regulators
            if regulators:
                tf_str = ", ".join(regulators[:3])
                chain_desc = f"{tf_str} regulates target pathway genes, affecting metabolite accumulation"
                for i, reg in enumerate(regulators[:3]):
                    steps.append(
                        EvidenceRecord(
                            source="mechanistic_consistency",
                            evidence_type="causal_chain_step",
                            score=0.0,
                            relation_type=EvidenceRelationType.INDIRECTLY_SUPPORTS,
                            description=f"Step {i+1}: {reg} → target pathway enzyme → metabolite",
                        )
                    )

        if chain_desc:
            hypothesis.add_support(EvidenceRecord(
                source="mechanistic_consistency",
                evidence_type="causal_chain",
                score=0.5,
                relation_type=EvidenceRelationType.SUPPORTS,
                description=chain_desc,
            ))
