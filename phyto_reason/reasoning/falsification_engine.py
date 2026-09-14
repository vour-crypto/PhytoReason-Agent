"""
falsification_engine.py — 主动反驳假设。

核心原则: 系统必须学会攻击自己的假设。
不是找"为什么对"，而是找"为什么可能是错的"。

Layer 4: 仅被 hypothesis reasoning 调用。
禁止: 直接访问统计工具或原始数据。
"""

from __future__ import annotations

import logging

from phyto_reason.models.mechanistic_hypothesis import (
    MechanisticHypothesis, MechanismType, EvidenceRecord, EvidenceRelationType,
    FalsificationTest,
)

logger = logging.getLogger("falsification_engine")


class FalsificationEngine:
    """主动寻找反驳证据的引擎。

    每个假设生成后必须经过 falsification 检查。
    """

    def attack(
        self,
        hypothesis: MechanisticHypothesis,
        expression_data: dict | None = None,
        metabolite_data: dict | None = None,
        n_samples: int = 0,
        species: str = "",
        tissues: list[str] | None = None,
        gene_ids_in_data: list[str] | None = None,
    ) -> MechanisticHypothesis:
        """对假设执行系统性的反驳检查。

        v3.0: 必须检查 species/tissue/condition mismatch。
        输出 why_wrong 说明。

        Args:
            gene_ids_in_data: 表达数据中存在的基因 ID 列表（已处理的列表，非原始矩阵）
        """
        if not hypothesis.proposed_regulators:
            return hypothesis

        contradictions: list[EvidenceRecord] = []

        contradictions.extend(self._check_expression_contradiction(
            hypothesis, gene_ids_in_data or (list(expression_data) if expression_data else None)))
        contradictions.extend(self._check_sample_size(hypothesis, n_samples))
        contradictions.extend(self._check_single_evidence(hypothesis))
        contradictions.extend(self._check_tissue_mismatch(hypothesis, tissues))
        contradictions.extend(self._check_species_mismatch(hypothesis, species))
        contradictions.extend(self._check_causal_gap(hypothesis))

        for c in contradictions:
            hypothesis.add_contradiction(c)

        hypothesis.falsification_tests = self._propose_falsification_tests(hypothesis)

        n_support = len(hypothesis.supporting_evidence)
        n_contra = len(hypothesis.contradictory_evidence)
        if n_contra >= n_support and n_support > 0:
            hypothesis.uncertainty_level = "high"

        # Populate why_wrong
        hypothesis.why_wrong = self._build_why_wrong(hypothesis, contradictions)

        return hypothesis

    def _check_expression_contradiction(
        self,
        hypothesis: MechanisticHypothesis,
        gene_ids_in_data: list[str] | None,
    ) -> list[EvidenceRecord]:
        """检查表达数据中的反驳证据。

        仅接收已处理的基因 ID 列表，不接收原始矩阵。
        """
        results = []
        if not gene_ids_in_data:
            return results

        regulators = hypothesis.proposed_regulators
        available = [
            r for r in regulators
            if r in gene_ids_in_data
        ]

        if not available:
            results.append(EvidenceRecord(
                source="falsification_engine",
                evidence_type="expression_check",
                score=0.0,
                relation_type=EvidenceRelationType.CONTRADICTS,
                description="候选调控因子在表达数据中未找到",
            ))
        return results

    def _check_sample_size(
        self,
        hypothesis: MechanisticHypothesis,
        n_samples: int,
    ) -> list[EvidenceRecord]:
        """小样本导致结论不可靠（通过 INSUFFICIENT 表达统计效力不足，不归为矛盾）。"""
        if n_samples < 6 and hypothesis.proposed_regulators:
            return [EvidenceRecord(
                source="falsification_engine",
                evidence_type="sample_size_check",
                score=0.0,
                relation_type=EvidenceRelationType.INSUFFICIENT,
                description=f"小样本量 (n={n_samples})，调控推断统计效力不足",
                contradiction_severity="low",
            )]
        return []

    def _check_single_evidence(
        self,
        hypothesis: MechanisticHypothesis,
    ) -> list[EvidenceRecord]:
        """单一证据源的假设不可靠。"""
        if len(hypothesis.supporting_evidence) < 2:
            return [EvidenceRecord(
                source="falsification_engine",
                evidence_type="evidence_diversity_check",
                score=0.0,
                relation_type=EvidenceRelationType.CONTRADICTS,
                description=f"仅 {len(hypothesis.supporting_evidence)} 个证据源，结论不可靠",
            )]
        return []

    def _check_tissue_mismatch(
        self,
        hypothesis: MechanisticHypothesis,
        tissues: list[str] | None,
    ) -> list[EvidenceRecord]:
        """检查组织错配。"""
        if not tissues:
            return []
        results = []
        if hypothesis.mechanism_type == MechanismType.TRANSCRIPTIONAL:
            results.append(EvidenceRecord(
                source="falsification_engine",
                evidence_type="tissue_mismatch_check",
                score=0.0,
                relation_type=EvidenceRelationType.CONTEXT_DEPENDENT,
                description=f"当前数据涵盖组织: {', '.join(tissues[:3])}。"
                            f"缺少目标组织独立验证数据，组织特异性调控的可能性未被排除",
                contradiction_severity="low",
            ))
        return results

    def _check_species_mismatch(
        self,
        hypothesis: MechanisticHypothesis,
        species: str,
    ) -> list[EvidenceRecord]:
        """检查物种错配。"""
        if not species:
            return []
        results = []
        for ev in hypothesis.supporting_evidence:
            if ev.species_scope and ev.species_scope != species:
                results.append(EvidenceRecord(
                    source="falsification_engine",
                    evidence_type="species_mismatch_check",
                    score=0.0,
                    relation_type=EvidenceRelationType.SPECIES_LIMITED,
                    description=f"证据 '{ev.source}' 来自 {ev.species_scope}，"
                                f"当前分析物种为 {species}，跨物种推断存在不确定性",
                    contradiction_severity="medium",
                ))
        return results

    def _check_causal_gap(
        self,
        hypothesis: MechanisticHypothesis,
    ) -> list[EvidenceRecord]:
        """检查因果链中的已知断裂。"""
        results = []
        if hypothesis.causal_chain and hypothesis.causal_chain.known_gaps:
            for gap in hypothesis.causal_chain.known_gaps:
                results.append(EvidenceRecord(
                    source="falsification_engine",
                    evidence_type="causal_gap_check",
                    score=0.0,
                    relation_type=EvidenceRelationType.INSUFFICIENT,
                    description=f"因果链断裂: {gap}",
                    contradiction_severity="high",
                ))
        return results

    def _build_why_wrong(
        self,
        hypothesis: MechanisticHypothesis,
        contradictions: list[EvidenceRecord],
    ) -> str:
        """生成"为什么这个假设可能是错的"说明。"""
        if not contradictions:
            return ""
        parts = ["该假设存在以下潜在问题:"]
        for c in contradictions[:5]:
            sev = f"[{c.contradiction_severity or 'info'}]"
            parts.append(f"  {sev} {c.description[:100]}")
        parts.append("建议在验证前先解决上述问题。")
        return "\n".join(parts)

    def _propose_falsification_tests(
        self,
        hypothesis: MechanisticHypothesis,
    ) -> list[FalsificationTest]:
        """生成可执行的证伪测试。"""
        tests = []
        if hypothesis.mechanism_type == MechanismType.TRANSCRIPTIONAL:
            tests.append(FalsificationTest(
                test_name="TF_knockout_prediction",
                prediction="敲除候选TF后目标代谢物应下降",
                required_data="CRISPR/RNAi 敲除株的代谢物数据",
                contradictory_prediction="代谢物无变化",
                feasibility="low",
            ))
            tests.append(FalsificationTest(
                test_name="temporal_consistency",
                prediction="TF表达变化应先于代谢物积累",
                required_data="时序表达数据",
                contradictory_prediction="代谢物变化先于TF或无时序关系",
                feasibility="medium",
            ))
        return tests
