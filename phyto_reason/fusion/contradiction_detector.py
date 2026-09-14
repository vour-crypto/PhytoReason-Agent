"""
contradiction_detector.py — 证据矛盾检测器。

识别不同证据维度之间的相互矛盾，生成 ContradictionFlag。
规则均为 deterministic。
"""

from __future__ import annotations

from phyto_reason.models.enums import BiologicalPlausibility
from phyto_reason.fusion.fusion_result import (
    ContradictionFlag, ContradictionLevel, FusionEvidence,
)


class ContradictionDetector:
    """证据矛盾检测器。

    检测以下矛盾模式:
      1. high_correlation_no_motif    — 高相关但无motif支持
      2. high_prior_wrong_tissue      — 强先验但组织不匹配
      3. pathway_no_upstream          — 通路分高但上游酶缺失
      4. strong_ortholog_wrong_path   — 强同源但通路不匹配
      5. high_module_low_corr         — 高模块归属但低相关性
      6. plausible_but_no_evidence    — 生物学合理但无数据支持
    """

    @staticmethod
    def detect(
        correlation_score: float = 0.0,
        motif_score: float = 0.0,
        pathway_score: float = 0.0,
        tf_prior_score: float = 0.0,
        tissue_score: float = 0.0,
        ortholog_score: float = 0.0,
        module_membership: float = 0.0,
        biological_plausibility: BiologicalPlausibility | str | None = None,
        evidence_list: list[FusionEvidence] | None = None,
    ) -> list[ContradictionFlag]:
        """检测所有矛盾模式。

        Args:
            各证据维度的分数 (0.0-1.0)
            biological_plausibility: 生物学可信度
            evidence_list: 已处理的证据列表 (可选)

        Returns:
            发现的 ContradictionFlag 列表
        """
        flags: list[ContradictionFlag] = []

        # 1. high_correlation_no_motif
        flags.extend(ContradictionDetector._rule_correlation_no_motif(
            correlation_score, motif_score
        ))

        # 2. high_prior_wrong_tissue
        flags.extend(ContradictionDetector._rule_prior_wrong_tissue(
            tf_prior_score, tissue_score
        ))

        # 3. pathway_no_upstream
        flags.extend(ContradictionDetector._rule_pathway_no_upstream(
            pathway_score, correlation_score, module_membership
        ))

        # 4. strong_ortholog_wrong_path
        flags.extend(ContradictionDetector._rule_ortholog_wrong_path(
            ortholog_score, pathway_score, tf_prior_score
        ))

        # 5. high_module_low_corr
        flags.extend(ContradictionDetector._rule_module_no_corr(
            module_membership, correlation_score
        ))

        # 6. plausible_but_no_evidence
        flags.extend(ContradictionDetector._rule_plausible_no_evidence(
            biological_plausibility, evidence_list
        ))

        return flags

    # ── Rules ──────────────────────────────────────────────

    @staticmethod
    def _rule_correlation_no_motif(
        corr: float, motif: float,
    ) -> list[ContradictionFlag]:
        """高相关但无motif → TF可能只是共表达, 并非直接调控。"""
        if corr >= 0.70 and motif < 0.20:
            return [ContradictionFlag(
                level=ContradictionLevel.HIGH,
                pattern="high_correlation_no_motif",
                description="高相关性但缺乏启动子motif结合证据，"
                           "TF-代谢物关联可能为间接共表达而非直接调控。",
                evidence_a=f"correlation={corr:.2f}",
                evidence_b=f"motif={motif:.2f}",
                penalty=0.20,
            )]
        return []

    @staticmethod
    def _rule_prior_wrong_tissue(
        prior: float, tissue: float,
    ) -> list[ContradictionFlag]:
        """强先验但组织不匹配 → TF在该组织中可能不发挥功能。"""
        if prior >= 0.70 and tissue < 0.25:
            return [ContradictionFlag(
                level=ContradictionLevel.MEDIUM,
                pattern="high_prior_wrong_tissue",
                description="TF家族先验支持强，但组织表达模式与代谢物积累组织不一致。",
                evidence_a=f"tf_prior={prior:.2f}",
                evidence_b=f"tissue={tissue:.2f}",
                penalty=0.15,
            )]
        return []

    @staticmethod
    def _rule_pathway_no_upstream(
        pathway: float, corr: float, module: float,
    ) -> list[ContradictionFlag]:
        """通路分高但无共表达支持 → 通路推理可能来自先验而非数据。"""
        if pathway >= 0.60 and corr < 0.30 and module < 0.30:
            return [ContradictionFlag(
                level=ContradictionLevel.MEDIUM,
                pattern="pathway_no_upstream",
                description="通路一致性评分高但缺乏共表达支持，"
                           "推理结果可能过度依赖先验知识。",
                evidence_a=f"pathway={pathway:.2f}",
                evidence_b=f"corr={corr:.2f} module={module:.2f}",
                penalty=0.10,
            )]
        return []

    @staticmethod
    def _rule_ortholog_wrong_path(
        ortholog: float, pathway: float, prior: float,
    ) -> list[ContradictionFlag]:
        """强同源但通路不匹配 → 同源TF在目标物种中可能调控不同通路。"""
        if ortholog >= 0.70 and pathway < 0.30 and prior < 0.30:
            return [ContradictionFlag(
                level=ContradictionLevel.LOW,
                pattern="strong_ortholog_wrong_path",
                description="拟南芥同源证据强但在当前物种中通路支持弱，"
                           "提示可能的功能分化。",
                evidence_a=f"ortholog={ortholog:.2f}",
                evidence_b=f"pathway={pathway:.2f} prior={prior:.2f}",
                penalty=0.05,
            )]
        return []

    @staticmethod
    def _rule_module_no_corr(
        module: float, corr: float,
    ) -> list[ContradictionFlag]:
        """高模块归属但低相关性 → 基因在模块中但与代谢物无直接关联。"""
        if module >= 0.70 and corr < 0.20:
            return [ContradictionFlag(
                level=ContradictionLevel.LOW,
                pattern="high_module_low_corr",
                description="高模块归属但TF与目标代谢物无显著相关性。",
                evidence_a=f"module_membership={module:.2f}",
                evidence_b=f"correlation={corr:.2f}",
                penalty=0.05,
            )]
        return []

    @staticmethod
    def _rule_plausible_no_evidence(
        plausibility: BiologicalPlausibility | str | None,
        evidence_list: list[FusionEvidence] | None,
    ) -> list[ContradictionFlag]:
        """生物学合理但无数据支持。"""
        if plausibility is None:
            return []

        plaus_str = plausibility.value if isinstance(plausibility, BiologicalPlausibility) else str(plausibility)

        if plaus_str in ("plausible",) and evidence_list:
            data_sources = [e for e in evidence_list if e.raw_score >= 0.30]
            if len(data_sources) < 2:
                return [ContradictionFlag(
                    level=ContradictionLevel.MEDIUM,
                    pattern="plausible_but_no_evidence",
                    description="生物学合理性判断为plausible但缺少足够的数据维度支持。",
                    evidence_a=f"plausibility={plaus_str}",
                    evidence_b=f"data_sources={len(data_sources)}",
                    penalty=0.10,
                )]
        return []
