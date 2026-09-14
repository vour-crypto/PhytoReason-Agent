"""
tissue_reasoner.py — 组织特异性表达推理。

评估候选TF的表达模式与目标代谢物的积累组织是否一致。
"""

from __future__ import annotations

import math

from phyto_reason.reasoning.reasoning_result import PriorEvidence, ReasoningResult


# ── 常见植物组织的代谢物偏好知识 ──────────────────────────────

TISSUE_METABOLITE_BIAS: dict[str, list[str]] = {
    "root": [
        "alkaloid", "nicotine", "pyridine", "tropane",
        "isoquinoline", "benzylisoquinoline",
    ],
    "leaf": [
        "flavonoid", "flavonol", "anthocyanin",
        "terpenoid", "monoterpene", "diterpene",
    ],
    "stem": ["lignin", "cellulose", "terpenoid"],
    "flower": ["anthocyanin", "flavonoid", "terpenoid"],
    "fruit": ["anthocyanin", "flavonoid", "carotenoid"],
    "seed": ["flavonoid", "proanthocyanidin", "lignin"],
    "bark": ["lignin", "terpenoid", "phenolic"],
    "tuber": ["alkaloid", "terpenoid"],
}


class TissueReasoner:
    """组织特异性推理器。

    支持两种模式:
      1. 基于先验组织-代谢物偏好 (无表达数据时备用)
      2. 基于实际表达数据计算组织特异性指数 (Tau)
    """

    @staticmethod
    def get_bias_for_metabolite(metabolite: str) -> list[str]:
        """获取与指定代谢物相关的组织列表。"""
        m = metabolite.lower()
        matches: list[tuple[str, str]] = []
        for tissue, metabolites in TISSUE_METABOLITE_BIAS.items():
            for meta in metabolites:
                if meta in m or m in meta:
                    matches.append((tissue, meta))
        return list(set(t for t, _ in matches))

    @staticmethod
    def reason_by_prior(
        target_metabolite: str,
        tf_tissue: str | None = None,
    ) -> ReasoningResult:
        """基于先验知识的组织-代谢物偏好推理。

        Args:
            target_metabolite: 目标代谢物名称
            tf_tissue: TF高表达的组织 (如 "root", "leaf")

        Returns:
            ReasoningResult
        """
        related_tissues = TissueReasoner.get_bias_for_metabolite(target_metabolite)

        if not related_tissues:
            return ReasoningResult(
                score=0.30, label="low",
                explanation=f"代谢物 '{target_metabolite}' 无已知组织积累偏好。",
            )

        details = [
            f"代谢物 '{target_metabolite}' 相关组织: {', '.join(related_tissues)}",
        ]

        if tf_tissue:
            match = tf_tissue.lower() in related_tissues
            score = 0.85 if match else 0.25
            if match:
                details.append(
                    f"TF表达组织 '{tf_tissue}' 与代谢物积累组织一致 → 组织共定位证据强"
                )
                explanation = (
                    f"TF在 '{tf_tissue}' 中高表达，与已知的 '{target_metabolite}' "
                    f"积累组织一致，提示存在组织特异性调控关系。"
                )
            else:
                details.append(
                    f"TF表达组织 '{tf_tissue}' 与代谢物积累组织 '{related_tissues[0]}' "
                    f"不一致 → 组织共定位证据弱"
                )
                explanation = (
                    f"TF在 '{tf_tissue}' 中表达，但 '{target_metabolite}' 通常在 "
                    f"'{', '.join(related_tissues)}' 中积累，组织共定位证据不足。"
                )
        else:
            score = 0.50
            explanation = (
                f"'{target_metabolite}' 已知在 '{', '.join(related_tissues)}' 中积累。"
                "需要TF组织表达数据进行共定位分析。"
            )

        evidence = PriorEvidence(
            strength="moderate" if score >= 0.5 else "weak",
            score=score,
            description=f"组织-代谢物偏好 (相关组织: {', '.join(related_tissues)})",
            pmids=[],
        )

        return ReasoningResult(
            score=score,
            explanation=explanation,
            details=details,
            evidence_list=[evidence],
        )

    @staticmethod
    def reason_by_expression(
        target_metabolite: str,
        tf_expression: dict[str, float] | None = None,
    ) -> ReasoningResult:
        """基于表达数据的组织特异性推理 (Tau指数)。

        Args:
            target_metabolite: 目标代谢物名称
            tf_expression: 字典 {组织名: 表达量}, None=无数据

        Returns:
            ReasoningResult
        """
        if not tf_expression or len(tf_expression) < 2:
            return TissueReasoner.reason_by_prior(target_metabolite)

        # 计算 Tau 组织特异性指数
        max_expr = max(tf_expression.values())
        if max_expr <= 0:
            return ReasoningResult(
                score=0.0, label="none",
                explanation="表达数据无效 (最大表达量≤0)。",
            )

        n = len(tf_expression)
        tau_sum = 0.0
        for val in tf_expression.values():
            normalized = val / max_expr
            tau_sum += 1.0 - normalized
        tau = tau_sum / (n - 1) if n > 1 else 0.0

        # 找高表达组织
        sorted_tissues = sorted(
            tf_expression.items(), key=lambda x: x[1], reverse=True
        )
        top_tissue = sorted_tissues[0][0]

        related_tissues = TissueReasoner.get_bias_for_metabolite(target_metabolite)
        tissue_match = top_tissue.lower() in related_tissues if related_tissues else False

        # Tau: 0=广泛表达, 1=完全特异性
        if tau >= 0.80 and tissue_match:
            score = 0.90
            explanation = (
                f"TF在 '{top_tissue}' 中高度特异性表达 (Tau={tau:.2f})，"
                f"与代谢物积累组织一致，组织共定位证据强。"
            )
        elif tau >= 0.50 and tissue_match:
            score = 0.65
            explanation = (
                f"TF有一定组织特异性 (Tau={tau:.2f})，且在 '{top_tissue}' 中"
                f"与代谢物积累组织一致。"
            )
        elif tau >= 0.80:
            score = 0.40
            explanation = (
                f"TF组织特异性高 (Tau={tau:.2f}) 但高表达组织 '{top_tissue}' "
                f"与代谢物积累组织不匹配。"
            )
        else:
            score = 0.20
            explanation = (
                f"TF广泛表达 (Tau={tau:.2f})，无组织特异性调控特征。"
            )

        details = [
            f"组织特异性指数 (Tau): {tau:.3f}",
            f"最高表达组织: {top_tissue} ({sorted_tissues[0][1]:.2f})",
            f"代谢物相关组织: {related_tissues or '未知'}",
            f"组织匹配: {'是' if tissue_match else '否'}",
        ]

        evidence = PriorEvidence(
            strength="strong" if score >= 0.7 else "moderate" if score >= 0.4 else "weak",
            score=score,
            description=f"Tau={tau:.3f}, top_tissue={top_tissue}",
            pmids=[],
        )

        return ReasoningResult(
            score=score,
            explanation=explanation,
            details=details,
            evidence_list=[evidence],
        )
