"""
mechanism_classifier.py — 代谢变化机制分类器。

判断当前代谢变化更可能属于哪种机制类型。
这是输入到 TF-reasoning 之前的必要门控，
防止系统在所有情况下默认假设"转录调控"。

可能的机制类型:
  - TRANSCRIPTIONAL: TF 表达变化导致酶基因变化
  - ENZYMATIC: 酶活性变化但不涉及转录
  - TRANSPORT: 转运蛋白介导的重分布
  - DEGRADATION: 降解速率变化
  - STRESS_RESPONSE: 胁迫诱导的一般性重编程
  - DEVELOPMENTAL: 发育程序性变化
  - COMPOSITIONAL: 细胞组成变化（如细胞类型比例）
  - PRECURSOR_LIMITATION: 前体供应不足
  - TECHNICAL_ARTIFACT: 技术假象
"""

from __future__ import annotations

import logging

from phyto_reason.models.mechanistic_hypothesis import MechanismType
from phyto_reason.ontology.metabolite_metadata import resolve_ontology

logger = logging.getLogger("mechanism_classifier")


class MechanismClassificationResult:
    """机制分类结果。"""
    def __init__(self) -> None:
        self.probabilities: dict[str, float] = {}
        self.primary: MechanismType = MechanismType.UNKNOWN
        self.supporting_evidence: list[str] = []
        self.uncertainty: str = "high"

    def to_dict(self) -> dict:
        return {
            "primary": self.primary.value,
            "probabilities": dict(sorted(
                self.probabilities.items(), key=lambda x: x[1], reverse=True,
            )),
            "uncertainty": self.uncertainty,
        }


class MechanismClassifier:
    """基于多维度特征的机制分类器。

    使用规则而非 ML，保证可解释性。
    """

    # Negative-first priority order for mechanism classification.
    # Non-transcriptional mechanisms are evaluated first;
    # TRANSCRIPTIONAL is only considered as the last option.
    # This prevents the system from defaulting to "transcriptional regulation"
    # for every observed metabolite change.
    _DEFAULT_PENALTY = 0.05  # 每个非转录机制的基础置信度衰减

    def classify(
        self,
        expression_data: dict | None = None,
        metabolite_data: dict | None = None,
        target_metabolite: str = "",
        cv: float = 0.0,
        n_samples: int = 0,
        has_promoter: bool = False,
    ) -> MechanismClassificationResult:
        """执行机制分类。

        默认假设: 代谢物变化不是转录调控。
        按优先级排除:
          1. TECHNICAL_ARTIFACT（技术假象）
          2. CELL_COMPOSITION_CHANGE（细胞组成变化）
          3. STRESS_RESPONSE（胁迫响应）
          4. DEVELOPMENTAL（发育程序）
          5. PRECURSOR_LIMITATION / ENZYMATIC / TRANSPORT（中等可能）
          最后: TRANSCRIPTIONAL（特异性调控）
        """
        result = MechanismClassificationResult()
        scores: dict[str, float] = {}
        reasons: list[str] = []

        onto = resolve_ontology(target_metabolite)

        # ── 第一优先级: 排除技术假象 ─────────────────────
        if cv < 0.1 and n_samples >= 3:
            scores[MechanismType.TECHNICAL_ARTIFACT.value] = 0.7
            reasons.append("代谢物在样本间几乎无变异 (CV<0.1)，可能为技术假象")
        elif cv < 0.2 and n_samples >= 3:
            scores[MechanismType.TECHNICAL_ARTIFACT.value] = 0.4
            reasons.append(f"代谢物变异较小 (CV={cv:.2f})，技术假象或生物学信号不足")

        # ── 第二优先级: 细胞组成变化 ────────────────────
        # 无直接检测方法，保持低调检测度
        scores[MechanismType.COMPOSITIONAL.value] = 0.15
        reasons.append("细胞组成变化无法排除（缺乏细胞类型分辨率数据）")

        # ── 第三优先级: 胁迫响应 ────────────────────────
        if onto and onto.stress_associations:
            stress_score = min(len(onto.stress_associations) * 0.15, 0.5)
            scores[MechanismType.STRESS_RESPONSE.value] = stress_score
            reasons.append(f"代谢物 '{target_metabolite}' 与胁迫响应相关: {', '.join(onto.stress_associations[:2])}")
        else:
            scores[MechanismType.STRESS_RESPONSE.value] = 0.1
            reasons.append("代谢物与胁迫无已知关联（环境因素仍需考虑）")

        # ── 第四优先级: 发育程序 ────────────────────────
        if onto and "developmental" in onto.ecological_roles:
            scores[MechanismType.DEVELOPMENTAL.value] = 0.4
            reasons.append("代谢物已知有发育阶段特异性表达")
        else:
            scores[MechanismType.DEVELOPMENTAL.value] = 0.1
            reasons.append("发育阶段数据缺失，无法排除发育调控")

        # ── 第五优先级: 前体限制 / 酶活性 / 转运 ─────────
        if onto and onto.precursor_pathways:
            scores[MechanismType.PRECURSOR_LIMITATION.value] = 0.3
            reasons.append(f"前体途径 ({', '.join(onto.precursor_pathways[:2])}) 可能限速")
        if onto and onto.transport_modes and onto.transport_modes[0] != "unknown":
            scores[MechanismType.TRANSPORT.value] = 0.3
            reasons.append(f"代谢物有已知转运机制 ({onto.transport_modes[0]})，转运重分布可能")
        scores[MechanismType.ENZYMATIC.value] = 0.2
        reasons.append("酶活性调控无法从表达数据推断")

        # ── 最后优先级: 转录调控 ────────────────────────
        # 有非转录机制存在合理证据时，衰减转录调控置信度
        non_transcriptional_max = max(
            scores.get(k, 0) for k in scores
            if k != MechanismType.TRANSCRIPTIONAL.value
        ) if scores else 0.0

        transcriptional_score = 0.15  # 基础分（低）
        if expression_data:
            transcriptional_score += 0.1
            reasons.append("有表达数据，转录调控理论上可检测")
        if has_promoter:
            transcriptional_score += 0.1
            reasons.append("有启动子数据，motif 检测可行")
        if expression_data and has_promoter:
            transcriptional_score += 0.1
            reasons.append("数据和启动子均存在，为转录调控分析提供了条件")

        # 非转录机制有合理证据(≥0.30)时, 衰减转录调控置信度
        if non_transcriptional_max >= 0.30:
            transcriptional_score *= 0.6
            reasons.append("非转录机制存在合理证据，转录调控置信度衰减")
        elif non_transcriptional_max >= 0.15:
            transcriptional_score *= 0.8
            reasons.append("非转录机制不能被排除，转录调控置信度小幅衰减")

        if expression_data:
            scores[MechanismType.TRANSCRIPTIONAL.value] = min(transcriptional_score, 0.5)
        else:
            scores[MechanismType.TRANSCRIPTIONAL.value] = 0.05
            reasons.append("无表达数据，无法评估转录调控")

        # ── 默认: 信息不足 ──────────────────────────────
        if not scores:
            scores[MechanismType.INSUFFICIENT_DATA.value] = 0.9
            reasons.append("信息不足以判断机制类型")

        # ── 归一化 ──────────────────────────────────────
        total = max(sum(scores.values()), 1.0)
        result.probabilities = {k: round(v / total, 3) for k, v in scores.items()}
        result.primary = MechanismType(max(result.probabilities, key=result.probabilities.get))
        result.supporting_evidence = reasons

        # ── 不确定性评估 ────────────────────────────────
        top_prob = max(result.probabilities.values())
        if top_prob < 0.3:
            result.uncertainty = "high"
        elif top_prob < 0.5:
            result.uncertainty = "moderate"
        else:
            result.uncertainty = "low"

        logger.info(
            f"Mechanism classification: primary={result.primary.value}, "
            f"uncertainty={result.uncertainty}, "
            f"top_probs={dict(list(result.probabilities.items())[:4])}"
        )

        return result
