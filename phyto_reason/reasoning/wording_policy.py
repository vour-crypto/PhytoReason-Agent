"""
wording_policy.py — 科研措辞控制策略。

防止系统在证据不足时输出 "regulates"、"controls" 等因果性措辞。

原则:
  - 仅 experimental validation → 允许 causal verbs
  - 仅 correlation → 只能 associative verbs
  - multi-evidence → candidate regulator
"""

from __future__ import annotations

from phyto_reason.fusion.fusion_result import UnifiedFusionResult


# ── 动词强度等级 ──────────────────────────────────────────

CAUSAL_VERBS = [
    "regulates", "controls", "activates", "represses",
    "is a master regulator of", "directly binds",
]

ASSOCIATIVE_VERBS = [
    "is associated with", "correlates with",
    "is co-expressed with", "shows correlation with",
]

CANDIDATE_VERBS = [
    "is a candidate regulator of",
    "may be involved in",
    "potentially linked to",
]

UNCERTAINTY_PREFIXES = [
    "Based on expression correlation data,",
    "Current evidence suggests that",
    "Limited to available data,",
]


def select_verb(
    fusion_result: UnifiedFusionResult | None = None,
    has_experimental_validation: bool = False,
    n_independent_evidences: int = 0,
    max_correlation: float = 0.0,
    has_motif: bool = False,
    has_literature: bool = False,
) -> str:
    """根据证据强度选择合适的动词。

    Args:
        fusion_result: 融合结果
        has_experimental_validation: 是否有湿实验验证
        n_independent_evidences: 独立证据源数量
        max_correlation: 最大相关系数
        has_motif: 是否有 motif 证据
        has_literature: 是否有文献支持

    Returns:
        合适的动词短语
    """
    if fusion_result is not None:
        n_independent_evidences = fusion_result.hierarchy.n_independent_sources
        max_correlation = max(
            (e.normalized_score for e in fusion_result.evidence_list
             if e.source == "correlation"),
            default=0.0,
        )
        has_motif = any(e.source in ("motif", "motif_real") and e.normalized_score >= 0.30
                        for e in fusion_result.evidence_list)
        has_literature = any(e.source == "literature" and e.normalized_score >= 0.30
                             for e in fusion_result.evidence_list)

    if has_experimental_validation:
        return CAUSAL_VERBS[0]

    if n_independent_evidences >= 3 and has_motif and has_literature:
        return CANDIDATE_VERBS[0]

    if n_independent_evidences >= 2 and max_correlation >= 0.70:
        return ASSOCIATIVE_VERBS[0]

    return ASSOCIATIVE_VERBS[1]


def format_hypothesis_statement(
    tf_id: str,
    target_enzymes: list[str],
    pathway_name: str,
    metabolite: str,
    fusion_result: UnifiedFusionResult | None = None,
) -> str:
    """生成符合措辞规范的假设语句。

    输出格式:
      "Current evidence suggests that [TF] is associated with [enzyme]"
      而非 "[TF] regulates [enzyme]"
    """
    verb = select_verb(fusion_result)
    prefix = ""

    if verb in CAUSAL_VERBS:
        prefix = ""
    elif verb in ASSOCIATIVE_VERBS:
        prefix = "Current evidence suggests that "
    else:
        prefix = "Limited to available data, "

    enzyme_str = ", ".join(target_enzymes[:3]) if target_enzymes else pathway_name

    if prefix:
        return f"{prefix}{tf_id} {verb} {enzyme_str} in the context of {metabolite} biosynthesis"
    return f"{tf_id} {verb} {enzyme_str} in {metabolite} biosynthesis"


# ── v3.0 置信度标签映射 ──────────────────────────────────
LEVEL_LABELS = {
    "high": "insufficient evidence",
    "moderate": "weakly supported",
    "low": "moderately plausible",
    "gold": "plausible",
    "silver": "weakly supported",
    "weak": "insufficient evidence",
    "insufficient": "insufficient evidence",
}


def format_confidence_label(uncertainty_level: str) -> str:
    """v3.0: 将 uncertainty_level 映射为证据支持的措辞。
    Gold/Silver → plausible / weakly supported
    """
    return LEVEL_LABELS.get(uncertainty_level.lower(), "insufficient evidence")


def validate_v3_output_text(text: str) -> list[str]:
    """v3.0 措辞检查 — 禁止 causal certainty 和 overclaim。"""
    violations: list[str] = []
    forbidden = [
        "identified regulator", "confirmed", "proved",
        "master regulator", "strongly proves", "demonstrates causality",
    ]
    for phrase in forbidden:
        if phrase in text.lower():
            violations.append(f"禁止措辞 '{phrase}' — 非实验验证结论不得使用")
    for verb in CAUSAL_VERBS:
        if verb in text.lower():
            violations.append(f"因果动词 '{verb}' 不应用于未经实验验证的结论")
    return violations


def build_v3_caveat_section(
    hypotheses: list | None = None,
    n_samples: int = 0,
    n_excluded: int = 0,
) -> list[str]:
    """v3.0 局限性声明 — 基于 MechanisticHypothesis 而非 UnifiedFusionResult。"""
    caveats: list[str] = []
    caveats.append("基于静态表达数据的相关性推断，无法确定因果关系")
    if n_samples < 8:
        caveats.append(f"小样本量 (n={n_samples})，统计效力有限")
    if hypotheses:
        n_ev = max(len(h.supporting_evidence) for h in hypotheses) if hypotheses else 0
        if n_ev < 3:
            caveats.append(f"所有假设的证据源均少于 3 个独立来源，结论可靠性受限")
    if n_excluded > 0:
        caveats.append(f"{n_excluded} 个候选因证据不足或矛盾被排除，剩余候选的可靠性需独立验证")
    return caveats


def validate_output_text(text: str) -> list[str]:
    """检查输出文本是否包含不允许的措辞。

    Returns:
        违反规则的措辞列表（空列表 = 通过）
    """
    violations: list[str] = []
    for verb in CAUSAL_VERBS:
        if verb in text.lower():
            violations.append(f"高因果强度动词 '{verb}' 不应用于未经实验验证的结论")
    return violations


def build_caveat_section(
    fusion_result: UnifiedFusionResult | None = None,
    n_samples: int = 0,
    has_deg: bool = False,
    has_perturbation: bool = False,
) -> list[str]:
    """构建结论的局限性和不确定性声明。"""
    caveats: list[str] = []

    if not has_perturbation:
        caveats.append("基于静态表达数据的相关性推断，无法确定因果关系")

    if n_samples < 8:
        caveats.append(f"小样本量 (n={n_samples})，统计效力有限")

    if fusion_result:
        n_ev = fusion_result.hierarchy.n_independent_sources
        if n_ev < 3:
            caveats.append(f"证据源单一 (仅 {n_ev} 个独立来源)，结论可靠性受限")

    if not has_deg:
        caveats.append("未经过差异表达分析和多重假设检验校正")

    return caveats
