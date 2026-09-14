"""
hypothesis_report.py — v3.0 hypothesis report formatter (pure functions).

仅做格式化, 不做 scientific reasoning。
禁止:
  - 调用统计工具
  - 访问原始矩阵
  - 生成新 hypothesis
"""

from __future__ import annotations

from phyto_reason.models.mechanistic_hypothesis import (
    MechanisticHypothesis, MechanismType, CausalChain, FalsificationTest,
)
from phyto_reason.reasoning.wording_policy import (
    format_confidence_label, build_v3_caveat_section, validate_v3_output_text,
)


# ── 机制类型映射 ──────────────────────────────────────────

MECHANISM_LABELS = {
    MechanismType.TRANSCRIPTIONAL: "transcriptional regulation",
    MechanismType.ENZYMATIC: "post-transcriptional / enzymatic regulation",
    MechanismType.TRANSPORT: "transport-mediated redistribution",
    MechanismType.DEGRADATION: "altered degradation rate",
    MechanismType.STRESS_RESPONSE: "stress-induced metabolic redistribution",
    MechanismType.DEVELOPMENTAL: "developmental program shift",
    MechanismType.COMPOSITIONAL: "cell composition change",
    MechanismType.PRECURSOR_LIMITATION: "precursor limitation",
    MechanismType.TECHNICAL_ARTIFACT: "potential technical artifact",
    MechanismType.INSUFFICIENT_DATA: "insufficient data",
    MechanismType.MIXED: "mixed mechanisms",
    MechanismType.UNKNOWN: "unknown mechanism",
}


def _format_evidence_list(
    evidence_list: list,
    label: str,
    max_items: int = 5,
) -> list[str]:
    """格式化证据列表。"""
    if not evidence_list:
        return [f"  No {label} evidence identified."]
    lines = []
    for ev in evidence_list[:max_items]:
        sev = f" [{ev.contradiction_severity}]" if ev.contradiction_severity else ""
        source = ev.source or ev.evidence_type
        lines.append(f"  - {source}: {ev.description[:120]}{sev}")
    if len(evidence_list) > max_items:
        lines.append(f"  ... and {len(evidence_list) - max_items} more")
    return lines


def _format_causal_chain(chain: CausalChain) -> list[str]:
    """格式化因果链。"""
    lines = []
    if chain.chain_description:
        lines.append(f"  Chain: {chain.chain_description}")
    if chain.known_gaps:
        for gap in chain.known_gaps:
            lines.append(f"  Gap: {gap}")
    if chain.chain_steps:
        for step in chain.chain_steps:
            status = "testable" if step.testable else "untested"
            lines.append(f"  Step {step.position}: {step.subject} {step.predicate} {step.object} [{status}]")
    return lines


def _format_falsification_tests(tests: list[FalsificationTest]) -> list[str]:
    """格式化证伪测试。"""
    if not tests:
        return ["  No falsification tests proposed."]
    lines = []
    for t in tests:
        lines.append(f"  - {t.test_name}: \"{t.prediction[:80]}\"")
        if t.required_data:
            lines.append(f"    Required data: {t.required_data}")
        lines.append(f"    If hypothesis is wrong: {t.contradictory_prediction[:80]}")
    return lines


def format_hypothesis_section(
    h: MechanisticHypothesis,
    rank: int,
    total: int,
) -> list[str]:
    """格式化单个假设的全部内容。"""
    lines = []
    label = MECHANISM_LABELS.get(h.mechanism_type, str(h.mechanism_type.value))
    confidence_label = format_confidence_label(h.uncertainty_level)

    lines.append(f"### Hypothesis {h.id}: {h.title}")
    lines.append(f"  Type: {label} | Confidence: {confidence_label} (uncertainty={h.uncertainty_level})")
    if h.upstream_signal:
        lines.append(f"  Upstream signal: {h.upstream_signal}")
    lines.append("")

    # Causal chain
    if h.causal_chain and h.causal_chain.chain_description:
        lines.append("  Causal chain:")
        lines.extend(_format_causal_chain(h.causal_chain))
        lines.append("")

    # Supporting evidence
    lines.append("  Supporting evidence:")
    lines.extend(_format_evidence_list(h.supporting_evidence, "supporting", max_items=5))
    lines.append("")

    # Contradictory evidence
    lines.append("  Contradictory evidence:")
    lines.extend(_format_evidence_list(h.contradictory_evidence, "contradictory", max_items=5))
    lines.append("")

    # Missing evidence
    lines.append("  Missing evidence:")
    if h.missing_evidence:
        for m in h.missing_evidence[:5]:
            lines.append(f"  - {m}")
    else:
        lines.append("  - None explicitly identified.")
    lines.append("")

    # Why wrong
    if h.why_wrong:
        for line in h.why_wrong.split("\n"):
            lines.append(f"  {line}")
        lines.append("")

    # Falsification tests
    lines.append("  Falsification tests:")
    lines.extend(_format_falsification_tests(h.falsification_tests))
    lines.append("")

    # Validation experiments
    if h.validation_experiments:
        lines.append("  Validation experiments:")
        for v in h.validation_experiments:
            lines.append(f"  - {v}")
        lines.append("")

    # Next best evidence
    if h.next_best_evidence:
        lines.append("  Next best evidence to acquire:")
        for n in h.next_best_evidence:
            lines.append(f"  - {n}")
        lines.append("")

    return lines


def format_competition_context(
    ranked: list[MechanisticHypothesis],
    idx: int,
) -> str:
    """解释为什么一个假设排名高/低。"""
    if idx >= len(ranked):
        return ""
    h = ranked[idx]
    n_support = len(h.supporting_evidence)
    n_contra = len(h.contradictory_evidence)
    n_missing = len(h.missing_evidence)
    total = n_support + n_contra + n_missing
    ratio = n_support / max(total, 1)

    if idx == 0:
        if ratio >= 0.5:
            return (f"Top-ranked: {n_support} supporting vs {n_contra} contradictory "
                    f"evidence sources (ratio={ratio:.2f})")
        return (f"Top-ranked despite low ratio ({ratio:.2f}): "
                f"other hypotheses have even less evidence")
    if ratio < 0.3:
        return f"Weak support ratio ({ratio:.2f}): {n_contra} contradictory items reduce confidence"
    return (f"Rank {idx + 1}/{len(ranked)}: {n_support} supporting, "
            f"{n_contra} contradictory, {n_missing} missing")


def build_full_report(
    competing_hypotheses: list[MechanisticHypothesis] | None = None,
    target_metabolite: str = "",
    pathway_name: str = "",
    n_samples: int = 0,
    n_excluded: int = 0,
    mechanism_probabilities: dict | None = None,
) -> list[str]:
    """构建完整 v3.0 hypothesis report。

    返回行列表，调用者负责 join/写入。
    """
    lines: list[str] = []
    hypotheses = competing_hypotheses or []

    # A. Executive Summary
    lines.append("# Executive Summary")
    if not hypotheses:
        lines.append("No reliable mechanistic hypotheses could be formed from available data.")
        if pathway_name:
            lines.append(f"Target pathway '{pathway_name}' may not be active under current conditions.")
        lines.append("Recommendation: confirm pathway activity with targeted metabolomics or perturbation data.")
        lines.append("")
        return lines

    top = hypotheses[0]
    top_label = MECHANISM_LABELS.get(top.mechanism_type, str(top.mechanism_type.value))
    lines.append(f"Target metabolite: {target_metabolite or 'unspecified'}")
    lines.append(f"Pathway: {pathway_name or 'unknown'}")
    lines.append(f"Top hypothesis: {top.title} ({top_label}, confidence={format_confidence_label(top.uncertainty_level)})")
    lines.append(f"Competing hypotheses: {len(hypotheses)}")
    lines.append("")

    # B. Mechanism Classification Summary
    lines.append("## Mechanism Classification Summary")
    if mechanism_probabilities:
        for k, v in sorted(mechanism_probabilities.items(), key=lambda x: -x[1]):
            lines.append(f"  - {k}: {v:.1%}")
    lines.append("")

    # C. Competing Hypotheses
    lines.append("## Competing Hypotheses")
    lines.append("")

    for i, h in enumerate(hypotheses):
        context = format_competition_context(hypotheses, i)
        section = format_hypothesis_section(h, i + 1, len(hypotheses))
        lines.extend(section)
        lines.append(f"  Competition context: {context}")
        lines.append("---")
        lines.append("")

    # D. Why the Top Hypothesis May Still Be Wrong
    lines.append("## Why the Top Hypothesis May Still Be Wrong")
    if top.why_wrong:
        for line in top.why_wrong.split("\n"):
            lines.append(f"  {line}")
    else:
        lines.append("  No specific contradictory evidence identified — but absence of evidence is not evidence of absence.")
    lines.append("")

    # E. Missing Critical Evidence
    lines.append("## Missing Critical Evidence")
    all_missing = set()
    for h in hypotheses:
        all_missing.update(h.missing_evidence)
    if all_missing:
        for m in all_missing:
            lines.append(f"  - {m}")
    else:
        lines.append("  No specific missing evidence cataloged.")
    lines.append("")

    # F. Recommended Validation Experiments
    lines.append("## Recommended Validation Experiments")
    all_experiments = set()
    for h in hypotheses:
        all_experiments.update(h.validation_experiments)
    if all_experiments:
        for v in all_experiments:
            lines.append(f"  - {v}")
    else:
        lines.append("  - Y1H/Dual-LUC: validate TF-promoter binding for the top candidate regulator")
        lines.append("  - CRISPR/Cas9 or RNAi knockdown: observe metabolite phenotype")
        lines.append("  - Temporal expression analysis: confirm ordering of regulatory events")
    lines.append("")

    # G. Final Uncertainty Statement
    lines.append("## Final Uncertainty Statement")
    caveats = build_v3_caveat_section(
        hypotheses=hypotheses,
        n_samples=n_samples,
        n_excluded=n_excluded,
    )
    for caveat in caveats:
        lines.append(f"  - {caveat}")
    lines.append("")

    return lines
