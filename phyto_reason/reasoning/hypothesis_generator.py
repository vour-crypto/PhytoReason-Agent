"""
hypothesis_generator.py — 科研假设生成器。
"""

from __future__ import annotations

from phyto_reason.fusion.fusion_result import UnifiedFusionResult


class HypothesisGenerator:
    """科研假设生成器。"""

    @staticmethod
    def generate(fusion_result: UnifiedFusionResult, target_metabolite: str = "", target_pathway: str = "") -> str:
        gene_id = fusion_result.target_gene_id
        calibrated = fusion_result.calibrated_score
        level = fusion_result.confidence_level
        level_str = level.value if hasattr(level, 'value') else str(level)

        supporting = []
        evidence_names = {
            "correlation": "表达相关性", "module_membership": "共表达模块归属",
            "motif": "启动子Motif结合", "pathway": "通路一致性",
            "tf_prior": "TF家族先验", "tissue": "组织特异性", "ortholog": "同源保守性",
        }
        for ev in fusion_result.evidence_list:
            name = evidence_names.get(ev.source, ev.source)
            if ev.normalized_score >= 0.40:
                supporting.append(f"{name} ({ev.normalized_score:.2f})")

        contradictions = fusion_result.contradictions
        contradiction_info = ""
        if contradictions:
            high = [c for c in contradictions if hasattr(c, 'level') and c.level.value == "high"]
            if high:
                contradiction_info = f" [注意: 存在 {len(high)} 项高级别证据矛盾]"

        parts = [f"假设: {gene_id}"]
        if supporting:
            parts.append(f"证据支持: {'; '.join(supporting[:5])}")
        if target_metabolite:
            parts.append(f"推测: {gene_id} 可能参与调控 {target_metabolite} 的生物合成")
        parts.append(f"置信度: {level_str} (融合分={calibrated:.3f}){contradiction_info}")

        if calibrated >= 0.70:
            parts.append("建议: 优先进行湿实验验证 (CRISPR/Cas9 敲除或过表达)")
        elif calibrated >= 0.40:
            parts.append("建议: 可进行 Y1H/Dual-LUC 实验验证结合")
        else:
            parts.append("建议: 需更多数据支持后再进行验证")

        return "\n".join(parts)

    @staticmethod
    def generate_summary(fusion_results: list[UnifiedFusionResult], target_metabolite: str = "") -> str:
        if not fusion_results:
            return "未找到候选TF，无法生成假设。"
        lines = [f"## 生物学假设摘要 (目标: {target_metabolite})", ""]
        ranked = sorted(fusion_results, key=lambda r: r.calibrated_score, reverse=True)
        for i, result in enumerate(ranked[:5], 1):
            lines.append(f"### Top-{i}: {result.target_gene_id}")
            lines.append(f"**融合分**: {result.calibrated_score:.3f} | **级别**: {result.confidence_level}")
            top_evs = sorted(result.evidence_list, key=lambda e: e.normalized_score, reverse=True)[:3]
            for ev in top_evs:
                lines.append(f"  - {ev.source}: {ev.normalized_score:.2f}")
            lines.append("")
        return "\n".join(lines)
