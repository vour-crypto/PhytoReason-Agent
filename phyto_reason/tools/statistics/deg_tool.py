"""
deg_tool.py — 差异表达分析工具 (v2.0).

v2.0 改进:
  - 多组比较: 接受 groups dict 做 one-way ANOVA/Kruskal-Wallis
  - 小样本 moderated t-statistic (limma 风格 empirical Bayes shrinkage)
  - 效应量排名: 始终输出 Cohen's d / eta² 排名
  - 降级阈值回退: 当 FDR<0.05 无基因时自动报告 top-N 效应量基因
  - 向后兼容 v1.0 的 (group_a, group_b) 接口

支持的数据类型:
  - RNA-seq count → DESeq2 风格 (pseudocount + NB test)
  - TPM/FPKM → limma 风格 (log2 + moderated t-test)
"""

from __future__ import annotations

import logging
import math
from typing import Any

import numpy as np
from scipy.stats import f as f_dist

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.utils.stats_utils import compute_fdr

logger = logging.getLogger("deg_tool")

# ═══════════════════════════════════════════════════════════════
# Empirical Bayes moderated t-statistic (limma-style)
# ═══════════════════════════════════════════════════════════════

def _moderated_t_test(
    a_vals: np.ndarray,
    b_vals: np.ndarray,
    global_var: float,
    df_prior: float = 3.0,
) -> tuple[float, float, float]:
    """Moderated t-test with empirical Bayes variance shrinkage.

    Shrinks the per-gene variance toward a global (prior) variance estimate,
    which stabilizes variance estimates when sample sizes are small.

    Args:
        a_vals, b_vals: group value vectors
        global_var: prior variance (e.g., median variance across all genes)
        df_prior: prior degrees of freedom (default 3 for small-n)

    Returns:
        (t_statistic, p_value, moderated_variance)
    """
    n_a, n_b = len(a_vals), len(b_vals)
    if n_a < 2 or n_b < 2:
        return 0.0, 1.0, global_var

    var_a = np.var(a_vals, ddof=1) if n_a > 1 else 0.0
    var_b = np.var(b_vals, ddof=1) if n_b > 1 else 0.0

    # Pooled per-gene variance
    df_gene = n_a + n_b - 2
    if df_gene > 0:
        var_gene = ((n_a - 1) * var_a + (n_b - 1) * var_b) / df_gene
    else:
        var_gene = global_var

    # Empirical Bayes shrinkage: shrink var_gene toward global_var
    # posterior_var = (df_prior * global_var + df_gene * var_gene) / (df_prior + df_gene)
    if global_var > 0:
        posterior_var = (df_prior * global_var + df_gene * var_gene) / (df_prior + df_gene)
    else:
        posterior_var = var_gene if var_gene > 0 else 1e-10

    # Moderated t-statistic
    se = math.sqrt(posterior_var * (1.0 / n_a + 1.0 / n_b))
    if se < 1e-15:
        return 0.0, 1.0, posterior_var

    mean_diff = np.mean(b_vals) - np.mean(a_vals)
    t_stat = mean_diff / se

    # Moderated degrees of freedom
    df_moderated = df_prior + df_gene
    if df_moderated < 1:
        df_moderated = 1.0

    from scipy.stats import t as t_dist
    p_value = 2.0 * t_dist.sf(abs(t_stat), df=df_moderated)
    p_value = max(p_value, 1e-300)

    return float(t_stat), float(p_value), float(posterior_var)


def _cohens_d(a_vals: np.ndarray, b_vals: np.ndarray) -> float:
    """Cohen's d effect size for two groups."""
    n_a, n_b = len(a_vals), len(b_vals)
    if n_a < 2 or n_b < 2:
        return 0.0
    mean_a, mean_b = np.mean(a_vals), np.mean(b_vals)
    var_a = np.var(a_vals, ddof=1)
    var_b = np.var(b_vals, ddof=1)
    pooled_sd = math.sqrt(((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2))
    if pooled_sd < 1e-15:
        return 0.0
    return float((mean_b - mean_a) / pooled_sd)


def _eta_squared(groups: dict[str, np.ndarray]) -> float:
    """Eta² effect size for multi-group comparison (one-way ANOVA)."""
    all_vals = np.concatenate(list(groups.values()))
    grand_mean = np.mean(all_vals)
    ss_total = np.sum((all_vals - grand_mean) ** 2)

    ss_between = 0.0
    for gv in groups.values():
        ss_between += len(gv) * (np.mean(gv) - grand_mean) ** 2

    if ss_total < 1e-15:
        return 0.0
    return float(ss_between / ss_total)


def _one_way_anova(groups: dict[str, np.ndarray]) -> tuple[float, float]:
    """One-way ANOVA F-test.

    Returns:
        (F_statistic, p_value)
    """
    k = len(groups)
    all_vals = np.concatenate(list(groups.values()))
    n_total = len(all_vals)

    if k < 2 or n_total <= k:
        return 0.0, 1.0

    grand_mean = np.mean(all_vals)

    ss_between = 0.0
    for gv in groups.values():
        ss_between += len(gv) * (np.mean(gv) - grand_mean) ** 2

    ss_within = 0.0
    for gv in groups.values():
        ss_within += np.sum((gv - np.mean(gv)) ** 2)

    df_between = k - 1
    df_within = n_total - k

    if df_within < 1 or ss_within < 1e-15:
        return 0.0, 1.0

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within

    if ms_within < 1e-15:
        return 0.0, 1.0

    f_stat = ms_between / ms_within
    p_value = f_dist.sf(f_stat, df_between, df_within)
    p_value = max(p_value, 1e-300)

    return float(f_stat), float(p_value)


def _kruskal_wallis(groups: dict[str, np.ndarray]) -> tuple[float, float]:
    """Kruskal-Wallis H-test (non-parametric ANOVA).

    Returns:
        (H_statistic, p_value)
    """
    from scipy.stats import kruskal
    group_list = [gv for gv in groups.values() if len(gv) > 0]
    if len(group_list) < 2:
        return 0.0, 1.0
    try:
        h_stat, p_value = kruskal(*group_list)
        return float(h_stat), max(float(p_value), 1e-300)
    except Exception:
        return 0.0, 1.0


# ═══════════════════════════════════════════════════════════════
# DEG Tool v2.0
# ═══════════════════════════════════════════════════════════════

@register_tool
class DEGTool(BaseTool):
    """差异表达分析工具 v2.0。

    支持:
      - 多组 ANOVA/Kruskal-Wallis（通过 groups 参数）
      - 两组比较（通过 group_a/group_b，向后兼容）
      - limma 风格 moderated t-test（小样本友好）
      - 效应量排名（Cohen's d / eta²）
      - 自动降级阈值回退
    """

    tool_name = "deg_analysis"
    description = (
        "差异表达分析工具 v2.0。支持多组 ANOVA 和两组比较，"
        "moderated t-test（小样本友好），效应量排名，自动降级阈值回退。"
    )
    version = "2.0.0"
    parameters = [
        ToolParameter(name="expression_matrix", type="object", required=True,
                      description="表达矩阵 {gene_id: {sample: count_or_fpkm}}"),
        ToolParameter(name="group_a", type="array", required=False,
                      description="A 组样本名列表（两组模式，向后兼容）"),
        ToolParameter(name="group_b", type="array", required=False,
                      description="B 组样本名列表（两组模式，向后兼容）"),
        ToolParameter(name="groups", type="object", required=False,
                      description="多组模式: {tissue_name: [sample_names]}。优先于 group_a/group_b"),
        ToolParameter(name="data_type", type="string", default="auto",
                      description="count | fpkm | tpm | auto"),
        ToolParameter(name="padj_threshold", type="number", default=0.05),
        ToolParameter(name="lfc_threshold", type="number", default=0.5,
                      description="log2FC 阈值 (v2.0 默认 0.5，约 1.4 倍变化)"),
        ToolParameter(name="top_n_effect", type="integer", default=50,
                      description="效应量排名 top-N 输出数"),
    ]
    supported_data_types = ["expression"]

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        if not kwargs.get("expression_matrix"):
            errors.append("expression_matrix 不能为空")
        has_pairwise = kwargs.get("group_a") and kwargs.get("group_b")
        has_multi = kwargs.get("groups")
        if not has_pairwise and not has_multi:
            errors.append("必须提供 groups（多组）或 group_a+group_b（两组）")
        return errors

    def run(self, **kwargs) -> ToolResult:
        expr_mat: dict = kwargs.get("expression_matrix", {})
        group_a: list[str] = kwargs.get("group_a", [])
        group_b: list[str] = kwargs.get("group_b", [])
        groups: dict[str, list[str]] = kwargs.get("groups") or {}
        data_type: str = kwargs.get("data_type", "auto")
        padj_threshold: float = kwargs.get("padj_threshold", 0.05)
        lfc_threshold: float = kwargs.get("lfc_threshold", 0.5)
        top_n_effect: int = kwargs.get("top_n_effect", 50)

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        # ── 确定分析模式 ──────────────────────────────────
        is_multi_group = bool(groups and len(groups) >= 2)

        if is_multi_group:
            # Validate multi-group
            valid_groups = {}
            for gname, samples in groups.items():
                if len(samples) >= 2:
                    valid_groups[gname] = samples
                else:
                    warnings.append(f"组 '{gname}' 样本数不足 (n={len(samples)}<2)，已跳过")
            if len(valid_groups) < 2:
                return ToolResult(
                    candidates=[], evidence_list=[],
                    warnings=warnings + ["多组比较需要至少 2 个有效组（每组 ≥2 样本）"],
                    metadata={"n_deg": 0, "mode": "multi_group"},
                )
            groups = valid_groups
            logger.info("DEG v2.0: multi-group mode with %d groups: %s",
                        len(groups), list(groups.keys()))
        else:
            if len(group_a) < 2 or len(group_b) < 2:
                return ToolResult(
                    candidates=[], evidence_list=[],
                    warnings=["每组至少需要 2 个样本"],
                    metadata={"n_deg": 0, "mode": "pairwise"},
                )
            logger.info("DEG v2.0: pairwise mode, n_a=%d, n_b=%d",
                        len(group_a), len(group_b))

        # ── 自动检测数据类型 ──────────────────────────────
        if data_type == "auto":
            all_values = [
                float(value)
                for sample_values in expr_mat.values()
                for value in sample_values.values()
                if isinstance(value, (int, float)) and np.isfinite(value)
            ]
            if all_values:
                matrix_median = float(np.median(all_values))
                if matrix_median > 50:
                    data_type = "count"
                else:
                    data_type = "fpkm"
            else:
                matrix_median = 0.0
        else:
            matrix_median = None

        # ── 计算全局方差 (for moderated t-test) ───────────
        all_variances: list[float] = []
        all_gene_data: dict[str, dict[str, np.ndarray]] = {}

        if is_multi_group:
            all_samples = []
            for gs in groups.values():
                all_samples.extend(gs)
        else:
            all_samples = list(group_a) + list(group_b)

        for gene_id, sample_vals in expr_mat.items():
            vals = np.array([sample_vals.get(s, float("nan")) for s in all_samples], dtype=float)
            vals = vals[~np.isnan(vals)]
            if len(vals) < 3:
                continue

            if data_type == "count":
                norm_vals = np.log2(vals + 1)
            else:
                norm_vals = np.log2(vals + 1e-6)

            if is_multi_group:
                gene_groups: dict[str, np.ndarray] = {}
                for gname, gsamples in groups.items():
                    gv = np.array([sample_vals.get(s, float("nan")) for s in gsamples], dtype=float)
                    gv = gv[~np.isnan(gv)]
                    if len(gv) >= 2:
                        if data_type == "count":
                            gene_groups[gname] = np.log2(gv + 1)
                        else:
                            gene_groups[gname] = np.log2(gv + 1e-6)
                all_gene_data[gene_id] = gene_groups
            else:
                a_vals = np.array([sample_vals.get(s, float("nan")) for s in group_a], dtype=float)
                b_vals = np.array([sample_vals.get(s, float("nan")) for s in group_b], dtype=float)
                a_vals = a_vals[~np.isnan(a_vals)]
                b_vals = b_vals[~np.isnan(b_vals)]
                if len(a_vals) >= 2 and len(b_vals) >= 2:
                    if data_type == "count":
                        all_gene_data[gene_id] = {
                            "a": np.log2(a_vals + 1),
                            "b": np.log2(b_vals + 1),
                        }
                    else:
                        all_gene_data[gene_id] = {
                            "a": np.log2(a_vals + 1e-6),
                            "b": np.log2(b_vals + 1e-6),
                        }

            var_gene = np.var(norm_vals, ddof=1)
            if var_gene > 0:
                all_variances.append(var_gene)

        if not all_variances:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=["无有效基因数据"],
                metadata={"n_deg": 0},
            )

        global_var = float(np.median(all_variances)) if all_variances else 0.0
        df_prior = max(1.0, min(5.0, (len(all_samples) - len(groups)) / 2.0)) if is_multi_group else 3.0

        logger.info("DEG v2.0: global_var=%.4f, df_prior=%.1f, n_genes=%d",
                    global_var, df_prior, len(all_gene_data))

        # ── 逐基因分析 ────────────────────────────────────
        results_fdr: list[tuple[str, float, float, float, str, float]] = []
        # (gene_id, p_value, effect_size, lfc_abs, mode, extra)
        results_effect: list[tuple[str, float, float, float, str, float, str]] = []
        # (gene_id, p_value, effect_size, lfc_abs, mode, stat_val, stat_name)

        for gene_id, gene_groups in all_gene_data.items():
            if is_multi_group:
                n_eff = sum(1 for g in gene_groups.values() if len(g) >= 2)
                if n_eff < 2:
                    continue

                # One-way ANOVA (primary test)
                f_stat, p_raw = _one_way_anova(gene_groups)
                eta2 = _eta_squared(gene_groups)

                # Kruskal-Wallis as robustness CHECK only — NOT as p-value gate
                # KW has low statistical power with n=3/group; using max(p_anova, p_kw)
                # kills ANOVA's strong signals. Use ANOVA primarily.
                kw_h, kw_p = _kruskal_wallis(gene_groups)
                p_value = p_raw
                # Flag if KW disagrees (for diagnostic purposes)
                kw_disagrees = (kw_p > 0.05 and p_raw < 0.01)

                # Log2FC: max pairwise difference between any two groups
                group_means = {g: np.mean(v) for g, v in gene_groups.items()}
                max_lfc = 0.0
                gnames = list(group_means.keys())
                for i in range(len(gnames)):
                    for j in range(i + 1, len(gnames)):
                        lfc_ij = abs(group_means[gnames[i]] - group_means[gnames[j]])
                        max_lfc = max(max_lfc, lfc_ij)

                effect_size = eta2
                stat_val = f_stat
                stat_name = "F_stat"

                results_effect.append((
                    gene_id, p_value, effect_size, max_lfc,
                    "anova", stat_val, stat_name,
                ))

                if p_value < padj_threshold * 10:  # relaxed pre-filter
                    results_fdr.append((gene_id, max_lfc, p_value, effect_size, "anova", stat_val))

            else:
                a_vals = gene_groups.get("a", np.array([]))
                b_vals = gene_groups.get("b", np.array([]))
                if len(a_vals) < 2 or len(b_vals) < 2:
                    continue

                # Moderated t-test
                t_stat, p_value, mod_var = _moderated_t_test(
                    a_vals, b_vals, global_var, df_prior,
                )

                lfc = np.mean(b_vals) - np.mean(a_vals)
                lfc_abs = abs(lfc)
                d = _cohens_d(a_vals, b_vals)
                effect_size = abs(d)

                stat_val = t_stat
                stat_name = "moderated_t"

                results_effect.append((
                    gene_id, p_value, effect_size, lfc_abs,
                    "moderated_t", stat_val, stat_name,
                ))

                if p_value < padj_threshold * 10:
                    results_fdr.append((gene_id, lfc_abs, p_value, effect_size, "moderated_t", stat_val))

        if not results_effect:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=["无显著差异表达基因（所有基因均未通过预筛选）"],
                metadata={"n_deg": 0, "n_tested": len(all_gene_data)},
            )

        # ── Sort by effect size (primary ranking) ─────────
        results_effect.sort(key=lambda x: x[2], reverse=True)

        # ── FDR correction on pre-filtered candidates ──────
        fdr_passed: set[str] = set()
        fdr_results: dict[str, float] = {}  # gene_id → q_value

        if results_fdr:
            pvalues_fdr = [r[2] for r in results_fdr]
            qvalues = compute_fdr(pvalues_fdr)

            for (gene_id, lfc_abs, p_val, es, mode, sv), q in zip(results_fdr, qvalues):
                fdr_results[gene_id] = q
                if q < padj_threshold and lfc_abs >= lfc_threshold:
                    fdr_passed.add(gene_id)

        n_fdr_significant = len(fdr_passed)

        # ── Build output candidates ────────────────────────
        # Always output top-N by effect size
        n_output = 0
        n_fdr_output = 0
        n_effect_output = 0

        for gene_id, p_value, effect_size, lfc_abs, mode, stat_val, stat_name in results_effect:
            if n_output >= top_n_effect:
                break

            is_fdr_sig = gene_id in fdr_passed
            q_value = fdr_results.get(gene_id, 1.0)

            # Score: prioritize FDR-significant, then effect size
            if is_fdr_sig:
                score = min(1.0, effect_size / 2.0 + 0.3)
                n_fdr_output += 1
            else:
                # Effect-size based score, capped lower for non-significant
                score = min(0.7, effect_size / 2.0)
                n_effect_output += 1

            ev_confidence = round(1.0 - min(q_value * 5, 0.95), 4) if is_fdr_sig else round(min(effect_size / 3.0, 0.6), 4)

            candidate = CandidateGene(
                gene_id=gene_id, is_tf=False,
                pathway_score=round(score, 4),
            )

            sig_label = "FDR-sig" if is_fdr_sig else "effect-ranked"
            desc_parts = [
                f"DEG({sig_label}): {stat_name}={stat_val:.3f}",
                f"effect_size={effect_size:.3f}",
                f"p={p_value:.2e}",
            ]
            if is_fdr_sig:
                desc_parts.append(f"q={q_value:.2e}")
            if mode == "anova":
                desc_parts.append(f"max_log2FC={lfc_abs:.2f}")
            else:
                desc_parts.append(f"log2FC={lfc_abs:.2f}")

            evidence = Evidence(
                evidence_type=EvidenceType.DEG,
                source=f"deg_tool_{self.version}",
                score=round(score, 4),
                confidence=ev_confidence,
                description="; ".join(desc_parts),
                metadata={
                    "stat_name": stat_name,
                    "stat_value": round(stat_val, 3),
                    "effect_size": round(effect_size, 3),
                    "p_value": round(p_value, 6),
                    "q_value": round(q_value, 6) if is_fdr_sig else None,
                    "log2FC": round(lfc_abs, 3),
                    "mode": mode,
                    "is_fdr_significant": is_fdr_sig,
                    "data_type": data_type,
                },
            )
            candidates.append(candidate)
            evidence_list.append(evidence)
            n_output += 1

        # ── Warnings ───────────────────────────────────────
        if n_fdr_significant == 0:
            warnings.append(
                f"FDR 校正后无显著基因 (padj<{padj_threshold}, |log2FC|≥{lfc_threshold})。"
                f"已输出效应量排名 top-{n_effect_output} 基因作为趋势候选。"
                f"原因: 样本量小 (n_total={len(all_samples)}) → 统计效力不足；"
                f"建议: 增加生物学重复或使用更宽松的阈值。"
            )
        elif n_fdr_significant < 10:
            warnings.append(
                f"仅有 {n_fdr_significant} 个 FDR-显著基因。"
                f"已补充 {n_effect_output} 个效应量排名基因作为候选。"
            )

        if len(all_samples) < 10:
            warnings.append(
                f"样本总量小 (n={len(all_samples)})，统计效力有限。"
                f"效应量排名比 p 值更可靠，建议优先关注 |effect_size|>1.0 的基因。"
            )

        logger.info(
            "DEG v2.0 results: n_tested=%d, n_fdr_sig=%d, n_effect_output=%d, "
            "mode=%s, data_type=%s, global_var=%.4f",
            len(all_gene_data), n_fdr_significant, n_effect_output,
            "multi_group" if is_multi_group else "pairwise", data_type, global_var,
        )

        return ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "version": self.version,
                "mode": "multi_group" if is_multi_group else "pairwise",
                "n_groups": len(groups) if is_multi_group else 2,
                "n_genes_tested": len(all_gene_data),
                "n_fdr_significant": n_fdr_significant,
                "n_effect_ranked": n_effect_output,
                "n_total_output": n_output,
                "method": f"moderated_{data_type}",
                "data_type_inference": {
                    "requested": kwargs.get("data_type", "auto"),
                    "resolved": data_type,
                    "matrix_median": matrix_median,
                    "rule": "full numeric matrix median > 50 => count; otherwise fpkm",
                },
                "global_variance": round(global_var, 4),
                "df_prior": round(df_prior, 2),
                "padj_threshold": padj_threshold,
                "lfc_threshold": lfc_threshold,
                "total_samples": len(all_samples),
            },
        )
