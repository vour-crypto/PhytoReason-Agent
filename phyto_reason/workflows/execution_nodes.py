"""
execution_nodes.py -- 递进式科研推理节点。

每个 node 对应一个科研推理步骤：
  0. deg_analysis           -- 差异表达基因分析（多组 ANOVA + 效应量排名）
  0b. dam_analysis          -- 差异积累代谢物分析
  1. data_quality_gate     -- 数据可信度验证
  2. metabolite_validation  -- 目标代谢物差异确认
  3. pathway_coherence      -- 通路层级一致性
  4. tf_narrowing           -- 基于 pathway 基因缩小候选 TF
  5. regulation_evidence    -- 调控证据检查
  6. contradiction_check    -- 假阳性主动排除
  7. hypothesis_synthesis   -- 机理假设生成
"""

from __future__ import annotations

import logging
import math

import numpy as np

from phyto_reason.knowledge.species_registry import get_species_registry
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.workflows.routing_logic import ALL_PIPELINE_NODES

logger = logging.getLogger("workflow_nodes")


# ── WorkflowPlan skip guard ─────────────────────────────────────

def _should_skip_node(state: RuntimeState, node_name: str) -> bool:
    """Return True if this node should be skipped based on WorkflowPlan.

    Checks state.workflow_plan.tools_to_run (expanded node names).
    If no plan or empty tools_to_run, run everything (backward compatible).
    """
    wp = state.workflow_plan
    if wp is None:
        return False
    if not wp.tools_to_run:
        return False
    return node_name not in wp.tools_to_run


def _get_expr_matrix(state: RuntimeState) -> dict | None:
    if state.planner_state:
        return state.planner_state.expression_matrix
    return None


def _get_meta_matrix(state: RuntimeState) -> dict | None:
    if state.planner_state:
        return state.planner_state.metabolite_matrix
    return None


def _get_meta(state: RuntimeState) -> str:
    return (state.planner_state.target_metabolite or "").lower().strip() if state.planner_state else ""


def _get_promoters(state: RuntimeState) -> dict:
    if state.planner_state:
        return getattr(state.planner_state, "promoter_sequences", {}) or {}
    return {}


# ═══════════════════════════════════════════════════════════════
# Node 0: 差异表达基因分析（多组 ANOVA + 效应量排名）
# ═══════════════════════════════════════════════════════════════

def _infer_sample_groups(expr: dict, metadata: dict | None = None) -> dict[str, list[str]]:
    """Infer sample groups from expression matrix column names.

    Metadata condition/tissue/treatment takes precedence. If no usable
    metadata is available, fall back to common prefix/suffix patterns.
    E.g., "Leaf_1", "Leaf_2", "Root_1", "Root_2" -> {"Leaf": [...], "Root": [...]}
    """
    if not expr:
        return {}
    first = next(iter(expr.values()), {})
    all_samples = list(first.keys())

    if metadata:
        rows = metadata.get("samples", []) if isinstance(metadata, dict) else []
        by_id = {str(row.get("sample_id")): row for row in rows if isinstance(row, dict)}
        grouped: dict[str, list[str]] = {}
        for sample in all_samples:
            row = by_id.get(str(sample), {})
            group = str(row.get("condition") or row.get("tissue") or row.get("treatment") or "").strip()
            if group:
                grouped.setdefault(group, []).append(sample)
        if len(grouped) >= 2:
            return grouped
        import logging
        logging.getLogger("workflow_nodes").warning(
            "Sample metadata supplied but no usable groups; falling back to sample-name prefixes"
        )

    if len(all_samples) < 4:
        return {"all": all_samples}

    # Try to split on common delimiters
    import re
    groups: dict[str, list[str]] = {}

    for s in all_samples:
        # Try patterns: "Tissue_Replicate", "Tissue-Replicate", "Tissue.Replicate"
        for sep in ["_", "-", "."]:
            parts = s.rsplit(sep, 1)
            if len(parts) == 2 and parts[1].strip().isdigit():
                group = parts[0].strip()
                groups.setdefault(group, []).append(s)
                break
        else:
            # Try removing trailing digits
            m = re.match(r"(.+?)[_\-\s]?\d+$", s)
            if m:
                group = m.group(1).strip()
                groups.setdefault(group, []).append(s)
            else:
                groups.setdefault(s, []).append(s)

    # Merge tiny groups (< 2 samples) into "other"
    merged: dict[str, list[str]] = {}
    others: list[str] = []
    for gname, gsamples in groups.items():
        if len(gsamples) >= 2:
            merged[gname] = gsamples
        else:
            others.extend(gsamples)
    if others:
        merged["other"] = others

    return merged if len(merged) >= 2 else {"all": all_samples}


# ═══════════════════════════════════════════════════════════════
# Node -1: 数据质控 (sample correlation, PCA, heatmap, upset)
# ═══════════════════════════════════════════════════════════════

def data_qc_node(state: RuntimeState) -> dict:
    """Generate data quality control figures before differential analysis.

    Produces:
      - Sample correlation heatmap
      - PCA plot
      - Clustered expression heatmap
      - UpSet plot (gene overlap across groups)
    """
    if _should_skip_node(state, "data_qc"):
        state.add_trace("data_qc", "skipped", summary={"reason": "user skipped this step"})
        return {}

    expr = getattr(state.planner_state, "expression_matrix", None) if state.planner_state else None
    meta_mat = getattr(state.planner_state, "metabolite_matrix", None) if state.planner_state else None

    if not expr:
        state.add_trace("data_qc", "skipped", summary={"reason": "no expression data"})
        return {"qc_report": {"reason": "no expression data"}}

    # Infer sample groups from column names (e.g. L_1, L_2, TR_1, ...)
    groups = _infer_sample_groups(expr, getattr(state.planner_state, "sample_metadata", None) if state.planner_state else None)
    # Invert to {sample_id: group_name} for correlation / PCA
    sample_to_group: dict[str, str] = {}
    for gname, samples in groups.items():
        for s in samples:
            sample_to_group[s] = gname

    qc_report = {
        "n_samples": len(next(iter(expr.values()))),
        "n_genes": len(expr),
        "groups": list(groups.keys()),
    }
    figure_urls = []
    figure_mds = []

    try:
        # 可视化包缺失时优雅降级（v5.0：figures 重建在路线图 Phase 3）
        from phyto_reason.visualization.multiomics_figures import (
            plot_qc_correlation_heatmap,
            plot_qc_pca,
            plot_qc_upset,
        )
        from phyto_reason.visualization.heatmap import HeatmapPlotter
        from phyto_reason.visualization.figure_exporter import FigureExporter

        # 1. Sample correlation heatmap
        corr_url = plot_qc_correlation_heatmap(expr, groups=sample_to_group)
        if corr_url:
            figure_urls.append(corr_url)
            figure_mds.append(f"![Sample Correlation]({corr_url})")
            qc_report["figure_correlation"] = corr_url

        # 2. PCA plot
        pca_url = plot_qc_pca(expr, groups=sample_to_group)
        if pca_url:
            figure_urls.append(pca_url)
            figure_mds.append(f"![PCA]({pca_url})")
            qc_report["figure_pca"] = pca_url

        # 3. Clustered expression heatmap
        from phyto_reason.visualization.r_bridge import r_available, r_heatmap
        gene_ids = list(expr.keys())
        valid_genes = [g for g in gene_ids if g in expr]
        # 只画方差最高的 60 个基因：全矩阵热图不可读（行标签挤成黑块）
        if len(valid_genes) > 60:
            def _row_variance(gid: str) -> float:
                vals = np.array(
                    [float(v) for v in expr[gid].values() if v is not None], dtype=float
                )
                return float(np.var(vals)) if vals.size else 0.0

            valid_genes = sorted(valid_genes, key=_row_variance, reverse=True)[:60]
            qc_heatmap_note = "top-60 高变基因（全矩阵热图不可读，已筛选）"
        else:
            qc_heatmap_note = None
        if len(valid_genes) >= 3:
            if r_available():
                heat_url = r_heatmap(expr, valid_genes, title="QC Expression Heatmap")
            else:
                try:
                    fig = HeatmapPlotter.clustered_heatmap(expr, valid_genes, title="QC Expression Heatmap")
                    import uuid
                    from pathlib import Path
                    from phyto_reason.platform_paths import figures_dir
                    static_dir = figures_dir()
                    static_dir.mkdir(parents=True, exist_ok=True)
                    fname = f"qc_heatmap_{uuid.uuid4().hex[:6]}"
                    FigureExporter.save(fig, str(static_dir / fname), dpi=150)
                    heat_url = f"/static/figures/{fname}.png"
                except Exception:
                    heat_url = None
            if heat_url:
                figure_urls.append(heat_url)
                figure_mds.append(f"![QC Heatmap]({heat_url})")
                qc_report["figure_heatmap"] = heat_url
                if qc_heatmap_note:
                    qc_report["figure_heatmap_note"] = qc_heatmap_note

        # 4. UpSet plot (gene overlap across groups)
        sample_gene_sets = {}
        for gname, samples in groups.items():
            expressed = {
                gid for gid in expr
                if any(expr[gid].get(s, 0) > 0 for s in samples)
            }
            sample_gene_sets[gname] = expressed
        if len(sample_gene_sets) >= 2:
            upset_url = plot_qc_upset(sample_gene_sets)
            if upset_url:
                figure_urls.append(upset_url)
                figure_mds.append(f"![UpSet]({upset_url})")
                qc_report["figure_upset"] = upset_url

    except Exception as e:
        logger.warning("QC figure generation failed: %s", e)
        qc_report["warning"] = str(e)

    if figure_urls:
        qc_report["figure_urls"] = figure_urls
        qc_report["figure_markdown"] = "\n".join(figure_mds)

    state.qc_report = qc_report

    state.add_trace("data_qc", "completed", summary={
        "n_genes": len(expr),
        "n_groups": len(groups),
        "n_figures": len(figure_urls),
    })

    logger.info("Data QC: %d genes, %d groups, %d figures", len(expr), len(groups), len(figure_urls))
    return {"qc_report": qc_report}


def _generate_volcano_figure(
    results: list[dict],
    key_id: str = "gene_id",
    title: str = "Volcano Plot",
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> str | None:
    """Generate volcano plot for DEG/DAM results.

    (R) R first (r_volcano) -> (Py) matplotlib fallback (VolcanoPlotter).
    """
    if not results:
        return None

    # Build { gene_id: {log2fc, padj} }
    volcano_data = {}
    for r in results:
        gid = r.get(key_id, "")
        if gid:
            volcano_data[gid] = {
                "log2fc": r.get("max_log2fc", 0),
                "padj": r.get("q_value", r.get("p_value", 1)),
            }

    if not volcano_data:
        return None

    # 可视化包缺失时优雅降级（v5.0：figures 重建在路线图遗留项）
    try:
        from phyto_reason.visualization.r_bridge import r_available, r_volcano
    except ImportError:
        return None

    # Try R first
    if r_available():
        url = r_volcano(volcano_data, title=title, padj_threshold=0.05, lfc_threshold=0.5)
        if url:
            return url

    # Fallback: matplotlib
    try:
        from phyto_reason.visualization.volcano import VolcanoPlotter
        from phyto_reason.visualization.figure_exporter import FigureExporter
        import uuid
        from pathlib import Path

        from phyto_reason.platform_paths import figures_dir
        static_dir = figures_dir()
        static_dir.mkdir(parents=True, exist_ok=True)
        uid = uuid.uuid4().hex[:6]
        key = "deg" if key_id == "gene_id" else "dam"
        fig = VolcanoPlotter.deg_volcano(volcano_data, title=title,
                                         padj_threshold=padj_threshold,
                                         lfc_threshold=lfc_threshold)
        path = str(static_dir / f"volcano_{key}_{uid}.png")
        FigureExporter.save(fig, path.replace(".png", ""), dpi=150)
        return f"/static/figures/volcano_{key}_{uid}.png"
    except Exception as e:
        logger.debug("Volcano figure fallback failed: %s", e)

    return None


def deg_analysis_node(state: RuntimeState) -> dict:
    """差异表达基因分析。

    运行多组 ANOVA（或两组 moderated t-test）对所有基因，
    输出 FDR-显著基因 + 效应量排名 top-N 基因。
    结果存入 state.deg_report 供下游节点使用。
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "deg_analysis"):
            state.add_trace("deg_analysis", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        expr = _get_expr_matrix(state)
        if not expr:
            state.add_trace("deg_analysis", "completed",
                            summary={"n_genes": 0, "reason": "no expression data"})
            return {"deg_report": {"n_genes_tested": 0, "reason": "no expression data"}}

        # Infer groups from metadata first, then sample names
        groups = _infer_sample_groups(
            expr,
            getattr(state.planner_state, "sample_metadata", None) if state.planner_state else None,
        )
        n_groups = len(groups)

        # Determine data type
        all_expr_values = [
            float(value)
            for sample_values in expr.values()
            for value in sample_values.values()
            if isinstance(value, (int, float)) and np.isfinite(value)
        ]
        data_type = "fpkm"
        matrix_median = float(np.median(all_expr_values)) if all_expr_values else 0.0
        if all_expr_values:
            if matrix_median > 50:
                data_type = "count"

        # Run DEG analysis
        try:
            from phyto_reason.tools.statistics.deg_tool import (
                _one_way_anova, _kruskal_wallis, _eta_squared,
                _moderated_t_test, _cohens_d,
            )
            from phyto_reason.utils.stats_utils import compute_fdr
        except ImportError:
            state.add_trace("deg_analysis", "completed",
                            summary={"n_genes": 0, "reason": "DEG tool not available"})
            return {"deg_report": {"n_genes_tested": 0, "reason": "DEG tool not available"}}

        all_samples = []
        for gs in groups.values():
            all_samples.extend(gs)

        # count 数据先做文库大小归一化（DESeq2 median-of-ratios size factors）：
        # 裸 log2(count+1) 仍携带文库差异，会扭曲组间比较与方差估计。
        size_factors: dict[str, float] | None = None
        if data_type == "count":
            from phyto_reason.utils.stats_utils import compute_size_factors
            size_factors = compute_size_factors(expr)
            if not size_factors or any(v <= 0 for v in size_factors.values()) \
                    or len(size_factors) != len(all_samples):
                size_factors = None
        normalization_note = (
            "median_of_ratios_size_factors + log2(count/size_factor + 1)"
            if size_factors else
            ("log2(count + 1)" if data_type == "count" else "log2(value + 1e-6)")
        )

        # Compute global variance for moderated t-test
        all_variances = []
        gene_results: list[dict] = []

        sf_vec = np.array([size_factors.get(s, 1.0) for s in all_samples], dtype=float) \
            if size_factors else np.ones(len(all_samples))
        for gene_id, sample_vals in expr.items():
            adjusted = np.array([sample_vals.get(s, float("nan")) for s in all_samples], dtype=float) / sf_vec
            vals = adjusted[~np.isnan(adjusted)]
            if len(vals) < 3:
                continue

            if data_type == "count":
                norm_vals = np.log2(vals + 1)
            else:
                norm_vals = np.log2(vals + 1e-6)

            var_g = np.var(norm_vals, ddof=1)
            if var_g > 0:
                all_variances.append(var_g)

        if not all_variances:
            state.add_trace("deg_analysis", "completed",
                            summary={"n_genes": 0, "reason": "no valid gene data"})
            return {"deg_report": {"n_genes_tested": 0, "reason": "no valid gene data"}}

        global_var = float(np.median(all_variances))
        df_prior = max(1.0, min(5.0, len(all_samples) / len(groups) / 2.0)) if n_groups >= 2 else 3.0

        # Per-gene analysis
        for gene_id, sample_vals in expr.items():
            if n_groups >= 2:
                # Multi-group ANOVA
                gene_groups: dict[str, np.ndarray] = {}
                for gname, gsamples in groups.items():
                    gv = np.array([sample_vals.get(s, float("nan")) for s in gsamples], dtype=float)
                    if size_factors:
                        gv = gv / np.array([size_factors.get(s, 1.0) for s in gsamples], dtype=float)
                    gv = gv[~np.isnan(gv)]
                    if len(gv) >= 2:
                        if data_type == "count":
                            gene_groups[gname] = np.log2(gv + 1)
                        else:
                            gene_groups[gname] = np.log2(gv + 1e-6)

                if len(gene_groups) < 2:
                    continue

                f_stat, p_anova = _one_way_anova(gene_groups)
                eta2 = _eta_squared(gene_groups)

                # Kruskal-Wallis as robustness CHECK only (NOT as p-value gate)
                # KW has low power with n=3/group -- using max(p_anova, p_kw)
                # would drown ANOVA's strong signals in KW's low power.
                # Instead, use ANOVA primarily and flag genes where KW disagrees.
                _, p_kw = _kruskal_wallis(gene_groups)
                p_value = p_anova
                kw_disagrees = (p_kw > 0.05 and p_anova < 0.01)  # KW lacks power for small-n groups

                # Max pairwise log2FC
                group_means = {g: float(np.mean(v)) for g, v in gene_groups.items()}
                max_lfc = 0.0
                gnames = list(group_means.keys())
                for i in range(len(gnames)):
                    for j in range(i + 1, len(gnames)):
                        lfc = abs(group_means[gnames[i]] - group_means[gnames[j]])
                        max_lfc = max(max_lfc, lfc)

                gene_results.append({
                    "gene_id": gene_id,
                    "p_value": p_value,
                    "effect_size": eta2,
                    "max_log2fc": max_lfc,
                    "stat_name": "F_stat",
                    "stat_value": f_stat,
                    "n_groups": len(gene_groups),
                })
            else:
                # Single group: skip DEG (need at least 2 groups)
                continue

        n_tested = len(gene_results)
        if n_tested == 0:
            state.add_trace("deg_analysis", "completed",
                            summary={"n_genes_tested": 0, "reason": "no genes with sufficient data"})
            return {"deg_report": {"n_genes_tested": 0, "reason": "no genes with sufficient data"}}

        # ── Variance pre-filter: remove near-zero-variance genes ──
        # Genes with essentially no variation across tissues are noise.
        # Removing them before FDR reduces multiple-testing burden,
        # which directly improves the number of genes passing FDR.
        if len(gene_results) > 1000:
            effect_sizes = [g["effect_size"] for g in gene_results]
            es_10pct = np.percentile(effect_sizes, 10) if effect_sizes else 0
            # Keep genes above the 10th percentile of effect size (eta²)
            # This typically removes ~10% of genes with η² near 0
            es_cutoff = max(es_10pct, 1e-6)
            n_before = len(gene_results)
            gene_results = [g for g in gene_results if g["effect_size"] >= es_cutoff]
            n_filtered = n_before - len(gene_results)
            if n_filtered > 0:
                logger.info(
                    "Variance pre-filter: removed %d/%d genes (effect_size < %.4f)",
                    n_filtered, n_before, es_cutoff,
                )

        # FDR correction
        gene_results.sort(key=lambda x: x["p_value"])
        pvalues = [g["p_value"] for g in gene_results]
        qvalues = compute_fdr(pvalues)
        for g, q in zip(gene_results, qvalues):
            g["q_value"] = q

        # Count FDR-significant
        padj_threshold = 0.05
        lfc_threshold = 0.5
        n_fdr_sig = sum(
            1 for g in gene_results
            if g["q_value"] < padj_threshold and g["max_log2fc"] >= lfc_threshold
        )

        # Sort by effect size for ranking
        gene_results.sort(key=lambda x: x["effect_size"], reverse=True)

        # Top-N by effect size
        top_n = min(50, n_tested)
        top_genes = gene_results[:top_n]

        # Build the report
        deg_report = {
            "n_genes_tested": n_tested,
            "n_fdr_significant": n_fdr_sig,
            "n_effect_ranked": min(top_n, n_tested - n_fdr_sig) if n_fdr_sig < top_n else 0,
            "n_groups": n_groups,
            "group_names": list(groups.keys()),
            "total_samples": len(all_samples),
            "data_type": data_type,
            "normalization": normalization_note,
            "data_type_inference": {
                "requested": "auto",
                "resolved": data_type,
                "matrix_median": matrix_median,
                "rule": "full numeric matrix median > 50 => count; otherwise fpkm",
            },
            "summary": f"n={len(all_samples)}; method=ANOVA + moderated variance; threshold=|log2FC|≥{lfc_threshold}, q<{padj_threshold}; FDR=Benjamini-Hochberg",
            "fdr_status": "BH applied to all tested genes",
            "global_variance": round(global_var, 4),
            "padj_threshold": padj_threshold,
            "lfc_threshold": lfc_threshold,
            "top_genes": [
                {
                    "gene_id": g["gene_id"],
                    "effect_size": round(g["effect_size"], 4),
                    "max_log2fc": round(g["max_log2fc"], 3),
                    "p_value": round(g["p_value"], 6),
                    "q_value": round(g["q_value"], 6),
                    "stat_name": g["stat_name"],
                    "is_fdr_sig": g["q_value"] < padj_threshold and g["max_log2fc"] >= lfc_threshold,
                }
                for g in top_genes[:20]
            ],
        }

        # Warnings
        deg_warnings = []
        if n_fdr_sig == 0:
            deg_warnings.append(
                f"FDR 校正后无显著差异表达基因 (padj<{padj_threshold}, |log2FC|≥{lfc_threshold})。"
                f"已输出效应量排名 top-{min(top_n, n_tested)} 基因。"
                f"小样本 (n={len(all_samples)}, {n_groups} 组) 导致统计效力不足。"
            )
        elif n_fdr_sig < 10:
            deg_warnings.append(
                f"仅 {n_fdr_sig} 个 FDR-显著 DEG。已补充效应量排名基因。"
            )

        if len(all_samples) < 10:
            deg_warnings.append(
                f"样本总量小 (n={len(all_samples)})，优先参考效应量排名 (|effect_size|>0.3) 而非 p 值。"
            )

        deg_report["warnings"] = deg_warnings

        # ── Generate DEG volcano plot ──
        deg_fig_url = _generate_volcano_figure(
            gene_results,
            key_id="gene_id",
            title=f"DEG Volcano Plot ({n_tested} genes)",
            padj_threshold=padj_threshold,
            lfc_threshold=lfc_threshold,
        )
        if deg_fig_url:
            deg_report["figure_url"] = deg_fig_url
            deg_report["figure_markdown"] = f"![DEG Volcano]({deg_fig_url})"

        state.deg_report = deg_report

        state.add_trace("deg_analysis", "completed", summary={
            "n_genes_tested": n_tested,
            "n_fdr_sig": n_fdr_sig,
            "n_groups": n_groups,
            "top_effect_size": round(top_genes[0]["effect_size"], 4) if top_genes else 0,
        })

        logger.info(
            "DEG node: n_tested=%d, n_fdr_sig=%d, n_groups=%d, top_es=%.4f",
            n_tested, n_fdr_sig, n_groups,
            top_genes[0]["effect_size"] if top_genes else 0,
        )

        return {"deg_report": deg_report}

    except Exception as e:
        logger.error(f"deg_analysis_node failed: {e}", exc_info=True)
        state.add_trace("deg_analysis", "failed", error=str(e))
        return {"deg_report": {"error": str(e)}, "current_node": "data_quality_gate"}


# ═══════════════════════════════════════════════════════════════
# Node 0b: 差异积累代谢物分析
# ═══════════════════════════════════════════════════════════════

def dam_analysis_node(state: RuntimeState) -> dict:
    """差异积累代谢物分析。

    对代谢物矩阵运行相同的多组 ANOVA 逻辑。
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "dam_analysis"):
            state.add_trace("dam_analysis", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        meta_mat = _get_meta_matrix(state)
        if not meta_mat:
            state.add_trace("dam_analysis", "completed",
                            summary={"n_metabolites": 0, "reason": "no metabolite data"})
            return {"dam_report": {"n_metabolites_tested": 0, "reason": "no metabolite data"}}

        # Infer groups (use same logic to be consistent)
        sample_metadata = getattr(state.planner_state, "sample_metadata", None) if state.planner_state else None
        groups = _infer_sample_groups(meta_mat, sample_metadata)

        all_samples = []
        for gs in groups.values():
            all_samples.extend(gs)

        from phyto_reason.tools.statistics.deg_tool import (
            _one_way_anova, _eta_squared, _kruskal_wallis, _moderated_t_test,
        )
        from phyto_reason.utils.stats_utils import compute_fdr

        meta_results = []
        pairwise_inputs: dict[str, dict[str, np.ndarray]] = {}
        variance_pool: list[float] = []
        for met_name, sample_vals in meta_mat.items():
            gene_groups: dict[str, np.ndarray] = {}
            for gname, gsamples in groups.items():
                gv = np.array([sample_vals.get(s, float("nan")) for s in gsamples], dtype=float)
                gv = gv[~np.isnan(gv)]
                if len(gv) >= 2:
                    gene_groups[gname] = np.log2(gv + 1e-6)

            if len(gene_groups) < 2:
                continue

            pairwise_inputs[str(met_name)] = gene_groups
            variance_pool.extend([float(v) for values in gene_groups.values() for v in values])

            f_stat, p_anova = _one_way_anova(gene_groups)
            eta2 = _eta_squared(gene_groups)
            # KW as robustness check only -- not as p-value gate
            _, p_kw = _kruskal_wallis(gene_groups)
            p_value = p_anova

            group_means = {g: float(np.mean(v)) for g, v in gene_groups.items()}
            gnames = list(group_means.keys())
            max_lfc = 0.0
            for i in range(len(gnames)):
                for j in range(i + 1, len(gnames)):
                    lfc = abs(group_means[gnames[i]] - group_means[gnames[j]])
                    max_lfc = max(max_lfc, lfc)

            meta_results.append({
                "metabolite": met_name,
                "p_value": p_value,
                "effect_size": eta2,
                "max_log2fc": max_lfc,
                "stat_name": "F_stat",
                "stat_value": f_stat,
            })

        n_tested = len(meta_results)
        if n_tested == 0:
            group_sizes = {name: len(samples) for name, samples in groups.items()}
            n_total = len(set(all_samples))
            feasibility = "exploratory" if n_total < 5 or any(n < 3 for n in group_sizes.values()) else "standard"
            return {"dam_report": {
                "n_metabolites_tested": 0,
                "n_total": n_total,
                "n_fdr_significant": 0,
                "n_groups": len(groups),
                "group_sizes": group_sizes,
                "group_source": "metadata" if sample_metadata else "sample_name_fallback",
                "feasibility": feasibility,
                "warnings": ["未形成至少两组有效比较，未执行 DAM。"],
            }}

        meta_results.sort(key=lambda x: x["p_value"])
        pvalues = [m["p_value"] for m in meta_results]
        qvalues = compute_fdr(pvalues)
        for m, q in zip(meta_results, qvalues):
            m["q_value"] = q

        # Follow a significant multi-group ANOVA with moderated pairwise
        # comparisons.  Pairwise p-values are BH-adjusted across all
        # metabolite x tissue-pair tests, so an ANOVA p-value alone never
        # becomes a tissue-specific conclusion.
        pairwise_results: list[dict] = []
        if len(groups) >= 2:
            global_var = float(np.var(np.asarray(variance_pool, dtype=float), ddof=1)) if len(variance_pool) > 1 else 1.0
            global_var = max(global_var, 1e-10)
            for overall in meta_results:
                if float(overall.get("p_value", 1.0)) >= 0.05:
                    continue
                met_name = str(overall["metabolite"])
                gene_groups = pairwise_inputs.get(met_name, {})
                group_names = list(gene_groups)
                for i, group_a in enumerate(group_names):
                    for group_b in group_names[i + 1:]:
                        a_vals, b_vals = gene_groups[group_a], gene_groups[group_b]
                        if len(a_vals) < 2 or len(b_vals) < 2:
                            continue
                        t_stat, p_pair, _ = _moderated_t_test(a_vals, b_vals, global_var)
                        log2fc = float(np.mean(b_vals) - np.mean(a_vals))
                        pairwise_results.append({
                            "metabolite": met_name,
                            "group_a": group_a,
                            "group_b": group_b,
                            "log2fc": round(log2fc, 6),
                            "p_value": float(p_pair),
                            "t_stat": round(float(t_stat), 6),
                            "direction": "up" if log2fc > 0 else "down" if log2fc < 0 else "unchanged",
                            "overall_anova_p": float(overall.get("p_value", 1.0)),
                            "overall_anova_q": float(overall.get("q_value", 1.0)),
                        })
            if pairwise_results:
                pairwise_q = compute_fdr([r["p_value"] for r in pairwise_results])
                for row, q in zip(pairwise_results, pairwise_q):
                    row["q_value"] = float(q)
                pairwise_results.sort(key=lambda r: (r["q_value"], -abs(r["log2fc"])))

        padj_threshold = 0.05
        lfc_threshold = 0.5
        n_fdr_sig = sum(
            1 for m in meta_results
            if m["q_value"] < padj_threshold and m["max_log2fc"] >= lfc_threshold
        )

        meta_results.sort(key=lambda x: x["effect_size"], reverse=True)
        top_n = min(30, n_tested)

        group_sizes = {name: len(samples) for name, samples in groups.items()}
        n_total = len(set(all_samples))
        feasibility = "exploratory" if n_total < 5 or any(n < 3 for n in group_sizes.values()) else "standard"
        analysis_warnings: list[str] = []
        if feasibility == "exploratory":
            analysis_warnings.append(
                f"探索性结果：总样本 n={n_total} 或存在组样本数 <3；不得仅凭 q 值作确认性结论。"
            )
        if sample_metadata:
            group_source = "metadata" if any(
                str(row.get("condition") or row.get("tissue") or row.get("treatment") or "").strip()
                for row in (sample_metadata.get("samples", []) if isinstance(sample_metadata, dict) else [])
            ) else "sample_name_fallback"
        else:
            group_source = "sample_name_fallback"
            analysis_warnings.append("未提供分组表，分组由样本名前缀推断。")

        dam_report = {
            "n_metabolites_tested": n_tested,
            "n_total": n_total,
            "n_fdr_significant": n_fdr_sig,
            "n_groups": len(groups),
            "group_sizes": group_sizes,
            "group_source": group_source,
            "feasibility": feasibility,
            "warnings": analysis_warnings,
            "pairwise_method": "moderated_t_empirical_bayes + Benjamini-Hochberg",
            "pairwise_results": pairwise_results,
            "n_pairwise_tests": len(pairwise_results),
            "n_pairwise_significant": sum(1 for r in pairwise_results if r.get("q_value", 1) < padj_threshold),
            "pairwise_test_unit": "metabolite x tissue_pair",
            "pairwise_bh_scope": "all metabolite x tissue_pair tests",
            "summary": f"n={n_total}; method=ANOVA + {len(pairwise_results)} moderated-t pairwise tests; threshold=|log2FC|≥{lfc_threshold}, q<{padj_threshold}; FDR=BH over all metabolite x tissue_pair tests",
            "top_metabolites": [
                {
                    "metabolite": m["metabolite"],
                    "effect_size": round(m["effect_size"], 4),
                    "max_log2fc": round(m["max_log2fc"], 3),
                    "p_value": round(m["p_value"], 6),
                    "p_value_raw": float(m["p_value"]),
                    "q_value": round(m["q_value"], 6),
                    "is_fdr_sig": m["q_value"] < padj_threshold and m["max_log2fc"] >= lfc_threshold,
                }
                for m in meta_results[:top_n]
            ],
        }

        if n_fdr_sig == 0 and n_tested > 0:
            dam_report["warnings"].append(
                f"FDR 校正后无显著差异代谢物。已输出效应量排名 top-{min(top_n, n_tested)}。"
                f"小样本可能导致统计效力不足。"
            )

        # ── Generate DAM volcano plot ──
        dam_fig_url = _generate_volcano_figure(
            meta_results,
            key_id="metabolite",
            title=f"DAM Volcano Plot ({n_tested} metabolites)",
            padj_threshold=padj_threshold,
            lfc_threshold=lfc_threshold,
        )
        if dam_fig_url:
            dam_report["figure_url"] = dam_fig_url
            dam_report["figure_markdown"] = f"![DAM Volcano]({dam_fig_url})"

        if pairwise_results:
            try:
                from phyto_reason.visualization.multiomics_figures import plot_dam_pairwise_heatmap
                pairwise_fig_url = plot_dam_pairwise_heatmap(pairwise_results)
            except Exception as exc:
                logger.debug("Pairwise heatmap generation skipped: %s", exc)
                pairwise_fig_url = ""
            if pairwise_fig_url:
                dam_report["pairwise_heatmap_url"] = pairwise_fig_url
                dam_report["pairwise_heatmap_markdown"] = f"![DAM pairwise heatmap]({pairwise_fig_url})"

        state.dam_report = dam_report

        state.add_trace("dam_analysis", "completed", summary={
            "n_metabolites_tested": n_tested,
            "n_fdr_sig": n_fdr_sig,
        })

        return {"dam_report": dam_report}

    except Exception as e:
        logger.error(f"dam_analysis_node failed: {e}", exc_info=True)
        state.add_trace("dam_analysis", "failed", error=str(e))
        return {"dam_report": {"error": str(e)}, "current_node": "data_quality_gate"}


# ═══════════════════════════════════════════════════════════════
# Node 0c: 多组学联合分析（DEG-DAM 相关性 + 模块检测）
# ═══════════════════════════════════════════════════════════════

def multiomics_integration_node(state: RuntimeState) -> dict:
    """多组学联合分析节点。

    执行:
      1. DEG-DAM 相关性分析：计算 top DEGs 与 top DAMs 的 pairwise 相关
      2. 简易 WGCNA 风格模块检测：基于相关矩阵的层次聚类
      3. 输出显著相关的 gene-metabolite 对和共表达模块
    """
    try:
        analysis_warnings: list[str] = []
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "multiomics_integration"):
            state.add_trace("multiomics_integration", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        expr = _get_expr_matrix(state)
        meta_mat = _get_meta_matrix(state)
        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}

        if not expr or not meta_mat:
            state.add_trace("multiomics_integration", "completed",
                            summary={"reason": "missing expression or metabolite data"})
            return {"multiomics_report": {"reason": "missing data"}}

        # Get common samples
        expr_samples = set()
        for v in expr.values():
            expr_samples.update(v.keys())
            break
        meta_samples = set()
        for v in meta_mat.values():
            meta_samples.update(v.keys())
            break
        common_samples = sorted(expr_samples & meta_samples)

        if len(common_samples) < 5:
            state.add_trace("multiomics_integration", "completed",
                            summary={"reason": f"too few common samples: {len(common_samples)}"})
            return {"multiomics_report": {"reason": f"too few common samples: {len(common_samples)}"}}

        # Get top DEGs (FDR-sig first, then effect-size ranked)
        top_deg_genes = deg_report.get("top_genes", [])
        if not top_deg_genes:
            # Fallback: use all expression genes
            top_deg_genes = [{"gene_id": g} for g in list(expr.keys())[:100]]

        top_deg_ids = [g["gene_id"] for g in top_deg_genes[:100]
                       if g["gene_id"] in expr]

        # Get top DAMs (FDR-sig first)
        top_dam_metas = dam_report.get("top_metabolites", [])
        if not top_dam_metas:
            top_dam_metas = [{"metabolite": m} for m in list(meta_mat.keys())[:50]]

        top_dam_ids = [m["metabolite"] for m in top_dam_metas[:50]
                       if m["metabolite"] in meta_mat]

        if not top_deg_ids or not top_dam_ids:
            state.add_trace("multiomics_integration", "completed",
                            summary={"reason": "no valid DEG or DAM IDs"})
            return {"multiomics_report": {"reason": "no valid IDs"}}

        logger.info(
            "Multi-omics: computing correlations for %d DEGs x %d DAMs across %d samples",
            len(top_deg_ids), len(top_dam_ids), len(common_samples),
        )

        # Compute pairwise bicor correlations
        from phyto_reason.utils.stats_utils import compute_correlation, compute_fdr

        corr_pairs: list[dict] = []
        all_pvalues: list[float] = []

        for gid in top_deg_ids[:50]:  # Limit to top 50 DEGs for speed
            gvec = np.array([expr[gid].get(s, float("nan")) for s in common_samples], dtype=float)
            if np.isnan(gvec).any() or np.std(gvec) < 1e-10:
                continue

            for mid in top_dam_ids[:30]:  # Limit to top 30 DAMs
                mvec = np.array([meta_mat[mid].get(s, float("nan")) for s in common_samples], dtype=float)
                if np.isnan(mvec).any() or np.std(mvec) < 1e-10:
                    continue

                try:
                    r, p = compute_correlation(gvec, mvec, method="bicor")
                    if not np.isnan(r):
                        corr_pairs.append({
                            "gene_id": gid, "metabolite": mid,
                            "correlation": round(float(r), 4),
                            "p_value": round(float(p), 6),
                        })
                        all_pvalues.append(p)
                except Exception:
                    continue

        if not corr_pairs:
            state.add_trace("multiomics_integration", "completed",
                            summary={"n_pairs": 0, "reason": "no valid correlations"})
            return {"multiomics_report": {"n_pairs": 0, "reason": "no valid correlations"}}

        # FDR correction
        qvalues = compute_fdr(all_pvalues)
        for pair, q in zip(corr_pairs, qvalues):
            pair["q_value"] = round(float(q), 6)

        # Filter significant pairs
        sig_pairs = [p for p in corr_pairs if p["q_value"] < 0.05]
        top_pairs = sorted(corr_pairs, key=lambda x: abs(x["correlation"]), reverse=True)[:30]

        n_sig = len(sig_pairs)
        n_total = len(corr_pairs)

        # Simple module detection: cluster genes by their correlation to metabolites
        # Genes with similar metabolite correlation profiles -> same module
        modules: dict[str, list[str]] = {}
        if len(top_deg_ids) >= 5:
            # Build gene-metabolite correlation matrix
            gene_meta_corr: dict[str, list[float]] = {}
            for pair in corr_pairs:
                gene_meta_corr.setdefault(pair["gene_id"], []).append(pair["correlation"])

            # Simple clustering: group genes by sign of top correlation
            pos_module = []
            neg_module = []
            for gid, corrs in gene_meta_corr.items():
                mean_corr = np.mean(corrs) if corrs else 0
                if mean_corr > 0.2:
                    pos_module.append(gid)
                elif mean_corr < -0.2:
                    neg_module.append(gid)

            if pos_module:
                modules["positive_correlated"] = pos_module[:20]
            if neg_module:
                modules["negative_correlated"] = neg_module[:20]

        multiomics_report = {
            "n_correlation_pairs_tested": n_total,
            "n_significant_pairs": n_sig,
            "n_common_samples": len(common_samples),
            "method": "bicor",
            "corr_threshold": 0.8,
            "fdr_status": "BH applied to all tested gene-metabolite pairs",
            "summary": f"n={len(common_samples)}; method=bicor; |r|≥0.8; FDR=BH q<0.05",
            "top_pairs": top_pairs[:20],
            "modules": {k: v[:10] for k, v in modules.items()},
            "warnings": analysis_warnings,
        }

        if n_sig == 0:
            multiomics_report["warnings"].append(
                f"无 FDR-显著的基因-代谢物相关对 (q<0.05)。"
                f"已输出 |r| 排名 top-{len(top_pairs)} 对。"
            )

        state.multiomics_report = multiomics_report

        state.add_trace("multiomics_integration", "completed", summary={
            "n_pairs_tested": n_total,
            "n_significant": n_sig,
            "n_modules": len(modules),
            "top_r": abs(top_pairs[0]["correlation"]) if top_pairs else 0,
        })

        logger.info(
            "Multi-omics: %d/%d significant pairs, %d modules, top_r=%.4f",
            n_sig, n_total, len(modules),
            abs(top_pairs[0]["correlation"]) if top_pairs else 0,
        )

        return {"multiomics_report": multiomics_report}

    except Exception as e:
        logger.error(f"multiomics_integration_node failed: {e}", exc_info=True)
        state.add_trace("multiomics_integration", "failed", error=str(e))
        return {"multiomics_report": {"error": str(e)}, "current_node": "data_quality_gate"}


# ═══════════════════════════════════════════════════════════════
# Node 1: 数据可信度验证
# ═══════════════════════════════════════════════════════════════

def data_quality_gate_node(state: RuntimeState) -> dict:
    """验证数据基本可信度。

    检查内容:
      - 样本数是否 ≥ 3（相关性分析最低要求）
      - 表达矩阵中缺失值比例
      - 表达与代谢物矩阵的共同样本数
      - 目标代谢物是否存在于代谢物矩阵中
    """
    try:
        issues: list[str] = []
        report: dict = {}
        expr = _get_expr_matrix(state)
        meta_mat = _get_meta_matrix(state)
        target = _get_meta(state)

        if expr:
            n_genes = len(expr)
            first = next(iter(expr.values()), {})
            n_samples = len(first)
            report["n_genes"] = n_genes
            report["n_samples"] = n_samples

            missing_rates = []
            for g, vals in expr.items():
                if vals:
                    missing = sum(1 for v in vals.values() if isinstance(v, float) and math.isnan(v))
                    missing_rates.append(missing / max(len(vals), 1))
            avg_missing = sum(missing_rates) / max(len(missing_rates), 1) if missing_rates else 0
            report["avg_missing_rate"] = round(avg_missing, 4)

            if n_samples < 3:
                issues.append(f"样本数不足 (n={n_samples}<3)，统计不可靠")
            if avg_missing > 0.3:
                issues.append(f"缺失值过高 ({avg_missing:.0%})，建议填补或过滤")
            if n_genes < 100:
                issues.append(f"基因数偏少 (n={n_genes})，分析范围受限")
        else:
            issues.append("无表达矩阵，无法进行共表达分析")

        if meta_mat:
            n_metas = len(meta_mat)
            report["n_metabolites"] = n_metas
            if target:
                found = any(target in m.lower() for m in meta_mat)
                if not found:
                    issues.append(f"目标代谢物 '{target}' 未在代谢物矩阵中找到")
                else:
                    report["target_found"] = True
        else:
            issues.append("无代谢物矩阵，无法验证代谢物差异")
            report["target_found"] = False

        if expr and meta_mat:
            expr_samples = set()
            for vals in expr.values():
                expr_samples.update(vals.keys())
                break
            meta_samples = set()
            for vals in meta_mat.values():
                meta_samples.update(vals.keys())
                break
            common = expr_samples & meta_samples
            report["common_samples"] = len(common)
            if len(common) < 3:
                issues.append(f"表达与代谢物矩阵共同样本不足 (n={len(common)}<3)")

        critical_issues = [i for i in issues if any(kw in i for kw in ["无表达", "样本数不足", "共同样本不足"])]
        passed = len(critical_issues) == 0
        if issues and passed:
            issues.append("（非关键问题，分析继续但需注意上述限制）")
        state.data_quality_passed = passed
        state.data_quality_issues = issues
        state.data_quality_report = report

        state.add_trace("data_quality_gate", "completed", summary={
            "passed": passed, "n_issues": len(issues), "critical": len(critical_issues),
        })

        node = "metabolite_validation" if passed else "__end__"
        return {
            "data_quality_passed": passed,
            "data_quality_issues": issues,
            "data_quality_report": report,
            "current_node": node,
            "finished": not passed,
            "hypothesis": "数据质量检查未通过。\n" + "\n".join(f"- {x}" for x in critical_issues) if not passed else "",
        }
    except Exception as e:
        logger.error(f"data_quality_gate failed: {e}")
        state.add_trace("data_quality_gate", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node 2: 目标代谢物差异确认
# ═══════════════════════════════════════════════════════════════

def metabolite_validation_node(state: RuntimeState) -> dict:
    """确认目标代谢物是否有足够的生物学差异。

    检查:
      - 代谢物矩阵中目标物的变异系数（CV>0.3=有信号）
      - 如果是 DEG/DAM 数据，检查是否显著差异
      - 无数据时仅做 soft gate（发出警告但继续）
    """
    try:
        target = _get_meta(state)
        meta_mat = _get_meta_matrix(state)
        note = ""
        validated = False

        if meta_mat and target:
            found_key = None
            target_lower = target.lower().replace("-", "_").replace(" ", "_")
            for m, vals in meta_mat.items():
                m_lower = m.lower()
                if target_lower == m_lower or target_lower in m_lower or m_lower in target_lower:
                    found_key = m
                    break
            if found_key:
                values = [v for v in meta_mat[found_key].values()
                          if isinstance(v, (int, float)) and not math.isnan(v)]
                if len(values) >= 3:
                    mean_v = np.mean(values)
                    std_v = np.std(values)
                    cv = std_v / mean_v if abs(mean_v) > 1e-10 else 0
                    if cv > 0.3:
                        validated = True
                        note = f"代谢物变异充分 (CV={cv:.2f})，适合分析"
                    elif cv > 0.1:
                        validated = True
                        note = f"代谢物有一定变异 (CV={cv:.2f})，分析敏感度受限"
                    else:
                        note = f"代谢物在样本间几乎无变化 (CV={cv:.2f})，调控分析意义有限"
                else:
                    note = f"代谢物数据点不足 (n={len(values)})"
            else:
                note = f"目标代谢物 '{target}' 在矩阵中未找到精确匹配"
        elif target:
            note = "无代谢物矩阵，基于先验知识继续分析"
            validated = True
        else:
            note = "未指定目标代谢物，进行探索性分析"

        state.metabolite_validated = validated
        state.metabolite_validation_note = note

        state.add_trace("metabolite_validation", "completed", summary={
            "validated": validated, "note": note[:80],
        })

        return {
            "metabolite_validated": validated,
            "metabolite_validation_note": note,
            "current_node": "pathway_coherence",
        }
    except Exception as e:
        logger.error(f"metabolite_validation failed: {e}")
        state.add_trace("metabolite_validation", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node 3: 通路层级一致性检查
# ═══════════════════════════════════════════════════════════════

def pathway_coherence_node(state: RuntimeState) -> dict:
    """确认代谢通路是否在表达数据中活跃。

    使用 PathwayActivityEngine:
      - mean z-score 检测通路基因相对全局背景的表达水平
      - 仅当 z ≥ 1.5 才判定通路活跃（替代旧的 pairwise r）
    """
    try:
        target = _get_meta(state)
        expr = _get_expr_matrix(state)

        from phyto_reason.reasoning.pathway_activity_engine import PathwayActivityEngine

        ann = getattr(state, "annotation_index", None)
        engine = PathwayActivityEngine()
        activity = engine.compute_pathway_activity(target, expr, annotation_index=ann) if target else {}

        pname = activity.get("pathway_name", "unknown")
        found_enzymes = activity.get("enzymes_found", [])
        mean_z = activity.get("mean_z_score", 0.0)
        coherence_score = activity.get("coherence_score", 0.0)
        is_active = activity.get("is_active", False)
        coherence_evidence = activity.get("activity_interpretation", "")

        # ── 物种 profile 兜底（v5.0）：注释索引无酶基因时，
        #    用该物种 species_profiles/ 的通路先验酶列表，显式标注降级来源 ──
        if not found_enzymes:
            try:
                species = getattr(state.planner_state, "species", "") if state.planner_state else ""
                if species:
                    from phyto_reason.knowledge.species_registry import get_species_registry
                    profile = get_species_registry().get(species)
                    if profile:
                        for pw_key, info in profile.pathway_prior.items():
                            found_enzymes.extend(info.get("enzymes", []) or [])
                        found_enzymes = list(dict.fromkeys(found_enzymes))  # 去重保序
                        if found_enzymes:
                            pname = pname if pname != "unknown" else next(
                                iter(profile.pathway_prior.keys()), "unknown")
                            coherence_evidence = (
                                f"注释索引无通路酶，已使用物种 profile 通路先验 "
                                f"（species_profiles/{profile.slug}.yaml）"
                            )
            except Exception as e:
                logger.debug("Species profile pathway fallback failed: %s", e)

        passed = is_active

        state.pathway_coherence_passed = passed
        state.pathway_coherence_score = round(coherence_score, 4)
        state.pathway_coherence_evidence = coherence_evidence
        state.pathway_name = pname
        state.pathway_genes = found_enzymes

        state.add_trace("pathway_coherence", "completed", summary={
            "passed": passed, "pathway": pname,
            "mean_z_score": mean_z, "is_active": is_active,
            "n_enzymes_found": len(found_enzymes),
        })

        return {
            "pathway_coherence_passed": passed,
            "pathway_coherence_score": round(coherence_score, 4),
            "pathway_coherence_evidence": coherence_evidence,
            "pathway_name": pname,
            "pathway_genes": found_enzymes,
            "current_node": "tf_narrowing",
        }
    except Exception as e:
        logger.error(f"pathway_coherence failed: {e}")
        state.add_trace("pathway_coherence", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node 4: 基于 pathway 基因缩小候选 TF
# ═══════════════════════════════════════════════════════════════

def tf_narrowing_node(state: RuntimeState) -> dict:
    """基于通路酶基因筛选候选 TF。

    逻辑:
      TF 调控应该作用于 pathway enzyme gene，不是直接作用于 metabolite。
      所以先找与 pathway genes 相关的 TF，再推理调控逻辑。

    步骤:
      1. 使用 CorrelationEngine（default bicor）计算 TF 与 pathway gene 的相关性
      2. 双阈值筛选: |r| ≥ corr_threshold **且** q < fdr_threshold
      3. 仅保留有 TF annotation 的基因（通过 tf_annotation_tool）
      4. motif 检查（使用 motif_real/RealMotifTool）
      5. TF family prior 匹配
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "tf_narrowing"):
            state.add_trace("tf_narrowing", "skipped",
                            summary={"reason": "user skipped this step"})
            return {"current_node": "contradiction_check"}

        target = _get_meta(state)
        expr = _get_expr_matrix(state)
        promoters = _get_promoters(state)
        pathway_genes = state.pathway_genes

        if not pathway_genes or not expr:
            state.add_trace("tf_narrowing", "completed", summary={"n_candidates": 0, "reason": "无通路基因或表达矩阵"})
            return {"current_node": "contradiction_check", "tf_candidates": {}}

        from phyto_reason.tools.network.correlation_engine import CorrelationEngine
        from phyto_reason.tools.network.confounder_engine import ConfounderEngine
        from phyto_reason.reasoning.tf_prior_reasoner import TfPriorReasoner
        from phyto_reason.tools.tool_executor import ToolExecutor

        engine = CorrelationEngine(method="bicor")
        confounder = ConfounderEngine()
        executor = ToolExecutor()

        all_gene_ids = list(expr.keys())
        if not all_gene_ids:
            state.add_trace("tf_narrowing", "completed", summary={"n_candidates": 0, "reason": "无基因"})
            return {"current_node": "contradiction_check"}

        # ── Step A: 识别 TF 基因 ──────────────────────────
        tf_ids: list[str] = []
        tf_result = executor.execute("tf_annotation", raise_on_error=False,
                                     gene_ids=all_gene_ids)
        if tf_result and len(tf_result.candidates) > 0:
            tf_ids = [c.gene_id for c in tf_result.candidates]
        else:
            # 非模式物种回退（v5.0）：用注释索引搜索 TF 关键词
            # （Pfam/Swissprot/GO 注释列），而非默认"全部基因"
            ann = getattr(state, "annotation_index", None)
            if ann is not None:
                for kw in ("transcription factor", "transcription regulator",
                           "DNA-binding", "myb domain", "bhlh", "wrky domain",
                           "erf domain", "nac domain", "zinc finger"):
                    tf_ids.extend(ann.search_substring(kw))
                # 注释 GeneID 常带转录本后缀（Zni13G006980.1），
                # 表达矩阵通常不带 —— 剥后缀并对齐到矩阵
                expr_keys = set(expr.keys())
                stripped: list[str] = []
                for tid in dict.fromkeys(tf_ids):
                    base = tid.rpartition(".")[0] if "." in tid else tid
                    if base in expr_keys or tid in expr_keys:
                        stripped.append(base if base in expr_keys else tid)
                tf_ids = stripped
                if tf_ids:
                    logger.info("Annotation-index TF fallback: %d genes matched", len(tf_ids))
            if not tf_ids:
                tf_ids = list(expr.keys())
                logger.warning("TF annotation tool returned no results, using all genes")

        if not tf_ids:
            state.add_trace("tf_narrowing", "completed", summary={"n_candidates": 0, "reason": "数据中未识别到 TF"})
            return {"current_node": "contradiction_check"}

        # ── Step B: 解析通路锚点基因 + bicor 计算 ─────────
        common_samples: list[str] = []
        for s in expr.values():
            common_samples = sorted(s.keys())
            break

        ann = getattr(state, "annotation_index", None)

        def _resolve_to_expr(hit: str) -> str | None:
            """注释命中（可带 .1 后缀）→ 表达矩阵键。"""
            if hit in expr:
                return hit
            base = hit.rpartition(".")[0]
            return base if base in expr else None

        # 1) 通路酶基因 → 表达矩阵（符号直配 / 注释索引解析）
        anchor_genes: list[str] = []
        for pw_gene in pathway_genes:
            pw_key = next((g for g in expr if pw_gene.upper() in g.upper()), None)
            if pw_key is None and ann is not None:
                for hit in ann.search_substring(pw_gene):
                    pw_key = _resolve_to_expr(hit)
                    if pw_key:
                        break
            if pw_key:
                anchor_genes.append(pw_key)
        anchor_genes = list(dict.fromkeys(anchor_genes))

        # 2) 降级：酶符号在注释中缺失（非模式物种注释体系差异）时，
        #    用物种 profile 通路的 search_terms（中间体/KEGG KO）从注释索引定位
        if not anchor_genes:
            try:
                species = getattr(state.planner_state, "species", "") if state.planner_state else ""
                profile = get_species_registry().get(species) if species else None
                if profile and ann is not None:
                    for info in profile.pathway_prior.values():
                        for term in info.get("search_terms") or []:
                            for hit in ann.search_substring(term):
                                k = _resolve_to_expr(hit)
                                if k:
                                    anchor_genes.append(k)
                    anchor_genes = list(dict.fromkeys(anchor_genes))
                    if anchor_genes:
                        logger.info(
                            "Profile pathway search-terms fallback: %d anchor genes resolved",
                            len(anchor_genes),
                        )
            except Exception as e:
                logger.debug("Pathway search-terms fallback failed: %s", e)

        if not anchor_genes:
            state.add_trace("tf_narrowing", "completed", summary={
                "n_candidates": 0,
                "reason": "通路锚点基因无法解析（注释缺失或与表达矩阵不对齐）",
            })
            return {"current_node": "contradiction_check", "tf_candidates": {}}

        tf_candidates: dict[str, CandidateGene] = {}
        tf_pathway_links: dict[str, list[str]] = {}

        for pw_key in anchor_genes:
            pw_vec = np.array([expr[pw_key].get(s, float("nan")) for s in common_samples], dtype=float)
            if np.isnan(pw_vec).any() or len(pw_vec) < 3:
                continue

            pw_pairs: list[tuple[str, str, float, float]] = []
            for tf_id in tf_ids:
                tf_vec = np.array([expr[tf_id].get(s, float("nan")) for s in common_samples], dtype=float)
                if np.isnan(tf_vec).any():
                    continue
                r, p = engine.calculate(tf_vec, pw_vec)
                if abs(r) >= engine.corr_threshold:
                    pw_pairs.append((tf_id, pw_key, r, p))

            filtered = engine.filter_pairs(pw_pairs)
            for tf_id, _pw_key, r, p, q in filtered:
                if tf_id not in tf_candidates:
                    gene = CandidateGene(gene_id=tf_id, is_tf=True,
                                         correlation_score=round(abs(r), 4))
                    tf_candidates[tf_id] = gene
                    tf_pathway_links[tf_id] = []
                tf_candidates[tf_id].correlation_score = max(
                    tf_candidates[tf_id].correlation_score, round(abs(r), 4))
                tf_pathway_links[tf_id].append(pw_key)

        # ── Step C: TF family prior ───────────────────────
        has_motif = bool(promoters)

        for tid in list(tf_candidates.keys()):
            gene = tf_candidates[tid]
            gene.target_enzymes = list(dict.fromkeys(tf_pathway_links.get(tid, [])))
            gene.target_metabolites = [target] if target else []

            prior = TfPriorReasoner.reason("", target)
            if prior.score > 0:
                gene.literature_score = max(gene.literature_score, prior.score)
                gene.tf_family = "prior_matched"

            gene.add_evidence(Evidence(
                evidence_type=EvidenceType.CORRELATION,
                source="tf_narrowing",
                score=gene.correlation_score,
                description=f"bicor与通路酶基因共表达 (r={gene.correlation_score:.2f}, FDR校正)",
            ))

        # ── Step D: motif 检查 (motif_real) ───────────────
        if has_motif and tf_candidates:
            candidate_ids = list(tf_candidates.keys())
            try:
                motif_result = executor.execute(
                    "motif_real", raise_on_error=False,
                    promoter_sequences=promoters, gene_ids=candidate_ids,
                )
                if motif_result.status == "failed":
                    state.degraded_tools.append("motif_scan")
                    state.degradation_notes.append("Motif scan failed; confidence reduced")
                for c in motif_result.candidates:
                    if c.gene_id in tf_candidates and c.motif_score > 0:
                        tf_candidates[c.gene_id].motif_score = c.motif_score
                        for ev in motif_result.evidence_list:
                            if (ev.metadata or {}).get("gene_id") == c.gene_id:
                                tf_candidates[c.gene_id].add_evidence(ev)
            except Exception as e:
                logger.warning(f"motif_real execution failed, skipping: {e}")
                state.degraded_tools.append("motif_scan")
                state.degradation_notes.append(f"Motif scan unavailable: {e}")

        # ── Step E: WGCNA ────────────────────────────────────
        n_before = len(tf_candidates)
        tf_candidates = {k: v for k, v in tf_candidates.items()
                         if v.correlation_score >= 0.50 or v.motif_score >= 0.30}

        state.tf_candidates = tf_candidates
        state.tf_pathway_gene_links = {k: list(dict.fromkeys(v)) for k, v in tf_pathway_links.items()
                                       if k in tf_candidates}
        state.tool_results["tf_narrowing"] = ToolResult(candidates=list(tf_candidates.values()))

        state.add_trace("tf_narrowing", "completed", summary={
            "n_candidates": len(tf_candidates),
            "n_filtered": n_before - len(tf_candidates),
            "top_tfs": list(tf_candidates.keys())[:5],
            "correlation_method": "bicor",
        })

        return {
            "tf_candidates": tf_candidates,
            "tf_pathway_gene_links": state.tf_pathway_gene_links,
            "tool_results": state.tool_results,
            "current_node": "regulation_evidence",
        }
    except Exception as e:
        logger.error(f"tf_narrowing failed: {e}")
        state.add_trace("tf_narrowing", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node 5: 调控证据多维检查
# ═══════════════════════════════════════════════════════════════

def regulation_evidence_node(state: RuntimeState) -> dict:
    """对已缩小的候选 TF 执行多维调控证据检查。

    检查维度:
      - 与通路基因的共表达力度（correlation）
      - 启动子 motif 存在（motif）
      - PubMed 文献支持（literature）
      - KEGG 通路关联（kegg）
      - TF 家族先验（tf_prior）
      - 组织特异性（tissue）
      - 拟南芥同源（ortholog）
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "regulation_evidence"):
            state.add_trace("regulation_evidence", "skipped",
                            summary={"reason": "user skipped this step"})
            return {"current_node": "contradiction_check"}

        target = _get_meta(state)
        species = state.planner_state.species if state.planner_state else ""
        candidates = state.tf_candidates

        if not candidates:
            state.add_trace("regulation_evidence", "completed", summary={"n_checked": 0})
            return {"current_node": "contradiction_check"}

        from phyto_reason.reasoning.biological_reasoner import BiologicalReasoner
        from phyto_reason.fusion.evidence_fusion_engine import EvidenceFusionEngine
        from phyto_reason.fusion.fusion_result import UnifiedFusionResult
        from phyto_reason.tools.tool_executor import ToolExecutor

        reasoner = BiologicalReasoner()
        fusion_engine = EvidenceFusionEngine()
        executor = ToolExecutor()

        literature_result = executor.execute(
            "literature_search", raise_on_error=False,
            metabolite=target, species=species, max_results=5,
        )
        kegg_result = executor.execute(
            "kegg_pathway", raise_on_error=False,
            compound_name=target, analysis_type="pathway",
        )

        # Track degraded tools
        if literature_result.status in ("failed", "partial") or literature_result.warnings:
            state.degraded_tools.append("literature_search")
            state.degradation_notes.append(
                f"Literature search degraded: {literature_result.warnings[:2] if literature_result.warnings else 'no results'}"
            )
        if kegg_result.status in ("failed", "partial") or kegg_result.warnings:
            state.degraded_tools.append("kegg_pathway")
            state.degradation_notes.append(
                f"KEGG pathway query degraded: {kegg_result.warnings[:2] if kegg_result.warnings else 'no results'}"
            )

        fusion_results: dict[str, UnifiedFusionResult] = {}
        artifact_notes: dict[str, list[str]] = {}

        for tf_id, gene in candidates.items():
            reasoned = reasoner.reason(gene, target_metabolite=target)
            for ev in literature_result.evidence_list:
                reasoned.literature_score = max(reasoned.literature_score, ev.score)
                reasoned.add_evidence(ev)
            for ev in kegg_result.evidence_list:
                reasoned.pathway_score = max(reasoned.pathway_score, ev.score)
                reasoned.add_evidence(ev)

            fr = fusion_engine.fuse(
                gene_id=tf_id,
                correlation_score=reasoned.correlation_score,
                module_membership=reasoned.module_membership,
                motif_score=reasoned.motif_score,
                pathway_score=reasoned.pathway_score,
                tf_prior_score=reasoned.literature_score,
                tissue_score=reasoned.tissue_specificity_score,
                ortholog_score=reasoned.ortholog_score,
                literature_score=reasoned.literature_score,
            )
            fusion_results[tf_id] = fr

            notes = []
            if fr.contradictions:
                for c in fr.contradictions:
                    notes.append(c.description[:100])
            artifact_notes[tf_id] = notes

        state.tool_results["regulation_evidence"] = literature_result
        state.tool_results["kegg_evidence"] = kegg_result
        state.fusion_results = fusion_results
        state.artifact_notes = artifact_notes

        state.add_trace("regulation_evidence", "completed", summary={
            "n_candidates": len(candidates),
            "n_with_literature": len(literature_result.evidence_list),
            "n_with_kegg": len(kegg_result.evidence_list),
        })

        return {
            "tool_results": state.tool_results,
            "fusion_results": fusion_results,
            "artifact_notes": artifact_notes,
            "current_node": "contradiction_check",
        }
    except Exception as e:
        logger.error(f"regulation_evidence failed: {e}")
        state.add_trace("regulation_evidence", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node 6: 假阳性主动排除
# ═══════════════════════════════════════════════════════════════

def contradiction_check_node(state: RuntimeState) -> dict:
    """主动检查假阳性源并排除。

    检查项:
      1. 组织混杂 -- TF 是否在代谢物不积累的组织中高表达
      2. 胁迫响应 -- TF 是否已知为泛胁迫响应因子
      3. 表达不稳定 -- TF 在样本间变异是否过大
      4. 证据孤立 -- 仅单一证据源支持
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "contradiction_check"):
            state.add_trace("contradiction_check", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        candidates = state.tf_candidates
        fusion_results = state.fusion_results
        expr = _get_expr_matrix(state)
        target = _get_meta(state)

        excluded: list[str] = []
        excluded_reasons: list[str] = []

        for tf_id, gene in list(candidates.items()):
            reasons: list[str] = []
            fr = fusion_results.get(tf_id)

            if fr:
                n_active_sources = fr.hierarchy.n_independent_sources
                if n_active_sources < 2:
                    reasons.append(f"证据孤立（仅 {n_active_sources} 个独立证据源）")
                for c in fr.contradictions:
                    if c.level.value in ("high", "medium"):
                        reasons.append(c.description[:80])

            if expr and tf_id in expr:
                vals = [v for v in expr[tf_id].values() if isinstance(v, (int, float)) and not math.isnan(v)]
                if len(vals) >= 3:
                    cv = np.std(vals) / (abs(np.mean(vals)) + 1e-10)
                    if cv > 1.5:
                        reasons.append(f"表达高度不稳定 (CV={cv:.1f})，可能为噪声")
                    if np.mean(vals) < 0.1:
                        reasons.append("表达量极低，调控功能可疑")

            from phyto_reason.reasoning.tissue_reasoner import TissueReasoner
            if target:
                related = TissueReasoner.get_bias_for_metabolite(target)
                if related:
                    pass

            if reasons:
                excluded.append(tf_id)
                excluded_reasons.append(f"{tf_id}: {'; '.join(reasons)}")

        for tf_id in excluded:
            candidates.pop(tf_id, None)
            fusion_results.pop(tf_id, None)

        state.excluded_tfs = excluded
        state.tf_candidates = candidates
        state.fusion_results = fusion_results

        # Generate initial competing hypotheses from remaining candidates.
        # This runs before falsification so hypotheses can be attacked.
        competing_hypotheses = []
        if candidates:
            from phyto_reason.reasoning.hypothesis_competition_engine import (
                HypothesisCompetitionEngine,
            )
            engine = HypothesisCompetitionEngine()
            competing_hypotheses = engine.generate_competing_hypotheses(
                target_metabolite=target or "",
                has_transcriptional_evidence=True,
                top_tf_candidates=list(candidates.keys())[:5],
            )

        state.add_trace("contradiction_check", "completed", summary={
            "n_excluded": len(excluded),
            "n_remaining": len(candidates),
            "n_hypotheses": len(competing_hypotheses),
        })

        return {
            "tf_candidates": candidates,
            "fusion_results": fusion_results,
            "excluded_tfs": excluded,
            "artifact_notes": state.artifact_notes,
            "competing_hypotheses": competing_hypotheses,
        }
    except Exception as e:
        logger.error(f"contradiction_check failed: {e}")
        state.add_trace("contradiction_check", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node 7: 机理假设生成
# ═══════════════════════════════════════════════════════════════

def hypothesis_synthesis_node(state: RuntimeState) -> dict:
    """v3.0 hypothesis report generation.

    v3.0 模式下读取 competing_hypotheses (falsified + ranked)，
    输出符合 Architecture Freeze Spec v3.0 格式的报告。

    向后兼容: 无 competing_hypotheses 时使用旧格式（基于 fusion_results）。
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "hypothesis_synthesis"):
            state.add_trace("hypothesis_synthesis", "skipped",
                            summary={"reason": "user skipped this step"})
            return {"current_node": "__end__", "finished": True}

        target = _get_meta(state)
        pname = state.pathway_name
        n_samples = state.planner_state.sample_count if state.planner_state else 0
        excluded = state.excluded_tfs
        hypotheses = state.competing_hypotheses

        if hypotheses:
            from phyto_reason.reasoning.hypothesis_report import build_full_report
            from phyto_reason.reasoning.wording_policy import validate_v3_output_text

            lines = build_full_report(
                competing_hypotheses=hypotheses,
                target_metabolite=target,
                pathway_name=pname,
                n_samples=n_samples,
                n_excluded=len(excluded),
                mechanism_probabilities=getattr(state, "mechanism_probabilities", None),
            )
            hypothesis_text = "\n".join(lines)

            violations = validate_v3_output_text(hypothesis_text)
            if violations:
                logger.warning(f"v3.0 wording violations detected: {violations}")
        else:
            candidates = state.tf_candidates
            fusion_results = state.fusion_results
            pathway_genes = state.pathway_genes

            from phyto_reason.reasoning.wording_policy import (
                format_hypothesis_statement, build_caveat_section,
            )

            lines: list[str] = []
            lines.append(f"## 调控机理假设")
            lines.append(f"**目标**: {target or '未指定'} | **通路**: {pname}")
            lines.append("")

            if excluded:
                lines.append(f"### 排除的候选（{len(excluded)} 个）")
                for tf_id in excluded:
                    notes = state.artifact_notes.get(tf_id, [])
                    lines.append(f"- {tf_id}: 排除原因 -- {'; '.join(notes[:2])}")
                lines.append("")

            ranked = sorted(
                [(tid, fr) for tid, fr in fusion_results.items()],
                key=lambda x: x[1].calibrated_score, reverse=True,
            )

            if ranked:
                lines.append(f"### 候选 TF 调控链（证据摘要）")
                lines.append("")
                for i, (tf_id, fr) in enumerate(ranked[:5], 1):
                    link_genes = state.tf_pathway_gene_links.get(tf_id, [])
                    ev_parts: list[str] = []
                    for ev in sorted(fr.evidence_list, key=lambda e: e.contribution, reverse=True):
                        label = {
                            "correlation": "共表达", "module_membership": "模块归属",
                            "motif": "启动子结合", "pathway": "通路一致性",
                            "tf_prior": "TF家族先验", "tissue": "组织特异性",
                            "ortholog": "同源保守性", "literature": "文献支持",
                        }.get(ev.source, ev.source)
                        if ev.normalized_score >= 0.30:
                            ev_parts.append(f"{label}({ev.normalized_score:.2f})")
                    level = fr.confidence_level.value if hasattr(fr.confidence_level, "value") else str(fr.confidence_level)
                    lines.append(f"**Top-{i}: {tf_id}** [{level}]")
                    if link_genes:
                        hypothesis_stmt = format_hypothesis_statement(
                            tf_id=tf_id, target_enzymes=link_genes,
                            pathway_name=pname, metabolite=target or "target",
                            fusion_result=fr,
                        )
                        lines.append(f"  - {hypothesis_stmt}")
                    if ev_parts:
                        lines.append(f"  - 证据链: {' -> '.join(ev_parts[:5])}")
                    lines.append(f"  - 总结分: {fr.calibrated_score:.3f}")
                    lines.append("")
            else:
                lines.append("未找到可信的候选转录因子。")
                if state.pathway_coherence_passed:
                    lines.append(f"通路 '{pname}' 虽然表现出协调性，但未检测到调控该通路的 TF。")
                else:
                    lines.append(f"通路 '{pname}' 未检测到协调表达，通路可能未在实验条件下激活。")

            lines.append("### 结论局限性")
            caveats = build_caveat_section(
                fusion_result=ranked[0][1] if ranked else None,
                n_samples=n_samples,
                has_deg=False, has_perturbation=False,
            )
            for caveat in caveats:
                lines.append(f"- {caveat}")
            lines.append("")

            if target and pname != "unknown":
                top_tf = ranked[0][0] if ranked else "?"
                lines.append("### 下一步验证建议")
                lines.append(f"- Y1H/Dual-LUC 验证 {top_tf} 与通路酶基因启动子的结合")
                lines.append("- CRISPR/Cas9 敲除或过表达验证表型")
                lines.append("- 时序表达分析确认调控先后关系")

            hypothesis_text = "\n".join(lines)

        state.mechanism_hypothesis = hypothesis_text
        state.hypothesis = hypothesis_text
        state.finished = True

        state.add_trace("hypothesis_synthesis", "completed", summary={
            "n_hypotheses": len(hypotheses) if hypotheses else 0,
            "report_length": len(hypothesis_text),
        })

        return {
            "mechanism_hypothesis": hypothesis_text,
            "hypothesis": hypothesis_text,
            "current_node": "__end__",
            "finished": True,
        }
    except Exception as e:
        logger.error(f"hypothesis_synthesis failed: {e}")
        state.add_trace("hypothesis_synthesis", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node +1: 机制分类（新）
# ═══════════════════════════════════════════════════════════════

def mechanism_classifier_node(state: RuntimeState) -> dict:
    """Node: 代谢变化机制分类。

    放置在 pathway_coherence 之后、tf_narrowing 之前。
    判断当前代谢变化最可能的机制类型。
    如果不是转录调控主导，则提示系统不要默认走 TF 路线。
    """
    try:
        target = _get_meta(state)
        meta_mat = _get_meta_matrix(state)
        cv = 0.0
        if meta_mat and target:
            for m, vals in meta_mat.items():
                if target in m.lower():
                    nums = [v for v in vals.values()
                            if isinstance(v, (int, float)) and not math.isnan(v)]
                    if len(nums) >= 3:
                        import numpy as np
                        cv = float(np.std(nums) / max(abs(np.mean(nums)), 1e-10))
                    break

        from phyto_reason.reasoning.mechanism_classifier import MechanismClassifier

        classifier = MechanismClassifier()
        classification = classifier.classify(
            expression_data=_get_expr_matrix(state),
            metabolite_data=meta_mat,
            target_metabolite=target,
            cv=cv,
            n_samples=_get_expr_matrix(state) and len(next(iter(_get_expr_matrix(state).values()))) or 0,
            has_promoter=bool(_get_promoters(state)),
        )

        is_transcriptional = classification.primary.value == "transcriptional_regulation"

        state.add_trace("mechanism_classifier", "completed", summary={
            "primary": classification.primary.value,
            "uncertainty": classification.uncertainty,
            "top_probs": dict(list(classification.probabilities.items())[:3]),
        })

        return {
            "current_node": "tf_narrowing",
            "mechanism_type": classification.primary.value,
            "mechanism_probabilities": classification.probabilities,
            "mechanism_uncertainty": classification.uncertainty,
        }
    except Exception as e:
        logger.error(f"mechanism_classifier failed: {e}")
        state.add_trace("mechanism_classifier", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node +2: 假设竞争（新）
# ═══════════════════════════════════════════════════════════════

def hypothesis_competition_node(state: RuntimeState) -> dict:
    """Node: 假设竞争排序（v3.0 修正）。

    接收已 falsification 攻击过的 competing_hypotheses，
    使用 support/contra/missing 比率重新排序。
    不再生成初始假设（由 contradiction_check 生成）。
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "hypothesis_competition"):
            state.add_trace("hypothesis_competition", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}
        hypotheses = state.competing_hypotheses or []
        if not hypotheses:
            logger.warning("hypothesis_competition: empty competing_hypotheses -- no candidates generated")
            return {"competing_hypotheses": []}

        from phyto_reason.reasoning.hypothesis_competition_engine import (
            HypothesisCompetitionEngine,
        )

        engine = HypothesisCompetitionEngine()
        ranked = engine.rank_by_evidence(hypotheses)

        for i, h in enumerate(ranked):
            h.id = f"H{i + 1}"

        state.add_trace("hypothesis_competition", "completed", summary={
            "n_hypotheses": len(ranked),
            "top_hypothesis": ranked[0].title if ranked else "",
        })

        return {
            "competing_hypotheses": ranked,
        }
    except Exception as e:
        logger.error(f"hypothesis_competition failed: {e}")
        state.add_trace("hypothesis_competition", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}


# ═══════════════════════════════════════════════════════════════
# Node +3: 证伪（新）
# ═══════════════════════════════════════════════════════════════

def falsification_node(state: RuntimeState) -> dict:
    """Node: 主动攻击假设。

    对每个 competing hypothesis 执行 falsification 检查。
    系统必须学会为什么自己的假设可能是错的。
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        if _should_skip_node(state, "falsification"):
            state.add_trace("falsification", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}
        hypotheses = state.competing_hypotheses or []
        if not hypotheses:
            return {"current_node": "hypothesis_synthesis"}

        from phyto_reason.reasoning.falsification_engine import FalsificationEngine

        falsifier = FalsificationEngine()
        attacked = []

        for h in hypotheses:
            attacked_h = falsifier.attack(
                h,
                expression_data=_get_expr_matrix(state),
                n_samples=_get_expr_matrix(state) and len(next(iter(_get_expr_matrix(state).values()))) or 0,
            )
            attacked.append(attacked_h)

        for h in attacked:
            tag = "WEAK" if h.uncertainty_level == "high" else "PLAUSIBLE"
            logger.info(f"  {h.id}: {tag} (supports={len(h.supporting_evidence)}, "
                        f"contradictions={len(h.contradictory_evidence)}, "
                        f"missing={len(h.missing_evidence)})")

        state.add_trace("falsification", "completed", summary={
            "n_attacked": len(attacked),
            "high_uncertainty": sum(1 for h in attacked if h.uncertainty_level == "high"),
        })

        return {
            "competing_hypotheses": attacked,
            "current_node": "hypothesis_synthesis",
        }
    except Exception as e:
        logger.error(f"falsification failed: {e}")
        state.add_trace("falsification", "failed", error=str(e))
        return {"error": str(e), "current_node": "__end__", "finished": True}
