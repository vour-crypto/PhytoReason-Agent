"""
wgcna.py — 加权基因共表达网络分析 (WGCNA) 节点。

自实现 WGCNA 核心算法（参考 Langfelder & Horvath 2008 BMC Bioinformatics）：
  1. Soft-thresholding power 选择 (scale-free topology fit)
  2. Signed adjacency matrix (a_ij = |(1+cor)/2|^β)
  3. Topological Overlap Matrix (TOM)
  4. 层次聚类 + 动态树切割
  5. 模块特征基因 (Module Eigengene) 计算
  6. 模块-性状关联分析

对标高分植物多组学文章中的 WGCNA 标准分析流程。
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("workflow_nodes.wgcna")


# ═══════════════════════════════════════════════════════════════
# Core WGCNA algorithms
# ═══════════════════════════════════════════════════════════════

def _correlation_matrix(data: np.ndarray, method: str = "pearson") -> np.ndarray:
    """Compute gene-gene correlation matrix.

    Args:
        data: (n_genes, n_samples) expression matrix
        method: pearson or spearman

    Returns:
        (n_genes, n_genes) correlation matrix
    """
    if method == "spearman":
        from scipy.stats import rankdata
        data = np.apply_along_axis(rankdata, 1, data)
    return np.corrcoef(data)


def _pick_soft_threshold(
    cor_matrix: np.ndarray,
    candidate_powers: list[int] | None = None,
    min_r2: float = 0.80,
) -> dict:
    """Pick the optimal soft-thresholding power β for scale-free topology.

    Tests each candidate power by computing the R² of the linear fit
    between log(p(k)) and log(k), where p(k) is the degree distribution
    of the weighted network.

    Args:
        cor_matrix: (n_genes, n_genes) correlation matrix
        candidate_powers: powers to test (default 1..20)
        min_r2: minimum R² to accept a power

    Returns:
        {power, r2, mean_connectivity, is_acceptable}
    """
    if candidate_powers is None:
        candidate_powers = list(range(1, 21))

    n_genes = cor_matrix.shape[0]
    results: list[dict] = []

    for power in candidate_powers:
        # Signed adjacency: a_ij = |0.5*(1+cor)|^power
        adj = np.abs(0.5 * (1.0 + cor_matrix)) ** power
        np.fill_diagonal(adj, 0)

        # Connectivity (degree) for each gene
        connectivity = np.sum(adj, axis=1)
        # Remove genes with zero connectivity
        pos_conn = connectivity[connectivity > 0]
        if len(pos_conn) < 10:
            results.append({
                "power": power, "r2": 0.0, "mean_connectivity": 0.0,
                "slope": 0.0, "is_acceptable": False,
            })
            continue

        # Degree distribution
        # Bin connectivity to get frequency
        log_k = np.log10(pos_conn)
        # Use histogram for p(k)
        n_bins = max(10, min(50, len(pos_conn) // 5))
        counts, bin_edges = np.histogram(log_k, bins=n_bins)
        # Use bin centers with count > 0
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        valid = counts > 0
        log_pk = np.log10(counts[valid] / len(pos_conn))
        log_k_bins = bin_centers[valid]

        if len(log_pk) < 3:
            results.append({
                "power": power, "r2": 0.0, "mean_connectivity": float(np.mean(connectivity)),
                "slope": 0.0, "is_acceptable": False,
            })
            continue

        # Linear regression: log(p(k)) ~ log(k)
        # R² = 1 - SS_res / SS_tot
        A = np.vstack([log_k_bins, np.ones_like(log_k_bins)]).T
        slope, intercept = np.linalg.lstsq(A, log_pk, rcond=None)[0]
        predicted = slope * log_k_bins + intercept
        ss_res = np.sum((log_pk - predicted) ** 2)
        ss_tot = np.sum((log_pk - np.mean(log_pk)) ** 2)
        r2 = 1.0 - ss_res / (ss_tot + 1e-15)
        r2 = max(0.0, min(1.0, r2))

        results.append({
            "power": power,
            "r2": round(float(r2), 4),
            "mean_connectivity": round(float(np.mean(connectivity)), 4),
            "slope": round(float(slope), 3),
            "is_acceptable": r2 >= min_r2 and float(slope) < 0,
        })

    # Find best power: first acceptable one with r2 >= min_r2 AND negative slope
    acceptable = [r for r in results if r["is_acceptable"]]
    if acceptable:
        # Pick the one with highest r2 among acceptable
        best = max(acceptable, key=lambda x: x["r2"])
    else:
        # If none acceptable, pick the one with highest r2
        best = max(results, key=lambda x: x["r2"])

    return {
        "best_power": best["power"],
        "best_r2": best["r2"],
        "is_satisfactory": best["is_acceptable"],
        "all_results": results,
    }


def _compute_tom(adj: np.ndarray) -> np.ndarray:
    """Compute Topological Overlap Matrix (TOM).

    TOM_ij = (Σ_k a_ik * a_kj + a_ij) / (min(k_i, k_j) + 1 - a_ij)

    where k_i = Σ_u a_iu (connectivity of node i).

    Args:
        adj: (n, n) adjacency matrix

    Returns:
        (n, n) TOM matrix
    """
    n = adj.shape[0]
    # Connectivity vector
    k = np.sum(adj, axis=1)

    # Σ_k a_ik * a_kj for all i,j: just A @ A
    # This is O(n³) so for large n we need to be efficient
    # For n ≤ 2000 this is fine
    a_sq = adj @ adj  # (n, n)

    tom = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                tom[i, j] = 1.0
            else:
                num = a_sq[i, j] + adj[i, j]
                den = min(k[i], k[j]) + 1.0 - adj[i, j]
                if den > 0:
                    tom[i, j] = num / den
                    tom[j, i] = tom[i, j]

    return tom


def _dynamic_tree_cut(
    dissimilarity: np.ndarray,
    min_module_size: int = 10,
    deep_split: int = 2,
) -> np.ndarray:
    """Simplified dynamic tree cut for module detection.

    Uses hierarchical clustering + adaptive height cutting.

    Args:
        dissimilarity: (n, n) dissimilarity matrix (1 - TOM)
        min_module_size: minimum genes per module
        deep_split: sensitivity (0-4, higher = more modules)

    Returns:
        (n,) array of module labels (0 = grey/unassigned)
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    n = dissimilarity.shape[0]

    # Condensed distance matrix
    dist_vec = squareform(dissimilarity, checks=False)

    # Hierarchical clustering (average linkage)
    Z = linkage(dist_vec, method="average")

    # Dynamic height threshold based on deep_split
    # deep_split: 0 = few modules, 4 = many modules
    heights = Z[:, 2]
    if len(heights) < 2:
        return np.zeros(n, dtype=int)

    # Adaptive cutoff based on percentile of merge heights
    percentile_cut = max(55, 95 - deep_split * 12)
    cut_height = np.percentile(heights, percentile_cut)

    labels = fcluster(Z, t=cut_height, criterion="distance")

    # Relabel as 0, 1, 2... (0 = grey)
    unique_labels = sorted(set(labels))
    label_map = {old: i + 1 for i, old in enumerate(unique_labels)}
    relabeled = np.array([label_map.get(l, 0) for l in labels], dtype=int)

    # Merge small modules into grey
    for label in list(label_map.values()):
        mask = relabeled == label
        if mask.sum() < min_module_size:
            relabeled[mask] = 0

    # Re-number to be contiguous (skip 0)
    kept = sorted(set(relabeled) - {0})
    final_map = {old: i + 1 for i, old in enumerate(kept)}
    final_map[0] = 0
    relabeled = np.array([final_map[l] for l in relabeled], dtype=int)

    return relabeled


def _module_eigengenes(
    expr: np.ndarray,
    module_labels: np.ndarray,
) -> dict[int, np.ndarray]:
    """Compute module eigengene (first principal component) for each module.

    Args:
        expr: (n_genes, n_samples) expression matrix
        module_labels: (n_genes,) module assignments (0 = grey)

    Returns:
        dict: {module_id: eigengene_vector (n_samples,)}
    """
    from numpy.linalg import svd

    eigengenes: dict[int, np.ndarray] = {}
    n_samples = expr.shape[1]

    for module_id in sorted(set(module_labels)):
        if module_id == 0:
            continue

        mask = module_labels == module_id
        module_expr = expr[mask, :]  # (n_mod_genes, n_samples)

        if module_expr.shape[0] < 2:
            # Single gene module: eigengene = expression
            eigengenes[module_id] = module_expr[0, :]
            continue

        # Center the data
        centered = module_expr - module_expr.mean(axis=1, keepdims=True)

        # SVD for first PC
        try:
            U, s, Vt = svd(centered, full_matrices=False)
            # First PC signs: choose direction to have positive correlation with most genes
            pc1 = Vt[0, :]  # (n_samples,)
            # Ensure positive mean correlation
            corrs = np.array([np.corrcoef(pc1, row)[0, 1] for row in centered])
            if np.mean(corrs) < 0:
                pc1 = -pc1
            eigengenes[module_id] = pc1
        except Exception:
            # Fallback: mean expression
            eigengenes[module_id] = module_expr.mean(axis=0)

    return eigengenes


def _module_trait_correlation(
    eigengenes: dict[int, np.ndarray],
    traits: np.ndarray,
    trait_names: list[str],
    module_sizes: dict[int, int],
) -> list[dict]:
    """Compute module eigengene - trait correlations.

    Args:
        eigengenes: {module_id: eigengene (n_samples,)}
        traits: (n_samples, n_traits) trait matrix
        trait_names: names of each trait column
        module_sizes: {module_id: n_genes}

    Returns:
        [{module_id, module_size, trait, correlation, p_value}, ...]
    """
    results = []
    for module_id, me in eigengenes.items():
        for j, trait_name in enumerate(trait_names):
            t = traits[:, j]
            # Remove any NaN samples
            valid = ~np.isnan(me) & ~np.isnan(t)
            if valid.sum() < 3:
                continue
            r = np.corrcoef(me[valid], t[valid])[0, 1]
            if np.isnan(r):
                continue
            # P-value via t-distribution
            n = valid.sum()
            t_stat = r * np.sqrt((n - 2) / (1 - r * r + 1e-15))
            from scipy.stats import t as t_dist
            p = 2.0 * (1.0 - t_dist.cdf(abs(t_stat), df=n - 2))

            results.append({
                "module_id": module_id,
                "module_size": module_sizes.get(module_id, 0),
                "trait": trait_name,
                "correlation": round(float(r), 4),
                "p_value": round(float(p), 6),
            })

    # FDR correction
    if results:
        from phyto_reason.utils.stats_utils import compute_fdr
        pvalues = [r["p_value"] for r in results]
        qvalues = compute_fdr(pvalues)
        for r, q in zip(results, qvalues):
            r["q_value"] = round(float(q), 6)

    return results


def _run_wgcna(
    expr: dict,
    gene_ids: list[str],
    common_samples: list[str],
    meta_mat: dict | None = None,
    meta_ids: list[str] | None = None,
    min_module_size: int = 10,
    deep_split: int = 2,
) -> dict:
    """Run full WGCNA pipeline.

    Args:
        expr: expression matrix {gene_id: {sample: value}}
        gene_ids: gene IDs to include (top DEGs)
        common_samples: common sample names
        meta_mat: metabolite matrix for trait correlation (optional)
        meta_ids: metabolite IDs for traits (optional)
        min_module_size: min genes per module
        deep_split: tree cutting sensitivity

    Returns:
        WGCNA report dict
    """
    n_genes = len(gene_ids)
    n_samples = len(common_samples)

    if n_genes < 20:
        return {"reason": f"too few genes ({n_genes}) for WGCNA", "n_genes": n_genes}

    # Build expression matrix
    expr_matrix = np.zeros((n_genes, n_samples))
    valid_genes: list[str] = []
    for i, gid in enumerate(gene_ids):
        if gid in expr:
            vals = [expr[gid].get(s, np.nan) for s in common_samples]
            if np.sum(~np.isnan(vals)) >= n_samples * 0.8:
                # Fill missing with gene mean
                arr = np.array(vals)
                arr[np.isnan(arr)] = np.nanmean(arr)
                expr_matrix[len(valid_genes), :] = arr
                valid_genes.append(gid)

    if len(valid_genes) < 20:
        return {"reason": f"too few valid genes ({len(valid_genes)})", "n_genes": len(valid_genes)}

    expr_matrix = expr_matrix[:len(valid_genes), :]
    n_genes = len(valid_genes)

    # Log-transform (WGCNA standard)；count 数据先做文库大小归一化
    positive_vals = expr_matrix[expr_matrix > 0]
    count_like = positive_vals.size > 0 and float(np.median(positive_vals)) > 50.0
    if count_like:
        from phyto_reason.utils.stats_utils import compute_size_factors_array
        size_factors = compute_size_factors_array(expr_matrix)
        expr_log = np.log2(expr_matrix / size_factors[None, :] + 1.0)
        normalization = "median_of_ratios_size_factors + log2(count/size_factor + 1)"
    elif positive_vals.size > 0 and float(positive_vals.max()) <= 30.0:
        # 已是对数标度（如公共数据 log2 值）——跳过二次 log2
        expr_log = expr_matrix
        normalization = "already_log_scale_pass_through"
    else:
        expr_log = np.log2(expr_matrix + 1.0)
        normalization = "log2(value + 1)"
    logger.info("WGCNA normalization: %s", normalization)

    # Step 1: Gene-gene correlation
    logger.info("WGCNA step 1: computing %d x %d correlation matrix", n_genes, n_genes)
    cor_mat = _correlation_matrix(expr_log, method="pearson")

    # Step 2: Soft threshold selection
    logger.info("WGCNA step 2: picking soft threshold power")
    threshold_result = _pick_soft_threshold(cor_mat, min_r2=0.70)
    power = threshold_result["best_power"]
    is_satisfactory = threshold_result["is_satisfactory"]

    # Step 3: Adjacency matrix
    logger.info("WGCNA step 3: computing adjacency (power=%d)", power)
    adj = np.abs(0.5 * (1.0 + cor_mat)) ** power
    np.fill_diagonal(adj, 0)

    # Step 4: TOM
    logger.info("WGCNA step 4: computing TOM (%d genes)", n_genes)
    tom = _compute_tom(adj)
    dissimilarity = 1.0 - tom

    # Step 5: Hierarchical clustering + module detection
    logger.info("WGCNA step 5: module detection")
    module_labels = _dynamic_tree_cut(dissimilarity, min_module_size=min_module_size, deep_split=deep_split)

    n_modules = len(set(module_labels) - {0})

    # Step 6: Module eigengenes
    logger.info("WGCNA step 6: computing module eigengenes")
    eigengenes = _module_eigengenes(expr_log, module_labels)

    # Module sizes
    module_sizes: dict[int, int] = {}
    for label in set(module_labels):
        if label != 0:
            module_sizes[label] = int(np.sum(module_labels == label))
    grey_count = int(np.sum(module_labels == 0))

    # Module gene lists
    module_genes: dict[int, list[str]] = {}
    for label in sorted(module_sizes.keys()):
        mask = module_labels == label
        module_genes[label] = [valid_genes[i] for i in range(n_genes) if mask[i]]

    # Step 7: Module-trait correlations
    trait_correlations: list[dict] = []
    if meta_mat and meta_ids and common_samples:
        # Build trait matrix from top metabolite levels
        n_traits = min(len(meta_ids), 10)
        traits = np.zeros((n_samples, n_traits))
        trait_names = []
        for j, mid in enumerate(meta_ids[:n_traits]):
            if mid in meta_mat:
                t_vals = [meta_mat[mid].get(s, np.nan) for s in common_samples]
                arr = np.array(t_vals)
                arr[np.isnan(arr)] = np.nanmean(arr)
                traits[:, j] = arr
            else:
                traits[:, j] = 0.0
            trait_names.append(mid)

        trait_correlations = _module_trait_correlation(
            eigengenes, traits, trait_names, module_sizes,
        )

    # Step 8: Hub genes (top intramodular connectivity)
    hub_genes: dict[int, list[dict]] = {}
    kIM = adj  # Intramodular connectivity (use adjacency as proxy)
    for label in sorted(module_sizes.keys()):
        mask = module_labels == label
        if mask.sum() == 0:
            continue
        # Mean connectivity within module
        mod_kIM = np.mean(kIM[mask, :][:, mask], axis=1)
        top_idx = np.argsort(mod_kIM)[-10:][::-1]  # Top 10
        hub_genes[label] = [
            {"gene_id": valid_genes[i], "kIM": round(float(mod_kIM[idx]), 4)}
            for idx, i in enumerate(top_idx)
        ]

    # Build report
    wgcna_report = {
        "n_genes_input": len(gene_ids),
        "n_genes_analyzed": n_genes,
        "n_samples": n_samples,
        "soft_power": power,
        "scale_free_r2": threshold_result["best_r2"],
        "scale_free_satisfactory": is_satisfactory,
        "n_modules": n_modules,
        "n_grey_genes": grey_count,
        "min_module_size": min_module_size,
        "module_sizes": {str(k): v for k, v in sorted(module_sizes.items(), key=lambda x: -x[1])},
        "module_genes": {str(k): v[:30] for k, v in module_genes.items()},
        "module_eigengenes": {str(k): v[:10].tolist() for k, v in list(eigengenes.items())[:10]},
        "hub_genes": {str(k): v for k, v in list(hub_genes.items())[:10]},
        "trait_correlations": trait_correlations[:50],
        "power_selection_results": threshold_result["all_results"][:10],
        "normalization": normalization,
        "warnings": [],
    }

    if not is_satisfactory:
        wgcna_report["warnings"].append(
            f"Scale-free topology fit R²={threshold_result['best_r2']:.3f} < 0.80. "
            f"Network may not be truly scale-free. Consider using more samples (n={n_samples}) "
            f"or different gene selection."
        )

    if n_modules == 0:
        wgcna_report["warnings"].append(
            "No modules detected. Try reducing min_module_size or increasing deep_split."
        )

    return wgcna_report


# ═══════════════════════════════════════════════════════════════
# Node
# ═══════════════════════════════════════════════════════════════

def wgcna_node(state: RuntimeState) -> dict:
    """WGCNA 共表达网络分析节点。

    执行完整的 WGCNA 流程:
      1. 软阈值选择 (scale-free topology)
      2. 加权邻接矩阵 (signed, power β)
      3. Topological Overlap Matrix (TOM)
      4. 层次聚类 + 动态模块检测
      5. 模块特征基因 (PC1) 计算
      6. 模块-代谢物性状关联
      7. Hub 基因鉴定 (intramodular connectivity)

    Returns:
        dict with wgcna_report
    """
    try:
        from phyto_reason.workflows.execution_nodes import _should_skip_node
        if _should_skip_node(state, "wgcna"):
            state.add_trace("wgcna", "skipped", summary={"reason": "user skipped this step"})
            return {}

        expr = None
        meta_mat = None
        if state.planner_state:
            expr = state.planner_state.expression_matrix
            meta_mat = state.planner_state.metabolite_matrix

        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}

        if not expr:
            state.add_trace("wgcna", "completed", summary={"reason": "no expression data"})
            return {"wgcna_report": {"reason": "no expression data"}}

        # Get sample list
        samples = set()
        for v in expr.values():
            samples.update(v.keys())
            break
        common_samples = sorted(samples)

        if len(common_samples) < 8:
            state.add_trace("wgcna", "completed",
                            summary={"reason": f"too few samples for WGCNA: {len(common_samples)} (need ≥8)"})
            return {"wgcna_report": {
                "reason": f"WGCNA requires ≥8 samples for reliable modules (got {len(common_samples)})",
                "n_samples": len(common_samples),
            }}

        # Select genes: FDR-sig DEGs first, then top effect-size
        gene_ids: list[str] = []
        if deg_report.get("top_genes"):
            deg_ids = [g["gene_id"] for g in deg_report["top_genes"]
                       if g["gene_id"] in expr]
            # Prioritize FDR-sig
            fdr_sig = [g["gene_id"] for g in deg_report["top_genes"]
                       if g.get("is_fdr_sig") and g["gene_id"] in expr]
            gene_ids = fdr_sig[:500] + [g for g in deg_ids if g not in fdr_sig][:500]
            # DEG reports intentionally expose a compact top list. WGCNA needs
            # a wider variable set to form modules, so complete it from the
            # expression matrix while retaining DEG priority.
            if len(gene_ids) < 100:
                gene_ids.extend(g for g in expr if g not in gene_ids)
        else:
            gene_ids = list(expr.keys())[:1000]

        # Limit to ~500 for performance (TOM is O(n³))
        gene_ids = gene_ids[:500]

        # Get metabolite trait IDs
        meta_ids: list[str] = []
        if dam_report.get("top_metabolites"):
            meta_ids = [m["metabolite"] for m in dam_report["top_metabolites"][:15]
                        if m["metabolite"] in (meta_mat or {})]

        logger.info("WGCNA: %d genes x %d samples, %d traits",
                    len(gene_ids), len(common_samples), len(meta_ids))

        wgcna_report = _run_wgcna(
            expr=expr,
            gene_ids=gene_ids,
            common_samples=common_samples,
            meta_mat=meta_mat,
            meta_ids=meta_ids,
            min_module_size=10,
            deep_split=2,
        )

        if "reason" in wgcna_report and "n_modules" not in wgcna_report:
            state.add_trace("wgcna", "completed", summary=wgcna_report)
            state.wgcna_report = wgcna_report
            return {"wgcna_report": wgcna_report}

        n_modules = wgcna_report.get("n_modules", 0)
        n_hub = len(wgcna_report.get("hub_genes", {}))
        n_sig_trait = sum(
            1 for t in wgcna_report.get("trait_correlations", [])
            if t.get("q_value", 1.0) < 0.05
        )

        # Interpretation
        if n_modules >= 3 and n_sig_trait > 0:
            interpretation = (
                f"WGCNA identified {n_modules} co-expression modules. "
                f"{n_sig_trait} significant module-trait associations found, "
                f"suggesting coordinated transcriptional programs linked to metabolite variation."
            )
        elif n_modules >= 3:
            interpretation = (
                f"WGCNA identified {n_modules} co-expression modules. "
                f"No significant module-trait associations at FDR<0.05. "
                f"Module membership may still be useful for functional annotation."
            )
        else:
            interpretation = (
                f"WGCNA detected {n_modules} modules. "
                f"Limited modular structure suggests weak co-expression patterns. "
                f"More samples or different gene selection may improve results."
            )

        wgcna_report["interpretation"] = interpretation

        # Generate WGCNA figures
        figure_urls: list[str] = []
        figure_mds: list[str] = []
        try:
            from phyto_reason.visualization.multiomics_figures import (
                plot_wgcna_scale_free, plot_wgcna_module_sizes,
            )
            # Scale-free topology plot
            power_results = wgcna_report.get("power_selection_results", [])
            sft_url = plot_wgcna_scale_free(
                power_results,
                best_power=wgcna_report.get("soft_power"),
                title="WGCNA Scale-Free Topology Fit",
            )
            if sft_url:
                figure_urls.append(sft_url)
                figure_mds.append(f"![WGCNA Scale-Free Topology]({sft_url})")

            # Module sizes bar chart
            mod_sizes = {str(k): int(v) for k, v in wgcna_report.get("module_sizes", {}).items()}
            if mod_sizes:
                mod_url = plot_wgcna_module_sizes(
                    mod_sizes,
                    grey_count=wgcna_report.get("n_grey_genes", 0),
                    title="WGCNA Module Sizes",
                )
                if mod_url:
                    figure_urls.append(mod_url)
                    figure_mds.append(f"![WGCNA Modules]({mod_url})")
        except Exception as e:
            logger.debug("WGCNA figure generation skipped: %s", e)

        if figure_urls:
            wgcna_report["figure_urls"] = figure_urls
            wgcna_report["figure_markdown"] = "\n".join(figure_mds)

        state.wgcna_report = wgcna_report

        state.add_trace("wgcna", "completed", summary={
            "n_genes": wgcna_report.get("n_genes_analyzed", 0),
            "n_modules": n_modules,
            "n_sig_trait_correlations": n_sig_trait,
            "soft_power": wgcna_report.get("soft_power", 0),
            "scale_free_r2": wgcna_report.get("scale_free_r2", 0),
        })

        logger.info("WGCNA: %d modules, %d sig trait corrs, power=%d",
                    n_modules, n_sig_trait, wgcna_report.get("soft_power", 0))

        return {"wgcna_report": wgcna_report}

    except Exception as e:
        logger.error("wgcna_node failed: %s", e, exc_info=True)
        state.add_trace("wgcna", "failed", error=str(e))
        return {"wgcna_report": {"error": str(e)}}
