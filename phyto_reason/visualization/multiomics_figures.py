from __future__ import annotations

import numpy as np

from ._common import PALETTE, gradient_cmap, new_figure, save_figure


def _save(fig, prefix: str) -> str:
    return save_figure(fig, prefix)


def _sample_matrix(matrix: dict):
    names = list(next(iter(matrix.values())).keys()) if matrix else []
    rows = [[float(v.get(s, 0)) for s in names] for v in matrix.values()]
    return names, np.asarray(rows, dtype=float)


def plot_qc_correlation_heatmap(expr: dict, groups=None) -> str:
    names, values = _sample_matrix(expr)
    if values.size == 0:
        return ""
    corr = np.corrcoef(values.T)
    fig, ax = new_figure("Sample Correlation", width=6.5, height=5.5)
    image = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdYlBu_r")
    ax.set_xticks(range(len(names)), names, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(names)), names, fontsize=8)
    fig.colorbar(image, ax=ax, label="Pearson r", shrink=0.8)
    return _save(fig, "qc_correlation")


def plot_qc_pca(expr: dict, groups=None) -> str:
    names, values = _sample_matrix(expr)
    if values.shape[1] < 2:
        return ""
    centered = values - values.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    coords = (centered.T @ (centered @ vt[:2].T)) if vt.shape[0] >= 2 else centered.T
    coords = coords[:, :2]
    fig, ax = new_figure("PCA of Samples")
    labels = [groups.get(n, "all") if isinstance(groups, dict) else "all" for n in names]
    unique = list(dict.fromkeys(labels))
    for i, label in enumerate(unique):
        idx = [j for j, x in enumerate(labels) if x == label]
        ax.scatter(coords[idx, 0], coords[idx, 1], s=42, label=label, color=[PALETTE["blue"], PALETTE["teal"], PALETTE["orange"]][i % 3])
    for i, name in enumerate(names):
        ax.annotate(name, (coords[i, 0], coords[i, 1]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    if len(unique) > 1: ax.legend(frameon=False)
    return _save(fig, "qc_pca")


def plot_qc_upset(sample_gene_sets: dict) -> str:
    labels = list(sample_gene_sets)
    counts = [len(sample_gene_sets[k]) for k in labels]
    fig, ax = new_figure("Expressed Features by Group")
    ax.bar(labels, counts, color=PALETTE["teal"])
    ax.set_ylabel("number of features")
    ax.tick_params(axis="x", rotation=35)
    return _save(fig, "qc_features")


def plot_wgcna_scale_free(power_results, best_power=None, title="WGCNA Scale-Free Topology Fit") -> str:
    if not power_results:
        return ""
    fig, ax = new_figure(title)
    x = [r.get("power", r.get("soft_power", 0)) for r in power_results]
    y = [r.get("r2", r.get("r_squared", r.get("scale_free_r2", 0))) for r in power_results]
    ax.plot(x, y, marker="o", color=PALETTE["blue"])
    if best_power is not None: ax.axvline(best_power, color=PALETTE["orange"], linestyle="--")
    ax.set_xlabel("soft-threshold power"); ax.set_ylabel("scale-free fit R²")
    return _save(fig, "wgcna_scale_free")


def plot_wgcna_module_sizes(module_sizes: dict, grey_count=0, title="WGCNA Module Sizes") -> str:
    if not module_sizes:
        return ""
    fig, ax = new_figure(title)
    labels, values = list(module_sizes), list(module_sizes.values())
    ax.bar(labels, values, color=PALETTE["teal"])
    if grey_count: ax.bar(["grey"], [grey_count], color=PALETTE["muted"])
    ax.set_ylabel("genes"); ax.tick_params(axis="x", rotation=40)
    return _save(fig, "wgcna_modules")


def plot_enrichment_bars(results, title="Pathway Enrichment") -> str:
    """论文气泡样式富集图：x=-log10(q)，点大小=富集倍数，颜色=显著性渐变。"""
    rows = results[:12] if isinstance(results, list) else []
    if not rows: return ""
    labels = [str(r.get("pathway_name", r.get("name", "pathway")))[:42] for r in rows]
    q = [max(float(r.get("q_value", r.get("p_value", 1)) or 1), 1e-300) for r in rows]
    x = [-np.log10(v) for v in q]
    sizes = []
    for r in rows:
        try:
            fold = float(r.get("fold_enrichment", r.get("fold", 1)) or 1)
        except (TypeError, ValueError):
            fold = 1.0
        sizes.append(min(max(fold, 1.0), 30.0) * 12 + 30)

    fig, ax = new_figure(title, width=8, height=max(4.5, len(rows) * 0.42))
    scatter = ax.scatter(x, range(len(labels)), s=sizes,
                         c=x, cmap=gradient_cmap(), edgecolors="black",
                         linewidths=0.4, alpha=0.9, zorder=3)
    ax.set_yticks(range(len(labels)), labels, fontsize=10)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.5)
    ax.set_xlabel("-log10 adjusted p-value")
    colorbar = fig.colorbar(scatter, ax=ax, shrink=0.8)
    colorbar.set_label("-log10(q)", fontsize=10)
    # 富集倍数尺寸图例说明
    ax.text(0.99, 0.02, "bubble size = fold enrichment",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="grey")
    return _save(fig, "pathway_enrichment")


def plot_dam_pairwise_heatmap(pairwise_results, title="DAM Pairwise Tissue Comparisons") -> str:
    """Summarize tissue-pair log2FC values for ANOVA-supported metabolites."""
    rows = list(pairwise_results or [])
    if not rows:
        return ""
    metabolites = list(dict.fromkeys(str(r.get("metabolite", "")) for r in rows))[:40]
    pairs = list(dict.fromkeys(
        f"{r.get('group_a', '')} vs {r.get('group_b', '')}" for r in rows
    ))
    if not metabolites or not pairs:
        return ""
    values = np.full((len(metabolites), len(pairs)), np.nan, dtype=float)
    row_index = {name: i for i, name in enumerate(metabolites)}
    col_index = {name: i for i, name in enumerate(pairs)}
    for r in rows:
        met = str(r.get("metabolite", ""))
        pair = f"{r.get('group_a', '')} vs {r.get('group_b', '')}"
        if met in row_index and pair in col_index:
            values[row_index[met], col_index[pair]] = float(r.get("log2fc", 0) or 0)
    # Cluster rows/columns when scipy is available; NaN means an untested
    # combination and is treated as zero only for ordering.
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        fill = np.nan_to_num(values, nan=0.0)
        if fill.shape[0] > 1:
            values = values[leaves_list(linkage(fill, method="average"))]
            metabolites = [metabolites[i] for i in leaves_list(linkage(fill, method="average"))]
        if fill.shape[1] > 1:
            order = leaves_list(linkage(fill.T, method="average"))
            values = values[:, order]
            pairs = [pairs[i] for i in order]
    except Exception:
        pass
    fig, ax = new_figure(title, width=max(7.0, len(pairs) * 0.7), height=max(4.5, len(metabolites) * 0.22))
    max_abs = float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 1.0
    max_abs = max(max_abs, 1e-6)
    image = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-max_abs, vmax=max_abs)
    ax.set_xticks(range(len(pairs)), pairs, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(metabolites)), metabolites, fontsize=7)
    fig.colorbar(image, ax=ax, label="log2FC (group B - group A)", shrink=0.8)
    return _save(fig, "dam_pairwise_heatmap")


def plot_correlation_network(nodes, edges, hubs=None, title="Correlation Network") -> str:
    if not nodes: return ""
    fig, ax = new_figure(title, width=8, height=6)
    n = len(nodes); angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = {str(node.get("id", node.get("name", i))): (np.cos(a), np.sin(a)) for i, (node, a) in enumerate(zip(nodes, angles))}
    for edge in edges or []:
        a, b = str(edge.get("source", "")), str(edge.get("target", ""))
        if a in pos and b in pos: ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#CBD5E1", linewidth=0.8, zorder=1)
    hub_set = {str(x) for x in (hubs or [])}
    for key, (x, y) in pos.items():
        ax.scatter(x, y, s=75 if key in hub_set else 38, color=PALETTE["orange"] if key in hub_set else PALETTE["blue"], zorder=2)
        ax.text(x, y, key[:16], fontsize=7, ha="center", va="bottom")
    ax.set_axis_off()
    return _save(fig, "correlation_network")


def plot_o2pls_scores(scores, title="O2PLS Joint Component Scores") -> str:
    if not scores: return ""
    arr = np.asarray(scores, dtype=float)
    if arr.ndim == 1: arr = arr.reshape(-1, 1)
    fig, ax = new_figure(title)
    x = arr[:, 0]; y = arr[:, 1] if arr.shape[1] > 1 else np.zeros(len(arr))
    ax.scatter(x, y, color=PALETTE["blue"], s=38)
    ax.set_xlabel("component 1"); ax.set_ylabel("component 2")
    return _save(fig, "o2pls_scores")


def plot_o2pls_loadings(top_x, top_y, title="O2PLS Variable Loadings") -> str:
    rows = []
    for item in list(top_x or [])[:8]: rows.append((str(item.get("name", item.get("gene_id", "X"))), float(item.get("loading", item.get("score", 0))), "X"))
    for item in list(top_y or [])[:8]: rows.append((str(item.get("name", item.get("metabolite", "Y"))), float(item.get("loading", item.get("score", 0))), "Y"))
    if not rows: return ""
    fig, ax = new_figure(title, width=8, height=max(4.5, len(rows) * 0.28))
    labels = [r[0][:30] for r in rows]; vals = [r[1] for r in rows]; colors = [PALETTE["blue"] if r[2] == "X" else PALETTE["orange"] for r in rows]
    ax.barh(labels[::-1], vals[::-1], color=colors[::-1]); ax.set_xlabel("loading")
    return _save(fig, "o2pls_loadings")
