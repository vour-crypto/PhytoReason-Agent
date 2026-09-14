"""
correlation_network.py — Spearman 相关网络节点。

基于 DEGs 和 DAMs 的 Spearman 秩相关构建基因-代谢物关联网络。
使用 NetworkX 进行图分析：hub 检测、模块发现、中心性计算。

升级替代 multiomics_integration_node 中的简易 bicor + 简易模块检测。
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("workflow_nodes.correlation_network")


def _spearman_corr(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """计算 Spearman 秩相关系数及 p-value。

    Args:
        x, y: 等长数值向量

    Returns:
        (rho, p_value)
    """
    n = len(x)
    if n < 3:
        return 0.0, 1.0

    from scipy.stats import spearmanr
    rho, p = spearmanr(x, y)
    # Handle NaN
    if np.isnan(rho):
        return 0.0, 1.0
    return float(rho), float(p)


def _build_correlation_network(
    expr: dict,
    meta_mat: dict,
    gene_ids: list[str],
    meta_ids: list[str],
    common_samples: list[str],
    min_abs_rho: float = 0.8,
) -> dict:
    """构建 Spearman 相关网络。

    Args:
        expr: 表达矩阵
        meta_mat: 代谢物矩阵
        gene_ids: 基因 ID 列表 (top DEGs)
        meta_ids: 代谢物 ID 列表 (top DAMs)
        common_samples: 共同样本名列表
        min_abs_rho: 最小绝对相关系数阈值

    Returns:
        network dict with nodes, edges, metrics
    """
    from phyto_reason.utils.stats_utils import compute_fdr

    # Compute all pairwise Spearman correlations
    edges: list[dict] = []
    all_pvalues: list[float] = []

    for gid in gene_ids:
        gvec = np.array([expr[gid].get(s, float("nan")) for s in common_samples], dtype=float)
        gvec = gvec[~np.isnan(gvec)]
        if len(gvec) < 3 or np.std(gvec) < 1e-10:
            continue

        for mid in meta_ids:
            mvec = np.array([meta_mat[mid].get(s, float("nan")) for s in common_samples], dtype=float)
            mvec = mvec[~np.isnan(mvec)]
            if len(mvec) < 3 or np.std(mvec) < 1e-10:
                continue

            # Use common non-nan indices
            g_vec = np.array([expr[gid].get(s, float("nan")) for s in common_samples], dtype=float)
            m_vec = np.array([meta_mat[mid].get(s, float("nan")) for s in common_samples], dtype=float)
            valid = ~np.isnan(g_vec) & ~np.isnan(m_vec)
            if valid.sum() < 3:
                continue

            rho, p = _spearman_corr(g_vec[valid], m_vec[valid])
            if abs(rho) >= min_abs_rho:
                edges.append({
                    "source": gid,
                    "target": mid,
                    "rho": round(float(rho), 4),
                    "p_value": round(float(p), 6),
                })
                all_pvalues.append(p)

    if not edges:
        return {"n_nodes": len(gene_ids) + len(meta_ids), "n_edges": 0, "reason": "no significant correlations"}

    # FDR correction
    qvalues = compute_fdr(all_pvalues)
    for e, q in zip(edges, qvalues):
        e["q_value"] = round(float(q), 6)

    # Strict default gate: effect size plus nominal p-value and BH-FDR.
    sig_edges = [
        e for e in edges
        if abs(e["rho"]) >= min_abs_rho
        and e["p_value"] < 0.05
        and e["q_value"] < 0.05
    ]

    # Build NetworkX graph
    try:
        import networkx as nx

        G = nx.Graph()

        # Add nodes with type
        for gid in gene_ids:
            G.add_node(gid, node_type="gene")
        for mid in meta_ids:
            G.add_node(mid, node_type="metabolite")

        for e in sig_edges:
            G.add_edge(e["source"], e["target"], weight=abs(e["rho"]), rho=e["rho"], q_value=e["q_value"])

        # ── Network metrics ──────────────────────────────
        # Degree centrality
        degree_cent = nx.degree_centrality(G)

        # Betweenness centrality (may be slow for large graphs)
        try:
            betweenness = nx.betweenness_centrality(G, weight="weight", normalized=True)
        except Exception:
            betweenness = {n: 0.0 for n in G.nodes()}

        # Clustering coefficient
        clustering = nx.clustering(G, weight="weight")

        # Hub detection: top degree + high betweenness
        hubs = []
        for node in G.nodes():
            deg = G.degree(node)
            if deg >= 3:  # At least 3 connections to be considered a hub
                hubs.append({
                    "node": node,
                    "type": G.nodes[node].get("node_type", "unknown"),
                    "degree": deg,
                    "degree_centrality": round(degree_cent.get(node, 0), 4),
                    "betweenness_centrality": round(betweenness.get(node, 0), 4),
                    "clustering_coefficient": round(clustering.get(node, 0), 4),
                    "neighbors": list(G.neighbors(node))[:15],
                })

        hubs.sort(key=lambda h: h["degree_centrality"], reverse=True)

        # Community detection (Louvain if available, else simple connected components)
        modules: list[dict] = []
        try:
            from networkx.algorithms.community import greedy_modularity_communities
            communities = greedy_modularity_communities(G, weight="weight")
            for i, comm in enumerate(communities):
                if len(comm) >= 2:
                    modules.append({
                        "module_id": f"M{i + 1}",
                        "size": len(comm),
                        "nodes": sorted(comm)[:30],  # Truncate for readability
                    })
        except Exception:
            # Fallback: connected components
            components = list(nx.connected_components(G))
            for i, comp in enumerate(components):
                if len(comp) >= 2:
                    modules.append({
                        "module_id": f"C{i + 1}",
                        "size": len(comp),
                        "nodes": sorted(comp)[:30],
                    })

        modules.sort(key=lambda m: m["size"], reverse=True)

        # ── Assemble network data ─────────────────────────
        network_data = {
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "n_sig_edges": len(sig_edges),
            "correlation_method": "spearman",
            "min_abs_rho": min_abs_rho,
            "p_threshold": 0.05,
            "fdr_status": "BH applied; edges require p<0.05 and q<0.05",
            "nodes": [
                {
                    "id": n,
                    "type": G.nodes[n].get("node_type", "unknown"),
                    "degree": G.degree(n),
                    "degree_centrality": round(degree_cent.get(n, 0), 4),
                    "betweenness_centrality": round(betweenness.get(n, 0), 4),
                    "clustering_coefficient": round(clustering.get(n, 0), 4),
                }
                for n in G.nodes()
            ],
            "edges": sig_edges[:200],  # Truncate for transport
            "hubs": hubs[:20],
            "modules": modules[:10],
            "network_density": round(nx.density(G), 4),
            "avg_clustering": round(np.mean(list(clustering.values())), 4) if clustering else 0,
        }

        return network_data

    except ImportError:
        logger.warning("networkx not available, returning edge list only")
        # Fallback without networkx
        return {
            "n_nodes": len(gene_ids) + len(meta_ids),
            "n_edges": len(sig_edges),
            "n_sig_edges": len(sig_edges),
            "correlation_method": "spearman",
            "min_abs_rho": min_abs_rho,
            "p_threshold": 0.05,
            "fdr_status": "BH applied; edges require p<0.05 and q<0.05",
            "edges": sig_edges[:200],
            "hubs": [],
            "modules": [],
            "fallback": "networkx_unavailable",
        }


def correlation_network_node(state: RuntimeState) -> dict:
    """Spearman 相关网络分析节点。

    使用 Spearman 秩相关 (对异常值更鲁棒) 构建基因-代谢物关联网络，
    并用 NetworkX 进行图分析和 hub 检测。

    需要:
      - state.deg_report (提供 top gene IDs)
      - state.dam_report (提供 top metabolite IDs)
      - 表达矩阵和代谢物矩阵

    返回:
      dict with correlation_network_report
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        from phyto_reason.workflows.execution_nodes import _should_skip_node
        if _should_skip_node(state, "correlation_network"):
            state.add_trace("correlation_network", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        expr = None
        meta_mat = None
        if state.planner_state:
            expr = state.planner_state.expression_matrix
            meta_mat = state.planner_state.metabolite_matrix

        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}

        if not expr or not meta_mat:
            state.add_trace("correlation_network", "completed",
                            summary={"reason": "missing expression or metabolite data"})
            return {"correlation_network_report": {"reason": "missing data"}}

        # Get common samples
        expr_samples: set[str] = set()
        for v in expr.values():
            expr_samples.update(v.keys())
            break
        meta_samples: set[str] = set()
        for v in meta_mat.values():
            meta_samples.update(v.keys())
            break
        common_samples = sorted(expr_samples & meta_samples)

        if len(common_samples) < 5:
            state.add_trace("correlation_network", "completed",
                            summary={"reason": f"too few common samples: {len(common_samples)}"})
            return {"correlation_network_report": {"reason": f"insufficient samples: {len(common_samples)}"}}

        # Get DEG gene IDs
        if deg_report.get("top_genes"):
            deg_ids = [g["gene_id"] for g in deg_report["top_genes"][:100]
                       if g["gene_id"] in expr]
        else:
            deg_ids = list(expr.keys())[:100]

        # Get DAM metabolite IDs
        if dam_report.get("top_metabolites"):
            dam_ids = [m["metabolite"] for m in dam_report["top_metabolites"][:50]
                       if m["metabolite"] in meta_mat]
        else:
            dam_ids = list(meta_mat.keys())[:50]

        if not deg_ids or not dam_ids:
            state.add_trace("correlation_network", "completed",
                            summary={"reason": "no valid DEG or DAM IDs"})
            return {"correlation_network_report": {"reason": "no valid IDs"}}

        logger.info(
            "Correlation network: %d genes × %d metabolites × %d samples",
            len(deg_ids), len(dam_ids), len(common_samples),
        )

        # Limit for performance
        max_genes = min(len(deg_ids), 80)
        max_metas = min(len(dam_ids), 40)
        network = _build_correlation_network(
            expr, meta_mat,
            deg_ids[:max_genes], dam_ids[:max_metas],
            common_samples,
            min_abs_rho=0.8,
        )

        # Interpretation
        n_edges = network.get("n_edges", 0)
        n_hubs = len(network.get("hubs", []))
        n_modules = len(network.get("modules", []))

        interpretation = ""
        if n_edges >= 50:
            interpretation = (
                f"Dense correlation network ({n_edges} edges). "
                f"Strong interconnectivity between transcriptome and metabolome. "
                f"{n_hubs} hub nodes identified."
            )
        elif n_edges >= 10:
            interpretation = (
                f"Moderate correlation network ({n_edges} edges). "
                f"Some gene-metabolite associations detected. "
                f"{n_hubs} hub nodes may be key regulatory points."
            )
        elif n_edges > 0:
            interpretation = (
                f"Sparse correlation network ({n_edges} edges). "
                f"Limited gene-metabolite associations suggest weak transcriptional "
                f"coordination or insufficient sample size (n={len(common_samples)})."
            )
        else:
            interpretation = (
                "No significant correlations found. Check data quality or "
                "relax correlation threshold."
            )

        correlation_network_report = {
            **network,
            "n_common_samples": len(common_samples),
            "n_genes_input": len(deg_ids[:max_genes]),
            "n_metabolites_input": len(dam_ids[:max_metas]),
            "interpretation": interpretation,
        }

        # Generate network visualization
        figure_url = None
        try:
            from phyto_reason.visualization.multiomics_figures import plot_correlation_network
            target_title = (state.planner_state.target_metabolite or "metabolite").strip() if state.planner_state else "metabolite"
            figure_url = plot_correlation_network(
                nodes=network.get("nodes", []),
                edges=network.get("edges", []),
                hubs=network.get("hubs", []),
                title=f"Gene-Metabolite Correlation Network — {target_title}",
            )
        except Exception as e:
            logger.debug("Network figure generation skipped: %s", e)
        if figure_url:
            correlation_network_report["figure_url"] = figure_url
            correlation_network_report["figure_markdown"] = f"![Correlation Network]({figure_url})"

        state.correlation_network_report = correlation_network_report

        state.add_trace("correlation_network", "completed", summary={
            "n_edges": n_edges,
            "n_hubs": n_hubs,
            "n_modules": n_modules,
            "density": network.get("network_density", 0),
        })

        logger.info(
            "Correlation network: %d edges, %d hubs, %d modules, density=%.4f",
            n_edges, n_hubs, n_modules, network.get("network_density", 0),
        )

        return {"correlation_network_report": correlation_network_report}

    except Exception as e:
        logger.error("correlation_network_node failed: %s", e, exc_info=True)
        state.add_trace("correlation_network", "failed", error=str(e))
        return {"correlation_network_report": {"error": str(e)}}
