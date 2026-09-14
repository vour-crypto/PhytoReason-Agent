"""
quadrant_plot.py — 四象限/九象限图节点。

对 DEG-DAM 按 log2FC 方向分类，生成四象限或九象限散点图。
对标真实植物多组学文章的标准可视化。

四象限模式:
  Q1 (+,+): 基因上调 + 代谢物上调 → 协同上调
  Q2 (-,+): 基因下调 + 代谢物上调 → 反向变化
  Q3 (-,-): 基因下调 + 代谢物下调 → 协同下调
  Q4 (+,-): 基因上调 + 代谢物下调 → 反向变化

九象限模式:
  加入 log2FC 阈值 (±T)，产生中心"无变化"区域
"""

from __future__ import annotations

import base64
import io
import logging
from typing import Any

import numpy as np

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("workflow_nodes.quadrant_plot")


def _classify_quadrant(gene_lfc: float, meta_lfc: float, threshold: float = 0.0) -> str:
    """分类 gene-metabolite 对的象限。

    Args:
        gene_lfc: 基因 log2 fold-change
        meta_lfc: 代谢物 log2 fold-change
        threshold: log2FC 阈值 (0 = 四象限, >0 = 九象限)

    Returns:
        象限标签: Q1-Q4 (四象限) 或 Q1-Q9 (九象限)
    """
    if threshold > 0:
        # Nine-quadrant classification
        if gene_lfc > threshold:
            if meta_lfc > threshold:
                return "Q1"
            elif meta_lfc < -threshold:
                return "Q4"
            else:
                return "Q5"
        elif gene_lfc < -threshold:
            if meta_lfc > threshold:
                return "Q2"
            elif meta_lfc < -threshold:
                return "Q3"
            else:
                return "Q7"
        else:
            if meta_lfc > threshold:
                return "Q8"
            elif meta_lfc < -threshold:
                return "Q6"
            else:
                return "Q9"
    else:
        # Four-quadrant classification
        if gene_lfc >= 0:
            if meta_lfc >= 0:
                return "Q1"  # gene up, meta up → coordinated up
            else:
                return "Q4"  # gene up, meta down → inverse
        else:
            if meta_lfc >= 0:
                return "Q2"  # gene down, meta up → inverse
            else:
                return "Q3"  # gene down, meta down → coordinated down


def _generate_quadrant_plot(
    pairs: list[dict],
    threshold: float = 0.0,
    gene_label: str = "Gene log2FC",
    meta_label: str = "Metabolite log2FC",
    title: str = "",
) -> str | None:
    """生成象限散点图并返回 base64 PNG 字符串。

    Args:
        pairs: [{gene_id, metabolite, gene_log2fc, meta_log2fc, quadrant}, ...]
        threshold: log2FC 阈值
        gene_label, meta_label: 轴标签
        title: 图标题

    Returns:
        base64-encoded PNG string, or None if matplotlib unavailable
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if not pairs:
            return None

        gene_lfcs = [p["gene_log2fc"] for p in pairs]
        meta_lfcs = [p["meta_log2fc"] for p in pairs]
        quadrants = [p.get("quadrant", "Q9") for p in pairs]

        # Color scheme
        quadrant_colors = {
            "Q1": "#e74c3c",   # red — coordinated up
            "Q2": "#e67e22",   # orange
            "Q3": "#27ae60",   # green — coordinated down
            "Q4": "#3498db",   # blue
            "Q5": "#95a5a6",   # gray — mid
            "Q6": "#95a5a6",
            "Q7": "#95a5a6",
            "Q8": "#95a5a6",
            "Q9": "#bdc3c7",   # light gray — no change
        }

        fig, ax = plt.subplots(figsize=(8, 7))

        # Plot each quadrant with its own color
        for q, color in quadrant_colors.items():
            mask = [quad == q for quad in quadrants]
            if any(mask):
                qx = [gene_lfcs[i] for i, m in enumerate(mask) if m]
                qy = [meta_lfcs[i] for i, m in enumerate(mask) if m]
                ax.scatter(qx, qy, c=color, label=q, alpha=0.7, s=40, edgecolors="white", linewidth=0.5)

        # Add reference lines
        max_abs = max(
            max(abs(x) for x in gene_lfcs) if gene_lfcs else 1,
            max(abs(y) for y in meta_lfcs) if meta_lfcs else 1,
        ) * 1.15
        ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
        ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

        if threshold > 0:
            ax.axhline(y=threshold, color="gray", linestyle=":", linewidth=0.5, alpha=0.3)
            ax.axhline(y=-threshold, color="gray", linestyle=":", linewidth=0.5, alpha=0.3)
            ax.axvline(x=threshold, color="gray", linestyle=":", linewidth=0.5, alpha=0.3)
            ax.axvline(x=-threshold, color="gray", linestyle=":", linewidth=0.5, alpha=0.3)

        ax.set_xlim(-max_abs, max_abs)
        ax.set_ylim(-max_abs, max_abs)
        ax.set_xlabel(gene_label, fontsize=12)
        ax.set_ylabel(meta_label, fontsize=12)
        ax.set_title(title or "DEG-DAM Quadrant Plot", fontsize=13, fontweight="bold")

        # Quadrant label annotations
        if threshold == 0:
            hw = max_abs * 0.5
            for pos, label in [((hw, hw), "Q1\nSyn↑"), ((-hw, hw), "Q2\nInv"),
                               ((-hw, -hw), "Q3\nSyn↓"), ((hw, -hw), "Q4\nInv")]:
                ax.annotate(label, xy=pos, fontsize=10, fontweight="bold",
                            ha="center", va="center", alpha=0.3)

        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
        ax.grid(True, alpha=0.2)
        plt.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")
        return f"data:image/png;base64,{img_base64}"

    except ImportError:
        logger.warning("matplotlib not available, skipping quadrant plot generation")
        return None
    except Exception as e:
        logger.error("Failed to generate quadrant plot: %s", e)
        return None


def _infer_log2fc(deg_report: dict, gene_id: str) -> float:
    """从 DEG report 推断基因的 log2FC。"""
    for g in deg_report.get("top_genes", []):
        if g.get("gene_id") == gene_id:
            # top_genes 使用的是 max_log2fc
            lfc = g.get("max_log2fc", 0)
            if lfc != 0:
                return float(lfc)
    return 0.0


def quadrant_plot_node(state: RuntimeState) -> dict:
    """四象限/九象限图节点。

    使用 DEG 和 DAM 的 log2FC 数据，将 gene-metabolite 对
    分类到象限中。可选生成 matplotlib 图。

    需要:
      - state.deg_report (含 top_genes with log2FC)
      - state.dam_report (含 top_metabolites with log2FC)
      - state.multiomics_report (含 gene-metabolite 相关对)

    返回:
      dict with quadrant_plot_report
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        from phyto_reason.workflows.execution_nodes import _should_skip_node
        if _should_skip_node(state, "quadrant_plot"):
            state.add_trace("quadrant_plot", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}
        multiomics_report = getattr(state, "multiomics_report", None) or {}

        if not deg_report or not dam_report:
            state.add_trace("quadrant_plot", "completed",
                            summary={"reason": "missing DEG or DAM data"})
            return {"quadrant_plot_report": {"reason": "missing DEG or DAM data"}}

        # Get DEG log2FC data
        deg_fc: dict[str, float] = {}
        for g in deg_report.get("top_genes", []):
            gid = g.get("gene_id", "")
            lfc = g.get("max_log2fc", 0)
            if gid and lfc != 0:
                deg_fc[gid] = float(lfc)

        # Get DAM log2FC data
        dam_fc: dict[str, float] = {}
        for m in dam_report.get("top_metabolites", []):
            mid = m.get("metabolite", "")
            lfc = m.get("max_log2fc", 0)
            if mid and lfc != 0:
                dam_fc[mid] = float(lfc)

        if not deg_fc or not dam_fc:
            state.add_trace("quadrant_plot", "completed",
                            summary={"reason": "no log2FC data available"})
            return {"quadrant_plot_report": {"reason": "no log2FC data"}}

        # If we have multiomics correlation pairs, use those for quadrant classification
        pairs: list[dict] = []
        if multiomics_report and multiomics_report.get("top_pairs"):
            for p in multiomics_report["top_pairs"]:
                # Use the same default evidence gate as the correlation network.
                # A pair without a measured correlation is not treated as an edge.
                corr = float(p.get("correlation", 0) or 0)
                p_value = float(p.get("p_value", 1) or 1)
                if (corr == 0 and p.get("p_value") is not None) or abs(corr) < 0.8 or p_value >= 0.05:
                    continue
                gid = p.get("gene_id", "")
                mid = p.get("metabolite", "")
                g_lfc = deg_fc.get(gid, _infer_log2fc(deg_report, gid))
                m_lfc = dam_fc.get(mid, 0.0)
                if g_lfc == 0 and m_lfc == 0:
                    continue
                quadrant = _classify_quadrant(g_lfc, m_lfc)
                pairs.append({
                    "gene_id": gid,
                    "metabolite": mid,
                    "gene_log2fc": g_lfc,
                    "meta_log2fc": m_lfc,
                    "correlation": p.get("correlation", 0),
                    "quadrant": quadrant,
                })
        else:
            # Use all DEG × DAM combinations (limited)
            for gid, g_lfc in list(deg_fc.items())[:50]:
                for mid, m_lfc in list(dam_fc.items())[:30]:
                    if g_lfc == 0 and m_lfc == 0:
                        continue
                    quadrant = _classify_quadrant(g_lfc, m_lfc)
                    pairs.append({
                        "gene_id": gid,
                        "metabolite": mid,
                        "gene_log2fc": g_lfc,
                        "meta_log2fc": m_lfc,
                        "correlation": 0,
                        "quadrant": quadrant,
                    })

        if not pairs:
            state.add_trace("quadrant_plot", "completed",
                            summary={"reason": "no valid pairs to classify"})
            return {"quadrant_plot_report": {"reason": "no valid pairs"}}

        # Count quadrants
        quadrant_counts: dict[str, int] = {}
        quadrant_genes: dict[str, list[str]] = {}
        for p in pairs:
            q = p["quadrant"]
            quadrant_counts[q] = quadrant_counts.get(q, 0) + 1
            quadrant_genes.setdefault(q, []).append(
                f"{p['gene_id']}-{p['metabolite']}"
            )

        # Generate plot
        target = (state.planner_state.target_metabolite or "metabolite").strip() if state.planner_state else "metabolite"
        plot_b64 = _generate_quadrant_plot(
            pairs,
            threshold=0.0,  # Four-quadrant by default
            gene_label="Gene log2FC",
            meta_label="Metabolite log2FC",
            title=f"DEG × DAM Quadrant Classification ({target})",
        )

        # Build interpretation
        n_total = len(pairs)
        n_syn_up = quadrant_counts.get("Q1", 0)  # Coordinated up
        n_syn_down = quadrant_counts.get("Q3", 0)  # Coordinated down
        n_inv = quadrant_counts.get("Q2", 0) + quadrant_counts.get("Q4", 0)
        syn_ratio = (n_syn_up + n_syn_down) / max(n_total, 1)

        interpretation = ""
        if syn_ratio >= 0.6:
            interpretation = (
                f"Strong coordinated regulation: {syn_ratio:.0%} of gene-metabolite pairs "
                f"show consistent direction (Q1+Q3). Suggests transcriptional control "
                f"of metabolite accumulation."
            )
        elif syn_ratio >= 0.4:
            interpretation = (
                f"Moderate coordination: {syn_ratio:.0%} of pairs are directionally consistent. "
                f"May involve both transcriptional and post-transcriptional regulation."
            )
        else:
            interpretation = (
                f"Weak directional coordination ({syn_ratio:.0%}). "
                f"Metabolite changes may be largely independent of transcript-level regulation."
            )

        quadrant_plot_report = {
            "n_pairs": n_total,
            "mode": "four_quadrant",
            "correlation_threshold": 0.8,
            "correlation_gate": "|r|>=0.8 and p<0.05 when using measured correlations",
            "quadrant_counts": quadrant_counts,
            "Q1_coordinated_up": n_syn_up,
            "Q2_gene_down_meta_up": quadrant_counts.get("Q2", 0),
            "Q3_coordinated_down": n_syn_down,
            "Q4_gene_up_meta_down": quadrant_counts.get("Q4", 0),
            "synchronicity_ratio": round(syn_ratio, 3),
            "interpretation": interpretation,
            "plot_base64": plot_b64,
            "top_coordinated_up": quadrant_genes.get("Q1", [])[:20],
            "top_coordinated_down": quadrant_genes.get("Q3", [])[:20],
        }

        # Save figure to disk for SSE figure events
        if plot_b64:
            file_url = None
            # ── Try R first ──
            try:
                from phyto_reason.visualization.r_bridge import r_available, r_quadrant
                if r_available():
                    pts = {}
                    for p in pairs:
                        key = f"{p['gene_id']}|{p['metabolite']}"
                        pts[key] = {
                            "x": p.get("gene_log2fc", 0),
                            "y": p.get("meta_log2fc", 0),
                            "label": f"{p['gene_id']}:{p['metabolite']}",
                            "quadrant": p.get("quadrant", "Q0"),
                        }
                    r_url = r_quadrant(pts, title=f"DEG × DAM Quadrant Classification ({target})",
                                       x_label="Gene log2FC", y_label="Metabolite log2FC")
                    if r_url:
                        file_url = r_url
            except Exception as e:
                logger.debug("R quadrant plot skipped: %s", e)

            # ── Fallback: matplotlib ──
            if not file_url:
                try:
                    from phyto_reason.visualization.multiomics_figures import _save
                    import matplotlib
                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt
                    import io as _io
                    img_bytes = base64.b64decode(plot_b64.split(",")[1] if "," in plot_b64 else plot_b64)
                    fig, ax = plt.subplots(figsize=(8, 7))
                    from PIL import Image
                    img = Image.open(_io.BytesIO(img_bytes))
                    ax.imshow(img)
                    ax.axis("off")
                    file_url = _save(fig, "quadrant")
                except Exception as e:
                    logger.debug("Quadrant figure save skipped: %s", e)

            if file_url:
                quadrant_plot_report["figure_url"] = file_url
                quadrant_plot_report["figure_markdown"] = f"![Quadrant Plot]({file_url})"

        state.quadrant_plot_report = quadrant_plot_report

        state.add_trace("quadrant_plot", "completed", summary={
            "n_pairs": n_total,
            "syn_ratio": round(syn_ratio, 3),
            "q1": n_syn_up,
            "q3": n_syn_down,
        })

        logger.info(
            "Quadrant plot: %d pairs, sync=%.3f, Q1=%d, Q3=%d",
            n_total, syn_ratio, n_syn_up, n_syn_down,
        )

        return {"quadrant_plot_report": quadrant_plot_report}

    except Exception as e:
        logger.error("quadrant_plot_node failed: %s", e, exc_info=True)
        state.add_trace("quadrant_plot", "failed", error=str(e))
        return {"quadrant_plot_report": {"error": str(e)}}
