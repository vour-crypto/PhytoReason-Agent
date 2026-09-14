from __future__ import annotations

import numpy as np

from ._common import GRADIENT_HIGH, PALETTE, new_figure


class VolcanoPlotter:
    """论文样式火山图：显著上调红 / 显著下调蓝 / 非显著灰，黑色虚线阈值。"""

    @staticmethod
    def deg_volcano(data: dict, title: str = "Differential Expression",
                    padj_threshold: float = 0.05, lfc_threshold: float = 1.0):
        fig, ax = new_figure(title)
        up_x, up_y, down_x, down_y, ns_x, ns_y = [], [], [], [], [], []
        top_labels: list[tuple[float, float, str]] = []
        for name, item in data.items():
            lfc = float(item.get("log2fc", item.get("effect_size", 0)) or 0)
            padj = max(float(item.get("padj", item.get("q_value", 1)) or 1), 1e-300)
            y = -np.log10(padj)
            significant = padj < padj_threshold and abs(lfc) >= lfc_threshold
            if significant and lfc >= 0:
                up_x.append(lfc); up_y.append(y)
            elif significant:
                down_x.append(lfc); down_y.append(y)
            else:
                ns_x.append(lfc); ns_y.append(y)
            top_labels.append((y, lfc, name))

        # 非显著在下层，显著在上层
        ax.scatter(ns_x, ns_y, c=PALETTE["muted"], s=18, alpha=0.6, edgecolors="none")
        if down_x:
            ax.scatter(down_x, down_y, c=PALETTE["blue"], s=26, alpha=0.85,
                       edgecolors="none", label=f"Down ({len(down_x)})")
        if up_x:
            ax.scatter(up_x, up_y, c=PALETTE["red"], s=26, alpha=0.85,
                       edgecolors="none", label=f"Up ({len(up_x)})")

        # 论文样式：黑色虚线阈值线
        ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=0.5)
        ax.axvline(-lfc_threshold, color="black", linestyle="--", linewidth=0.5)
        ax.axvline(lfc_threshold, color="black", linestyle="--", linewidth=0.5)

        # Top5 显著基因标注
        top_labels.sort(reverse=True)
        for y_value, lfc, name in top_labels[:5]:
            if y_value > 1.3:
                ax.annotate(name, (lfc, y_value), fontsize=7,
                            xytext=(3, 3), textcoords="offset points")

        if up_x or down_x:
            ax.legend(loc="upper right", fontsize=10)
        ax.set_xlabel("log2(Fold Change)")
        ax.set_ylabel("-log10(padj)")
        return fig
