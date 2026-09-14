from __future__ import annotations

import numpy as np

from ._common import PALETTE, new_figure, save_figure


class HeatmapPlotter:
    @staticmethod
    def clustered_heatmap(matrix: dict, feature_ids: list[str] | None = None,
                          title: str = "Expression Heatmap"):
        feature_ids = feature_ids or list(matrix)
        feature_ids = [f for f in feature_ids if f in matrix]
        samples = list(next(iter(matrix.values())).keys()) if feature_ids else []
        values = np.array([[float(matrix[f].get(s, 0)) for s in samples] for f in feature_ids])
        if values.size == 0:
            raise ValueError("No numeric values available for heatmap")
        values = (values - values.mean(axis=1, keepdims=True)) / (values.std(axis=1, keepdims=True) + 1e-9)
        fig, ax = new_figure(title, width=7.5, height=max(4.0, min(9.0, 2.5 + len(feature_ids) * 0.16)))
        image = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_xticks(range(len(samples)), samples, rotation=45, ha="right", fontsize=8)
        if len(feature_ids) <= 50:
            ax.set_yticks(range(len(feature_ids)), feature_ids, fontsize=8)
        else:
            ax.set_yticks([])  # 行数过多时标签不可读，隐藏（图例交给调用方筛选）
        fig.colorbar(image, ax=ax, label="row z-score", shrink=0.8)
        return fig
