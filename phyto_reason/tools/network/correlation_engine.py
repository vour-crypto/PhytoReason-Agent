"""
correlation_engine.py — 统一相关计算引擎。

支持：
  - pearson
  - spearman
  - bicor（默认，优先）
  - partial correlation（调用 confounder_engine）
"""

from __future__ import annotations

from typing import Any

import numpy as np

from phyto_reason.utils.stats_utils import compute_correlation, compute_fdr, pearson_pvalue


class CorrelationEngine:
    """统一相关计算引擎。

    用法:
        engine = CorrelationEngine(method="bicor")
        r, p = engine.calculate(vec_a, vec_b)
        results = engine.batch_calculate(matrix_a, matrix_b)
    """

    def __init__(
        self,
        method: str = "bicor",
        corr_threshold: float = 0.8,
        fdr_threshold: float = 0.05,
        min_samples: int = 3,
    ) -> None:
        self.method = method
        self.corr_threshold = corr_threshold
        self.fdr_threshold = fdr_threshold
        self.min_samples = min_samples

    def calculate(
        self,
        x: np.ndarray,
        y: np.ndarray,
    ) -> tuple[float, float]:
        """计算两个向量间的相关系数和 p-value。"""
        return compute_correlation(x, y, method=self.method)

    def filter_pairs(
        self,
        pairs: list[tuple[str, str, float, float]],
    ) -> list[tuple[str, str, float, float, float]]:
        """对配对列表做 FDR 校正 + 双阈值筛选。

        Args:
            pairs: [(id_a, id_b, r, p), ...]

        Returns:
            [(id_a, id_b, r, p, q), ...]
            仅保留满足 r ≥ corr_threshold 且 q < fdr_threshold 的项
        """
        if not pairs:
            return []

        pvalues = [p for _, _, _, p in pairs]
        qvalues = compute_fdr(pvalues)

        filtered = []
        for (id_a, id_b, r, p), q in zip(pairs, qvalues):
            if abs(r) >= self.corr_threshold and q < self.fdr_threshold:
                filtered.append((id_a, id_b, r, p, q))
        return filtered

    def batch_calculate(
        self,
        vec_dict_a: dict[str, np.ndarray],
        vec_dict_b: dict[str, np.ndarray],
    ) -> list[tuple[str, str, float, float, float]]:
        """批量计算两组向量间的 pairwise correlation。

        Args:
            vec_dict_a: {id: np.array}
            vec_dict_b: {id: np.array}

        Returns:
            FDR 校正后的 [(id_a, id_b, r, p, q), ...]
        """
        pairs: list[tuple[str, str, float, float]] = []
        for id_a, vec_a in vec_dict_a.items():
            if len(vec_a) < self.min_samples:
                continue
            for id_b, vec_b in vec_dict_b.items():
                if len(vec_b) < self.min_samples:
                    continue
                if len(vec_a) != len(vec_b):
                    continue
                r, p = self.calculate(vec_a, vec_b)
                if abs(r) >= self.corr_threshold:
                    pairs.append((id_a, id_b, r, p))

        return self.filter_pairs(pairs)

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "corr_threshold": self.corr_threshold,
            "fdr_threshold": self.fdr_threshold,
            "min_samples": self.min_samples,
        }
