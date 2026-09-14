"""
confounder_engine.py — 混杂因子校正引擎。

当存在 tissue / sample_group / batch metadata 时，
自动构建 covariate matrix 并对 TF-enzyme correlation
做 partial correlation 校正。

减少组织特异表达导致的伪相关。
"""

from __future__ import annotations

from typing import Any

import numpy as np


class ConfounderEngine:
    """混杂因子校正引擎。

    用法:
        engine = ConfounderEngine()
        r, p = engine.partial_corr(x, y, covariates)
    """

    def __init__(self, min_samples: int = 5) -> None:
        self.min_samples = min_samples

    def partial_corr(
        self,
        x: np.ndarray,
        y: np.ndarray,
        covariates: np.ndarray | None = None,
    ) -> tuple[float, float]:
        """计算偏相关系数（控制协变量后的 conditional correlation）。

        Args:
            x: 基因表达向量
            y: 目标向量（另一基因或代谢物）
            covariates: (n_samples, n_covariates) 协变量矩阵

        Returns:
            (partial_r, p_value)

        References:
            https://en.wikipedia.org/wiki/Partial_correlation
        """
        n = len(x)
        if n < self.min_samples:
            return 0.0, 1.0

        if covariates is None or covariates.shape[1] == 0:
            from phyto_reason.utils.stats_utils import compute_correlation
            return compute_correlation(x, y, method="bicor")

        x_resid = self._regress_out(x, covariates)
        y_resid = self._regress_out(y, covariates)

        from phyto_reason.utils.stats_utils import compute_fdr, pearson_pvalue
        from scipy.stats import pearsonr

        r, _ = pearsonr(x_resid, y_resid)
        dof = n - covariates.shape[1] - 2
        if dof < 1:
            return 0.0, 1.0

        t_stat = r * np.sqrt(dof / max(1 - r * r, 1e-15))
        from scipy.stats import t as t_dist
        p = 2.0 * (1.0 - t_dist.cdf(abs(t_stat), df=dof))
        return float(r), float(p)

    def build_covariate_matrix(
        self,
        metadata: dict[str, list[str]],
        samples: list[str],
    ) -> np.ndarray:
        """从 metadata 构建协变量矩阵。

        Args:
            metadata: {tissue: [sample1, sample2, ...], ...}
                      或 {batch: [sample1, ...], ...}
            samples: 所有样本名的有序列表

        Returns:
            (n_samples, n_covariates) 的 one-hot 矩阵
        """
        if not metadata:
            return np.empty((len(samples), 0))

        covariates = []
        for category, members in metadata.items():
            vec = np.array([1.0 if s in members else 0.0 for s in samples])
            if np.std(vec) > 0:
                covariates.append(vec)

        if not covariates:
            return np.empty((len(samples), 0))

        return np.column_stack(covariates)

    @staticmethod
    def _regress_out(y: np.ndarray, covariates: np.ndarray) -> np.ndarray:
        """从 y 中回归掉协变量的影响，返回残差。"""
        X = np.column_stack([np.ones(len(y)), covariates])
        try:
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            return y - X @ beta
        except np.linalg.LinAlgError:
            return y - np.mean(y)

    @staticmethod
    def detect_tissue_from_samples(
        sample_names: list[str],
        tissue_keywords: dict[str, list[str]] | None = None,
    ) -> dict[str, list[str]]:
        """从样本名中自动检测组织信息。

        通过样本名关键词匹配推断组织归属。

        Args:
            sample_names: 样本名列表
            tissue_keywords: {tissue: [keyword1, keyword2, ...]} 自定义映射
                             默认内置常见植物组织关键词

        Returns:
            {tissue: [sample1, sample2, ...], ...}
        """
        default_keywords = {
            "root": ["root", "radicle", "rhiz"],
            "leaf": ["leaf", "leave", "foliar"],
            "stem": ["stem", "culm", "stalk"],
            "flower": ["flower", "blossom", "floral"],
            "fruit": ["fruit", "berry", "pod", "seed"],
            "bark": ["bark", "cortex"],
            "callus": ["callus", "calli"],
        }
        keywords = tissue_keywords or default_keywords

        tissue_map: dict[str, list[str]] = {}
        for sample in sample_names:
            s_lower = sample.lower()
            for tissue, kws in keywords.items():
                if any(kw in s_lower for kw in kws):
                    tissue_map.setdefault(tissue, []).append(sample)
                    break

        return tissue_map
