"""
stats_utils.py — 统计工具函数。

所有 correlation-based 筛选必须经过 FDR correction。
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np


def compute_fdr(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction。

    Args:
        pvalues: 原始 p-value 列表

    Returns:
        校正后的 q-value 列表（长度与输入相同）
    """
    n = len(pvalues)
    if n == 0:
        return []
    if n == 1:
        return [min(pvalues[0], 1.0)]

    ranked = sorted([(p, i) for i, p in enumerate(pvalues)])
    qvals = [0.0] * n
    prev = 1.0
    for rank, (p, idx) in enumerate(reversed(ranked), 1):
        q = p * n / (n - rank + 1)
        q = min(q, prev, 1.0)
        qvals[idx] = q
        prev = q
    return qvals


def pearson_pvalue(r: float, n: int) -> float:
    """Pearson correlation 的 p-value（t-distribution）。"""
    if n < 3 or abs(r) >= 1.0:
        return 0.0 if abs(r) >= 1.0 else 1.0
    t_stat = r * math.sqrt((n - 2) / (1 - r * r + 1e-15))
    from scipy.stats import t as t_dist
    return 2.0 * (1.0 - t_dist.cdf(abs(t_stat), df=n - 2))


def bicor(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Biweight midcorrelation（bicor）。

    参考 WGCNA 的 bicor 实现。
    对 outlier 鲁棒，适合小样本组学数据。

    Returns:
        (bicor_r, p_value)
    """
    n = len(x)
    if n < 3:
        return 0.0, 1.0

    x_med = np.median(x)
    y_med = np.median(y)
    x_mad = np.median(np.abs(x - x_med)) + 1e-10
    y_mad = np.median(np.abs(y - y_med)) + 1e-10

    u = (x - x_med) / (9.0 * x_mad)
    v = (y - y_med) / (9.0 * y_mad)

    w_x = (1 - u ** 2) ** 2
    w_y = (1 - v ** 2) ** 2
    w_x = np.where(np.abs(u) < 1, w_x, 0.0)
    w_y = np.where(np.abs(v) < 1, w_y, 0.0)

    x_centered = x - x_med
    y_centered = y - y_med

    num = np.sum(w_x * x_centered * w_y * y_centered)
    den = math.sqrt(np.sum(w_x * x_centered ** 2) * np.sum(w_y * y_centered ** 2)) + 1e-15
    r = num / den
    r = max(-1.0, min(1.0, r))

    p = pearson_pvalue(r, n)
    return r, p


def compute_correlation(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "bicor",
) -> tuple[float, float]:
    """统一相关计算接口。

    Args:
        x, y: 等长向量
        method: pearson | spearman | bicor

    Returns:
        (r, p_value)
    """
    n = len(x)
    if n < 3:
        return 0.0, 1.0

    if method == "pearson":
        from scipy.stats import pearsonr
        r, p = pearsonr(x, y)
        return float(r), float(p)

    elif method == "spearman":
        from scipy.stats import spearmanr
        r, p = spearmanr(x, y)
        return float(r), float(p)

    elif method == "bicor":
        return bicor(x, y)

    raise ValueError(f"Unknown correlation method: {method}")


def compute_size_factors(counts: dict[str, dict[str, float]]) -> dict[str, float]:
    """DESeq2 式 median-of-ratios 文库大小因子。

    counts: gene -> sample -> count（原始计数）。
    返回 sample -> size factor（>0，均值≈1）。仅使用全样本均 >0 的基因；
    若无法估计（如全基因非全样本覆盖），退化为全 1（调用方按未归一处理）。
    """
    samples: list[str] = []
    for sample_vals in counts.values():
        for sample in sample_vals:
            if sample not in samples:
                samples.append(sample)
    if not samples:
        return {}
    ratios: dict[str, list[float]] = {s: [] for s in samples}
    for sample_vals in counts.values():
        vals: list[float] = []
        ok = True
        for s in samples:
            v = float(sample_vals.get(s, float("nan")))
            if math.isnan(v) or v <= 0:
                ok = False
                break
            vals.append(v)
        if not ok:
            continue
        gmean = math.exp(sum(math.log(v) for v in vals) / len(vals))
        if gmean <= 0:
            continue
        for s, v in zip(samples, vals):
            ratios[s].append(v / gmean)

    raw = {}
    for s in samples:
        rs = ratios[s]
        raw[s] = float(np.median(rs)) if rs else 1.0
    positive = [f for f in raw.values() if f > 0]
    if not positive:
        return {s: 1.0 for s in samples}
    mean_f = sum(positive) / len(positive)
    return {s: (raw[s] / mean_f) if raw[s] > 0 else 1.0 for s in samples}


def compute_size_factors_array(mat: "np.ndarray") -> "np.ndarray":
    """矩阵版 median-of-ratios：输入 基因×样本，返回每列 size factor（长度=样本数）。

    仅用全样本 >0 的行估计；无可用行时返回全 1。
    """
    mat = np.asarray(mat, dtype=float)
    if mat.ndim != 2 or mat.shape[1] == 0:
        return np.ones(mat.shape[1] if mat.ndim == 2 else 0)
    positive = mat > 0
    complete_rows = positive.all(axis=1)
    sub = mat[complete_rows]
    if sub.shape[0] == 0:
        return np.ones(mat.shape[1])
    # 每基因跨样本几何均值作参考（DESeq2 median-of-ratios），再取每样本比值中位数
    with np.errstate(divide="ignore"):
        log_geo = np.log(sub).mean(axis=1)
    ref = np.exp(log_geo)
    ref[ref <= 0] = np.nan
    ratio = sub / ref[:, None]
    factors = np.median(ratio, axis=0)
    positive_factors = factors[factors > 0]
    mean_f = float(positive_factors.mean()) if positive_factors.size else 1.0
    factors = np.where(factors > 0, factors / mean_f, 1.0)
    return factors
