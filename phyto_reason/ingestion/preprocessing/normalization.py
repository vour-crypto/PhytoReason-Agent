"""
normalization.py — 数据标准化处理器。

支持:
  - TPM
  - log2(TPM + 1)
  - z-score
  - quantile normalization
  - min-max scaling
"""

from __future__ import annotations

import math
from typing import Any


class Normalizer:
    """数据标准化器。"""

    @staticmethod
    def log2_transform(matrix: dict[str, dict[str, float]],
                        pseudo_count: float = 1.0) -> dict[str, dict[str, float]]:
        """log2(x + pseudo_count) 变换。"""
        result: dict[str, dict[str, float]] = {}
        for gid, vals in matrix.items():
            result[gid] = {
                s: round(math.log2(max(v, 0) + pseudo_count), 6)
                for s, v in vals.items()
            }
        return result

    @staticmethod
    def z_score(matrix: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
        """z-score 标准化 (每行)。"""
        result: dict[str, dict[str, float]] = {}
        for gid, vals in matrix.items():
            values = list(vals.values())
            mean = sum(values) / max(len(values), 1)
            var = sum((v - mean) ** 2 for v in values) / max(len(values), 1)
            std = math.sqrt(max(var, 1e-10))

            result[gid] = {
                s: round((v - mean) / std, 6)
                for s, v in vals.items()
            }
        return result

    @staticmethod
    def min_max(matrix: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
        """min-max 归一化到 [0, 1]。"""
        result: dict[str, dict[str, float]] = {}
        for gid, vals in matrix.items():
            values = list(vals.values())
            vmin = min(values)
            vmax = max(values)
            span = max(vmax - vmin, 1e-10)

            result[gid] = {
                s: round((v - vmin) / span, 6)
                for s, v in vals.items()
            }
        return result

    @staticmethod
    def quantile(matrix: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
        """分位数标准化 (所有样本分布相同)。"""
        samples = list({s for vals in matrix.values() for s in vals})
        if not samples:
            return matrix

        # 提取样本向量
        sample_vectors: dict[str, list[float]] = {s: [] for s in samples}
        for vals in matrix.values():
            for s in samples:
                if s in vals:
                    sample_vectors[s].append(vals[s])

        if not sample_vectors or not all(sample_vectors.values()):
            return matrix

        # 排序并计算参考分位
        sorted_vals = {s: sorted(v) for s, v in sample_vectors.items()}
        n = min(len(v) for v in sorted_vals.values()) if sorted_vals else 0
        if n == 0:
            return matrix

        reference = []
        for i in range(n):
            ref_i = sum(sorted_vals[s][i] for s in samples) / len(samples)
            reference.append(ref_i)

        result: dict[str, dict[str, float]] = {}
        for gid, vals in matrix.items():
            sorted_gene = sorted(vals.items(), key=lambda x: x[1])
            result[gid] = {}
            for (sample, _), ref_val in zip(sorted_gene, reference):
                result[gid][sample] = round(ref_val, 6)

            # Ensure all samples are in result
            for s in samples:
                if s not in result[gid]:
                    result[gid][s] = round(sum(reference) / max(n, 1), 6)

        return result
