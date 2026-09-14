"""
sample_aligner.py — 样本名自动对齐器。

处理:
  - 命名不一致: "Root_1A" vs "root-1a"
  - 大小写: "Root" vs "root"
  - 分隔符: "_" vs "-" vs "."
  - 多余前缀/后缀
"""

from __future__ import annotations

import re
from typing import Any


class SampleAligner:
    """样本名对齐器。"""

    @staticmethod
    def normalize(name: str) -> str:
        """归一化样本名: lowercase + separator normalization。"""
        n = name.lower().strip()
        n = re.sub(r"[\s_\-\.]+", "_", n)
        n = n.rstrip("_")
        return n

    @staticmethod
    def build_alignment_map(
        expr_samples: list[str],
        meta_samples: list[str],
    ) -> dict[str, str]:
        """构建表达 → 代谢物 样本名映射。

        Returns:
            {expression_sample_name: metabolite_sample_name}
        """
        normalized_expr = {SampleAligner.normalize(s): s for s in expr_samples}
        normalized_meta = {SampleAligner.normalize(s): s for s in meta_samples}

        mapping: dict[str, str] = {}
        for norm_name, expr_name in normalized_expr.items():
            if norm_name in normalized_meta:
                mapping[expr_name] = normalized_meta[norm_name]

        return mapping

    @staticmethod
    def align_matrices(
        expr_matrix: dict[str, dict[str, float]],
        meta_matrix: dict[str, dict[str, float]],
    ) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], list[str]]:
        """对齐两个矩阵的样本名。

        Returns:
            (aligned_expr, aligned_meta, warnings)
        """
        expr_samples = set()
        for vals in expr_matrix.values():
            expr_samples.update(vals.keys())

        meta_samples = set()
        for vals in meta_matrix.values():
            meta_samples.update(vals.keys())

        warnings = []
        mapping = SampleAligner.build_alignment_map(
            list(expr_samples), list(meta_samples)
        )

        if not mapping:
            warnings.append("No common samples found between matrices")
            return expr_matrix, meta_matrix, warnings

        common_expr = list(mapping.keys())
        common_meta = list(mapping.values())

        aligned_expr: dict[str, dict[str, float]] = {}
        for gid, vals in expr_matrix.items():
            aligned_expr[gid] = {
                s: vals[s] for s in common_expr if s in vals
            }

        aligned_meta: dict[str, dict[str, float]] = {}
        for mid, vals in meta_matrix.items():
            aligned_meta[mid] = {
                s: vals[s] for s in common_meta if s in vals
            }

        if len(common_expr) < len(expr_samples):
            warnings.append(f"Aligned {len(common_expr)}/{len(expr_samples)} expression samples")
        if len(common_meta) < len(meta_samples):
            warnings.append(f"Aligned {len(common_meta)}/{len(meta_samples)} metabolite samples")

        return aligned_expr, aligned_meta, warnings
