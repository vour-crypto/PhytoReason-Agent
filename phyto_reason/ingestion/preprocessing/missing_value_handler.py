"""
missing_value_handler.py — 缺失值处理器。

支持:
  - NA / NaN / None
  - 空字符串
  - 无限值 (inf)
  - 零值过滤
  - 低表达过滤
"""

from __future__ import annotations

import math
from typing import Any

from pydantic import BaseModel, Field


class MissingValueReport(BaseModel):
    replaced: int = 0
    removed_genes: int = 0
    removed_metabolites: int = 0
    low_count_genes: int = 0
    method: str = ""


class MissingValueHandler:
    """缺失值处理器。"""

    @staticmethod
    def is_missing(val: Any) -> bool:
        if val is None:
            return True
        if isinstance(val, float):
            if math.isnan(val):
                return True
            if math.isinf(val):
                return True
        if isinstance(val, str):
            if val.strip().lower() in ("", "na", "nan", "null", "none", "-", "inf"):
                return True
        return False

    @staticmethod
    def report(matrix: dict[str, dict[str, float]]) -> MissingValueReport:
        """返回缺失值统计。"""
        n_missing = 0
        for vals in matrix.values():
            for v in vals.values():
                if MissingValueHandler.is_missing(v):
                    n_missing += 1
        return MissingValueReport(replaced=n_missing)

    @staticmethod
    def replace_with_mean(matrix: dict[str, dict[str, float]],
                           axis: str = "row") -> dict[str, dict[str, float]]:
        """用均值填充缺失值。"""
        result: dict[str, dict[str, float]] = {}
        replaced = 0

        for gid, vals in matrix.items():
            clean_vals = [v for v in vals.values()
                          if not MissingValueHandler.is_missing(v)]
            mean = sum(clean_vals) / len(clean_vals) if clean_vals else 0.0

            result[gid] = {}
            for sample, v in vals.items():
                if MissingValueHandler.is_missing(v):
                    result[gid][sample] = round(mean, 6)
                    replaced += 1
                else:
                    result[gid][sample] = v

        return result

    @staticmethod
    def replace_with_zero(matrix: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
        """用 0 填充缺失值。"""
        result: dict[str, dict[str, float]] = {}
        replaced = 0

        for gid, vals in matrix.items():
            result[gid] = {}
            for sample, v in vals.items():
                if MissingValueHandler.is_missing(v):
                    result[gid][sample] = 0.0
                    replaced += 1
                else:
                    result[gid][sample] = v

        return result

    @staticmethod
    def filter_low_expression(matrix: dict[str, dict[str, float]],
                               min_mean: float = 1.0,
                               min_present: float = 0.5) -> dict[str, dict[str, float]]:
        """过滤低表达基因。"""
        result: dict[str, dict[str, float]] = {}
        removed = 0

        for gid, vals in matrix.items():
            n_samples = len(vals)
            n_present = sum(1 for v in vals.values()
                            if not MissingValueHandler.is_missing(v) and v > 0)
            present_ratio = n_present / max(n_samples, 1)
            mean_val = sum(v for v in vals.values()
                           if not MissingValueHandler.is_missing(v)) / max(n_present, 1)

            if present_ratio >= min_present and mean_val >= min_mean:
                result[gid] = vals
            else:
                removed += 1

        return result

    @staticmethod
    def remove_na_genes(matrix: dict[str, dict[str, float]],
                         max_na_ratio: float = 0.8) -> dict[str, dict[str, float]]:
        """移除缺失比例过高的基因。"""
        result: dict[str, dict[str, float]] = {}
        removed = 0

        for gid, vals in matrix.items():
            n_na = sum(1 for v in vals.values() if MissingValueHandler.is_missing(v))
            if n_na / max(len(vals), 1) <= max_na_ratio:
                result[gid] = vals
            else:
                removed += 1

        return result
