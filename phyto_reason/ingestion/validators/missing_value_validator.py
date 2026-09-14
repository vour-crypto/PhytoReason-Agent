"""
missing_value_validator.py — 缺失值验证（纯函数）。

仅做 warning/masking/reporting。
禁止自动 imputation。
"""

from __future__ import annotations

import math


class MissingValueReport:
    def __init__(self):
        self.feature_missingness: dict[str, float] = {}
        self.sample_missingness: dict[str, float] = {}
        self.overall_rate: float = 0.0
        self.high_missing_features: list[str] = []
        self.high_missing_samples: list[str] = []
        self.warnings: list[str] = []
        self.passed: bool = True

    def to_dict(self) -> dict:
        return {
            "overall_missing_rate": round(self.overall_rate, 4),
            "n_high_missing_features": len(self.high_missing_features),
            "n_high_missing_samples": len(self.high_missing_samples),
            "n_warnings": len(self.warnings),
            "passed": self.passed,
        }


def validate_missing_values(
    matrix: dict[str, dict[str, float]] | None = None,
    feature_axis: bool = True,
    max_feature_missing_rate: float = 0.5,
    max_sample_missing_rate: float = 0.3,
) -> MissingValueReport:
    """检查矩阵中的缺失值。"""
    report = MissingValueReport()

    if matrix is None or len(matrix) == 0:
        report.warnings.append("Empty matrix — no missing value analysis possible")
        report.passed = False
        return report

    all_vals = []
    total_cells = 0
    missing_cells = 0

    for feature_id, sample_vals in matrix.items():
        n_missing = sum(1 for v in sample_vals.values() if isinstance(v, float) and math.isnan(v))
        n_total = len(sample_vals)
        rate = n_missing / max(n_total, 1)
        report.feature_missingness[feature_id] = round(rate, 4)
        all_vals.extend(sample_vals.values())
        total_cells += n_total
        missing_cells += n_missing

        if rate > max_feature_missing_rate:
            report.high_missing_features.append(feature_id)

    report.overall_rate = missing_cells / max(total_cells, 1)

    sample_ids = set()
    for fv in matrix.values():
        sample_ids.update(fv.keys())
    for sid in sorted(sample_ids):
        n_miss = sum(1 for fv in matrix.values()
                     if isinstance(fv.get(sid), float) and math.isnan(fv.get(sid, float("nan"))))
        n_total = len(matrix)
        rate = n_miss / max(n_total, 1)
        report.sample_missingness[sid] = round(rate, 4)
        if rate > max_sample_missing_rate:
            report.high_missing_samples.append(sid)

    if report.high_missing_features:
        report.warnings.append(
            f"Features with >{max_feature_missing_rate:.0%} missing: {len(report.high_missing_features)}"
        )
    if report.high_missing_samples:
        report.warnings.append(
            f"Samples with >{max_sample_missing_rate:.0%} missing: {len(report.high_missing_samples)}"
        )

    return report
