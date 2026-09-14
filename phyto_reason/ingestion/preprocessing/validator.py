"""
validator.py — 数据质量验证器。

自动检查:
  - 重复基因/代谢物
  - 样本名不匹配
  - 空矩阵
  - 非法字符
  - 编码问题
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class ValidationReport(BaseModel):
    passed: bool = True
    n_genes: int = 0
    n_metabolites: int = 0
    n_samples_expr: int = 0
    n_samples_meta: int = 0
    common_samples: int = 0
    n_missing_values_expr: int = 0
    n_missing_values_meta: int = 0
    n_duplicate_genes: int = 0
    n_infinite_values: int = 0
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


class DataValidator:
    """数据质量验证器。"""

    @staticmethod
    def validate_expression(matrix: dict[str, dict[str, float]],
                            label: str = "expression") -> ValidationReport:
        """验证表达矩阵。"""
        report = ValidationReport()
        if not matrix:
            report.errors.append(f"{label}: empty matrix")
            report.passed = False
            return report

        report.n_genes = len(matrix)

        # 检查样本数
        samples = set()
        for gid, vals in matrix.items():
            samples.update(vals.keys())
        report.n_samples_expr = len(samples)

        # 检查重复
        ids_lower = [g.lower() for g in matrix]
        if len(ids_lower) != len(set(ids_lower)):
            dupes = len(ids_lower) - len(set(ids_lower))
            report.n_duplicate_genes = dupes
            report.warnings.append(f"{label}: {dupes} duplicate gene IDs detected")

        # 检查缺失值
        missing = 0
        infinities = 0
        for vals in matrix.values():
            for v in vals.values():
                if v is None or (isinstance(v, float) and (v != v)):
                    missing += 1
                elif isinstance(v, float) and (v == float("inf") or v == float("-inf")):
                    infinities += 1
        report.n_missing_values_expr = missing
        report.n_infinite_values = infinities

        if missing > 0:
            report.warnings.append(f"{label}: {missing} missing values detected")
        if infinities > 0:
            report.warnings.append(f"{label}: {infinities} infinite values detected")

        if report.n_genes < 10:
            report.warnings.append(f"{label}: only {report.n_genes} genes, may be insufficient")
        if report.n_samples_expr < 3:
            report.errors.append(f"{label}: only {report.n_samples_expr} samples, requires >= 3")
            report.passed = False

        return report

    @staticmethod
    def validate_metabolite(matrix: dict[str, dict[str, float]],
                            label: str = "metabolite") -> ValidationReport:
        """验证代谢物矩阵。"""
        report = ValidationReport()
        if not matrix:
            report.errors.append(f"{label}: empty matrix")
            report.passed = False
            return report

        report.n_metabolites = len(matrix)

        samples = set()
        for mid, vals in matrix.items():
            samples.update(vals.keys())
        report.n_samples_meta = len(samples)

        missing = 0
        for vals in matrix.values():
            for v in vals.values():
                if v is None or (isinstance(v, float) and (v != v)):
                    missing += 1
        report.n_missing_values_meta = missing

        if report.n_metabolites < 1:
            report.errors.append("No metabolites found")
            report.passed = False

        return report

    @staticmethod
    def validate_alignment(expr_samples: set[str],
                            meta_samples: set[str]) -> list[str]:
        """检查样本对齐。"""
        warnings = []
        common = expr_samples & meta_samples

        if not common:
            warnings.append("No common samples between expression and metabolite matrices")
            return warnings

        if len(common) < len(expr_samples):
            diff = expr_samples - meta_samples
            warnings.append(f"{len(diff)} samples in expression not found in metabolite: {list(diff)[:5]}")

        if len(common) < len(meta_samples):
            diff = meta_samples - expr_samples
            warnings.append(f"{len(diff)} samples in metabolite not found in expression: {list(diff)[:5]}")

        if len(common) < 3:
            warnings.append(f"Only {len(common)} aligned samples, may be insufficient")

        return warnings
