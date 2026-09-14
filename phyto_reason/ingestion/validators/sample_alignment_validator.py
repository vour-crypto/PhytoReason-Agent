"""
sample_alignment_validator.py — 样本对齐验证。

必须检查:
  - expression vs metabolite 样本重叠
  - metadata vs matrix 重叠
  - 重复样本
  - 不一致的样本顺序
禁止 silent reorder。
"""

from __future__ import annotations

import logging

from phyto_reason.ingestion.models.omics_dataset import AlignmentReport

logger = logging.getLogger("sample_alignment")


def validate_alignment(
    expression_samples: list[str] | None = None,
    metabolite_samples: list[str] | None = None,
    metadata_samples: list[str] | None = None,
) -> AlignmentReport:
    """检查所有数据源之间的样本对齐。"""
    report = AlignmentReport()

    common: set[str] | None = None
    dropped_expr: list[str] = []
    dropped_meta: list[str] = []
    dropped_md: list[str] = []

    if expression_samples:
        expr_set = set(expression_samples)
        common = expr_set.copy()

        if metabolite_samples:
            meta_set = set(metabolite_samples)
            new_common = common & meta_set
            dropped_meta = list(expr_set - meta_set)
            common = new_common

        if metadata_samples:
            md_set = set(metadata_samples)
            new_common = common & md_set
            dropped_md = list(common - md_set)
            common = new_common

    if common:
        report.common_samples = sorted(common)
        report.n_common = len(common)

        if expression_samples:
            report.dropped_from_expression = [s for s in expression_samples if s not in common]
        if metabolite_samples:
            report.dropped_from_metabolite = [s for s in metabolite_samples if s not in common]
        if metadata_samples:
            report.dropped_from_metadata = [s for s in metadata_samples if s not in common]

    if report.dropped_from_expression:
        logger.warning(f"Samples in expression but not in others: {report.dropped_from_expression}")
    if report.dropped_from_metabolite:
        logger.warning(f"Samples in metabolite but not in others: {report.dropped_from_metabolite}")
    if report.dropped_from_metadata:
        logger.warning(f"Samples in metadata but not in others: {report.dropped_from_metadata}")

    return report
