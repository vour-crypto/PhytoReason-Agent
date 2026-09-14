"""
preprocessing_pipeline.py — 完整预处理流水线编排器。

整合:
  - 文件加载
  - 格式检测
  - 样本对齐
  - 缺失值处理
  - 标准化
  - 验证
  - 日志
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from phyto_reason.ingestion.preprocessing.matrix_loader import MatrixLoader
from phyto_reason.ingestion.preprocessing.sample_aligner import SampleAligner
from phyto_reason.ingestion.preprocessing.missing_value_handler import MissingValueHandler
from phyto_reason.ingestion.preprocessing.normalization import Normalizer
from phyto_reason.ingestion.preprocessing.validator import DataValidator, ValidationReport

logger = logging.getLogger("preprocessing_pipeline")


class PreprocessingResult:
    """预处理结果。"""
    def __init__(self, expr_matrix: dict, meta_matrix: dict,
                 promoter_seqs: dict | None, validation: ValidationReport,
                 warnings: list[str], metadata: dict) -> None:
        self.expr_matrix = expr_matrix
        self.meta_matrix = meta_matrix
        self.promoter_seqs = promoter_seqs
        self.validation = validation
        self.warnings = warnings
        self.metadata = metadata

    def to_summary(self) -> dict:
        return {
            "genes": len(self.expr_matrix),
            "metabolites": len(self.meta_matrix),
            "samples": self.validation.common_samples,
            "has_promoters": self.promoter_seqs is not None,
            "warnings": len(self.warnings),
            "passed": self.validation.passed,
        }


class PreprocessingPipeline:
    """完整预处理流水线。"""

    def __init__(self) -> None:
        self.loader = MatrixLoader()
        self.aligner = SampleAligner()
        self.mvh = MissingValueHandler()
        self.normalizer = Normalizer()
        self.validator = DataValidator()

    def run(self, expr_path: str | Path, meta_path: str | Path,
            promoter_path: str | Path | None = None,
            normalize_method: str = "log2",
            fill_method: str = "mean",
            min_mean_expr: float = 0.5) -> PreprocessingResult:
        """执行完整预处理。"""
        warnings: list[str] = []
        metadata: dict[str, Any] = {}

        logger.info(f"Loading expression: {expr_path}")
        expr = self.loader.load_expression_matrix(expr_path)
        logger.info(f"  -> {len(expr)} genes loaded")

        logger.info(f"Loading metabolite: {meta_path}")
        meta = self.loader.load_metabolite_matrix(meta_path)
        logger.info(f"  -> {len(meta)} metabolites loaded")

        promoters = None
        if promoter_path:
            logger.info(f"Loading promoters: {promoter_path}")
            promoters = self.loader.load_promoter_fasta(promoter_path)
            logger.info(f"  -> {len(promoters)} sequences loaded")

        # Sample alignment
        logger.info("Aligning sample names...")
        expr, meta, align_warnings = self.aligner.align_matrices(expr, meta)
        warnings.extend(align_warnings)

        common_samples = set()
        for vals in expr.values():
            common_samples.update(vals.keys())
        metadata["n_samples"] = len(common_samples)

        # Missing values
        logger.info("Handling missing values...")
        expr_report = self.mvh.report(expr)
        if expr_report.replaced > 0:
            logger.info(f"  {expr_report.replaced} missing values in expression")
            if fill_method == "mean":
                expr = self.mvh.replace_with_mean(expr)
            else:
                expr = self.mvh.replace_with_zero(expr)

        meta_report = self.mvh.report(meta)
        if meta_report.replaced > 0:
            logger.info(f"  {meta_report.replaced} missing values in metabolite")
            meta = self.mvh.replace_with_mean(meta)

        # Low expression filtering
        n_before = len(expr)
        expr = self.mvh.filter_low_expression(expr, min_mean=min_mean_expr)
        n_removed = n_before - len(expr)
        if n_removed > 0:
            logger.info(f"  Filtered {n_removed} low-expression genes")
            warnings.append(f"Removed {n_removed} low-expression genes")

        # Normalization
        logger.info(f"Normalizing ({normalize_method})...")
        if normalize_method == "log2":
            expr = self.normalizer.log2_transform(expr)
        elif normalize_method == "zscore":
            expr = self.normalizer.z_score(expr)
        elif normalize_method == "minmax":
            expr = self.normalizer.min_max(expr)
        metadata["normalization"] = normalize_method

        # Validation
        logger.info("Validating...")
        val = DataValidator.validate_expression(expr)
        val.common_samples = len(common_samples)
        warnings.extend(val.warnings)
        warnings.extend(val.errors)

        metadata["n_genes"] = len(expr)
        metadata["n_metabolites"] = len(meta)

        logger.info(f"Preprocessing complete: {len(expr)} genes, "
                     f"{len(meta)} metabolites, {len(common_samples)} samples")

        return PreprocessingResult(
            expr_matrix=expr, meta_matrix=meta,
            promoter_seqs=promoters,
            validation=val, warnings=warnings, metadata=metadata,
        )
