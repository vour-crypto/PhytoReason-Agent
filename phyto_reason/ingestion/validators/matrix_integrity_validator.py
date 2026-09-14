"""
matrix_integrity_validator.py — 矩阵完整性验证。
"""

from __future__ import annotations

import pandas as pd

from phyto_reason.ingestion.models.omics_dataset import QCReport


def validate_matrix_integrity(
    matrix: dict[str, dict[str, float]] | None = None,
    name: str = "matrix",
) -> QCReport:
    """检查矩阵的完整性。"""
    report = QCReport()

    if matrix is None or len(matrix) == 0:
        report.sample_count_warning = True
        report.malformed_warning = f"{name}: empty or None"
        return report

    df = pd.DataFrame.from_dict(matrix, orient="index")
    n_features, n_samples = df.shape

    if n_samples < 3:
        report.sample_count_warning = True
        report.malformed_warning = f"{name}: only {n_samples} samples (< 3)"

    zero_var = df.columns[df.std() == 0].tolist() if df.shape[1] > 0 else []
    report.zero_variance_features = [str(c) for c in zero_var]

    high_miss = df.columns[df.isna().mean() > 0.5].tolist() if df.shape[1] > 0 else []
    report.high_missing_features = [str(c) for c in high_miss]

    dup_features = df.index[df.index.duplicated()].tolist()
    report.duplicated_features = [str(g) for g in dup_features]

    report.all_passed = not any([
        report.sample_count_warning,
        len(report.zero_variance_features) > 0,
        len(report.high_missing_features) > 5,
        len(report.malformed_warning) > 0,
    ])

    return report


def detect_batch_effect(matrix=None, metadata=None, name="matrix"):
    import numpy as np
    from sklearn.decomposition import PCA
    result = {}
    result["batch_score"] = 0.0
    result["confounded"] = False
    result["warnings"] = []
    if matrix is None or len(matrix) < 3:
        return result
    if metadata is None:
        return result
    has_batch = any("batch" in str(c).lower() for c in metadata.columns_present)
    if not has_batch:
        return result
    try:
        df = pd.DataFrame.from_dict(matrix, orient="index").T.fillna(0)
        pca = PCA(n_components=min(5, df.shape[1], df.shape[0]))
        pca.fit(df.values)
        var = float(pca.explained_variance_ratio_[0])
        result["batch_score"] = round(var, 3)
    except Exception as e:
        result["warnings"].append(str(e))
    return result
