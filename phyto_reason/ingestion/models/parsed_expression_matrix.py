"""
parsed_expression_matrix.py — 表达矩阵的规范化内部表示。

纯数据结构。不包含 biological interpretation。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class ParsedExpressionMatrix(BaseModel):
    """规范化的表达矩阵。"""
    matrix: dict[str, dict[str, float]] = Field(default_factory=dict)

    feature_ids: list[str] = Field(default_factory=list)
    sample_ids: list[str] = Field(default_factory=list)

    n_features: int = 0
    n_samples: int = 0

    missing_rate: float = 0.0
    missing_features: list[str] = Field(default_factory=list)
    missing_samples: list[str] = Field(default_factory=list)

    normalization_status: str = "unknown"
    data_type: str = "unknown"

    species: str = ""
    source_file: str = ""
    warnings: list[str] = Field(default_factory=list)
    provenance: dict = Field(default_factory=dict)

    def to_workflow_dict(self) -> dict[str, dict[str, float]]:
        """转换为 WorkflowRunner.run() 期望的 dict 格式。"""
        return self.matrix
