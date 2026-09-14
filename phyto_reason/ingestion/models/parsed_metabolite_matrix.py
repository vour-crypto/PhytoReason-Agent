"""
parsed_metabolite_matrix.py — 代谢物矩阵的规范化内部表示。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class ParsedMetaboliteMatrix(BaseModel):
    """规范化的代谢物矩阵。"""
    matrix: dict[str, dict[str, float]] = Field(default_factory=dict)

    feature_ids: list[str] = Field(default_factory=list)
    sample_ids: list[str] = Field(default_factory=list)
    original_names: list[str] = Field(default_factory=list)

    n_features: int = 0
    n_samples: int = 0

    missing_rate: float = 0.0
    missing_features: list[str] = Field(default_factory=list)
    missing_samples: list[str] = Field(default_factory=list)

    normalized: bool = False
    source_file: str = ""
    warnings: list[str] = Field(default_factory=list)
    provenance: dict = Field(default_factory=dict)

    # Column classification details for user feedback
    classification: dict = Field(default_factory=dict)

    def to_workflow_dict(self) -> dict[str, dict[str, float]]:
        return self.matrix
