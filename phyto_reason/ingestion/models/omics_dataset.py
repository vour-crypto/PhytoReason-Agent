"""
omics_dataset.py — 统一的多组学数据集。

这是 ingestion 层的最终输出。
包含 表达矩阵 + 代谢物矩阵 + 元数据 + 对齐报告 + QC 报告。
不含任何 biological reasoning。
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from phyto_reason.ingestion.models.parsed_expression_matrix import ParsedExpressionMatrix
from phyto_reason.ingestion.models.parsed_metabolite_matrix import ParsedMetaboliteMatrix
from phyto_reason.ingestion.models.sample_metadata_model import SampleMetadataTable


class AlignmentReport(BaseModel):
    """样本对齐报告。"""
    common_samples: list[str] = Field(default_factory=list)
    n_common: int = 0
    dropped_from_expression: list[str] = Field(default_factory=list)
    dropped_from_metabolite: list[str] = Field(default_factory=list)
    dropped_from_metadata: list[str] = Field(default_factory=list)


class QCReport(BaseModel):
    """质量控制报告。"""
    sample_count_warning: bool = False
    zero_variance_features: list[str] = Field(default_factory=list)
    high_missing_features: list[str] = Field(default_factory=list)
    duplicated_features: list[str] = Field(default_factory=list)
    malformed_warning: str = ""
    all_passed: bool = False


class OmicsDataset(BaseModel):
    """统一的多组学数据集 — ingestion 层最终输出。"""
    expression: ParsedExpressionMatrix | None = None
    metabolite: ParsedMetaboliteMatrix | None = None
    metadata: SampleMetadataTable | None = None

    species: str = ""
    target_metabolite: str = ""

    alignment: AlignmentReport = Field(default_factory=AlignmentReport)
    qc_report: QCReport = Field(default_factory=QCReport)

    warnings: list[str] = Field(default_factory=list)
    source_files: dict[str, str] = Field(default_factory=dict)

    def to_workflow_inputs(self) -> dict:
        """转换为 WorkflowRunner.run() 可接受的参数。"""
        return {
            "expression_matrix": self.expression.to_workflow_dict() if self.expression else None,
            "metabolite_matrix": self.metabolite.to_workflow_dict() if self.metabolite else None,
            "species": self.species,
            "target_metabolite": self.target_metabolite,
            "has_expression": self.expression is not None,
            "has_metabolite": self.metabolite is not None,
            "sample_count": self.alignment.n_common,
        }

    def is_usable(self) -> bool:
        """数据是否足以启动 workflow。"""
        return self.expression is not None and self.alignment.n_common >= 3
