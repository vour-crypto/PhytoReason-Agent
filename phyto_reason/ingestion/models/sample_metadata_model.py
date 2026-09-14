"""
sample_metadata_model.py — 样本元数据的规范化表示。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class SampleMetadata(BaseModel):
    """样本级元数据。"""
    sample_id: str
    condition: str = ""
    treatment: str = ""
    tissue: str = ""
    genotype: str = ""
    replicate: str = ""
    batch: str = ""
    extra: dict = Field(default_factory=dict)


class SampleMetadataTable(BaseModel):
    """完整的样本元数据表。"""
    samples: list[SampleMetadata] = Field(default_factory=list)
    sample_ids: list[str] = Field(default_factory=list)

    columns_present: list[str] = Field(default_factory=list)
    columns_missing: list[str] = Field(default_factory=list)

    source_file: str = ""
    warnings: list[str] = Field(default_factory=list)

    def has_column(self, col: str) -> bool:
        return col in self.columns_present

    def get_condition_map(self) -> dict[str, str]:
        return {s.sample_id: s.condition for s in self.samples}

    def get_tissue_map(self) -> dict[str, str]:
        return {s.sample_id: s.tissue for s in self.samples if s.tissue}
