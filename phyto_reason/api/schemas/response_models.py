"""
response_models.py — API 响应体 Pydantic 模型。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class CandidateInfo(BaseModel):
    gene_id: str = ""
    tf_family: str | None = None
    confidence_score: float = 0.0
    confidence_level: str = ""
    biological_plausibility: str = ""
    n_evidences: int = 0
    top_scores: dict[str, float] = Field(default_factory=dict)


class AnalysisResponse(BaseModel):
    success: bool = True
    session_id: str = ""
    workflow_summary: dict = Field(default_factory=dict)
    top_candidates: list[CandidateInfo] = Field(default_factory=list)
    confidence: dict = Field(default_factory=dict)
    hypothesis: str = ""
    suggested_validation: list[str] = Field(default_factory=list)
    execution_trace: list[dict] = Field(default_factory=list)
    error: str | None = None


class ErrorResponse(BaseModel):
    success: bool = False
    error: str = ""
    detail: str = ""
