"""
request_models.py — API 请求体 Pydantic 模型。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class UploadResponse(BaseModel):
    session_id: str = ""
    file_id: str = ""
    filename: str = ""
    size_bytes: int = 0
    content_type: str = ""


class AnalysisRequest(BaseModel):
    species: str = ""
    target_metabolite: str = ""
    target_pathway: str = ""
    research_question: str = ""
    sample_count: int = Field(default=0, ge=0)
    has_expression: bool = False
    has_metabolite: bool = False
    has_annotation: bool = False
    has_promoter: bool = False
    has_itak: bool = False
    upload_session_id: str = ""
    data_path: str = ""


class HealthResponse(BaseModel):
    status: str = "ok"
    version: str = "4.0.0"
    service: str = "plantomics-agent"
