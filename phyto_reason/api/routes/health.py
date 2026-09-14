"""
health.py — GET /health 健康检查端点。
"""

from __future__ import annotations

from fastapi import APIRouter

from phyto_reason.api.schemas.request_models import HealthResponse

router = APIRouter(tags=["health"])


@router.get("/health", response_model=HealthResponse)
async def health_check():
    """系统健康检查。"""
    return HealthResponse(status="ok")
