"""
analysis.py — 分析执行、图表列表与报告导出端点。
"""

from __future__ import annotations

from fastapi import APIRouter
from pydantic import BaseModel, Field

from phyto_reason.api.schemas.request_models import AnalysisRequest
from phyto_reason.api.schemas.response_models import AnalysisResponse
from phyto_reason.api.services.analysis_service import AnalysisService

router = APIRouter(tags=["analysis"])

_service = AnalysisService()


class ReportExportRequest(BaseModel):
    session_id: str


@router.post("/run_analysis", response_model=AnalysisResponse)
async def run_analysis(request: AnalysisRequest):
    """执行完整分析流程。"""
    return _service.run_analysis(request)


@router.delete("/figures")
async def clear_figures():
    """清空图库：删除用户图表目录下全部已生成图表（PNG/SVG，均可再生）。"""
    import logging

    from phyto_reason.platform_paths import figures_dir

    root = figures_dir()
    deleted = 0
    for pattern in ("*.png", "*.svg"):
        for path in root.glob(pattern):
            try:
                path.unlink()
                deleted += 1
            except OSError:
                continue
    logging.getLogger("api.analysis").info("Cleared %d figures from %s", deleted, root)
    return {"deleted": deleted, "directory": str(root)}


@router.get("/figures")
async def list_figures(limit: int = 50, since: float | None = None):
    """列出用户图表目录：PNG 预览 + 同名 SVG 矢量下载（按修改时间倒序）。

    `since`（Unix 秒）用于会话隔离：前端默认只取本次会话产生的图表；
    传 since 则只返回该时间点之后修改的文件。
    """
    import logging

    from phyto_reason.platform_paths import figures_dir
    from phyto_reason.reports.analysis_report import FIGURE_CAPTIONS, _figure_prefix

    root = figures_dir()
    items = []
    try:
        pngs = sorted(root.glob("*.png"), key=lambda p: p.stat().st_mtime, reverse=True)
    except OSError:
        pngs = []
    for path in pngs:
        try:
            if since is not None and path.stat().st_mtime < since:
                continue
        except OSError:
            continue
        try:
            svg = path.with_suffix(".svg")
            prefix = _figure_prefix(path.stem)
            caption_en, caption_zh, _howto = FIGURE_CAPTIONS.get(prefix, ("", "", ""))
            items.append({
                "name": path.stem,
                "png": f"/static/figures/{path.name}",
                "svg": f"/static/figures/{svg.name}" if svg.exists() else None,
                "size_bytes": path.stat().st_size,
                "caption_en": caption_en,
                "caption_zh": caption_zh,
            })
        except OSError:
            continue
        if len(items) >= max(1, min(limit, 200)):
            break
    logging.getLogger("api.analysis").info("Listed %d figures from %s (since=%s)", len(items), root, since)
    return {"figures": items, "directory": str(root)}


class ReportExportRequest(BaseModel):
    session_id: str


@router.post("/report/export")
async def export_report(request: ReportExportRequest):
    """导出会话最近一次管线运行的 SCI 式分析报告（Markdown）。

    分节只覆盖实际完成的节点；图表按会话注册表编号并引用。
    """
    from fastapi import HTTPException

    from phyto_reason.agent.conversation_store import store
    from phyto_reason.reports.analysis_report import render_markdown

    session = store.get(request.session_id)
    if session is None or not getattr(session, "analysis_digest", None):
        raise HTTPException(
            status_code=409,
            detail="该会话尚无分析结果：请先在分析页运行一次分析，再导出报告。",
        )
    digest = session.analysis_digest
    markdown = render_markdown(digest)
    species = (digest.get("species") or "report").replace(" ", "_")[:40]
    target = (digest.get("target_metabolite") or "").replace(" ", "_")[:30]
    filename = f"analysis_report_{species}" + (f"_{target}" if target else "") + ".md"
    return {"markdown": markdown, "filename": filename,
            "n_sections": len(digest.get("sections", [])),
            "figures": digest.get("figure_urls", [])}
