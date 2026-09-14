"""
app.py — FastAPI 应用入口。
"""

from __future__ import annotations
import logging
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from phyto_reason.api.routes.health import router as health_router
from phyto_reason.api.routes.upload import router as upload_router
from phyto_reason.api.routes.analysis import router as analysis_router
from phyto_reason.api.routes.chat import router as chat_router
from phyto_reason.api.routes.settings import router as settings_router
from phyto_reason.platform_paths import figures_dir

# ── Logging ──────────────────────────────────────────────────
from phyto_reason.utils.logger import setup_logging as _setup_logging
_log = _setup_logging("phyto_reason", level=logging.INFO)

# Also configure uvicorn loggers to output through our handlers
for _uv_name in ("uvicorn", "uvicorn.access", "uvicorn.error"):
    _uv_log = logging.getLogger(_uv_name)
    _uv_log.handlers.clear()
    _uv_log.propagate = False  # don't bubble up to root
    _uv_log.setLevel(logging.INFO)
    # Use the same handlers as our root plantomics logger
    for _h in _log.handlers:
        _uv_log.addHandler(_h)

_log.info("PhytoReason-Agent v5.0 starting — logging configured")

app = FastAPI(
    title="PhytoReason-Agent API",
    description="多物种药用植物次生代谢调控科研推理 Agent",
    version="5.0.0",
)

# ── CORS ────────────────────────────────────────────────────
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── API routes ──────────────────────────────────────────────
app.include_router(health_router)
app.include_router(upload_router)
app.include_router(analysis_router)
app.include_router(chat_router)
app.include_router(settings_router)

# ── Frontend static ─────────────────────────────────────────
_STATIC_DIR = Path(__file__).resolve().parent / "static"
_STATIC_DIR.mkdir(parents=True, exist_ok=True)


@app.get("/static/figures/{figure_path:path}", include_in_schema=False)
async def serve_user_figure(figure_path: str):
    """按请求动态解析图表目录——目录覆盖保存后立即生效，无需重启。"""
    from fastapi import HTTPException

    root = figures_dir().resolve()
    target = (root / figure_path).resolve()
    if not str(target).startswith(str(root)) or not target.is_file():
        raise HTTPException(status_code=404, detail="figure not found")
    return FileResponse(target)


app.mount("/static", StaticFiles(directory=str(_STATIC_DIR)), name="static")


@app.get("/")
async def serve_frontend():
    """Serve the main frontend page."""
    index_path = _STATIC_DIR / "index.html"
    if index_path.exists():
        return FileResponse(index_path)
    return {"message": "PhytoReason-Agent API is running. Frontend not built."}


def run():
    """使用 uvicorn 启动服务。"""
    import uvicorn
    uvicorn.run(
        "phyto_reason.api.app:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        log_level="info",
        log_config=None,  # disable uvicorn log config, use our own
    )
    # Set uvicorn loggers to propagate to root
    logging.getLogger("uvicorn").handlers.clear()
    logging.getLogger("uvicorn.access").handlers.clear()
    logging.getLogger("uvicorn").addHandler(logging.StreamHandler())
    logging.getLogger("uvicorn.access").addHandler(logging.StreamHandler())


if __name__ == "__main__":
    run()
