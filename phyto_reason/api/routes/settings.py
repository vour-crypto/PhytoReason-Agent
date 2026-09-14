"""LLM provider settings and data-directory/backup endpoints.

Phase 6.4 Step 1：
- GET  /settings/data        数据目录与各子路径盘点
- POST /settings/data/open   在系统文件管理器中打开目录
- POST /settings/backup/export  导出用户数据备份（zip）
- POST /settings/backup/import  导入备份（zip-slip 校验 + 1GB 解压上限 + .bak 覆盖保护）
"""

from __future__ import annotations

import io
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path

from fastapi import APIRouter, File, UploadFile
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from phyto_reason.config.llm_config import probe_config, public_config, save_config
from phyto_reason.platform_paths import (
    cache_dir,
    config_path,
    default_backup_dir,
    dir_overrides,
    figures_dir,
    log_dir,
    memory_db_path,
    output_dir,
    report_dir,
    sessions_db_path,
    set_dir_overrides,
    user_data_dir,
    workbench_db_path,
)

router = APIRouter(prefix="/settings", tags=["settings"])

# 任务书 Step 1：导入解压总量上限（防 zip bomb）
MAX_IMPORT_BYTES = 1024**3
# 备份排除项：缓存可再生、backups 防嵌套、临时探测文件
_BACKUP_EXCLUDE_DIRS = {"cache", "backups", "__pycache__"}
_BACKUP_EXCLUDE_NAMES = {".write_test", ".migration_v1.complete"}


class LLMSettingsRequest(BaseModel):
    provider: str = Field(default="deepseek")
    base_url: str = Field(default="https://api.deepseek.com/v1")
    model: str = Field(default="deepseek-chat")
    api_key: str = Field(default="")


class OpenDirRequest(BaseModel):
    target: str = Field(default="data")  # data | figures | outputs | logs | cache


class ExportBackupRequest(BaseModel):
    target_dir: str | None = Field(default=None, description="备份落盘目录；缺省为用户数据目录下 backups/")


class DirectoriesRequest(BaseModel):
    figures_dir: str = Field(default="", description="图表保存目录覆盖；空=默认")
    reports_dir: str = Field(default="", description="报告保存目录覆盖；空=默认")


@router.get("/llm")
async def get_llm_settings():
    return public_config()


@router.post("/llm")
async def update_llm_settings(request: LLMSettingsRequest):
    try:
        save_config(request.model_dump())
        return public_config()
    except ValueError as exc:
        return JSONResponse(status_code=400, content={"error": str(exc)})


@router.post("/llm/probe")
async def probe_llm_settings():
    return probe_config()


# ── 数据目录 ────────────────────────────────────────────────

def _dir_size(path: Path, exclude: set[str] | None = None) -> int:
    exclude = exclude or set()
    total = 0
    for dirpath, dirnames, filenames in os.walk(path):
        dirnames[:] = [d for d in dirnames if d not in exclude]
        for name in filenames:
            try:
                total += (Path(dirpath) / name).stat().st_size
            except OSError:
                continue
    return total


@router.get("/data")
async def get_data_info():
    root = user_data_dir()

    def entry(label: str, path: Path) -> dict:
        return {"label": label, "path": str(path), "exists": path.exists()}

    return {
        "data_dir": str(root),
        "usage_bytes": _dir_size(root, exclude=_BACKUP_EXCLUDE_DIRS),
        "paths": [
            entry("数据目录", root),
            entry("会话数据库", sessions_db_path()),
            entry("记忆库", memory_db_path()),
            entry("桌面工作台库", workbench_db_path()),
            entry("LLM 配置", config_path()),
            entry("图表目录", figures_dir()),
            entry("用户产物", output_dir()),
            entry("报告目录", report_dir()),
            entry("缓存目录", cache_dir()),
            entry("日志目录", log_dir()),
        ],
    }


def _resolve_target(target: str) -> Path:
    mapping = {
        "data": user_data_dir,
        "figures": figures_dir,
        "outputs": output_dir,
        "logs": log_dir,
        "cache": cache_dir,
    }
    if target not in mapping:
        raise ValueError("target 只允许 data/figures/outputs/logs/cache")
    return mapping[target]()


# ── 产物目录自定义 ──────────────────────────────────────────

@router.get("/directories")
async def get_directories():
    return {
        "figures_dir": dir_overrides().get("figures_dir", ""),
        "reports_dir": dir_overrides().get("reports_dir", ""),
        "figures_default": str(figures_dir()),
        "reports_default": str(report_dir()),
        "backup_default": str(default_backup_dir()),
    }


@router.post("/directories")
async def update_directories(request: DirectoriesRequest):
    try:
        saved = set_dir_overrides(request.model_dump())
    except ValueError as exc:
        return JSONResponse(status_code=400, content={"error": str(exc)})
    return {"saved": saved}


def _native_directory_picker() -> str | None:
    """系统原生"选择文件夹"对话框（tkinter，stdlib 自带，无需额外依赖）。

    在独立线程中运行（starlette threadpool），用自己的隐藏 Tk root，
    不与运行中的 GUI 工具箱冲突；用户取消返回 None。
    """
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:  # pragma: no cover - 少数精简 Python 无 tkinter
        return None
    try:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            chosen = filedialog.askdirectory(title="选择保存目录")
        finally:
            root.destroy()
        return chosen or None
    except Exception:
        return None


@router.post("/directories/picker")
async def pick_directory():
    from starlette.concurrency import run_in_threadpool

    chosen = await run_in_threadpool(_native_directory_picker)
    if not chosen:
        return JSONResponse(
            status_code=501,
            content={"error": "本机无法打开原生目录选择对话框，请手动输入路径"},
        )
    return {"picked": True, "path": str(Path(chosen))}


def _open_in_file_manager(path: Path) -> None:
    if sys.platform == "win32":
        os.startfile(str(path))  # noqa: S606 - 本地桌面应用的有意行为
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


@router.post("/data/open")
async def open_data_dir(request: OpenDirRequest):
    try:
        path = _resolve_target(request.target)
    except ValueError as exc:
        return JSONResponse(status_code=400, content={"error": str(exc)})
    if not path.exists():
        return JSONResponse(status_code=400, content={"error": f"目录不存在: {path}"})
    _open_in_file_manager(path)
    return {"opened": True, "path": str(path)}


# ── 备份 ────────────────────────────────────────────────────

def _sqlite_snapshot(source: Path, destination: Path) -> bool:
    """用 SQLite backup API 做一致性快照；非 SQLite 文件返回 False。"""
    try:
        with sqlite3.connect(source) as conn, sqlite3.connect(destination) as dest:
            conn.backup(dest)
        return True
    except sqlite3.Error:
        return False


def _iter_backup_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if rel.parts and rel.parts[0] in _BACKUP_EXCLUDE_DIRS:
            continue
        if path.name in _BACKUP_EXCLUDE_NAMES or path.suffix == ".zip":
            continue
        yield path, rel


@router.post("/backup/export")
async def export_backup(request: ExportBackupRequest):
    root = user_data_dir()
    if request.target_dir:
        target = Path(request.target_dir).expanduser()
        if not target.is_absolute():
            return JSONResponse(status_code=400, content={"error": "target_dir 必须是绝对路径"})
        if not target.is_dir():
            return JSONResponse(status_code=400, content={"error": f"目标目录不存在: {target}"})
    else:
        target = default_backup_dir()
        target.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_path = target / f"phytoreason_backup_{stamp}.zip"
    manifest = {"created_at": stamp, "data_dir": str(root), "files": []}

    with tempfile.TemporaryDirectory() as staging_name:
        staging = Path(staging_name)
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for path, rel in _iter_backup_files(root):
                arcname = str(rel).replace("\\", "/")
                if path.suffix in {".db", ".sqlite", ".sqlite3"}:
                    snapshot = staging / rel
                    snapshot.parent.mkdir(parents=True, exist_ok=True)
                    if not _sqlite_snapshot(path, snapshot):
                        shutil.copy2(path, snapshot)
                    zf.write(snapshot, arcname=arcname)
                else:
                    zf.write(path, arcname=arcname)
                manifest["files"].append(arcname)
            zf.writestr("backup_manifest.json", json.dumps(manifest, ensure_ascii=False, indent=2))

    return {"path": str(zip_path), "size_bytes": zip_path.stat().st_size, "files": len(manifest["files"])}


def _validate_zip_entries(zf: zipfile.ZipFile, root: Path) -> dict[str, Path]:
    """安全路径校验：拒绝绝对路径、盘符、`..` 与越界解析（zip-slip）。"""
    mapping: dict[str, Path] = {}
    for info in zf.infolist():
        if info.is_dir():
            continue
        name = info.filename.replace("\\", "/")
        if not name or name.startswith("/") or ":" in name.split("/")[0]:
            raise ValueError(f"拒绝绝对路径条目: {info.filename}")
        parts = [part for part in name.split("/") if part not in ("", ".")]
        if any(part == ".." for part in parts):
            raise ValueError(f"拒绝相对路径条目: {info.filename}")
        resolved = (root / Path(*parts)).resolve()
        if not str(resolved).startswith(str(root.resolve())):
            raise ValueError(f"条目越界: {info.filename}")
        mapping[name] = resolved
    return mapping


@router.post("/backup/import")
async def import_backup(file: UploadFile = File(...)):
    payload = await file.read()
    if len(payload) > MAX_IMPORT_BYTES:
        return JSONResponse(status_code=400, content={"error": "备份文件超过 1GB 上传上限"})
    root = user_data_dir().resolve()

    try:
        zf = zipfile.ZipFile(io.BytesIO(payload))
    except zipfile.BadZipFile:
        return JSONResponse(status_code=400, content={"error": "不是有效的 zip 备份文件"})

    try:
        targets = _validate_zip_entries(zf, root)
    except ValueError as exc:
        zf.close()
        return JSONResponse(status_code=400, content={"error": str(exc)})

    total_uncompressed = sum(info.file_size for info in zf.infolist() if not info.is_dir())
    if total_uncompressed > MAX_IMPORT_BYTES:
        zf.close()
        return JSONResponse(status_code=400, content={"error": f"解压后总量 {total_uncompressed} 字节超过 1GB 上限"})

    restored = 0
    backed_up: list[str] = []
    try:
        for name, target in targets.items():
            if name == "backup_manifest.json":
                continue
            if target.exists():
                bak = Path(str(target) + ".bak")
                if bak.exists():
                    bak.unlink()  # 只保留一代 .bak，旧 .bak 让位于本次
                target.rename(bak)
                backed_up.append(str(target.relative_to(root)))
            target.parent.mkdir(parents=True, exist_ok=True)
            with zf.open(name) as source, open(target, "wb") as dest:
                shutil.copyfileobj(source, dest)
            restored += 1
    except OSError as exc:
        return JSONResponse(status_code=500, content={"error": f"导入失败: {exc}", "restored": restored})
    finally:
        zf.close()

    return {"restored": restored, "backed_up": backed_up, "data_dir": str(root)}
