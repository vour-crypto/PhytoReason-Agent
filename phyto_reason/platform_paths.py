"""Per-user runtime paths and one-time migration for PhytoReason."""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

# 允许用户自定义的产物目录（directories.json，绝对路径覆盖默认位置）
_DIR_OVERRIDE_KEYS = ("figures_dir", "reports_dir")


def user_data_dir() -> Path:
    override = os.getenv("PHYTOREASON_DATA_DIR")
    if override:
        root = Path(override).expanduser()
    elif sys.platform == "win32":
        root = Path(os.getenv("LOCALAPPDATA", Path.home() / "AppData" / "Local")) / "PhytoReason"
    elif sys.platform == "darwin":
        root = Path.home() / "Library" / "Application Support" / "PhytoReason"
    else:
        root = Path(os.getenv("XDG_DATA_HOME", Path.home() / ".local" / "share")) / "PhytoReason"
    candidates = [root, Path.home() / ".phytoreason", Path(tempfile.gettempdir()) / "PhytoReason"]
    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            probe = candidate / ".write_test"
            probe.write_text("ok", encoding="ascii")
            probe.unlink()
            root = candidate
            break
        except OSError:
            continue
    else:
        raise OSError("无法创建 PhytoReason 用户数据目录")
    _migrate_once(root)
    return root


def directories_config_path() -> Path:
    return user_data_dir() / "directories.json"


def dir_overrides() -> dict[str, str]:
    """当前生效的目录覆盖（仅含已设置的键）。"""
    path = directories_config_path()
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(value, dict):
        return {}
    return {k: str(v) for k, v in value.items() if k in _DIR_OVERRIDE_KEYS and v}


def set_dir_overrides(payload: dict) -> dict[str, str]:
    """保存产物目录覆盖；空值=清除该键的覆盖，回落默认位置。

    校验绝对路径且可创建；figures_dir 的变更需重启应用生效（静态挂载在启动时建立）。
    """
    cleaned: dict[str, str] = {}
    for key in _DIR_OVERRIDE_KEYS:
        raw = str(payload.get(key, "") or "").strip()
        if not raw:
            continue
        target = Path(raw).expanduser()
        if not target.is_absolute():
            raise ValueError(f"{key} 必须是绝对路径: {raw}")
        try:
            target.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise ValueError(f"{key} 目录无法创建: {raw}（{exc}）") from exc
        cleaned[key] = str(target)
    path = directories_config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cleaned, ensure_ascii=False, indent=2), encoding="utf-8")
    return cleaned


def _override_or_default(key: str, default: Path) -> Path:
    override = dir_overrides().get(key, "")
    if not override:
        return default
    path = Path(override)
    try:
        path.mkdir(parents=True, exist_ok=True)
        return path
    except OSError:
        return default  # 自定义目录不可用时回落默认位置


def figures_dir() -> Path:
    return _override_or_default("figures_dir", user_data_dir() / "figures")


def sessions_db_path() -> Path:
    path = user_data_dir() / "data" / "sessions.db"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def memory_db_path() -> Path:
    """Agent 记忆库路径。

    记忆表目前与会话同库（sessions.db）——本函数用于把路径合同
    显式化并可测试；将来拆库时只需改此处，调用方不动。
    """
    return sessions_db_path()


def workbench_db_path() -> Path:
    path = user_data_dir() / "data" / "workbench.db"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def config_path() -> Path:
    return user_data_dir() / "config.json"


def cache_dir() -> Path:
    path = user_data_dir() / "cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def output_dir() -> Path:
    path = user_data_dir() / "data" / "outputs"
    path.mkdir(parents=True, exist_ok=True)
    return path


def report_dir() -> Path:
    return _override_or_default("reports_dir", output_dir() / "reports")


def log_dir() -> Path:
    path = user_data_dir() / "logs"
    path.mkdir(parents=True, exist_ok=True)
    return path


def copy_if_missing(source: Path, target: Path) -> None:
    """迁移原语：目标已存在则绝不覆盖（copy-if-missing 语义）。

    公开给测试使用，保证迁移幂等性有可断言的合同。
    """
    if not source.exists() or target.exists():
        return
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        return
    if source.is_dir():
        try:
            shutil.copytree(source, target)
        except OSError:
            return
    else:
        try:
            shutil.copy2(source, target)
        except OSError:
            return


def _migrate_once(root: Path) -> None:
    marker = root / ".migration_v1.complete"
    if marker.exists():
        return
    package_root = Path(__file__).resolve().parent
    old_home = Path.home() / ".phytoreason"
    copy_if_missing(old_home / "workbench.db", root / "data" / "workbench.db")
    copy_if_missing(package_root.parent / "data" / "sessions.db", root / "data" / "sessions.db")
    copy_if_missing(package_root / "api" / "static" / "figures", root / "figures")
    copy_if_missing(package_root.parent / "logs", root / "logs")
    copy_if_missing(package_root.parent / "outputs", root / "cache")
    try:
        marker.write_text("migration_v1_complete\n", encoding="utf-8")
    except OSError:
        pass


def default_backup_dir() -> Path:
    """备份导出默认目录：优先系统"下载"文件夹（用户最熟悉），否则用户数据目录/backups。"""
    downloads = Path.home() / "Downloads"
    try:
        if downloads.is_dir():
            return downloads
    except OSError:
        pass
    return user_data_dir() / "backups"
