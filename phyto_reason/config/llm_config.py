"""Persistent, masked multi-provider LLM configuration."""

from __future__ import annotations

import json
import os
import stat
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

from phyto_reason.platform_paths import config_path

PROVIDERS = {"deepseek", "openai", "custom"}
DEFAULTS = {"provider": "deepseek", "base_url": "https://api.deepseek.com/v1", "model": "deepseek-chat", "api_key": ""}


def _read_file() -> dict:
    path = config_path()
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def file_config() -> dict:
    """Raw config.json content (no env merge). Public accessor for callers
    that must implement their own precedence, e.g. the desktop shell."""
    return _read_file()


def load_config() -> dict:
    file_config = _read_file()
    result = dict(DEFAULTS)
    result.update({
        "api_key": os.getenv("OPENAI_API_KEY", DEFAULTS["api_key"]),
        "base_url": os.getenv("OPENAI_BASE_URL", DEFAULTS["base_url"]),
        "model": os.getenv("LLM_MODEL", DEFAULTS["model"]),
    })
    result.update({key: value for key, value in file_config.items() if value is not None})
    return result


def config_fingerprint() -> tuple[int, int]:
    path = config_path()
    try:
        info = path.stat()
        return info.st_mtime_ns, info.st_size
    except OSError:
        return 0, 0


def validate_config(payload: dict) -> dict:
    provider = str(payload.get("provider", "deepseek")).lower()
    if provider not in PROVIDERS:
        raise ValueError("provider 必须是 deepseek、openai 或 custom")
    base_url = str(payload.get("base_url", DEFAULTS["base_url"])).strip().rstrip("/")
    parsed = urlparse(base_url)
    if parsed.scheme not in {"https", "http"} or not parsed.netloc:
        raise ValueError("Base URL 只允许 https 或 localhost 的 http")
    if parsed.scheme == "http" and parsed.hostname not in {"localhost", "127.0.0.1", "::1"}:
        raise ValueError("http Base URL 只允许 localhost")
    model = str(payload.get("model", "")).strip()
    if not model:
        raise ValueError("模型名不能为空")
    return {"provider": provider, "base_url": base_url, "model": model, "api_key": str(payload.get("api_key", ""))}


def save_config(payload: dict) -> dict:
    current = _read_file()
    incoming = dict(payload)
    if not incoming.get("api_key") and current.get("api_key"):
        incoming["api_key"] = current["api_key"]
    normalized = validate_config(incoming)
    for key in ("probe", "updated_at"):
        if key in current:
            normalized[key] = current[key]
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(normalized, ensure_ascii=False, indent=2), encoding="utf-8")
    try:
        path.chmod(stat.S_IRUSR | stat.S_IWUSR)
    except OSError:
        pass
    return normalized


def mask_key(value: str) -> str:
    if not value:
        return ""
    return value[:3] + "***" + value[-4:] if len(value) > 7 else "***"


def public_config() -> dict:
    value = load_config()
    return {
        "provider": value["provider"],
        "base_url": value["base_url"],
        "model": value["model"],
        "api_key": mask_key(value.get("api_key", "")),
        "probe": value.get("probe"),
    }


def record_probe(result: dict) -> None:
    """Persist the probe result without ever writing env-derived secrets.

    只有 config.json 中已存在的内容会被保留；环境变量里的 key/base_url
    绝不因此落盘（Phase 6.4 Step 2 挂账修复）。
    """
    value = _read_file()
    value["probe"] = {**result, "timestamp": datetime.now(timezone.utc).isoformat()}
    path = config_path()
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    try:
        path.chmod(stat.S_IRUSR | stat.S_IWUSR)
    except OSError:
        pass


def probe_config() -> dict:
    """Probe authentication and function calling without exposing the key."""
    value = load_config()
    result = {"connectivity": False, "supports_tools": False, "error": None, "provider": value["provider"]}
    if not value.get("api_key"):
        result["error"] = "未配置 API key"
        record_probe(result)
        return result
    try:
        from openai import OpenAI
        client = OpenAI(api_key=value["api_key"], base_url=value["base_url"], timeout=10.0)
        client.models.list()
        result["connectivity"] = True
        try:
            client.chat.completions.create(
                model=value["model"],
                messages=[{"role": "user", "content": "ping"}],
                tools=[{"type": "function", "function": {"name": "probe", "description": "连通性探测", "parameters": {"type": "object", "properties": {}}}}],
                max_tokens=1,
                timeout=10.0,
            )
            result["supports_tools"] = True
        except Exception as exc:
            result["error"] = "连通成功，但工具调用探测失败: " + type(exc).__name__
    except Exception as exc:
        result["error"] = type(exc).__name__ + ": " + str(exc)[:200]
    record_probe(result)
    return result
