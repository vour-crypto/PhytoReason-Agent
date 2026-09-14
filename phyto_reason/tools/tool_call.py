"""
tool_call.py — 统一工具调用协议对象。
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, Field


class ToolCall(BaseModel):
    tool_name: str
    arguments: dict[str, Any] = Field(default_factory=dict)
    caller: str = "unknown"
    call_id: str = Field(default_factory=lambda: datetime.now().strftime("call_%Y%m%d_%H%M%S_%f"))
    timestamp: str = Field(default_factory=lambda: datetime.now().isoformat())

    def to_log(self) -> dict:
        return {
            "call_id": self.call_id,
            "tool": self.tool_name,
            "caller": self.caller,
            "n_args": len(self.arguments),
            "arg_keys": list(self.arguments.keys()),
        }
