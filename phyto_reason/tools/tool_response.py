"""
tool_response.py — 统一工具调用响应对象。
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.execution_trace import CallRecord


class ToolResponse(BaseModel):
    success: bool = False
    tool_name: str = ""
    result: ToolResult | None = None
    content: str = ""
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)
    duration: float = 0.0
    trace: CallRecord | None = None

    def to_llm_message(self) -> str:
        lines = [f"工具: {self.tool_name}"]
        if self.success:
            lines.append(f"状态: 成功 ({self.duration:.2f}s)")
            if self.content:
                lines.append(self.content)
            if self.result:
                summary = self.result.to_summary()
                lines.append(f"输出: {summary}")
        else:
            lines.append(f"状态: 失败")
            if self.errors:
                lines.extend([f"  - {e}" for e in self.errors])
        if self.warnings:
            lines.append(f"警告: {'; '.join(self.warnings)}")
        return "\n".join(lines)

    def to_summary(self) -> dict:
        return {
            "success": self.success,
            "tool": self.tool_name,
            "duration": round(self.duration, 3),
            "n_warnings": len(self.warnings),
            "n_errors": len(self.errors),
            "has_result": self.result is not None,
        }
