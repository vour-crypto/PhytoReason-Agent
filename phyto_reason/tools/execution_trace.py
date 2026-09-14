"""
execution_trace.py — 工具调用跟踪记录。
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, Field


class CallRecord(BaseModel):
    call_id: str = ""
    tool_name: str = ""
    caller: str = ""
    arguments: dict = Field(default_factory=dict)
    started_at: str = ""
    finished_at: str = ""
    duration: float = 0.0
    success: bool = False
    result_summary: dict = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)
    metadata: dict = Field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "call_id": self.call_id,
            "tool": self.tool_name,
            "caller": self.caller,
            "duration": round(self.duration, 3),
            "success": self.success,
            "n_warnings": len(self.warnings),
            "n_errors": len(self.errors),
        }


class ExecutionTrace:
    """工具调用追踪器。"""

    def __init__(self) -> None:
        self.records: list[CallRecord] = []

    def start_call(self, call_id: str, tool_name: str, caller: str = "",
                   arguments: dict | None = None) -> CallRecord:
        record = CallRecord(
            call_id=call_id,
            tool_name=tool_name,
            caller=caller,
            arguments=arguments or {},
            started_at=datetime.now().isoformat(),
        )
        self.records.append(record)
        return record

    def finish_call(self, record: CallRecord, success: bool,
                    result_summary: dict | None = None,
                    warnings: list[str] | None = None,
                    errors: list[str] | None = None) -> CallRecord:
        record.finished_at = datetime.now().isoformat()
        from datetime import datetime as dt
        try:
            start = dt.fromisoformat(record.started_at)
            end = dt.fromisoformat(record.finished_at)
            record.duration = (end - start).total_seconds()
        except Exception:
            record.duration = 0.0
        record.success = success
        record.result_summary = result_summary or {}
        record.warnings = warnings or []
        record.errors = errors or []
        return record

    def get_history(self, tool_name: str | None = None, limit: int = 10) -> list[CallRecord]:
        filtered = self.records if tool_name is None else [
            r for r in self.records if r.tool_name == tool_name
        ]
        return filtered[-limit:]

    def summary(self) -> dict:
        n_total = len(self.records)
        n_success = sum(1 for r in self.records if r.success)
        return {
            "total_calls": n_total,
            "successful": n_success,
            "failed": n_total - n_success,
            "tools_used": list(set(r.tool_name for r in self.records)),
        }
