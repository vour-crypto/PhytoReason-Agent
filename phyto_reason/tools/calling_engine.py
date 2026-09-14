"""
calling_engine.py — 统一工具调用引擎。

整合 validation → routing → execution → tracing。
"""

from __future__ import annotations

import logging
import time
from datetime import datetime
from typing import Any

from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.tool_call import ToolCall
from phyto_reason.tools.tool_response import ToolResponse
from phyto_reason.tools.argument_validator import ArgumentValidator
from phyto_reason.tools.tool_router import ToolRouter
from phyto_reason.tools.execution_trace import CallRecord, ExecutionTrace
from phyto_reason.tools.tool_executor import ToolExecutor as OldToolExecutor

logger = logging.getLogger("calling_engine")


class CallingEngine:
    """统一工具调用引擎。

    调用流程:
      1. 解析 ToolCall
      2. 验证参数
      3. 路由到 tool
      4. 执行
      5. 追踪
      6. 返回 ToolResponse
    """

    def __init__(self, registry=None) -> None:
        self.router = ToolRouter(registry)
        self.validator = ArgumentValidator()
        self.executor = OldToolExecutor(registry)
        self.trace = ExecutionTrace()

    def execute(self, call: ToolCall) -> ToolResponse:
        """执行一个 ToolCall。"""
        start = time.time()
        record = self.trace.start_call(
            call_id=call.call_id,
            tool_name=call.tool_name,
            caller=call.caller,
            arguments=call.arguments,
        )

        # 1. 路由到 tool
        tool = self.router.resolve(call.tool_name)
        if tool is None:
            errors = [f"未找到工具: '{call.tool_name}'. 可用: {self.router.list_available()}"]
            record = self.trace.finish_call(record, success=False, errors=errors)
            return ToolResponse(
                success=False, tool_name=call.tool_name,
                errors=errors, duration=time.time() - start, trace=record,
            )

        # 2. 参数验证
        validation = self.validator.validate(tool, call.arguments)
        if not validation.valid:
            record = self.trace.finish_call(
                record, success=False,
                errors=validation.errors, warnings=validation.warnings,
            )
            return ToolResponse(
                success=False, tool_name=call.tool_name,
                errors=validation.errors, warnings=validation.warnings,
                duration=time.time() - start, trace=record,
            )

        # 3. 执行
        try:
            result: ToolResult = tool.run(**call.arguments)
            elapsed = time.time() - start
            content = self._build_content(tool.tool_name, result)
            record = self.trace.finish_call(
                record, success=True,
                result_summary=result.to_summary(),
                warnings=result.warnings,
            )
            return ToolResponse(
                success=True, tool_name=call.tool_name,
                result=result, content=content,
                warnings=result.warnings,
                duration=elapsed, trace=record,
            )
        except Exception as e:
            elapsed = time.time() - start
            err_msg = f"执行失败: {type(e).__name__}: {e}"
            logger.error(err_msg)
            record = self.trace.finish_call(record, success=False, errors=[err_msg])
            return ToolResponse(
                success=False, tool_name=call.tool_name,
                errors=[err_msg], duration=elapsed, trace=record,
            )

    def execute_py(self, tool_name: str, caller: str = "python",
                   **kwargs) -> ToolResponse:
        """Python 直接调用快捷方法。"""
        call = ToolCall(tool_name=tool_name, arguments=kwargs, caller=caller)
        return self.execute(call)

    def get_trace_history(self, tool_name: str | None = None,
                          limit: int = 10) -> list[CallRecord]:
        return self.trace.get_history(tool_name, limit)

    def summary(self) -> dict:
        return self.trace.summary()

    @staticmethod
    def _build_content(tool_name: str, result: ToolResult) -> str:
        """从 ToolResult 构建人类可读的内容。"""
        parts = [f"工具 '{tool_name}' 执行完成。"]
        summary = result.to_summary()
        parts.append(f"  候选: {summary['n_candidates']} 个")
        parts.append(f"  证据: {summary['n_evidences']} 条")
        if result.warnings:
            parts.append(f"  警告: {'; '.join(result.warnings[:3])}")
        return "\n".join(parts)
