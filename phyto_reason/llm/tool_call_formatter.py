"""
tool_call_formatter.py — LLM tool calling 格式化器。

将 BaseTool 的 function schema 转换为 LLM-compatible format，
并解析 LLM 返回的 tool calls。
"""

from __future__ import annotations

import json
from typing import Any

from phyto_reason.tools.tool_registry import TOOL_REGISTRY


class ToolCallFormatter:
    """Tool calling 格式化器。"""

    @staticmethod
    def get_llm_tools(tool_names: list[str] | None = None) -> list[dict]:
        """获取 LLM-compatible tool schemas。"""
        if tool_names:
            tools = []
            for name in tool_names:
                tool = TOOL_REGISTRY.get(name)
                if tool:
                    tools.append(ToolCallFormatter._to_llm_schema(tool))
            return tools

        return [
            ToolCallFormatter._to_llm_schema(t)
            for t_name in TOOL_REGISTRY.list_tools()
            if (t := TOOL_REGISTRY.get(t_name)) is not None
        ]

    @staticmethod
    def _to_llm_schema(tool) -> dict:
        """将 BaseTool 转换为 OpenAI tool schema。"""
        params = tool.parameters or []
        properties = {}
        required = []
        for p in params:
            properties[p.name] = {
                "type": p.type,
                "description": p.description,
            }
            if p.required:
                required.append(p.name)
        return {
            "type": "function",
            "function": {
                "name": tool.tool_name,
                "description": tool.description,
                "parameters": {
                    "type": "object",
                    "properties": properties,
                    "required": required,
                },
            },
        }

    @staticmethod
    def parse_tool_call(tool_call: dict) -> tuple[str | None, dict]:
        """解析 LLM 返回的 tool call。

        Args:
            tool_call: LLM tool_call dict

        Returns:
            (tool_name, arguments)
        """
        func = tool_call.get("function", tool_call)
        name = func.get("name")
        args = func.get("arguments", {})
        if isinstance(args, str):
            try:
                args = json.loads(args)
            except json.JSONDecodeError:
                args = {"_raw": args}
        return name, args

    @staticmethod
    def format_tool_results(tool_name: str, result) -> str:
        """将 ToolResult 格式化为 LLM-readable string。"""
        summary = result.to_summary() if hasattr(result, 'to_summary') else {}
        return json.dumps({
            "tool": tool_name,
            "candidates_count": summary.get("n_candidates", 0),
            "evidences_count": summary.get("n_evidences", 0),
            "warnings": result.warnings if hasattr(result, 'warnings') else [],
            "details": summary,
        }, ensure_ascii=False)
