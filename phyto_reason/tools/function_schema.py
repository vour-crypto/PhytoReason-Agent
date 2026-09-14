"""
function_schema.py — OpenAI-compatible function schema。

从 BaseTool 自动生成 OpenAI/DeepSeek function calling schema。
"""

from __future__ import annotations

from phyto_reason.tools.tool_registry import TOOL_REGISTRY


class FunctionSchema:
    """工具函数 schema 生成器。"""

    TYPE_MAP = {
        "string": "string",
        "number": "number",
        "integer": "integer",
        "boolean": "boolean",
        "array": "array",
        "object": "object",
    }

    @staticmethod
    def from_tool(tool) -> dict:
        """从 BaseTool 实例生成 OpenAI-compatible function schema。"""
        params = tool.parameters or []
        properties = {}
        required = []

        for p in params:
            json_type = FunctionSchema.TYPE_MAP.get(p.type, "string")
            prop = {"type": json_type, "description": p.description}
            if p.default is not None:
                pass
            properties[p.name] = prop
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
    def from_tool_simple(tool) -> dict:
        """简化的 schema (不嵌套 function 字段)。"""
        params = tool.parameters or []
        properties = {}
        required = []
        for p in params:
            properties[p.name] = {"type": p.type, "description": p.description}
            if p.required:
                required.append(p.name)
        return {
            "name": tool.tool_name,
            "description": tool.description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": required,
            },
        }


def generate_all_schemas() -> list[dict]:
    """生成所有已注册工具的 OpenAI-compatible schemas。"""
    tools = [
        TOOL_REGISTRY.get(name) for name in TOOL_REGISTRY.list_tools()
        if TOOL_REGISTRY.get(name) is not None
    ]
    return [FunctionSchema.from_tool(t) for t in tools]


def generate_all_simple() -> list[dict]:
    """生成简化版本的所有工具 schema。"""
    tools = [
        TOOL_REGISTRY.get(name) for name in TOOL_REGISTRY.list_tools()
        if TOOL_REGISTRY.get(name) is not None
    ]
    return [FunctionSchema.from_tool_simple(t) for t in tools]
