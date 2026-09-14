"""
tool_registry.py — 工具注册中心。

所有工具在模块加载时自动注册到 TOOL_REGISTRY。
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from phyto_reason.tools.base_tool import BaseTool


class ToolRegistry:
    """工具注册中心，管理所有注册的工具。"""

    def __init__(self) -> None:
        self._tools: dict[str, BaseTool] = {}

    def register(self, tool: BaseTool) -> None:
        """注册一个工具实例。"""
        if not tool.tool_name:
            raise ValueError(f"Tool must have a non-empty tool_name: {type(tool).__name__}")
        if tool.tool_name in self._tools:
            from warnings import warn
            warn(f"Overwriting existing tool: {tool.tool_name}")
        self._tools[tool.tool_name] = tool

    def get(self, tool_name: str) -> BaseTool | None:
        """通过名称获取工具。"""
        return self._tools.get(tool_name)

    def list_tools(self) -> list[str]:
        """列出所有已注册的工具名称。"""
        return sorted(self._tools.keys())

    def list_by_data_type(self, data_type: str) -> list[BaseTool]:
        """筛选支持指定数据类型的工具。"""
        return [
            t for t in self._tools.values()
            if data_type in t.supported_data_types
        ]

    @property
    def count(self) -> int:
        return len(self._tools)

    def get_tool_schemas(self) -> list[dict]:
        """返回所有工具的 function calling schema。"""
        return [t.get_parameter_schema() for t in self._tools.values()]


# 全局单例
TOOL_REGISTRY = ToolRegistry()


def register_tool(tool_cls):
    """装饰器: 实例化工具类并注册到全局 registry。"""
    instance = tool_cls()
    TOOL_REGISTRY.register(instance)
    return tool_cls
