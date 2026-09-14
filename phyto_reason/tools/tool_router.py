"""
tool_router.py — 动态工具路由器。

从 registry 中查找工具，支持别名和错误处理。
"""

from __future__ import annotations

from phyto_reason.tools.tool_registry import TOOL_REGISTRY


TOOL_ALIASES: dict[str, str] = {
    "correlation": "correlation",
    "corr": "correlation",
    "pearson": "correlation",
    "wgcna": "wgcna",
    "network": "wgcna",
    "coexpression": "wgcna",
    "motif": "motif_scan",
    "motif_scan": "motif_scan",
    "promoter": "motif_scan",
}


class ToolRouter:
    """工具路由器。"""

    def __init__(self, registry=None) -> None:
        self.registry = registry or TOOL_REGISTRY

    def resolve(self, tool_name: str):
        """将 tool_name 解析为工具实例。支持别名。"""
        resolved = TOOL_ALIASES.get(tool_name.lower(), tool_name)
        tool = self.registry.get(resolved)
        return tool

    def list_available(self) -> list[str]:
        """列出所有已注册的工具名称 (含别名)。"""
        return self.registry.list_tools()

    def list_with_aliases(self) -> dict[str, str]:
        """列出所有工具及其别名。"""
        result: dict[str, str] = {}
        for tool_name in self.registry.list_tools():
            result[tool_name] = tool_name
        for alias, target in TOOL_ALIASES.items():
            if alias != target:
                result[alias] = target
        return result
