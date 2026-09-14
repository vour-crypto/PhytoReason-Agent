"""
base_tool.py — 统一科研工具抽象基类。

所有分析工具必须继承 BaseTool 并实现 run() 方法。
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from phyto_reason.models.tool_result import ToolResult


class ToolParameter:
    """工具输入参数定义。"""
    def __init__(self, name: str, type: str = "string", description: str = "",
                 required: bool = False, default: Any = None) -> None:
        self.name = name
        self.type = type
        self.description = description
        self.required = required
        self.default = default


class BaseTool(ABC):
    tool_name: str = ""
    description: str = ""
    version: str = "1.0.0"

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if not hasattr(cls, "parameters") or cls.parameters is None:
            cls.parameters = []
        if not hasattr(cls, "supported_data_types") or cls.supported_data_types is None:
            cls.supported_data_types = []

    @abstractmethod
    def validate_input(self, **kwargs) -> list[str]:
        """验证输入参数，返回错误信息列表。空列表=验证通过。"""
        ...

    @abstractmethod
    def run(self, **kwargs) -> ToolResult:
        """执行工具分析，返回统一 ToolResult。"""
        ...

    def summarize_result(self, result: ToolResult) -> dict:
        """生成结果摘要。"""
        return result.to_summary()

    def get_parameter_schema(self) -> dict:
        """返回 OpenAI/Claude function calling 兼容的 schema。"""
        properties = {}
        required = []
        params = self.parameters or []
        for p in params:
            properties[p.name] = {
                "type": p.type,
                "description": p.description,
            }
            if p.required:
                required.append(p.name)
        return {
            "name": self.tool_name,
            "description": self.description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": required,
            },
        }
