"""
tool_executor.py — 工具执行引擎。

统一负责:
  - 工具查找
  - 输入验证
  - 执行
  - 异常处理
  - 日志记录
"""

from __future__ import annotations

import logging
import traceback
from datetime import datetime
from typing import Any

from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.tool_registry import TOOL_REGISTRY

logger = logging.getLogger("tool_executor")


class ToolExecutionError(Exception):
    """工具执行错误。"""
    pass


class ToolExecutor:
    """工具执行器。"""

    def __init__(self, registry: Any | None = None) -> None:
        self.registry = registry or TOOL_REGISTRY

    def execute(
        self,
        tool_name: str,
        raise_on_error: bool = False,
        log_execution: bool = True,
        **kwargs,
    ) -> ToolResult:
        """查找并执行指定工具。

        Args:
            tool_name: 工具名称 (在 registry 中注册的名称)
            raise_on_error: 是否在失败时抛出异常
            log_execution: 是否记录执行日志
            **kwargs: 传递给工具 run() 的参数

        Returns:
            ToolResult: 工具执行结果 (失败时返回含错误信息的 ToolResult)
        """
        tool = self.registry.get(tool_name)
        if tool is None:
            msg = f"未找到工具: '{tool_name}'. 可用工具: {self.registry.list_tools()}"
            if log_execution:
                logger.error(msg)
            if raise_on_error:
                raise ToolExecutionError(msg)
            return ToolResult(
                candidates=[],
                evidence_list=[],
                warnings=[msg],
                metadata={"tool_name": tool_name, "status": "not_found"},
            )

        # 输入验证
        errors = tool.validate_input(**kwargs)
        if errors:
            msg = f"工具 '{tool_name}' 输入验证失败: {'; '.join(errors)}"
            if log_execution:
                logger.warning(msg)
            if raise_on_error:
                raise ToolExecutionError(msg)
            return ToolResult(
                candidates=[],
                evidence_list=[],
                warnings=[msg],
                metadata={"tool_name": tool_name, "status": "validation_failed"},
            )

        # 缓存检查
        from phyto_reason.utils.cache import tool_cache, tool_key

        cache_ttl = _get_tool_cache_ttl(tool_name)
        ck = tool_key(tool_name, **kwargs)
        if cache_ttl > 0:
            cached = tool_cache.get(ck)
            if cached is not None:
                if log_execution:
                    logger.info(f"工具 '{tool_name}' 命中缓存, 命中率={tool_cache.stats()['hit_rate']:.0%}")
                return cached

        # 执行
        if log_execution:
            logger.info(f"执行工具: {tool_name}, 参数: {kwargs}")

        try:
            result = tool.run(**kwargs)
            result.metadata["tool_name"] = tool_name
            result.metadata["status"] = "completed"
            result.metadata["tool_version"] = tool.version

            # 缓存结果
            if cache_ttl > 0 and result.status == "success":
                tool_cache.set(ck, result, ttl=cache_ttl)

            if log_execution:
                logger.info(f"工具 '{tool_name}' 执行完成: {result.to_summary()}")
            return result

        except Exception as e:
            tb = traceback.format_exc()
            msg = f"工具 '{tool_name}' 执行失败: {e}"
            if log_execution:
                logger.error(f"{msg}\n{tb}")
            if raise_on_error:
                raise ToolExecutionError(msg) from e
            return ToolResult(
                candidates=[],
                evidence_list=[],
                warnings=[msg],
                metadata={
                    "tool_name": tool_name,
                    "status": "failed",
                    "error": str(e),
                },
            )


def _get_tool_cache_ttl(tool_name: str) -> int:
    """Return cache TTL in seconds for a given tool.

    0 = disable caching for this tool.
    """
    ttl_map: dict[str, int] = {
        # API calls — cache 1 hour
        "kegg_pathway": 3600,
        "literature_search": 3600,
        # Heavy computation — cache 30 min
        "correlation": 1800,
        "wgcna": 1800,
        "wgcna_real": 1800,
        "motif_scan": 1800,
        "motif_real": 1800,
        "deg_analysis": 1800,
        "differential_analysis": 1800,
        # Lightweight / semi-static — cache 1 hour
        "tf_annotation": 3600,
    }
    return ttl_map.get(tool_name, 0)
