"""phyto_reason.tools — 统一科研工具系统。"""

# ── 工具抽象层 ──────────────────────────────────────────
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import TOOL_REGISTRY, ToolRegistry, register_tool
from phyto_reason.tools.tool_executor import ToolExecutor, ToolExecutionError

# ── Function Calling Protocol (展平自 tools/calling/) ────
from phyto_reason.tools.tool_call import ToolCall
from phyto_reason.tools.tool_response import ToolResponse
from phyto_reason.tools.function_schema import FunctionSchema, generate_all_schemas
from phyto_reason.tools.argument_validator import ArgumentValidator, ValidationResult
from phyto_reason.tools.tool_router import ToolRouter
from phyto_reason.tools.execution_trace import ExecutionTrace, CallRecord
from phyto_reason.tools.calling_engine import CallingEngine

# ── 活跃查询工具 ─────────────────────────────────────────
from phyto_reason.tools.literature import PubMedSearch
from phyto_reason.tools.rag import RAGEngine

# 注册新工具 (import 触发装饰器)
from phyto_reason.tools.network import correlation_tool   # noqa: F401
from phyto_reason.tools.network import wgcna_tool         # noqa: F401
from phyto_reason.tools.motif import motif_tool            # noqa: F401
import phyto_reason.tools.network.wgcna_r              # noqa: F401
import phyto_reason.tools.motif.fimo_wrapper           # noqa: F401
import phyto_reason.tools.pathway.kegg_client          # noqa: F401
from phyto_reason.tools.literature.pubmed_tool import PubMedTool  # noqa: F401
from phyto_reason.tools.annotation.tf_annotation_tool import TFAnnotationTool  # noqa: F401
from phyto_reason.tools.statistics.deg_tool import DEGTool  # noqa: F401
from phyto_reason.tools.statistics.differential_tool import DifferentialTool  # noqa: F401
import phyto_reason.tools.phase2_tools  # noqa: F401
import phyto_reason.tools.spectra_tools  # noqa: F401

__all__ = [
    # 抽象层
    "BaseTool", "ToolParameter", "TOOL_REGISTRY", "ToolRegistry",
    "register_tool", "ToolExecutor", "ToolExecutionError",
    # Function Calling
    "ToolCall", "ToolResponse", "FunctionSchema", "generate_all_schemas",
    "ArgumentValidator", "ValidationResult", "ToolRouter",
    "ExecutionTrace", "CallRecord", "CallingEngine",
    # 查询工具
    "PubMedSearch", "RAGEngine",
]
