"""
routing_logic.py -- 科学门控路由 (v3.1).

所有路由函数集中于此。
graph 不包含任何科学逻辑 -- 仅做 orchestration。

路由函数回答：
  - 当前验证了什么？
  - 是否能继续下一步？
  - 如果不能，应该输出什么？

v3.1: + 分析步骤选择 (STEP_TO_NODES, STEP_OPTIONS, expand_analysis_selection)
"""

from __future__ import annotations

from langgraph.graph import END

from phyto_reason.workflows.runtime_state import RuntimeState


# ── 分析步骤选择: 用户可见步骤 -> 图节点映射 ──────────────────────

STEP_TO_NODES: dict[str, set[str]] = {
    "data_qc": {"data_qc"},
    "DEG": {"deg_analysis"},
    "DAM": {"dam_analysis"},
    "multiomics": {"deg_analysis", "dam_analysis", "multiomics_integration",
                   "joint_enrichment", "quadrant_plot", "correlation_network",
                   "wgcna", "o2pls"},
    "joint_enrichment": {"joint_enrichment"},
    "quadrant_plot": {"quadrant_plot"},
    "correlation_network": {"correlation_network"},
    "wgcna": {"wgcna"},
    "o2pls": {"o2pls"},
    "tf_narrowing": {"tf_narrowing"},
    "hypothesis": {"regulation_evidence", "contradiction_check", "falsification",
                   "hypothesis_competition", "hypothesis_synthesis"},
}

"""始终运行的 4 个路由关键节点 -- 保证图的条件路由不中断。"""
ALWAYS_RUN_NODES: set[str] = {
    "data_quality_gate", "metabolite_validation", "pathway_coherence",
    "mechanism_classifier",
}

"""全部 19 个管线节点。"""
ALL_PIPELINE_NODES: set[str] = ALWAYS_RUN_NODES | {
    "data_qc", "deg_analysis", "dam_analysis", "multiomics_integration",
    "joint_enrichment", "quadrant_plot", "correlation_network",
    "wgcna", "o2pls",
    "tf_narrowing", "regulation_evidence", "contradiction_check",
    "falsification", "hypothesis_competition", "hypothesis_synthesis",
}

STEP_OPTIONS: list[dict] = [
    {"id": "data_qc", "label": "数据质控 (QC)", "label_en": "Data Quality Control",
     "description": "样本相关性热图 + PCA + 聚类热图 + UpSet 基因重叠图"},
    {"id": "DEG", "label": "差异表达分析 (DEG)", "label_en": "Differential Expression Analysis",
     "description": "鉴定不同组织/处理间差异表达基因，含 limma 风格 moderated t-test + 效应量排名"},
    {"id": "DAM", "label": "差异代谢物分析 (DAM)", "label_en": "Differential Metabolite Analysis",
     "description": "鉴定差异积累代谢物，ANOVA + 效应量排名"},
    {"id": "multiomics", "label": "多组学整合", "label_en": "Multi-omics Integration",
     "description": "DEG-DAM 相关性分析 + 联合 KEGG 富集 + 四象限图 + Spearman 相关网络"},
    {"id": "joint_enrichment", "label": "联合 KEGG 富集", "label_en": "Joint KEGG Enrichment",
     "description": "超几何检验检测 DEG/DAM 在 KEGG 通路中的联合富集"},
    {"id": "quadrant_plot", "label": "四象限图", "label_en": "Quadrant Plot",
     "description": "DEG-DAM log2FC 四象限/九象限分类散点图"},
    {"id": "correlation_network", "label": "Spearman 相关网络", "label_en": "Spearman Correlation Network",
     "description": "基于 Spearman 秩相关的基因-代谢物关联网络 + hub 检测 + 模块发现"},
    {"id": "wgcna", "label": "WGCNA 共表达网络", "label_en": "WGCNA Co-expression Network",
     "description": "加权基因共表达网络分析：软阈值 + TOM + 模块检测 + 模块-性状关联"},
    {"id": "o2pls", "label": "O2PLS 多变量整合", "label_en": "O2PLS Integration",
     "description": "双向正交偏最小二乘：转录组-代谢组联合分解 + 正交信号校正 + VIP 评分"},
    {"id": "tf_narrowing", "label": "TF 候选鉴定", "label_en": "TF Candidate Identification",
     "description": "基于相关性、基序、文献等多维证据筛选候选转录因子"},
    {"id": "hypothesis", "label": "假说合成", "label_en": "Hypothesis Synthesis",
     "description": "生成并排名竞争性机制假说 + 证伪检验 + 验证实验建议"},
]


def expand_analysis_selection(tools_to_run: list[str]) -> list[str]:
    """将用户可见步骤 ID 展开为图节点名称列表。

    如果 tools_to_run 为空或包含 "all"，返回全部 14 个节点。
    否则：始终运行的 4 个节点 + 选中步骤对应的节点。

    返回去重后的有序列表（始终运行节点在前，选中节点在后）。
    """
    if not tools_to_run:
        return sorted(ALL_PIPELINE_NODES)
    if "all" in tools_to_run:
        return sorted(ALL_PIPELINE_NODES)

    expanded: set[str] = set(ALWAYS_RUN_NODES)
    for step_id in tools_to_run:
        nodes = STEP_TO_NODES.get(step_id)
        if nodes:
            expanded |= nodes

    # 有序输出：始终运行节点在前，选中节点在后
    result = sorted(ALWAYS_RUN_NODES)
    for node in sorted(expanded - ALWAYS_RUN_NODES):
        result.append(node)
    return result


# ── 图路由函数 ──────────────────────────────────────────────


def gate_after_quality(state: RuntimeState) -> str:
    """数据质量门控。

    验证: 数据基本可信（样本数、缺失率、共同样本）
    通过: -> metabolite_validation
    失败: -> END
    """
    if state.data_quality_passed:
        return "metabolite_validation"
    return END


def route_after_classification(state: RuntimeState) -> str:
    """Mechanism-based branching (v4.0).

    Routes to different evidence acquisition branches based on mechanism_type:

      TRANSCRIPTIONAL  -> tf_narrowing (TF co-expression + motif)
      TRANSPORT        -> transport_evidence (transporter gene scan)
      STRESS           -> stress_evidence (stress marker activation)
      ENZYMATIC        -> tf_narrowing (enzyme genes are transcriptional targets)
      UNKNOWN / MIXED  -> tf_narrowing (default: most common path)

    All branches merge back at regulation_evidence.
    """
    mechanism_type = getattr(state, "mechanism_type", "") or ""

    if mechanism_type == "transport_redistribution":
        return "transport_evidence"
    if mechanism_type == "stress_induced_redistribution":
        return "stress_evidence"
    # transcriptional_regulation, enzymatic_regulation, unknown, mixed, etc.
    return "tf_narrowing"
