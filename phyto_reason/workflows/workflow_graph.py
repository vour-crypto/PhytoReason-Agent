"""
workflow_graph.py -- 科研推理状态机 (v4.7).

图结构（20 节点，分支 DAG）:

                            deg_analysis ─── dam_analysis
                                  │
                          multiomics_integration
                                  │
                     joint_enrichment ── quadrant_plot
                                  │
                          correlation_network
                                  │
                                wgcna
                                  │
                                o2pls
                                  │
                            data_quality_gate
                                  │ PASS
                                  ▼
                         metabolite_validation
                                  │
                                  ▼
                          pathway_coherence
                                  │
                                  ▼
                        mechanism_classifier
                        /        |         \
                TRANSCRIPTIONAL  TRANSPORT  STRESS
                      │              │         │
                tf_narrowing  transport_ev  stress_ev
                      │              │         │
                      └──────────────┼─────────┘
                                     │
                            regulation_evidence  <- 融合点
                                     │
                            contradiction_check
                                     │
                               falsification
                                     │
                          hypothesis_competition
                                     │
                           hypothesis_synthesis ──-> END

v4.7 新增: joint_enrichment, quadrant_plot, correlation_network, wgcna, o2pls (5 个多组学模块)

路由集中度:
  - 所有条件路由在 routing_logic.py 中
  - graph 只做 orchestration, 不包含科学逻辑
  - node 不包含 next_node 逻辑
"""

from __future__ import annotations

from langgraph.graph import StateGraph, END

from phyto_reason.workflows.runtime_state import RuntimeState
from phyto_reason.workflows.execution_nodes import (
    data_qc_node,
    deg_analysis_node,
    dam_analysis_node,
    multiomics_integration_node,
    data_quality_gate_node,
    metabolite_validation_node,
    pathway_coherence_node,
    mechanism_classifier_node,
    tf_narrowing_node,
    regulation_evidence_node,
    contradiction_check_node,
    hypothesis_competition_node,
    falsification_node,
    hypothesis_synthesis_node,
)
from phyto_reason.workflows.nodes.transport_evidence import transport_evidence_node
from phyto_reason.workflows.nodes.stress_evidence import stress_evidence_node
from phyto_reason.workflows.nodes.joint_enrichment import joint_enrichment_node
from phyto_reason.workflows.nodes.quadrant_plot import quadrant_plot_node
from phyto_reason.workflows.nodes.correlation_network import correlation_network_node
from phyto_reason.workflows.nodes.wgcna import wgcna_node
from phyto_reason.workflows.nodes.o2pls import o2pls_node
from phyto_reason.workflows.routing_logic import (
    gate_after_quality,
    route_after_classification,
)


def _node_with_trace(fn):
    """节点包装：LangGraph 给每个节点发状态副本，节点内原地写入的
    execution_trace / completed_nodes / failures 会随副本丢弃——
    这里在节点返回 dict 中带上副本上的累计 trace，随 channel 合并传播。

    注: RuntimeState 字段为 replace 语义，每次返回完整累计 trace 即可。
    """

    def wrapped(state):
        result = fn(state) or {}
        if getattr(state, "execution_trace", None):
            result["execution_trace"] = state.execution_trace
        if getattr(state, "completed_nodes", None):
            result["completed_nodes"] = state.completed_nodes
        if getattr(state, "failures", None):
            result["failures"] = state.failures
        return result

    return wrapped


def build_workflow() -> StateGraph:
    workflow = StateGraph(RuntimeState)

    # ── 注册所有节点（统一包 _node_with_trace）───────────
    workflow.add_node("data_qc", _node_with_trace(data_qc_node))
    workflow.add_node("deg_analysis", _node_with_trace(deg_analysis_node))
    workflow.add_node("dam_analysis", _node_with_trace(dam_analysis_node))
    workflow.add_node("multiomics_integration", _node_with_trace(multiomics_integration_node))
    workflow.add_node("joint_enrichment", _node_with_trace(joint_enrichment_node))
    workflow.add_node("quadrant_plot", _node_with_trace(quadrant_plot_node))
    workflow.add_node("correlation_network", _node_with_trace(correlation_network_node))
    workflow.add_node("wgcna", _node_with_trace(wgcna_node))
    workflow.add_node("o2pls", _node_with_trace(o2pls_node))
    workflow.add_node("data_quality_gate", _node_with_trace(data_quality_gate_node))
    workflow.add_node("metabolite_validation", _node_with_trace(metabolite_validation_node))
    workflow.add_node("pathway_coherence", _node_with_trace(pathway_coherence_node))
    workflow.add_node("mechanism_classifier", _node_with_trace(mechanism_classifier_node))
    workflow.add_node("tf_narrowing", _node_with_trace(tf_narrowing_node))
    workflow.add_node("transport_evidence", _node_with_trace(transport_evidence_node))
    workflow.add_node("stress_evidence", _node_with_trace(stress_evidence_node))
    workflow.add_node("regulation_evidence", _node_with_trace(regulation_evidence_node))
    workflow.add_node("contradiction_check", _node_with_trace(contradiction_check_node))
    workflow.add_node("hypothesis_competition", _node_with_trace(hypothesis_competition_node))
    workflow.add_node("falsification", _node_with_trace(falsification_node))
    workflow.add_node("hypothesis_synthesis", _node_with_trace(hypothesis_synthesis_node))

    workflow.set_entry_point("data_qc")

    # ── QC -> DEG -> DAM -> Multi-omics -> ... -> O2PLS -> quality_gate ──
    workflow.add_edge("data_qc", "deg_analysis")
    workflow.add_edge("deg_analysis", "dam_analysis")
    workflow.add_edge("dam_analysis", "multiomics_integration")
    workflow.add_edge("multiomics_integration", "joint_enrichment")
    workflow.add_edge("joint_enrichment", "quadrant_plot")
    workflow.add_edge("quadrant_plot", "correlation_network")
    workflow.add_edge("correlation_network", "wgcna")
    workflow.add_edge("wgcna", "o2pls")
    workflow.add_edge("o2pls", "data_quality_gate")

    workflow.add_edge("deg_analysis", "dam_analysis")

    # ── 条件路由 ─────────────────────────────────────
    workflow.add_conditional_edges(
        "data_quality_gate",
        gate_after_quality,
        {"metabolite_validation": "metabolite_validation", END: END},
    )

    # mechanism_classifier -> 3 条分支
    workflow.add_conditional_edges(
        "mechanism_classifier",
        route_after_classification,
        {
            "tf_narrowing": "tf_narrowing",
            "transport_evidence": "transport_evidence",
            "stress_evidence": "stress_evidence",
        },
    )

    # ── 无条件边 ─────────────────────────────────────
    workflow.add_edge("metabolite_validation", "pathway_coherence")
    workflow.add_edge("pathway_coherence", "mechanism_classifier")

    # All 3 branches merge at regulation_evidence
    workflow.add_edge("tf_narrowing", "regulation_evidence")
    workflow.add_edge("transport_evidence", "regulation_evidence")
    workflow.add_edge("stress_evidence", "regulation_evidence")

    workflow.add_edge("regulation_evidence", "contradiction_check")
    workflow.add_edge("contradiction_check", "falsification")
    workflow.add_edge("falsification", "hypothesis_competition")
    workflow.add_edge("hypothesis_competition", "hypothesis_synthesis")
    workflow.add_edge("hypothesis_synthesis", END)

    return workflow


def compile_workflow():
    graph = build_workflow()
    return graph.compile()
