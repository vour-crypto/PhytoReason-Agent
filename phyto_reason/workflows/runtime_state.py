"""
runtime_state.py — LangGraph 运行时状态。
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, Field

from phyto_reason.models.planner_state import PlannerState
from phyto_reason.models.workflow_plan import WorkflowPlan
from phyto_reason.reasoning.confidence_evaluator import EvaluationResult
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence
from phyto_reason.fusion.fusion_result import UnifiedFusionResult
from phyto_reason.models.mechanistic_hypothesis import MechanisticHypothesis


class ExecutionStep(BaseModel):
    node: str = ""
    status: str = "pending"
    started_at: str = ""
    finished_at: str = ""
    summary: dict = Field(default_factory=dict)
    error: str | None = None
    retry_count: int = 0


class RuntimeState(BaseModel):
    session_id: str = Field(default_factory=lambda: datetime.now().strftime("PA-%Y%m%d-%H%M%S"))
    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())

    planner_state: PlannerState | None = None
    workflow_plan: WorkflowPlan | None = None
    tool_results: dict[str, ToolResult] = Field(default_factory=dict)
    fusion_results: dict[str, UnifiedFusionResult] = Field(default_factory=dict)
    evaluation: EvaluationResult | None = None
    hypothesis: str = ""

    execution_trace: list[ExecutionStep] = Field(default_factory=list)
    failures: list[dict] = Field(default_factory=list)

    current_node: str = "data_quality_gate"
    completed_nodes: list[str] = Field(default_factory=list)
    error: str | None = None
    finished: bool = False

    # ── Gate states ─────────────────────────────────────
    data_quality_passed: bool = False
    data_quality_issues: list[str] = Field(default_factory=list)
    data_quality_report: dict = Field(default_factory=dict)

    metabolite_validated: bool = False
    metabolite_validation_note: str = ""

    pathway_coherence_passed: bool = False
    pathway_coherence_score: float = 0.0
    pathway_coherence_evidence: str = ""
    pathway_name: str = ""
    pathway_genes: list[str] = Field(default_factory=list)

    # ── TF narrowing results ────────────────────────────
    tf_candidates: dict[str, CandidateGene] = Field(default_factory=dict)
    tf_pathway_gene_links: dict[str, list[str]] = Field(default_factory=dict)

    # ── Contradiction / artifact ────────────────────────
    excluded_tfs: list[str] = Field(default_factory=list)
    artifact_notes: dict[str, list[str]] = Field(default_factory=dict)

    # ── DEG/DAM analysis results ────────────────────────
    deg_report: dict = Field(default_factory=dict)
    # {n_genes_tested, n_fdr_sig, n_effect_ranked, top_genes: [{gene_id, effect_size, p_value, q_value}], warnings}
    dam_report: dict = Field(default_factory=dict)
    # {n_metabolites_tested, n_fdr_sig, n_effect_ranked, top_metabolites: [...]}
    qc_report: dict = Field(default_factory=dict)
    # {n_samples, n_genes, groups, figure_urls, figure_markdown, ...}

    multiomics_report: dict = Field(default_factory=dict)
    # {n_correlation_pairs_tested, n_significant_pairs, top_pairs, modules, warnings}

    # ── v4.7 Multi-omics extension reports ────────────────
    joint_enrichment_report: dict = Field(default_factory=dict)
    # {n_pathways_tested, n_significant_pathways, results: [...], significant: [...]}
    quadrant_plot_report: dict = Field(default_factory=dict)
    # {n_pairs, quadrant_counts, synchronicity_ratio, plot_base64, interpretation}
    correlation_network_report: dict = Field(default_factory=dict)
    # {n_nodes, n_edges, nodes: [...], edges: [...], hubs: [...], modules: [...]}

    wgcna_report: dict = Field(default_factory=dict)
    # {n_genes_analyzed, n_modules, soft_power, module_sizes, hub_genes, trait_correlations}

    o2pls_report: dict = Field(default_factory=dict)
    # {n_joint_components, x_joint_variance_explained, y_joint_variance_explained,
    #  top_x_variables, top_y_variables, joint_scores}

    # ── v3.0 Mechanism classification ──────────────────
    mechanism_type: str = ""
    mechanism_probabilities: dict = Field(default_factory=dict)
    mechanism_uncertainty: str = "high"

    # ── v3.0 Competing hypotheses ──────────────────────
    competing_hypotheses: list[MechanisticHypothesis] = Field(default_factory=list)

    # ── Ingestion (Phase 5) ────────────────────────────
    ingestion_warnings: list[str] = Field(default_factory=list)
    parsing_errors: list[str] = Field(default_factory=list)
    unmapped_entities: dict[str, list[str]] = Field(default_factory=dict)
    sample_alignment_report: dict = Field(default_factory=dict)
    annotation_index: Any = None

    # ── Graceful degradation ────────────────────────────
    degraded_tools: list[str] = Field(default_factory=list)
    degradation_notes: list[str] = Field(default_factory=list)

    # ── Progress streaming ──────────────────────────────
    progress_messages: list[str] = Field(default_factory=list)

    # ── Final hypothesis ────────────────────────────────
    mechanism_hypothesis: str = ""

    def add_trace(self, node: str, status: str, summary: dict | None = None,
                   error: str | None = None, retry_count: int = 0) -> None:
        step = ExecutionStep(
            node=node, status=status, summary=summary or {},
            error=error, retry_count=retry_count,
            started_at=datetime.now().isoformat(),
            finished_at=datetime.now().isoformat(),
        )
        self.execution_trace.append(step)
        if status == "completed":
            self.completed_nodes.append(node)
        elif status == "failed":
            self.failures.append({"node": node, "error": error})
            self.error = error
        # Also emit progress
        self.emit_progress(node, status)

    def emit_progress(self, node: str, status: str, message: str = "") -> None:
        """Emit a progress event for streaming feedback."""
        ts = datetime.now().strftime("%H:%M:%S")
        emoji = {"completed": "[OK]", "failed": "[FAIL]", "running": "[...]"}.get(status, "")
        msg = f"{ts} {emoji} {node}" + (f": {message}" if message else "")
        self.progress_messages.append(msg)

    def to_summary(self) -> dict:
        return {
            "session_id": self.session_id,
            "current_node": self.current_node,
            "n_completed": len(self.completed_nodes),
            "n_failures": len(self.failures),
            "data_quality_passed": self.data_quality_passed,
            "metabolite_validated": self.metabolite_validated,
            "pathway_coherence_passed": self.pathway_coherence_passed,
            "n_tf_candidates": len(self.tf_candidates),
            "n_excluded_tfs": len(self.excluded_tfs),
            "n_ingestion_warnings": len(self.ingestion_warnings),
            "has_hypothesis": len(self.mechanism_hypothesis) > 0,
            "finished": self.finished,
            "error": self.error is not None,
            "deg": {
                "n_tested": self.deg_report.get("n_genes_tested", 0),
                "n_fdr_sig": self.deg_report.get("n_fdr_significant", 0),
            } if self.deg_report else {},
            "dam": {
                "n_tested": self.dam_report.get("n_metabolites_tested", 0),
                "n_fdr_sig": self.dam_report.get("n_fdr_significant", 0),
            } if self.dam_report else {},
            "multiomics": {
                "n_pairs": self.multiomics_report.get("n_correlation_pairs_tested", 0),
                "n_sig": self.multiomics_report.get("n_significant_pairs", 0),
            } if self.multiomics_report else {},
            "joint_enrichment": {
                "n_pathways": self.joint_enrichment_report.get("n_pathways_tested", 0),
                "n_sig": self.joint_enrichment_report.get("n_significant_pathways", 0),
            } if self.joint_enrichment_report else {},
            "quadrant_plot": {
                "n_pairs": self.quadrant_plot_report.get("n_pairs", 0),
                "sync_ratio": self.quadrant_plot_report.get("synchronicity_ratio", 0),
            } if self.quadrant_plot_report else {},
            "correlation_network": {
                "n_edges": self.correlation_network_report.get("n_edges", 0),
                "n_hubs": len(self.correlation_network_report.get("hubs", [])),
                "n_modules": len(self.correlation_network_report.get("modules", [])),
            } if self.correlation_network_report else {},
            "wgcna": {
                "n_modules": self.wgcna_report.get("n_modules", 0),
                "soft_power": self.wgcna_report.get("soft_power", 0),
            } if self.wgcna_report else {},
            "o2pls": {
                "n_joint": self.o2pls_report.get("n_joint_components", 0),
                "x_var": self.o2pls_report.get("x_joint_variance_explained", 0),
            } if self.o2pls_report else {},
        }
