"""AgentState — unified state schema for all agent operations.
All tools read from and write to this state. No exceptions."""
from __future__ import annotations
from pydantic import BaseModel, Field
from typing import Any
from phyto_reason.models.mechanistic_hypothesis import MechanisticHypothesis
from phyto_reason.models.candidate_gene import CandidateGene


class AgentState(BaseModel):
    """Standardized state across all pipelines.

    Every pipeline (differential, regulation, exploration, QC)
    reads from and writes to this state.
    """
    # ── Data ───────────────────────────────────
    has_expression: bool = False
    has_metabolite: bool = False
    expression_matrix: Any = None     # dict[str, dict[str, float]]
    metabolite_matrix: Any = None
    annotation_index: Any = None      # AnnotationIndex

    # ── Samples ─────────────────────────────────
    sample_count: int = 0
    sample_ids: list[str] = Field(default_factory=list)
    group_a: list[str] = Field(default_factory=list)
    group_b: list[str] = Field(default_factory=list)
    species: str = ""
    target_metabolite: str = ""

    # ── QC / Data quality ──────────────────────
    data_quality_passed: bool = False
    data_quality_issues: list[str] = Field(default_factory=list)

    # ── Differential analysis ──────────────────
    deg_candidates: list[CandidateGene] = Field(default_factory=list)
    volcano_data: list[dict] = Field(default_factory=list)       # [{gene, lfc, p, q, np10}]
    pca_scores: dict = Field(default_factory=dict)               # {samples, pc1, pc2, var1, var2}

    # ── Pathway analysis ───────────────────────
    pathway_name: str = ""
    pathway_coherence_passed: bool = False
    pathway_genes: list[str] = Field(default_factory=list)
    pathway_enzymes_found: list[str] = Field(default_factory=list)

    # ── TF regulation (existing) ───────────────
    tf_candidates: dict[str, CandidateGene] = Field(default_factory=dict)
    fusion_results: Any = None

    # ── Hypothesis (existing) ──────────────────
    competing_hypotheses: list[MechanisticHypothesis] = Field(default_factory=list)
    falsification_results: list[str] = Field(default_factory=list)

    # ── Report ─────────────────────────────────
    hypothesis: str = ""
    ingestion_warnings: list[str] = Field(default_factory=list)

    # ── Execution history ──────────────────────
    executed_tools: list[str] = Field(default_factory=list)
    task: str = ""          # DIFFERENTIAL | REGULATION | QC | EXPLAIN | CHAT

    def to_summary(self) -> dict:
        return {
            "task": self.task,
            "samples": self.sample_count,
            "data_quality": self.data_quality_passed,
            "pathway": self.pathway_name if self.pathway_coherence_passed else None,
            "deg_count": len(self.deg_candidates),
            "tf_count": len(self.tf_candidates),
            "hypotheses": len(self.competing_hypotheses),
        }
