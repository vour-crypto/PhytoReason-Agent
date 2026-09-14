"""
GlobalState — 轻量级全局状态。

不存储大型矩阵（表达矩阵、代谢物矩阵等）。
仅保存:
  - workflow decisions
  - candidate objects
  - evidence
  - tool output summaries
  - lightweight metadata
"""

from __future__ import annotations

from datetime import datetime
from pydantic import BaseModel, Field

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.candidate_tf import CandidateTF
from phyto_reason.models.enums import WorkflowStage
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.workflow_decision import WorkflowDecision
from phyto_reason.models.tool_result import ToolResult


class DataQualityReport(BaseModel):
    has_expression_data: bool = False
    has_metabolite_data: bool = False
    has_annotation: bool = False
    has_itak: bool = False
    sample_count: int = 0
    common_samples: int = 0
    species: str = ""
    target_metabolite: str = ""
    warnings: list[str] = Field(default_factory=list)
    passed: bool = False


class GlobalState(BaseModel):
    session_id: str = Field(default_factory=lambda: datetime.now().strftime("PA-%Y%m%d-%H%M%S"))
    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())
    updated_at: str = Field(default_factory=lambda: datetime.now().isoformat())

    research_question: str = ""
    species: str = ""
    target_metabolite: str = ""
    data_quality: DataQualityReport = Field(default_factory=DataQualityReport)

    workflow_decisions: list[WorkflowDecision] = Field(default_factory=list)

    candidates: dict[str, CandidateTF] = Field(default_factory=dict)
    candidate_genes: dict[str, CandidateGene] = Field(default_factory=dict)

    evidence_pool: dict[str, list[Evidence]] = Field(default_factory=dict)

    tool_outputs: dict[str, ToolResult] = Field(default_factory=dict)

    current_stage: WorkflowStage = WorkflowStage.RESEARCH_PLANNING
    completed_stages: list[WorkflowStage] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)

    def add_candidate_tf(self, tf: CandidateTF) -> None:
        self.candidates[tf.tf_id] = tf
        self._touch()

    def add_candidate_gene(self, gene: CandidateGene) -> None:
        self.candidate_genes[gene.gene_id] = gene
        self._touch()

    def add_workflow_decision(self, decision: WorkflowDecision) -> None:
        self.workflow_decisions.append(decision)
        self._touch()

    def add_evidence(self, evidence: Evidence) -> None:
        etype = evidence.evidence_type
        if isinstance(etype, EvidenceType):
            etype = etype.value
        if etype not in self.evidence_pool:
            self.evidence_pool[etype] = []
        self.evidence_pool[etype].append(evidence)
        self._touch()

    def add_tool_output(self, tool_name: str, result: ToolResult) -> None:
        self.tool_outputs[tool_name] = result
        self._touch()

    def complete_stage(self, stage: WorkflowStage) -> None:
        if stage not in self.completed_stages:
            self.completed_stages.append(stage)
        self.current_stage = stage
        self._touch()

    def add_error(self, error: str) -> None:
        self.errors.append(error)
        self._touch()

    def to_summary(self) -> dict:
        stage = self.current_stage.value if isinstance(self.current_stage, WorkflowStage) else self.current_stage
        return {
            "session_id": self.session_id,
            "species": self.species,
            "target_metabolite": self.target_metabolite,
            "current_stage": stage,
            "n_completed_stages": len(self.completed_stages),
            "n_candidate_tfs": len(self.candidates),
            "n_candidate_genes": len(self.candidate_genes),
            "n_evidence_types": len(self.evidence_pool),
            "n_tool_outputs": len(self.tool_outputs),
            "data_quality_passed": self.data_quality.passed,
            "n_errors": len(self.errors),
        }

    def _touch(self) -> None:
        self.updated_at = datetime.now().isoformat()


