"""
planner_state.py — PI-Agent 规划状态模型。
"""

from __future__ import annotations

from datetime import datetime
from pydantic import BaseModel, Field


class PlannerDecision(BaseModel):
    stage: str = ""
    decision: str = ""
    rationale: str = ""
    warnings: list[str] = Field(default_factory=list)
    alternatives: list[str] = Field(default_factory=list)
    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())


class PlannerState(BaseModel):
    research_question: str = ""
    species: str = ""
    target_metabolite: str = ""
    target_pathway: str = ""

    sample_count: int = 0
    has_expression: bool = False
    has_metabolite: bool = False
    has_annotation: bool = False
    has_promoter: bool = False
    has_itak: bool = False
    missing_data_ratio: float = 0.0

    decisions: list[PlannerDecision] = Field(default_factory=list)
    selected_tools: list[str] = Field(default_factory=list)
    skipped_analyses: list[str] = Field(default_factory=list)

    confidence_assessment: str = "not_assessed"
    overall_feasibility: str = "unknown"
    warnings: list[str] = Field(default_factory=list)

    final_hypothesis: str = ""
    top_candidates: list[str] = Field(default_factory=list)

    # Real data (loaded from uploaded files or file paths)
    data_path: str = ""
    expression_matrix: dict | None = None
    metabolite_matrix: dict | None = None
    promoter_sequences: dict | None = None
    # Sample-level design metadata: {sample_id: {condition, tissue, ...}}
    sample_metadata: dict | None = None

    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())
    updated_at: str = Field(default_factory=lambda: datetime.now().isoformat())

    def add_decision(self, stage: str, decision: str, rationale: str,
                     warnings: list[str] | None = None,
                     alternatives: list[str] | None = None) -> None:
        self.decisions.append(PlannerDecision(
            stage=stage, decision=decision, rationale=rationale,
            warnings=warnings or [], alternatives=alternatives or [],
        ))
        self._touch()

    def add_warning(self, warning: str) -> None:
        if warning not in self.warnings:
            self.warnings.append(warning)
            self._touch()

    def _touch(self) -> None:
        self.updated_at = datetime.now().isoformat()

    def to_summary(self) -> dict:
        return {
            "question": self.research_question,
            "species": self.species,
            "target": self.target_metabolite,
            "sample_count": self.sample_count,
            "n_decisions": len(self.decisions),
            "selected_tools": self.selected_tools,
            "skipped": self.skipped_analyses,
            "feasibility": self.overall_feasibility,
            "warnings": len(self.warnings),
            "n_top_candidates": len(self.top_candidates),
        }
