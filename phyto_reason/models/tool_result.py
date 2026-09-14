"""ToolResult — standard output for all tools."""
from __future__ import annotations
from datetime import datetime
from pydantic import BaseModel, Field
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.candidate_tf import CandidateTF
from phyto_reason.models.evidence import Evidence


class ToolResult(BaseModel):
    """Standardized output for ALL tools.

    Every tool must return this. No exceptions.
    """
    status: str = "success"                  # success | partial | failed | skipped
    candidates: list[CandidateGene | CandidateTF] = Field(default_factory=list)
    evidence_list: list[Evidence] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    plots: list[dict] = Field(default_factory=list)      # {"type":"png","data":"base64..."}
    metadata: dict = Field(default_factory=dict)
    next_suggestions: list[str] = Field(default_factory=list)
    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())

    def merge(self, other: ToolResult) -> ToolResult:
        self.candidates.extend(other.candidates)
        self.evidence_list.extend(other.evidence_list)
        self.warnings.extend(other.warnings)
        self.plots.extend(other.plots)
        self.metadata.update(other.metadata)
        self.next_suggestions.extend(other.next_suggestions)
        if other.status == "failed":
            self.status = "partial"
        return self

    def to_summary(self) -> dict:
        return {
            "status": self.status,
            "n_candidates": len(self.candidates),
            "n_evidences": len(self.evidence_list),
            "n_warnings": len(self.warnings),
            "n_plots": len(self.plots),
            "n_suggestions": len(self.next_suggestions),
            "metadata_keys": list(self.metadata.keys()),
        }
