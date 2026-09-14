from __future__ import annotations

from datetime import datetime
from pydantic import BaseModel, Field

from phyto_reason.models.enums import WorkflowStage


class WorkflowDecision(BaseModel):
    stage: WorkflowStage = WorkflowStage.RESEARCH_PLANNING
    selected_workflow: str = ""
    sample_size: int = 0
    chosen_methods: list[str] = Field(default_factory=list)
    rejected_methods: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    reason: str = ""
    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())

    def to_summary(self) -> dict:
        stage_val = self.stage.value if isinstance(self.stage, WorkflowStage) else self.stage
        return {
            "stage": stage_val,
            "workflow": self.selected_workflow,
            "sample_size": self.sample_size,
            "chosen": self.chosen_methods,
            "rejected": self.rejected_methods,
            "n_warnings": len(self.warnings),
        }


