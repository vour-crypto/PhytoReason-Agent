"""
workflow_plan.py — WorkflowPlan model.

Moved from brain/workflow_selector.py (deprecated).
"""

from __future__ import annotations

from pydantic import BaseModel


class WorkflowPlan(BaseModel):
    """Workflow execution plan. Used by RuntimeState and workflow graph."""
    selected_workflow: str = ""
    tools_to_run: list[str] = []           # expanded graph-node names (deg_analysis, dam_analysis, ...)
    tools_to_skip: list[str] = []
    original_selection: list[str] = []      # user-facing step IDs (DEG, DAM, multiomics, ...)
    reasoning: str = ""
    estimated_duration: str = ""
