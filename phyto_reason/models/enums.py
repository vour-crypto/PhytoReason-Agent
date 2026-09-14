from __future__ import annotations

from enum import Enum


class ConfidenceLevel(str, Enum):
    GOLD = "Gold"
    SILVER = "Silver"
    WEAK = "Weak"
    NONE = "None"


class BiologicalPlausibility(str, Enum):
    PLAUSIBLE = "plausible"
    WEAK = "weak"
    UNLIKELY = "unlikely"
    UNKNOWN = "unknown"


class WorkflowStage(str, Enum):
    RESEARCH_PLANNING = "research_planning"
    DATA_QUALITY_CHECK = "data_quality_check"
    WORKFLOW_SELECTION = "workflow_selection"
    OMICS_ANALYSIS = "omics_analysis"
    BIOLOGICAL_REASONING = "biological_reasoning"
    EVIDENCE_FUSION = "evidence_fusion"
    HYPOTHESIS_GENERATION = "hypothesis_generation"
    REPORT_GENERATION = "report_generation"
    COMPLETED = "completed"
    FAILED = "failed"
