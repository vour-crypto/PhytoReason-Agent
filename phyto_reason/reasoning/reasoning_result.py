"""
reasoning_result.py — 统一的推理结果类型。

所有 reasoning 子模块的输出均使用 ReasoningResult。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class PriorEvidence(BaseModel):
    strength: str = "unknown"        # "strong" | "moderate" | "weak" | "unknown"
    score: float = 0.0               # 0.0-1.0
    description: str = ""
    pmids: list[str] = Field(default_factory=list)
    species_scope: str = "general"   # "general" | "solanaceae" | "fabaceae" | ...


class ReasoningResult(BaseModel):
    score: float = Field(default=0.0, ge=0.0, le=1.0)
    label: str = "not_assessed"      # "high" | "moderate" | "low" | "none" | "not_assessed"
    explanation: str = ""
    details: list[str] = Field(default_factory=list)
    evidence_list: list[PriorEvidence] = Field(default_factory=list)

    def update_label(self) -> None:
        if self.score >= 0.70:
            self.label = "high"
        elif self.score >= 0.40:
            self.label = "moderate"
        elif self.score >= 0.10:
            self.label = "low"
        else:
            self.label = "none"
