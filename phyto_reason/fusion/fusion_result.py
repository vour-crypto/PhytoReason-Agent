"""
fusion_result.py — 证据融合统一输出类型。
"""

from __future__ import annotations

from enum import Enum
from pydantic import BaseModel, Field

from phyto_reason.models.enums import ConfidenceLevel
from phyto_reason.models.evidence import EvidenceType


class FusionEvidence(BaseModel):
    """单个证据项在融合过程中的记录。"""
    source: str = ""
    evidence_type: EvidenceType | None = None
    raw_score: float = Field(default=0.0, ge=0.0, le=1.0)
    normalized_score: float = Field(default=0.0, ge=0.0, le=1.0)
    weight: float = Field(default=0.0, ge=0.0, le=1.0)
    contribution: float = Field(default=0.0, ge=0.0, le=1.0)
    explanation: str = ""

    def to_summary(self) -> str:
        return (
            f"[{self.source}] raw={self.raw_score:.2f} norm={self.normalized_score:.2f} "
            f"wt={self.weight:.2f} contrib={self.contribution:.2f}"
        )


class ContradictionLevel(str, Enum):
    NONE = "none"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class ContradictionFlag(BaseModel):
    """证据矛盾记录。"""
    level: ContradictionLevel = ContradictionLevel.NONE
    pattern: str = ""
    description: str = ""
    evidence_a: str = ""
    evidence_b: str = ""
    penalty: float = Field(default=0.0, ge=0.0, le=1.0)

    def to_summary(self) -> str:
        return (
            f"[{self.level.value}] {self.pattern}: "
            f"{self.evidence_a} vs {self.evidence_b} "
            f"(penalty={self.penalty:.2f})"
        )


class EvidenceHierarchy(BaseModel):
    """证据重要性层次。"""
    top_sources: list[str] = Field(default_factory=list)
    convergence_score: float = Field(default=0.0, ge=0.0, le=1.0)
    diversity_score: float = Field(default=0.0, ge=0.0, le=1.0)
    n_independent_sources: int = 0

    def to_summary(self) -> str:
        return (
            f"top={self.top_sources[:3]} convergence={self.convergence_score:.2f} "
            f"diversity={self.diversity_score:.2f} n_src={self.n_independent_sources}"
        )


class UnifiedFusionResult(BaseModel):
    """统一融合输出。"""
    target_gene_id: str = ""

    raw_fused_score: float = Field(default=0.0, ge=0.0, le=1.0)
    calibrated_score: float = Field(default=0.0, ge=0.0, le=1.0)
    confidence_level: ConfidenceLevel = ConfidenceLevel.NONE

    evidence_list: list[FusionEvidence] = Field(default_factory=list)
    contradictions: list[ContradictionFlag] = Field(default_factory=list)
    hierarchy: EvidenceHierarchy = Field(default_factory=EvidenceHierarchy)

    explanation: str = ""

    def to_summary(self) -> dict:
        level = self.confidence_level
        if isinstance(level, ConfidenceLevel):
            level = level.value
        return {
            "gene_id": self.target_gene_id,
            "raw_score": round(self.raw_fused_score, 3),
            "calibrated_score": round(self.calibrated_score, 3),
            "confidence_level": level,
            "n_evidences": len(self.evidence_list),
            "n_contradictions": len(self.contradictions),
            "convergence": round(self.hierarchy.convergence_score, 3),
            "diversity": self.hierarchy.n_independent_sources,
        }
