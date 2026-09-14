from __future__ import annotations

from pydantic import BaseModel, Field

from phyto_reason.models.enums import ConfidenceLevel
from phyto_reason.models.evidence import Evidence, EvidenceType


class CandidateTF(BaseModel):
    tf_id: str
    tf_family: str | None = None

    target_genes: list[str] = Field(default_factory=list)
    target_enzymes: list[str] = Field(default_factory=list)
    target_pathways: list[str] = Field(default_factory=list)

    evidence_list: list[Evidence] = Field(default_factory=list)

    confidence_level: ConfidenceLevel = ConfidenceLevel.NONE
    hypothesis: str | None = None

    def add_evidence(self, evidence: Evidence) -> None:
        self.evidence_list.append(evidence)

    def num_evidence_types(self) -> int:
        seen: set[str] = set()
        for ev in self.evidence_list:
            etype = ev.evidence_type
            seen.add(etype.value if isinstance(etype, EvidenceType) else etype)
        return len(seen)

    def compute_confidence_level(self) -> ConfidenceLevel:
        support_count = self.num_evidence_types()
        if support_count >= 3:
            self.confidence_level = ConfidenceLevel.GOLD
        elif support_count >= 2:
            self.confidence_level = ConfidenceLevel.SILVER
        elif support_count >= 1:
            self.confidence_level = ConfidenceLevel.WEAK
        else:
            self.confidence_level = ConfidenceLevel.NONE
        return self.confidence_level

    def to_summary(self) -> dict:
        level = self.confidence_level.value if isinstance(self.confidence_level, ConfidenceLevel) else self.confidence_level
        return {
            "tf_id": self.tf_id,
            "tf_family": self.tf_family,
            "n_target_genes": len(self.target_genes),
            "n_target_enzymes": len(self.target_enzymes),
            "n_target_pathways": len(self.target_pathways),
            "n_evidence_items": len(self.evidence_list),
            "confidence_level": level,
            "has_hypothesis": self.hypothesis is not None,
        }


