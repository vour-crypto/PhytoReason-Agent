from __future__ import annotations

from pydantic import BaseModel, Field

from phyto_reason.models.enums import BiologicalPlausibility, ConfidenceLevel
from phyto_reason.models.evidence import Evidence


class CandidateGene(BaseModel):
    gene_id: str
    gene_name: str | None = None
    species: str = ""

    module_id: str | None = None
    is_tf: bool = False
    tf_family: str | None = None
    is_synthetic: bool = False

    correlation_score: float = Field(default=0.0, ge=0.0, le=1.0)
    module_membership: float = Field(default=0.0, ge=0.0, le=1.0)
    motif_score: float = Field(default=0.0, ge=0.0, le=1.0)
    pathway_score: float = Field(default=0.0, ge=0.0, le=1.0)
    tissue_specificity_score: float = Field(default=0.0, ge=0.0, le=1.0)
    literature_score: float = Field(default=0.0, ge=0.0, le=1.0)
    ortholog_score: float = Field(default=0.0, ge=0.0, le=1.0)

    biological_plausibility: BiologicalPlausibility = BiologicalPlausibility.UNKNOWN
    confidence_score: float = Field(default=0.0, ge=0.0, le=1.0)
    confidence_level: ConfidenceLevel = ConfidenceLevel.NONE

    evidence_list: list[Evidence] = Field(default_factory=list)

    target_metabolites: list[str] = Field(default_factory=list)
    target_enzymes: list[str] = Field(default_factory=list)
    suggested_validation: list[str] = Field(default_factory=list)

    def add_evidence(self, evidence: Evidence) -> None:
        self.evidence_list.append(evidence)

    def compute_confidence_score(
        self,
        weights: dict[str, float] | None = None,
    ) -> float:
        if weights is None:
            weights = {
                "correlation": 0.15,
                "module_membership": 0.20,
                "motif": 0.10,
                "pathway": 0.15,
                "tissue": 0.05,
                "literature": 0.15,
                "ortholog": 0.05,
                "plausibility": 0.15,
            }

        score = 0.0
        score += self.correlation_score * weights.get("correlation", 0.0)
        score += self.module_membership * weights.get("module_membership", 0.0)
        score += self.motif_score * weights.get("motif", 0.0)
        score += self.pathway_score * weights.get("pathway", 0.0)
        score += self.tissue_specificity_score * weights.get("tissue", 0.0)
        score += self.literature_score * weights.get("literature", 0.0)
        score += self.ortholog_score * weights.get("ortholog", 0.0)

        plausibility_map = {
            BiologicalPlausibility.PLAUSIBLE: 1.0,
            BiologicalPlausibility.WEAK: 0.5,
            BiologicalPlausibility.UNLIKELY: 0.0,
            BiologicalPlausibility.UNKNOWN: 0.0,
        }
        score += plausibility_map.get(self.biological_plausibility, 0.0) * weights.get("plausibility", 0.0)

        self.confidence_score = min(score, 1.0)
        self._assign_level()
        return self.confidence_score

    def _assign_level(self) -> None:
        if self.confidence_score >= 0.70:
            self.confidence_level = ConfidenceLevel.GOLD
        elif self.confidence_score >= 0.40:
            self.confidence_level = ConfidenceLevel.SILVER
        else:
            self.confidence_level = ConfidenceLevel.WEAK

    def to_summary(self) -> dict:
        level = self.confidence_level.value if isinstance(self.confidence_level, ConfidenceLevel) else self.confidence_level
        bio = self.biological_plausibility.value if isinstance(self.biological_plausibility, BiologicalPlausibility) else self.biological_plausibility
        return {
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "tf_family": self.tf_family,
            "confidence_score": round(self.confidence_score, 3),
            "confidence_level": level,
            "biological_plausibility": bio,
            "n_evidences": len(self.evidence_list),
            "n_target_metabolites": len(self.target_metabolites),
        }


