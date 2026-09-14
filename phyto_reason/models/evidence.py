from __future__ import annotations

from datetime import datetime
from enum import Enum
from pydantic import BaseModel, Field


class EvidenceType(str, Enum):
    WGCNA_MODULE_MEMBERSHIP = "wgcna_module_membership"
    CORRELATION = "correlation"
    MOTIF_BINDING = "motif_binding"
    PATHWAY_CONSISTENCY = "pathway_consistency"
    TF_FAMILY_PRIOR = "tf_family_prior"
    LITERATURE = "literature"
    TISSUE_SPECIFICITY = "tissue_specificity"
    ORTHOLOG = "ortholog"
    ENRICHMENT = "enrichment"
    ANNOTATION = "annotation"
    EXPRESSION = "expression"
    DEG = "differential_expression"


class EvidenceRelationType(str, Enum):
    SUPPORTS = "supports"
    CONTRADICTS = "contradicts"
    CONDITIONAL = "conditionally_supports"
    INDIRECT = "indirectly_supports"
    SPECIES_LIMITED = "species_limited_evidence"
    TEMPORAL = "temporal_evidence"
    CONTEXT_DEPENDENT = "context_dependent"
    INFERRED = "inferred_by_homology"
    INSUFFICIENT = "insufficient_to_conclude"


class Evidence(BaseModel):
    evidence_type: EvidenceType
    source: str = ""
    score: float = Field(default=0.0, ge=0.0, le=1.0)
    confidence: float = Field(default=1.0, ge=0.0, le=1.0)
    description: str = ""
    metadata: dict = Field(default_factory=dict)
    created_at: str = Field(default_factory=lambda: datetime.now().isoformat())
    relation_type: EvidenceRelationType = EvidenceRelationType.SUPPORTS
    species_scope: str = ""
    tissue_scope: str = ""
    condition_scope: str = ""

    def to_summary(self) -> str:
        return (
            f"[{self.evidence_type.value}] {self.source} "
            f"score={self.score:.2f} conf={self.confidence:.2f} | {self.description[:100]}"
        )
