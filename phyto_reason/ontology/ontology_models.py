"""
ontology_models.py — 本体数据模型。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class PathwayNode(BaseModel):
    node_id: str = ""
    name: str = ""
    parent: str | None = None
    children: list[str] = Field(default_factory=list)
    kegg_ids: list[str] = Field(default_factory=list)
    enzymes: list[str] = Field(default_factory=list)
    bottleneck_enzymes: list[str] = Field(default_factory=list)


class MetaboliteNode(BaseModel):
    node_id: str = ""
    name: str = ""
    parent: str | None = None
    children: list[str] = Field(default_factory=list)
    aliases: list[str] = Field(default_factory=list)
    kegg_id: str | None = None
    pathway: str | None = None
    pathway_ids: list[str] = Field(default_factory=list)
    is_leaf: bool = False


class TFPriorEntry(BaseModel):
    tf_family: str = ""
    target_pathway: str = ""
    strength: str = "weak"
    score: float = Field(default=0.0, ge=0.0, le=1.0)
    pmids: list[str] = Field(default_factory=list)
    description: str = ""
    inheritable: bool = True


class OntologyResult(BaseModel):
    original_term: str = ""
    resolved_node: str = ""
    node_type: str = ""  # "metabolite" | "pathway" | "tf_family"
    aliases: list[str] = Field(default_factory=list)
    parent_nodes: list[str] = Field(default_factory=list)
    child_nodes: list[str] = Field(default_factory=list)
    pathway_ids: list[str] = Field(default_factory=list)
    kegg_id: str | None = None
    inheritance_chain: list[str] = Field(default_factory=list)
    tf_prior_scores: dict[str, float] = Field(default_factory=dict)
    error: str | None = None

    def to_summary(self) -> str:
        chain = " → ".join(self.inheritance_chain) if self.inheritance_chain else self.resolved_node
        return f"{self.original_term} → {chain} (parents={self.parent_nodes})"
