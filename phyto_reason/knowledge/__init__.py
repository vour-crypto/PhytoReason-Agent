"""knowledge — Zero-data knowledge layer for PhytoReason-Agent.

Enables Layer 0 (no uploaded data) through public database queries
and built-in biological prior knowledge.

Modules:
    pubmed_retriever:  Structured PubMed search with smart query building
    kegg_browser:      KEGG pathway/compound browsing
    tf_knowledge_base: TF family ↔ metabolite regulatory relationships
    species_knowledge: Cross-species ortholog inference
    species_registry:  Species profile registry + degrade chain (v5.0)
    compound_profiles: Metabolite-class MS/MS diagnostic rules (v5.0)
"""

from phyto_reason.knowledge.pubmed_retriever import (
    PubMedRetriever, LiteratureResult,
)
from phyto_reason.knowledge.kegg_browser import (
    KEGGBrowser, PathwayInfo,
)
from phyto_reason.knowledge.tf_knowledge_base import (
    TFKnowledgeBase, TFMetaboliteRelation,
)
from phyto_reason.knowledge.species_knowledge import (
    SpeciesKnowledge, OrthologInference,
)
from phyto_reason.knowledge.ortholog_finder import (
    OrthologFinder, OrthologHit, OrthologResult, find_orthologs,
)
from phyto_reason.knowledge.species_registry import (
    SpeciesRegistry, SpeciesProfile, MarkerMetabolite, get_species_registry,
)
from phyto_reason.knowledge.compound_profiles import (
    load_compound_profile,
    list_compound_classes,
    find_diagnostic_rules,
    get_class_info,
)

__all__ = [
    "PubMedRetriever", "LiteratureResult",
    "KEGGBrowser", "PathwayInfo",
    "TFKnowledgeBase", "TFMetaboliteRelation",
    "SpeciesKnowledge", "OrthologInference",
    "OrthologFinder", "OrthologHit", "OrthologResult", "find_orthologs",
    "SpeciesRegistry", "SpeciesProfile", "MarkerMetabolite", "get_species_registry",
    "load_compound_profile", "list_compound_classes",
    "find_diagnostic_rules", "get_class_info",
]
