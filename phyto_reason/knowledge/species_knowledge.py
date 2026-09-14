"""
species_knowledge.py — Cross-species knowledge and ortholog inference.

Enables Layer 1 (cross-species mode): given a model species' known
TF → metabolite relationship, infer whether a target species may have
similar regulation.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class OrthologInference:
    """Cross-species regulatory inference result."""
    source_species: str
    target_species: str
    tf_family: str
    metabolite: str
    confidence: str = "low"       # "high" | "medium" | "low" | "speculative"
    reasoning: str = ""
    source_pmids: list[str] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)


@dataclass
class SpeciesKnowledgeResult:
    query_species: str
    query_metabolite: str = ""
    known_pathways: list[str] = field(default_factory=list)
    known_metabolites: list[str] = field(default_factory=list)
    ortholog_inferences: list[OrthologInference] = field(default_factory=list)
    note: str = ""


# ── Species-Pathway-Metabolite Knowledge ──────────────────────
# Maps species → known secondary metabolites → associated pathways
SPECIES_KNOWLEDGE: dict[str, dict] = {
    "coptis": {
        "common_name": "Coptis (Huanglian)",
        "metabolites": ["berberine", "coptisine", "palmatine", "jateorhizine"],
        "pathways": ["isoquinoline alkaloid biosynthesis"],
        "tf_families": ["bHLH", "WRKY", "MYB", "ERF"],
        "model_species": "Coptis japonica, Coptis chinensis",
        "refs": ["25296257"],
    },
    "arabidopsis": {
        "common_name": "Arabidopsis thaliana",
        "metabolites": ["anthocyanin", "flavonol", "glucosinolate", "lignin"],
        "pathways": ["flavonoid biosynthesis", "glucosinolate biosynthesis"],
        "tf_families": ["MYB", "bHLH", "WD40", "WRKY", "NAC", "ERF", "bZIP"],
        "model_species": "Arabidopsis thaliana",
        "refs": ["25228336", "26847441"],
    },
    "tobacco": {
        "common_name": "Tobacco (Nicotiana)",
        "metabolites": ["nicotine", "nornicotine", "anatabine", "anabasine"],
        "pathways": ["pyridine alkaloid biosynthesis"],
        "tf_families": ["ERF", "bHLH", "MYC2", "WRKY"],
        "model_species": "Nicotiana tabacum, Nicotiana benthamiana",
        "refs": ["17419843", "19033552"],
    },
    "catharanthus": {
        "common_name": "Madagascar periwinkle",
        "metabolites": ["vincristine", "vinblastine", "catharanthine", "vindoline"],
        "pathways": ["terpenoid indole alkaloid biosynthesis"],
        "tf_families": ["ERF", "bHLH", "MYB", "WRKY"],
        "model_species": "Catharanthus roseus",
        "refs": ["17419843", "25296257"],
    },
    "artemisia": {
        "common_name": "Artemisia annua (Sweet wormwood)",
        "metabolites": ["artemisinin", "arteannuin B"],
        "pathways": ["sesquiterpenoid biosynthesis"],
        "tf_families": ["ERF", "WRKY", "MYB", "bZIP"],
        "model_species": "Artemisia annua",
        "refs": ["28417071", "27207470"],
    },
    "rice": {
        "common_name": "Rice (Oryza sativa)",
        "metabolites": ["flavonoid", "anthocyanin", "diterpenoid phytoalexin"],
        "pathways": ["flavonoid biosynthesis", "diterpenoid biosynthesis"],
        "tf_families": ["MYB", "bHLH", "WRKY", "bZIP", "NAC"],
        "model_species": "Oryza sativa",
        "refs": ["31530398"],
    },
    "tomato": {
        "common_name": "Tomato (Solanum lycopersicum)",
        "metabolites": ["anthocyanin", "flavonoid", "tomatine"],
        "pathways": ["flavonoid biosynthesis", "steroidal glycoalkaloid biosynthesis"],
        "tf_families": ["MYB", "bHLH", "ERF", "WRKY"],
        "model_species": "Solanum lycopersicum",
        "refs": ["26847441"],
    },
    "grape": {
        "common_name": "Grape (Vitis vinifera)",
        "metabolites": ["anthocyanin", "resveratrol", "flavonoid", "tannin"],
        "pathways": ["flavonoid biosynthesis", "stilbenoid biosynthesis"],
        "tf_families": ["MYB", "bHLH", "WD40", "WRKY"],
        "model_species": "Vitis vinifera",
        "refs": ["26847441"],
    },
    "maize": {
        "common_name": "Maize (Zea mays)",
        "metabolites": ["anthocyanin", "flavonoid", "benzoxazinoid"],
        "pathways": ["flavonoid biosynthesis", "benzoxazinoid biosynthesis"],
        "tf_families": ["MYB", "bHLH", "bZIP"],
        "model_species": "Zea mays",
        "refs": ["31530398"],
    },
}

# ── Conservation rules: when TF family X regulates metabolite Y in species A,
#   how likely is the same relationship in species B?
CONSERVATION_RULES = [
    # (tf_family, metabolite_class, conservation_level, note)
    ("MYB", "anthocyanin", "high", "MBW complex is deeply conserved in angiosperms"),
    ("MYB", "flavonoid", "high", "R2R3-MYB flavonoid regulation conserved in all land plants"),
    ("bHLH", "anthocyanin", "high", "MBW complex component, deeply conserved"),
    ("WD40", "anthocyanin", "high", "MBW complex component, deeply conserved"),
    ("ERF", "alkaloid", "medium", "JA-responsive ERFs regulate alkaloids but species-specific"),
    ("WRKY", "alkaloid", "medium", "WRKY-alkaloid link present in multiple families"),
    ("bHLH", "alkaloid", "low", "Species-specific bHLH-alkaloid interactions"),
    ("NAC", "lignin", "high", "NAC secondary wall regulation deeply conserved in vascular plants"),
    ("MYB", "lignin", "high", "MYB58/63-type lignin regulation conserved in angiosperms"),
    ("ERF", "terpenoid", "low", "Species-specific ERF-terpenoid interactions"),
]


class SpeciesKnowledge:
    """Cross-species knowledge and ortholog inference engine."""

    def __init__(self) -> None:
        self._species = SPECIES_KNOWLEDGE
        self._conservation = CONSERVATION_RULES

    def lookup(self, species: str) -> dict | None:
        """Look up a species in the knowledge base."""
        species_lower = species.lower().strip()
        for key, info in self._species.items():
            if key in species_lower or species_lower in key:
                return info
        # Try common name match
        for key, info in self._species.items():
            common = info.get("common_name", "").lower()
            if species_lower in common or key in species_lower:
                return info
        return None

    def infer_regulation(
        self,
        target_species: str,
        metabolite: str = "",
        tf_family: str = "",
        source_species: str = "",
    ) -> SpeciesKnowledgeResult:
        """Infer whether a TF family likely regulates a metabolite in a target species.

        Uses conservation rules and known species-specific knowledge.
        """
        species_info = self.lookup(target_species)
        known_mets = species_info.get("metabolites", []) if species_info else []
        known_paths = species_info.get("pathways", []) if species_info else []
        known_tfs = species_info.get("tf_families", []) if species_info else []

        inferences: list[OrthologInference] = []

        # 1. If TF family is known in this species, it's a direct match
        if tf_family and species_info:
            tf_lower = tf_family.lower()
            if any(tf_lower in kt.lower() for kt in known_tfs):
                inferences.append(OrthologInference(
                    source_species=source_species or target_species,
                    target_species=target_species,
                    tf_family=tf_family,
                    metabolite=metabolite,
                    confidence="high",
                    reasoning=f"{tf_family} TFs are known in {target_species}",
                    caveats=["Specific target genes need experimental validation"],
                ))

        # 2. Check conservation rules
        if tf_family and metabolite:
            tf_lower = tf_family.lower()
            met_lower = metabolite.lower()
            for rule_tf, rule_meta, level, note in self._conservation:
                if tf_lower in rule_tf.lower() and met_lower in rule_meta.lower():
                    inferences.append(OrthologInference(
                        source_species=source_species or "model species",
                        target_species=target_species,
                        tf_family=tf_family,
                        metabolite=metabolite,
                        confidence=level,
                        reasoning=note,
                        caveats=[
                            "Cross-species inference should be validated experimentally",
                            "Promoter motif conservation should be checked",
                        ],
                    ))

        # 3. Build note
        note = ""
        if species_info:
            note = (
                f"{species_info['common_name']}: known metabolites include "
                f"{', '.join(known_mets[:5])}. "
                f"TF families known: {', '.join(known_tfs)}."
            )
        else:
            note = (
                f"No species-specific knowledge for '{target_species}'. "
                f"Using general conservation rules only."
            )

        if not inferences:
            note += (
                f" No direct regulatory inference available for {tf_family} → {metabolite} "
                f"in {target_species}."
            )

        return SpeciesKnowledgeResult(
            query_species=target_species,
            query_metabolite=metabolite,
            known_pathways=known_paths,
            known_metabolites=known_mets,
            ortholog_inferences=inferences,
            note=note,
        )

    def cross_species(
        self,
        source_species: str,
        target_species: str,
        tf_family: str = "",
        metabolite: str = "",
    ) -> SpeciesKnowledgeResult:
        """Full cross-species inference: source → target."""
        source_info = self.lookup(source_species)
        target_info = self.lookup(target_species)

        known_mets = target_info.get("metabolites", []) if target_info else []
        known_paths = target_info.get("pathways", []) if target_info else []

        inferences: list[OrthologInference] = []

        # Check if both species share TF families
        if source_info and target_info:
            shared_tfs = set(
                tf.lower() for tf in source_info.get("tf_families", [])
            ) & set(
                tf.lower() for tf in target_info.get("tf_families", [])
            )
            for stf in shared_tfs:
                inferences.append(OrthologInference(
                    source_species=source_species,
                    target_species=target_species,
                    tf_family=stf.upper(),
                    metabolite=metabolite or "shared metabolites",
                    confidence="medium" if stf in ["myb", "bhlh", "nac"] else "low",
                    reasoning=(
                        f"{stf.upper()} TFs exist in both {source_species} and {target_species}. "
                        f"Cross-species functional conservation is possible but requires validation."
                    ),
                    caveats=[
                        "Functional validation needed",
                        "Check promoter motif conservation",
                    ],
                ))

        return SpeciesKnowledgeResult(
            query_species=target_species,
            query_metabolite=metabolite,
            known_pathways=known_paths,
            known_metabolites=known_mets,
            ortholog_inferences=inferences,
            note=(
                f"Cross-species inference: {source_species} → {target_species}. "
                f"Found {len(inferences)} potential conserved regulatory relationships."
            ) if inferences else f"No conserved regulatory patterns found between {source_species} and {target_species}.",
        )

    def list_species(self) -> list[str]:
        """List all species with knowledge base entries."""
        return sorted(self._species.keys())
