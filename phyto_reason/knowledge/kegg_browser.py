"""
kegg_browser.py — KEGG database browser for zero-data mode.

Provides structured pathway and compound lookups without requiring user data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger("kegg_browser")


@dataclass
class PathwayInfo:
    """Structured KEGG pathway information."""
    pathway_id: str = ""
    pathway_name: str = ""
    pathway_class: str = ""
    compound_name: str = ""
    enzyme_genes: list[str] = field(default_factory=list)
    enzyme_names: list[str] = field(default_factory=list)
    related_pathways: list[str] = field(default_factory=list)
    source_url: str = ""


@dataclass
class KEGGBrowseResult:
    query: str
    pathways: list[PathwayInfo] = field(default_factory=list)
    note: str = ""


# Built-in pathway knowledge (fallback when KEGG API unavailable)
# Covers common plant secondary metabolism pathways
BUILTIN_PATHWAYS: dict[str, PathwayInfo] = {
    "berberine": PathwayInfo(
        pathway_id="map00950",
        pathway_name="Isoquinoline alkaloid biosynthesis",
        pathway_class="Alkaloid biosynthesis",
        compound_name="berberine",
        enzyme_genes=[
            "BBE", "CYP80G2", "SMT", "CNMT", "6OMT", "4'OMT",
            "CYP719A1", "TNMT", "CYP80B1", "NMCH",
        ],
        enzyme_names=[
            "berberine bridge enzyme",
            "corytuberine synthase",
            "scoulerine 9-O-methyltransferase",
            "coclaurine N-methyltransferase",
            "norcoclaurine 6-O-methyltransferase",
            "3'-hydroxy-N-methylcoclaurine 4'-O-methyltransferase",
        ],
        source_url="https://www.genome.jp/pathway/map00950",
    ),
    "nicotine": PathwayInfo(
        pathway_id="map00960",
        pathway_name="Tropane, piperidine and pyridine alkaloid biosynthesis",
        pathway_class="Alkaloid biosynthesis",
        compound_name="nicotine",
        enzyme_genes=[
            "PMT", "MPO", "QPT", "A622", "BBL",
            "ODC", "ADC", "DAO",
        ],
        enzyme_names=[
            "putrescine N-methyltransferase",
            "N-methylputrescine oxidase",
            "quinolinate phosphoribosyltransferase",
        ],
        source_url="https://www.genome.jp/pathway/map00960",
    ),
    "anthocyanin": PathwayInfo(
        pathway_id="map00942",
        pathway_name="Anthocyanin biosynthesis",
        pathway_class="Flavonoid biosynthesis",
        compound_name="anthocyanin",
        enzyme_genes=[
            "CHS", "CHI", "F3H", "F3'H", "F3'5'H",
            "DFR", "ANS", "UFGT", "OMT",
        ],
        enzyme_names=[
            "chalcone synthase",
            "chalcone isomerase",
            "flavanone 3-hydroxylase",
            "dihydroflavonol 4-reductase",
            "anthocyanidin synthase",
        ],
        source_url="https://www.genome.jp/pathway/map00942",
    ),
    "flavonoid": PathwayInfo(
        pathway_id="map00941",
        pathway_name="Flavonoid biosynthesis",
        pathway_class="Phenylpropanoid biosynthesis",
        compound_name="flavonoid",
        enzyme_genes=[
            "PAL", "C4H", "4CL", "CHS", "CHI", "F3H",
            "FNS", "FLS", "DFR", "LAR", "ANR",
        ],
        enzyme_names=[
            "phenylalanine ammonia-lyase",
            "cinnamate 4-hydroxylase",
            "4-coumarate:CoA ligase",
        ],
        source_url="https://www.genome.jp/pathway/map00941",
    ),
    "artemisinin": PathwayInfo(
        pathway_id="map00909",
        pathway_name="Sesquiterpenoid and triterpenoid biosynthesis",
        pathway_class="Terpenoid biosynthesis",
        compound_name="artemisinin",
        enzyme_genes=[
            "ADS", "CYP71AV1", "DBR2", "ALDH1", "HMGR",
            "FPS", "SQS",
        ],
        enzyme_names=[
            "amorpha-4,11-diene synthase",
            "amorpha-4,11-diene monooxygenase",
            "artemisinic aldehyde reductase",
        ],
        source_url="https://www.genome.jp/pathway/map00909",
    ),
    "phenylpropanoid": PathwayInfo(
        pathway_id="map00940",
        pathway_name="Phenylpropanoid biosynthesis",
        pathway_class="General metabolism",
        compound_name="phenylpropanoid",
        enzyme_genes=["PAL", "C4H", "4CL", "CCR", "CAD", "COMT", "CCoAOMT"],
        enzyme_names=[
            "phenylalanine ammonia-lyase",
            "cinnamate 4-hydroxylase",
            "4-coumarate:CoA ligase",
        ],
        source_url="https://www.genome.jp/pathway/map00940",
    ),
}


class KEGGBrowser:
    """KEGG database browser with built-in pathway knowledge fallback."""

    def __init__(self) -> None:
        self._builtin = BUILTIN_PATHWAYS

    def browse(self, compound_name: str) -> KEGGBrowseResult:
        """Look up pathway information for a compound.

        Tries KEGG API first, falls back to built-in knowledge.
        """
        compound_lower = compound_name.lower().strip()

        # 1. Try KEGG API
        api_result = self._query_kegg_api(compound_name)
        if api_result and api_result.pathways:
            return api_result

        # 2. Fall back to built-in knowledge
        builtin = self._lookup_builtin(compound_lower)
        if builtin:
            return KEGGBrowseResult(
                query=compound_name,
                pathways=[builtin],
                note=f"Built-in knowledge for {compound_name} (KEGG API unavailable or no results)",
            )

        # 3. Try fuzzy match
        fuzzy = self._fuzzy_match(compound_lower)
        if fuzzy:
            return KEGGBrowseResult(
                query=compound_name,
                pathways=[fuzzy],
                note=f"Fuzzy match: {fuzzy.pathway_name} (no exact match for '{compound_name}')",
            )

        return KEGGBrowseResult(
            query=compound_name,
            note=f"No pathway information found for '{compound_name}'. Try a different metabolite name.",
        )

    def list_known_compounds(self) -> list[str]:
        """List all compounds with built-in pathway knowledge."""
        return sorted(self._builtin.keys())

    def _lookup_builtin(self, name: str) -> PathwayInfo | None:
        """Look up a compound in the built-in knowledge base."""
        name_clean = name.replace("-", "_").replace(" ", "_")
        for key, info in self._builtin.items():
            if key in name_clean or name_clean in key:
                return info
        return None

    def _fuzzy_match(self, name: str) -> PathwayInfo | None:
        """Fuzzy match compound name to pathway class."""
        class_keywords = {
            "alkaloid": "berberine",
            "isoquinoline": "berberine",
            "flavonoid": "flavonoid",
            "anthocyanin": "anthocyanin",
            "terpenoid": "artemisinin",
            "sesquiterpene": "artemisinin",
            "phenylpropanoid": "phenylpropanoid",
            "nicotine": "nicotine",
            "tropane": "nicotine",
        }
        for kw, ref in class_keywords.items():
            if kw in name:
                return self._builtin.get(ref)
        return None

    def _query_kegg_api(self, compound_name: str) -> KEGGBrowseResult | None:
        """Try to query KEGG API directly."""
        try:
            from phyto_reason.tools.pathway.kegg_client import query_kegg

            result = query_kegg(compound_name)
            if not result:
                return None

            pathways = []
            if isinstance(result, dict):
                pi = PathwayInfo(
                    pathway_name=result.get("pathway", ""),
                    pathway_class=result.get("pathway_class", ""),
                    compound_name=compound_name,
                    enzyme_genes=result.get("enzymes", []),
                )
                if pi.pathway_name:
                    pathways.append(pi)
            elif isinstance(result, str):
                pathways.append(PathwayInfo(
                    pathway_name=result[:200],
                    compound_name=compound_name,
                ))

            return KEGGBrowseResult(query=compound_name, pathways=pathways)

        except Exception as e:
            logger.debug(f"KEGG API unavailable: {e}")
            return None
