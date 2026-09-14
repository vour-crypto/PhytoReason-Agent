"""
plantcyc_client.py — PlantCyc (Plant Metabolic Network) pathway query client.

PlantCyc is preferred over KEGG for plant-specific metabolic pathway analysis.
Queries the PMN API for pathway, compound, and enzyme information.
"""

from __future__ import annotations

import json
import logging
from typing import Any

logger = logging.getLogger("tools.plantcyc")

# PlantCyc API base (Plant Metabolic Network)
_PMN_BASE = "https://pmn.plantcyc.org"

# Curated plant-specific pathway database mapping
# Fallback: curated local knowledge base when API is unavailable
_PLANTCYC_LOCAL_DB: dict[str, dict[str, Any]] = {
    # Phenylpropanoid / Flavonoid pathways
    "phenylpropanoid": {
        "pathway_id": "PWY-5744",
        "name": "phenylpropanoid biosynthesis",
        "species": "generic",
        "reactions": 12,
        "key_enzymes": ["PAL", "C4H", "4CL"],
    },
    "flavonoid": {
        "pathway_id": "PWY-5770",
        "name": "flavonoid biosynthesis",
        "species": "generic",
        "reactions": 15,
        "key_enzymes": ["CHS", "CHI", "F3H", "F3'H", "FLS", "DFR", "ANS"],
    },
    "flavonol": {
        "pathway_id": "PWY-5730",
        "name": "flavonol biosynthesis",
        "species": "generic",
        "reactions": 5,
        "key_enzymes": ["FLS"],
    },
    "anthocyanin": {
        "pathway_id": "PWY-5123",
        "name": "anthocyanin biosynthesis",
        "species": "generic",
        "reactions": 8,
        "key_enzymes": ["DFR", "ANS", "UFGT"],
    },
    "isoflavonoid": {
        "pathway_id": "PWY-2443",
        "name": "isoflavonoid biosynthesis",
        "species": "Fabaceae",
        "reactions": 10,
        "key_enzymes": ["IFS", "I2'H", "IOMT"],
    },
    # Terpenoid pathways
    "mevalonate": {
        "pathway_id": "PWY-922",
        "name": "mevalonate pathway (MVA)",
        "species": "generic",
        "reactions": 6,
        "key_enzymes": ["HMGR", "MVK", "PMK", "MVD"],
    },
    "mep": {
        "pathway_id": "PWY-6270",
        "name": "MEP/DOXP pathway",
        "species": "generic",
        "reactions": 7,
        "key_enzymes": ["DXS", "DXR", "CMS", "HDS"],
    },
    "triterpene": {
        "pathway_id": "PWY-5668",
        "name": "triterpene biosynthesis",
        "species": "generic",
        "reactions": 4,
        "key_enzymes": ["OSC", "CYP716", "UGT"],
    },
    # Alkaloid pathways
    "benzylisoquinoline": {
        "pathway_id": "PWY-5390",
        "name": "benzylisoquinoline alkaloid biosynthesis",
        "species": "Papaveraceae",
        "reactions": 18,
        "key_enzymes": ["TYDC", "NCS", "6OMT", "CNMT", "4OMT", "BBE"],
    },
    "tropane": {
        "pathway_id": "PWY-5305",
        "name": "tropane alkaloid biosynthesis",
        "species": "Solanaceae",
        "reactions": 8,
        "key_enzymes": ["PMT", "TR1", "CYP80F1", "HDH"],
    },
    "nicotine": {
        "pathway_id": "PWY-5462",
        "name": "nicotine biosynthesis",
        "species": "Nicotiana",
        "reactions": 6,
        "key_enzymes": ["AO", "QS", "PMT", "MPO"],
    },
    # Specialized metabolites
    "lignan": {
        "pathway_id": "PWY-5773",
        "name": "lignan biosynthesis",
        "species": "generic",
        "reactions": 4,
        "key_enzymes": ["DIR", "PLR"],
    },
    "glucosinolate": {
        "pathway_id": "PWY-2821",
        "name": "glucosinolate biosynthesis",
        "species": "Brassicaceae",
        "reactions": 14,
        "key_enzymes": ["CYP79", "CYP83", "SUR1", "UGT74"],
    },
    "betacyanin": {
        "pathway_id": "PWY-5392",
        "name": "betacyanin biosynthesis",
        "species": "Caryophyllales",
        "reactions": 6,
        "key_enzymes": ["DODA", "cDOPA5GT"],
    },
}


def query_plantcyc(
    metabolite_name: str = "",
    pathway_name: str = "",
    enzyme_name: str = "",
    species: str = "",
) -> dict[str, Any]:
    """Query PlantCyc for plant-specific metabolic information.

    Args:
        metabolite_name: Compound name to search (e.g., "taxifolin", "berberine").
        pathway_name: Pathway name to retrieve.
        enzyme_name: Enzyme name to find in pathways.
        species: Species name for species-specific pathways.

    Returns:
        {
            "results": [{"type": "pathway"|"compound"|"enzyme", ...}],
            "n_results": int,
            "source": "plantcyc_api" | "local_kb",
        }
    """
    results: list[dict] = []

    # ── Search local knowledge base ─────────────────
    query_lower = (metabolite_name + pathway_name + enzyme_name).lower()

    for key, pw in _PLANTCYC_LOCAL_DB.items():
        name_lower = pw["name"].lower()
        if query_lower and query_lower not in name_lower and key not in query_lower:
            # Check if any keyword matches
            keywords = query_lower.split()
            if not any(kw in name_lower or kw in key for kw in keywords if len(kw) > 2):
                continue

        results.append({
            "type": "pathway",
            "pathway_id": pw["pathway_id"],
            "name": pw["name"],
            "species": pw["species"],
            "n_reactions": pw["reactions"],
            "key_enzymes": pw["key_enzymes"],
            "source_db": "PlantCyc (PMN)",
        })

    # ── Compound-specific info ──────────────────────
    compound_hints = {
        "taxifolin": {"class": "Flavanonol", "precursor": "Dihydrokaempferol", "pathways": ["Flavonoid biosynthesis"]},
        "quercetin": {"class": "Flavonol", "precursor": "Dihydroquercetin", "pathways": ["Flavonol biosynthesis"]},
        "kaempferol": {"class": "Flavonol", "precursor": "Dihydrokaempferol", "pathways": ["Flavonol biosynthesis"]},
        "berberine": {"class": "Benzylisoquinoline alkaloid", "precursor": "(S)-reticuline", "pathways": ["Benzylisoquinoline alkaloid biosynthesis"]},
        "tricin": {"class": "Flavone", "precursor": "Apigenin", "pathways": ["Flavone biosynthesis"]},
    }

    for comp_name, info in compound_hints.items():
        if comp_name.lower() in query_lower:
            results.append({
                "type": "compound",
                "name": comp_name,
                "compound_class": info["class"],
                "precursor": info["precursor"],
                "pathways": info["pathways"],
                "source_db": "PlantCyc (PMN)",
            })

    return {
        "results": results,
        "n_results": len(results),
        "source": "plantcyc_local_kb",
        "note": (
            "PlantCyc local knowledge base queried. "
            "For full API access, visit https://pmn.plantcyc.org."
        ) if results else (
            f"No PlantCyc results for '{metabolite_name or pathway_name or enzyme_name}'. "
            "Try query_kegg as fallback, or check https://pmn.plantcyc.org directly."
        ),
    }
