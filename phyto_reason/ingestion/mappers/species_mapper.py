"""
species_mapper.py — 物种名标准化（纯函数）。
"""

from __future__ import annotations


class SpeciesMappingResult:
    def __init__(self):
        self.canonical: str = ""
        self.confidence: float = 1.0
        self.warnings: list[str] = []


SPECIES_MAP = {
    "\u4e24\u9762\u9488": "zanthoxylum_nitidum",
    "\u4e24\u9762\u9488(zanthoxylum nitidum)": "zanthoxylum_nitidum",
    "zanthoxylum nitidum": "zanthoxylum_nitidum",
    "arabidopsis": "arabidopsis_thaliana",
    "arabidopsis thaliana": "arabidopsis_thaliana",
    "ath": "arabidopsis_thaliana",
    "rice": "oryza_sativa",
    "oryza sativa": "oryza_sativa",
    "osa": "oryza_sativa",
    "tomato": "solanum_lycopersicum",
    "solanum lycopersicum": "solanum_lycopersicum",
    "sly": "solanum_lycopersicum",
    "tobacco": "nicotiana_tabacum",
    "nicotiana tabacum": "nicotiana_tabacum",
    "nta": "nicotiana_tabacum",
    "maize": "zea_mays",
    "zea mays": "zea_mays",
    "zma": "zea_mays",
    "soybean": "glycine_max",
    "glycine max": "glycine_max",
    "gma": "glycine_max",
    "grape": "vitis_vinifera",
    "vitis vinifera": "vitis_vinifera",
    "vvi": "vitis_vinifera",
}


def map_species(species: str) -> SpeciesMappingResult:
    """标准化物种名。"""
    result = SpeciesMappingResult()
    key = species.strip().lower().replace("-", "_")
    canonical = SPECIES_MAP.get(key)
    if canonical:
        result.canonical = canonical
        result.confidence = 1.0
    else:
        result.canonical = key
        result.confidence = 0.5
        result.warnings.append(f"Species '{species}' not in canonical list — using raw name")
    return result
