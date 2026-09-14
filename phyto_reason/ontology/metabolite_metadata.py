"""
metabolite_ontology.py — 代谢物本体的结构化表示。

reasoning 必须基于 ontology 而不是字符串匹配。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class MetaboliteOntology(BaseModel):
    """代谢物本体 — 替代简单的字符串比较。"""
    metabolite: str

    superclass: str = ""
    subclass: str = ""

    precursor_pathways: list[str] = Field(default_factory=list)

    competing_branches: list[str] = Field(default_factory=list)

    ecological_roles: list[str] = Field(default_factory=list)

    stress_associations: list[str] = Field(default_factory=list)

    tissue_accumulation: list[str] = Field(default_factory=list)

    transport_modes: list[str] = Field(default_factory=list)

    known_tf_families: list[str] = Field(default_factory=list)


# 内置本体库 — 可扩展
BUILTIN_ONTOLOGY: dict[str, MetaboliteOntology] = {
    "berberine": MetaboliteOntology(
        metabolite="berberine",
        superclass="alkaloid",
        subclass="benzylisoquinoline",
        precursor_pathways=["tyrosine", "dopamine"],
        competing_branches=["sanguinarine", "morphine"],
        ecological_roles=["antimicrobial", "antifeedant"],
        stress_associations=["pathogen", "wounding"],
        tissue_accumulation=["root", "rhizome"],
        transport_modes=["abc_transporter", "vesicle"],
        known_tf_families=["WRKY", "MYB", "bHLH", "ERF"],
    ),
    "nicotine": MetaboliteOntology(
        metabolite="nicotine",
        superclass="alkaloid",
        subclass="pyridine",
        precursor_pathways=["nicotinic_acid", "putrescine"],
        competing_branches=["nornicotine", "anatabine"],
        ecological_roles=["antifeedant", "defense"],
        stress_associations=["herbivory", "wounding", "JA"],
        tissue_accumulation=["root", "leaf"],
        transport_modes=["abc_transporter", "xylem"],
        known_tf_families=["ERF", "MYC2", "bHLH", "WRKY"],
    ),
    "anthocyanin": MetaboliteOntology(
        metabolite="anthocyanin",
        superclass="flavonoid",
        subclass="anthocyanin",
        precursor_pathways=["phenylalanine", "malonyl_coa"],
        competing_branches=["proanthocyanidin", "flavonol"],
        ecological_roles=["pigmentation", "UV_protection", "pollination"],
        stress_associations=["light", "cold", "nutrient"],
        tissue_accumulation=["leaf", "fruit", "flower", "seed"],
        transport_modes=["mrp_transporter", "vesicle"],
        known_tf_families=["MYB", "bHLH", "WD40", "WRKY"],
    ),
    "artemisinin": MetaboliteOntology(
        metabolite="artemisinin",
        superclass="terpenoid",
        subclass="sesquiterpene_lactone",
        precursor_pathways=["farnesyl_diphosphate", "MVA", "MEP"],
        competing_branches=["sterol", "triterpenoid"],
        ecological_roles=["antimalarial", "defense"],
        stress_associations=["oxidative", "pathogen"],
        tissue_accumulation=["leaf_trichome", "gland"],
        transport_modes=["unknown"],
        known_tf_families=["ERF", "WRKY", "MYB", "bZIP"],
    ),
}


def resolve_ontology(metabolite_name: str) -> MetaboliteOntology | None:
    """按名称解析代谢物本体，支持模糊匹配。"""
    name_lower = metabolite_name.lower().replace("-", "_").replace(" ", "_")
    for key, onto in BUILTIN_ONTOLOGY.items():
        if key in name_lower or name_lower in key:
            return onto
    for key, onto in BUILTIN_ONTOLOGY.items():
        if onto.superclass in name_lower or onto.subclass in name_lower:
            return onto
    return None


def get_competing_hypotheses_metabolites(metabolite: str) -> list[str]:
    """基于本体生成竞争假设的代谢物目标列表。"""
    onto = resolve_ontology(metabolite)
    if onto is None:
        return []
    return onto.competing_branches + [onto.superclass]
