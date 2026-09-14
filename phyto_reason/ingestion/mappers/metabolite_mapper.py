"""
metabolite_mapper.py — 代谢物名映射与标准化（纯函数）。

支持: KEGG / HMDB / MetaCyc hook。
"""

from __future__ import annotations

import re


class MetaboliteMappingResult:
    def __init__(self):
        self.mapped: dict[str, str] = {}
        self.unmapped: list[str] = []
        self.kegg_hits: dict[str, str] = {}
        self.hmdb_hits: dict[str, str] = {}
        self.warnings: list[str] = []


# 内置的同义词映射 — 扩展自 metabolite_ontology.py
SYNONYM_MAP: dict[str, str] = {
    "berberine": "berberine_alkaloid",
    "ber": "berberine_alkaloid",
    "nicotine": "nicotine_alkaloid",
    "nic": "nicotine_alkaloid",
    "anthocyanin": "anthocyanin_flavonoid",
    "cyanidin": "anthocyanin_flavonoid",
    "delphinidin": "anthocyanin_flavonoid",
    "artemisinin": "artemisinin_terpenoid",
    "arte": "artemisinin_terpenoid",
    "kaempferol": "kaempferol_flavonoid",
    "quercetin": "quercetin_flavonoid",
    "naringenin": "naringenin_flavonoid",
    "juglone": "juglone_naphthoquinone",
    "camptothecin": "camptothecin_alkaloid",
    "vinblastine": "vinblastine_alkaloid",
}

PATTERN_MAP: list[tuple[re.Pattern, str]] = [
    (re.compile(r"(isoquinoline|berberine|benzylisoquinoline)", re.IGNORECASE), "alkaloid"),
    (re.compile(r"(anthocyanin|flavonol|flavone|flavonoid)", re.IGNORECASE), "flavonoid"),
    (re.compile(r"(terpenoid|terpene|sesquiterpene)", re.IGNORECASE), "terpenoid"),
    (re.compile(r"(alkaloid)", re.IGNORECASE), "alkaloid"),
    (re.compile(r"(lignin|lignan)", re.IGNORECASE), "lignin"),
    (re.compile(r"(phenolic|phenylpropanoid)", re.IGNORECASE), "phenylpropanoid"),
]


def normalize_name(metabolite_name: str) -> str:
    """标准化代谢物名。"""
    name = str(metabolite_name).strip().lower()
    name = re.sub(r"[_\- ]+", "_", name)
    return name


def map_metabolites(
    metabolite_names: list[str],
) -> MetaboliteMappingResult:
    """映射代谢物名到规范名称。"""
    result = MetaboliteMappingResult()

    for mname in metabolite_names:
        norm = normalize_name(mname)
        if not norm:
            result.unmapped.append(mname)
            continue

        if norm in SYNONYM_MAP:
            result.mapped[mname] = SYNONYM_MAP[norm]
            result.kegg_hits[mname] = f"synonym:{SYNONYM_MAP[norm]}"
            continue

        mapped = False
        for pattern, category in PATTERN_MAP:
            if pattern.search(norm):
                canonical = f"{norm}_{category}" if category not in norm else norm
                result.mapped[mname] = canonical
                result.kegg_hits[mname] = f"pattern:{category}"
                mapped = True
                break

        if not mapped:
            # Keep original, mark as uncertain
            result.mapped[mname] = norm
            result.warnings.append(f"No annotation for '{mname}' — using raw name")

    return result
