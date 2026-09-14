"""
synonym_mapper.py — 代谢物同义词/别名映射。

支持:
  - 通用名 ↔ 学名
  - KEGG 交叉引用
  - 大小写归一化
  - 多语言别名
"""

from __future__ import annotations


SYNONYM_MAP: dict[str, list[str]] = {
    # 黄酮类
    "quercetin": ["quercetin", "Quercetin", "3,3',4',5,7-pentahydroxyflavone",
                   "quercetin_aglycone", "C00389", "Q4951"],
    "kaempferol": ["kaempferol", "Kaempferol", "3,4',5,7-tetrahydroxyflavone",
                    "kampferol", "C05903"],
    "myricetin": ["myricetin", "Myricetin", "3,3',4',5,5',7-hexahydroxyflavone"],
    "naringenin": ["naringenin", "Naringenin", "4',5,7-trihydroxyflavanone",
                    "C00509"],
    "apigenin": ["apigenin", "Apigenin", "4',5,7-trihydroxyflavone"],
    "luteolin": ["luteolin", "Luteolin", "3',4',5,7-tetrahydroxyflavone"],
    "cyanidin": ["cyanidin", "Cyanidin", "C05905"],
    "delphinidin": ["delphinidin", "Delphinidin"],
    "catechin": ["catechin", "Catechin", "C00114"],
    "epicatechin": ["epicatechin", "Epicatechin"],
    "proanthocyanidin": ["proanthocyanidin", "condensed tannin", "PAs"],

    # 生物碱类
    "berberine": ["berberine", "Berberine", "C00721", "BBR"],
    "coptisine": ["coptisine", "Coptisine"],
    "palmatine": ["palmatine", "Palmatine"],
    "magnoflorine": ["magnoflorine", "Magnoflorine"],
    "nicotine": ["nicotine", "Nicotine", "C00738"],
    "morphine": ["morphine", "Morphine", "C00876"],
    "codeine": ["codeine", "Codeine", "C00477"],
    "caffeine": ["caffeine", "Caffeine", "C00148"],
    "capsaicin": ["capsaicin", "Capsaicin", "C00781"],
    "sanguinarine": ["sanguinarine", "Sanguinarine"],
    "chelerythrine": ["chelerythrine", "Chelerythrine"],
    "strychnine": ["strychnine", "Strychnine"],
    "quinine": ["quinine", "Quinine", "C00366"],

    # 萜类
    "artemisinin": ["artemisinin", "Artemisinin", "artemisinine", "C06412"],
    "tanshinone": ["tanshinone", "Tanshinone", "tanshinone_iia"],
    "ginsenoside": ["ginsenoside", "Ginsenoside"],
    "menthol": ["menthol", "Menthol", "C00083"],
    "limonene": ["limonene", "Limonene", "C00114"],
    "beta_carotene": ["beta_carotene", "β-carotene", "beta-carotene", "C02094"],
    "lycopene": ["lycopene", "Lycopene", "C02020"],
    "taxol": ["taxol", "paclitaxel", "Paclitaxel", "C00940"],

    # 酚酸类
    "salicylic_acid": ["salicylic acid", "salicylate", "SA", "C00805"],
    "caffeic_acid": ["caffeic acid", "caffeate", "C00190"],
    "p_coumaric_acid": ["p-coumaric acid", "4-coumaric acid", "C00811"],
    "ferulic_acid": ["ferulic acid", "ferulate", "C00194"],
    "sinapic_acid": ["sinapic acid", "sinapate", "C00195"],

    # 通路/类别名
    "flavonoid": ["flavonoid", "flavonoids", "flavonoid compound",
                   "flavonoid biosynthesis"],
    "flavonol": ["flavonol", "flavonols", "flavonol glycoside"],
    "anthocyanin": ["anthocyanin", "anthocyanins", "anthocyanidin"],
    "isoflavonoid": ["isoflavonoid", "isoflavone", "isoflavones"],
    "alkaloid": ["alkaloid", "alkaloids", "alkaloid compound"],
    "benzylisoquinoline": ["benzylisoquinoline", "benzylisoquinoline alkaloid", "BIA"],
    "terpenoid": ["terpenoid", "terpenoids", "terpene", "terpenes",
                   "isoprenoid", "isoprenoids"],
    "monoterpene": ["monoterpene", "monoterpenoid", "C10 terpene"],
    "sesquiterpene": ["sesquiterpene", "sesquiterpenoid", "C15 terpene"],
    "diterpene": ["diterpene", "diterpenoid", "C20 terpene"],
    "triterpene": ["triterpene", "triterpenoid", "C30 terpene"],
    "lignin": ["lignin", "lignins", "lignification"],
    "phenylpropanoid": ["phenylpropanoid", "phenylpropanoids",
                         "phenylpropanoid biosynthesis"],
    "phenolic": ["phenolic", "phenolics", "phenolic compound",
                  "phenolic acid"],
    "coumarin": ["coumarin", "coumarins"],
    "lignan": ["lignan", "lignans"],
    "defense": ["defense", "defence", "defense response",
                 "plant defense", "defense metabolite", "phytoalexin"],
}


class SynonymMapper:
    """同义词映射器。"""

    @staticmethod
    def normalize(term: str) -> str:
        """归一化: lowercase + strip + hyphen→underscore."""
        result = term.strip().lower()
        result = result.replace("β-", "beta_").replace("β", "beta")
        result = result.replace("-", "_").replace(" ", "_")
        return result

    @staticmethod
    def resolve(term: str) -> str | None:
        """将任意别名解析为规范 ID。"""
        normalized = SynonymMapper.normalize(term)
        for canonical, aliases in SYNONYM_MAP.items():
            for alias in aliases:
                if SynonymMapper.normalize(alias) == normalized:
                    return canonical
        return None

    @staticmethod
    def get_aliases(canonical: str) -> list[str]:
        """获取规范 ID 的所有别名。"""
        return SYNONYM_MAP.get(canonical, [])

    @staticmethod
    def search(partial: str) -> list[str]:
        """模糊搜索 — 返回所有包含 partial 的规范 ID。"""
        normalized = SynonymMapper.normalize(partial)
        matches = []
        for canonical, aliases in SYNONYM_MAP.items():
            if normalized in canonical.lower():
                matches.append(canonical)
            else:
                for alias in aliases:
                    if normalized in alias.lower():
                        matches.append(canonical)
                        break
        return matches
