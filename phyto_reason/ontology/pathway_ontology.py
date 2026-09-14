"""
pathway_ontology.py — 通路本体。

定义 pathway → sub-pathway → enzyme 的层级关系。
"""

from __future__ import annotations

from phyto_reason.ontology.ontology_models import PathwayNode


PATHWAY_TREE: dict[str, PathwayNode] = {
    "secondary_metabolism": PathwayNode(
        node_id="secondary_metabolism", name="Secondary metabolism",
        children=["phenylpropanoid_biosynthesis", "alkaloid_biosynthesis",
                   "terpenoid_biosynthesis"],
    ),
    "phenylpropanoid_biosynthesis": PathwayNode(
        node_id="phenylpropanoid_biosynthesis", name="Phenylpropanoid biosynthesis",
        parent="secondary_metabolism",
        children=["flavonoid_biosynthesis", "lignin_biosynthesis",
                   "coumarin_biosynthesis", "phenolic_biosynthesis"],
        kegg_ids=["map00940"],
        enzymes=["PAL", "C4H", "4CL"],
        bottleneck_enzymes=["PAL", "C4H"],
    ),
    "flavonoid_biosynthesis": PathwayNode(
        node_id="flavonoid_biosynthesis", name="Flavonoid biosynthesis",
        parent="phenylpropanoid_biosynthesis",
        children=["flavonol_biosynthesis", "anthocyanin_biosynthesis",
                   "isoflavonoid_biosynthesis"],
        kegg_ids=["map00941"],
        enzymes=["CHS", "CHI", "F3H", "F3'H", "FLS"],
        bottleneck_enzymes=["CHS"],
    ),
    "flavonol_biosynthesis": PathwayNode(
        node_id="flavonol_biosynthesis", name="Flavonol biosynthesis",
        parent="flavonoid_biosynthesis",
        enzymes=["CHS", "CHI", "F3H", "F3'H", "FLS"],
        bottleneck_enzymes=["FLS"],
    ),
    "anthocyanin_biosynthesis": PathwayNode(
        node_id="anthocyanin_biosynthesis", name="Anthocyanin biosynthesis",
        parent="flavonoid_biosynthesis",
        enzymes=["CHS", "CHI", "F3H", "DFR", "ANS", "UFGT"],
        bottleneck_enzymes=["DFR", "ANS"],
    ),
    "lignin_biosynthesis": PathwayNode(
        node_id="lignin_biosynthesis", name="Lignin biosynthesis",
        parent="phenylpropanoid_biosynthesis",
        enzymes=["PAL", "C4H", "4CL", "HCT", "CCR", "CAD"],
        bottleneck_enzymes=["CCR", "CAD"],
    ),
    "alkaloid_biosynthesis": PathwayNode(
        node_id="alkaloid_biosynthesis", name="Alkaloid biosynthesis",
        parent="secondary_metabolism",
        children=["benzylisoquinoline_biosynthesis", "indole_alkaloid_biosynthesis"],
        kegg_ids=["map00950", "map00960"],
        enzymes=["TYDC", "NCS", "CYP450"],
        bottleneck_enzymes=["NCS", "TYDC"],
    ),
    "benzylisoquinoline_biosynthesis": PathwayNode(
        node_id="benzylisoquinoline_biosynthesis", name="Benzylisoquinoline alkaloid biosynthesis",
        parent="alkaloid_biosynthesis",
        enzymes=["TYDC", "NCS", "BBE", "CYP80B", "SDR", "OMT"],
        bottleneck_enzymes=["NCS", "CYP80B"],
    ),
    "terpenoid_biosynthesis": PathwayNode(
        node_id="terpenoid_biosynthesis", name="Terpenoid biosynthesis",
        parent="secondary_metabolism",
        children=["sesquiterpene_biosynthesis", "diterpene_biosynthesis",
                   "monoterpene_biosynthesis", "triterpene_biosynthesis"],
        kegg_ids=["map00900", "map00909"],
        enzymes=["HMGR", "DXS", "DXR", "TPS"],
        bottleneck_enzymes=["HMGR", "TPS"],
    ),
    "sesquiterpene_biosynthesis": PathwayNode(
        node_id="sesquiterpene_biosynthesis", name="Sesquiterpene biosynthesis",
        parent="terpenoid_biosynthesis",
        enzymes=["HMGR", "FPPS", "TPS"],
        bottleneck_enzymes=["TPS"],
    ),
    "diterpene_biosynthesis": PathwayNode(
        node_id="diterpene_biosynthesis", name="Diterpene biosynthesis",
        parent="terpenoid_biosynthesis",
        enzymes=["DXS", "DXR", "GGPPS", "TPS"],
        bottleneck_enzymes=["GGPPS", "TPS"],
    ),
}

# 别名映射: 常见 pathway 名称 → node_id
PATHWAY_ALIASES: dict[str, str] = {
    "flavonoid": "flavonoid_biosynthesis",
    "flavonoid biosynthesis": "flavonoid_biosynthesis",
    "flavonol": "flavonol_biosynthesis",
    "flavonol biosynthesis": "flavonol_biosynthesis",
    "anthocyanin": "anthocyanin_biosynthesis",
    "anthocyanin biosynthesis": "anthocyanin_biosynthesis",
    "alkaloid": "alkaloid_biosynthesis",
    "alkaloid biosynthesis": "alkaloid_biosynthesis",
    "terpenoid": "terpenoid_biosynthesis",
    "terpenoid biosynthesis": "terpenoid_biosynthesis",
    "phenylpropanoid": "phenylpropanoid_biosynthesis",
    "lignin": "lignin_biosynthesis",
    "lignin biosynthesis": "lignin_biosynthesis",
    "benzylisoquinoline": "benzylisoquinoline_biosynthesis",
    "benzylisoquinoline alkaloid": "benzylisoquinoline_biosynthesis",
}


class PathwayOntology:
    """通路本体查询。"""

    @staticmethod
    def get_node(node_id: str) -> PathwayNode | None:
        return PATHWAY_TREE.get(node_id)

    @staticmethod
    def resolve_alias(name: str) -> str | None:
        return PATHWAY_ALIASES.get(name.lower().strip())

    @staticmethod
    def get_parent_pathways(node_id: str) -> list[str]:
        parents = []
        current = PATHWAY_TREE.get(node_id)
        while current and current.parent:
            parents.append(current.parent)
            current = PATHWAY_TREE.get(current.parent)
        return parents

    @staticmethod
    def get_all_enzymes(node_id: str) -> list[str]:
        """获取通路 + 所有子通路的酶。"""
        enzymes = []
        node = PATHWAY_TREE.get(node_id)
        if node:
            enzymes.extend(node.enzymes)
            for child in node.children:
                child_enzymes = PathwayOntology.get_all_enzymes(child)
                for e in child_enzymes:
                    if e not in enzymes:
                        enzymes.append(e)
        return enzymes
