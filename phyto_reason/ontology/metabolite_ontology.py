"""
metabolite_ontology.py — 代谢物本体层级。

定义化合物 → 亚类 → 大类 → 通路的完整继承链。
"""

from __future__ import annotations

from phyto_reason.ontology.ontology_models import MetaboliteNode


# ── 代谢物本体层级 ─────────────────────────────────────────
# 结构: node_id → MetaboliteNode

METABOLITE_TREE: dict[str, MetaboliteNode] = {
    # ════════════════════════════════════════════════════════
    # Level 0: 最上级
    # ════════════════════════════════════════════════════════
    "secondary_metabolite": MetaboliteNode(
        node_id="secondary_metabolite", name="secondary metabolite",
        children=["phenylpropanoid", "alkaloid", "terpenoid"],
    ),
    "defense_metabolite": MetaboliteNode(
        node_id="defense_metabolite", name="defense metabolite",
        children=["alkaloid", "terpenoid", "phenolic", "phytoalexin"],
    ),

    # ════════════════════════════════════════════════════════
    # Level 1: 主要代谢物大类
    # ════════════════════════════════════════════════════════
    "phenylpropanoid": MetaboliteNode(
        node_id="phenylpropanoid", name="phenylpropanoid",
        parent="secondary_metabolite",
        children=["flavonoid", "lignin", "coumarin", "phenolic", "lignan"],
        pathway_ids=["map00940", "map00941", "map00942"],
    ),
    "alkaloid": MetaboliteNode(
        node_id="alkaloid", name="alkaloid",
        parent="secondary_metabolite",
        children=["benzylisoquinoline", "indole_alkaloid", "tropane_alkaloid",
                   "purine_alkaloid", "steroidal_alkaloid"],
        pathway_ids=["map00950", "map00960"],
    ),
    "terpenoid": MetaboliteNode(
        node_id="terpenoid", name="terpenoid",
        parent="secondary_metabolite",
        children=["monoterpene", "sesquiterpene", "diterpene",
                   "triterpene", "tetraterpene", "polyterpene"],
        pathway_ids=["map00900", "map00901", "map00902", "map00903",
                      "map00904", "map00905", "map00906"],
    ),

    # ════════════════════════════════════════════════════════
    # Level 2: 亚类
    # ════════════════════════════════════════════════════════
    "flavonoid": MetaboliteNode(
        node_id="flavonoid", name="flavonoid",
        parent="phenylpropanoid",
        children=["flavonol", "anthocyanin", "flavanone", "flavone",
                   "isoflavonoid", "flavanol", "chalcone"],
        pathway_ids=["map00941"],
    ),
    "flavonol": MetaboliteNode(
        node_id="flavonol", name="flavonol",
        parent="flavonoid",
        children=["quercetin", "kaempferol", "myricetin"],
    ),
    "anthocyanin": MetaboliteNode(
        node_id="anthocyanin", name="anthocyanin",
        parent="flavonoid",
        children=["cyanidin", "delphinidin", "pelargonidin", "peonidin"],
    ),
    "flavanone": MetaboliteNode(
        node_id="flavanone", name="flavanone",
        parent="flavonoid",
        children=["naringenin"],
    ),
    "flavone": MetaboliteNode(
        node_id="flavone", name="flavone",
        parent="flavonoid",
        children=["apigenin", "luteolin"],
    ),
    "flavanol": MetaboliteNode(
        node_id="flavanol", name="flavanol",
        parent="flavonoid",
        children=["catechin", "epicatechin"],
    ),
    "lignin": MetaboliteNode(
        node_id="lignin", name="lignin",
        parent="phenylpropanoid",
        children=["g_lignin", "s_lignin", "h_lignin"],
        pathway_ids=["map00940", "map00941"],
    ),
    "coumarin": MetaboliteNode(
        node_id="coumarin", name="coumarin",
        parent="phenylpropanoid",
    ),
    "phenolic": MetaboliteNode(
        node_id="phenolic", name="phenolic compound",
        parent="phenylpropanoid",
        children=["salicylic_acid", "caffeic_acid", "p_coumaric_acid",
                   "ferulic_acid", "sinapic_acid"],
    ),
    "benzylisoquinoline": MetaboliteNode(
        node_id="benzylisoquinoline", name="benzylisoquinoline alkaloid",
        parent="alkaloid",
        children=["berberine", "coptisine", "palmatine", "magnoflorine",
                   "morphine", "codeine", "sanguinarine", "chelerythrine",
                   "papaverine", "noscapine"],
    ),
    "indole_alkaloid": MetaboliteNode(
        node_id="indole_alkaloid", name="indole alkaloid",
        parent="alkaloid",
        children=["strychnine", "quinine", "vinblastine", "vincristine"],
    ),
    "tropane_alkaloid": MetaboliteNode(
        node_id="tropane_alkaloid", name="tropane alkaloid",
        parent="alkaloid",
        children=["atropine", "scopolamine", "cocaine"],
    ),
    "purine_alkaloid": MetaboliteNode(
        node_id="purine_alkaloid", name="purine alkaloid",
        parent="alkaloid",
        children=["caffeine", "theobromine"],
    ),
    "sesquiterpene": MetaboliteNode(
        node_id="sesquiterpene", name="sesquiterpene",
        parent="terpenoid",
        children=["artemisinin", "capsaicin"],
    ),
    "diterpene": MetaboliteNode(
        node_id="diterpene", name="diterpene",
        parent="terpenoid",
        children=["tanshinone", "taxol", "gibberellin"],
    ),
    "triterpene": MetaboliteNode(
        node_id="triterpene", name="triterpene",
        parent="terpenoid",
        children=["ginsenoside", "saponin"],
    ),
    "monoterpene": MetaboliteNode(
        node_id="monoterpene", name="monoterpene",
        parent="terpenoid",
        children=["menthol", "limonene", "linalool"],
    ),
    "tetraterpene": MetaboliteNode(
        node_id="tetraterpene", name="tetraterpene",
        parent="terpenoid",
        children=["beta_carotene", "lycopene", "lutein"],
    ),

    # ════════════════════════════════════════════════════════
    # Level 3: 具体化合物 (leaf nodes)
    # ════════════════════════════════════════════════════════
    "quercetin": MetaboliteNode(
        node_id="quercetin", name="quercetin",
        parent="flavonol", is_leaf=True, kegg_id="C00389",
        aliases=["3,3',4',5,7-pentahydroxyflavone", "quercetin aglycone"],
    ),
    "kaempferol": MetaboliteNode(
        node_id="kaempferol", name="kaempferol",
        parent="flavonol", is_leaf=True, kegg_id="C05903",
    ),
    "myricetin": MetaboliteNode(
        node_id="myricetin", name="myricetin",
        parent="flavonol", is_leaf=True,
    ),
    "naringenin": MetaboliteNode(
        node_id="naringenin", name="naringenin",
        parent="flavanone", is_leaf=True, kegg_id="C00509",
    ),
    "apigenin": MetaboliteNode(
        node_id="apigenin", name="apigenin",
        parent="flavone", is_leaf=True,
    ),
    "luteolin": MetaboliteNode(
        node_id="luteolin", name="luteolin",
        parent="flavone", is_leaf=True,
    ),
    "cyanidin": MetaboliteNode(
        node_id="cyanidin", name="cyanidin",
        parent="anthocyanin", is_leaf=True, kegg_id="C05905",
    ),
    "delphinidin": MetaboliteNode(
        node_id="delphinidin", name="delphinidin",
        parent="anthocyanin", is_leaf=True,
    ),
    "catechin": MetaboliteNode(
        node_id="catechin", name="catechin",
        parent="flavanol", is_leaf=True, kegg_id="C00114",
    ),
    "berberine": MetaboliteNode(
        node_id="berberine", name="berberine",
        parent="benzylisoquinoline", is_leaf=True, kegg_id="C00721",
    ),
    "coptisine": MetaboliteNode(
        node_id="coptisine", name="coptisine",
        parent="benzylisoquinoline", is_leaf=True,
    ),
    "palmatine": MetaboliteNode(
        node_id="palmatine", name="palmatine",
        parent="benzylisoquinoline", is_leaf=True,
    ),
    "magnoflorine": MetaboliteNode(
        node_id="magnoflorine", name="magnoflorine",
        parent="benzylisoquinoline", is_leaf=True,
    ),
    "morphine": MetaboliteNode(
        node_id="morphine", name="morphine",
        parent="benzylisoquinoline", is_leaf=True, kegg_id="C00876",
    ),
    "sanguinarine": MetaboliteNode(
        node_id="sanguinarine", name="sanguinarine",
        parent="benzylisoquinoline", is_leaf=True,
    ),
    "nicotine": MetaboliteNode(
        node_id="nicotine", name="nicotine",
        parent="pyridine_alkaloid", is_leaf=True, kegg_id="C00738",
    ),
    "caffeine": MetaboliteNode(
        node_id="caffeine", name="caffeine",
        parent="purine_alkaloid", is_leaf=True, kegg_id="C00148",
    ),
    "artemisinin": MetaboliteNode(
        node_id="artemisinin", name="artemisinin",
        parent="sesquiterpene", is_leaf=True, kegg_id="C06412",
    ),
    "tanshinone": MetaboliteNode(
        node_id="tanshinone", name="tanshinone",
        parent="diterpene", is_leaf=True,
    ),
    "ginsenoside": MetaboliteNode(
        node_id="ginsenoside", name="ginsenoside",
        parent="triterpene", is_leaf=True,
    ),
    "menthol": MetaboliteNode(
        node_id="menthol", name="menthol",
        parent="monoterpene", is_leaf=True, kegg_id="C00083",
    ),
    "limonene": MetaboliteNode(
        node_id="limonene", name="limonene",
        parent="monoterpene", is_leaf=True,
    ),
    "beta_carotene": MetaboliteNode(
        node_id="beta_carotene", name="beta-carotene",
        parent="tetraterpene", is_leaf=True, kegg_id="C02094",
    ),
    "lycopene": MetaboliteNode(
        node_id="lycopene", name="lycopene",
        parent="tetraterpene", is_leaf=True, kegg_id="C02020",
    ),
    "taxol": MetaboliteNode(
        node_id="taxol", name="taxol (paclitaxel)",
        parent="diterpene", is_leaf=True, kegg_id="C00940",
    ),
    "capsaicin": MetaboliteNode(
        node_id="capsaicin", name="capsaicin",
        parent="sesquiterpene", is_leaf=True, kegg_id="C00781",
    ),
    "salicylic_acid": MetaboliteNode(
        node_id="salicylic_acid", name="salicylic acid",
        parent="phenolic", is_leaf=True, kegg_id="C00805",
    ),
    "caffeic_acid": MetaboliteNode(
        node_id="caffeic_acid", name="caffeic acid",
        parent="phenolic", is_leaf=True, kegg_id="C00190",
    ),
}


class MetaboliteOntology:
    """代谢物本体 — 化合物层级查询。"""

    @staticmethod
    def get_node(node_id: str) -> MetaboliteNode | None:
        return METABOLITE_TREE.get(node_id)

    @staticmethod
    def get_inheritance_chain(node_id: str) -> list[str]:
        """从 leaf 到 root 的继承链。"""
        chain = [node_id]
        current = METABOLITE_TREE.get(node_id)
        while current and current.parent:
            chain.append(current.parent)
            current = METABOLITE_TREE.get(current.parent)
        return chain

    @staticmethod
    def get_all_parents(node_id: str) -> list[str]:
        """所有 parent nodes。"""
        chain = MetaboliteOntology.get_inheritance_chain(node_id)
        return chain[1:] if len(chain) > 1 else []

    @staticmethod
    def get_children(node_id: str) -> list[str]:
        node = METABOLITE_TREE.get(node_id)
        return node.children if node else []

    @staticmethod
    def get_leaves(parent_id: str) -> list[str]:
        """获取某个类别下所有 leaf 化合物。"""
        leaves = []
        node = METABOLITE_TREE.get(parent_id)
        if node and node.is_leaf:
            return [parent_id]
        if node:
            for child in node.children:
                leaves.extend(MetaboliteOntology.get_leaves(child))
        return leaves

    @staticmethod
    def get_pathways(node_id: str) -> list[str]:
        """获取从 leaf 向上继承的所有 pathway IDs。"""
        node = METABOLITE_TREE.get(node_id)
        if node and node.pathway_ids:
            return node.pathway_ids
        parents = MetaboliteOntology.get_all_parents(node_id)
        for p in parents:
            p_node = METABOLITE_TREE.get(p)
            if p_node and p_node.pathway_ids:
                return p_node.pathway_ids
        return []

    @staticmethod
    def list_all() -> list[str]:
        return list(METABOLITE_TREE.keys())
