"""
ontology_resolver.py — 统一本体解析器。

输入任意代谢物/通路/TF名称，自动:
  1. 别名解析 (SynonymMapper)
  2. 代谢物层级定位 (MetaboliteOntology)
  3. 通路映射 (PathwayOntology)
  4. TF先验匹配 (TFPriorOntology)
  5. 返回完整 OntologyResult
"""

from __future__ import annotations

from phyto_reason.ontology.ontology_models import OntologyResult
from phyto_reason.ontology.synonym_mapper import SynonymMapper
from phyto_reason.ontology.metabolite_ontology import MetaboliteOntology
from phyto_reason.ontology.pathway_ontology import PathwayOntology
from phyto_reason.ontology.tf_prior_ontology import TFPriorOntology


class OntologyResolver:
    """统一本体解析器。"""

    def __init__(self) -> None:
        self.synonyms = SynonymMapper()
        self.metabolites = MetaboliteOntology()
        self.pathways = PathwayOntology()
        self.tf_prior = TFPriorOntology()

    def resolve(self, term: str) -> OntologyResult:
        """解析任意生物学术语。

        Args:
            term: 代谢物名、通路名、TF家族等

        Returns:
            OntologyResult — 包含完整的继承链、parents、pathways、TF prior scores
        """
        term_stripped = term.strip()

        # 0. Check for pathway-like terms FIRST (before synonym/metabolite lookup)
        if ("biosynthesis" in term_stripped.lower()
                or "pathway" in term_stripped.lower()
                or "metabolism" in term_stripped.lower()):
            from phyto_reason.ontology.pathway_ontology import PATHWAY_ALIASES
            pw_id = self.pathways.resolve_alias(term_stripped)
            if not pw_id:
                for alias_key, node_id in PATHWAY_ALIASES.items():
                    if term_stripped.lower() in alias_key or alias_key in term_stripped.lower():
                        pw_id = node_id
                        break
            if pw_id:
                node = self.pathways.get_node(pw_id)
                if node:
                    parents = self.pathways.get_parent_pathways(pw_id)
                    enzymes = self.pathways.get_all_enzymes(pw_id)
                    return OntologyResult(
                        original_term=term_stripped, resolved_node=pw_id,
                        node_type="pathway", parent_nodes=parents,
                        child_nodes=node.children, pathway_ids=node.kegg_ids,
                        inheritance_chain=[pw_id] + parents,
                    )

        # 1. 同义词解析
        canonical = self.synonyms.resolve(term_stripped)

        # 2. 代谢物查找
        if canonical and self.metabolites.get_node(canonical):
            node = self.metabolites.get_node(canonical)
            inheritance = self.metabolites.get_inheritance_chain(canonical)
            parents = self.metabolites.get_all_parents(canonical)
            children = self.metabolites.get_children(canonical)
            pathways = self.metabolites.get_pathways(canonical)

            # 计算 TF prior scores for all known TF families
            tf_scores = {}
            for family in self.tf_prior.list_families():
                score = self.tf_prior.get_best_score(family, canonical, parents)
                if score > 0:
                    tf_scores[family] = score

            return OntologyResult(
                original_term=term_stripped,
                resolved_node=canonical,
                node_type="metabolite",
                aliases=self.synonyms.get_aliases(canonical),
                parent_nodes=parents,
                child_nodes=children,
                pathway_ids=pathways,
                kegg_id=node.kegg_id,
                inheritance_chain=inheritance,
                tf_prior_scores=tf_scores,
            )

        # 3. 通路查找 (general)
        pw_id = self.pathways.resolve_alias(term_stripped)
        if pw_id:
            node = self.pathways.get_node(pw_id)
            if node:
                parents = self.pathways.get_parent_pathways(pw_id)
                enzymes = self.pathways.get_all_enzymes(pw_id)
                return OntologyResult(
                    original_term=term_stripped,
                    resolved_node=pw_id,
                    node_type="pathway",
                    parent_nodes=parents,
                    child_nodes=node.children,
                    pathway_ids=node.kegg_ids,
                    inheritance_chain=[pw_id] + parents,
                )

        # 4. TF family 查找
        family_upper = term_stripped.upper()
        if family_upper in self.tf_prior.list_families() or \
           family_upper.replace("-", "/") in self.tf_prior.list_families() or \
           family_upper.replace("/", "-") in self.tf_prior.list_families():
            return OntologyResult(
                original_term=term_stripped,
                resolved_node=family_upper,
                node_type="tf_family",
            )

        # 5. 模糊搜索 (fallback)
        fuzzy = self.synonyms.search(term_stripped)
        if fuzzy:
            result = self.resolve(fuzzy[0])
            result.error = "fuzzy_match"
            return result

        return OntologyResult(
            original_term=term_stripped,
            resolved_node=term_stripped,
            node_type="unknown",
            error=f"无法解析: '{term}'",
        )

    def resolve_with_prior(self, term: str, tf_family: str) -> OntologyResult:
        """解析并添加特定 TF 家族的 prior 评分。"""
        result = self.resolve(term)

        if result.node_type == "metabolite":
            score = self.tf_prior.get_best_score(
                tf_family, result.resolved_node, result.parent_nodes
            )
            result.tf_prior_scores[tf_family] = score

        return result

    def tf_prior_for_metabolite(self, term: str) -> dict[str, float]:
        """返回所有 TF 家族对该代谢物的 prior 评分。"""
        result = self.resolve(term)
        return result.tf_prior_scores
