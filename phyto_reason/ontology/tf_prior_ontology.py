"""
tf_prior_ontology.py — TF 先验本体。

TF family ↔ pathway 映射，支持 ontology 继承。
当查询 MYB ↔ quercetin 时:
  quercetin → flavonol → flavonoid (ontology 继承)
  MYB ↔ flavonoid → strong (prior 继承)
"""

from __future__ import annotations

from phyto_reason.ontology.ontology_models import TFPriorEntry
from phyto_reason.knowledge.tf_knowledge_base import CANONICAL_KNOWLEDGE, CLASS_ALIASES


# ── TF 先验本体 (映射到 ontology node_id) ────────────────
# v4.0: 数据源统一为 canonical knowledge base
# target 使用 metabolite_ontology 中的 node_id
# 子类自动继承父类的 TF prior


def _build_prior_ontology() -> dict[str, dict[str, TFPriorEntry]]:
    """从 canonical knowledge base 构建 TF_PRIOR_ONTOLOGY。

    转换 TFMetaboliteRelation → TFPriorEntry。
    """
    result: dict[str, dict[str, TFPriorEntry]] = {}
    for rel in CANONICAL_KNOWLEDGE:
        tf = rel.tf_family
        target = rel.metabolite_class if rel.metabolite_class != "defense_metabolite" else "defense_metabolite"
        if tf not in result:
            result[tf] = {}
        if target in result[tf] and result[tf][target].score >= rel.score:
            continue
        result[tf][target] = TFPriorEntry(
            tf_family=tf,
            target_pathway=target,
            strength=rel.strength if rel.strength != "inferred" else "weak",
            score=rel.score,
            pmids=rel.pmids,
            description=rel.description,
            inheritable=rel.strength != "inferred",
        )
    # CLASS_ALIASES entries
    for alias, cls in CLASS_ALIASES.items():
        for tf, metas in result.items():
            if cls in metas and alias not in metas:
                entry = metas[cls]
                result[tf][alias] = TFPriorEntry(
                    tf_family=tf,
                    target_pathway=alias,
                    strength=entry.strength,
                    score=entry.score,
                    pmids=entry.pmids,
                    description=entry.description,
                    inheritable=entry.inheritable,
                )
    return result


TF_PRIOR_ONTOLOGY: dict[str, dict[str, TFPriorEntry]] = _build_prior_ontology()


# Legacy hardcoded data removed in v4.0. See knowledge/tf_knowledge_base.py.
_EMPTY_ONTOLOGY: dict[str, dict[str, TFPriorEntry]] = {}

class TFPriorOntology:
    """TF 先验本体 — 支持 ontology-aware 继承查询。"""

    @staticmethod
    def get_direct(tf_family: str, target: str) -> TFPriorEntry | None:
        """直接查询 TF ↔ target 的先验。"""
        family_data = TF_PRIOR_ONTOLOGY.get(tf_family)
        if family_data:
            return family_data.get(target)
        return None

    @staticmethod
    def get_inherited(tf_family: str, target: str,
                      parent_targets: list[str]) -> list[tuple[str, TFPriorEntry, str]]:
        """查询 TF ↔ target 及其 ontology 父类的先验。

        Returns:
            [(matched_target, entry, inheritance_type), ...]
            其中 inheritance_type: "direct" | "inherited_from_parent"
        """
        results: list[tuple[str, TFPriorEntry, str]] = []

        direct = TFPriorOntology.get_direct(tf_family, target)
        if direct:
            results.append((target, direct, "direct"))

        for parent in parent_targets:
            inherited = TFPriorOntology.get_direct(tf_family, parent)
            if inherited and inherited.inheritable:
                results.append((parent, inherited, f"inherited_from_{parent}"))

        return results

    @staticmethod
    def get_best_score(tf_family: str, target: str,
                       parent_targets: list[str] | None = None) -> float:
        """获取最佳先验分数 (含继承)。"""
        entries = TFPriorOntology.get_inherited(
            tf_family, target, parent_targets or []
        )
        if not entries:
            return 0.0
        return max(e[1].score for e in entries)

    @staticmethod
    def get_best_entry(tf_family: str, target: str,
                       parent_targets: list[str] | None = None) -> TFPriorEntry | None:
        entries = TFPriorOntology.get_inherited(
            tf_family, target, parent_targets or []
        )
        if not entries:
            return None
        best = max(entries, key=lambda e: e[1].score)
        return best[1]

    @staticmethod
    def list_families() -> list[str]:
        return list(TF_PRIOR_ONTOLOGY.keys())
