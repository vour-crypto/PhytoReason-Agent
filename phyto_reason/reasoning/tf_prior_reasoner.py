"""
tf_prior_reasoner.py — TF家族 ↔ 代谢通路 先验知识推理。

内置来自已发表植物生物学研究的先验知识映射表。
所有规则均为 deterministic，不调用 LLM。
"""

from __future__ import annotations

from phyto_reason.reasoning.reasoning_result import PriorEvidence, ReasoningResult
from phyto_reason.knowledge.tf_knowledge_base import CANONICAL_KNOWLEDGE, CLASS_ALIASES


# ── 先验知识库 ──────────────────────────────────────────────────
# 结构: tf_family → target_metabolite → PriorEvidence
# v4.0: 数据源统一为 canonical knowledge base (knowledge/tf_knowledge_base.py)


def _build_prior_knowledge() -> dict[str, dict[str, PriorEvidence]]:
    """从 canonical knowledge base 构建 PRIOR_KNOWLEDGE dict。

    转换 TFMetaboliteRelation → PriorEvidence 并组织为 dict[TF][meta] 结构。
    """
    result: dict[str, dict[str, PriorEvidence]] = {}
    for rel in CANONICAL_KNOWLEDGE:
        tf = rel.tf_family
        meta = rel.metabolite_class
        if tf not in result:
            result[tf] = {}
        # 如果已有条目，保留分数更高的
        if meta in result[tf] and result[tf][meta].score >= rel.score:
            continue
        result[tf][meta] = PriorEvidence(
            strength=rel.strength if rel.strength != "inferred" else "weak",
            score=rel.score,
            description=rel.description,
            pmids=rel.pmids,
            species_scope=rel.species_scope,
        )
    # 为 specific_metabolite 也创建条目（用 CLASS_ALIASES 做反向映射）
    for rel in CANONICAL_KNOWLEDGE:
        if rel.specific_metabolite and rel.specific_metabolite not in result.get(rel.tf_family, {}):
            tf_dict = result.setdefault(rel.tf_family, {})
            tf_dict[rel.specific_metabolite] = PriorEvidence(
                strength=rel.strength if rel.strength != "inferred" else "weak",
                score=rel.score,
                description=rel.description,
                pmids=rel.pmids,
                species_scope=rel.species_scope,
            )
    # 为 CLASS_ALIASES 创建条目
    for alias, cls in CLASS_ALIASES.items():
        for tf, metas in result.items():
            if cls in metas and alias not in metas:
                metas[alias] = metas[cls]
    return result


PRIOR_KNOWLEDGE: dict[str, dict[str, PriorEvidence]] = _build_prior_knowledge()


# Legacy hardcoded data removed in v4.0. See knowledge/tf_knowledge_base.py.
# Remaining empty dict preserved for direct access pattern compatibility.
_EMPTY_PRIOR: dict[str, dict[str, PriorEvidence]] = {}

class TfPriorReasoner:
    """TF家族先验知识推理器。

    基于已发表的调控关系知识库，评估给定TF家族调控目标代谢通路的先验可能性。
    纯规则引擎，无LLM调用。
    """

    @staticmethod
    def _resolve_family(tf_family: str) -> str | None:
        """尝试多种格式匹配TF家族名。"""
        candidates = [tf_family, tf_family.upper(), tf_family.lower(), tf_family.capitalize()]
        for c in candidates:
            if c in PRIOR_KNOWLEDGE:
                return c
        # handle "/" vs "-" differences
        for c in candidates:
            alt = c.replace("/", "-")
            if alt in PRIOR_KNOWLEDGE:
                return alt
            alt = c.replace("-", "/")
            if alt in PRIOR_KNOWLEDGE:
                return alt
        return None

    @staticmethod
    def get_supported_metabolites(tf_family: str) -> list[str]:
        """获取给定TF家族已知调控的代谢物类别列表。"""
        key = TfPriorReasoner._resolve_family(tf_family)
        if key is None:
            return []
        entry = PRIOR_KNOWLEDGE.get(key, {})
        sorted_items = sorted(entry.items(), key=lambda x: x[1].score, reverse=True)
        return [m for m, _ in sorted_items]

    @staticmethod
    def reason(
        tf_family: str | None,
        target_metabolite: str,
    ) -> ReasoningResult:
        """评估TF家族调控目标代谢物的先验可能性。

        Args:
            tf_family: TF家族名称 (如 "MYB", "WRKY", "bHLH")
            target_metabolite: 目标代谢物名称或类别

        Returns:
            ReasoningResult: 包含先验分数、解释和文献证据
        """
        if not tf_family:
            return ReasoningResult(
                score=0.0, label="none",
                explanation="未提供TF家族信息，无法评估先验知识。",
            )

        family_key = TfPriorReasoner._resolve_family(tf_family)
        if family_key is None:
            return ReasoningResult(
                score=0.0, label="none",
                explanation=f"知识库中未找到TF家族 '{tf_family}' 的调控记录。",
            )
        normalized = target_metabolite.lower().strip()

        # 尝试精确匹配代谢物类别
        candidates = []
        family_data = PRIOR_KNOWLEDGE.get(family_key, {})

        for meta_class, evidence in family_data.items():
            if meta_class in normalized or normalized in meta_class:
                candidates.append((evidence.score, meta_class, evidence))
                break

        # 无精确匹配 → 检查部分匹配
        if not candidates:
            for meta_class, evidence in family_data.items():
                if any(word in normalized for word in meta_class.split("_")) or \
                   any(word in meta_class for word in normalized.split()):
                    candidates.append((evidence.score, meta_class, evidence))

        # 仍然无匹配 → 使用该家族的最佳general匹配
        if not candidates:
            family_data = PRIOR_KNOWLEDGE.get(family_key, {})
            if family_data:
                best = max(family_data.values(), key=lambda x: x.score)
                candidates.append((best.score * 0.5, "general", best))
            else:
                return ReasoningResult(
                    score=0.0, label="none",
                    explanation=f"知识库中未找到TF家族 '{tf_family}' 的调控记录。",
                )

        # 取最佳匹配
        best_score, best_class, best_evidence = max(candidates, key=lambda x: x[0])

        details = [
            f"TF家族 {tf_family} 调控 {best_class}：强度 {best_evidence.strength}",
        ]
        if best_evidence.pmids:
            details.append(f"文献支持: PMID {' '.join(best_evidence.pmids)}")
        if best_evidence.species_scope != "general":
            details.append(f"物种范围: {best_evidence.species_scope}")

        return ReasoningResult(
            score=best_score,
            label=best_evidence.strength,
            explanation=best_evidence.description,
            details=details,
            evidence_list=[best_evidence],
        )

    @staticmethod
    def batch_reason(
        tf_family: str | None,
        target_metabolites: list[str],
    ) -> ReasoningResult:
        """批量评估多个目标代谢物的先验证据。

        对所有代谢物取最高分。
        """
        if not target_metabolites:
            return ReasoningResult(
                score=0.0, label="none",
                explanation="未提供目标代谢物列表。",
            )

        best: ReasoningResult | None = None
        for meta in target_metabolites:
            result = TfPriorReasoner.reason(tf_family, meta)
            if best is None or result.score > best.score:
                best = result

        if best is None:
            return ReasoningResult(
                score=0.0, label="none",
                explanation=f"TF家族 '{tf_family}' 与任何目标代谢物均无已知调控关系。",
            )

        return best

    # ── Ontology-aware reasoning ──────────────────────────

    @staticmethod
    def reason_ontology(
        tf_family: str | None,
        target_metabolite: str,
    ) -> ReasoningResult:
        """Ontology-aware TF prior reasoning。

        当直接字符串匹配失败时, 通过生物学本体继承链
        (化合物 → 亚类 → 大类) 继承先验知识。

        例如:
          MYB ↔ quercetin
          → 解析: quercetin → flavonol → flavonoid
          → 继承: MYB ↔ flavonoid (strong score=0.90)
        """
        if not tf_family:
            return ReasoningResult(score=0.0, label="none",
                explanation="未提供TF家族信息，无法评估先验知识。")

        # Try direct matching first — only return if it's a strong exact match
        direct = TfPriorReasoner.reason(tf_family, target_metabolite)
        if direct.score >= 0.70:
            return direct

        # Ontology resolution
        try:
            from phyto_reason.ontology.ontology_resolver import OntologyResolver
            resolver = OntologyResolver()
            resolved = resolver.resolve_with_prior(target_metabolite, tf_family)

            if resolved.error and resolved.error != "fuzzy_match":
                return direct

            inheritance = resolved.inheritance_chain
            score = resolved.tf_prior_scores.get(tf_family, 0.0)

            if score > 0:
                chain_str = " → ".join(inheritance)
                explanation = (
                    f"通过本体继承链 '{chain_str}' 间接匹配: "
                    f"TF家族 {tf_family} 已知调控代谢物大类 "
                    f"'{resolved.parent_nodes[0] if resolved.parent_nodes else '?'}', "
                    f"该大类包含目标代谢物 {target_metabolite} ({resolved.resolved_node}). "
                )
                return ReasoningResult(
                    score=score,
                    label="strong" if score >= 0.7 else "moderate" if score >= 0.4 else "weak",
                    explanation=explanation,
                    details=[f"本体继承链: {chain_str}",
                              f"匹配大类: {resolved.parent_nodes}",
                              f"TF先验分数: {score:.2f}"],
                )
        except ImportError:
            pass

        return direct
