"""
ortholog_reasoner.py — 同源基因推理（进化保守性评估）。

通过拟南芥同源基因的已知调控关系，评估目标物种中
候选TF的调控功能保守性。
"""

from __future__ import annotations

from phyto_reason.reasoning.reasoning_result import PriorEvidence, ReasoningResult


# ── 拟南芥中已知的TF-通路调控关系 ──────────────────────────
# 来源: TAIR, 文献挖掘
# 结构: tf_family → target_pathway → (score, description, at_genes)

ARABIDOPSIS_REGULON: dict[str, dict[str, dict]] = {
    "MYB": {
        "flavonoid": {
            "score": 0.90,
            "description": "AtMYB11/12/111 (R2R3-MYB) activate flavonol biosynthesis",
            "at_genes": ["AT3G62610", "AT2G47460", "AT5G49330"],
        },
        "anthocyanin": {
            "score": 0.95,
            "description": "AtMYB75/PAP1, AtMYB90/PAP2, AtMYB113/114 activate anthocyanin",
            "at_genes": ["AT1G56650", "AT1G66390", "AT1G66380"],
        },
        "lignin": {
            "score": 0.85,
            "description": "AtMYB58/63/85 specifically activate lignin biosynthesis",
            "at_genes": ["AT1G16490", "AT1G79180", "AT4G22680"],
        },
        "phenylpropanoid": {
            "score": 0.70,
            "description": "AtMYB4/7/32 regulate phenylpropanoid pathway",
            "at_genes": ["AT4G38620", "AT2G16720", "AT4G34990"],
        },
    },
    "bHLH": {
        "anthocyanin": {
            "score": 0.90,
            "description": "AtTT8, AtGL3, AtEGL3 form MBW complex for anthocyanin",
            "at_genes": ["AT4G09820", "AT5G41315", "AT1G63650"],
        },
        "flavonoid": {
            "score": 0.50,
            "description": "bHLH partners in MBW complex, indirect regulation",
            "at_genes": ["AT4G09820", "AT5G41315"],
        },
    },
    "WRKY": {
        "alkaloid": {
            "score": 0.40,
            "description": "WRKYs regulate defense-related alkaloids in Arabidopsis",
            "at_genes": ["AT2G23320", "AT4G23810"],
        },
        "flavonoid": {
            "score": 0.30,
            "description": "AtWRKY23/54 weakly affect flavonoid accumulation under stress",
            "at_genes": ["AT2G47260", "AT4G23810"],
        },
    },
    "AP2-ERF": {
        "terpenoid": {
            "score": 0.50,
            "description": "AtERF/ORA family regulates terpenoid-related defense genes",
            "at_genes": ["AT1G06160", "AT5G61600"],
        },
    },
    "NAC": {
        "lignin": {
            "score": 0.90,
            "description": "AtNST1/2/3, AtVND6/7 master regulators of secondary wall formation",
            "at_genes": ["AT2G46770", "AT3G61910", "AT1G32770"],
        },
    },
    "bZIP": {
        "flavonoid": {
            "score": 0.60,
            "description": "AtHY5 (bZIP) activates flavonoid genes under light",
            "at_genes": ["AT5G11260"],
        },
    },
}


class OrthologReasoner:
    """同源基因推理器。

    基于拟南芥中已表征的TF-通路调控关系，推断目标物种中
    同源TF的调控功能保守性。
    """

    @staticmethod
    def _resolve_family(tf_family: str) -> str | None:
        """尝试多种格式匹配TF家族名。"""
        candidates = [tf_family, tf_family.upper(), tf_family.lower(), tf_family.capitalize()]
        for c in candidates:
            if c in ARABIDOPSIS_REGULON:
                return c
        for c in candidates:
            alt = c.replace("/", "-")
            if alt in ARABIDOPSIS_REGULON:
                return alt
            alt = c.replace("-", "/")
            if alt in ARABIDOPSIS_REGULON:
                return alt
        return None

    @staticmethod
    def get_arabidopsis_evidence(tf_family: str, pathway: str) -> dict | None:
        """获取拟南芥中特定TF家族调控特定通路的证据。"""
        key = OrthologReasoner._resolve_family(tf_family)
        if key is None:
            return None
        family_data = ARABIDOPSIS_REGULON.get(key, {})
        # 精确匹配
        if pathway in family_data:
            return family_data[pathway]
        # 模糊匹配
        for pw, data in family_data.items():
            if pw in pathway or pathway in pw:
                return data
        return None

    @staticmethod
    def reason(
        tf_family: str | None,
        target_metabolite: str,
        has_ortholog: bool | None = None,
        ortholog_identity: float | None = None,
    ) -> ReasoningResult:
        """评估同源调控保守性。

        Args:
            tf_family: TF家族名称
            target_metabolite: 目标代谢物
            has_ortholog: 是否找到拟南芥同源基因, None=未知
            ortholog_identity: 同源基因的序列一致性 (0-100), None=未知

        Returns:
            ReasoningResult
        """
        if not tf_family:
            return ReasoningResult(
                score=0.0, label="none",
                explanation="未提供TF家族，无法进行同源推理。",
            )

        # 查询拟南芥已知调控关系
        pathway_key = target_metabolite.lower()
        evidence = OrthologReasoner.get_arabidopsis_evidence(tf_family, pathway_key)

        if evidence is None:
            key = OrthologReasoner._resolve_family(tf_family)
            family_data = ARABIDOPSIS_REGULON.get(key, {}) if key else {}
            if family_data:
                best_pathway, best_evidence = max(
                    family_data.items(), key=lambda x: x[1]["score"]
                )
                score = best_evidence["score"] * 0.5
                evidence_text = (
                    f"拟南芥{best_evidence['description']}. "
                    f"目标代谢物 '{target_metabolite}' 不一致, 但家族保守性部分支持."
                )
                return ReasoningResult(
                    score=score, label="low",
                    explanation=evidence_text,
                    details=[
                        f"目标物 '{target_metabolite}' ≠ 已知通路 '{best_pathway}'",
                        f"家族保守性: at_genes={', '.join(best_evidence['at_genes'])}",
                    ],
                )
            else:
                return ReasoningResult(
                    score=0.0, label="none",
                    explanation=f"TF家族 '{tf_family}' 在拟南芥中无已知调控关系。",
                )

        # 找到匹配
        base_score = evidence["score"]

        # 根据同源信息调整
        if has_ortholog is False:
            base_score *= 0.5
        elif has_ortholog is True and ortholog_identity is not None:
            if ortholog_identity >= 70:
                base_score *= 1.0
            elif ortholog_identity >= 50:
                base_score *= 0.8
            else:
                base_score *= 0.6

        details = [
            f"拟南芥同源基因: {', '.join(evidence['at_genes'])}",
            evidence["description"],
        ]
        if ortholog_identity is not None:
            details.append(f"序列一致性: {ortholog_identity:.1f}%")
        if has_ortholog is not None:
            details.append(f"同源基因已找到: {'是' if has_ortholog else '否'}")

        arab_evidence = PriorEvidence(
            strength="strong" if base_score >= 0.7 else "moderate" if base_score >= 0.4 else "weak",
            score=base_score,
            description=evidence["description"],
            pmids=[],  # 可扩展
        )

        return ReasoningResult(
            score=min(base_score, 1.0),
            explanation=evidence["description"],
            details=details,
            evidence_list=[arab_evidence],
        )
