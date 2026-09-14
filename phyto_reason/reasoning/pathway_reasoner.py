"""
pathway_reasoner.py — 代谢通路一致性推理。

评估候选基因在目标代谢通路中的调控逻辑一致性:
  - 上游酶共表达一致性 (upstream enzyme coherence)
  - 通路分支竞争 (branch competition)
  - 限速酶调控 (bottleneck inference)
"""

from __future__ import annotations

from phyto_reason.reasoning.reasoning_result import PriorEvidence, ReasoningResult

# ── 通路知识库 ──────────────────────────────────────────────────
# 定义关键代谢通路的结构: 酶名 → 是否是限速酶

PATHWAY_KNOWLEDGE: dict[str, dict] = {
    "flavonoid": {
        "name": "Flavonoid biosynthesis",
        "enzymes": [
            "PAL", "C4H", "4CL", "CHS", "CHI",
            "F3H", "F3'H", "FLS", "DFR", "ANS", "UFGT",
        ],
        "bottleneck_enzymes": ["CHS", "DFR", "ANS"],
        "upstream_enzymes": ["PAL", "C4H", "4CL"],
        "branch_points": {
            "F3H": "Flavonol vs anthocyanin branch",
            "DFR": "Anthocyanin branch entry",
        },
    },
    "alkaloid": {
        "name": "Alkaloid biosynthesis",
        "enzymes": [
            "TYDC", "NCS", "STS", "CYP80B", "BBE",
            "SDR", "OMT", "NMCH", "CYP450",
        ],
        "bottleneck_enzymes": ["NCS", "TYDC", "CYP80B"],
        "upstream_enzymes": ["TYDC", "NCS"],
        "branch_points": {
            "NCS": "Benzylisoquinoline branch entry",
            "CYP80B": "Berberine vs morphinan branch",
        },
    },
    "terpenoid": {
        "name": "Terpenoid biosynthesis",
        "enzymes": [
            "HMGR", "HMGS", "MK", "PMK", "MVD",
            "DXS", "DXR", "MCT", "CMK", "MDS", "HDS", "HDR",
            "GPPS", "FPPS", "GGPPS", "TPS",
        ],
        "bottleneck_enzymes": ["HMGR", "DXS", "TPS"],
        "upstream_enzymes": ["HMGR", "DXS", "DXR"],
        "branch_points": {
            "FPPS": "Sesquiterpene vs triterpene branch",
            "GGPPS": "Diterpene branch entry",
            "TPS": "Terpene scaffold diversification",
        },
    },
    "phenylpropanoid": {
        "name": "Phenylpropanoid biosynthesis",
        "enzymes": [
            "PAL", "C4H", "4CL", "HCT", "C3'H",
            "CCoAOMT", "CCR", "CAD",
        ],
        "bottleneck_enzymes": ["PAL", "C4H"],
        "upstream_enzymes": ["PAL", "C4H", "4CL"],
        "branch_points": {
            "4CL": "Flavonoid vs lignin vs coumarin branch",
        },
    },
    "lignin": {
        "name": "Lignin biosynthesis",
        "enzymes": [
            "PAL", "C4H", "4CL", "HCT", "C3'H",
            "CCoAOMT", "CCR", "CAD", "F5H", "COMT", "LAC", "PRX",
        ],
        "bottleneck_enzymes": ["CCR", "CAD"],
        "upstream_enzymes": ["PAL", "C4H", "4CL"],
        "branch_points": {
            "HCT": "G vs S lignin branch",
            "F5H": "S lignin branch entry",
        },
    },
}


def _detect_pathway_class(metabolite: str) -> str | None:
    """根据代谢物名称或类别推测其所属通路。"""
    m = metabolite.lower()
    for pathway, info in PATHWAY_KNOWLEDGE.items():
        if pathway in m:
            return pathway
    # fuzzy matching
    keyword_map: dict[str, str] = {
        "flavonoid": "flavonoid", "anthocyanin": "flavonoid",
        "flavone": "flavonoid", "flavonol": "flavonoid",
        "alkaloid": "alkaloid", "isoquinoline": "alkaloid",
        "terpene": "terpenoid", "terpenoid": "terpenoid",
        "lignin": "lignin", "phenylpropanoid": "phenylpropanoid",
        "phenolic": "phenylpropanoid",
    }
    for keyword, path in keyword_map.items():
        if keyword in m:
            return path
    return None


class PathwayReasoner:
    """代谢通路一致性推理器。

    评估候选基因在代谢通路中的调控位置、上下游一致性和限速步骤。
    """

    @staticmethod
    def get_pathway_info(pathway: str) -> dict | None:
        """获取指定通路的结构信息。"""
        return PATHWAY_KNOWLEDGE.get(pathway)

    @staticmethod
    def reason(
        tf_family: str | None,
        target_metabolite: str,
        upstream_genes_found: list[str] | None = None,
        coexpression_consistency: float | None = None,
        bottleneck_hit: bool | None = None,
    ) -> ReasoningResult:
        """评估候选TF对目标代谢通路的调控一致性。

        Args:
            tf_family: TF家族
            target_metabolite: 目标代谢物
            upstream_genes_found: 在DEG/共表达中发现的通路上游酶基因列表
            coexpression_consistency: 上游酶与TF的共表达一致性 (0-1), None=未知
            bottleneck_hit: TF是否调控限速酶基因, None=未知

        Returns:
            ReasoningResult
        """
        pathway = _detect_pathway_class(target_metabolite)
        if not pathway:
            return ReasoningResult(
                score=0.30, label="low",
                explanation=f"无法将 '{target_metabolite}' 映射到已知代谢通路。"
                            "使用默认通路一致性分数。",
            )

        info = PATHWAY_KNOWLEDGE[pathway]
        details = [
            f"目标代谢物映射到通路: {info['name']}",
            f"通路包含 {len(info['enzymes'])} 个已知酶基因",
        ]

        upstream_score = 0.0
        bottleneck_score = 0.0

        # 上游酶一致性评估
        if upstream_genes_found is not None:
            upstream = info["upstream_enzymes"]
            found_upstream = [g for g in upstream_genes_found if g in upstream]
            if upstream:
                ratio = len(found_upstream) / len(upstream)
                upstream_score = min(ratio, 1.0)
                details.append(
                    f"上游酶检测: 在DEG中发现 {len(found_upstream)}/{len(upstream)} "
                    f"个上游酶 ({found_upstream})"
                )
            else:
                upstream_score = 0.5  # 无为中性
        else:
            upstream_score = 0.5  # 无数据时取中性

        # 共表达一致性
        if coexpression_consistency is not None:
            upstream_score = (upstream_score + coexpression_consistency) / 2
            details.append(f"共表达一致性: {coexpression_consistency:.2f}")

        # 限速酶调控评估
        if bottleneck_hit is not None:
            bottleneck_score = 0.8 if bottleneck_hit else 0.3
            label = "调控" if bottleneck_hit else "未调控"
            details.append(f"限速酶: {label}通路限速酶 "
                           f"({', '.join(info['bottleneck_enzymes'])})")
        else:
            bottleneck_score = 0.5  # 无数据时取中性

        # 通路分支竞争
        branch_info = ""
        if info.get("branch_points"):
            branches = list(info["branch_points"].values())
            branch_info = f"通路分支点: {'; '.join(branches)}"
            details.append(branch_info)

        final_score = 0.4 * upstream_score + 0.3 * bottleneck_score + 0.3 * 0.5

        explanation_parts = [
            f"通路 '{info['name']}' 一致性分析:",
            f"  - 上游酶一致性: {upstream_score:.2f}",
            f"  - 限速酶调控: {bottleneck_score:.2f}",
            f"  - 综合通路分数: {final_score:.2f}",
        ]
        if branch_info:
            explanation_parts.append(f"  - {branch_info}")

        evidence = PriorEvidence(
            strength="moderate" if final_score >= 0.5 else "weak",
            score=final_score,
            description=info["name"],
            pmids=[],
        )

        return ReasoningResult(
            score=final_score,
            explanation="\n".join(explanation_parts),
            details=details,
            evidence_list=[evidence],
        )

    @staticmethod
    def _resolve_via_ontology(metabolite: str) -> str | None:
        """通过 ontology 解析代谢物的通路。"""
        try:
            from phyto_reason.ontology.ontology_resolver import OntologyResolver
            resolver = OntologyResolver()
            resolved = resolver.resolve(metabolite)
            if resolved.node_type == "metabolite" and resolved.pathway_ids:
                for parent in resolved.parent_nodes:
                    pathway_key = _detect_pathway_class(parent)
                    if pathway_key:
                        return pathway_key
                return _detect_pathway_class(resolved.resolved_node)
        except ImportError:
            pass
        return None

    @staticmethod
    def reason_ontology(
        tf_family: str | None,
        target_metabolite: str,
        upstream_genes_found: list[str] | None = None,
        coexpression_consistency: float | None = None,
        bottleneck_hit: bool | None = None,
    ) -> ReasoningResult:
        """Ontology-aware 通路一致性推理。

        先尝试 ontology 解析，失败时回退到原始推理。
        """
        # Try direct mapping first
        direct_pathway = _detect_pathway_class(target_metabolite)
        if not direct_pathway:
            ontology_pathway = PathwayReasoner._resolve_via_ontology(target_metabolite)
            if ontology_pathway:
                return PathwayReasoner.reason(
                    tf_family, ontology_pathway,
                    upstream_genes_found, coexpression_consistency, bottleneck_hit,
                )
        return PathwayReasoner.reason(
            tf_family, target_metabolite,
            upstream_genes_found, coexpression_consistency, bottleneck_hit,
        )

    @staticmethod
    def reason_bottleneck(
        tf_family: str | None,
        target_metabolite: str,
        target_enzymes: list[str],
    ) -> ReasoningResult:
        """评估TF是否调控限速酶。"""
        pathway = _detect_pathway_class(target_metabolite)
        if not pathway:
            return ReasoningResult(
                score=0.0, label="none",
                explanation="无法识别通路, 无法评估限速酶调控。",
            )

        info = PATHWAY_KNOWLEDGE[pathway]
        bottlenecks = info["bottleneck_enzymes"]
        hits = [e for e in target_enzymes if e.upper() in [b.upper() for b in bottlenecks]]

        if hits:
            return ReasoningResult(
                score=0.80, label="high",
                explanation=f"TF可能调控限速酶: {', '.join(hits)}",
                details=[f"在 {info['name']} 通路中调控限速酶"],
            )
        else:
            return ReasoningResult(
                score=0.30, label="low",
                explanation=f"未发现TF调控通路限速酶 ({', '.join(bottlenecks)})",
                details=["未命中限速酶"],
            )
