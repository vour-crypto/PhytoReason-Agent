"""
pubmed_tool.py — PubMed 文献检索工具 (BaseTool 封装)。
"""

from __future__ import annotations

import logging

from phyto_reason.models.tool_result import ToolResult
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.tools.literature.pubmed_client import PubMedSearch

logger = logging.getLogger("pubmed_tool")


@register_tool
class PubMedTool(BaseTool):
    tool_name = "literature_search"
    description = "PubMed 文献检索，基于 TF 家族、代谢物名自动构建查询。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="tf_family", type="string",
                      description="TF 家族名称 (如 MYB, WRKY)", default=""),
        ToolParameter(name="metabolite", type="string",
                      description="目标代谢物名称", default=""),
        ToolParameter(name="gene_id", type="string",
                      description="候选基因 ID", default=""),
        ToolParameter(name="species", type="string",
                      description="物种名称", default=""),
        ToolParameter(name="max_results", type="integer",
                      description="最多返回几篇文献", default=5),
    ]
    supported_data_types = ["literature"]

    def __init__(self) -> None:
        self.client = PubMedSearch()

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        if not kwargs.get("tf_family") and not kwargs.get("metabolite") and not kwargs.get("gene_id"):
            errors.append("至少需要提供 tf_family、metabolite 或 gene_id 之一")
        return errors

    def run(self, **kwargs) -> ToolResult:
        tf_family: str = kwargs.get("tf_family", "") or ""
        metabolite: str = kwargs.get("metabolite", "") or ""
        gene_id: str = kwargs.get("gene_id", "") or ""
        species: str = kwargs.get("species", "") or ""
        max_results: int = kwargs.get("max_results", 5)

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        results = self.client.smart_search(
            tf_gene_id=gene_id, metabolite=metabolite,
            tf_family=tf_family, species=species,
            max_results=max_results,
        )

        if not results and (tf_family or metabolite):
            query = " ".join(filter(None, [tf_family, metabolite, "biosynthesis", "plant"]))
            if query:
                results = self.client.search(query, max_results=max_results)

        if results:
            n_pmids = len(results)
            score = min(n_pmids / 3, 1.0) * 0.8

            pmids_list = [r.get("PMID", "") for r in results]

            if gene_id:
                candidates.append(CandidateGene(
                    gene_id=gene_id, is_tf=bool(tf_family), tf_family=tf_family or None,
                    literature_score=round(score, 4),
                ))

            evidence_list.append(Evidence(
                evidence_type=EvidenceType.LITERATURE,
                source=f"pubmed_tool_{self.version}",
                score=round(score, 4),
                confidence=round(score, 4),
                description=f"PubMed: {n_pmids} articles for {tf_family or ''} {metabolite or ''}",
                metadata={
                    "tf_family": tf_family,
                    "metabolite": metabolite,
                    "n_articles": n_pmids,
                    "pmids": pmids_list,
                    "titles": [r.get("Title", "")[:100] for r in results[:5]],
                },
            ))
        else:
            warnings.append(f"PubMed 未找到相关文献")

        return ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "tf_family": tf_family,
                "metabolite": metabolite,
                "n_articles": len(results),
                "count": len(results),
                "records": [{
                    "title": r.get("Title", r.get("title", "")),
                    "abstract": r.get("Abstract", r.get("abstract", "")),
                    "source": r.get("Source", r.get("source", "PubMed")),
                } for r in results[:10]],
                "query_type": "smart_search",
            },
        )
