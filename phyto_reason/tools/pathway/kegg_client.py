"""
kegg_client.py — KEGG REST API 真实数据集成。

支持:
  - 查询化合物-通路映射
  - 查询通路-酶映射
  - 查询通路层次结构
  - 查询相关代谢物
  - 搜索化合物

使用 KEGG REST API (https://rest.kegg.jp/)。
无需本地数据库。
"""

from __future__ import annotations

import json
import logging
import time
import os
from urllib.parse import quote
from pathlib import Path
from typing import Any

import requests

from phyto_reason.reasoning.reasoning_result import PriorEvidence, ReasoningResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.candidate_gene import CandidateGene

logger = logging.getLogger("kegg_client")

KEGG_BASE = "https://rest.kegg.jp"
REQUEST_INTERVAL = 0.3


class KEGGClient:
    """KEGG REST API 客户端。"""

    def __init__(self, cache_dir: str | None = None) -> None:
        self._cache: dict[str, Any] = {}
        self._cache_dir = Path(cache_dir) if cache_dir else None

    # ── 核心 API 调用 ──────────────────────────────────────

    def _get(self, endpoint: str) -> str | None:
        cache_key = endpoint
        if cache_key in self._cache:
            return self._cache[cache_key]

        url = f"{KEGG_BASE}/{endpoint}"
        try:
            time.sleep(REQUEST_INTERVAL)
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            text = resp.text
            self._cache[cache_key] = text
            return text
        except Exception as e:
            logger.warning(f"KEGG API error: {url}: {e}")
            return None

    def _get_json(self, endpoint: str) -> dict | None:
        text = self._get(endpoint)
        if text is None:
            return None
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return None

    # ── 通路查询 ───────────────────────────────────────────

    def find_compound(self, keyword: str, max_results: int = 10) -> list[dict]:
        """搜索 KEGG 化合物。"""
        text = self._get(f"find/compound/{keyword}")
        if not text:
            return []
        results = []
        for line in text.strip().split("\n")[:max_results]:
            parts = line.split("\t")
            if len(parts) >= 2:
                results.append({"id": parts[0].strip(), "name": parts[1].strip()})
        return results

    def get_compound_pathways(self, compound_id: str) -> list[dict]:
        """查询化合物参与的通路。"""
        text = self._get(f"link/pathway/{compound_id}")
        if not text:
            return []
        pathways = []
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pw_id = parts[0].strip()
                pw_name = parts[1].strip()
                pathways.append({"id": pw_id, "name": pw_name})
        return pathways

    def get_pathway_enzymes(self, pathway_id: str) -> list[dict]:
        """查询通路中的酶。"""
        text = self._get(f"link/ec/{pathway_id}")
        if not text:
            return []
        enzymes = []
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                enzymes.append({"ec": parts[0].strip(), "name": parts[1].strip()})
        return enzymes

    def get_pathway_genes(self, pathway_id: str, organism: str = "ath") -> list[dict]:
        """查询通路中的基因 (指定物种)。"""
        endpoint = f"link/genes/{pathway_id}" if not organism else f"link/{organism}/{pathway_id}"
        text = self._get(endpoint)
        if not text:
            return []
        genes = []
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                genes.append({"gene_id": parts[0].strip(), "name": parts[1].strip()})
        return genes

    def get_compound_detail(self, compound_id: str) -> dict:
        """获取化合物详细信息。"""
        text = self._get(f"get/{compound_id}")
        if not text:
            return {}
        detail: dict[str, Any] = {"id": compound_id}
        current_key = None
        for line in text.split("\n"):
            if line.startswith(" ") and current_key:
                detail.setdefault(current_key, "").append(line.strip())
            elif line.strip() and not line.startswith(" "):
                if "  " in line:
                    key, value = line.split("  ", 1)
                    current_key = key.strip()
                    detail[current_key] = [value.strip()]
        for k in list(detail.keys()):
            if isinstance(detail.get(k), list):
                detail[k] = "\n".join(detail[k])
        return detail

    def search_compound_by_name(self, name: str) -> list[dict]:
        """按名称搜索化合物，模糊匹配。"""
        results = self.find_compound(name)
        if not results:
            results = self.find_compound(name.lower())
        return results

    def get_pathway_hierarchy(self, pathway_id: str) -> list[dict]:
        """获取通路的层次分类。"""
        detail = {}
        text = self._get(f"get/{pathway_id}")
        if not text:
            return []
        hierarchy = []
        for line in text.split("\n"):
            if "Module" in line:
                hierarchy.append({"type": "module", "info": line.strip()})
        return hierarchy


@register_tool
class KEGGTool(BaseTool):
    """KEGG 通路分析工具。"""
    tool_name = "kegg_pathway"
    description = "查询 KEGG 化合物-通路-酶映射关系。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="compound_name", type="string",
                      description="化合物名称 (如 berberine, flavonoid)", required=True),
        ToolParameter(name="analysis_type", type="string",
                      description="分析类型: pathway|enzymes|detail|all", default="all"),
    ]
    supported_data_types = ["metabolite", "pathway"]

    def __init__(self) -> None:
        self.client = KEGGClient()

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        if not kwargs.get("compound_name"):
            errors.append("compound_name 不能为空")
        return errors

    def run(self, **kwargs) -> ToolResult:
        compound_name: str = kwargs.get("compound_name", "")
        analysis_type: str = kwargs.get("analysis_type", "all")

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        compounds = self.client.search_compound_by_name(compound_name)
        if not compounds:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=[f"KEGG 中未找到化合物: '{compound_name}'"],
                metadata={"compound": compound_name, "found": False},
            )

        all_pathways: list[dict] = []
        all_enzymes: list[dict] = []

        for cpd in compounds[:3]:
            cpd_id = cpd["id"]
            if analysis_type in ("pathway", "all"):
                pathways = self.client.get_compound_pathways(cpd_id)
                all_pathways.extend(pathways)
            if analysis_type in ("enzymes", "all"):
                for pw in self.client.get_compound_pathways(cpd_id)[:5]:
                    enzymes = self.client.get_pathway_enzymes(pw["id"])
                    all_enzymes.extend(enzymes)

        fallback_used = False
        if not all_pathways and compounds:
            try:
                direct = requests.get(
                    f"{KEGG_BASE}/find/pathway/{quote(compound_name)}",
                    timeout=15,
                )
                direct.raise_for_status()
                for line in direct.text.splitlines():
                    if "\t" in line:
                        pathway_id, description = line.split("\t", 1)
                        all_pathways.append({"pathway_id": pathway_id, "name": description.strip()})
                fallback_used = bool(all_pathways)
            except Exception as exc:
                logger.debug("KEGG pathway-name fallback failed: %s", exc)
        if fallback_used:
            warnings.append(f"化合物无直接通路链接，已返回名称匹配的通路 {len(all_pathways)} 条")
        elif not all_pathways:
            warnings.append(
                f"KEGG 找到 {len(compounds)} 个化合物，但未找到关联通路；"
                "可尝试更具体的化合物名称。"
            )

        evidence_list.append(Evidence(
            evidence_type=EvidenceType.PATHWAY_CONSISTENCY,
            source=f"kegg_tool_{self.version}",
            score=min(len(all_pathways) / 5, 1.0),
            description=f"KEGG compound={compound_name}: {len(all_pathways)} pathways, {len(all_enzymes)} enzymes",
            metadata={
                "compound": compound_name,
                "pathways": all_pathways[:10],
                "enzymes": all_enzymes[:20],
            },
        ))

        return ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "compound_searched": compound_name,
                "compounds_found": len(compounds),
                "pathways_found": len(all_pathways),
                "count": len(all_pathways),
                "enzymes_found": len(all_enzymes),
                "pathway_names": [p["name"] for p in all_pathways[:5]],
            },
        )
