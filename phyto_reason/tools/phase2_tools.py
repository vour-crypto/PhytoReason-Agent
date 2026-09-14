"""Phase 2 public literature and plant-resource tools (v2 per API)."""
from __future__ import annotations

import hashlib
import os
import sqlite3
import threading
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import requests
from tenacity import retry, stop_after_attempt, wait_exponential

from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.utils.cache import knowledge_cache

_LAST_REQUEST = 0.0
_RATE_LOCK = threading.Lock()
_ATOM = "{http://www.w3.org/2005/Atom}"
_PLANTTF_DB = Path(__file__).resolve().parents[1] / "data" / "planttf.db"


def _throttle(interval: float = 0.2) -> None:
    global _LAST_REQUEST
    with _RATE_LOCK:
        delay = interval - (time.monotonic() - _LAST_REQUEST)
        if delay > 0:
            time.sleep(delay)
        _LAST_REQUEST = time.monotonic()


def _headers() -> dict[str, str]:
    return {"User-Agent": "PhytoReason/5.1"}


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=0.4, max=3), reraise=True)
def _get_json(url: str, params: dict[str, Any], headers: dict[str, str] | None = None) -> Any:
    key = "p2v2:" + hashlib.md5((url + repr(sorted(params.items()))).encode()).hexdigest()
    cached = knowledge_cache.get(key)
    if cached is not None:
        return cached
    _throttle()
    response = requests.get(url, params=params, timeout=15, headers=headers or _headers())
    response.raise_for_status()
    value = response.json()
    knowledge_cache.set(key, value, ttl=86400)
    return value


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=0.4, max=3), reraise=True)
def _get_xml(url: str, params: dict[str, Any]) -> ET.Element:
    key = "p2v2xml:" + hashlib.md5((url + repr(sorted(params.items()))).encode()).hexdigest()
    cached = knowledge_cache.get(key)
    if cached is not None:
        return ET.fromstring(cached)
    _throttle()
    response = requests.get(url, params=params, timeout=15, headers=_headers())
    response.raise_for_status()
    knowledge_cache.set(key, response.text, ttl=86400)
    return ET.fromstring(response.text)


def _result(source: str, query: str, records: list[dict[str, Any]], error: str = "") -> ToolResult:
    score = min(len(records) / 5, 1.0)
    evidence = Evidence(evidence_type=EvidenceType.LITERATURE, source=source, score=score, confidence=score,
                        description=f"{source}: {len(records)} records for {query}", metadata={"query": query, "records": records[:10]})
    return ToolResult(status="partial" if error else "success", evidence_list=[evidence], warnings=[error] if error else [],
                      metadata={"source": source, "query": query, "count": len(records)})


class _SearchTool(BaseTool):
    parameters = [ToolParameter("query", "string", "Literature or resource query", required=True), ToolParameter("max_results", "integer", "Maximum records", default=10)]
    supported_data_types = ["literature"]
    endpoint = ""
    source = ""

    def validate_input(self, **kwargs) -> list[str]:
        return [] if str(kwargs.get("query", "")).strip() else ["query 不能为空"]

    def build_params(self, query: str, limit: int) -> dict[str, Any]:
        raise NotImplementedError

    def parse_records(self, data: Any) -> list[dict[str, Any]]:
        raise NotImplementedError

    def run(self, **kwargs) -> ToolResult:
        query = str(kwargs.get("query", "")).strip()
        if not query:
            return _result(self.source, "", [], "query 不能为空")
        limit = min(max(1, int(kwargs.get("max_results", 10))), 20)
        try:
            data = _get_json(self.endpoint, self.build_params(query, limit))
            return _result(self.source, query, self.parse_records(data))
        except Exception as exc:
            return _result(self.source, query, [], f"{self.source} 请求失败: {exc}")


@register_tool
class EuropePMCSearchTool(_SearchTool):
    tool_name, description, version, source = "europe_pmc_search", "Europe PMC 公开文献检索（含摘要）", "1.1.0", "Europe PMC"
    endpoint = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    def build_params(self, query, limit): return {"query": query, "format": "json", "pageSize": limit, "resultType": "core"}
    def parse_records(self, data): return list(data.get("resultList", {}).get("result", []))


@register_tool
class OpenAlexSearchTool(_SearchTool):
    tool_name, description, version, source = "openalex_search", "OpenAlex 学术文献检索", "1.1.0", "OpenAlex"
    endpoint = "https://api.openalex.org/works"
    def build_params(self, query, limit): return {"search": query, "per_page": limit}
    def parse_records(self, data): return list(data.get("results", []))


@register_tool
class CrossrefSearchTool(_SearchTool):
    tool_name, description, version, source = "crossref_search", "Crossref DOI 元数据检索", "1.1.0", "Crossref"
    endpoint = "https://api.crossref.org/works"
    def build_params(self, query, limit): return {"query": query, "rows": limit}
    def parse_records(self, data): return list(data.get("message", {}).get("items", []))


@register_tool
class BioRxivSearchTool(_SearchTool):
    tool_name, description, version, source = "biorxiv_search", "bioRxiv 预印本检索（经 Europe PMC SRC:PPR 索引）", "1.1.0", "bioRxiv（经 Europe PMC 索引）"
    endpoint = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    def build_params(self, query, limit): return {"query": f"({query}) AND (SRC:PPR)", "format": "json", "pageSize": limit, "resultType": "core"}
    def parse_records(self, data): return list(data.get("resultList", {}).get("result", []))


@register_tool
class ArxivSearchTool(_SearchTool):
    tool_name, description, version, source = "arxiv_search", "arXiv 预印本检索（Atom XML 解析）", "1.1.0", "arXiv"
    endpoint = "http://export.arxiv.org/api/query"
    def run(self, **kwargs) -> ToolResult:
        query = str(kwargs.get("query", "")).strip()
        if not query: return _result(self.source, "", [], "query 不能为空")
        try:
            root = _get_xml(self.endpoint, {"search_query": f"all:{query}", "max_results": min(max(1, int(kwargs.get("max_results", 10))), 20)})
            records = []
            for entry in root.findall(f"{_ATOM}entry"):
                def text(name):
                    node = entry.find(f"{_ATOM}{name}"); return (node.text or "").strip() if node is not None else ""
                records.append({"title": text("title"), "summary": text("summary")[:300], "published": text("published"), "id": text("id")})
            return _result(self.source, query, records)
        except Exception as exc: return _result(self.source, query, [], f"arXiv 请求失败: {exc}")


@register_tool
class SemanticScholarSearchTool(_SearchTool):
    tool_name, description, version, source = "semantic_scholar_search", "Semantic Scholar 文献检索", "1.1.0", "Semantic Scholar"
    endpoint = "https://api.semanticscholar.org/graph/v1/paper/search"
    def request_headers(self):
        key = os.getenv("SEMANTIC_SCHOLAR_API_KEY", "").strip()
        return {**_headers(), **({"x-api-key": key} if key else {})}
    def build_params(self, query, limit): return {"query": query, "limit": limit, "fields": "title,abstract,year,venue,externalIds"}
    def parse_records(self, data): return list(data.get("data", []))
    def run(self, **kwargs):
        query = str(kwargs.get("query", "")).strip()
        if not query: return _result(self.source, "", [], "query 不能为空")
        limit = min(max(1, int(kwargs.get("max_results", 10))), 20)
        try:
            return _result(self.source, query, self.parse_records(_get_json(self.endpoint, self.build_params(query, limit), self.request_headers())))
        except Exception as exc:
            # Anonymous Semantic Scholar quotas can return 429. OpenAlex is
            # a compatible public scholarly index, so preserve availability
            # while making the fallback explicit in warnings and source data.
            try:
                fallback = _get_json("https://api.openalex.org/works", {"search": query, "per_page": limit})
                records = list(fallback.get("results", [])) if isinstance(fallback, dict) else []
                result = _result(self.source, query, records)
                result.status = "partial"
                if result.evidence_list:
                    result.evidence_list[0].source = "OpenAlex（Semantic Scholar 降级）"
                    result.evidence_list[0].metadata["actual_source"] = "OpenAlex"
                result.warnings.append(f"Semantic Scholar 请求失败（{exc}），已降级为 OpenAlex 记录")
                result.metadata["fallback_source"] = "OpenAlex"
                return result
            except Exception:
                return _result(self.source, query, [], f"{self.source} 请求失败: {exc}")


class _PlantTool(_SearchTool):
    supported_data_types = ["annotation", "pathway"]


@register_tool
class PlantTFDBLookupTool(_PlantTool):
    tool_name, description, version, source = "planttfdb_lookup", "PlantTFDB 转录因子查询（本地批量数据；未构建时诚实降级）", "1.1.0", "PlantTFDB（本地批量数据）"
    def run(self, **kwargs):
        query = str(kwargs.get("query", "")).strip()
        if not query: return _result(self.source, "", [], "query 不能为空")
        if not _PLANTTF_DB.exists(): return _result(self.source, query, [], "PlantTFDB 无公开 REST API 且本地索引未构建")
        with sqlite3.connect(_PLANTTF_DB) as conn:
            conn.row_factory = sqlite3.Row; like = f"%{query}%"; rows = conn.execute("SELECT species,gene_id,family,description FROM tf WHERE gene_id LIKE ? OR family LIKE ? OR description LIKE ? LIMIT 20", (like, like, like)).fetchall()
        return _result(self.source, query, [dict(row) for row in rows])


@register_tool
class JASPARMotifLookupTool(_PlantTool):
    tool_name, description, version, source = "jaspar_motif_lookup", "JASPAR 植物转录因子 motif 查询", "1.1.0", "JASPAR"
    endpoint = "https://jaspar.genereg.net/api/v1/matrix/"
    def build_params(self, query, limit): return {"search": query, "tax_group": "plants", "page_size": limit}
    def parse_records(self, data):
        rows = data.get("results", []) if isinstance(data, dict) else list(data)
        fields = ("matrix_id", "name", "tax_group", "species", "collection", "class")
        records = []
        for row in rows:
            detail = {}
            detail_url = row.get("url")
            if detail_url:
                try:
                    detail = _get_json(detail_url, {})
                except Exception:
                    detail = {}
            records.append({key: detail.get(key, row.get(key)) for key in fields})
        return records


@register_tool
class EnsemblPlantsLookupTool(_PlantTool):
    tool_name, description, version, source = "ensembl_plants_lookup", "Ensembl Plants 同源基因查询（query=stable_id；species 可选）", "1.2.0", "Ensembl Plants"
    base_url = "https://rest.ensembl.org/homology/id/"
    parameters = [ToolParameter("query", "string", "Ensembl stable_id", required=True), ToolParameter("species", "string", "Ensembl species name", default="arabidopsis_thaliana")]
    def run(self, **kwargs):
        stable_id = str(kwargs.get("query", "")).strip()
        if not stable_id: return _result(self.source, "", [], "query 应为 Ensembl stable_id")
        species = str(kwargs.get("species", "")).strip() or "arabidopsis_thaliana"
        try:
            data = _get_json(f"{self.base_url}{species}/{stable_id}", {"compara": "plants", "type": "orthologues"}, {**_headers(), "Content-Type": "application/json", "Accept": "application/json"})
            homologies = data.get("data", [{}])[0].get("homologies", [])
            records = [{"type": h.get("type"), "source_id": h.get("source_id"), "target_id": h.get("target_id"), "target_species": h.get("target_species")} for h in homologies]
            return _result(self.source, f"{species}/{stable_id}", records)
        except Exception as exc: return _result(self.source, f"{species}/{stable_id}", [], f"Ensembl Plants 请求失败: {exc}")
