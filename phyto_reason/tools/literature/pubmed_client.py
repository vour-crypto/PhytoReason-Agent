"""
PubMedLiterature — NCBI E-utilities 文献检索工具。

升级版：
  - 基础 search：自由 query 检索
  - smart_search：基于 TF-代谢物关联对自动构建精准查询
    "ZniTF1 (MYB) AND Berberine biosynthesis"
"""

from __future__ import annotations

import json
import os
import time
from typing import Any
from urllib.parse import quote

import pandas as pd
import requests

from phyto_reason.config.settings import Settings

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


class PubMedSearch:
    """PubMed 检索工具，支持智能 query 构建。"""

    def __init__(self, settings: Settings | None = None) -> None:
        self.settings = settings or Settings()
        self.cache_path = self.settings.cache_dir / ".pubmed_cache.json"
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self._cache: dict[str, list[dict]] = {}
        self._load_cache()

    # ── 缓存 ──────────────────────────────────────────────────

    def _load_cache(self) -> None:
        if self.cache_path.exists():
            try:
                self._cache = json.loads(self.cache_path.read_text(encoding="utf-8"))
            except Exception:
                self._cache = {}

    def _save_cache(self) -> None:
        self.cache_path.write_text(
            json.dumps(self._cache, ensure_ascii=False), encoding="utf-8"
        )

    # ── API 调用 ──────────────────────────────────────────────

    def _eutils_get(self, url: str, params: dict) -> dict:
        full_url = url + "?" + "&".join(
            f"{k}={quote(str(v))}" for k, v in params.items()
        )
        for attempt in range(self.settings.pubmed_max_retries):
            try:
                resp = requests.get(
                    full_url,
                    timeout=self.settings.pubmed_timeout,
                    headers={"User-Agent": self.settings.ncbi_tool},
                )
                resp.raise_for_status()
                return resp.json()
            except Exception:
                if attempt < self.settings.pubmed_max_retries - 1:
                    time.sleep(1.5 ** attempt)
        return {}

    def search(self, query: str, max_results: int = 5) -> list[dict]:
        """通用 PubMed 检索（自由 query）。"""
        cache_key = f"pm::{query}::{max_results}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        time.sleep(self.settings.pubmed_request_interval)
        id_data = self._eutils_get(ESEARCH, {
            "db": "pubmed", "term": query, "retmax": max_results,
            "retmode": "json", "email": self.settings.ncbi_email,
            "tool": self.settings.ncbi_tool,
        })
        ids = id_data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            self._cache[cache_key] = []
            return []

        summary_data = self._eutils_get(ESUMMARY, {
            "db": "pubmed", "id": ",".join(ids), "retmode": "json",
            "email": self.settings.ncbi_email,
            "tool": self.settings.ncbi_tool,
        })
        results = []
        for pmid in ids:
            entry = summary_data.get("result", {}).get(pmid, {})
            if not entry:
                continue
            doi = next(
                (a["value"] for a in entry.get("articleids", [])
                 if a.get("idtype") == "doi"), "",
            )
            results.append({
                "PMID": pmid,
                "Title": entry.get("title", ""),
                "Source": entry.get("source", ""),
                "PubDate": entry.get("pubdate", ""),
                "DOI": doi,
            })

        self._cache[cache_key] = results
        self._save_cache()
        return results

    # ── 智能 query 构建 ────────────────────────────────────

    def smart_search(
        self,
        tf_gene_id: str = "",
        metabolite: str = "",
        tf_family: str = "",
        metabolite_class: str = "",
        species: str = "",
        max_results: int = 5,
    ) -> list[dict]:
        """
        基于关联对自动构建精准 PubMed 查询。

        query 格式:
          "TF_ID (TF_FAMILY) AND metabolite biosynthesis plant"
        例:
          "Zni22G003130 (WRKY) AND p-coumaroylagmatine biosynthesis"
        """
        parts = []
        if tf_gene_id:
            # Expand gene ID with synonyms for better PubMed recall
            try:
                from phyto_reason.knowledge.gene_resolver import expand_gene_query
                expanded = expand_gene_query(tf_gene_id, species)
                parts.append(expanded if expanded != tf_gene_id else tf_gene_id)
            except Exception:
                parts.append(tf_gene_id)
        if tf_family:
            parts.append(f"({tf_family})")
        if metabolite:
            parts.append(f"\"{metabolite}\" biosynthesis")
        if metabolite_class:
            parts.append(metabolite_class)
        if species:
            parts.append(species)
        else:
            sp = self.settings.species_rules.get("species", {}).get("name", "")
            if sp:
                parts.append(sp)

        query = " AND ".join(parts) if parts else ""
        if not query:
            return []

        return self.search(query, max_results=max_results)

    # ── 便捷快捷方式 ─────────────────────────────────────────

    def search_tf_family(self, family: str) -> list[dict]:
        sp = self.settings.species_rules.get("species", {}).get("name", "plant")
        queries = [
            f"{family} transcription factor {sp}",
            f"{family} transcription factor plant secondary metabolism",
        ]
        seen: set[str] = set()
        results: list[dict] = []
        for q in queries:
            for r in self.search(q):
                pmid = r.get("PMID", "")
                if pmid and pmid not in seen:
                    seen.add(pmid)
                    results.append(r)
        return results

    # ── 格式化 ──────────────────────────────────────────────

    def format_results(self, results: list[dict], max_items: int = 5) -> str:
        if not results:
            return "（PubMed 未找到相关文献）"
        lines = [f"## PubMed 检索结果（共 {len(results)} 篇）\n"]
        for r in results[:max_items]:
            lines.append(f"- **{r.get('Title', '?')}**")
            lines.append(f"  {r.get('Source', '')} ({r.get('PubDate', '')})")
            lines.append(f"  PMID: {r.get('PMID', '')} | DOI: {r.get('DOI', '')}")
            lines.append("")
        return "\n".join(lines)
