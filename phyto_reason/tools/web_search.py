"""
web_search.py — Broad search tool using NCBI E-utilities across multiple databases.

Queries PubMed, PMC (full-text), Gene, and Nucleotide databases simultaneously,
providing a much broader search surface than PubMed alone. Falls back to
MyGene.info for gene-specific queries.

All databases are free and already accessible in this network environment.
"""

from __future__ import annotations

import concurrent.futures
import logging
import re
import time
from typing import Any

import requests

from phyto_reason.utils.cache import knowledge_cache

logger = logging.getLogger("web_search")

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "phyto_reason"
EMAIL = "plantomics@example.com"

# Databases to search simultaneously
SEARCH_DATABASES = [
    ("pubmed", "PubMed abstracts"),
    ("pmc", "PubMed Central full-text"),
    ("gene", "NCBI Gene database"),
    ("nucleotide", "Nucleotide sequences + annotations"),
]


def _search_ncbi_db(db: str, query: str, retmax: int = 5) -> dict:
    """Search a single NCBI database. Returns (db, label, id_list, count)."""
    try:
        resp = requests.get(
            NCBI_BASE + "esearch.fcgi",
            params={
                "db": db,
                "term": query,
                "retmax": retmax,
                "retmode": "json",
                "tool": TOOL_NAME,
                "email": EMAIL,
            },
            timeout=10,
        )
        if not resp.ok:
            return {"db": db, "label": "", "ids": [], "count": 0}
        data = resp.json()
        result = data.get("esearchresult", {})
        return {
            "db": db,
            "label": "",
            "ids": [str(i) for i in result.get("idlist", [])],
            "count": int(result.get("count", 0)),
        }
    except Exception:
        return {"db": db, "label": "", "ids": [], "count": 0}


def _fetch_pubmed_summaries(pmids: list[str]) -> list[dict]:
    """Fetch titles for PubMed/PMC IDs."""
    if not pmids:
        return []
    try:
        resp = requests.get(
            NCBI_BASE + "esummary.fcgi",
            params={
                "db": "pubmed",
                "id": ",".join(pmids[:5]),
                "retmode": "json",
                "tool": TOOL_NAME,
                "email": EMAIL,
            },
            timeout=10,
        )
        if not resp.ok:
            return []
        data = resp.json()
        results = []
        for pid in pmids[:5]:
            rec = data.get("result", {}).get(pid, {})
            if rec:
                title = rec.get("title", "")
                source = rec.get("source", "")
                pubdate = rec.get("pubdate", "")
                results.append({"title": title, "source": source, "date": pubdate, "id": pid})
        return results
    except Exception:
        return []


def _search_mygene(query: str) -> list[dict]:
    """Search MyGene.info as a gene-specific fallback."""
    try:
        resp = requests.get(
            "https://mygene.info/v3/query",
            params={
                "q": query,
                "fields": "symbol,name,taxid,alias,summary",
                "size": 5,
            },
            timeout=10,
            headers={"User-Agent": f"{TOOL_NAME}/4.5"},
        )
        if not resp.ok:
            return []
        hits = resp.json().get("hits", [])
        results = []
        for h in hits[:5]:
            results.append({
                "title": h.get("symbol", "") + " — " + (h.get("name", "")),
                "source": f"MyGene.info (taxid: {h.get('taxid', '')})",
                "date": "",
                "id": str(h.get("_id", "")),
                "aliases": h.get("alias", []),
                "summary": h.get("summary", "")[:200],
            })
        return results
    except Exception:
        return []


def _search_semantic_scholar(query: str, max_results: int = 5) -> list[dict]:
    """Search Semantic Scholar API for academic papers.

    Semantic Scholar is a free, open-access academic search engine covering
    all scientific disciplines. It complements NCBI PubMed with broader coverage,
    citation data, and often more recent papers.

    Used as a supplementary source when NCBI results are sparse.
    """
    try:
        resp = requests.get(
            "https://api.semanticscholar.org/graph/v1/paper/search",
            params={
                "query": query,
                "limit": max_results,
                "fields": "title,abstract,year,authors,externalIds,url,publicationVenue",
            },
            timeout=12,
            headers={"User-Agent": f"{TOOL_NAME}/4.5"},
        )
        if not resp.ok:
            logger.debug("Semantic Scholar returned %d: %s", resp.status_code, resp.text[:100])
            return []

        data = resp.json()
        papers = data.get("data", [])

        results = []
        for paper in papers[:max_results]:
            authors = paper.get("authors", [])
            author_names = [a.get("name", "") for a in authors[:5]]

            external_ids = paper.get("externalIds", {})
            doi = external_ids.get("DOI", "")
            pmid = external_ids.get("PubMed", "")

            results.append({
                "title": paper.get("title", ""),
                "abstract": (paper.get("abstract") or "")[:300],
                "year": paper.get("year", ""),
                "authors": author_names,
                "venue": (paper.get("publicationVenue") or {}).get("name", ""),
                "url": paper.get("url", ""),
                "doi": doi,
                "pmid": pmid,
                "id": paper.get("paperId", ""),
                "source": "Semantic Scholar",
            })

        logger.info("Semantic Scholar: %d results for '%s'", len(results), query[:60])
        return results

    except Exception as e:
        logger.debug("Semantic Scholar search failed: %s", e)
        return []


class WebSearchTool:
    """Broad-spectrum search across NCBI databases + MyGene.info + Semantic Scholar.

    Multi-source search with progressive enrichment:
      1. NCBI E-utilities (PubMed, PMC, Gene, Nucleotide) — in parallel
      2. Semantic Scholar API — academic papers, complements NCBI
      3. MyGene.info — gene-specific lookup

    All sources are free and academic-focused. No DuckDuckGo dependency.

    Usage:
        tool = WebSearchTool()
        result = tool.search("Os01g0884300 NAC rice transcription factor")
        print(tool.format_for_llm(result))
    """

    def __init__(self, timeout: int = 10):
        self.timeout = timeout
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": f"{TOOL_NAME}/4.5"})

    def search(self, query: str, max_results: int = 5) -> dict[str, Any]:
        """Search across multiple sources in parallel.

        Returns:
            dict with:
                - query: original query
                - db_results: list of per-database NCBI results
                - summaries: PubMed/PMC article summaries
                - semantic_scholar_hits: Semantic Scholar academic papers
                - mygene_hits: MyGene.info gene matches
                - total_found: total count across NCBI databases
        """
        if not query or not query.strip():
            return {
                "query": query,
                "db_results": [],
                "summaries": [],
                "semantic_scholar_hits": [],
                "mygene_hits": [],
                "total_found": 0,
            }

        query = query.strip()
        cache_key = f"web_search:{query}:{max_results}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        db_results: list[dict] = []
        all_ids: set[str] = set()

        # Search NCBI databases + Semantic Scholar in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            # NCBI databases
            ncbi_futures = {
                executor.submit(_search_ncbi_db, db, query, max_results): (db, label)
                for db, label in SEARCH_DATABASES
            }
            # Semantic Scholar (runs in parallel with NCBI)
            sem_scholar_future = executor.submit(
                _search_semantic_scholar, query, max_results
            )

            for future in concurrent.futures.as_completed(ncbi_futures):
                try:
                    result = future.result(timeout=10)
                    for db_name, label in SEARCH_DATABASES:
                        if result["db"] == db_name:
                            result["label"] = label
                            break
                    if result["ids"]:
                        db_results.append(result)
                        all_ids.update(result["ids"])
                except Exception as e:
                    logger.debug("NCBI db search failed: %s", e)

            # Collect Semantic Scholar results
            semantic_scholar_hits: list[dict] = []
            try:
                semantic_scholar_hits = sem_scholar_future.result(timeout=15)
            except Exception as e:
                logger.debug("Semantic Scholar search failed: %s", e)

        # Sort by count descending
        db_results.sort(key=lambda r: r["count"], reverse=True)

        # Fetch summaries for PubMed/PMC results
        summaries: list[dict] = []
        if all_ids:
            pmids = [i for i in list(all_ids)[:5] if len(i) <= 12]
            if pmids:
                summaries = _fetch_pubmed_summaries(pmids)
                time.sleep(0.35)

        # MyGene.info search (gene-specific)
        mygene_hits = _search_mygene(query)

        total_found = sum(r["count"] for r in db_results)

        result = {
            "query": query,
            "db_results": db_results,
            "summaries": summaries,
            "semantic_scholar_hits": semantic_scholar_hits,
            "mygene_hits": mygene_hits,
            "total_found": total_found,
        }

        knowledge_cache.set(cache_key, result, ttl=3600)
        logger.info(
            "Web search '%s': %d total across %d NCBI DBs, %d Semantic Scholar, %d MyGene hits",
            query[:60], total_found, len(db_results),
            len(semantic_scholar_hits), len(mygene_hits),
        )
        return result

    def format_for_llm(self, result: dict[str, Any]) -> str:
        """Format search results for LLM consumption."""
        query = result.get("query", "")
        lines = [f"Web search results for '{query}':", ""]

        db_results = result.get("db_results", [])
        summaries = result.get("summaries", [])
        sem_scholar_hits = result.get("semantic_scholar_hits", [])
        mygene_hits = result.get("mygene_hits", [])

        if not db_results and not mygene_hits and not sem_scholar_hits:
            lines.append("No results found across NCBI databases, Semantic Scholar, or MyGene.info.")
            return "\n".join(lines)

        # NCBI database summary
        if db_results:
            total = result.get("total_found", 0)
            lines.append(f"**NCBI** ({total} total across {len(db_results)} databases):")
            for db_r in db_results[:4]:
                icon = {"pubmed": "[PubMed]", "pmc": "[PMC]", "gene": "[Gene]", "nucleotide": "[Nucleotide]"}.get(db_r["db"], "")
                lines.append(f"  {icon} {db_r['label']}: {db_r['count']} results")
            lines.append("")

        # NCBI article summaries
        if summaries:
            lines.append("**Top articles (PubMed)**:")
            for i, s in enumerate(summaries, 1):
                lines.append(f"{i}. {s['title'][:200]}")
                lines.append(f"   {s.get('source', '')} ({s.get('date', '')}) — ID: {s['id']}")
                lines.append("")
        elif db_results:
            lines.append("(No article summaries available — try searching PubMed directly)")
            lines.append("")

        # Semantic Scholar results
        if sem_scholar_hits:
            lines.append("**Academic papers (Semantic Scholar)**:")
            for i, h in enumerate(sem_scholar_hits, 1):
                authors = ", ".join(h.get("authors", [])[:3])
                year = h.get("year", "")
                lines.append(f"{i}. {h['title'][:200]}")
                if authors:
                    lines.append(f"   Authors: {authors} ({year})")
                if h.get("abstract"):
                    lines.append(f"   Abstract: {h['abstract'][:250]}")
                if h.get("venue"):
                    lines.append(f"   Published in: {h['venue']}")
                if h.get("pmid"):
                    lines.append(f"   PMID: {h['pmid']}")
                if h.get("doi"):
                    lines.append(f"   DOI: {h['doi']}")
                if h.get("url"):
                    lines.append(f"   URL: {h['url']}")
                lines.append("")

        # MyGene.info hits
        if mygene_hits:
            lines.append("**Gene matches (MyGene.info)**:")
            for i, h in enumerate(mygene_hits, 1):
                lines.append(f"{i}. {h['title'][:200]}")
                if h.get("aliases"):
                    lines.append(f"   Aliases: {', '.join(h['aliases'][:6])}")
                if h.get("summary"):
                    lines.append(f"   {h['summary'][:200]}")
                lines.append("")

        return "\n".join(lines)


def web_search(query: str, max_results: int = 5) -> str:
    """Search NCBI + MyGene.info and return formatted text for LLM."""
    tool = WebSearchTool()
    result = tool.search(query, max_results=max_results)
    return tool.format_for_llm(result)
