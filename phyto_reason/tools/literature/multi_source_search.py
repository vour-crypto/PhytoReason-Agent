"""
multi_source_search.py — Multi-source academic literature search engine.

Aggregates results from 6 sources in parallel:
  1. PubMed (NCBI E-utilities)
  2. Semantic Scholar (academic papers with citations)
  3. Europe PMC (PubMed + PMC + Agricola + preprints)
  4. OpenAlex (free open index of 250M+ works — aggregates WoS/Scopus/Crossref)
  5. arXiv (preprints — q-bio, q-genomics, bioinformatics)
  6. Local RAG (previous analysis results, reports, cached data)

All sources are queried simultaneously via ThreadPoolExecutor. Each source
degrades gracefully — failure of one does not block results from others.

Usage:
    engine = MultiSourceLiteratureSearch()
    result = engine.search("berberine biosynthesis MYB Coptis", max_results=5)
    print(result)  # formatted text ready for LLM consumption
"""

from __future__ import annotations

import concurrent.futures
import logging
import re
from typing import Any

from phyto_reason.utils.cache import knowledge_cache

logger = logging.getLogger("multi_source_search")


# ── Source-specific search functions ──────────────────────────

def _search_pubmed(query: str, max_results: int) -> list[dict]:
    """Search PubMed via existing PubMedSearch client."""
    try:
        from phyto_reason.tools.literature.pubmed_client import PubMedSearch

        client = PubMedSearch()
        return client.search(query, max_results=max_results)
    except Exception as e:
        logger.debug("PubMed search failed in multi-source: %s", e)
        return []


def _search_semantic_scholar(query: str, max_results: int) -> list[dict]:
    """Search Semantic Scholar (reuse from web_search module)."""
    try:
        from phyto_reason.tools.web_search import _search_semantic_scholar
        return _search_semantic_scholar(query, max_results)
    except Exception as e:
        logger.debug("Semantic Scholar failed in multi-source: %s", e)
        return []


def _search_europe_pmc(query: str, max_results: int) -> list[dict]:
    """Search Europe PMC."""
    try:
        from phyto_reason.tools.literature.europe_pmc import _search_europe_pmc
        return _search_europe_pmc(query, max_results)
    except Exception as e:
        logger.debug("Europe PMC failed in multi-source: %s", e)
        return []


def _search_local_rag(query: str, max_results: int) -> list[dict]:
    """Search local RAG document store."""
    try:
        from phyto_reason.tools.rag import get_rag_engine

        engine = get_rag_engine()
        results = engine.search(query, top_k=max_results)
        return results
    except Exception as e:
        logger.debug("Local RAG failed in multi-source: %s", e)
        return []


def _search_openalex(query: str, max_results: int) -> list[dict]:
    """Search OpenAlex — free, open index of 250M+ scholarly works.

    OpenAlex aggregates data from Crossref, PubMed, ORCID, ROR, Unpaywall, etc.
    It includes coverage of Web of Science and Scopus metadata via their
    institutional partnerships. No API key required.

    API: https://docs.openalex.org/api
    """
    try:
        import requests

        resp = requests.get(
            "https://api.openalex.org/works",
            params={
                "search": query,
                "per_page": min(max_results, 10),
                "sort": "relevance_score:desc",
            },
            timeout=12,
            headers={"User-Agent": "plantomics-agent/4.5", "Accept": "application/json"},
        )

        if not resp.ok:
            logger.debug("OpenAlex returned %d: %s", resp.status_code, resp.text[:200])
            return []

        data = resp.json()
        entries = data.get("results", [])

        results: list[dict] = []
        for entry in entries:
            # Reconstruct abstract from inverted index
            abstract = ""
            inverted = entry.get("abstract_inverted_index")
            if inverted:
                words = []
                for word, positions in inverted.items():
                    for pos in positions:
                        words.append((pos, word))
                words.sort(key=lambda x: x[0])
                abstract = " ".join(w for _, w in words)[:300]

            # Authors
            authors = []
            for a in entry.get("authorships", []):
                name = a.get("author", {}).get("display_name", "")
                if name:
                    authors.append(name)

            # Journal / venue
            loc = entry.get("primary_location") or {}
            source = loc.get("source") or {}
            journal = source.get("display_name", "") or ""
            oa_status = source.get("open_access", {}).get("oa_status", "")

            # OA URL
            oa = entry.get("open_access") or {}
            oa_url = oa.get("oa_url", "") or ""

            # PMID is inside ids dict as a URL
            ids = entry.get("ids") or {}
            pmid_raw = ids.get("pmid", "") or ""
            pmid = pmid_raw.replace("https://pubmed.ncbi.nlm.nih.gov/", "")

            results.append({
                "title": (entry.get("title", "") or ""),
                "authors": authors[:5],
                "year": str(entry.get("publication_year", "") or ""),
                "journal": journal,
                "doi": (entry.get("doi", "") or "").replace("https://doi.org/", ""),
                "pmid": pmid,
                "abstract": abstract,
                "source": "OpenAlex",
                "cited_count": int(entry.get("cited_by_count", 0)),
                "open_access": oa_status,
                "oa_url": oa_url,
                "type": entry.get("type", ""),
                "id": entry.get("id", ""),
            })

        logger.info(
            "OpenAlex: %d results for '%s'", len(results), query[:60]
        )
        return results

    except Exception as e:
        logger.debug("OpenAlex search failed: %s", e)
        return []


def _search_arxiv(query: str, max_results: int) -> list[dict]:
    """Search arXiv preprints via the arXiv API.

    arXiv hosts preprints in quantitative biology (q-bio), including
    plant science, genomics, and bioinformatics.

    API: https://info.arxiv.org/help/api/index.html
    """
    try:
        import urllib.request
        import urllib.parse
        import xml.etree.ElementTree as ET

        # arXiv API uses its own query syntax: replace empty spaces with AND
        # and wrap phrases in quotes
        arxiv_query = query.strip()
        if len(arxiv_query.split()) > 1:
            # Group into a phrase query for better precision
            arxiv_query = "(" + " AND ".join(arxiv_query.split()[:6]) + ")"

        url = (
            "https://export.arxiv.org/api/query?"
            + urllib.parse.urlencode({
                "search_query": f"all:{arxiv_query}",
                "start": 0,
                "max_results": min(max_results, 10),
                "sortBy": "relevance",
                "sortOrder": "descending",
            })
        )

        req = urllib.request.Request(
            url,
            headers={"User-Agent": "plantomics-agent/4.5"},
        )

        with urllib.request.urlopen(req, timeout=12) as resp:
            xml_data = resp.read().decode("utf-8")

        # arXiv API returns Atom XML
        ns = {
            "atom": "http://www.w3.org/2005/Atom",
            "arxiv": "http://arxiv.org/schemas/atom",
        }

        root = ET.fromstring(xml_data)
        entries = root.findall("atom:entry", ns)

        results: list[dict] = []
        for entry in entries:
            title = (entry.findtext("atom:title", "", ns) or "").replace("\n", " ").strip()

            # Authors
            authors = []
            for author_elem in entry.findall("atom:author", ns):
                name = author_elem.findtext("atom:name", "", ns)
                if name:
                    authors.append(name)

            # Published year
            published = entry.findtext("atom:published", "", ns)
            year = published[:4] if published else ""

            # arXiv ID
            arxiv_id = entry.findtext("atom:id", "", ns)
            arxiv_id = arxiv_id.replace("http://arxiv.org/abs/", "") if arxiv_id else ""

            # DOI
            doi = ""
            for link in entry.findall("atom:link", ns):
                href = link.get("href", "")
                if "doi.org" in href:
                    doi = href.replace("https://doi.org/", "").replace("http://dx.doi.org/", "")
                    break

            # Categories
            cats = [c.get("term", "") for c in entry.findall("atom:category", ns)]

            # Summary
            summary = (entry.findtext("atom:summary", "", ns) or "").replace("\n", " ").strip()[:300]

            results.append({
                "title": title,
                "authors": authors[:5],
                "year": year,
                "journal": "arXiv",
                "doi": doi,
                "pmid": "",
                "abstract": summary,
                "source": "arXiv",
                "arxiv_id": arxiv_id,
                "categories": cats[:5],
            })

        logger.info(
            "arXiv: %d results for '%s'", len(results), query[:60]
        )
        return results

    except Exception as e:
        logger.debug("arXiv search failed: %s", e)
        return []


# ── Formatting helpers ────────────────────────────────────────

def _format_pubmed(results: list[dict]) -> str:
    if not results:
        return ""
    lines = [f"### PubMed ({len(results)} results)", ""]
    for i, r in enumerate(results, 1):
        title = r.get("Title", r.get("title", "?"))
        source = r.get("Source", r.get("journal", r.get("source", "")))
        date = r.get("PubDate", r.get("year", r.get("date", "")))
        pmid = r.get("PMID", r.get("pmid", ""))
        doi = r.get("DOI", r.get("doi", ""))

        lines.append(f"{i}. **{title[:200]}**")
        cite = f"   *{source}* ({date})"
        lines.append(cite)
        id_line = f"   PMID: {pmid}"
        if doi:
            id_line += f" | DOI: {doi}"
        lines.append(id_line)
        lines.append("")
    return "\n".join(lines)


def _format_semantic_scholar(results: list[dict]) -> str:
    if not results:
        return ""
    lines = [f"### Semantic Scholar ({len(results)} results)", ""]
    for i, r in enumerate(results, 1):
        title = r.get("title", "?")
        authors = ", ".join(r.get("authors", [])[:3])
        year = r.get("year", "")
        venue = r.get("venue", "")
        doi = r.get("doi", "")
        pmid = r.get("pmid", "")
        abstract = r.get("abstract", "")

        lines.append(f"{i}. **{title[:200]}**")
        if authors:
            cite = f"   {authors} ({year})"
            if venue:
                cite += f" — *{venue}*"
            lines.append(cite)
        if pmid:
            lines.append(f"   PMID: {pmid}")
        if doi:
            lines.append(f"   DOI: {doi}")
        if abstract:
            lines.append(f"   > {abstract[:250]}")
        lines.append("")
    return "\n".join(lines)


def _format_europe_pmc(results: list[dict]) -> str:
    if not results:
        return ""
    lines = [f"### Europe PMC ({len(results)} results)", ""]
    for i, r in enumerate(results, 1):
        title = r.get("title", "?")
        authors = ", ".join(r.get("authors", [])[:3])
        year = r.get("year", "")
        journal = r.get("journal", "")
        pmid = r.get("pmid", "")
        doi = r.get("doi", "")
        abstract = r.get("abstract", "")
        src_type = r.get("source_type", "")

        lines.append(f"{i}. **{title[:200]}**")
        if authors:
            cite = f"   {authors} ({year})"
            if journal:
                cite += f" — *{journal}*"
            if src_type and src_type != "MED":
                cite += f" [{src_type}]"
            lines.append(cite)
        if pmid:
            id_line = f"   PMID: {pmid}"
            if doi:
                id_line += f" | DOI: {doi}"
            lines.append(id_line)
        elif doi:
            lines.append(f"   DOI: {doi}")
        if abstract:
            lines.append(f"   > {abstract[:250]}")
        lines.append("")
    return "\n".join(lines)


def _format_openalex(results: list[dict]) -> str:
    if not results:
        return ""
    lines = [f"### OpenAlex ({len(results)} results)", ""]
    for i, r in enumerate(results, 1):
        title = r.get("title", "?")
        authors = ", ".join(r.get("authors", [])[:3])
        year = r.get("year", "")
        journal = r.get("journal", "")
        doi = r.get("doi", "")
        pmid = r.get("pmid", "")
        abstract = r.get("abstract", "")
        cited = r.get("cited_count", 0)
        oa = r.get("open_access", "")
        c = "🔓" if oa else "🔒"

        lines.append(f"{i}. **{title[:200]}**")
        if authors:
            cite = f"   {authors} ({year})"
            if journal:
                cite += f" — *{journal}*"
            lines.append(cite)
        id_line_parts = []
        if pmid:
            id_line_parts.append(f"PMID: {pmid}")
        if doi:
            id_line_parts.append(f"DOI: {doi}")
        if id_line_parts:
            lines.append(f"   {' | '.join(id_line_parts)}")
        if cited > 0:
            lines.append(f"   Cited: {cited} times {c}")
        if abstract:
            lines.append(f"   > {abstract[:250]}")
        lines.append("")
    return "\n".join(lines)


def _format_arxiv(results: list[dict]) -> str:
    if not results:
        return ""
    lines = [f"### arXiv ({len(results)} preprints)", ""]
    for i, r in enumerate(results, 1):
        title = r.get("title", "?")
        authors = ", ".join(r.get("authors", [])[:3])
        year = r.get("year", "")
        arxiv_id = r.get("arxiv_id", "")
        doi = r.get("doi", "")
        categories = r.get("categories", [])
        abstract = r.get("abstract", "")

        lines.append(f"{i}. **{title[:200]}**")
        if authors:
            lines.append(f"   {authors} ({year})")
        id_line = f"   arXiv: {arxiv_id}"
        if doi:
            id_line += f" | DOI: {doi}"
        lines.append(id_line)
        if categories:
            lines.append(f"   Categories: {', '.join(categories[:4])}")
        if abstract:
            lines.append(f"   > {abstract[:250]}")
        lines.append("")
    return "\n".join(lines)


def _format_rag(results: list[dict]) -> str:
    if not results:
        return ""
    lines = [f"### 📂 Local Knowledge Base ({len(results)} matches from previous analyses)", ""]
    for i, r in enumerate(results, 1):
        meta = r.get("metadata", {})
        src = meta.get("source", "?")
        section = meta.get("section", "")
        col = meta.get("column", "")
        score = r.get("score", 0)

        loc = src
        if section:
            loc += f" → {section}"
        if col:
            loc += f" | col: {col}"

        lines.append(f"{i}. [{loc}] (relevance: {score:.2f})")
        text = r.get("text", "")[:300]
        lines.append(f"   {text}")
        lines.append("")
    return "\n".join(lines)


# ── Reference list ────────────────────────────────────────────

def _normalize_ref(r: dict) -> dict:
    """Normalize a result dict from any source into a uniform reference shape."""
    title = (r.get("Title") or r.get("title") or "").replace("<i>", "*").replace("</i>", "*").strip()
    authors_list = r.get("authors") or []
    if not authors_list:
        author_str = (r.get("authorString") or r.get("authors", "") or "")
        if isinstance(author_str, str):
            authors_list = [a.strip() for a in author_str.rstrip(".").split(",") if a.strip()]
    authors = ", ".join(authors_list[:5]) if authors_list else ""

    year = r.get("year") or r.get("PubDate") or r.get("publication_year") or ""
    if isinstance(year, int):
        year = str(year)
    journal = (r.get("journal") or r.get("Source") or r.get("source") or
               r.get("venue") or r.get("Journal") or "")
    # Skip vol/issue numbers masquerading as journal names
    if re.fullmatch(r"v\d+(\(\d+\))?", journal.strip()):
        journal = ""
    doi = r.get("doi") or r.get("DOI") or ""
    pmid = r.get("pmid") or r.get("PMID") or ""
    url = r.get("url") or r.get("oa_url") or ""
    arxiv_id = r.get("arxiv_id") or ""
    source_name = r.get("source", "?")

    return {
        "title": title, "authors": authors, "year": year, "journal": journal,
        "doi": doi, "pmid": pmid, "url": url, "arxiv_id": arxiv_id,
        "source_name": source_name,
    }


def _format_references(all_results: dict[str, list[dict]]) -> str:
    """Build a unified numbered reference list from all source results.

    Deduplicates by PMID/DOI — if the same paper appears in multiple
    sources, only the first occurrence is listed.
    """
    seen_dois: set[str] = set()
    seen_pmids: set[str] = set()
    seen_titles: set[str] = set()
    refs: list[dict] = []
    # Track dedup keys → index in refs for merging on duplicate
    ref_by_key: dict[str, int] = {}

    for source_key in ["pubmed", "semantic_scholar", "europe_pmc", "openalex", "arxiv"]:
        for r in all_results.get(source_key, []):
            ref = _normalize_ref(r)
            doi = ref["doi"].lower().strip()
            pmid = ref["pmid"].strip()
            title_key = re.sub(r"[<>*/\[\]\(\)&\"]+", "", ref["title"].lower()).strip().rstrip(".")[:80]

            # Determine dedup key (prefer DOI > PMID > title)
            dedup_key = doi or pmid or title_key
            if not dedup_key or not ref["title"]:
                continue

            if dedup_key in ref_by_key:
                # Merge: keep the entry with more complete data
                existing_idx = ref_by_key[dedup_key]
                existing = refs[existing_idx]
                for field in ["authors", "journal", "doi", "pmid", "url"]:
                    if not existing.get(field) and ref.get(field):
                        existing[field] = ref[field]
                continue

            # Store ALL possible dedup keys pointing to this reference
            idx = len(refs)
            refs.append(ref)
            if doi:
                ref_by_key[doi] = idx
            if pmid:
                ref_by_key[pmid] = idx
            if title_key and title_key not in ref_by_key:
                ref_by_key[title_key] = idx

    if not refs:
        return ""

    lines = ["", "## 📚 References", ""]
    for i, r in enumerate(refs, 1):
        parts = []
        if r["authors"]:
            parts.append(f"{r['authors']}")
            if r["year"]:
                parts[-1] += f" ({r['year']})"
        elif r["year"]:
            parts.append(f"({r['year']})")
        parts.append(f"*{r['title']}*")
        if r["journal"]:
            parts.append(f"*{r['journal']}*")
        lines.append(f"{i:2d}. {' '.join(parts)}")

        id_parts = []
        if r["pmid"]:
            id_parts.append(f"PMID: {r['pmid']}")
        if r["doi"]:
            id_parts.append(f"DOI: {r['doi']}")
        if r["arxiv_id"]:
            id_parts.append(f"arXiv: {r['arxiv_id']}")
        if r["url"]:
            id_parts.append(f"URL: {r['url']}")
        if id_parts:
            lines.append(f"      {' | '.join(id_parts)}")
        lines.append("")

    return "\n".join(lines)


# ── Main engine ───────────────────────────────────────────────

class MultiSourceLiteratureSearch:
    """Unified literature search across PubMed + Semantic Scholar + Europe PMC + OpenAlex + arXiv + RAG.

    All sources are queried in parallel. Each source degrades independently —
    if one fails, results from others are still returned.

    Usage:
        engine = MultiSourceLiteratureSearch()
        formatted_text = engine.search("berberine MYB biosynthesis Coptis", max_results=5)
    """

    def __init__(self) -> None:
        self._max_workers = 6

    def search(self, query: str, max_results: int = 5) -> str:
        """Search all academic sources + local RAG in parallel.

        Args:
            query: Search query string.
            max_results: Max results per source (capped at 10).

        Returns:
            Formatted markdown text ready for LLM consumption, with all
            source results clearly labeled in separate sections.
        """
        if not query or not query.strip():
            return "Error: search query is required."

        query = query.strip()
        max_r = min(max_results, 10)

        # Check cache
        cache_key = f"multi_lit:{query}:{max_r}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        # ── Parallel search across all 6 sources ──────────────
        pubmed_results: list[dict] = []
        sem_scholar_results: list[dict] = []
        epmc_results: list[dict] = []
        openalex_results: list[dict] = []
        arxiv_results: list[dict] = []
        rag_results: list[dict] = []

        with concurrent.futures.ThreadPoolExecutor(max_workers=self._max_workers) as executor:
            futures = {
                executor.submit(_search_pubmed, query, max_r): "pubmed",
                executor.submit(_search_semantic_scholar, query, max_r): "semantic_scholar",
                executor.submit(_search_europe_pmc, query, max_r): "europe_pmc",
                executor.submit(_search_openalex, query, max_r): "openalex",
                executor.submit(_search_arxiv, query, max_r): "arxiv",
                executor.submit(_search_local_rag, query, max_r): "rag",
            }

            for future in concurrent.futures.as_completed(futures):
                source = futures[future]
                try:
                    result = future.result(timeout=15)
                    if source == "pubmed":
                        pubmed_results = result or []
                    elif source == "semantic_scholar":
                        sem_scholar_results = result or []
                    elif source == "europe_pmc":
                        epmc_results = result or []
                    elif source == "openalex":
                        openalex_results = result or []
                    elif source == "arxiv":
                        arxiv_results = result or []
                    elif source == "rag":
                        rag_results = result or []
                except concurrent.futures.TimeoutError:
                    logger.debug("Source '%s' timed out during multi-source search", source)
                except Exception as e:
                    logger.debug("Source '%s' failed: %s", source, e)

        # ── Format results ────────────────────────────────────
        sections = []
        total_sources = 0

        # Header
        sections.append(f"## Literature Search Results for: '{query}'")
        sections.append("*Queried: PubMed, Semantic Scholar, Europe PMC, OpenAlex, arXiv, and Local Knowledge Base*")
        sections.append("")

        # PubMed
        if pubmed_results:
            sections.append(_format_pubmed(pubmed_results))
            total_sources += 1

        # Semantic Scholar
        if sem_scholar_results:
            sections.append(_format_semantic_scholar(sem_scholar_results))
            total_sources += 1

        # Europe PMC
        if epmc_results:
            sections.append(_format_europe_pmc(epmc_results))
            total_sources += 1

        # OpenAlex
        if openalex_results:
            sections.append(_format_openalex(openalex_results))
            total_sources += 1

        # arXiv
        if arxiv_results:
            sections.append(_format_arxiv(arxiv_results))
            total_sources += 1

        # Local RAG
        if rag_results:
            sections.append(_format_rag(rag_results))
            total_sources += 1

        if total_sources == 0:
            sections.append(
                "⚠️ No results found across PubMed, Semantic Scholar, Europe PMC, "
                "OpenAlex, arXiv, or the local knowledge base. Try broadening your "
                "query or using more general search terms."
            )

        # Unified reference list (deduplicated across all sources)
        all_raw = {
            "pubmed": pubmed_results,
            "semantic_scholar": sem_scholar_results,
            "europe_pmc": epmc_results,
            "openalex": openalex_results,
            "arxiv": arxiv_results,
        }
        refs_section = _format_references(all_raw)
        if refs_section:
            sections.append(refs_section)

        # Footer with guidance
        sections.append("---")
        sections.append(
            "*Results aggregated from multiple databases. Each source may "
            "have overlapping coverage. PMIDs and DOIs are provided for "
            "cross-referencing.*"
        )

        result = "\n".join(sections)
        knowledge_cache.set(cache_key, result, ttl=1800)  # cache 30 min
        return result

    def search_structured(
        self, query: str, max_results: int = 5
    ) -> dict[str, list[dict]]:
        """Search all sources and return structured results (no formatting).

        Useful for programmatic access. Returns dict with keys:
            pubmed, semantic_scholar, europe_pmc, openalex, arxiv, rag
        """
        max_r = min(max_results, 10)
        return {
            "query": query,
            "pubmed": _search_pubmed(query, max_r),
            "semantic_scholar": _search_semantic_scholar(query, max_r),
            "europe_pmc": _search_europe_pmc(query, max_r),
            "openalex": _search_openalex(query, max_r),
            "arxiv": _search_arxiv(query, max_r),
            "rag": _search_local_rag(query, max_r),
        }


# ── Convenience function ──────────────────────────────────────

def multi_source_literature_search(query: str, max_results: int = 5) -> str:
    """Convenience function: search all literature sources and return formatted text."""
    engine = MultiSourceLiteratureSearch()
    return engine.search(query, max_results)
