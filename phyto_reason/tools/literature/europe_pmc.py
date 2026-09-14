"""
europe_pmc.py — Europe PMC literature search client.

Europe PMC is a free, open-access academic search engine that indexes:
  - PubMed / MEDLINE
  - PubMed Central full-text
  - Agricola (USDA agricultural literature — great for plant science)
  - Preprints (bioRxiv, Research Square, etc.)
  - Patents and more

Its coverage of plant science literature is excellent and often broader
than PubMed alone. Zero API key required.

API docs: https://europepmc.org/RestfulWebService
"""

from __future__ import annotations

import logging
from typing import Any

import requests

logger = logging.getLogger("europe_pmc")

# Base URL for the search endpoint (lite result type = fast, no full text parsing)
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Query wrapper: include all article types (MED=PubMed, PMC=full text, AGR=Agricola, PPR=Preprints)
# We don't restrict source types at query time — instead we filter/annotate in the results
QUERY_PREFIX = ""


def _search_europe_pmc(query: str, max_results: int = 5) -> list[dict[str, Any]]:
    """Search Europe PMC and return structured results.

    Returns:
        List of dicts with keys:
            title, authors, year, journal, pmid, doi, abstract, source,
            cited_count, first_publication_date
    """
    if not query or not query.strip():
        return []

    # Build query
    full_query = query.strip()
    if QUERY_PREFIX:
        full_query = f"{QUERY_PREFIX} AND ({query.strip()})"

    try:
        resp = requests.get(
            EPMC_SEARCH,
            params={
                "query": full_query,
                "resultType": "core",
                "format": "json",
                "pageSize": min(max_results, 20),
            },
            timeout=12,
            headers={"User-Agent": "plantomics-agent/4.5"},
        )

        if not resp.ok:
            logger.debug(
                "Europe PMC returned %d: %s", resp.status_code, resp.text[:200]
            )
            return []

        data = resp.json()
        entries = data.get("resultList", {}).get("result", [])

        results: list[dict[str, Any]] = []
        for entry in entries:
            # Extract author names
            authors = []
            author_string = entry.get("authorString", "")
            if author_string:
                authors = [a.strip() for a in author_string.split(",") if a.strip()]

            # Extract citation info — Europe PMC doesn't return journalTitle in
            # core result type; try journalInfo and bookOrReportDetails
            source_type = entry.get("source", "")  # MED, PMC, AGR, PPR
            journal = ""
            journal_info = entry.get("journalInfo", {}) or {}
            if journal_info.get("volume"):
                journal = f"v{journal_info['volume']}"
                if journal_info.get("issue"):
                    journal = f"{journal}({journal_info['issue']})"

            # For book/report entries
            book_details = entry.get("bookOrReportDetails", {}) or {}
            if not journal and book_details.get("title"):
                journal = book_details["title"]

            doi = entry.get("doi", "") or ""
            pmid = ""
            pmcid = entry.get("pmcid", "") or ""

            # For MED (PubMed) entries, the 'id' field IS the PubMed ID
            # For PMC entries, 'id' is the PMC ID
            eid = str(entry.get("id", ""))
            if source_type == "MED" and eid and eid.isdigit():
                pmid = eid
            elif source_type == "PMC":
                pmcid = pmcid or eid

            # Year: prefer pubYear (MED), fall back to firstPublicationDate
            year = entry.get("pubYear", "") or ""
            if not year:
                first_date = entry.get("firstPublicationDate", "")
                if first_date and len(first_date) >= 4:
                    year = first_date[:4]

            # Abstract
            abstract = (entry.get("abstractText", "") or "")[:300]

            results.append({
                "title": (entry.get("title", "") or "").strip(),
                "authors": authors[:5],
                "year": year,
                "journal": journal,
                "pmid": pmid,
                "pmcid": pmcid,
                "doi": doi,
                "abstract": abstract,
                "source": "Europe PMC",
                "cited_count": int(entry.get("citedByCount", 0)),
                "first_publication_date": entry.get("firstPublicationDate", ""),
                "source_type": source_type,
            })

        logger.info(
            "Europe PMC: %d results for '%s'", len(results), query[:60]
        )
        return results

    except requests.exceptions.Timeout:
        logger.debug("Europe PMC request timed out for: %s", query[:60])
        return []
    except Exception as e:
        logger.debug("Europe PMC search failed: %s", e)
        return []


class EuropePMCSearch:
    """Europe PMC search client with caching.

    Usage:
        client = EuropePMCSearch()
        results = client.search("berberine biosynthesis MYB plant")
        formatted = client.format_for_llm(results)
    """

    def __init__(self) -> None:
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": "plantomics-agent/4.5"})

    def search(self, query: str, max_results: int = 5) -> list[dict[str, Any]]:
        """Search Europe PMC for academic literature."""
        return _search_europe_pmc(query, max_results)

    def format_for_llm(self, results: list[dict]) -> str:
        """Format Europe PMC results for LLM consumption."""
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
            cited = r.get("cited_count", 0)
            src_type = r.get("source_type", "")

            lines.append(f"{i}. **{title[:200]}**")
            if authors:
                cite_info = f"{authors} ({year})"
                if journal:
                    cite_info += f" — *{journal}*"
                lines.append(f"   {cite_info}")
            if pmid:
                line = f"   PMID: {pmid}"
                if doi:
                    line += f" | DOI: {doi}"
                lines.append(line)
            elif doi:
                lines.append(f"   DOI: {doi}")
            if cited > 0:
                lines.append(f"   Cited by: {cited}")
            abstract = r.get("abstract", "")
            if abstract:
                lines.append(f"   > {abstract[:250]}")
            lines.append("")

        return "\n".join(lines)
