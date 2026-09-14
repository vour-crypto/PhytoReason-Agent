"""
pubmed_retriever.py — Structured PubMed retrieval for zero-data mode.

Smart query building: automatically constructs TF + metabolite + species
queries from user intent.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger("pubmed_retriever")


@dataclass
class LiteratureResult:
    """Structured literature search result."""
    pmid: str
    title: str = ""
    snippet: str = ""
    authors: str = ""
    year: str = ""
    journal: str = ""
    relevance: str = ""  # "high" | "medium" | "low"
    key_finding: str = ""


@dataclass
class LiteratureSearchResult:
    query: str
    total_hits: int
    results: list[LiteratureResult] = field(default_factory=list)
    query_note: str = ""


class PubMedRetriever:
    """Structured PubMed retriever with plant-science-aware query templates.

    Builds precise queries from user intent rather than passing raw strings.
    """

    # Query templates for common plant biology questions
    TEMPLATES = {
        "tf_regulation": (
            '("{tf_family}"[All Fields] OR "transcription factor"[All Fields]) '
            'AND ("{metabolite}"[All Fields] OR "{pathway}"[All Fields]) '
            'AND ("biosynthesis"[All Fields] OR "regulation"[All Fields] OR "metabolism"[All Fields])'
        ),
        "species_tf": (
            '("{species}"[All Fields]) AND ("transcription factor"[All Fields] OR "{tf_family}"[All Fields]) '
            'AND ("{metabolite}"[All Fields] OR "secondary metabolism"[All Fields])'
        ),
        "pathway_enzymes": (
            '("{metabolite}"[All Fields] OR "{pathway}"[All Fields]) '
            'AND ("biosynthesis"[All Fields] OR "pathway"[All Fields]) '
            'AND ("gene"[All Fields] OR "enzyme"[All Fields])'
        ),
        "general_metabolite": (
            '"{metabolite}"[All Fields] AND '
            '("regulation"[All Fields] OR "biosynthesis"[All Fields] OR "transcription"[All Fields])'
        ),
    }

    def __init__(self) -> None:
        self._cache: dict[str, LiteratureSearchResult] = {}

    def search(
        self,
        query: str = "",
        metabolite: str = "",
        tf_family: str = "",
        species: str = "",
        pathway: str = "",
        template: str = "general_metabolite",
        max_results: int = 10,
        use_cache: bool = True,
    ) -> LiteratureSearchResult:
        """Search PubMed with smart query construction.

        Precedence: explicit query > template-based construction.
        """
        if query:
            search_query = query
        else:
            search_query = self._build_query(
                metabolite=metabolite, tf_family=tf_family,
                species=species, pathway=pathway, template=template,
            )

        cache_key = f"{search_query}:{max_results}"
        if use_cache and cache_key in self._cache:
            return self._cache[cache_key]

        result = self._execute_search(search_query, max_results)
        if use_cache:
            self._cache[cache_key] = result
        return result

    def smart_search(
        self, topic: str, max_results: int = 10
    ) -> LiteratureSearchResult:
        """High-level: deduce intent and search.

        'topic' is a natural-language topic like:
        "berberine regulation by MYB in Coptis"
        """
        # Simple keyword extraction (LLM does better, this is fallback)
        topic_lower = topic.lower()
        query_parts = [topic]

        # If the topic mentions a metabolite, add focused terms
        if any(w in topic_lower for w in ["regulation", "调控", "regulate"]):
            query_parts.append("(regulation OR biosynthesis OR metabolism)")
        if any(w in topic_lower for w in ["biosynthesis", "合成", "pathway", "通路"]):
            query_parts.append("(biosynthesis OR pathway)")

        search_query = " AND ".join(f'("{p.strip()}"[All Fields])' for p in topic.split() if len(p) > 2)
        if len(search_query) < 20:
            search_query = topic

        return self._execute_search(search_query, max_results)

    def _build_query(
        self,
        metabolite: str = "",
        tf_family: str = "",
        species: str = "",
        pathway: str = "",
        template: str = "general_metabolite",
    ) -> str:
        """Build a PubMed query from structured parameters."""
        tmpl = self.TEMPLATES.get(template, self.TEMPLATES["general_metabolite"])
        return tmpl.format(
            metabolite=metabolite or "secondary metabolite",
            tf_family=tf_family or "transcription factor",
            species=species or "plant",
            pathway=pathway or "biosynthesis",
        )

    def _execute_search(
        self, query: str, max_results: int = 10
    ) -> LiteratureSearchResult:
        """Execute PubMed search via the existing PubMedSearch client."""
        try:
            from phyto_reason.tools.literature.pubmed_client import PubMedSearch

            client = PubMedSearch()
            raw = client.search(query, max_results=min(max_results, 20))

            results: list[LiteratureResult] = []
            for item in raw if isinstance(raw, list) else []:
                if isinstance(item, dict):
                    results.append(LiteratureResult(
                        pmid=str(item.get("pmid", item.get("uid", ""))),
                        title=item.get("title", "")[:200],
                        snippet=item.get("snippet", item.get("abstract", ""))[:300],
                        authors=item.get("authors", ""),
                        year=item.get("year", item.get("pubdate", "")),
                        journal=item.get("journal", item.get("source", "")),
                    ))

            total = len(results)
            return LiteratureSearchResult(
                query=query,
                total_hits=total,
                results=results[:max_results],
                query_note=(
                    f"PubMed returned {total} results for: {query[:100]}"
                    if total > 0 else f"No PubMed results for: {query[:100]}"
                ),
            )

        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            return LiteratureSearchResult(
                query=query,
                total_hits=0,
                query_note=f"PubMed search unavailable: {e}",
            )
