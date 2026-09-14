"""
public_expression_client.py — Public expression database integration (v4.5).

Integrates with:
  - BAR (Bio-Analytic Resource) — bar.utoronto.ca
    * AGI → gene annotation
    * Tissue expression data
    * Expressologs (cross-species expression similarity)
    * Protein-protein interactions
  - eFP Browser — expression pictographs (image URLs)

Used to enrich ortholog results with public expression context.
All queries cached for responsiveness.

Reference: Sullivan et al. (2024), NAR 53(D1):D1576-D1586. PMID: 39441075
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from typing import Any

import requests

from phyto_reason.utils.cache import knowledge_cache

logger = logging.getLogger("public_expression")

BAR_BASE = "https://bar.utoronto.ca/webservices/"
TOOL_NAME = "phyto_reason"
TIMEOUT = 10  # seconds


# ── Data models ────────────────────────────────────────────────

@dataclass
class TissueExpression:
    """Expression level in a specific tissue/condition."""
    tissue: str = ""
    category: str = ""               # e.g., "root", "leaf", "flower"
    expression_level: float = 0.0    # normalized expression value
    is_expressed: bool = False       # above detection threshold


@dataclass
class GeneExpressionProfile:
    """Full expression profile for a gene across tissues."""
    gene_id: str = ""
    gene_symbol: str = ""
    species: str = ""
    data_source: str = ""            # "BAR", "eFP", etc.
    tissues: list[TissueExpression] = field(default_factory=list)
    top_expressed_tissues: list[str] = field(default_factory=list)
    note: str = ""


@dataclass
class ExpressologHit:
    """An expressolog — a gene with similar expression pattern across species."""
    query_gene: str = ""
    target_gene: str = ""
    target_species: str = ""
    expression_correlation: float = 0.0  # 0-1
    data_source: str = ""


@dataclass
class ExpressologResult:
    """Expressolog search results."""
    query_gene: str = ""
    hits: list[ExpressologHit] = field(default_factory=list)
    note: str = ""


# ── BAR Client ─────────────────────────────────────────────────

class BARClient:
    """Client for BAR (Bio-Analytic Resource for Plant Biology) web services.

    Provides gene annotation, tissue expression, expressologs, and PPI data.
    All methods cache results and fail gracefully.
    """

    def __init__(self, base_url: str = BAR_BASE):
        self.base_url = base_url
        self._session = requests.Session()
        self._session.headers.update({
            "User-Agent": f"{TOOL_NAME}/4.5",
        })

    # ── Gene annotation ────────────────────────────────────────

    def get_gene_annotation(self, gene_id: str) -> dict:
        """Fetch gene annotation from BAR (AGI → annotation).

        Returns dict with keys: agi, annotation, symbol, aliases, etc.
        """
        cache_key = f"bar_annot_{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached:
            return json.loads(cached) if isinstance(cached, str) else cached

        try:
            url = self.base_url + "agiToAnnot.php"
            resp = self._session.get(
                url,
                params={"request": json.dumps({"agi": gene_id})},
                timeout=TIMEOUT,
            )
            if resp.ok:
                data = resp.json()
                knowledge_cache.set(cache_key, json.dumps(data), ttl=86400 * 7)
                return data
        except Exception as e:
            logger.warning("BAR gene annotation failed for %s: %s", gene_id, e)

        return {}

    # ── Tissue expression ──────────────────────────────────────

    def get_tissue_expression(self, gene_id: str) -> GeneExpressionProfile:
        """Get tissue-specific expression profile for a gene.

        Uses BAR expression signal endpoint.
        Falls back to empty profile if unavailable.
        """
        cache_key = f"bar_tissue_expr_{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached:
            if isinstance(cached, str):
                return GeneExpressionProfile(**json.loads(cached))
            return cached

        profile = GeneExpressionProfile(gene_id=gene_id, data_source="BAR")

        try:
            # Try BAR expression signal endpoint
            url = self.base_url + "get_expression_data.php"
            resp = self._session.get(
                url,
                params={"request": json.dumps([{"gene": gene_id}])},
                timeout=TIMEOUT,
            )

            if resp.ok:
                data = resp.json()
                tissues = self._parse_tissue_expression(data, gene_id)
                profile.tissues = tissues
                profile.top_expressed_tissues = [
                    t.tissue for t in sorted(tissues, key=lambda x: x.expression_level, reverse=True)[:5]
                ]
                profile.note = f"Data from BAR ({len(tissues)} tissues)"
            else:
                profile.note = "BAR expression data unavailable (HTTP error)"

        except requests.Timeout:
            profile.note = "BAR expression data timed out"
        except Exception as e:
            profile.note = f"BAR expression data error: {e}"
            logger.warning("BAR tissue expression failed for %s: %s", gene_id, e)

        # Cache and return
        try:
            knowledge_cache.set(cache_key, json.dumps(profile.__dict__), ttl=86400)
        except Exception:
            pass

        return profile

    def _parse_tissue_expression(self, data: Any, gene_id: str) -> list[TissueExpression]:
        """Parse BAR expression data into TissueExpression list."""
        results: list[TissueExpression] = []
        try:
            if isinstance(data, list):
                for item in data:
                    if isinstance(item, dict):
                        for key, value in item.items():
                            if key not in ("gene", "agi", "request"):
                                try:
                                    level = float(value) if value is not None else 0.0
                                    results.append(TissueExpression(
                                        tissue=key,
                                        category=self._categorize_tissue(key),
                                        expression_level=level,
                                        is_expressed=level > 1.0,
                                    ))
                                except (ValueError, TypeError):
                                    pass
            elif isinstance(data, dict):
                for key, value in data.items():
                    if key not in ("gene", "agi", "request", "status"):
                        try:
                            level = float(value) if value is not None else 0.0
                            results.append(TissueExpression(
                                tissue=key,
                                category=self._categorize_tissue(key),
                                expression_level=level,
                                is_expressed=level > 1.0,
                            ))
                        except (ValueError, TypeError):
                            pass
        except Exception as e:
            logger.debug("Failed to parse tissue expression for %s: %s", gene_id, e)

        return results

    @staticmethod
    def _categorize_tissue(tissue_name: str) -> str:
        """Categorize tissue into broad category."""
        t = tissue_name.lower()
        if any(k in t for k in ("root", "radicle", "lateral root")):
            return "root"
        elif any(k in t for k in ("leaf", "rosette", "cotyledon", "cauline", "shoot")):
            return "shoot"
        elif any(k in t for k in ("flower", "petal", "sepal", "stamen", "carpel", "pollen", "anther", "inflorescence")):
            return "flower"
        elif any(k in t for k in ("seed", "embryo", "endosperm", "silique", "fruit")):
            return "seed"
        elif any(k in t for k in ("stem", "hypocotyl", "internode", "node")):
            return "stem"
        elif any(k in t for k in ("stress", "drought", "salt", "cold", "heat", "aba", "wound")):
            return "stress"
        elif any(k in t for k in ("guard", "stomata")):
            return "guard_cell"
        return "other"

    # ── Expressologs ───────────────────────────────────────────

    def get_expressologs(self, gene_id: str) -> ExpressologResult:
        """Find cross-species expressologs for a gene.

        Expressologs are genes with similar expression patterns across species,
        identified via orthology + co-expression conservation.
        """
        cache_key = f"bar_expressolog_{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached:
            if isinstance(cached, str):
                data = json.loads(cached)
                return ExpressologResult(
                    query_gene=data.get("query_gene", gene_id),
                    hits=[ExpressologHit(**h) for h in data.get("hits", [])],
                    note=data.get("note", ""),
                )
            return cached

        result = ExpressologResult(query_gene=gene_id)

        try:
            url = self.base_url + "get_expressologs.php"
            resp = self._session.get(
                url,
                params={"request": json.dumps([{"gene": gene_id}])},
                timeout=TIMEOUT,
            )

            if resp.ok:
                data = resp.json()
                hits = self._parse_expressologs(data, gene_id)
                result.hits = hits
                result.note = f"BAR expressologs: {len(hits)} hits"
            else:
                result.note = "BAR expressologs unavailable"

        except requests.Timeout:
            result.note = "BAR expressologs timed out"
        except Exception as e:
            result.note = f"BAR expressologs error: {e}"
            logger.warning("BAR expressologs failed for %s: %s", gene_id, e)

        # Cache
        try:
            cache_data = {
                "query_gene": result.query_gene,
                "hits": [h.__dict__ for h in result.hits],
                "note": result.note,
            }
            knowledge_cache.set(cache_key, json.dumps(cache_data), ttl=86400)
        except Exception:
            pass

        return result

    def _parse_expressologs(self, data: Any, gene_id: str) -> list[ExpressologHit]:
        """Parse BAR expressolog response."""
        hits: list[ExpressologHit] = []
        try:
            items = data if isinstance(data, list) else [data]
            for item in items:
                if isinstance(item, dict):
                    # Extract expressolog data
                    for key, value in item.items():
                        if key in ("gene", "agi", "request"):
                            continue
                        if isinstance(value, dict):
                            target_species = value.get("species", "")
                            corr = float(value.get("correlation", 0) or 0)
                            target_gene = value.get("gene", key)
                            hits.append(ExpressologHit(
                                query_gene=gene_id,
                                target_gene=str(target_gene),
                                target_species=str(target_species),
                                expression_correlation=corr,
                                data_source="BAR",
                            ))
                        elif isinstance(value, list):
                            for v in value:
                                if isinstance(v, dict):
                                    hits.append(ExpressologHit(
                                        query_gene=gene_id,
                                        target_gene=str(v.get("gene", "")),
                                        target_species=str(v.get("species", "")),
                                        expression_correlation=float(v.get("correlation", 0) or 0),
                                        data_source="BAR",
                                    ))
        except Exception as e:
            logger.debug("Failed to parse expressologs for %s: %s", gene_id, e)

        # Sort by correlation descending
        hits.sort(key=lambda x: x.expression_correlation, reverse=True)
        return hits[:10]  # top 10

    # ── PPI (protein-protein interactions) ────────────────────

    def get_interactions(self, gene_id: str) -> list[dict]:
        """Get protein-protein interactions for a gene from BAR."""
        cache_key = f"bar_ppi_{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached:
            return json.loads(cached) if isinstance(cached, str) else cached

        try:
            url = self.base_url + "aiv/get_interactions.php"
            resp = self._session.get(
                url,
                params={"request": json.dumps([{"agi": gene_id}])},
                timeout=TIMEOUT,
            )
            if resp.ok:
                data = resp.json()
                knowledge_cache.set(cache_key, json.dumps(data), ttl=86400 * 7)
                return data
        except Exception as e:
            logger.warning("BAR PPI failed for %s: %s", gene_id, e)

        return []

    # ── Convenience: full gene context ─────────────────────────

    def get_gene_context(self, gene_id: str) -> dict:
        """Get comprehensive gene context from public databases.

        Returns dict with:
          - annotation: gene annotation
          - tissue_expression: top expressed tissues
          - expressologs: cross-species expression homologs
          - ppi_count: number of known protein interactions
        """
        # Run independent queries in sequence (BAR rate-limits parallel requests)
        annotation = self.get_gene_annotation(gene_id)
        tissue_profile = self.get_tissue_expression(gene_id)
        expressologs = self.get_expressologs(gene_id)
        ppi = self.get_interactions(gene_id)

        return {
            "gene_id": gene_id,
            "symbol": annotation.get("symbol", "") if isinstance(annotation, dict) else "",
            "annotation": annotation.get("annotation", "") if isinstance(annotation, dict) else "",
            "top_tissues": tissue_profile.top_expressed_tissues,
            "tissue_count": len(tissue_profile.tissues),
            "expressolog_count": len(expressologs.hits),
            "expressologs": [
                {"gene": h.target_gene, "species": h.target_species, "correlation": h.expression_correlation}
                for h in expressologs.hits[:5]
            ],
            "ppi_count": len(ppi),
            "data_sources": ["BAR"],
            "note": tissue_profile.note,
        }


# ── Module-level singleton ─────────────────────────────────────

_bar_client: BARClient | None = None


def get_bar_client() -> BARClient:
    """Get or create BAR client singleton."""
    global _bar_client
    if _bar_client is None:
        _bar_client = BARClient()
    return _bar_client


# ── Convenience functions ──────────────────────────────────────

def get_arabidopsis_expression(gene_id: str) -> GeneExpressionProfile:
    """Get tissue expression profile for an Arabidopsis gene from BAR."""
    return get_bar_client().get_tissue_expression(gene_id)


def get_expressologs(gene_id: str) -> ExpressologResult:
    """Get cross-species expressologs for a gene from BAR."""
    return get_bar_client().get_expressologs(gene_id)


def get_gene_public_context(gene_id: str) -> dict:
    """Get comprehensive public database context for a gene."""
    return get_bar_client().get_gene_context(gene_id)
