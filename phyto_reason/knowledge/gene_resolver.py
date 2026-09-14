"""
gene_resolver.py — Universal gene ID resolution via MyGene.info + NCBI.

Multi-level dynamic resolution pipeline:
  Level 1: MyGene.info REST API (covers ~20,000 species, any ID format)
  Level 2: NCBI E-utilities (all species in NCBI Gene)
  Level 3: Pattern-based fallback

No hardcoded gene data — all results come from live database queries.
All queries cached for performance.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from typing import Any

import requests

from phyto_reason.utils.cache import knowledge_cache

logger = logging.getLogger("gene_resolver")

# ═══════════════════════════════════════════════════════════════
# Data model
# ═══════════════════════════════════════════════════════════════


@dataclass
class GeneInfo:
    """Resolved gene information from database queries."""
    input_id: str = ""
    primary_symbol: str = ""          # e.g. "OsNAC6" (preferred common name)
    synonyms: list[str] = field(default_factory=list)  # all known aliases
    description: str = ""             # functional description
    species: str = ""                 # e.g. "Oryza sativa"
    ncbi_gene_id: str = ""            # NCBI Gene ID (numeric string)
    taxid: str = ""                   # NCBI Taxonomy ID
    tf_family: str = ""               # e.g. "NAC" if inferrable
    source: str = ""                  # "mygene" | "ncbi" | "pattern"


# ═══════════════════════════════════════════════════════════════
# ID format detection
# ═══════════════════════════════════════════════════════════════

_RAPDB_PATTERN = re.compile(r'^Os\d{2}g\d{7}$', re.IGNORECASE)
_MSU_PATTERN = re.compile(r'^LOC_Os\d{2}g\d{5,6}$', re.IGNORECASE)
_TAIR_PATTERN = re.compile(r'^AT[1-5MC]G\d{5}(\.[0-9]+)?$', re.IGNORECASE)
_ENSEMBL_PLANT_PATTERN = re.compile(
    r'^(Zm|Solyc|Glyma\.|Bradi|GSVIV|Sobic|Traes|GRMZM|AC|evm\.)\w+$', re.IGNORECASE
)
_NCBI_GENE_PATTERN = re.compile(r'^\d{4,}$')


def detect_id_format(gene_id: str) -> str:
    """Detect the naming system of a gene ID.

    Returns one of: "rapdb", "msu", "tair", "ensembl_plant", "ncbi_gene", "symbol", "unknown"
    """
    if not gene_id or not gene_id.strip():
        return "unknown"
    gene_id = gene_id.strip()
    if _RAPDB_PATTERN.match(gene_id):
        return "rapdb"
    if _MSU_PATTERN.match(gene_id):
        return "msu"
    if _TAIR_PATTERN.match(gene_id):
        return "tair"
    if _NCBI_GENE_PATTERN.match(gene_id):
        return "ncbi_gene"
    if _ENSEMBL_PLANT_PATTERN.match(gene_id):
        return "ensembl_plant"
    # Check for gene symbol patterns (e.g., OsNAC6, WRKY45, MYB75)
    if re.search(r'^(Os|AT|GRMZM|LOC_|Zm)', gene_id, re.IGNORECASE):
        return "symbol"
    return "symbol"


# ═══════════════════════════════════════════════════════════════
# MyGene.info client
# ═══════════════════════════════════════════════════════════════


class MyGeneClient:
    """Client for MyGene.info REST API (https://mygene.info).

    MyGene.info is a free, public gene annotation service covering
    ~20,000 species. It accepts any gene ID format (RAP-DB, TAIR, Ensembl,
    NCBI Gene ID, gene symbol) and returns normalized gene information.

    Rate limit: ~10 queries/second (generous for our use case).
    """

    BASE_URL = "https://mygene.info/v3"
    TIMEOUT = 10

    def __init__(self, timeout: int = 10):
        self.timeout = timeout
        self._session = requests.Session()
        self._session.headers.update({
            "User-Agent": "phyto_reason/4.4",
            "Accept": "application/json",
        })

    def query(self, gene_id: str) -> dict | None:
        """Query MyGene.info for a gene ID.

        Returns the top hit as a dict, or None if not found.
        """
        if not gene_id or not gene_id.strip():
            return None

        gene_id = gene_id.strip()

        # Check cache
        cache_key = f"mygene:query:{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached if cached != "__none__" else None

        try:
            params = {
                "q": gene_id,
                "fields": "symbol,name,taxid,alias,summary,entrezgene,other_names",
            }
            resp = self._session.get(
                f"{self.BASE_URL}/query",
                params=params,
                timeout=self.timeout,
            )

            if resp.status_code == 429:
                time.sleep(1)
                resp = self._session.get(
                    f"{self.BASE_URL}/query",
                    params=params,
                    timeout=self.timeout,
                )

            if not resp.ok:
                logger.debug("MyGene.info query failed for %s: HTTP %s", gene_id, resp.status_code)
                knowledge_cache.set(cache_key, "__none__", ttl=3600)
                return None

            data = resp.json()
            hits = data.get("hits", [])
            if not hits:
                knowledge_cache.set(cache_key, "__none__", ttl=86400)
                return None

            hit = hits[0]
            result = {
                "_id": str(hit.get("_id", "")),
                "symbol": hit.get("symbol", ""),
                "name": hit.get("name", ""),
                "taxid": str(hit.get("taxid", "")),
                "alias": hit.get("alias", []),
                "summary": hit.get("summary", ""),
                "entrezgene": str(hit.get("entrezgene", hit.get("_id", ""))),
                "other_names": hit.get("other_names", []),
            }

            knowledge_cache.set(cache_key, result, ttl=86400)
            logger.info("MyGene.info hit: %s → %s (NCBI:%s)", gene_id, result["symbol"], result["entrezgene"])
            return result

        except requests.Timeout:
            logger.debug("MyGene.info timeout for %s", gene_id)
            return None
        except Exception as e:
            logger.debug("MyGene.info error for %s: %s", gene_id, e)
            return None

    def get_by_ncbi_id(self, ncbi_gene_id: str) -> dict | None:
        """Get gene info by NCBI Gene ID."""
        if not ncbi_gene_id:
            return None

        cache_key = f"mygene:gene:{ncbi_gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached if cached != "__none__" else None

        try:
            resp = self._session.get(
                f"{self.BASE_URL}/gene/{ncbi_gene_id}",
                params={"fields": "symbol,name,taxid,alias,summary,other_names"},
                timeout=self.timeout,
            )

            if not resp.ok:
                knowledge_cache.set(cache_key, "__none__", ttl=3600)
                return None

            data = resp.json()
            result = {
                "_id": str(data.get("_id", ncbi_gene_id)),
                "symbol": data.get("symbol", ""),
                "name": data.get("name", ""),
                "taxid": str(data.get("taxid", "")),
                "alias": data.get("alias", []),
                "summary": data.get("summary", ""),
                "entrezgene": str(ncbi_gene_id),
                "other_names": data.get("other_names", []),
            }

            knowledge_cache.set(cache_key, result, ttl=86400)
            return result

        except Exception as e:
            logger.debug("MyGene.info gene/%s error: %s", ncbi_gene_id, e)
            return None


# ═══════════════════════════════════════════════════════════════
# NCBI E-utilities helpers (gene search + summary)
# ═══════════════════════════════════════════════════════════════


def _ncbi_esearch(query: str, retmax: int = 3) -> list[str]:
    """Search NCBI Gene database, return list of NCBI Gene IDs."""
    cache_key = f"gene_resolver:ncbi:search:{query}"
    cached = knowledge_cache.get(cache_key)
    if cached is not None:
        return cached

    try:
        resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={
                "db": "gene",
                "term": query,
                "retmax": retmax,
                "retmode": "json",
                "tool": "phyto_reason",
                "email": "plantomics@example.com",
            },
            timeout=10,
        )
        data = resp.json()
        ids = [str(i) for i in data.get("esearchresult", {}).get("idlist", [])]
        knowledge_cache.set(cache_key, ids, ttl=86400)
        return ids
    except Exception as e:
        logger.debug("NCBI esearch failed: %s", e)
        return []


def _ncbi_esummary(ncbi_gene_id: str) -> dict | None:
    """Get gene summary from NCBI esummary."""
    cache_key = f"gene_resolver:ncbi:summary:{ncbi_gene_id}"
    cached = knowledge_cache.get(cache_key)
    if cached is not None:
        return cached if cached != "__none__" else None

    try:
        resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={
                "db": "gene",
                "id": ncbi_gene_id,
                "retmode": "json",
                "tool": "phyto_reason",
                "email": "plantomics@example.com",
            },
            timeout=10,
        )
        data = resp.json()
        gene_record = data.get("result", {}).get(str(ncbi_gene_id), {})
        if not gene_record or "error" in gene_record:
            knowledge_cache.set(cache_key, "__none__", ttl=86400)
            return None

        result = {
            "name": gene_record.get("name", ""),
            "description": gene_record.get("description", ""),
            "taxid": str(gene_record.get("taxid", "")),
        }
        knowledge_cache.set(cache_key, result, ttl=86400)
        return result
    except Exception as e:
        logger.debug("NCBI esummary failed for %s: %s", ncbi_gene_id, e)
        return None


# ═══════════════════════════════════════════════════════════════
# TF family inference
# ═══════════════════════════════════════════════════════════════

_TF_FAMILIES = [
    "NAC", "WRKY", "MYB", "BHLH", "ERF", "BZIP", "WD40", "AP2",
    "DREB", "HD-ZIP", "TCP", "GRAS", "SPL", "ARF", "GATA", "DOF", "NF-Y",
    "MADS", "HSF", "LOB", "C2H2", "Trihelix",
]

# TF patterns that can appear as substrings in gene symbols
# Use case-insensitive substring match with boundary check
# (e.g., OsbHLH007, NAC100_7, OsMADS23, OsWRKY45)
_TF_FAMILIES_UPPER = [fam.upper() for fam in _TF_FAMILIES]


def _infer_tf_family(symbol: str, description: str, aliases: list[str]) -> str:
    """Infer TF family from gene symbol, description, and aliases."""
    raw_text = f"{symbol} {' '.join(aliases)} {description}"
    upper_text = raw_text.upper()

    # Sort by length descending so longer names match first (HD-ZIP before ZIP)
    indexed_fams = sorted(enumerate(_TF_FAMILIES), key=lambda x: len(x[1]), reverse=True)

    for orig_idx, fam in indexed_fams:
        fam_upper = fam.upper()
        idx = upper_text.find(fam_upper)
        if idx < 0:
            continue

        # Boundary check on original-case text.
        # The family name must be preceded by one of:
        #   - Start of text
        #   - Non-alpha char (space, digit, underscore, dash)
        #   - Lowercase letter (camelCase: OsNAC6, OsbHLH007)
        #   - Up to 2 uppercase prefix letters followed by non-alpha or lowercase
        #     (species prefixes: AT→MYB75, Os→WRKY45, Zm→NAC100)
        before_ok = False
        if idx == 0:
            before_ok = True
        else:
            for lookback in range(1, min(4, idx + 1)):
                ch = raw_text[idx - lookback]
                if not ch.isalpha() or ch.islower():
                    before_ok = True
                    break

        after_ok = (
            (idx + len(fam) >= len(raw_text))
            or not raw_text[idx + len(fam)].isalpha()
        )
        if before_ok and after_ok:
            return fam

    return ""


# ═══════════════════════════════════════════════════════════════
# Species → taxid mapping (for NCBI fallback)
# ═══════════════════════════════════════════════════════════════

_PLANT_TAXIDS: dict[str, str] = {
    "arabidopsis thaliana": "3702",
    "arabidopsis": "3702",
    "oryza sativa": "4530",
    "rice": "4530",
    "nicotiana tabacum": "4097",
    "tobacco": "4097",
    "zea mays": "4577",
    "maize": "4577",
    "corn": "4577",
    "solanum lycopersicum": "4081",
    "tomato": "4081",
    "vitis vinifera": "29760",
    "grape": "29760",
    "glycine max": "3847",
    "soybean": "3847",
}


def _resolve_taxid(species: str) -> str:
    """Map species name to NCBI Taxonomy ID."""
    if not species:
        return ""
    species_lower = species.lower().strip()
    if species_lower in _PLANT_TAXIDS:
        return _PLANT_TAXIDS[species_lower]
    for known, taxid in _PLANT_TAXIDS.items():
        if species_lower in known or known in species_lower:
            return taxid
    return ""


# ═══════════════════════════════════════════════════════════════
# GeneResolver — 3-level dynamic resolution pipeline
# ═══════════════════════════════════════════════════════════════


class GeneResolver:
    """Universal gene ID resolver using MyGene.info + NCBI + pattern matching.

    All gene information comes from live database queries — no hardcoded gene data.

    Usage:
        resolver = GeneResolver()
        info = resolver.resolve("Os01g0884300", species="rice")
        # → GeneInfo(symbol="OsNAC6", ncbi_gene_id="4325006", tf_family="NAC", ...)

        info = resolver.resolve("AT1G56650", species="arabidopsis")
        # → GeneInfo(symbol="PAP1", ncbi_gene_id="842120", tf_family="MYB", ...)
    """

    def __init__(self):
        self._mygene = MyGeneClient()

    # ═════════════════════════════════════════════════════════
    # Public API
    # ═════════════════════════════════════════════════════════

    def resolve(self, gene_id: str, species: str = "") -> GeneInfo:
        """Resolve a gene ID using the 3-level pipeline.

        Args:
            gene_id: Any gene identifier.
            species: Optional species hint.

        Returns:
            GeneInfo with symbol, synonyms, NCBI Gene ID, and TF family if known.
        """
        if not gene_id or not gene_id.strip():
            return GeneInfo(input_id=gene_id or "", source="pattern")

        gene_id = gene_id.strip()
        taxid_hint = _resolve_taxid(species)

        # ── Level 1: MyGene.info ──────────────────────────
        result = self._resolve_via_mygene(gene_id)
        if result:
            return result

        # ── Level 2: NCBI E-utilities ─────────────────────
        result = self._resolve_via_ncbi(gene_id, taxid_hint)
        if result:
            return result

        # Try MyGene.info again with species-adjusted search
        # (MyGene might need different query format for some IDs)
        if species:
            # Try searching by species-adjusted query
            alt_query = f"{gene_id} {species}"
            result = self._resolve_via_mygene(alt_query)
            if result and result.primary_symbol != alt_query:
                return result

        # ── Level 3: Pattern-based fallback ───────────────
        fmt = detect_id_format(gene_id)
        return GeneInfo(
            input_id=gene_id,
            primary_symbol=gene_id,
            synonyms=[gene_id],
            description="",
            species=species,
            source="pattern",
        )

    def expand_search_query(self, gene_id: str, species: str = "") -> str:
        """Build an expanded PubMed/KEGG search query with all known aliases.

        Example: "Os01g0884300" → "(Os01g0884300 OR OsNAC6 OR SNAC2 OR NAC6)"
        """
        info = self.resolve(gene_id, species)
        # Collect all unique names: input + symbol + synonyms
        all_names = [info.input_id, info.primary_symbol] + info.synonyms
        # Deduplicate, remove numeric-only entries, limit to 5
        seen = set()
        unique = []
        for name in all_names:
            name = name.strip()
            if name and name not in seen and not name.isdigit():
                seen.add(name)
                unique.append(name)
                if len(unique) >= 5:
                    break

        if len(unique) <= 1:
            return gene_id
        return "(" + " OR ".join(unique) + ")"

    def get_tf_family(self, gene_id: str, species: str = "") -> str:
        """Get the TF family for a gene if known."""
        info = self.resolve(gene_id, species)
        return info.tf_family

    # ═════════════════════════════════════════════════════════
    # Level 1: MyGene.info
    # ═════════════════════════════════════════════════════════

    def _resolve_via_mygene(self, gene_id: str) -> GeneInfo | None:
        """Resolve gene via MyGene.info REST API."""
        hit = self._mygene.query(gene_id)
        if not hit:
            return None

        symbol = hit.get("symbol", "")

        # Reject garbage results: single-char symbols or obviously wrong mappings
        if symbol and len(symbol.strip()) <= 1 and not symbol.isdigit():
            logger.debug("MyGene.info returned garbage symbol '%s' for %s, rejecting", symbol, gene_id)
            return None

        all_aliases: list[str] = list(hit.get("alias", []))
        other_names: list[str] = list(hit.get("other_names", []))
        name = hit.get("name", "")
        ncbi_id = hit.get("entrezgene", hit.get("_id", ""))
        taxid = hit.get("taxid", "")
        summary = hit.get("summary", "")

        # Build comprehensive synonym list
        synonyms = [gene_id]
        if symbol and symbol != gene_id:
            synonyms.append(symbol)
        for alias in all_aliases:
            if alias and alias not in synonyms:
                synonyms.append(alias)
        for oname in other_names:
            if oname and oname not in synonyms:
                synonyms.append(oname)

        # Prefer the most "human-readable" alias as primary symbol
        primary = symbol
        if not primary or primary.startswith("LOC"):
            # Score each alias for desirability as primary name
            def _alias_score(alias: str) -> int:
                """Higher score = better primary symbol candidate.
                Returns -1000 for entries that are clearly not gene symbols."""
                # Reject entries that are clearly descriptions, not gene symbols
                if len(alias) > 30:
                    return -1000  # descriptions are long
                if " " in alias:
                    return -1000  # gene symbols don't have spaces
                if alias.lower().startswith("domain-containing"):
                    return -1000  # part of description
                if alias.lower() in ("protein", "putative", "uncharacterized", "predicted"):
                    return -1000

                s = 0
                # Penalize LOC-prefixed names (locus tags)
                if alias.startswith("LOC") or alias.startswith("OsJ_"):
                    s -= 10
                # Penalize names with many digits (locus-like)
                digits = sum(c.isdigit() for c in alias)
                if digits > 3:
                    s -= 5
                # Bonus for containing letters (gene-like, not pure numbers)
                if re.search(r'[A-Za-z]', alias):
                    s += 3
                # Bonus for known organism prefixes
                if re.match(r'^(Os|AT|Zm|GRMZM|Solyc|Glyma)', alias):
                    s += 2
                # Big bonus for known TF family names in the alias
                if _infer_tf_family(alias, "", []):
                    s += 5
                # Prefer shorter names (easier to read)
                if 3 <= len(alias) <= 12:
                    s += 2
                return s

            best_alias = None
            best_score = -100
            for alias in all_aliases + other_names:
                score = _alias_score(alias)
                if score > best_score:
                    best_score = score
                    best_alias = alias
            if best_alias:
                primary = best_alias

        description = summary or name or ""

        # Infer TF family
        tf_family = _infer_tf_family(symbol, description, all_aliases + other_names)

        return GeneInfo(
            input_id=gene_id,
            primary_symbol=primary or symbol,
            synonyms=synonyms,
            description=description,
            species="",
            ncbi_gene_id=ncbi_id,
            taxid=taxid,
            tf_family=tf_family,
            source="mygene",
        )

    # ═════════════════════════════════════════════════════════
    # Level 2: NCBI E-utilities
    # ═════════════════════════════════════════════════════════

    def _resolve_via_ncbi(self, gene_id: str, taxid: str = "") -> GeneInfo | None:
        """Resolve gene via NCBI E-utilities as fallback."""
        # Build search query
        if taxid:
            query = f"{gene_id}[All Fields] AND txid{taxid}[Organism]"
        else:
            query = f"{gene_id}[All Fields]"

        ids = _ncbi_esearch(query)
        if not ids and taxid:
            # Try without organism constraint
            ids = _ncbi_esearch(f"{gene_id}[All Fields]")

        if not ids:
            return None

        ncbi_id = ids[0]
        summary = _ncbi_esummary(ncbi_id)
        if not summary:
            return None

        symbol = summary.get("name", gene_id)
        desc = summary.get("description", "")
        ncbi_taxid = summary.get("taxid", taxid)

        # Try MyGene.info for richer alias data using the NCBI Gene ID
        mygene_hit = self._mygene.get_by_ncbi_id(ncbi_id)
        if mygene_hit:
            aliases: list[str] = [gene_id, symbol]
            for a in mygene_hit.get("alias", []):
                if a and a not in aliases:
                    aliases.append(a)
            for a in mygene_hit.get("other_names", []):
                if a and a not in aliases:
                    aliases.append(a)
            mg_desc = mygene_hit.get("summary", "") or mygene_hit.get("name", "")
            if mg_desc:
                desc = mg_desc
            tf_family = _infer_tf_family(
                symbol, desc, mygene_hit.get("alias", []) + mygene_hit.get("other_names", []),
            )
        else:
            aliases = [gene_id, symbol]
            tf_family = _infer_tf_family(symbol, desc, [])

        return GeneInfo(
            input_id=gene_id,
            primary_symbol=symbol,
            synonyms=aliases,
            description=desc,
            species="",
            ncbi_gene_id=ncbi_id,
            taxid=ncbi_taxid,
            tf_family=tf_family,
            source="ncbi",
        )


# ═══════════════════════════════════════════════════════════════
# Module-level convenience functions
# ═══════════════════════════════════════════════════════════════

_default_resolver: GeneResolver | None = None


def get_gene_resolver() -> GeneResolver:
    """Get or create the default GeneResolver instance."""
    global _default_resolver
    if _default_resolver is None:
        _default_resolver = GeneResolver()
    return _default_resolver


def resolve_gene(gene_id: str, species: str = "") -> GeneInfo:
    """Convenience function to resolve a gene ID using live database queries."""
    return get_gene_resolver().resolve(gene_id, species)


def expand_gene_query(gene_id: str, species: str = "") -> str:
    """Build an expanded search query with all known gene aliases."""
    return get_gene_resolver().expand_search_query(gene_id, species)
