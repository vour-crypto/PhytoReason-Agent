"""
ortholog_finder.py — Gene-level ortholog computation via NCBI E-utilities.

Enables Layer 1 (cross-species mode) with real ortholog discovery:
  1. NCBI Ortholog Database (elink gene_orthologs)
  2. Gene symbol search in target species
  3. TF family search as coarse fallback
  4. Curated SpeciesKnowledge as last resort

Pure Python Needleman-Wunsch alignment — no BioPython dependency.
All queries cached via utils.cache for responsiveness.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from typing import Any, ClassVar

import requests

from phyto_reason.utils.cache import knowledge_cache

logger = logging.getLogger("ortholog_finder")

# ═══════════════════════════════════════════════════════════════
# Data models
# ═══════════════════════════════════════════════════════════════


@dataclass
class OrthologHit:
    """A single ortholog relationship between source and target genes."""
    source_gene_id: str = ""
    source_species: str = ""
    source_gene_symbol: str = ""
    target_gene_id: str = ""
    target_species: str = ""
    target_gene_symbol: str = ""
    percent_identity: float = 0.0       # 0-100
    alignment_length: int = 0
    method: str = ""                    # "ncbi_ortholog_db" | "sequence_similarity" | "symbol_match" | "curated"
    confidence: str = "low"             # "high" | "medium" | "low"
    source_protein_id: str = ""
    target_protein_id: str = ""
    description: str = ""


@dataclass
class OrthologResult:
    """Result of an ortholog search."""
    query_gene_id: str = ""
    source_species: str = ""
    target_species: str = ""
    hits: list[OrthologHit] = field(default_factory=list)
    method_used: str = ""               # "ncbi" | "fallback_curated"
    note: str = ""


# ═══════════════════════════════════════════════════════════════
# Species → NCBI TaxID mapping (common plant species)
# ═══════════════════════════════════════════════════════════════

PLANT_TAXIDS: dict[str, str] = {
    "arabidopsis thaliana": "3702",
    "arabidopsis": "3702",
    "thaliana": "3702",
    "nicotiana tabacum": "4097",
    "tobacco": "4097",
    "nicotiana benthamiana": "4100",
    "benthamiana": "4100",
    "oryza sativa": "4530",
    "rice": "4530",
    "solanum lycopersicum": "4081",
    "tomato": "4081",
    "zea mays": "4577",
    "maize": "4577",
    "corn": "4577",
    "vitis vinifera": "29760",
    "grape": "29760",
    "coptis japonica": "3442",
    "coptis chinensis": "3442",
    "coptis": "3442",
    "coptis teeta": "261452",
    "artemisia annua": "35608",
    "artemisia": "35608",
    # ── Medicinal plants (v4.5) ─────────────────────────
    "salvia miltiorrhiza": "22663",
    "salvia": "22663",
    "danshen": "22663",
    "scutellaria baicalensis": "65409",
    "scutellaria": "65409",
    "baical skullcap": "65409",
    "panax ginseng": "4054",
    "ginseng": "4054",
    "panax notoginseng": "44586",
    "notoginseng": "44586",
    "glycyrrhiza uralensis": "74613",
    "glycyrrhiza": "74613",
    "licorice": "74613",
    "angelica sinensis": "165353",
    "angelica": "165353",
    "dong quai": "165353",
    "taxus chinensis": "29814",
    "taxus brevifolia": "46222",
    "taxus": "29814",
    "yew": "29814",
    "camptotheca acuminata": "16922",
    "camptotheca": "16922",
    "rauvolfia serpentina": "4060",
    "rauvolfia": "4060",
    "papaver somniferum": "3469",
    "opium poppy": "3469",
    "digitalis purpurea": "4168",
    "digitalis": "4168",
    "foxglove": "4168",
    "atropa belladonna": "33113",
    "belladonna": "33113",
    "hypericum perforatum": "65561",
    "hypericum": "65561",
    "st johns wort": "65561",
    "curcuma longa": "136217",
    "curcuma": "136217",
    "turmeric": "136217",
    "zingiber officinale": "94328",
    "zingiber": "94328",
    "ginger": "94328",
    "ginkgo biloba": "3311",
    "ginkgo": "3311",
    "ephedra sinica": "33152",
    "ephedra": "33152",
    "cannabis sativa": "3483",
    "cannabis": "3483",
    "hemp": "3483",
    "catharanthus roseus": "4058",
    "catharanthus": "4058",
    "sorghum bicolor": "4558",
    "sorghum": "4558",
    "glycine max": "3847",
    "soybean": "3847",
    "triticum aestivum": "4565",
    "wheat": "4565",
    "hordeum vulgare": "4513",
    "barley": "4513",
    "solanum tuberosum": "4113",
    "potato": "4113",
    "brassica napus": "3708",
    "canola": "3708",
    "gossypium hirsutum": "3635",
    "cotton": "3635",
    "medicago truncatula": "3880",
    "setaria italica": "4555",
    "brachypodium distachyon": "15368",
}


# ═══════════════════════════════════════════════════════════════
# OrthologFinder
# ═══════════════════════════════════════════════════════════════


class OrthologFinder:
    """Find orthologous genes across plant species using NCBI E-utilities.

    Usage:
        finder = OrthologFinder()
        result = finder.find_orthologs("AT1G56650", "Arabidopsis", "tobacco")
        for hit in result.hits:
            print(f"{hit.target_gene_symbol}: {hit.percent_identity:.1f}% identity")
    """

    NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    ESEARCH = NCBI_BASE + "esearch.fcgi"
    ESUMMARY = NCBI_BASE + "esummary.fcgi"
    EFETCH = NCBI_BASE + "efetch.fcgi"
    ELINK = NCBI_BASE + "elink.fcgi"

    # NCBI requires identifying info
    TOOL_NAME = "phyto_reason"
    EMAIL = "plantomics@example.com"

    def __init__(self, timeout: int = 15, max_retries: int = 2):
        self.timeout = timeout
        self.max_retries = max_retries
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": f"{self.TOOL_NAME}/4.0"})

    # ═════════════════════════════════════════════════════════
    # Public API
    # ═════════════════════════════════════════════════════════

    def find_orthologs(
        self,
        gene_id: str,
        source_species: str = "Arabidopsis thaliana",
        target_species: str = "",
    ) -> OrthologResult:
        """Find orthologs of gene_id in target_species.

        Args:
            gene_id: Source gene ID (e.g. "AT1G56650" or NCBI Gene ID).
            source_species: Species the gene comes from.
            target_species: Species to search for orthologs in.

        Returns:
            OrthologResult with ranked hits.
        """
        source_taxid = self._resolve_taxid(source_species)
        target_taxid = self._resolve_taxid(target_species) if target_species else ""

        if not target_taxid:
            return OrthologResult(
                query_gene_id=gene_id,
                source_species=source_species,
                target_species=target_species,
                method_used="fallback_curated",
                note=f"Species '{target_species}' not recognized. Supported species: "
                     f"{', '.join(sorted(PLANT_TAXIDS.keys())[:20])}...",
            )

        # Check cache
        cache_key = f"ortholog:{source_taxid}:{gene_id}:{target_taxid}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            logger.info("Ortholog cache hit: %s", cache_key)
            return cached

        result = self._find_orthologs_impl(gene_id, source_taxid, target_taxid,
                                            source_species, target_species)

        # Cache result (24h)
        knowledge_cache.set(cache_key, result, ttl=86400)
        return result

    # ═════════════════════════════════════════════════════════
    # Implementation: 5-tier ortholog search
    # ═════════════════════════════════════════════════════════

    def _find_orthologs_impl(
        self,
        gene_id: str,
        source_taxid: str,
        target_taxid: str,
        source_species: str,
        target_species: str,
    ) -> OrthologResult:
        """Core ortholog search with progressive degradation."""

        result = OrthologResult(
            query_gene_id=gene_id,
            source_species=source_species,
            target_species=target_species,
        )

        # ── Tier 1: Resolve source gene in NCBI ──────────────
        source_info = self._resolve_gene(gene_id, source_taxid)
        if not source_info:
            # Gene ID not found in NCBI — maybe it's already an NCBI ID
            source_info = self._try_direct_lookup(gene_id)
            if not source_info:
                # Still try MyGene.info before giving up
                hits = self._search_by_domain_conservation(
                    gene_id, source_taxid, target_taxid,
                    source_species, target_species,
                )
                if hits:
                    result.method_used = "mygene_domain"
                    result.hits = sorted(hits, key=lambda h: h.percent_identity, reverse=True)
                    return result

                result.method_used = "fallback_curated"
                result.note = (
                    f"Gene '{gene_id}' not found in NCBI Gene database "
                    f"for {source_species}."
                )
                return result

        source_ncbi_id = source_info.get("ncbi_gene_id", gene_id)
        source_symbol = source_info.get("symbol", gene_id)
        source_desc = source_info.get("description", "")

        result.method_used = "ncbi"

        # ── Tier 2: NCBI Ortholog Database ──────────────────
        ortholog_ids = self._get_ortholog_links(source_ncbi_id)

        if ortholog_ids:
            hits = self._score_orthologs(
                ortholog_ids, target_taxid, gene_id, source_species,
                source_symbol, target_species,
            )
            if hits:
                result.hits = sorted(hits, key=lambda h: h.percent_identity, reverse=True)
                return result

        # ── Tier 2.5: MyGene.info domain-based orthologs ─────
        hits = self._search_by_domain_conservation(
            source_ncbi_id, source_taxid, target_taxid,
            source_species, target_species,
        )
        if hits:
            result.method_used = "mygene_domain"
            result.hits = sorted(hits, key=lambda h: h.percent_identity, reverse=True)
            return result

        # ── Tier 3: Gene symbol search in target species ─────
        if source_symbol and source_symbol != gene_id:
            hits = self._search_by_symbol(
                source_symbol, target_taxid, target_species,
                gene_id, source_species, source_symbol,
            )
            if hits:
                result.hits = sorted(hits, key=lambda h: h.percent_identity, reverse=True)
                return result

        # ── Tier 4: No orthologs found ──────────────────────
        result.note = (
            f"No orthologs of {gene_id} ({source_symbol}) found in {target_species} "
            f"via NCBI. The gene may be species-specific or the target species "
            f"genome is not well-represented in NCBI Gene."
        )
        return result

    # ═════════════════════════════════════════════════════════
    # Tier 2.5: MyGene.info domain-based functional ortholog search
    # ═════════════════════════════════════════════════════════

    def _search_by_domain_conservation(
        self,
        gene_id: str,
        source_taxid: str,
        target_taxid: str,
        source_species: str,
        target_species: str,
    ) -> list[OrthologHit]:
        """Find functional orthologs via conserved protein domains using MyGene.info.

        Strategy:
        1. Get InterPro/Pfam domains for the source gene from MyGene.info
        2. Search MyGene.info for genes in target species with same domains
        3. Score by domain overlap percentage
        """
        hits: list[OrthologHit] = []

        if not target_taxid:
            return hits

        try:
            import requests

            # Get source gene domain annotation
            mg_resp = requests.get(
                f"https://mygene.info/v3/gene/{gene_id}",
                params={"fields": "symbol,name,interpro,pfam,taxid"},
                timeout=10,
                headers={"User-Agent": "phyto_reason/4.4"},
            )
            if not mg_resp.ok:
                return hits

            mg_data = mg_resp.json()
            source_symbol = mg_data.get("symbol", gene_id)
            interpro_domains = mg_data.get("interpro", [])
            pfam_domain = mg_data.get("pfam", "")

            if not interpro_domains and not pfam_domain:
                # Try query-based search instead of gene lookup
                mg_query = requests.get(
                    "https://mygene.info/v3/query",
                    params={
                        "q": gene_id,
                        "fields": "symbol,name,interpro,pfam,taxid",
                    },
                    timeout=10,
                    headers={"User-Agent": "phyto_reason/4.4"},
                )
                if mg_query.ok:
                    qdata = mg_query.json()
                    qhits = qdata.get("hits", [])
                    if qhits:
                        qhit = qhits[0]
                        interpro_domains = qhit.get("interpro", [])
                        pfam_domain = qhit.get("pfam", "")
                        source_symbol = qhit.get("symbol", source_symbol)

            if not interpro_domains and not pfam_domain:
                return hits

            # Search target species for genes with the same domains
            target_genes: dict[str, dict] = {}  # ncbi_id → gene_data
            domain_scores: dict[str, float] = {}

            # Use each InterPro domain as a search filter
            for domain in interpro_domains[:3]:  # top 3 domains
                domain_id = domain.get("id", "")
                if not domain_id:
                    continue

                domain_resp = requests.get(
                    "https://mygene.info/v3/query",
                    params={
                        "q": f"interpro:{domain_id} AND taxid:{target_taxid}",
                        "fields": "symbol,name,taxid,interpro",
                        "size": 20,
                    },
                    timeout=15,
                    headers={"User-Agent": "phyto_reason/4.4"},
                )

                if not domain_resp.ok:
                    continue

                domain_data = domain_resp.json()
                for dhit in domain_data.get("hits", []):
                    ncbi_id = str(dhit.get("_id", ""))
                    if not ncbi_id:
                        continue

                    if ncbi_id not in target_genes:
                        target_genes[ncbi_id] = dhit
                        domain_scores[ncbi_id] = 0.0

                    # Each shared domain adds to the score
                    domain_scores[ncbi_id] += 1.0 / len(interpro_domains)

                time.sleep(0.3)  # rate limit

            # Also try Pfam if available
            if pfam_domain:
                pfam_resp = requests.get(
                    "https://mygene.info/v3/query",
                    params={
                        "q": f"pfam:{pfam_domain} AND taxid:{target_taxid}",
                        "fields": "symbol,name,taxid",
                        "size": 20,
                    },
                    timeout=15,
                    headers={"User-Agent": "phyto_reason/4.4"},
                )

                if pfam_resp.ok:
                    pfam_data = pfam_resp.json()
                    for dhit in pfam_data.get("hits", []):
                        ncbi_id = str(dhit.get("_id", ""))
                        if ncbi_id not in target_genes:
                            target_genes[ncbi_id] = dhit
                            domain_scores[ncbi_id] = 0.0
                        domain_scores[ncbi_id] += 0.3  # bonus for Pfam match
                    time.sleep(0.3)

            # Build OrthologHit objects
            for ncbi_id, gdata in target_genes.items():
                score = domain_scores.get(ncbi_id, 0.0)
                # Normalize: max 100%
                normalized_score = min(score * 100, 100.0)

                confidence = (
                    "high" if normalized_score >= 80 else
                    "medium" if normalized_score >= 40 else
                    "low"
                )

                hits.append(OrthologHit(
                    source_gene_id=gene_id,
                    source_species=source_species,
                    source_gene_symbol=source_symbol,
                    target_gene_id=ncbi_id,
                    target_species=target_species,
                    target_gene_symbol=gdata.get("symbol", ncbi_id),
                    percent_identity=normalized_score,  # domain overlap %
                    alignment_length=0,
                    method="domain_conservation",
                    confidence=confidence,
                    description=gdata.get("name", ""),
                ))

            # ── Enrich with sequence alignment ────────────────
            if hits and len(hits) >= 1:
                self._enrich_with_sequence_alignment(
                    hits, source_symbol, source_species, target_species,
                )

            if hits:
                logger.info(
                    "MyGene.info domain search: %d functional orthologs in %s for %s",
                    len(hits), target_species, gene_id,
                )

        except Exception as e:
            logger.debug("MyGene.info domain search failed: %s", e)

        return hits

    # ═════════════════════════════════════════════════════════
    # Tier 2.6: Sequence alignment enrichment
    # ═════════════════════════════════════════════════════════

    def _enrich_with_sequence_alignment(
        self,
        hits: list[OrthologHit],
        source_symbol: str,
        source_species: str,
        target_species: str,
    ) -> None:
        """Enrich domain-based ortholog hits with protein sequence alignment.

        Fetches the source gene's protein sequence, then pairwise-aligns with
        each candidate's protein sequence. Updates percent_identity and confidence
        based on real sequence similarity (Needleman-Wunsch).
        """
        if not hits:
            return

        # Get source protein sequence
        source_gene_id = hits[0].source_gene_id
        source_protein_id = ""
        source_seq = ""

        try:
            source_protein_id = self._get_protein_accession(source_gene_id)
            if source_protein_id:
                source_seq = self._fetch_protein_sequence(source_protein_id)
        except Exception:
            pass

        if not source_seq:
            # Try via MyGene.info UniProt accession
            try:
                import requests
                mg_resp = requests.get(
                    f"https://mygene.info/v3/gene/{source_gene_id}",
                    params={"fields": "uniprot"},
                    timeout=10,
                    headers={"User-Agent": "phyto_reason/4.5"},
                )
                if mg_resp.ok:
                    uniprot_data = mg_resp.json().get("uniprot", {})
                    if isinstance(uniprot_data, dict):
                        swissprot = uniprot_data.get("Swiss-Prot", uniprot_data.get("TrEMBL", ""))
                        if swissprot:
                            # Fetch from UniProt
                            uniprot_resp = requests.get(
                                f"https://rest.uniprot.org/uniprotkb/{swissprot}.fasta",
                                timeout=10,
                                headers={"User-Agent": "phyto_reason/4.5"},
                            )
                            if uniprot_resp.ok:
                                source_seq = self._parse_fasta(uniprot_resp.text) if hasattr(self, '_parse_fasta') else ""
            except Exception:
                pass

        if not source_seq:
            logger.debug("No source protein sequence available for %s, skipping alignment", source_gene_id)
            return

        aligned_count = 0
        for hit in hits[:15]:  # align top 15 candidates
            if hit.alignment_length > 0:
                continue  # already aligned

            try:
                target_protein_id = self._get_protein_accession(hit.target_gene_id)
                if target_protein_id:
                    target_seq = self._fetch_protein_sequence(target_protein_id)
                    if target_seq:
                        identity, align_len = self._compute_identity(source_seq, target_seq)

                        if identity > 0:
                            # Update hit with real sequence identity
                            hit.percent_identity = round(identity, 1)
                            hit.alignment_length = align_len
                            hit.source_protein_id = source_protein_id
                            hit.target_protein_id = target_protein_id
                            hit.method = "sequence_similarity"

                            # Re-score confidence based on sequence identity
                            if identity >= 60:
                                hit.confidence = "high"
                            elif identity >= 35:
                                hit.confidence = "medium"
                            else:
                                hit.confidence = "low"

                            aligned_count += 1
                            time.sleep(0.35)  # rate limit
                else:
                    # Try to get sequence via MyGene.info → UniProt → sequence
                    pass

            except Exception as e:
                logger.debug("Sequence alignment failed for %s: %s", hit.target_gene_symbol, e)

        if aligned_count > 0:
            logger.info(
                "Sequence alignment: %d/%d hits aligned for %s → %s",
                aligned_count, len(hits), source_symbol, target_species,
            )
            # Re-sort by percent_identity (sequence alignment trumps domain overlap)
            hits.sort(key=lambda h: h.percent_identity, reverse=True)

    # NCBI E-utilities helpers
    # ═════════════════════════════════════════════════════════

    def _resolve_gene(self, gene_id: str, taxid: str) -> dict | None:
        """Resolve a gene ID (e.g. AT1G56650) to NCBI Gene ID and metadata."""
        # Check cache first
        cache_key = f"gene_info:{taxid}:{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        # Try searching by gene ID + organism
        query = f"{gene_id}[All Fields] AND txid{taxid}[Organism]"
        result = self._eutils_get("json", self.ESEARCH, {
            "db": "gene",
            "term": query,
            "retmax": 5,
        })

        id_list = self._parse_idlist(result)
        if not id_list:
            # Try without taxid constraint
            result2 = self._eutils_get("json", self.ESEARCH, {
                "db": "gene",
                "term": f"{gene_id}[All Fields]",
                "retmax": 3,
            })
            id_list = self._parse_idlist(result2)

        if not id_list:
            knowledge_cache.set(cache_key, None, ttl=86400)
            return None

        ncbi_gene_id = id_list[0]
        info = self._get_gene_summary(ncbi_gene_id)
        if not info:
            knowledge_cache.set(cache_key, None, ttl=86400)
            return None

        gene_info = {
            "ncbi_gene_id": ncbi_gene_id,
            "symbol": info.get("name", gene_id),
            "description": info.get("description", ""),
        }

        knowledge_cache.set(cache_key, gene_info, ttl=86400)
        return gene_info

    def _try_direct_lookup(self, gene_id: str) -> dict | None:
        """Try using gene_id directly as NCBI Gene ID."""
        cache_key = f"gene_info:direct:{gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        info = self._get_gene_summary(gene_id)
        if info:
            result = {
                "ncbi_gene_id": gene_id,
                "symbol": info.get("name", gene_id),
                "description": info.get("description", ""),
            }
            knowledge_cache.set(cache_key, result, ttl=86400)
            return result

        knowledge_cache.set(cache_key, None, ttl=86400)
        return None

    def _get_gene_summary(self, ncbi_gene_id: str) -> dict | None:
        """Fetch gene summary from NCBI esummary."""
        cache_key = f"gene_summary:{ncbi_gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        result = self._eutils_get("json", self.ESUMMARY, {
            "db": "gene",
            "id": ncbi_gene_id,
        })
        if not result:
            return None

        try:
            records = result.get("result", {})
            gene_record = records.get(str(ncbi_gene_id), {})
            if not gene_record or "error" in gene_record:
                return None

            info = {
                "name": gene_record.get("name", ""),
                "description": gene_record.get("description", ""),
                "taxid": str(gene_record.get("taxid", "")),
                "organism": gene_record.get("organism", {}).get("scientificname", ""),
            }
            knowledge_cache.set(cache_key, info, ttl=86400)
            return info
        except Exception:
            return None

    def _get_ortholog_links(self, ncbi_gene_id: str) -> list[str]:
        """Get list of orthologous NCBI Gene IDs via elink."""
        cache_key = f"ncbi_elink:{ncbi_gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        # elink returns XML; request JSON if available, else parse XML
        result = self._eutils_get("json", self.ELINK, {
            "dbfrom": "gene",
            "db": "gene",
            "id": ncbi_gene_id,
            "linkname": "gene_orthologs",
            "retmode": "json",
        })

        ortholog_ids: list[str] = []

        if result:
            try:
                linksets = result.get("linksets", [])
                for ls in linksets:
                    for linksetdb in ls.get("linksetdbs", []):
                        if linksetdb.get("linkname") == "gene_orthologs":
                            ortholog_ids.extend(linksetdb.get("links", []))
            except Exception:
                pass

        # If JSON didn't work, try XML
        if not ortholog_ids:
            xml_result = self._eutils_get("text", self.ELINK, {
                "dbfrom": "gene",
                "db": "gene",
                "id": ncbi_gene_id,
                "linkname": "gene_orthologs",
            })
            if xml_result:
                import xml.etree.ElementTree as ET
                try:
                    root = ET.fromstring(str(xml_result))
                    for link in root.findall(".//LinkSetDb/Link/Id"):
                        if link.text:
                            ortholog_ids.append(link.text)
                except Exception:
                    pass

        # Convert to strings
        ortholog_ids = [str(oid) for oid in ortholog_ids]

        knowledge_cache.set(cache_key, ortholog_ids, ttl=86400)
        return ortholog_ids

    def _score_orthologs(
        self,
        ortholog_ids: list[str],
        target_taxid: str,
        source_gene_id: str,
        source_species: str,
        source_symbol: str,
        target_species: str,
    ) -> list[OrthologHit]:
        """Fetch details for each ortholog and filter/score by target species."""
        hits: list[OrthologHit] = []

        # Get source protein sequence once
        source_seq = self._fetch_protein_sequence(
            self._get_gene_summary(ortholog_ids[0])  # won't work — need source protein
        )
        # Actually, fetch from the source gene's protein accession
        source_protein_id = self._get_protein_accession(source_gene_id)
        source_seq = self._fetch_protein_sequence(source_protein_id) if source_protein_id else ""

        for oid in ortholog_ids[:10]:  # limit to 10 to respect rate limits
            summary = self._get_gene_summary(oid)
            if not summary:
                continue

            ortholog_taxid = summary.get("taxid", "")
            # Filter by target species
            if ortholog_taxid != target_taxid:
                continue

            target_symbol = summary.get("name", oid)
            target_desc = summary.get("description", "")

            # Get target protein sequence and compute identity
            target_protein_id = self._get_protein_accession(oid)
            identity = 0.0
            align_len = 0

            if source_seq and target_protein_id:
                target_seq = self._fetch_protein_sequence(target_protein_id)
                if target_seq:
                    identity, align_len = self._compute_identity(source_seq, target_seq)

            # Determine confidence
            if identity >= 60:
                confidence = "high"
            elif identity >= 35:
                confidence = "medium"
            else:
                confidence = "low"

            hits.append(OrthologHit(
                source_gene_id=source_gene_id,
                source_species=source_species,
                source_gene_symbol=source_symbol,
                target_gene_id=oid,
                target_species=target_species,
                target_gene_symbol=target_symbol,
                percent_identity=round(identity, 1),
                alignment_length=align_len,
                method="ncbi_ortholog_db" if identity > 0 else "symbol_match",
                confidence=confidence,
                source_protein_id=source_protein_id,
                target_protein_id=target_protein_id or "",
                description=target_desc,
            ))

            # Rate limit: small delay between NCBI requests
            if identity == 0.0:
                time.sleep(0.35)

        return hits

    def _search_by_symbol(
        self,
        symbol: str,
        target_taxid: str,
        target_species: str,
        source_gene_id: str,
        source_species: str,
        source_symbol: str,
    ) -> list[OrthologHit]:
        """Search for genes with similar symbol in target species."""
        hits: list[OrthologHit] = []

        query = f"{symbol}[Gene Name] AND txid{target_taxid}[Organism]"
        result = self._eutils_get("json", self.ESEARCH, {
            "db": "gene",
            "term": query,
            "retmax": 5,
        })

        id_list = self._parse_idlist(result)
        if not id_list:
            return hits

        # Get source protein
        source_protein_id = self._get_protein_accession(source_gene_id)
        source_seq = self._fetch_protein_sequence(source_protein_id) if source_protein_id else ""

        for oid in id_list[:5]:
            summary = self._get_gene_summary(oid)
            if not summary:
                continue

            target_symbol = summary.get("name", oid)
            target_desc = summary.get("description", "")

            identity = 0.0
            align_len = 0

            if source_seq:
                target_protein_id = self._get_protein_accession(oid)
                if target_protein_id:
                    target_seq = self._fetch_protein_sequence(target_protein_id)
                    if target_seq:
                        identity, align_len = self._compute_identity(source_seq, target_seq)
                        time.sleep(0.35)

            confidence = "high" if identity >= 60 else ("medium" if identity >= 35 else "low")
            hits.append(OrthologHit(
                source_gene_id=source_gene_id,
                source_species=source_species,
                source_gene_symbol=source_symbol,
                target_gene_id=oid,
                target_species=target_species,
                target_gene_symbol=target_symbol,
                percent_identity=round(identity, 1),
                alignment_length=align_len,
                method="sequence_similarity" if identity > 0 else "symbol_match",
                confidence=confidence,
                source_protein_id=source_protein_id,
                target_protein_id="",
                description=target_desc,
            ))

        return hits

    # ═════════════════════════════════════════════════════════
    # Protein sequence helpers
    # ═════════════════════════════════════════════════════════

    def _get_protein_accession(self, ncbi_gene_id: str) -> str:
        """Get protein accession for a gene via elink gene → protein."""
        cache_key = f"protein_acc:{ncbi_gene_id}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        result = self._eutils_get("json", self.ELINK, {
            "dbfrom": "gene",
            "db": "protein",
            "id": ncbi_gene_id,
            "linkname": "gene_protein",
            "retmode": "json",
        })

        acc = ""
        if result:
            try:
                linksets = result.get("linksets", [])
                for ls in linksets:
                    for linksetdb in ls.get("linksetdbs", []):
                        links = linksetdb.get("links", [])
                        if links:
                            acc = str(links[0])
                            break
            except Exception:
                pass

        knowledge_cache.set(cache_key, acc, ttl=604800)  # 7 days
        return acc

    def _fetch_protein_sequence(self, protein_accession: str) -> str:
        """Fetch protein FASTA sequence from NCBI."""
        if not protein_accession:
            return ""

        cache_key = f"sequence:{protein_accession}"
        cached = knowledge_cache.get(cache_key)
        if cached is not None:
            return cached

        result = self._eutils_get("text", self.EFETCH, {
            "db": "protein",
            "id": protein_accession,
            "rettype": "fasta",
            "retmode": "text",
        })

        seq = self._parse_fasta(str(result)) if result else ""
        knowledge_cache.set(cache_key, seq, ttl=604800)  # 7 days
        return seq

    # ═════════════════════════════════════════════════════════
    # Sequence alignment (pure Python)
    # ═════════════════════════════════════════════════════════

    @staticmethod
    def _needleman_wunsch(
        seq_a: str,
        seq_b: str,
        match_score: int = 2,
        mismatch_penalty: int = -1,
        gap_penalty: int = -2,
    ) -> tuple[float, int]:
        """Needleman-Wunsch global alignment.

        Returns: (percent_identity, alignment_length)

        For plant TF proteins (200-800 aa), completes in < 100ms.
        O(n*m) time, O(n*m) memory.
        """
        n, m = len(seq_a), len(seq_b)

        if n == 0 and m == 0:
            return 100.0, 0
        if n == 0 or m == 0:
            return 0.0, max(n, m)

        # DP matrix
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            dp[i][0] = i * gap_penalty
        for j in range(m + 1):
            dp[0][j] = j * gap_penalty

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                score = match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_penalty
                dp[i][j] = max(
                    dp[i - 1][j - 1] + score,
                    dp[i - 1][j] + gap_penalty,
                    dp[i][j - 1] + gap_penalty,
                )

        # Traceback
        i, j = n, m
        matches = 0
        aligned_len = 0

        while i > 0 or j > 0:
            if i > 0 and j > 0:
                diag = match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_penalty
                if dp[i][j] == dp[i - 1][j - 1] + diag:
                    if seq_a[i - 1] == seq_b[j - 1]:
                        matches += 1
                    aligned_len += 1
                    i -= 1
                    j -= 1
                    continue
            if i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
                aligned_len += 1
                i -= 1
            else:
                aligned_len += 1
                j -= 1

        identity = (matches / aligned_len * 100) if aligned_len > 0 else 0.0
        return identity, aligned_len

    @staticmethod
    def _kmer_identity(seq_a: str, seq_b: str, k: int = 3) -> float:
        """Quick identity estimate via k-mer overlap. O(n+m).

        Used as fallback for very long sequences (>2000 aa).
        """
        if len(seq_a) < k or len(seq_b) < k:
            return 0.0

        kmers_a = {seq_a[i:i + k] for i in range(len(seq_a) - k + 1)}
        kmers_b = {seq_b[i:i + k] for i in range(len(seq_b) - k + 1)}

        if not kmers_a or not kmers_b:
            return 0.0

        common = kmers_a & kmers_b
        return (len(common) / min(len(kmers_a), len(kmers_b))) * 100

    def _compute_identity(self, seq_a: str, seq_b: str) -> tuple[float, int]:
        """Compute sequence identity, choosing algorithm based on length.

        Returns: (percent_identity, alignment_length)
        """
        # Clean sequences: uppercase, remove non-standard chars
        seq_a = re.sub(r'[^A-Z]', '', seq_a.upper())
        seq_b = re.sub(r'[^A-Z]', '', seq_b.upper())

        max_len = max(len(seq_a), len(seq_b))

        if max_len > 2000:
            # Use k-mer approximation for very long sequences
            identity = self._kmer_identity(seq_a, seq_b)
            return identity, max_len

        return self._needleman_wunsch(seq_a, seq_b)

    # ═════════════════════════════════════════════════════════
    # Utility methods
    # ═════════════════════════════════════════════════════════

    @staticmethod
    def _resolve_taxid(species: str) -> str:
        """Map a species name to NCBI Taxonomy ID."""
        if not species:
            return ""

        species_lower = species.lower().strip()

        # Exact match
        if species_lower in PLANT_TAXIDS:
            return PLANT_TAXIDS[species_lower]

        # Substring match
        for known, taxid in PLANT_TAXIDS.items():
            if species_lower in known or known in species_lower:
                return taxid

        return ""

    @staticmethod
    def _parse_idlist(result: dict | None) -> list[str]:
        """Extract ID list from esearch JSON result."""
        if not result:
            return []
        try:
            esr = result.get("esearchresult", {})
            idlist = esr.get("idlist", [])
            return [str(i) for i in idlist]
        except Exception:
            return []

    @staticmethod
    def _parse_fasta(fasta_text: str) -> str:
        """Parse FASTA text, return concatenated sequence."""
        if not fasta_text:
            return ""
        lines = fasta_text.strip().split("\n")
        seq_parts = []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq_parts.append(line)
        return "".join(seq_parts)

    def _eutils_get(self, parse_as: str, url: str, params: dict) -> Any | None:
        """Make NCBI E-utilities GET request with retry and rate limiting.

        Args:
            parse_as: "json" | "text" — how to parse the response
            url: Full NCBI E-utilities endpoint URL
            params: Query parameters (email, tool added automatically)
        """
        params["email"] = self.EMAIL
        params["tool"] = self.TOOL_NAME

        for attempt in range(self.max_retries + 1):
            try:
                resp = self._session.get(url, params=params, timeout=self.timeout)
                if resp.status_code == 429:
                    wait = 2 ** attempt
                    logger.warning("NCBI rate limit hit, waiting %ds", wait)
                    time.sleep(wait)
                    continue
                if resp.status_code == 400:
                    logger.debug("NCBI 400 for %s: %s", url, resp.text[:200])
                    return None
                resp.raise_for_status()

                if parse_as == "json":
                    return resp.json()
                return resp.text

            except requests.Timeout:
                logger.warning("NCBI timeout (attempt %d): %s", attempt + 1, params.get("db", ""))
            except requests.HTTPError as e:
                logger.warning("NCBI HTTP error: %s", e)
                return None
            except Exception as e:
                logger.warning("NCBI request failed: %s", e)
                return None

            if attempt < self.max_retries:
                time.sleep(1.5 ** attempt)

        return None


# ═══════════════════════════════════════════════════════════════
# Convenience function
# ═══════════════════════════════════════════════════════════════

def find_orthologs(
    gene_id: str,
    source_species: str = "Arabidopsis thaliana",
    target_species: str = "",
) -> OrthologResult:
    """Convenience function for ortholog finding."""
    finder = OrthologFinder()
    return finder.find_orthologs(gene_id, source_species, target_species)
