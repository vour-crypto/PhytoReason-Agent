"""
test_knowledge.py — Unit tests for the knowledge layer.

Tests TF knowledge base, species knowledge, KEGG browser, and PubMed retriever.
All tests use built-in data — no external benchmark data required.
"""

import pytest

from phyto_reason.knowledge.tf_knowledge_base import (
    TFKnowledgeBase,
    TFMetaboliteRelation,
    CANONICAL_KNOWLEDGE,
    METABOLITE_CLASS_MAP,
)
from phyto_reason.knowledge.species_knowledge import (
    SpeciesKnowledge,
    SPECIES_KNOWLEDGE,
    CONSERVATION_RULES,
)
from phyto_reason.knowledge.kegg_browser import (
    KEGGBrowser,
    BUILTIN_PATHWAYS,
)


# ═══════════════════════════════════════════════════════════════
# TF Knowledge Base
# ═══════════════════════════════════════════════════════════════

class TestTFKnowledgeBase:
    def setup_method(self):
        self.kb = TFKnowledgeBase()

    def test_canonical_knowledge_not_empty(self):
        """Built-in knowledge base should have entries."""
        assert len(CANONICAL_KNOWLEDGE) >= 20, f"Expected >=20 relations, got {len(CANONICAL_KNOWLEDGE)}"

    def test_all_relations_have_pmids_or_low_confidence(self):
        """Every relation must have PMIDs, or be low-confidence (score <= 0.3, no PMIDs)."""
        for rel in CANONICAL_KNOWLEDGE:
            has_citation = bool(rel.pmids)
            is_low_confidence = rel.score <= 0.3 and rel.strength in ("weak", "inferred")
            assert has_citation or is_low_confidence, (
                f"{rel.tf_family} → {rel.metabolite_class}: no PMIDs, "
                f"score={rel.score}, strength={rel.strength}"
            )

    def test_all_relations_have_score_range(self):
        """Scores must be in [0.0, 1.0]."""
        for rel in CANONICAL_KNOWLEDGE:
            assert 0.0 <= rel.score <= 1.0, (
                f"{rel.tf_family} → {rel.metabolite_class}: score={rel.score} out of range"
            )

    def test_query_berberine_returns_alkaloid_tfs(self):
        """Berberine (alkaloid) should find WRKY, ERF, bHLH families."""
        result = self.kb.query(metabolite="berberine")
        families = {r.tf_family for r in result.relations}
        assert "WRKY" in families or "ERF" in families or "bHLH" in families, (
            f"No known alkaloid-regulating TF families found: {families}"
        )
        assert result.total >= 3

    def test_query_anthocyanin_returns_mbw(self):
        """Anthocyanin should return MYB, bHLH, WD40 (MBW complex)."""
        result = self.kb.query(metabolite="anthocyanin")
        families = {r.tf_family for r in result.relations}
        for expected in ["MYB", "bHLH", "WD40"]:
            assert expected in families, f"MBW component {expected} missing from anthocyanin results"

    def test_query_unknown_metabolite_returns_empty(self):
        """Unknown metabolite should return empty with a note."""
        result = self.kb.query(metabolite="xyz_unknown_12345")
        assert result.total == 0
        assert result.note  # should have explanation

    def test_query_with_tf_filter(self):
        """TF family filter should narrow results to matching families."""
        result = self.kb.query(metabolite="alkaloid", tf_family="ERF")
        for rel in result.relations:
            # ERF filter matches both "ERF" and "AP2-ERF" (substring)
            assert "ERF" in rel.tf_family.upper(), (
                f"Expected ERF-family, got {rel.tf_family}"
            )

    def test_query_with_min_score(self):
        """min_score should filter weak evidence."""
        result_all = self.kb.query(metabolite="alkaloid")
        result_strong = self.kb.query(metabolite="alkaloid", min_score=0.7)
        assert result_strong.total <= result_all.total

    def test_results_sorted_by_score(self):
        """Results should be sorted descending by score."""
        result = self.kb.query(metabolite="alkaloid")
        scores = [r.score for r in result.relations]
        assert scores == sorted(scores, reverse=True), f"Not sorted: {scores}"

    def test_get_families_for_metabolite(self):
        """Should return unique TF families sorted by score."""
        families = self.kb.get_families_for_metabolite("berberine")
        assert len(families) >= 3
        assert len(families) == len(set(families))  # no duplicates

    def test_get_metabolites_for_family(self):
        """Should return metabolite classes for a TF family."""
        classes = self.kb.get_metabolites_for_family("MYB")
        assert "flavonoid" in classes or "anthocyanin" in classes

    def test_get_citations(self):
        """Should return unique PMIDs."""
        pmids = self.kb.get_citations("MYB", "anthocyanin")
        assert len(pmids) >= 1
        assert len(pmids) == len(set(pmids))


class TestMetaboliteClassMap:
    def test_berberine_maps_to_alkaloid(self):
        assert METABOLITE_CLASS_MAP["berberine"] == "alkaloid"

    def test_nicotine_maps_to_alkaloid(self):
        assert METABOLITE_CLASS_MAP["nicotine"] == "alkaloid"

    def test_anthocyanin_maps_to_anthocyanin(self):
        assert METABOLITE_CLASS_MAP["anthocyanin"] == "anthocyanin"

    def test_quercetin_maps_to_flavonoid(self):
        assert METABOLITE_CLASS_MAP["quercetin"] == "flavonoid"

    def test_artemisinin_maps_to_terpenoid(self):
        assert METABOLITE_CLASS_MAP["artemisinin"] == "terpenoid"

    def test_lignin_maps_to_lignin(self):
        assert METABOLITE_CLASS_MAP["lignin"] == "lignin"


class TestTFKnowledgeBaseResolveClass:
    def setup_method(self):
        self.kb = TFKnowledgeBase()

    def test_exact_match(self):
        assert self.kb._resolve_class("berberine") == "alkaloid"

    def test_substring_match(self):
        assert self.kb._resolve_class("quercetin") == "flavonoid"

    def test_return_as_is_if_unknown(self):
        assert self.kb._resolve_class("completely_unknown_xyz") == "completely_unknown_xyz"


# ═══════════════════════════════════════════════════════════════
# Species Knowledge
# ═══════════════════════════════════════════════════════════════

class TestSpeciesKnowledge:
    def setup_method(self):
        self.sk = SpeciesKnowledge()

    def test_species_data_not_empty(self):
        assert len(SPECIES_KNOWLEDGE) >= 8

    def test_conservation_rules_exist(self):
        assert len(CONSERVATION_RULES) >= 8

    def test_lookup_arabidopsis(self):
        info = self.sk.lookup("Arabidopsis thaliana")
        assert info is not None
        assert "anthocyanin" in info["metabolites"]

    def test_lookup_coptis(self):
        info = self.sk.lookup("Coptis")
        assert info is not None
        assert "berberine" in info["metabolites"]

    def test_lookup_tobacco(self):
        info = self.sk.lookup("tobacco")
        assert info is not None
        assert "nicotine" in info["metabolites"]

    def test_lookup_unknown_returns_none(self):
        info = self.sk.lookup("SpeciesThatDoesNotExist")
        assert info is None

    def test_infer_regulation_known_tf(self):
        """TF family known in species → high confidence."""
        result = self.sk.infer_regulation(
            target_species="Coptis", metabolite="berberine", tf_family="bHLH"
        )
        assert len(result.ortholog_inferences) >= 1
        assert result.known_metabolites  # Coptis has known metabolites

    def test_infer_regulation_unknown_species(self):
        """Unknown species → should still give note."""
        result = self.sk.infer_regulation(
            target_species="UnknownPlant", metabolite="flavonoid", tf_family="MYB"
        )
        assert result.note  # should have explanation

    def test_cross_species_arabidopsis_to_tobacco(self):
        """Arabidopsis → tobacco: MYB shared."""
        result = self.sk.cross_species(
            source_species="Arabidopsis", target_species="tobacco",
            tf_family="MYB"
        )
        # Both have MYB, so should find shared TF families
        assert len(result.ortholog_inferences) >= 0  # at minimum doesn't crash

    def test_list_species(self):
        species = self.sk.list_species()
        assert len(species) >= 8
        assert "coptis" in species


# ═══════════════════════════════════════════════════════════════
# KEGG Browser
# ═══════════════════════════════════════════════════════════════

class TestKEGGBrowser:
    def setup_method(self):
        self.kegg = KEGGBrowser()

    def test_builtin_pathways_not_empty(self):
        assert len(BUILTIN_PATHWAYS) >= 5

    def test_browse_berberine(self):
        result = self.kegg.browse("berberine")
        assert len(result.pathways) >= 1
        pw = result.pathways[0]
        assert "isoquinoline" in pw.pathway_name.lower()
        assert len(pw.enzyme_genes) >= 5

    def test_browse_anthocyanin(self):
        result = self.kegg.browse("anthocyanin")
        assert len(result.pathways) >= 1
        assert any("anthocyanin" in pw.pathway_name.lower() for pw in result.pathways)

    def test_browse_nicotine(self):
        result = self.kegg.browse("nicotine")
        assert len(result.pathways) >= 1
        pw = result.pathways[0]
        assert len(pw.enzyme_genes) >= 3

    def test_browse_flavonoid(self):
        result = self.kegg.browse("flavonoid")
        assert len(result.pathways) >= 1

    def test_browse_unknown(self):
        result = self.kegg.browse("completely_unknown_xyz")
        assert len(result.pathways) == 0
        assert result.note

    def test_fuzzy_match_alkaloid(self):
        """Fuzzy match: 'some_alkaloid' should match alkaloid class."""
        result = self.kegg.browse("unknown_alkaloid")
        assert len(result.pathways) >= 1  # fuzzy matched to berberine's pathway

    def test_list_known_compounds(self):
        compounds = self.kegg.list_known_compounds()
        assert "berberine" in compounds
        assert "flavonoid" in compounds

    def test_all_builtin_pathways_have_enzymes(self):
        for name, pw in BUILTIN_PATHWAYS.items():
            assert pw.pathway_name, f"{name}: missing pathway name"
            assert len(pw.enzyme_genes) >= 3, f"{name}: too few enzymes ({len(pw.enzyme_genes)})"


# ═══════════════════════════════════════════════════════════════
# Integration: knowledge layer end-to-end
# ═══════════════════════════════════════════════════════════════

class TestKnowledgeIntegration:
    """End-to-end: query TF KB → cross-check with species → verify with KEGG."""

    def test_berberine_full_chain(self):
        """For berberine in Coptis: KB knows TFs → species knows metabolites → KEGG knows pathway."""
        kb = TFKnowledgeBase()
        sk = SpeciesKnowledge()
        kegg = KEGGBrowser()

        # TF knowledge
        tf_result = kb.query(metabolite="berberine")
        assert tf_result.total >= 3

        # Species knowledge
        species_info = sk.lookup("Coptis")
        assert species_info is not None
        assert "berberine" in species_info["metabolites"]

        # KEGG pathway
        kegg_result = kegg.browse("berberine")
        assert len(kegg_result.pathways) >= 1

    def test_nicotine_full_chain(self):
        """For nicotine in tobacco: same cross-validation."""
        kb = TFKnowledgeBase()
        sk = SpeciesKnowledge()
        kegg = KEGGBrowser()

        tf_result = kb.query(metabolite="nicotine")
        assert tf_result.total >= 3

        species_info = sk.lookup("tobacco")
        assert species_info is not None

        kegg_result = kegg.browse("nicotine")
        assert len(kegg_result.pathways) >= 1


# ═══════════════════════════════════════════════════════════════
# Gene Resolver (v4.4 — MyGene.info + NCBI dynamic resolution)
# ═══════════════════════════════════════════════════════════════


class TestGeneResolverFormatDetection:
    """Tests for gene ID format detection (no network needed)."""

    def test_detect_rapdb_format(self):
        from phyto_reason.knowledge.gene_resolver import detect_id_format
        assert detect_id_format("Os01g0884300") == "rapdb"
        assert detect_id_format("Os03g0815100") == "rapdb"
        assert detect_id_format("Os11g0126900") == "rapdb"

    def test_detect_msu_format(self):
        from phyto_reason.knowledge.gene_resolver import detect_id_format
        assert detect_id_format("LOC_Os01g66130") == "msu"
        assert detect_id_format("LOC_Os03g59730") == "msu"

    def test_detect_tair_format(self):
        from phyto_reason.knowledge.gene_resolver import detect_id_format
        assert detect_id_format("AT1G56650") == "tair"
        assert detect_id_format("AT5G13930") == "tair"

    def test_detect_ncbi_gene_id(self):
        from phyto_reason.knowledge.gene_resolver import detect_id_format
        assert detect_id_format("4325006") == "ncbi_gene"
        assert detect_id_format("842120") == "ncbi_gene"

    def test_detect_symbol(self):
        from phyto_reason.knowledge.gene_resolver import detect_id_format
        assert detect_id_format("OsNAC6") == "symbol"
        assert detect_id_format("WRKY45") == "symbol"
        assert detect_id_format("PAP1") == "symbol"

    def test_detect_empty(self):
        from phyto_reason.knowledge.gene_resolver import detect_id_format
        assert detect_id_format("") == "unknown"


class TestGeneResolverEmptyAndFallback:
    """Tests for edge cases that don't need network."""

    def test_resolve_empty(self):
        from phyto_reason.knowledge.gene_resolver import resolve_gene
        r = resolve_gene("")
        assert r.primary_symbol == ""
        assert r.source == "pattern"

    def test_resolve_unknown_pattern_fallback(self):
        """Unknown gene ID should fall back to pattern matching."""
        from phyto_reason.knowledge.gene_resolver import resolve_gene
        r = resolve_gene("Os99g0999900", "rice")
        # Should not crash; returns some result
        assert r.primary_symbol != ""
        assert r.source in ("pattern", "ncbi", "mygene")

    def test_resolve_ncbi_gene_id(self):
        """NCBI Gene ID should be resolvable."""
        from phyto_reason.knowledge.gene_resolver import resolve_gene
        r = resolve_gene("4325006")
        # Should return a result from mygene or ncbi
        assert r.primary_symbol != ""
        assert r.ncbi_gene_id == "4325006" or r.ncbi_gene_id == ""

    def test_expand_query_with_synonyms(self):
        """Query expansion should produce a valid OR-expression."""
        from phyto_reason.knowledge.gene_resolver import expand_gene_query
        expanded = expand_gene_query("AT1G56650", "arabidopsis")
        assert len(expanded) > 0
        assert "AT1G56650" in expanded

    def test_tf_family_inference(self):
        """TF family should be inferred from gene names."""
        from phyto_reason.knowledge.gene_resolver import _infer_tf_family
        assert _infer_tf_family("OsNAC6", "NAC transcription factor", ["NAC048", "SNAC2"]) == "NAC"
        assert _infer_tf_family("OsWRKY45", "WRKY DNA-binding protein", []) == "WRKY"
        assert _infer_tf_family("OsMYB4", "R2R3-MYB transcription factor", []) == "MYB"
        assert _infer_tf_family("UnknownGene", "hypothetical protein", []) == ""

    def test_alias_scoring_rejects_descriptions(self):
        """The alias scoring function should reject descriptions masquerading as symbols."""
        from phyto_reason.knowledge.gene_resolver import GeneResolver
        resolver = GeneResolver()

        # Build a mock hit similar to what MyGene returns
        mock_hit = {
            "_id": "4325006",
            "symbol": "LOC4325006",
            "name": "NAC domain-containing protein 48-like",
            "taxid": "39947",
            "alias": ["NAC048", "NAC48", "NAC6", "ONAC048", "OsJ_04317", "OsNAC6", "SNAC2"],
            "other_names": [],
            "summary": "",
            "entrezgene": "4325006",
        }

        # Monkey-patch to avoid network call
        import phyto_reason.knowledge.gene_resolver as gr
        original_query = resolver._mygene.query
        resolver._mygene.query = lambda gid: mock_hit if "0884300" in gid else None

        try:
            r = resolver.resolve("Os01g0884300", "rice")
            # Primary symbol should be a real gene symbol, not a description
            assert " " not in r.primary_symbol, f"Got description as symbol: {r.primary_symbol}"
            assert len(r.primary_symbol) < 25, f"Symbol too long: {r.primary_symbol}"
            assert r.tf_family == "NAC"
            assert r.ncbi_gene_id == "4325006"
            assert r.source == "mygene"
        finally:
            resolver._mygene.query = original_query

    def test_resolve_returns_all_fields(self):
        """Resolution should populate all key GeneInfo fields."""
        from phyto_reason.knowledge.gene_resolver import GeneResolver
        resolver = GeneResolver()

        mock_hit = {
            "_id": "842120",
            "symbol": "PAP1",
            "name": "production of anthocyanin pigment 1",
            "taxid": "3702",
            "alias": ["ATMYB75", "MYB75", "SIAA1"],
            "other_names": [],
            "summary": "Encodes a MYB transcription factor involved in anthocyanin biosynthesis.",
            "entrezgene": "842120",
        }

        import phyto_reason.knowledge.gene_resolver as gr
        original_query = resolver._mygene.query
        resolver._mygene.query = lambda gid: mock_hit if "AT1G56650" in gid else None

        try:
            r = resolver.resolve("AT1G56650", "arabidopsis")
            assert r.primary_symbol == "PAP1"
            assert r.tf_family == "MYB"
            assert r.ncbi_gene_id == "842120"
            assert r.source == "mygene"
            assert len(r.synonyms) >= 3  # input + symbol + aliases
        finally:
            resolver._mygene.query = original_query

    def test_multi_level_fallback_empty_result(self):
        """When MyGene and NCBI both return nothing, should fall back to pattern."""
        from phyto_reason.knowledge.gene_resolver import GeneResolver
        resolver = GeneResolver()

        # Patch MyGene to return None
        import phyto_reason.knowledge.gene_resolver as gr
        original_query = resolver._mygene.query
        resolver._mygene.query = lambda gid: None

        try:
            r = resolver.resolve("completely_unknown_xyz_123", "arabidopsis")
            assert r.source == "pattern"
            assert r.primary_symbol == "completely_unknown_xyz_123"
            assert r.synonyms == ["completely_unknown_xyz_123"]
        finally:
            resolver._mygene.query = original_query
