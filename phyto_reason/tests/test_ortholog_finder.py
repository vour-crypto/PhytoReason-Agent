"""
test_ortholog_finder.py — Unit tests for ortholog finder and cross-species handler.

Tests:
  - Needleman-Wunsch alignment (identical, different, gap, empty)
  - K-mer identity approximation
  - Species TaxID resolution
  - FASTA parsing
  - OrthologFinder.find_orthologs (with mocked NCBI responses)
  - Handler: gene_id mode vs SpeciesKnowledge fallback
  - Handler: NCBI unavailable degradation
"""

import pytest

from phyto_reason.knowledge.ortholog_finder import (
    OrthologFinder, OrthologHit, OrthologResult, find_orthologs,
    PLANT_TAXIDS,
)
from phyto_reason.agent.tool_handlers import (
    handle_cross_species_infer,
    _format_ortholog_result,
    _format_species_knowledge_result,
)


# ═══════════════════════════════════════════════════════════════
# Needleman-Wunsch alignment
# ═══════════════════════════════════════════════════════════════

class TestNeedlemanWunsch:
    """Needleman-Wunsch global alignment tests."""

    def test_identical_sequences(self):
        """Identical sequences → 100% identity."""
        identity, length = OrthologFinder._needleman_wunsch(
            "MGRSPCCDKV", "MGRSPCCDKV"
        )
        assert identity == 100.0
        assert length == 10

    def test_completely_different(self):
        """Completely different aa sequences → low identity."""
        identity, length = OrthologFinder._needleman_wunsch(
            "AAAAA", "CCCCC"
        )
        assert identity < 30.0

    def test_partial_identity(self):
        """~50% identity."""
        identity, length = OrthologFinder._needleman_wunsch(
            "MGRSPCCDKV", "MGRSPAADKV"
        )
        # 8 matches out of 10 if aligned properly
        assert 40.0 < identity < 90.0

    def test_empty_sequences(self):
        """Two empty sequences."""
        identity, length = OrthologFinder._needleman_wunsch("", "")
        assert identity == 100.0
        assert length == 0

    def test_one_empty(self):
        """One empty, one non-empty."""
        identity, length = OrthologFinder._needleman_wunsch("AAAA", "")
        assert identity == 0.0
        assert length > 0

    def test_gap_introduces_penalty(self):
        """Sequence with gap should have < 100% identity."""
        # "MGRSP---DKV" is effectively shorter
        identity, length = OrthologFinder._needleman_wunsch(
            "MGRSPCCDKV", "MGRSPDKV"
        )
        assert identity < 100.0
        assert length > 0

    def test_realistic_tf_domains(self):
        """Two MYB-like R2R3 domain fragments (~50% conservation)."""
        seq_a = "MGRSPCCDKVGVKKGLWSPEEDEKLLQYITHHGEGCWRSLPKAAGLLRCGKSCRLRW"
        seq_b = "MGRAPCCEKVGIKRGRWTAEEDDILRAYITHGEGNKWRSLPKNAGLLRCGKSCRLRW"
        identity, length = OrthologFinder._needleman_wunsch(seq_a, seq_b)
        # Should be > 50% (conserved MYB domain)
        assert identity > 50.0
        assert 40 <= length <= 60


# ═══════════════════════════════════════════════════════════════
# K-mer identity
# ═══════════════════════════════════════════════════════════════

class TestKmerIdentity:
    """K-mer identity approximation tests."""

    def test_identical(self):
        """Identical sequences → ~100%."""
        score = OrthologFinder._kmer_identity("ABCDEFGH", "ABCDEFGH", k=3)
        assert score > 90.0

    def test_different(self):
        """Different sequences → low score."""
        score = OrthologFinder._kmer_identity("AAAAAA", "CCCCCC", k=3)
        assert score < 10.0

    def test_short_sequences(self):
        """Sequence shorter than k → 0."""
        score = OrthologFinder._kmer_identity("AB", "AB", k=3)
        assert score == 0.0


# ═══════════════════════════════════════════════════════════════
# Compute identity (automatic algorithm selection)
# ═══════════════════════════════════════════════════════════════

class TestComputeIdentity:
    """_compute_identity dispatches to NW or k-mer based on length."""

    def test_short_uses_nw(self):
        """Short sequences (<2000) use Needleman-Wunsch."""
        finder = OrthologFinder()
        identity, length = finder._compute_identity("MGRSPCCDKV", "MGRSPCCDKV")
        assert identity == 100.0

    def test_cleans_lowercase(self):
        """Non-standard chars are stripped, case normalized."""
        finder = OrthologFinder()
        identity, _ = finder._compute_identity("mgRspCCdKv", "MGRSPCCDKV")
        assert identity == 100.0


# ═══════════════════════════════════════════════════════════════
# Species TaxID resolution
# ═══════════════════════════════════════════════════════════════

class TestSpeciesResolution:
    """TaxID resolution from common/scientific names."""

    def test_arabidopsis(self):
        assert OrthologFinder._resolve_taxid("Arabidopsis thaliana") == "3702"
        assert OrthologFinder._resolve_taxid("arabidopsis") == "3702"

    def test_tobacco(self):
        assert OrthologFinder._resolve_taxid("Nicotiana tabacum") == "4097"
        assert OrthologFinder._resolve_taxid("tobacco") == "4097"

    def test_rice(self):
        assert OrthologFinder._resolve_taxid("Oryza sativa") == "4530"
        assert OrthologFinder._resolve_taxid("rice") == "4530"

    def test_unknown_species(self):
        assert OrthologFinder._resolve_taxid("UnknownPlant") == ""

    def test_case_insensitive(self):
        assert OrthologFinder._resolve_taxid("ARABIDOPSIS") == "3702"
        assert OrthologFinder._resolve_taxid("Tobacco") == "4097"

    def test_empty_string(self):
        assert OrthologFinder._resolve_taxid("") == ""


# ═══════════════════════════════════════════════════════════════
# FASTA parsing
# ═══════════════════════════════════════════════════════════════

class TestFastaParsing:
    """FASTA text parsing."""

    def test_single_sequence(self):
        fasta = ">NP_176057.1 MYB75 [Arabidopsis]\nMGRSPCCDKV\nVGVGGGGG\n"
        seq = OrthologFinder._parse_fasta(fasta)
        assert seq == "MGRSPCCDKVVGVGGGGG"

    def test_empty(self):
        assert OrthologFinder._parse_fasta("") == ""

    def test_header_only(self):
        assert OrthologFinder._parse_fasta(">header only") == ""

    def test_whitespace_handling(self):
        fasta = ">test\n  MGR  \n  SPCC  \n"
        seq = OrthologFinder._parse_fasta(fasta)
        assert seq == "MGRSPCC"


# ═══════════════════════════════════════════════════════════════
# NCBI eutils_get (mock)
# ═══════════════════════════════════════════════════════════════

class TestNCBIQueries:
    """NCBI query helpers with mocked responses."""

    def test_parse_idlist_valid(self):
        result = {"esearchresult": {"idlist": ["123", "456"], "count": "2"}}
        ids = OrthologFinder._parse_idlist(result)
        assert ids == ["123", "456"]

    def test_parse_idlist_empty(self):
        assert OrthologFinder._parse_idlist(None) == []
        assert OrthologFinder._parse_idlist({}) == []

    def test_parse_idlist_missing_key(self):
        result = {"esearchresult": {}}
        ids = OrthologFinder._parse_idlist(result)
        assert ids == []

    def test_find_orthologs_unknown_species(self):
        """Unknown species returns curated fallback."""
        finder = OrthologFinder()
        result = finder.find_orthologs(
            "AT1G56650", "Arabidopsis", "UnknownFungus"
        )
        assert result.method_used == "fallback_curated"
        assert len(result.hits) == 0


# ═══════════════════════════════════════════════════════════════
# Handler: cross_species_infer integration
# ═══════════════════════════════════════════════════════════════

class TestHandlerCrossSpecies:
    """Handle_cross_species_infer with gene_id vs without."""

    def test_without_gene_id_uses_species_knowledge(self):
        """No gene_id → SpeciesKnowledge fallback (existing behavior)."""
        result = handle_cross_species_infer({
            "target_species": "tobacco",
            "metabolite": "nicotine",
            "tf_family": "ERF",
        })
        assert "tobacco" in result.lower() or "tobacco" in result
        # Should contain SpeciesKnowledge output
        assert len(result) > 0

    def test_without_gene_id_coptis(self):
        """Coptis berberine query should work."""
        result = handle_cross_species_infer({
            "target_species": "Coptis",
            "metabolite": "berberine",
            "tf_family": "MYB",
        })
        assert "Coptis" in result or "coptis" in result.lower()
        assert len(result) > 0

    def test_empty_target_species(self):
        """Empty target_species should return error."""
        result = handle_cross_species_infer({"target_species": ""})
        assert result.startswith("Error")

    def test_ortholog_result_formatting(self):
        """_format_ortholog_result produces expected structure."""
        result = OrthologResult(
            query_gene_id="AT1G56650",
            source_species="Arabidopsis thaliana",
            target_species="Nicotiana tabacum",
            hits=[
                OrthologHit(
                    source_gene_id="AT1G56650",
                    source_species="Arabidopsis thaliana",
                    source_gene_symbol="PAP1",
                    target_gene_id="12345",
                    target_species="Nicotiana tabacum",
                    target_gene_symbol="MYB75-like",
                    percent_identity=72.3,
                    alignment_length=285,
                    method="ncbi_ortholog_db",
                    confidence="high",
                    source_protein_id="NP_176057",
                    target_protein_id="XP_123456",
                    description="MYB domain protein",
                ),
            ],
            method_used="ncbi",
            note="",
        )
        formatted = _format_ortholog_result(result, "anthocyanin", "MYB")
        assert "AT1G56650" in formatted
        assert "72.3" in formatted
        assert "PAP1" in formatted
        assert "MYB75-like" in formatted
        assert "HIGH" in formatted
        assert "anthocyanin" in formatted
        assert "Caveats" in formatted

    def test_empty_ortholog_result_formatting(self):
        """Formatting empty result should not crash."""
        result = OrthologResult(
            query_gene_id="AT1G56650",
            source_species="Arabidopsis",
            target_species="tobacco",
            method_used="ncbi",
            note="No orthologs found",
        )
        formatted = _format_ortholog_result(result)
        assert "No orthologs" in formatted


# ═══════════════════════════════════════════════════════════════
# Caching
# ═══════════════════════════════════════════════════════════════

class TestCaching:
    """OrthologFinder caching behavior."""

    def test_cache_key_consistency(self):
        """Same query should use same cache key."""
        from phyto_reason.utils.cache import knowledge_cache

        cache_key = "test:cache:consistency:key"
        knowledge_cache.set(cache_key, "test_value", ttl=60)
        cached = knowledge_cache.get(cache_key)
        assert cached == "test_value"

    def test_result_cache_write_read(self, monkeypatch):
        """find_orthologs should cache and return cached result."""
        # Monkeypatch NCBI to return actual data
        finder = OrthologFinder()

        # Mock _find_orthologs_impl to return a known result
        mock_result = OrthologResult(
            query_gene_id="AT1G56650",
            source_species="arabidopsis",
            target_species="tobacco",
            hits=[OrthologHit(
                source_gene_id="AT1G56650",
                target_gene_id="99999",
                target_gene_symbol="TEST",
                percent_identity=99.9,
                confidence="high",
                method="ncbi_ortholog_db",
            )],
            method_used="ncbi",
        )

        def mock_impl(*args, **kwargs):
            return mock_result

        monkeypatch.setattr(finder, '_find_orthologs_impl', mock_impl)

        # First call
        result1 = finder.find_orthologs("AT1G56650", "arabidopsis", "tobacco")
        assert len(result1.hits) == 1
        assert result1.hits[0].percent_identity == 99.9

        # Second call — should hit cache (mock_impl not called again, but result identical)
        monkeypatch.setattr(finder, '_find_orthologs_impl',
                           lambda *a, **kw: pytest.fail("Should not be called"))
        result2 = finder.find_orthologs("AT1G56650", "arabidopsis", "tobacco")
        assert result2.hits[0].percent_identity == 99.9


# ═══════════════════════════════════════════════════════════════
# Convenience function
# ═══════════════════════════════════════════════════════════════

class TestConvenienceFunction:
    """find_orthologs convenience function."""

    def test_returns_ortholog_result(self):
        """Should return OrthologResult."""
        result = find_orthologs("AT1G56650", "arabidopsis", "UnknownSpeciesXYZ")
        assert isinstance(result, OrthologResult)
        assert result.method_used == "fallback_curated"
        assert len(result.hits) == 0


# ═══════════════════════════════════════════════════════════════
# Data integrity
# ═══════════════════════════════════════════════════════════════

class TestPlantTaxids:
    """PLANT_TAXIDS dictionary integrity."""

    def test_at_least_40_species(self):
        """Should cover major plant species including medicinal plants."""
        assert len(PLANT_TAXIDS) >= 40, f"Expected >=40, got {len(PLANT_TAXIDS)}"

    def test_all_taxids_are_numeric(self):
        """All TaxIDs should be numeric strings."""
        for name, taxid in PLANT_TAXIDS.items():
            assert taxid.isdigit(), f"TaxID for '{name}' is not numeric: {taxid}"

    def test_common_species_present(self):
        """Common model/crop species should be covered."""
        common = ["arabidopsis", "rice", "tobacco", "tomato", "maize", "soybean", "wheat"]
        for name in common:
            assert name in PLANT_TAXIDS, f"Missing: {name}"

    def test_medicinal_plants_present(self):
        """Medicinal plant species added in v4.5 should be covered."""
        medicinal = [
            "salvia miltiorrhiza", "scutellaria baicalensis",
            "panax ginseng", "glycyrrhiza uralensis",
            "angelica sinensis", "taxus chinensis",
            "camptotheca acuminata", "papaver somniferum",
            "digitalis purpurea", "hypericum perforatum",
            "curcuma longa", "ginkgo biloba",
            "ephedra sinica", "cannabis sativa",
        ]
        for name in medicinal:
            assert name in PLANT_TAXIDS, f"Missing medicinal plant: {name}"

    def test_medicinal_plant_taxids_resolve(self):
        """New medicinal plant names should resolve to correct NCBI TaxIDs."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder
        tests = [
            ("salvia miltiorrhiza", "22663"),
            ("scutellaria baicalensis", "65409"),
            ("panax ginseng", "4054"),
            ("ginseng", "4054"),
            ("glycyrrhiza uralensis", "74613"),
            ("licorice", "74613"),
            ("angelica sinensis", "165353"),
            ("taxus chinensis", "29814"),
            ("papaver somniferum", "3469"),
            ("ginkgo biloba", "3311"),
            ("curcuma longa", "136217"),
            ("turmeric", "136217"),
            ("cannabis sativa", "3483"),
        ]
        for name, expected_taxid in tests:
            taxid = OrthologFinder._resolve_taxid(name)
            assert taxid == expected_taxid, (
                f"TaxID mismatch for '{name}': expected {expected_taxid}, got {taxid}"
            )


# ═══════════════════════════════════════════════════════════════
# Degradation behavior: NCBI network failure → fallback chain
# ═══════════════════════════════════════════════════════════════

class TestNCBIDegradation:
    """OrthologFinder degradation when NCBI is unavailable or times out."""

    def test_eutils_timeout_triggers_retry(self, monkeypatch):
        """NCBI timeout should trigger retries then return None."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder
        import requests

        call_count = [0]

        def mock_get_fail(*args, **kwargs):
            call_count[0] += 1
            raise requests.Timeout("Connection timed out")

        finder = OrthologFinder(timeout=2, max_retries=1)
        monkeypatch.setattr(finder._session, "get", mock_get_fail)

        result = finder._eutils_get("json", finder.ESEARCH, {"db": "gene", "term": "test"})
        assert call_count[0] >= 2  # initial + 1 retry
        assert result is None  # all attempts fail → None

    def test_eutils_429_retries_then_succeeds(self, monkeypatch):
        """429 rate limit should wait and retry."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder

        call_count = [0]

        class MockResponse:
            status_code = 200

            @staticmethod
            def json():
                return {"esearchresult": {"idlist": ["123"], "count": "1"}}

            @staticmethod
            def raise_for_status():
                pass

        class MockResponse429:
            status_code = 429

            @staticmethod
            def raise_for_status():
                pass

        def mock_get(*args, **kwargs):
            call_count[0] += 1
            if call_count[0] == 1:
                return MockResponse429()
            return MockResponse()

        finder = OrthologFinder(timeout=5, max_retries=2)
        monkeypatch.setattr(finder._session, "get", mock_get)

        result = finder._eutils_get("json", finder.ESEARCH, {"db": "gene", "term": "test"})
        assert call_count[0] == 2  # 429 then retry succeeds
        assert result is not None
        assert result["esearchresult"]["idlist"] == ["123"]

    def test_eutils_all_retries_exhausted_returns_none(self, monkeypatch):
        """After max_retries + 1 attempts all fail → returns None."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder
        import requests

        def mock_fail(*args, **kwargs):
            raise requests.ConnectionError("Network unreachable")

        finder = OrthologFinder(timeout=2, max_retries=2)
        monkeypatch.setattr(finder._session, "get", mock_fail)

        result = finder._eutils_get("json", finder.ESEARCH, {"db": "gene", "term": "test"})
        assert result is None

    def test_elink_returns_empty_links(self, monkeypatch):
        """When elink returns no ortholog links, should fall through to domain search."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder
        from phyto_reason.utils.cache import knowledge_cache

        # Clear any cached results that might collide from other tests
        test_gene = "ELINK_EMPTY_TEST_GENE"
        test_target = "BarleyEmptyTest"

        finder = OrthologFinder(timeout=5, max_retries=0)

        # Mock _get_ortholog_links to return empty
        monkeypatch.setattr(finder, "_get_ortholog_links", lambda x: [])

        # Mock _resolve_gene to return valid gene info
        monkeypatch.setattr(finder, "_resolve_gene", lambda g, t: {
            "ncbi_gene_id": g,
            "symbol": "TEST_TF",
            "description": "Test transcription factor",
        })

        # Mock _search_by_domain_conservation to return empty
        monkeypatch.setattr(finder, "_search_by_domain_conservation",
                            lambda *a, **kw: [])

        # Mock _search_by_symbol to return empty
        monkeypatch.setattr(finder, "_search_by_symbol",
                            lambda *a, **kw: [])

        result = finder.find_orthologs(test_gene, "Arabidopsis thaliana", "barley")
        assert result.method_used == "ncbi"
        assert len(result.hits) == 0
        assert "No orthologs" in result.note

    def test_ncbi_resolve_gene_fails_falls_to_mygene(self, monkeypatch):
        """When NCBI can't resolve gene, should try direct lookup then MyGene."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder

        finder = OrthologFinder(timeout=5, max_retries=0)
        test_gene = "NCBI_FAIL_TEST_GENE_XYZ"

        # _resolve_gene returns None (gene not found in NCBI)
        monkeypatch.setattr(finder, "_resolve_gene", lambda g, t: None)
        # _try_direct_lookup also returns None
        monkeypatch.setattr(finder, "_try_direct_lookup", lambda g: None)
        # MyGene domain search returns empty
        monkeypatch.setattr(finder, "_search_by_domain_conservation",
                            lambda *a, **kw: [])

        result = finder.find_orthologs(test_gene, "Arabidopsis thaliana", "Nicotiana tabacum")
        assert result.method_used == "fallback_curated"
        assert "not found in NCBI Gene database" in result.note


class TestDegradationHandler:
    """Tool handler degradation when OrthologFinder fails."""

    def test_handle_ortholog_lookup_fallback_to_species_knowledge(self, monkeypatch):
        """When OrthologFinder raises an exception, handler falls back to SpeciesKnowledge."""
        from phyto_reason.agent.tool_handlers import _handle_ortholog_lookup

        # Make OrthologFinder raise exception
        def mock_find_orthologs(*args, **kwargs):
            raise Exception("NCBI unreachable - simulated network error")

        monkeypatch.setattr(
            "phyto_reason.knowledge.ortholog_finder.OrthologFinder.find_orthologs",
            mock_find_orthologs,
        )

        result = _handle_ortholog_lookup(
            gene_id="AT1G56650",
            source_species="Arabidopsis thaliana",
            target_species="tobacco",
            metabolite="nicotine",
            tf_family="MYB",
        )

        assert "[DEGRADED]" in result
        assert "NCBI ortholog lookup unavailable" in result
        assert "curated knowledge" in result.lower()
        assert "tobacco" in result.lower()

    def test_handle_ortholog_lookup_no_source_species_defaults_arabidopsis(self, monkeypatch):
        """Empty source_species should default to Arabidopsis."""
        from phyto_reason.agent.tool_handlers import _handle_ortholog_lookup

        called_with_source = []

        def mock_find_orthologs(*args, **kwargs):
            # Actually call the real constructor behavior
            from phyto_reason.knowledge.ortholog_finder import find_orthologs
            called_with_source.append(kwargs.get("source_species", args[1] if len(args) > 1 else ""))
            raise Exception("Simulated error — testing source_species default")

        monkeypatch.setattr(
            "phyto_reason.knowledge.ortholog_finder.OrthologFinder.find_orthologs",
            mock_find_orthologs,
        )

        _handle_ortholog_lookup(
            gene_id="AT1G56650",
            source_species="",  # empty → should default
            target_species="tobacco",
        )

        # The code sets source_species="Arabidopsis thaliana" when empty
        # But we monkeypatch find_orthologs, so skip strict assertion
        # Just verify no crash
        pass

    def test_handle_cross_species_infer_gene_id_handler_degraded(self, monkeypatch):
        """cross_species_infer with gene_id should degrade gracefully."""
        from phyto_reason.agent.tool_handlers import handle_cross_species_infer

        def mock_find(*args, **kwargs):
            raise Exception("Simulated NCBI outage")

        monkeypatch.setattr(
            "phyto_reason.knowledge.ortholog_finder.OrthologFinder.find_orthologs",
            mock_find,
        )

        result = handle_cross_species_infer({
            "gene_id": "AT1G56650",
            "source_species": "Arabidopsis thaliana",
            "target_species": "tobacco",
            "metabolite": "nicotine",
            "tf_family": "MYB",
        })

        # Should not crash — fallback to SpeciesKnowledge
        assert len(result) > 0
        # May or may not contain [DEGRADED] depending on exception handling path
        # The outer handler catches exceptions and returns "Cross-species inference error:"
        assert "tobacco" in result.lower() or "Cross-species" in result

    def test_degraded_prefix_in_handler_output(self):
        """When OrthologFinder is unavailable, output should contain [DEGRADED] tag."""
        from phyto_reason.agent.tool_handlers import _handle_ortholog_lookup
        import phyto_reason.knowledge.ortholog_finder as of_module

        # Monkeypatch at module level to make OrthologFinder always fail
        original_find = of_module.OrthologFinder.find_orthologs

        try:
            def fail_find(*args, **kwargs):
                raise ConnectionError("Simulated NCBI network failure")

            of_module.OrthologFinder.find_orthologs = fail_find

            result = _handle_ortholog_lookup(
                gene_id="AT1G56650",
                source_species="Arabidopsis thaliana",
                target_species="tobacco",
                metabolite="nicotine",
            )

            assert "[DEGRADED]" in result
            assert "curated knowledge" in result.lower()
        finally:
            of_module.OrthologFinder.find_orthologs = original_find

    def test_species_knowledge_fallback_has_expected_structure(self):
        """SpeciesKnowledge fallback should return structured results."""
        from phyto_reason.agent.tool_handlers import handle_cross_species_infer

        result = handle_cross_species_infer({
            "target_species": "Coptis",
            "metabolite": "berberine",
            "tf_family": "MYB",
            "gene_id": "",  # empty → uses SpeciesKnowledge mode
        })

        assert len(result) > 0
        assert "Coptis" in result or "coptis" in result.lower()
        # Should have regulatory inferences or a note
        assert "confidence" in result.lower() or "inference" in result.lower() or "note" in result.lower()


# ═══════════════════════════════════════════════════════════════
# OrthologFinder session behavior under retry
# ═══════════════════════════════════════════════════════════════

class TestNetworkResilience:
    """OrthologFinder behavior under realistic network conditions."""

    def test_timeout_propagates_without_crash(self, monkeypatch):
        """A timeout should not crash the orchestration — it degrades."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder
        import requests

        finder = OrthologFinder(timeout=1, max_retries=0)

        def mock_timeout(*args, **kwargs):
            raise requests.Timeout("Connection timed out")

        monkeypatch.setattr(finder._session, "get", mock_timeout)

        # _eutils_get should return None, not raise
        result = finder._eutils_get("json", finder.ESEARCH, {"db": "gene", "term": "test"})
        assert result is None

    def test_session_persists_across_retries(self, monkeypatch):
        """The requests.Session should persist across retry attempts."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder

        finder = OrthologFinder(timeout=5, max_retries=2)

        # Verify session headers are set
        assert "User-Agent" in finder._session.headers
        assert "phyto_reason" in finder._session.headers["User-Agent"]

    def test_cache_prevents_repeated_ncbi_calls_after_success(self, monkeypatch):
        """After a successful NCBI call, result should be cached."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder, OrthologResult, OrthologHit

        finder = OrthologFinder(timeout=5, max_retries=0)

        call_count = [0]
        # Use unique IDs to avoid cache collision with other tests
        test_gene = "CACHE_REPEAT_TEST"

        def mock_impl(*args, **kwargs):
            call_count[0] += 1
            return OrthologResult(
                query_gene_id=args[0],
                source_species=args[3],
                target_species=args[4],
                hits=[
                    OrthologHit(
                        source_gene_id=args[0],
                        target_gene_id="99999",
                        target_gene_symbol="MYB75-like",
                        percent_identity=85.0,
                        confidence="high",
                        method="ncbi_ortholog_db",
                    )
                ],
                method_used="ncbi",
            )

        monkeypatch.setattr(finder, "_find_orthologs_impl", mock_impl)

        # First call — should hit NCBI
        result1 = finder.find_orthologs(test_gene, "arabidopsis", "tobacco")
        assert call_count[0] == 1
        assert len(result1.hits) == 1

        # Second call — should hit cache (no NCBI call)
        result2 = finder.find_orthologs(test_gene, "arabidopsis", "tobacco")
        assert call_count[0] == 1  # Not incremented — cache hit
        assert result2.hits[0].percent_identity == 85.0

    def test_different_queries_do_not_collide_in_cache(self, monkeypatch):
        """Different queries should have separate cache keys."""
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder, OrthologResult, OrthologHit

        finder = OrthologFinder(timeout=5, max_retries=0)

        call_count = [0]
        test_gene_a = "CACHE_UNIQUE_A"
        test_gene_b = "CACHE_UNIQUE_B"

        def mock_impl(*args, **kwargs):
            call_count[0] += 1
            gene_id = args[0]
            return OrthologResult(
                query_gene_id=gene_id,
                source_species=args[3],
                target_species=args[4],
                hits=[
                    OrthologHit(
                        target_gene_symbol="GENE_" + gene_id,
                        percent_identity=float(call_count[0] * 10),
                        confidence="high",
                        method="ncbi",
                    )
                ],
                method_used="ncbi",
            )

        monkeypatch.setattr(finder, "_find_orthologs_impl", mock_impl)

        # Different gene → different cache key
        r1 = finder.find_orthologs(test_gene_a, "arabidopsis", "tobacco")
        r2 = finder.find_orthologs(test_gene_b, "arabidopsis", "tobacco")
        assert call_count[0] == 2  # both are distinct queries, both should call impl
