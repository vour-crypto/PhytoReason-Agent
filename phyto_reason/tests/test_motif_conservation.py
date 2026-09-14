"""
test_motif_conservation.py — Tests for promoter motif conservation analysis (v4.5).

Tests:
  - MotifScanner IUPAC pattern compilation
  - Motif scanning on synthetic promoters
  - Motif conservation between ortholog pairs
  - Family conservation rules
  - Edge cases (empty/very short sequences)
"""

from __future__ import annotations

import pytest

from phyto_reason.knowledge.motif_conservation import (
    MotifScanner,
    MotifHit,
    MotifScanResult,
    MotifConservationResult,
    FamilyConservationRule,
    PLANT_TF_MOTIFS,
    FAMILY_CONSERVATION_RULES,
    scan_promoter_motifs,
    compare_promoter_conservation,
    get_family_conservation,
    get_all_conserved_families,
)


# ── Fixtures ──────────────────────────────────────────────────

@pytest.fixture
def scanner():
    return MotifScanner()


@pytest.fixture
def sample_promoter():
    """Synthetic promoter with known motifs embedded."""
    return (
        "ATGCATGCACGTGATGC"        # G-box (CACGTG) at pos 8
        "TTGACTTTGACT"              # W-box (TTGACY) at pos 8 (relative)
        "ACCAACCATGC"               # MYB AC-element (ACCAACC)
        "AGCCGCCATGC"               # GCC-box (AGCCGCC)
        "AAAGAAAGAAAG"              # DOF core (AAAG)
        "TGASTCATGACG"              # bZIP GCN4-like (TGASTCA) / OCS (TGACG)
    )


@pytest.fixture
def divergent_promoter():
    """Promoter with only a few conserved motifs."""
    return (
        "NNNNNNNNNNCACGTGNNNNNNNNNN"   # Only G-box
        "NNNNNNNNNNTTGACYNNNNNNNNNN"   # Only W-box
    )


# ── Test IUPAC compilation ───────────────────────────────────

class TestMotifScannerCompilation:
    def test_basic_iupac(self, scanner):
        """Simple IUPAC patterns should compile."""
        pattern = scanner._compile("CACGTG")
        assert pattern.pattern == "CACGTG"  # plain sequence

    def test_iupac_w(self, scanner):
        """W = [AT], Y = [CT]."""
        pattern = scanner._compile("TTGACY")
        assert pattern.search("TTGACT")   # TTGAC + T  (Y=[CT] → C)
        assert pattern.search("TTGACC")   # TTGAC + C  (Y=[CT] → C)
        assert not pattern.search("TTGATT")  # TTGAT not TTGAC

    def test_iupac_n(self, scanner):
        """N = [ACGT]."""
        pattern = scanner._compile("CANNTG")
        assert pattern.search("CATATG")
        assert pattern.search("CAGCTG")
        assert pattern.search("CACCTG")

    def test_iupac_r(self, scanner):
        """R = [AG]."""
        pattern = scanner._compile("RCCGAC")
        assert pattern.search("ACCGAC")
        assert pattern.search("GCCGAC")
        assert not pattern.search("TCCGAC")

    def test_iupac_m(self, scanner):
        """M = [AC], Y = [CT], B = [CGT]."""
        pattern = scanner._compile("MACGYGB")
        assert pattern.search("AACGCGC")  # M=[AC] A C G Y=[CT] G B=[CGT]

    def test_case_insensitive(self, scanner):
        """Motif scanning should be case-insensitive."""
        pattern = scanner._compile("CACGTG")
        assert pattern.search("cacgtg")
        assert pattern.search("CACGTG")
        assert pattern.search("CacGtG")


# ── Test motif scanning ──────────────────────────────────────

class TestMotifScanning:
    def test_scan_finds_gbox(self, scanner):
        """Should find G-box (CACGTG) in promoter."""
        result = scanner.scan("NNNNCACGTGNNNN", "test_gene")
        gbox_hits = [h for h in result.hits if h.motif_name == "bHLH-Gbox"]
        assert len(gbox_hits) >= 1
        assert gbox_hits[0].tf_family == "bHLH"

    def test_scan_finds_wbox(self, scanner):
        """Should find W-box (TTGACY)."""
        result = scanner.scan("NNNTTGACTNNN", "test_gene")
        wbox_hits = [h for h in result.hits if h.motif_name == "WRKY-Wbox"]
        assert len(wbox_hits) >= 1
        assert wbox_hits[0].tf_family == "WRKY"

    def test_scan_empty_sequence(self, scanner):
        """Empty sequence should return empty result."""
        result = scanner.scan("")
        assert result.promoter_length == 0
        assert len(result.hits) == 0

    def test_scan_short_sequence(self, scanner):
        """Very short sequence should return empty result."""
        result = scanner.scan("AT")
        assert len(result.hits) == 0

    def test_scan_returns_gene_id(self, scanner):
        """MotifScanResult should preserve gene_id."""
        result = scanner.scan("CACGTG", "AT1G56650")
        assert result.gene_id == "AT1G56650"

    def test_scan_finds_multiple_families(self, scanner, sample_promoter):
        """Synthetic promoter should have hits from multiple TF families."""
        result = scanner.scan(sample_promoter)
        families = {h.tf_family for h in result.hits}
        # Should find at least: bHLH, WRKY, MYB, ERF, DOF, bZIP
        assert len(families) >= 4
        assert "bHLH" in families or any("bHLH" in h.motif_name for h in result.hits)
        assert "WRKY" in families or any("WRKY" in h.motif_name for h in result.hits)

    def test_scan_score_range(self, scanner):
        """Motif scores should be in [0, 1]."""
        result = scanner.scan("CACGTGTTGACTACCAACC", "test")
        for hit in result.hits:
            assert 0.0 <= hit.score <= 1.0, f"{hit.motif_name}: score={hit.score}"

    def test_scan_hit_positions(self, scanner):
        """Hit positions should be within promoter bounds."""
        seq = "NNNNCACGTGNNNN"
        result = scanner.scan(seq, "test")
        for hit in result.hits:
            assert 0 <= hit.start < len(seq)
            assert hit.start < hit.end <= len(seq)


# ── Test motif conservation ──────────────────────────────────

class TestMotifConservation:
    def test_identical_promoters(self, scanner):
        """Identical promoters should have high conservation score."""
        seq = "CACGTGTTGACTACCAACC"
        result = scanner.scan_conservation(seq, seq, "geneA", "geneB")
        assert result.conservation_score > 0.7
        assert len(result.shared_families) > 0

    def test_divergent_promoters(self, scanner):
        """Very different promoters should have low conservation score."""
        result = scanner.scan_conservation(
            "CACGTGTTGACTACCAACCAGCCGCCAAAG",
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
            "geneA", "geneB",
        )
        assert result.conservation_score < 0.5

    def test_conservation_shared_families(self, scanner):
        """Should correctly identify shared TF families."""
        # Both have G-box (bHLH) and W-box (WRKY)
        seq_a = "CACGTGTTGACT"
        seq_b = "CACGTGTTGACT"
        result = scanner.scan_conservation(seq_a, seq_b, "geneA", "geneB")
        assert "bHLH" in result.shared_families
        assert "WRKY" in result.shared_families

    def test_conservation_empty_source(self, scanner):
        """Empty source sequence should handle gracefully."""
        result = scanner.scan_conservation("", "CACGTG", "geneA", "geneB")
        assert result.conservation_score == 0.0

    def test_conservation_empty_both(self, scanner):
        """Both empty — should return zero score."""
        result = scanner.scan_conservation("", "", "geneA", "geneB")
        assert result.conservation_score == 0.0


# ── Test family conservation rules ───────────────────────────

class TestFamilyConservationRules:
    def test_all_major_families_have_rules(self):
        """All major plant TF families should have conservation rules."""
        expected = {"MYB", "bHLH", "WRKY", "NAC", "ERF", "bZIP", "DOF", "MADS", "HD-ZIP", "SBP"}
        for fam in expected:
            assert fam in FAMILY_CONSERVATION_RULES, f"Missing rule for {fam}"

    def test_get_family_conservation_known(self):
        """Should return rule for known family."""
        rule = get_family_conservation("MYB")
        assert rule is not None
        assert rule.conserved_in_plants is True
        assert rule.tf_family == "MYB"

    def test_get_family_conservation_unknown(self):
        """Should return None for unknown family."""
        rule = get_family_conservation("UNKNOWN_TF_FAMILY")
        assert rule is None

    def test_get_family_conservation_case_insensitive(self):
        """Should handle case differences."""
        rule = get_family_conservation("myb")
        assert rule is not None
        assert rule.tf_family == "MYB"

    def test_get_all_conserved_families(self):
        """Should return all conserved families."""
        families = get_all_conserved_families()
        assert len(families) >= 8
        assert "MYB" in families
        assert "WRKY" in families
        assert "bHLH" in families

    def test_conservation_rules_have_citations(self):
        """Each rule should have a literature citation (PMID)."""
        for fam, rule in FAMILY_CONSERVATION_RULES.items():
            assert "PMID" in rule.note, f"{fam} rule missing PMID"

    def test_low_divergence_families(self):
        """Families with low divergence should be the most deeply conserved ones."""
        low_div = [
            fam for fam, rule in FAMILY_CONSERVATION_RULES.items()
            if rule.motif_divergence == "low"
        ]
        # MYB, bHLH, WRKY, ERF, DOF, SBP should be low divergence
        assert "MYB" in low_div
        assert "bHLH" in low_div
        assert "WRKY" in low_div


# ── Test motif database ──────────────────────────────────────

class TestMotifDatabase:
    def test_all_motifs_have_pmids(self):
        """Every motif should have a literature citation."""
        for name, (consensus, family, pmid) in PLANT_TF_MOTIFS.items():
            assert pmid.startswith("PMID:"), f"{name} missing PMID"

    def test_motifs_cover_major_families(self):
        """Motifs should cover all major plant TF families."""
        families = {fam for _, (_, fam, _) in PLANT_TF_MOTIFS.items()}
        expected = {"MYB", "bHLH", "WRKY", "ERF", "NAC", "bZIP", "DOF", "MADS", "HD-ZIP", "SBP", "TCP"}
        for fam in expected:
            assert fam in families, f"Missing motifs for {fam}"

    def test_iupac_consensus_valid(self):
        """All IUPAC consensus sequences should use valid codes."""
        valid = set("ACGTRYSWKMBDHVN")
        for name, (consensus, _, _) in PLANT_TF_MOTIFS.items():
            for char in consensus.upper():
                assert char in valid, f"{name}: invalid IUPAC char '{char}'"

    def test_minimal_consensus_length(self):
        """Consensus should be at least 4 bp (avoid false positives)."""
        for name, (consensus, _, _) in PLANT_TF_MOTIFS.items():
            core = consensus.upper().replace("N", "")
            assert len(consensus) >= 4, f"{name} too short: {consensus}"


# ── Test convenience functions ────────────────────────────────

class TestConvenienceFunctions:
    def test_scan_promoter_motifs(self):
        result = scan_promoter_motifs("CACGTGTTGACT", "AT1G56650")
        assert isinstance(result, MotifScanResult)
        assert result.gene_id == "AT1G56650"
        assert len(result.hits) > 0

    def test_compare_promoter_conservation(self):
        result = compare_promoter_conservation(
            "CACGTGTTGACT", "CACGTGTTGACT",
            "AT1G56650", "Os01g0884300",
        )
        assert isinstance(result, MotifConservationResult)
        assert result.conservation_score > 0.5
        assert len(result.shared_families) > 0
