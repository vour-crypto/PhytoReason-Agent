"""
motif_conservation.py — Promoter motif conservation analysis (v4.5).

Provides:
  1. Built-in plant TF binding motifs (consensus from PlantTFDB / JASPAR / AthaMap)
  2. Simple PWM scanning (no FIMO dependency)
  3. Motif conservation scoring between ortholog promoters
  4. Integration with ortholog pipeline for cross-species validation

Pure Python implementation — no external bioinformatics tools required.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any, ClassVar

logger = logging.getLogger("motif_conservation")


# ═══════════════════════════════════════════════════════════════
# Data models
# ═══════════════════════════════════════════════════════════════

@dataclass
class MotifHit:
    """A single motif occurrence in a promoter."""
    motif_name: str = ""
    tf_family: str = ""
    consensus: str = ""
    start: int = 0       # 0-based position in promoter
    end: int = 0
    strand: str = "+"
    score: float = 0.0   # 0-1 normalized
    sequence: str = ""   # actual matched sequence


@dataclass
class MotifScanResult:
    """Result of scanning a promoter for motifs."""
    gene_id: str = ""
    promoter_length: int = 0
    hits: list[MotifHit] = field(default_factory=list)
    unique_families: int = 0
    note: str = ""


@dataclass
class MotifConservationResult:
    """Motif conservation between two orthologous promoters."""
    source_gene: str = ""
    target_gene: str = ""
    source_hits: list[MotifHit] = field(default_factory=list)
    target_hits: list[MotifHit] = field(default_factory=list)
    shared_families: list[str] = field(default_factory=list)
    conservation_score: float = 0.0  # 0-1
    note: str = ""


# ═══════════════════════════════════════════════════════════════
# Built-in plant TF binding motifs (IUPAC consensus)
# ═══════════════════════════════════════════════════════════════

# Key: motif name → (IUPAC consensus, TF family, PMID/reference)
PLANT_TF_MOTIFS: dict[str, tuple[str, str, str]] = {
    # ── MYB family ──────────────────────────────────────────
    "MYB-core":        ("YAACKG", "MYB", "PMID:12679534"),
    "MYB-AC-element":  ("ACCAACC", "MYB", "PMID:12068109"),
    "MYB-CCA1":        ("AAAAATCT", "MYB", "PMID:11500506"),
    "MYB-PHR1":        ("GNATATNC", "MYB", "PMID:14675444"),
    "MYB-WER":         ("ACCWAMC", "MYB", "PMID:11069708"),

    # ── bHLH family ─────────────────────────────────────────
    "bHLH-Gbox":       ("CACGTG", "bHLH", "PMID:14675432"),
    "bHLH-Ebox":       ("CANNTG", "bHLH", "PMID:14675432"),
    "bHLH-MYC2":       ("CACATG", "bHLH", "PMID:21799859"),

    # ── WRKY family ─────────────────────────────────────────
    "WRKY-Wbox":       ("TTGACY", "WRKY", "PMID:10758498"),
    "WRKY-WKbox":      ("TTTTCCAC", "WRKY", "PMID:10758498"),

    # ── ERF/AP2 family ──────────────────────────────────────
    "ERF-GCCbox":      ("AGCCGCC", "ERF", "PMID:10950870"),
    "ERF-DRE":         ("RCCGAC", "ERF", "PMID:9851978"),
    "ERF-CRT":         ("CCGAC", "ERF", "PMID:9851978"),

    # ── NAC family ──────────────────────────────────────────
    "NAC-core":        ("CACG", "NAC", "PMID:15805479"),
    "NAC-CBNAC":       ("NNNCTNNNNNNNNAAG", "NAC", "PMID:17582538"),
    "NAC-ANAC":        ("CATGTG", "NAC", "PMID:15805479"),

    # ── bZIP family ─────────────────────────────────────────
    "bZIP-ABRE":       ("MACGYGB", "bZIP", "PMID:11340185"),
    "bZIP-GCN4":       ("TGASTCA", "bZIP", "PMID:19141215"),
    "bZIP-OCS":        ("TGACG", "bZIP", "PMID:19141215"),

    # ── DOF family ──────────────────────────────────────────
    "DOF-core":        ("AAAG", "DOF", "PMID:12581302"),

    # ── MADS family ─────────────────────────────────────────
    "MADS-CArG":       ("CCWWWWWWGG", "MADS", "PMID:10852934"),  # CC(A/T)6GG

    # ── HD-ZIP family ───────────────────────────────────────
    "HDZIP-ATHB":      ("CAATWATTG", "HD-ZIP", "PMID:11597504"),

    # ── SBP family ──────────────────────────────────────────
    "SBP-core":        ("GTAC", "SBP", "PMID:15173568"),

    # ── TCP family ──────────────────────────────────────────
    "TCP-core":        ("GGNCCC", "TCP", "PMID:15304637"),
}

# IUPAC code mapping
IUPAC_MAP: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


# ═══════════════════════════════════════════════════════════════
# MotifScanner
# ═══════════════════════════════════════════════════════════════

class MotifScanner:
    """Scan promoter sequences for known plant TF binding motifs.

    Uses IUPAC consensus patterns — no FIMO/MEME dependency.
    """

    def __init__(self, motifs: dict | None = None):
        """
        Args:
            motifs: Optional custom motif dictionary.
                    Defaults to PLANT_TF_MOTIFS.
        """
        self.motifs = motifs or PLANT_TF_MOTIFS
        self._compiled: dict[str, re.Pattern] = {}

    def _compile(self, iupac: str) -> re.Pattern:
        """Convert IUPAC consensus to compiled regex."""
        if iupac in self._compiled:
            return self._compiled[iupac]

        pattern_parts = []
        for char in iupac.upper():
            mapped = IUPAC_MAP.get(char, char)
            pattern_parts.append(mapped)

        regex = re.compile("".join(pattern_parts), re.IGNORECASE)
        self._compiled[iupac] = regex
        return regex

    def scan(self, sequence: str, gene_id: str = "") -> MotifScanResult:
        """Scan a promoter sequence for all known motifs.

        Args:
            sequence: DNA sequence (promoter region, e.g. -1000bp to TSS)
            gene_id: Optional gene identifier for reporting

        Returns:
            MotifScanResult with all motif hits
        """
        clean_seq = sequence.upper().replace("\n", "").replace(" ", "")
        result = MotifScanResult(
            gene_id=gene_id,
            promoter_length=len(clean_seq),
        )

        if not clean_seq or len(clean_seq) < 4:
            result.note = "Sequence too short for motif scanning"
            return result

        families_seen = set()

        for motif_name, (consensus, tf_family, pmid) in self.motifs.items():
            try:
                pattern = self._compile(consensus)
                for match in pattern.finditer(clean_seq):
                    # Score: motif length / consensus length (simplified PWM)
                    hit_len = match.end() - match.start()
                    motif_len = len(consensus.replace("N", ""))
                    score = min(hit_len / max(motif_len, 1), 1.0)

                    result.hits.append(MotifHit(
                        motif_name=motif_name,
                        tf_family=tf_family,
                        consensus=consensus,
                        start=match.start(),
                        end=match.end(),
                        strand="+",
                        score=score,
                        sequence=match.group(),
                    ))
                    families_seen.add(tf_family)
            except re.error:
                logger.debug("Invalid motif pattern: %s (%s)", motif_name, consensus)

        result.unique_families = len(families_seen)
        result.note = (
            f"Found {len(result.hits)} motif hits across {result.unique_families} "
            f"TF families in {len(clean_seq)}bp promoter"
        )
        return result

    def scan_conservation(
        self,
        source_seq: str,
        target_seq: str,
        source_gene: str = "",
        target_gene: str = "",
    ) -> MotifConservationResult:
        """Compare motif content between two orthologous promoters.

        Returns MotifConservationResult with shared families and conservation score.
        """
        source_result = self.scan(source_seq, source_gene)
        target_result = self.scan(target_seq, target_gene)

        # Find shared TF families
        source_families = {h.tf_family for h in source_result.hits}
        target_families = {h.tf_family for h in target_result.hits}
        shared = sorted(source_families & target_families)

        # Conservation score: Jaccard-like + length-normalized hit count
        total_families = len(source_families | target_families)
        if total_families == 0:
            conservation_score = 0.0
        elif not source_result.hits or not target_result.hits:
            # One side has no motifs → no conservation
            conservation_score = 0.0
        else:
            jaccard = len(shared) / total_families
            # Bonus for shared motifs in similar count ratios
            ratio = min(
                len(source_result.hits) / max(len(target_result.hits), 1),
                len(target_result.hits) / max(len(source_result.hits), 1),
            )
            hit_ratio = 0.5 + 0.5 * ratio

            conservation_score = round(jaccard * 0.7 + hit_ratio * 0.3, 3)

        note = (
            f"Motif conservation: {len(shared)}/{total_families} shared families, "
            f"score={conservation_score:.2f}"
        )

        return MotifConservationResult(
            source_gene=source_gene,
            target_gene=target_gene,
            source_hits=source_result.hits,
            target_hits=target_result.hits,
            shared_families=shared,
            conservation_score=conservation_score,
            note=note,
        )


# ═══════════════════════════════════════════════════════════════
# Built-in per-family conservation rules (no promoters needed)
# ═══════════════════════════════════════════════════════════════

# When no promoter sequences are available, use curated conservation rules
# based on known TF family evolutionary conservation patterns.

@dataclass
class FamilyConservationRule:
    """Motif conservation rule for a TF family."""
    tf_family: str
    conserved_in_plants: bool = True       # Motif is conserved across angiosperms
    motif_divergence: str = "low"          # "low" | "medium" | "high"
    target_clades: list[str] = field(default_factory=list)  # Specific clades
    note: str = ""


# Curated rules based on literature (PMID cited)
FAMILY_CONSERVATION_RULES: dict[str, FamilyConservationRule] = {
    "MYB": FamilyConservationRule(
        tf_family="MYB",
        conserved_in_plants=True,
        motif_divergence="low",
        note="MYB binding sites (AC elements) are deeply conserved in land plants. PMID: 12068109",
    ),
    "bHLH": FamilyConservationRule(
        tf_family="bHLH",
        conserved_in_plants=True,
        motif_divergence="low",
        note="G-box / E-box motifs are conserved from algae to angiosperms. PMID: 14675432",
    ),
    "WRKY": FamilyConservationRule(
        tf_family="WRKY",
        conserved_in_plants=True,
        motif_divergence="low",
        note="W-box (TTGACY) is highly conserved across plant species. PMID: 10758498",
    ),
    "NAC": FamilyConservationRule(
        tf_family="NAC",
        conserved_in_plants=True,
        motif_divergence="medium",
        note="NAC binding shows family-level but not universal conservation. PMID: 17582538",
    ),
    "ERF": FamilyConservationRule(
        tf_family="ERF",
        conserved_in_plants=True,
        motif_divergence="low",
        note="GCC-box (AGCCGCC) and DRE/CRT elements conserved in plants. PMID: 10950870",
    ),
    "bZIP": FamilyConservationRule(
        tf_family="bZIP",
        conserved_in_plants=True,
        motif_divergence="medium",
        note="ABRE and GCN4-like motifs show clade-specific variations. PMID: 11340185",
    ),
    "DOF": FamilyConservationRule(
        tf_family="DOF",
        conserved_in_plants=True,
        motif_divergence="low",
        note="AAAG core motif is conserved in all plant DOF proteins. PMID: 12581302",
    ),
    "MADS": FamilyConservationRule(
        tf_family="MADS",
        conserved_in_plants=True,
        motif_divergence="medium",
        note="CArG box is conserved but flanking preferences vary. PMID: 10852934",
    ),
    "HD-ZIP": FamilyConservationRule(
        tf_family="HD-ZIP",
        conserved_in_plants=True,
        motif_divergence="medium",
        note="HD-ZIP binding shows lineage-specific variations. PMID: 11597504",
    ),
    "SBP": FamilyConservationRule(
        tf_family="SBP",
        conserved_in_plants=True,
        motif_divergence="low",
        note="SBP motif (GTAC) is conserved but flanking context varies. PMID: 15173568",
    ),
}


# ═══════════════════════════════════════════════════════════════
# Convenience functions
# ═══════════════════════════════════════════════════════════════

def scan_promoter_motifs(sequence: str, gene_id: str = "") -> MotifScanResult:
    """Scan a promoter sequence for known plant TF motifs."""
    scanner = MotifScanner()
    return scanner.scan(sequence, gene_id)


def compare_promoter_conservation(
    source_seq: str,
    target_seq: str,
    source_gene: str = "",
    target_gene: str = "",
) -> MotifConservationResult:
    """Compare motif conservation between two orthologous promoters."""
    scanner = MotifScanner()
    return scanner.scan_conservation(source_seq, target_seq, source_gene, target_gene)


def get_family_conservation(tf_family: str) -> FamilyConservationRule | None:
    """Get motif conservation rule for a TF family (no promoter needed)."""
    return FAMILY_CONSERVATION_RULES.get(tf_family.upper())


def get_all_conserved_families(target_clade: str = "") -> list[str]:
    """Get list of TF families with conserved motifs.

    Args:
        target_clade: Optional filter (e.g., "eudicots", "monocots", "angiosperms")
    """
    families = []
    for fam, rule in FAMILY_CONSERVATION_RULES.items():
        if rule.conserved_in_plants:
            if target_clade:
                if target_clade.lower() in [c.lower() for c in rule.target_clades]:
                    families.append(fam)
                elif not rule.target_clades:  # no specific clade = universal
                    families.append(fam)
            else:
                families.append(fam)
    return sorted(families)
