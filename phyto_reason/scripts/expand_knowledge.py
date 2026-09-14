"""
expand_knowledge.py — 从 PubMed 批量检索，扩展 TF-代谢物关系知识库。

用法:
    cd D:/plant_agent
    python phyto_reason/scripts/expand_knowledge.py           # 打印新增条目到终端
    python phyto_reason/scripts/expand_knowledge.py --append  # 自动追加到 knowledge/tf_knowledge_base.py

依赖:
    pip install requests  (如果已有 phyto_reason 环境则已安装)

NCBI E-utilities 限速:
    无 API key:  3 次/秒
    有 API key: 10 次/秒 (在 .env 中设置 NCBI_API_KEY)

输出格式与 CANONICAL_KNOWLEDGE 完全一致，可直接粘贴到 tf_knowledge_base.py。
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests

# ── Config ──────────────────────────────────────────────────
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")
RATE_LIMIT = 0.15 if NCBI_API_KEY else 0.35  # seconds between requests
MAX_RESULTS_PER_QUERY = 10

# ── TF families to expand ───────────────────────────────────
TF_FAMILIES = [
    "MYB", "bHLH", "WRKY", "ERF", "NAC", "bZIP", "WD40",
    "SPL", "GRAS", "ARF", "HD-ZIP", "MADS-box", "TCP", "C2H2",
    "HSF", "Dof", "GATA", "Trihelix", "LOB", "GeBP",
]

# ── Metabolite classes to query ─────────────────────────────
METABOLITE_CLASSES = [
    "alkaloid", "flavonoid", "anthocyanin", "terpenoid",
    "phenylpropanoid", "lignin", "glucosinolate", "phenolic",
]

# ── Specific metabolites for focused queries ────────────────
SPECIFIC_METABOLITES = [
    "berberine", "nicotine", "artemisinin", "taxol", "vinblastine",
    "caffeine", "morphine", "quercetin", "resveratrol", "tanshinone",
    "gossypol", "catharanthine", "solanine", "capsaicin",
]


@dataclass
class PubMedHit:
    """Raw PubMed search hit."""
    pmid: str
    title: str = ""
    year: str = ""
    journal: str = ""
    snippet: str = ""
    tf_families_found: list[str] = field(default_factory=list)
    metabolites_found: list[str] = field(default_factory=list)


@dataclass
class ExtractedRelation:
    """Extracted TF-metabolite relationship ready for knowledge base."""
    tf_family: str
    metabolite_class: str
    specific_metabolite: str = ""
    strength: str = "moderate"
    score: float = 0.5
    description: str = ""
    pmids: list[str] = field(default_factory=list)
    species_scope: str = "general"
    mechanism: str = ""
    known_examples: list[str] = field(default_factory=list)


class PubMedExpander:
    """Query PubMed and extract TF-metabolite regulatory relationships."""

    def __init__(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "PhytoReason-Agent/4.0 (scientific research; mailto:research@example.com)",
        })
        self.results: list[ExtractedRelation] = []
        self._call_count = 0

    def _rate_limit(self) -> None:
        """Respect NCBI rate limits."""
        time.sleep(RATE_LIMIT)
        self._call_count += 1
        if self._call_count % 20 == 0:
            print(f"  ... {self._call_count} API calls made, continuing ...", file=sys.stderr)

    def search(self, query: str, max_results: int = 10) -> list[str]:
        """Search PubMed and return list of PMIDs."""
        params: dict[str, Any] = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "sort": "relevance",
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            self._rate_limit()
            resp = self.session.get(ESEARCH, params=params, timeout=15)
            if resp.status_code != 200:
                return []
            data = resp.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            return id_list
        except Exception as e:
            print(f"  Search error: {e}", file=sys.stderr)
            return []

    def fetch_summaries(self, pmids: list[str]) -> list[PubMedHit]:
        """Fetch article summaries for a list of PMIDs."""
        if not pmids:
            return []

        params: dict[str, Any] = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
        }
        if NCBI_API_KEY:
            params["api_key"] = NCBI_API_KEY

        try:
            self._rate_limit()
            resp = self.session.get(ESUMMARY, params=params, timeout=15)
            if resp.status_code != 200:
                return []
            data = resp.json()
            results = data.get("result", {})
            hits: list[PubMedHit] = []
            for pmid in pmids:
                article = results.get(pmid, {})
                if article and "error" not in article:
                    hit = PubMedHit(
                        pmid=pmid,
                        title=article.get("title", ""),
                        year=article.get("pubdate", "")[:4],
                        journal=article.get("source", ""),
                    )
                    hits.append(hit)
            return hits
        except Exception as e:
            print(f"  Summary error: {e}", file=sys.stderr)
            return []

    def extract_tf_families(self, text: str) -> list[str]:
        """Extract TF family names from text."""
        found = []
        text_upper = text.upper()
        for tf in TF_FAMILIES:
            tf_upper = tf.upper().replace("-", " ")
            # Match as whole word or with common suffixes
            pattern = re.compile(
                rf'\b{re.escape(tf_upper)}\b'
                rf'|\b{re.escape(tf_upper)}S\b'
                rf'|\b{re.escape(tf_upper)}\s*TRANSCRIPTION\s*FACTOR',
                re.IGNORECASE,
            )
            if pattern.search(text):
                found.append(tf)
        return list(set(found))

    def extract_metabolites(self, text: str) -> list[str]:
        """Extract metabolite names from text."""
        found = []
        text_lower = text.lower()
        for meta in METABOLITE_CLASSES + SPECIFIC_METABOLITES:
            if meta in text_lower:
                found.append(meta)
        return list(set(found))

    def classify_metabolite_class(self, metabolite: str) -> str:
        """Map specific metabolite to its class."""
        mapping: dict[str, str] = {
            "berberine": "alkaloid", "nicotine": "alkaloid", "morphine": "alkaloid",
            "caffeine": "alkaloid", "vinblastine": "alkaloid", "catharanthine": "alkaloid",
            "solanine": "alkaloid", "capsaicin": "alkaloid",
            "quercetin": "flavonoid", "resveratrol": "phenylpropanoid",
            "artemisinin": "terpenoid", "taxol": "terpenoid", "tanshinone": "terpenoid",
            "gossypol": "terpenoid", "anthocyanin": "anthocyanin",
            "flavonoid": "flavonoid", "alkaloid": "alkaloid",
            "terpenoid": "terpenoid", "lignin": "lignin",
            "phenylpropanoid": "phenylpropanoid", "glucosinolate": "phenylpropanoid",
            "phenolic": "phenylpropanoid",
        }
        return mapping.get(metabolite, metabolite)

    def assess_strength(self, tf: str, meta_class: str, n_papers: int) -> tuple[str, float]:
        """Heuristic: estimate evidence strength from TF family and paper count."""
        # Well-studied TF-metabolite pairs
        strong_pairs = {
            ("MYB", "flavonoid"), ("MYB", "anthocyanin"), ("MYB", "lignin"),
            ("bHLH", "anthocyanin"), ("WD40", "anthocyanin"),
            ("WRKY", "alkaloid"), ("ERF", "alkaloid"), ("NAC", "lignin"),
            ("bHLH", "flavonoid"), ("bZIP", "flavonoid"),
        }
        moderate_pairs = {
            ("WRKY", "terpenoid"), ("ERF", "terpenoid"), ("bHLH", "alkaloid"),
            ("MYB", "alkaloid"), ("NAC", "anthocyanin"),
            ("bZIP", "alkaloid"), ("SPL", "flavonoid"), ("TCP", "flavonoid"),
            ("ARF", "flavonoid"), ("HD-ZIP", "phenylpropanoid"),
        }

        pair = (tf, meta_class)
        if pair in strong_pairs and n_papers >= 3:
            return "strong", 0.80
        elif pair in strong_pairs or (pair in moderate_pairs and n_papers >= 2):
            return "moderate", 0.55
        elif n_papers >= 2:
            return "moderate", 0.50
        elif n_papers >= 1:
            return "weak", 0.30
        else:
            return "weak", 0.20

    def run(self, dry_run: bool = False) -> list[ExtractedRelation]:
        """Run the full expansion pipeline.

        Queries each TF family × specific metabolite, extracts relationships,
        and returns structured ExtractedRelation objects.

        Args:
            dry_run: If True, only print planned queries without executing.
        """
        all_relations: dict[tuple[str, str], ExtractedRelation] = {}

        # Build query list
        queries: list[tuple[str, str, str]] = []
        for tf in TF_FAMILIES:
            for meta in SPECIFIC_METABOLITES[:8]:  # top 8 most common
                meta_class = self.classify_metabolite_class(meta)
                # Skip known well-covered pairs
                queries.append((tf, meta, meta_class))

        if dry_run:
            print(f"Planned queries: {len(queries)}")
            for tf, meta, mc in queries[:10]:
                print(f"  {tf} × {meta} ({mc})")
            return []

        print(f"Starting expansion: {len(queries)} queries")
        print(f"Rate limit: {RATE_LIMIT}s/request (key={'yes' if NCBI_API_KEY else 'no'})")
        print()

        for i, (tf, meta, meta_class) in enumerate(queries):
            query = f'("{tf}"[All Fields] OR "{tf} transcription factor"[All Fields]) AND ("{meta}"[All Fields]) AND ("regulation"[All Fields] OR "biosynthesis"[All Fields] OR "secondary metabolism"[All Fields])'

            print(f"[{i+1}/{len(queries)}] {tf} × {meta} ...", end=" ", flush=True)

            pmids = self.search(query, max_results=5)
            if not pmids:
                print("0 hits")
                continue

            hits = self.fetch_summaries(pmids)
            print(f"{len(hits)} papers")

            if not hits:
                continue

            # Build description from top hit
            top = hits[0]
            desc = (
                f"{tf} TFs are associated with {meta} ({meta_class}) regulation. "
                f"Top paper: {top.title[:180]}. "
                f"PMID: {top.pmid}."
            )
            if len(hits) >= 2:
                desc += f" Additional evidence: PMID:{hits[1].pmid}."

            strength, score = self.assess_strength(tf, meta_class, len(hits))

            key = (tf, meta_class)
            if key in all_relations and all_relations[key].score >= score:
                continue  # keep the stronger entry

            all_relations[key] = ExtractedRelation(
                tf_family=tf,
                metabolite_class=meta_class,
                specific_metabolite=meta if meta_class != meta else "",
                strength=strength,
                score=score,
                description=desc,
                pmids=[h.pmid for h in hits[:3]],
                species_scope="general",
            )

        self.results = sorted(all_relations.values(), key=lambda r: r.score, reverse=True)
        return self.results


# ═══════════════════════════════════════════════════════════════
# Output formatters
# ═══════════════════════════════════════════════════════════════

def format_as_python(results: list[ExtractedRelation]) -> str:
    """Format results as Python code to paste into CANONICAL_KNOWLEDGE."""
    lines = ["    # ════════════════════════════════════════════════════════════"]
    lines.append("    # Auto-expanded from PubMed (scripts/expand_knowledge.py)")
    lines.append("    # ════════════════════════════════════════════════════════════")

    for r in results:
        lines.append("    TFMetaboliteRelation(")
        lines.append(f'        tf_family="{r.tf_family}", metabolite_class="{r.metabolite_class}",')
        if r.specific_metabolite and r.specific_metabolite != r.metabolite_class:
            lines.append(f'        specific_metabolite="{r.specific_metabolite}",')
        lines.append(f'        strength="{r.strength}", score={r.score:.2f},')
        # Truncate description to fit nicely
        desc = r.description[:300].replace('"', '\\"')
        lines.append(f'        description="{desc}",')
        lines.append(f"        pmids={r.pmids},")
        lines.append(f'        species_scope="{r.species_scope}",')
        if r.mechanism:
            lines.append(f'        mechanism="{r.mechanism}",')
        if r.known_examples:
            lines.append(f"        known_examples={r.known_examples},")
        lines.append("    ),")

    return "\n".join(lines)


def format_as_report(results: list[ExtractedRelation]) -> str:
    """Format results as a human-readable report."""
    lines = [
        "=" * 60,
        f"PubMed Expansion Results — {len(results)} new relations found",
        "=" * 60,
        "",
    ]
    for r in results:
        lines.append(f"## {r.tf_family} → {r.metabolite_class} "
                     f"({r.strength}, {r.score:.0%})")
        if r.specific_metabolite:
            lines.append(f"  Specific: {r.specific_metabolite}")
        lines.append(f"  {r.description[:200]}")
        lines.append(f"  PMIDs: {', '.join(r.pmids)}")
        lines.append(f"  Species: {r.species_scope}")
        lines.append("")

    # Summary by TF family
    by_family: dict[str, int] = {}
    for r in results:
        by_family[r.tf_family] = by_family.get(r.tf_family, 0) + 1
    lines.append("---")
    lines.append("Summary by TF family:")
    for tf, count in sorted(by_family.items(), key=lambda x: -x[1]):
        lines.append(f"  {tf}: {count} relations")
    lines.append(f"  Total: {len(results)}")

    return "\n".join(lines)


def append_to_knowledge_base(results: list[ExtractedRelation], filepath: str) -> bool:
    """Append new relations to tf_knowledge_base.py CANONICAL_KNOWLEDGE list.

    Inserts before the closing ']' of CANONICAL_KNOWLEDGE.

    Returns True on success, False if the marker wasn't found.
    """
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        new_entries = []
        for r in results:
            entry_lines = [
                "    TFMetaboliteRelation(",
                f'        tf_family="{r.tf_family}", metabolite_class="{r.metabolite_class}",',
            ]
            if r.specific_metabolite and r.specific_metabolite != r.metabolite_class:
                entry_lines.append(f'        specific_metabolite="{r.specific_metabolite}",')
            entry_lines.append(f'        strength="{r.strength}", score={r.score:.2f},')
            desc = r.description[:300].replace('"', '\\"')
            entry_lines.append(f'        description="{desc}",')
            entry_lines.append(f"        pmids={r.pmids},")
            entry_lines.append(f'        species_scope="{r.species_scope}",')
            if r.known_examples:
                entry_lines.append(f"        known_examples={r.known_examples},")
            entry_lines.append("    ),")
            new_entries.append("\n".join(entry_lines))

        insertion = "\n" + "\n".join(new_entries) + "\n"

        # Insert before the closing ']' of CANONICAL_KNOWLEDGE
        marker = "\n]"
        if marker not in content:
            print("ERROR: Could not find closing ']' of CANONICAL_KNOWLEDGE.", file=sys.stderr)
            return False

        new_content = content.replace(marker, insertion + marker, 1)

        with open(filepath, "w", encoding="utf-8") as f:
            f.write(new_content)

        return True

    except Exception as e:
        print(f"ERROR: Failed to update knowledge base: {e}", file=sys.stderr)
        return False


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Expand TF-metabolite knowledge base from PubMed",
    )
    parser.add_argument(
        "--append", action="store_true",
        help="Auto-append results to knowledge/tf_knowledge_base.py",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show planned queries without executing",
    )
    parser.add_argument(
        "--format", choices=["python", "report"], default="python",
        help="Output format (default: python)",
    )
    parser.add_argument(
        "--limit", type=int, default=0,
        help="Limit number of queries (for testing)",
    )
    args = parser.parse_args()

    # Paths
    base_dir = Path(__file__).resolve().parent.parent
    kb_path = base_dir / "knowledge" / "tf_knowledge_base.py"

    if not kb_path.exists():
        print(f"ERROR: Knowledge base not found at {kb_path}", file=sys.stderr)
        sys.exit(1)

    expander = PubMedExpander()

    if args.dry_run:
        expander.run(dry_run=True)
        return

    results = expander.run()

    if not results:
        print("No new relations found. The existing knowledge base may already cover these queries.")
        return

    # Output
    if args.format == "python":
        print(format_as_python(results))
    else:
        print(format_as_report(results))

    # Append
    if args.append:
        ok = append_to_knowledge_base(results, str(kb_path))
        if ok:
            print(f"\n{'='*60}")
            print(f"Appended {len(results)} relations to {kb_path}")
            print(f"Run 'python -m phyto_reason.benchmark.scoring' to verify.")
        else:
            print("\nAuto-append failed. Copy the output above and paste manually.")

    print(f"\n{'='*60}")
    print(f"Total: {len(results)} new relations")
    print(f"API calls: {expander._call_count}")
    print(f"Paste the above entries into {kb_path} inside CANONICAL_KNOWLEDGE = [")


if __name__ == "__main__":
    main()
