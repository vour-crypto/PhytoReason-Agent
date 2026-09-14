"""Phase 2 smoke test: one real request per registered Phase 2 adapter."""
from __future__ import annotations

from phyto_reason.tools import TOOL_REGISTRY

CASES = [
    ("europe_pmc_search", {"query": "OsNAC6 salt stress", "max_results": 3}),
    ("openalex_search", {"query": "anthocyanin biosynthesis", "max_results": 3}),
    ("biorxiv_search", {"query": "rice salinity", "max_results": 3}),
    ("crossref_search", {"query": "tanshinone biosynthesis", "max_results": 3}),
    ("arxiv_search", {"query": "transformer protein", "max_results": 3}),
    ("semantic_scholar_search", {"query": "plant transcription factor", "max_results": 3}),
    ("planttfdb_lookup", {"query": "NAC"}),
    ("jaspar_motif_lookup", {"query": "MYB", "max_results": 3}),
    ("ensembl_plants_lookup", {"query": "AT1G01010", "species": "arabidopsis_thaliana"}),
    ("literature_search", {"metabolite": "salt stress rice", "max_results": 3}),
    ("kegg_pathway", {"compound_name": "gibberellin", "analysis_type": "pathway"}),
]


def main() -> None:
    for name, kwargs in CASES:
        tool = TOOL_REGISTRY.get(name)
        if tool is None:
            print(f"{name:26s} | 未注册 <== 需核实注册名")
            continue
        result = tool.run(**kwargs)
        records = result.evidence_list[0].metadata.get("records", []) if result.evidence_list else []
        first = str(records[0])[:70] if records else "-"
        warns = "; ".join(result.warnings) or "-"
        count = result.metadata.get("count", result.metadata.get("n_articles", result.metadata.get("pathways_found", 0)))
        print(f"{name:26s} | {result.status:7s} | count={count} | {first} | {warns}")
        if name == "jaspar_motif_lookup" and records:
            tax_group = records[0].get("tax_group")
            tg_str = ",".join(tax_group) if isinstance(tax_group, list) else str(tax_group)
            verdict = "OK" if "plants" in tg_str else "FAIL <== 过滤未生效，需排查"
            print(f"{'':26s} | JASPAR 过滤验证: tax_group={tg_str} · {verdict}")

    result = TOOL_REGISTRY.get("kegg_pathway").run(compound_name="diterpenoid", analysis_type="pathway")
    names = result.metadata.get("pathway_names", [])
    print(f"{'kegg_fallback':26s} | {result.status:7s} | count={result.metadata.get('count', '?')} | 命中: {names} | {result.warnings}")


if __name__ == "__main__":
    main()
