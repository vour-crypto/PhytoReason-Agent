# Phase 2 tool runtime notes

The seven public literature adapters are `literature_search` (PubMed),
`europe_pmc_search`, `openalex_search`, `biorxiv_search`, `crossref_search`,
`arxiv_search`, and `semantic_scholar_search`. The four general plant adapters
are `planttfdb_lookup`, `jaspar_motif_lookup`, `kegg_pathway`, and
`ensembl_plants_lookup`.

The six adapters in `phase2_tools.py` share a 0.2-second minimum request
interval, three attempts with exponential backoff, and a 24-hour in-memory TTL
cache keyed by endpoint and parameters. KEGG and PubMed retain their existing
client-specific rate limits, retries, and caches. Every adapter returns a
`ToolResult` with source, query, count, records, and a partial-failure warning
when an upstream service is unavailable.
