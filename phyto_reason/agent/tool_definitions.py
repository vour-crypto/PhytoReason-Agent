"""
tool_definitions.py -- DeepSeek function calling JSON schemas.

Defines the tools available to the LLM.
Schema format: OpenAI-compatible (used by DeepSeek).
"""

from __future__ import annotations

# ── Tool Definitions ──────────────────────────────────────────

TOOL_DEFINITIONS = [
    {
        "type": "function",
        "function": {
            "name": "search_literature",
            "description": "Search academic literature across MULTIPLE databases (PubMed + Semantic Scholar + Europe PMC + OpenAlex + arXiv) "
                           "AND the local knowledge base (RAG) -- all in parallel. "
                           "Returns results from all sources in a single response. "
                           "Use this as your PRIMARY tool for any literature or scientific knowledge query. "
                           "Covers: published papers, preprints, agricultural science (Agricola), "
                           "open-access metadata (OpenAlex aggregates WoS/Scopus/Crossref), "
                           "preprints (arXiv q-bio), and previous analysis results stored locally. "
                           "No need to call search_rag separately -- it's already included.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query combining metabolite name, species, and regulation keywords. "
                                       "Example: 'berberine biosynthesis transcription factor Coptis'",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Maximum number of results to return (default 5, max 10)",
                        "default": 5,
                    },
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "query_kegg",
            "description": "Query the KEGG database for metabolic pathway information. "
                           "Returns pathway names, enzyme genes, and pathway class for a given compound.",
            "parameters": {
                "type": "object",
                "properties": {
                    "compound_name": {
                        "type": "string",
                        "description": "The compound/metabolite name to query. "
                                       "Example: 'berberine', 'anthocyanin', 'flavonoid'",
                    },
                },
                "required": ["compound_name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_scientific_pipeline",
            "description": "Run the full scientific analysis pipeline on uploaded omics data. "
                           "Generates, falsifies, and ranks mechanistic hypotheses about what "
                           "regulates the target metabolite. "
                           "REQUIRES: user must have uploaded expression + metabolite data first. "
                           "If no data is uploaded, use search_literature and query_kegg instead.",
            "parameters": {
                "type": "object",
                "properties": {
                    "species": {
                        "type": "string",
                        "description": "Species name. Example: 'Arabidopsis thaliana', 'Coptis chinensis', 'Zanthoxylum nitidum'",
                    },
                    "target_metabolite": {
                        "type": "string",
                        "description": "The metabolite of interest. If left empty, the pipeline auto-detects "
                                       "the most differentially accumulated metabolite from the data. "
                                       "Example: 'berberine', 'anthocyanin', 'alkaloid'",
                    },
                    "target_pathway": {
                        "type": "string",
                        "description": "Optional: target pathway name if known. Example: 'isoquinoline alkaloid biosynthesis'",
                    },
                },
                "required": ["species"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "check_data_quality",
            "description": "Check the quality of uploaded data before running analysis. "
                           "Validates sample counts, missing values, and sample alignment.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Session ID for the current analysis session",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_knowledge_base",
            "description": "Search built-in TF-metabolite regulatory knowledge base. "
                           "Returns known TF families that regulate specific metabolites or metabolite classes, "
                           "with supporting citations (PMIDs). Use this when the user asks what TFs are "
                           "known to regulate a metabolite, especially in zero-data mode. "
                           "Covers: MYB, bHLH, WRKY, ERF, NAC, bZIP, WD40 families and their target metabolite classes.",
            "parameters": {
                "type": "object",
                "properties": {
                    "metabolite": {
                        "type": "string",
                        "description": "Metabolite name or class. Examples: 'berberine', 'alkaloid', 'anthocyanin', 'nicotine', 'flavonoid'",
                    },
                    "tf_family": {
                        "type": "string",
                        "description": "Optional: filter by TF family. Examples: 'MYB', 'bHLH', 'WRKY', 'ERF', 'NAC'",
                    },
                    "species": {
                        "type": "string",
                        "description": "Optional: species name for scope filtering. Examples: 'Arabidopsis', 'Coptis', 'tobacco'",
                    },
                },
                "required": ["metabolite"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "cross_species_infer",
            "description": "🔬 FIND ORTHOLOGS/HOMOLOGS ACROSS SPECIES -- use this FIRST for any cross-species gene query. "
                           "USE WHEN: (1) user asks about a gene in a non-model species (e.g., 'genes like AT1G56650 in tobacco'), "
                           "(2) user mentions orthologs/homologs/同源基因/直系同源, "
                           "(3) user compares two species at the gene level, "
                           "(4) user asks 'what regulates X in <non-model plant>?' -- "
                           "cross-species inference is the ONLY way to get candidate regulators for non-model species. "
                           "When gene_id is provided, queries NCBI databases for real orthologs with protein sequence identity scores. "
                           "Without gene_id, uses built-in TF family conservation rules with curated plant species knowledge.",
            "parameters": {
                "type": "object",
                "properties": {
                    "target_species": {
                        "type": "string",
                        "description": "The species you want to find orthologs in. Example: 'Nicotiana tabacum' or 'tobacco'",
                    },
                    "gene_id": {
                        "type": "string",
                        "description": "Optional: specific gene ID to find orthologs for. Example: 'AT1G56650' (TAIR locus), "
                                       "'836077' (NCBI Gene ID). When provided, real NCBI ortholog computation is used.",
                    },
                    "source_species": {
                        "type": "string",
                        "description": "Optional: species the query gene comes from. Default: 'Arabidopsis thaliana'. "
                                       "Example: 'Arabidopsis thaliana'",
                    },
                    "metabolite": {
                        "type": "string",
                        "description": "Optional: metabolite of interest for regulatory context. Example: 'anthocyanin'",
                    },
                    "tf_family": {
                        "type": "string",
                        "description": "Optional: TF family for regulatory context. Example: 'MYB', 'bHLH'",
                    },
                },
                "required": ["target_species"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "search_rag",
            "description": "Search the local document store for relevant information from previous "
                           "analysis results, cached reports, and exported data tables. "
                           "Use this for deep, focused retrieval from local files when "
                           "search_literature's built-in RAG section needs more detail. "
                           "Note: search_literature already includes RAG results automatically -- "
                           "use this tool only for deeper follow-up queries into local data.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query -- can be a gene name, metabolite, pathway, "
                                       "species name, or natural language question. "
                                       "Example: 'Coptis berberine biosynthesis genes'",
                    },
                    "top_k": {
                        "type": "integer",
                        "description": "Number of results to return (default 5, max 10)",
                        "default": 5,
                    },
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "web_search",
            "description": "Search the web for information about genes, metabolites, species, or any plant biology topic. "
                           "Use this when PubMed, KEGG, and the knowledge base return no results, or when you need "
                           "broad background information not available in structured databases. "
                           "Particularly useful for: resolving unknown gene IDs, finding gene functions, "
                           "checking species-specific information, and filling knowledge gaps.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query -- can include gene IDs, metabolite names, species, "
                                       "or natural language questions. Example: 'Os01g0884300 rice NAC transcription factor gene'",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Number of results to return (default 8, max 10)",
                        "default": 8,
                    },
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "query_public_expression",
            "description": "Query public expression databases (BAR/eFP) for tissue-specific gene expression data, "
                           "expressologs (cross-species expression conservation), and protein interactions. "
                           "USE WHEN: you need tissue expression context for candidate genes, want to check "
                           "if orthologs have similar expression patterns, or need expression-based validation. "
                           "Currently supports Arabidopsis genes (AGI format: AT1G56650).",
            "parameters": {
                "type": "object",
                "properties": {
                    "gene_id": {
                        "type": "string",
                        "description": "Gene ID to query. For Arabidopsis, use AGI format. "
                                       "Example: 'AT1G56650'. Other species gene IDs are resolved via BAR.",
                    },
                },
                "required": ["gene_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "recall_memory",
            "description": "Search the agent's cross-session memory stream for relevant past observations, "
                           "analysis results, user preferences, and insights. "
                           "The memory stream persists across sessions -- use this to recall what was "
                           "discussed or discovered in previous conversations. "
                           "Particularly useful for: remembering user preferences, recalling past "
                           "analysis results, and checking if a question was answered before.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "What to search for in memory. Can be a topic, species name, "
                                       "metabolite, gene family, or natural language query. "
                                       "Example: 'user preference about analysis detail level' or "
                                       "'previous findings about berberine in Coptis'",
                    },
                    "top_k": {
                        "type": "integer",
                        "description": "Number of memories to retrieve (default 5, max 10)",
                        "default": 5,
                    },
                },
                "required": ["query"],
            },
        },
    },

    # ── v5.0 Fine-grained pipeline node tools (ReAct agent) ──

    {
        "type": "function",
        "function": {
            "name": "run_deg_analysis",
            "description": "Run differential expression analysis on uploaded expression data. "
                           "Uses multi-group ANOVA with FDR correction and effect size (eta-squared) ranking. "
                           "Returns: number of genes tested, FDR-significant count, top DEGs ranked by effect size. "
                           "Also generates a volcano plot figure. "
                           "Requires: expression matrix with at least 2 sample groups. "
                           "Call check_data_quality first to understand the data before running this.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Current session ID",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_dam_analysis",
            "description": "Run differential accumulation analysis on uploaded metabolite data. "
                           "Identifies differentially accumulated metabolites between sample groups "
                           "using multi-group ANOVA with FDR correction; when ANOVA is significant, "
                           "runs moderated-t pairwise tissue comparisons with BH correction. "
                           "Group priority is condition > tissue > treatment, with metadata or sample-name fallback marked. "
                           "ANOVA alone is not a tissue-specific conclusion. "
                           "Requires: uploaded metabolite matrix with at least 2 groups.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Current session ID",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_multiomics",
            "description": "Run multi-omics integration analysis. Computes: "
                           "(1) pairwise gene-metabolite correlations (bicor + FDR), "
                           "(2) joint KEGG pathway enrichment, "
                           "(3) quadrant plot (DEG vs DAM log2FC), "
                           "(4) Spearman correlation network with hub detection. "
                           "Requires: DEG and DAM results available. "
                           "Call run_deg_analysis and run_dam_analysis first.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Current session ID",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_wgcna",
            "description": "Run WGCNA (Weighted Gene Co-expression Network Analysis). "
                           "Identifies co-expression modules and correlates them with traits. "
                           "Returns: number of modules, soft power, module sizes, hub genes. "
                           "Requires: expression matrix with 8+ samples.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Current session ID",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_tf_analysis",
            "description": "Identify transcription factor candidates regulating the target metabolite. "
                           "Combines: bicor correlation with pathway genes, promoter motif scanning, "
                           "literature evidence, TF family prior knowledge, and evidence fusion scoring. "
                           "Returns ranked TF candidates with multi-dimensional scores. "
                           "Call after DEG and multi-omics analyses are completed.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Current session ID",
                    },
                    "target_metabolite": {
                        "type": "string",
                        "description": "Target metabolite name. Auto-detected if empty.",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_hypothesis_synthesis",
            "description": "Generate ranked competing mechanistic hypotheses from all accumulated evidence. "
                           "Includes falsification checks, evidence strength comparison, "
                           "and proposed validation experiments. "
                           "Call this LAST as the final analysis step.",
            "parameters": {
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Current session ID",
                    },
                },
                "required": ["session_id"],
            },
        },
    },
    # ── 代谢组学工具 ──

    {
        "type": "function",
        "function": {
            "name": "annotate_ms2_spectrum",
            "description": "Annotate an unknown MS/MS spectrum (metabolite identification). "
                           "Input: path to spectrum file (MGF/CSV) OR inline precursor m/z + fragment ions. "
                           "Output: adduct inference, molecular formula candidates (CHNOPS + golden rules), "
                           "structural class via diagnostic fragment matching (compound_profiles, 31 classes), "
                           "candidate metabolites, and MSI confidence level (2/3/4; Level 1 needs standards). "
                           "Optional species context prioritizes profile marker metabolites.",
            "parameters": {
                "type": "object",
                "properties": {
                    "spectrum_path": {
                        "type": "string",
                        "description": "Path to MS/MS spectrum file (.mgf/.csv/.tsv; mzML needs pyteomics).",
                    },
                    "precursor_mz": {
                        "type": "number",
                        "description": "Precursor ion m/z (used when no file is given).",
                    },
                    "fragments": {
                        "type": "array",
                        "items": {"type": "number"},
                        "description": "Fragment ion m/z values (used when no file is given).",
                    },
                    "species": {
                        "type": "string",
                        "description": "Optional species name — profile marker metabolites are used as anchors.",
                    },
                    "target_metabolite": {
                        "type": "string",
                        "description": "Optional target metabolite to prioritize.",
                    },
                    "remote_query": {
                        "type": "string",
                        "description": "Optional compound name for MassBank/MoNA comparison; omitted means local-only MSI annotation.",
                    },
                    "tolerance_ppm": {
                        "type": "number",
                        "description": "Mass tolerance in ppm (default 10).",
                    },
                    "polarity": {
                        "type": "string",
                        "enum": ["", "positive", "negative"],
                        "description": "Ionization mode (default auto: filename pos/neg or MGF CHARGE).",
                    },
                },
                "required": [],
            },
        },
    },

    # ── 植物专用工具 ──

    {
        "type": "function",
        "function": {
            "name": "query_plantcyc",
            "description": "Query PlantCyc (Plant Metabolic Network) for plant-specific metabolic pathways, "
                           "compounds, and enzyme information. PREFERRED over KEGG for plant metabolism. "
                           "Contains curated pathways for specialized metabolism: flavonoids, alkaloids, "
                           "terpenoids, glucosinolates, etc.",
            "parameters": {
                "type": "object",
                "properties": {
                    "metabolite_name": {
                        "type": "string",
                        "description": "Metabolite/compound name to search (e.g., 'taxifolin', 'berberine').",
                    },
                    "pathway_name": {
                        "type": "string",
                        "description": "Pathway name to look up (e.g., 'flavonoid biosynthesis').",
                    },
                    "enzyme_name": {
                        "type": "string",
                        "description": "Enzyme name to find in pathways (e.g., 'CHS', 'PAL').",
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name for species-specific pathways.",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "massbank_spectrum_search",
            "description": "Search MassBank3 spectra by compound, precursor mass, ion mode, or peak list.",
            "parameters": {"type": "object", "properties": {
                "compound_name": {"type": "string"}, "precursor_mz": {"type": "number"},
                "tolerance_ppm": {"type": "number", "default": 10},
                "ion_mode": {"type": "string"}, "ms_type": {"type": "string"},
                "peaks": {"type": "array"}, "limit": {"type": "integer", "default": 10},
            }, "required": []},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "mona_spectrum_search",
            "description": "Search MoNA spectra using a non-empty Spring Filter query with pagination.",
            "parameters": {"type": "object", "properties": {
                "query": {"type": "string"}, "size": {"type": "integer", "default": 10},
                "page": {"type": "integer", "default": 0},
            }, "required": ["query"]},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "pubchem_compound_properties",
            "description": "Get PubChem compound properties or batch-resolve names to CAS/InChIKey, with explicit 30-minute circuit-breaker fallback and source/actual_source provenance.",
            "parameters": {"type": "object", "properties": {
                "name": {"type": "string"}, "cid": {"type": "string"},
                "properties": {"type": "array"},
                "names": {"type": "array", "items": {"type": "string"}, "description": "Optional batch of compound names."},
            }, "required": []},
        },
    },
]
