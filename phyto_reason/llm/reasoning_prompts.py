"""
reasoning_prompts.py — 科研推理 prompt 模板。

每个 prompt 包含 {placeholders}，由 PromptManager 填充。
"""

from __future__ import annotations

REASONING_PROMPTS: dict[str, str] = {
    "explain_candidate": (
        "Given the following evidence for transcription factor {gene_id}:\n\n"
        "Fusion Score: {fused_score}\n"
        "Confidence Level: {confidence_level}\n"
        "Contradictions: {contradictions}\n\n"
        "Evidence details:\n{evidence_table}\n\n"
        "Biological context:\n"
        "- Target metabolite: {target_metabolite}\n"
        "- TF family: {tf_family}\n"
        "- Ontology inheritance: {ontology_chain}\n\n"
        "Generate a scientific explanation of why this TF is a candidate.\n"
        "Include: confidence assessment, each evidence dimension's contribution, "
        "biological plausibility, and any caveats.\n"
        "If the confidence level is Gold, explain what makes it Gold.\n"
        "If Silver or Weak, explain what's missing.\n"
    ),
    "compare_candidates": (
        "Compare the following candidate TFs for regulating {target_metabolite}:\n\n{candidate_table}\n\n"
        "For each candidate, explain:\n"
        "1. What evidence supports it\n"
        "2. What evidence is missing or contradictory\n"
        "3. Why it ranks where it does\n\n"
        "Then provide a summary comparison and recommendation.\n"
    ),
    "generate_hypothesis": (
        "Based on the top candidate TFs for regulating {target_metabolite} "
        "in {species}, generate a structured biological hypothesis.\n\n"
        "Top candidates:\n{candidate_summaries}\n\n"
        "Available pathway knowledge:\n{pathway_info}\n\n"
        "Generate hypothesis addressing:\n"
        "1. Which TF(s) most likely regulate the pathway\n"
        "2. Predicted target enzymes\n"
        "3. Regulatory mechanism\n"
        "4. Confidence level per claim\n"
        "5. Key experiments to validate\n"
    ),
    "experimental_design": (
        "Given the following candidate TFs for {target_metabolite}:\n\n{candidate_summaries}\n\n"
        "Suggest specific wet-lab validation experiments. For each suggested experiment:\n"
        "- Which TF to test\n"
        "- Why this experiment is appropriate (based on specific evidence)\n"
        "- Expected outcome\n"
        "- Priority (high/medium/low)\n\n"
        "Available validation methods:\n- CRISPR/Cas9 knockout\n- Overexpression\n- Yeast One-Hybrid (Y1H)\n"
        "- Dual-Luciferase (Dual-LUC)\n- qRT-PCR\n- Electrophoretic Mobility Shift Assay (EMSA)\n"
        "- Chromatin Immunoprecipitation (ChIP-qPCR)\n"
    ),
    "troubleshoot": (
        "The user is reporting an issue with their analysis. Given:\n\n"
        "Error/Warning: {error_info}\n"
        "Current state: {state_summary}\n"
        "Data available: {data_info}\n\n"
        "Diagnose the issue and suggest solutions.\n"
        "Common issues:\n"
        "- Insufficient samples for WGCNA\n"
        "- Missing data matrices\n"
        "- No significant correlations found\n"
        "- Low candidate count\n"
        "- Contradictory evidence\n"
    ),
}
