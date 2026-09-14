"""
tool_handlers.py -- Tool execution handlers.

When the LLM calls a function, these handlers execute the actual logic
and return results formatted for the LLM to consume.
"""

from __future__ import annotations

import logging
from pathlib import Path
import re as _re
from typing import Any

from phyto_reason.models.workflow_plan import WorkflowPlan

logger = logging.getLogger("tool_handlers")


def handle_search_literature(arguments: dict) -> str:
    """Execute multi-source literature search (PubMed + Semantic Scholar + Europe PMC + RAG).

    All sources are queried in parallel. The local RAG document store is
    automatically included -- no need for a separate search_rag call.
    """
    query = arguments.get("query", "")
    max_results = min(arguments.get("max_results", 5), 10)

    if not query:
        return "Error: search query is required."

    # Resolve gene IDs in the query to expand with synonyms
    expanded_query = _expand_gene_ids_in_query(query)

    try:
        from phyto_reason.tools.literature.multi_source_search import (
            MultiSourceLiteratureSearch,
        )

        engine = MultiSourceLiteratureSearch()

        # Search with expanded query first
        result = engine.search(expanded_query, max_results=max_results)

        # If the expanded query returned nothing (all sources empty),
        # try the original query
        if "No results found across" in result and expanded_query != query:
            result = engine.search(query, max_results=max_results)

        return result

    except Exception as e:
        logger.error(f"Multi-source literature search failed: {e}")

        # Graceful degradation: fall back to PubMed only
        try:
            logger.info("Falling back to PubMed-only search")
            from phyto_reason.tools.literature.pubmed_client import PubMedSearch

            client = PubMedSearch()
            pm_results = client.search(query, max_results=max_results)
            if pm_results:
                lines = [
                    f"[Degraded -- PubMed only] Search results for '{query}':",
                    "",
                ]
                for i, r in enumerate(pm_results[:max_results], 1):
                    title = r.get("Title", r.get("title", "?"))
                    pmid = r.get("PMID", r.get("pmid", ""))
                    source = r.get("Source", r.get("source", ""))
                    date = r.get("PubDate", r.get("date", ""))
                    lines.append(f"{i}. {title[:200]}")
                    lines.append(f"   *{source}* ({date}) | PMID: {pmid}")
                    lines.append("")
                return "\n".join(lines)
            return f"Literature search unavailable. PubMed returned no results for: {query}"
        except Exception as e2:
            logger.error(f"PubMed fallback also failed: {e2}")
            return f"Literature search error: {e}"


def handle_query_kegg(arguments: dict) -> str:
    """Execute KEGG pathway query and return formatted results."""
    compound_name = arguments.get("compound_name", "")

    if not compound_name:
        return "Error: compound_name is required."

    try:
        from phyto_reason.knowledge.kegg_browser import KEGGBrowser

        browser = KEGGBrowser()
        result = browser.browse(compound_name)

        if not result.pathways:
            return result.note or f"No KEGG data found for '{compound_name}'"

        lines = [f"KEGG pathway information for '{compound_name}':", ""]
        for pw in result.pathways:
            lines.append(f"Pathway: {pw.pathway_name}")
            if pw.pathway_class:
                lines.append(f"Class: {pw.pathway_class}")
            if pw.enzyme_genes:
                lines.append(f"Enzyme genes ({len(pw.enzyme_genes)}): {', '.join(pw.enzyme_genes[:12])}")
            if pw.enzyme_names:
                lines.append(f"Enzymes: {', '.join(pw.enzyme_names[:5])}")
            if pw.source_url:
                lines.append(f"URL: {pw.source_url}")
            lines.append("")

        return "\n".join(lines)

    except Exception as e:
        logger.error(f"KEGG query failed: {e}")
        return f"KEGG query error: {e}"


def handle_query_plantcyc(arguments: dict) -> str:
    """Execute PlantCyc (PMN) query — plant-preferred pathway lookup."""
    from phyto_reason.tools.plantcyc_client import query_plantcyc

    try:
        result = query_plantcyc(
            metabolite_name=arguments.get("metabolite_name", ""),
            pathway_name=arguments.get("pathway_name", ""),
            enzyme_name=arguments.get("enzyme_name", ""),
            species=arguments.get("species", ""),
        )
        if not result.get("results"):
            return result.get("note") or "PlantCyc: 无匹配条目。"
        lines = [f"PlantCyc ({result.get('source', 'local_kb')}) — {result['n_results']} 条结果:", ""]
        for r in result["results"][:8]:
            kind = r.get("type", "?")
            lines.append(f"**{r.get('name', '')}** ({kind})")
            if r.get("pathway_id"):
                lines.append(f"   ID: {r['pathway_id']}")
            if r.get("key_enzymes"):
                lines.append(f"   关键酶: {', '.join(r['key_enzymes'][:8])}")
            if r.get("species"):
                lines.append(f"   物种: {r['species']}")
            if r.get("description"):
                lines.append(f"   {r['description'][:150]}")
            lines.append("")
        return "\n".join(lines)
    except Exception as e:
        logger.error(f"PlantCyc query failed: {e}")
        return f"PlantCyc query error: {e}"


def _format_annotation(result, label: str = "") -> str:
    """格式化注释结果为 LLM 可消费文本。"""
    from phyto_reason.metabolomics.annotator import AnnotationResult

    lines = []
    if label:
        lines.append(f"## {label}")
    lines.append(f"precursor {result.precursor_mz:.4f} → 加合物 {result.adduct}，"
                 f"中性质量 {result.neutral_mass:.4f} Da")
    valid = [f for f in result.formula_candidates if not f["violations"]]
    if valid:
        lines.append("分子式候选: " + ", ".join(
            f"{f['formula']} ({f['ppm']:.1f} ppm)" for f in valid[:5]
        ) + "（非唯一，需碎片/文献收敛）")
    else:
        lines.append("无通过 golden rules 的分子式候选")
    if result.candidates:
        lines.append(f"注释候选（{len(result.candidates)} 个，按分排序）:")
        for c in result.candidates:
            name = f" = {c.name}" if c.name else ""
            lines.append(
                f"  #{c.rank} [{c.class_name}]{name} — MSI Level {c.msi_level}, "
                f"score {c.score} ({c.rule_confidence})"
            )
            if c.formula:
                lines.append(f"     公式 {c.formula}; 通路 {c.pathway[:60]}")
            if c.matched_fragments:
                frag_str = ", ".join(
                    f"{m['observed_mz']:.4f}←{m['annotation']}" for m in c.matched_fragments[:4]
                )
                lines.append(f"     匹配碎片: {frag_str}")
    if result.profile_hits:
        lines.append("物种 profile 命中: " + ", ".join(
            f"{h['name']}（{h['ppm']} ppm, {h['source']}）" for h in result.profile_hits
        ))
    lines.append(result.note)
    return "\n".join(lines)


def handle_annotate_ms2_spectrum(arguments: dict) -> str:
    """MS² 谱图注释（未知物鉴定）——代谢组方向核心工具。"""
    from phyto_reason.metabolomics.spectrum_io import Spectrum, load_spectrum
    from phyto_reason.metabolomics.annotator import annotate_spectrum

    tolerance = float(arguments.get("tolerance_ppm", 10) or 10)
    species = arguments.get("species", "")
    target = arguments.get("target_metabolite", "")
    remote_query = arguments.get("remote_query", "") or target
    polarity = arguments.get("polarity", "")

    def annotate_one(spectrum, label: str = "") -> str:
        local_result = annotate_spectrum(
            spectrum, species=species, target_metabolite=target,
            tolerance_ppm=tolerance, polarity=polarity,
        )
        text = _format_annotation(local_result, label=label)
        if remote_query.strip():
            try:
                from phyto_reason.tools.spectra_tools import remote_annotation_fallback
                remote = remote_annotation_fallback(
                    remote_query.strip(), precursor_mz=float(spectrum.precursor_mz),
                    peaks=spectrum.fragments,
                )
                source = remote.metadata.get("source", "unknown")
                actual_source = remote.evidence_list[0].metadata.get("actual_source", "") if remote.evidence_list else ""
                text += (
                    f"\n\nRemote spectral comparison: source={source}"
                    f"{f' (actual={actual_source})' if actual_source else ''}; "
                    f"status={remote.status}; count={remote.metadata.get('count', 0)}"
                )
                for warning in remote.warnings:
                    text += f"\nWarning: {warning}"
            except Exception as exc:
                text += f"\nWarning: remote spectral comparison failed: {exc}"
        else:
            text += "\nRemote spectral comparison: skipped (no searchable compound name provided)."
        return text

    try:
        path = arguments.get("spectrum_path", "")
        if path:
            spectra = load_spectrum(path, polarity=polarity)
            if not spectra:
                return f"Error: 未从 {path} 解析到任何谱图"
            blocks = [
                annotate_one(
                    s,
                    label=f"谱图 {i + 1}（{s.scan_id or s.precursor_mz:.4f}）",
                )
                for i, s in enumerate(spectra[:5])
            ]
            if len(spectra) > 5:
                blocks.append(f"…（共 {len(spectra)} 张谱图，仅显示前 5 张）")
            return "\n\n".join(blocks)

        precursor_mz = arguments.get("precursor_mz", 0)
        fragments = arguments.get("fragments", [])
        if not precursor_mz or not fragments:
            return "Error: 需要 spectrum_path，或 (precursor_mz + fragments) 内联输入"
        spec = Spectrum(precursor_mz=float(precursor_mz),
                        fragments=[(float(f), 1.0) for f in fragments])
        return annotate_one(spec)

    except Exception as e:
        logger.error(f"MS2 annotation failed: {e}")
        return f"MS2 annotation error: {e}"


def _auto_detect_target_metabolite(metabolite_matrix: dict) -> str:
    """Auto-detect the most promising target metabolite from the data.

    Uses variance-based ranking: metabolites with the highest coefficient
    of variation across samples are likely differentially accumulated.
    Returns the name of the top-ranked metabolite.
    """
    import numpy as np

    best_metabolite = ""
    best_cv = 0.0

    try:
        for meta_name, sample_values in metabolite_matrix.items():
            if not sample_values or len(sample_values) < 2:
                continue

            # Convert to float array
            try:
                values = np.array([float(v) for v in sample_values.values()], dtype=np.float64)
            except (ValueError, TypeError):
                continue

            if len(values) < 2:
                continue

            mean_val = np.mean(values)
            if mean_val <= 0:
                continue  # skip near-zero metabolites

            std_val = np.std(values)
            cv = std_val / mean_val  # coefficient of variation

            if cv > best_cv:
                best_cv = cv
                best_metabolite = meta_name

    except Exception as e:
        logger.warning("Auto-detect target metabolite failed: %s", e)

    return best_metabolite


def handle_run_pipeline(
    arguments: dict,
    expression_matrix: dict | None = None,
    metabolite_matrix: dict | None = None,
    promoter_sequences: dict | None = None,
    species: str = "",
    target_metabolite: str = "",
    session=None,
    workflow_plan: WorkflowPlan | None = None,
) -> str:
    """Execute the scientific pipeline and return formatted results.

    If session is provided, stores structured candidate scores for downstream
    visualization (radar charts, network plots, etc.).
    """
    species = arguments.get("species", species)
    target_metabolite = arguments.get("target_metabolite", target_metabolite)
    target_pathway = arguments.get("target_pathway", "")

    # Check if data is available
    has_data = bool(expression_matrix)

    if not has_data:
        # Check if session has data in DB but not passed to this call
        session_has_data = False
        session_id = ""
        if session is not None:
            session_expr = getattr(session, "expression_matrix", None)
            session_has_data = bool(session_expr)
            session_id = getattr(session, "session_id", "")

        if session_has_data:
            # Data exists in session but wasn't passed to handler -- retry with session data
            logger.warning(
                "Pipeline called without expression_matrix but session %s has %d genes. "
                "Re-reading from session.",
                session_id, len(session_expr) if session_expr else 0,
            )
            expression_matrix = session_expr
            metabolite_matrix = getattr(session, "metabolite_matrix", None) or metabolite_matrix
            has_data = True
        else:
            return (
                f"⚠️ 没有可用的上传数据。\n\n"
                f"当前会话: {session_id}\n"
                f"表达矩阵: {'已加载' if expression_matrix else '未加载 (None)'}\n"
                f"代谢物矩阵: {'已加载' if metabolite_matrix else '未加载 (None)'}\n\n"
                f"请通过以下步骤上传数据:\n"
                f"1. 在左侧面板拖拽或点击上传 CSV/TSV/XLSX 表达矩阵文件\n"
                f"2. 上传代谢物矩阵文件（文件名包含 'metabolite'）\n"
                f"3. 输入物种名称和目标代谢物\n"
                f"4. 上传完成后，再发送分析请求\n\n"
                f"如果已经上传过数据，请检查:\n"
                f"- 是否刷新了页面（刷新会创建新会话，之前上传的数据无法访问）\n"
                f"- 文件格式是否正确（CSV/TSV/XLSX，基因行 x 样本列）\n"
                f"- 上传响应是否显示 '✅ done'\n\n"
                f"⚠️ No uploaded data available.\n"
                f"Please upload expression + metabolite data files first, "
                f"then request analysis. If you already uploaded data, "
                f"do NOT refresh the page -- refreshing creates a new session."
            )

    # Auto-detect target metabolite from the metabolite matrix if not specified
    if not target_metabolite and metabolite_matrix:
        target_metabolite = _auto_detect_target_metabolite(metabolite_matrix)
        if target_metabolite:
            logger.info("Auto-detected target metabolite: %s", target_metabolite)
        else:
            return (
                f"Data is available ({len(expression_matrix)} genes, {len(metabolite_matrix)} metabolites) "
                f"but no target metabolite was specified. Please tell me which metabolite or "
                f"compound class you are interested in. For example: "
                f"'analyze alkaloids' or 'find regulators of nitidine'."
            )

    try:
        from phyto_reason.workflows.workflow_runner import WorkflowRunner

        runner = WorkflowRunner()

        # Determine sample count from data
        sample_count = 0
        if expression_matrix:
            first_key = next(iter(expression_matrix), None)
            if first_key:
                sample_count = len(expression_matrix[first_key])

        state = runner.run(
            species=species,
            target_metabolite=target_metabolite,
            target_pathway=target_pathway,
            sample_count=sample_count,
            has_expression=bool(expression_matrix),
            has_metabolite=bool(metabolite_matrix),
            has_promoter=bool(promoter_sequences),
            expression_matrix=expression_matrix,
            metabolite_matrix=metabolite_matrix,
            promoter_sequences=promoter_sequences,
            workflow_plan=workflow_plan,
        )

        # Extract structured candidate scores for visualization
        if session is not None:
            _store_candidate_scores(state, session)

        # 报告支持：分节摘要 + 图表编号注册（仅覆盖实际完成的节点）
        if session is not None:
            try:
                from phyto_reason.reports.analysis_report import build_analysis_digest
                session.analysis_digest = build_analysis_digest(
                    state,
                    species=getattr(getattr(state, "planner_state", None), "species", "") or getattr(session, "species", ""),
                    target_metabolite=getattr(getattr(state, "planner_state", None), "target_metabolite", "") or getattr(session, "target_metabolite", ""),
                    figure_registry=getattr(session, "figure_captions", None),
                )
            except Exception as exc:
                logger.warning("analysis digest build failed: %s", exc)

        # Auto-refresh RAG index so new results are searchable
        try:
            from phyto_reason.tools.rag import refresh_rag_engine
            n = refresh_rag_engine()
            logger.info("RAG index refreshed after pipeline: %d chunks", n)
        except Exception:
            pass  # non-critical

        # Detect language from session history
        language = "cn"  # default Chinese
        if session is not None:
            from phyto_reason.agent.orchestrator import _detect_chinese
            all_user_msgs = " ".join(
                t.get("content", "") for t in (getattr(session, "history", []) or [])
                if t.get("role") == "user"
            )
            if not _detect_chinese(all_user_msgs):
                language = "en"

        result_text = _format_pipeline_result(state, language=language)
        if session is not None:
            # 假设池数据源：把竞争假设持久化到会话
            try:
                hypotheses = getattr(state, "competing_hypotheses", None) or []
                session.set_pipeline_result(
                    result_text, [h.model_dump() for h in hypotheses]
                )
            except Exception as exc:
                logger.warning("hypothesis persistence failed: %s", exc)
        return result_text

    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}", exc_info=True)
        return f"Pipeline error: {e}"


def _store_candidate_scores(state, session) -> None:

    """Produces: session.candidate_scores = {
        "WRKY1": {"correlation": 0.85, "motif": 0.6, "pathway": 0.7, ...},
        "MYB2": {"correlation": 0.72, ...},
        ...
    }
    """
    try:
        fusion = getattr(state, "fusion_results", None) or {}
        if not fusion:
            logger.info("No fusion_results in pipeline state -- skipping candidate scores")
            return

        scored: dict[str, dict[str, float]] = {}
        for gene_id, result in fusion.items():
            gene_scores: dict[str, float] = {}
            # Collect scores from individual evidence items
            for ev in result.evidence_list:
                source = ev.source.lower().replace(" ", "_")
                gene_scores[source] = round(ev.normalized_score, 4)
            # Add the overall calibrated score
            gene_scores["overall"] = round(result.calibrated_score, 4)
            scored[gene_id] = gene_scores

        session.candidate_scores = scored
        logger.info(
            "Stored candidate scores for %d genes (dims: %s)",
            len(scored),
            list(next(iter(scored.values())).keys()) if scored else "none",
        )

        # Also load DEG/DAM data for volcano plots
        _load_deg_dam_data(session, state=state)

    except Exception as e:
        logger.warning("Failed to extract candidate scores: %s", e)


def _load_deg_dam_data(session, state=None) -> None:
    """Load DEG and DAM data from pipeline state or CSV files into session.

    Priority: pipeline state.deg_report/dam_report over CSV files (outputs/step1/).
    Stores structured dicts in session.deg_data and session.dam_data for volcano plots.
    """
    # ── Try pipeline state first (modern LangGraph pipeline) ──
    if state is not None:
        deg_report = getattr(state, "deg_report", None)
        if deg_report and deg_report.get("top_genes"):
            deg_dict: dict[str, dict[str, float]] = {}
            for g in deg_report["top_genes"]:
                gid = g.get("gene_id", "")
                if gid:
                    deg_dict[gid] = {
                        "log2fc": g.get("max_log2fc", 0),
                        "padj": g.get("q_value", g.get("p_value", 1)),
                    }
            session.deg_data = deg_dict
            logger.info("Loaded DEG data from pipeline state: %d genes", len(deg_dict))

        dam_report = getattr(state, "dam_report", None)
        if dam_report and dam_report.get("top_metabolites"):
            dam_dict: dict[str, dict[str, float]] = {}
            for m in dam_report["top_metabolites"]:
                name = m.get("metabolite", "")
                if name:
                    dam_dict[name] = {
                        "log2fc": m.get("max_log2fc", 0),
                        "pvalue": m.get("q_value", m.get("p_value", 1)),
                    }
            session.dam_data = dam_dict
            logger.info("Loaded DAM data from pipeline state: %d metabolites", len(dam_dict))

        if session.deg_data or session.dam_data:
            return  # State data is sufficient

    # ── Fallback: CSV files (old pipeline) ──
    from pathlib import Path
    import pandas as pd

    project_root = Path(__file__).resolve().parent.parent.parent
    step1_dir = project_root / "outputs" / "step1"

    # DEG data
    deg_path = step1_dir / "DEGs_significant.csv"
    if deg_path.exists():
        try:
            df = pd.read_csv(deg_path)
            deg_dict: dict[str, dict[str, float]] = {}
            for _, row in df.iterrows():
                gid = str(row.get("GeneID", row.get("gene_id", "")))
                if not gid or gid == "nan":
                    continue
                lfc = float(row.get("log2FoldChange", row.get("log2fc", 0)))
                padj = float(row.get("padj", row.get("pvalue", 1)))
                deg_dict[gid] = {"log2fc": lfc, "padj": padj}
            session.deg_data = deg_dict
            logger.info("Loaded DEG data: %d genes from %s", len(deg_dict), deg_path)
        except Exception as e:
            logger.warning("Failed to load DEG data: %s", e)

    # ── DAM data ────────────────────────────────────────
    dam_path = step1_dir / "DAMs_significant.csv"
    if dam_path.exists():
        try:
            df = pd.read_csv(dam_path)
            dam_dict: dict[str, dict[str, float]] = {}
            for _, row in df.iterrows():
                name = str(row.get("Name", row.get("name", row.get("Metabolite", ""))))
                if not name or name == "nan":
                    continue
                lfc = float(row.get("log2FC", row.get("log2fc", 0)))
                pval = float(row.get("pvalue", row.get("padj", 1)))
                dam_dict[name] = {"log2fc": lfc, "pvalue": pval}
            session.dam_data = dam_dict
            logger.info("Loaded DAM data: %d metabolites from %s", len(dam_dict), dam_path)
        except Exception as e:
            logger.warning("Failed to load DAM data: %s", e)


def handle_search_knowledge_base(arguments: dict) -> str:
    """Query the built-in TF-metabolite knowledge base."""
    metabolite = arguments.get("metabolite", "")
    tf_family = arguments.get("tf_family", "")
    species = arguments.get("species", "")

    if not metabolite:
        return "Error: metabolite parameter is required."

    try:
        from phyto_reason.knowledge.tf_knowledge_base import TFKnowledgeBase
        from phyto_reason.knowledge.species_registry import get_species_registry

        # ── 物种知识段落（降级链）───────────────────────
        species_lines: list[str] = []
        if species:
            sk_result, level = get_species_registry().knowledge_for(
                species, metabolite=metabolite, tf_family=tf_family
            )
            species_lines.append(f"## 物种知识: {species}")
            if level == "profile":
                profile = get_species_registry().get(species)
                species_lines.append(f"知识层级: 本地 profile（置信度正常）")
                for m in profile.marker_metabolites:
                    mz = f" [M+H]+ {m.mz:.4f}" if m.mz else ""
                    species_lines.append(f"- {m.name}{mz}（来源 {m.source}）")
                if profile.pathway_prior:
                    for pw, info in profile.pathway_prior.items():
                        enz = ", ".join(info.get("enzymes", []) or []) or "未标注"
                        species_lines.append(f"- 通路 {pw}: 酶 [{enz}]")
                if profile.refs:
                    species_lines.append(f"- 参考: {', '.join(profile.refs[:2])}")
            elif level == "ortholog":
                species_lines.append("知识层级: 无本地 profile → 同源迁移（置信度降级）")
            else:
                species_lines.append("知识层级: 显式知识缺口（无公共数据）")
            species_lines.append("")

        kb = TFKnowledgeBase()
        result = kb.query(metabolite=metabolite, tf_family=tf_family, species=species)

        # ── MS/MS 诊断规则段落（compound_profiles）───────
        from phyto_reason.knowledge.compound_profiles import find_diagnostic_rules
        diag = find_diagnostic_rules(metabolite)
        diag_lines: list[str] = []
        if diag:
            diag_lines.append(f"## MS/MS 诊断规则（{diag['class_file']}）")
            diag_lines.append(f"结构类: {diag['class']}")
            if diag.get("formula") and diag.get("mass"):
                diag_lines.append(f"母核: {diag['formula']} (monoisotopic {diag['mass']:.4f})")
            if diag.get("adducts"):
                diag_lines.append(f"加合物: {', '.join(diag['adducts'])}")
            if diag.get("fragments"):
                diag_lines.append("诊断碎片:")
                for fr in diag["fragments"][:6]:
                    mz = fr.get("mz", "?")
                    ann = fr.get("annotation", "")
                    formula = fr.get("formula", "")
                    diag_lines.append(f"  - {mz:.4f} ({formula}) {ann}")
            diag_lines.append("")

        if not result.relations:
            return ("\n".join(species_lines + diag_lines)
                    + (result.note or f"No known TF regulators found for '{metabolite}'."))

        lines = species_lines + diag_lines + [f"Known TF regulators for '{metabolite}':", ""]

        # Group by TF family
        families: dict[str, list] = {}
        for r in result.relations:
            families.setdefault(r.tf_family, []).append(r)

        for family, rels in sorted(families.items(), key=lambda x: -max(r.score for r in x[1])):
            best = max(rels, key=lambda r: r.score)
            strength_emoji = {"strong": "🟢", "moderate": "🟡", "weak": "🔴", "inferred": "⚪"}
            emoji = strength_emoji.get(best.strength, "⚪")

            lines.append(f"{emoji} **{family}** -- {best.strength} ({best.score:.0%})")
            lines.append(f"   {best.description[:200]}")
            if best.known_examples:
                lines.append(f"   Examples: {', '.join(best.known_examples[:5])}")
            if best.pmids:
                lines.append(f"   PMIDs: {', '.join(best.pmids[:3])}")
            if best.species_scope != "general":
                lines.append(f"   Scope: {best.species_scope}")
            lines.append("")

        lines.append(f"---")
        lines.append(f"Found {result.total} relationships.")
        if result.note:
            lines.append(result.note)

        return "\n".join(lines)

    except Exception as e:
        logger.error(f"Knowledge base query failed: {e}")
        return f"Knowledge base query error: {e}"


def handle_cross_species_infer(arguments: dict) -> str:
    """Cross-species regulatory inference -- with real ortholog computation.

    Two modes:
      1. gene_id provided -> OrthologFinder (NCBI E-utilities + sequence alignment)
      2. no gene_id -> SpeciesKnowledge (curated TF family rules, backward compatible)
    """
    target_species = arguments.get("target_species", "")
    metabolite = arguments.get("metabolite", "")
    tf_family = arguments.get("tf_family", "")
    source_species = arguments.get("source_species", "")
    gene_id = arguments.get("gene_id", "")

    if not target_species:
        return "Error: target_species parameter is required."

    try:
        # ── Mode 1: Gene-level ortholog computation ─────────
        if gene_id:
            return _handle_ortholog_lookup(
                gene_id=gene_id,
                source_species=source_species,
                target_species=target_species,
                metabolite=metabolite,
                tf_family=tf_family,
            )

        # ── Mode 2: TF-family level（走物种降级链）────────
        # 降级链: 有 profile → 直接答；无 → 同源迁移；都无 → 显式缺口
        from phyto_reason.knowledge.species_registry import get_species_registry

        if source_species:
            # 显式 source → target 模式：沿用 cross_species（profile 优先）
            from phyto_reason.knowledge.species_knowledge import SpeciesKnowledge
            registry = get_species_registry()
            target_profile = registry.get(target_species)
            sk = SpeciesKnowledge()
            result = sk.cross_species(
                source_species=source_species,
                target_species=target_species,
                tf_family=tf_family,
                metabolite=metabolite,
            )
            level = "profile" if target_profile else "ortholog"
        else:
            result, level = get_species_registry().knowledge_for(
                target_species,
                metabolite=metabolite,
                tf_family=tf_family,
            )

        return _format_species_knowledge_result(result, target_species, level=level)

    except Exception as e:
        logger.error(f"Cross-species inference failed: {e}")
        return f"Cross-species inference error: {e}"


# ═══════════════════════════════════════════════════════════════
# Gene-level ortholog lookup (new)
# ═══════════════════════════════════════════════════════════════

def _handle_ortholog_lookup(
    gene_id: str,
    source_species: str,
    target_species: str,
    metabolite: str = "",
    tf_family: str = "",
) -> str:
    """Query NCBI for orthologs of gene_id in target_species."""
    if not source_species:
        source_species = "Arabidopsis thaliana"

    # Resolve gene ID to get primary symbol and NCBI Gene ID
    resolved_gene_id = gene_id
    gene_symbol = ""
    resolved_tf_family = tf_family
    try:
        from phyto_reason.knowledge.gene_resolver import resolve_gene
        gene_info = resolve_gene(gene_id, species=source_species)
        if gene_info and gene_info.primary_symbol and gene_info.primary_symbol != gene_id:
            gene_symbol = gene_info.primary_symbol
            if gene_info.ncbi_gene_id:
                resolved_gene_id = gene_info.ncbi_gene_id
            if gene_info.tf_family and not tf_family:
                resolved_tf_family = gene_info.tf_family
            logger.info(
                "Gene resolved: %s -> %s (NCBI:%s, TF:%s)",
                gene_id, gene_symbol, resolved_gene_id, resolved_tf_family,
            )
    except Exception as e:
        logger.debug("GeneResolver not available for ortholog lookup: %s", e)

    try:
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder

        finder = OrthologFinder()
        result = finder.find_orthologs(
            gene_id=resolved_gene_id,
            source_species=source_species,
            target_species=target_species,
        )

        return _format_ortholog_result(result, metabolite, resolved_tf_family)

    except Exception as e:
        logger.warning(
            "OrthologFinder failed, falling back to SpeciesKnowledge: %s", e
        )
        # Graceful degradation: fall back to curated knowledge
        from phyto_reason.knowledge.species_knowledge import SpeciesKnowledge

        sk = SpeciesKnowledge()
        sk_result = sk.infer_regulation(
            target_species=target_species,
            metabolite=metabolite,
            tf_family=resolved_tf_family,
        )
        return (
            f"[DEGRADED] NCBI ortholog lookup unavailable ({e}). "
            f"Using curated knowledge instead.\n\n"
            + _format_species_knowledge_result(sk_result, target_species)
        )


def _enrich_with_public_expression(lines: list[str], result) -> None:
    """Enrich ortholog result with public expression data (BAR/eFP).

    Fetches tissue expression for the top 2 ortholog hits if they are
    Arabidopsis genes. Runs with a 5-second timeout to avoid blocking.
    """
    if not result.hits:
        return

    try:
        from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutTimeout

        top_hits = result.hits[:2]

        def _fetch_context():
            """Fetch BAR context for ortholog hits."""
            from phyto_reason.knowledge.public_expression_client import get_bar_client
            bar = get_bar_client()
            contexts = []
            for hit in top_hits:
                ctx = bar.get_gene_context(hit.target_gene_id)
                contexts.append((hit, ctx))
            return contexts

        try:
            with ThreadPoolExecutor(max_workers=1) as pool:
                future = pool.submit(_fetch_context)
                contexts = future.result(timeout=5.0)
        except (FutTimeout, Exception):
            return  # BAR unavailable, skip enrichment

        if not contexts:
            return

        lines.append("── Public expression context (BAR/eFP) ──")
        lines.append("")

        for hit, ctx in contexts:
            symbol = ctx.get("symbol", "") or hit.target_gene_symbol
            top_tissues = ctx.get("top_tissues", [])
            expressolog_count = ctx.get("expressolog_count", 0)

            if top_tissues:
                tissue_str = ", ".join(top_tissues[:5])
                lines.append(
                    f"  {symbol}: expressed in {tissue_str}"
                )
                if expressolog_count > 0:
                    lines.append(
                        f"    -> {expressolog_count} expressolog(s) found "
                        f"(genes with correlated expression patterns)"
                    )
            elif expressolog_count > 0:
                lines.append(
                    f"  {symbol}: {expressolog_count} expressolog(s) found"
                )

        lines.append("")
        lines.append("  Source: BAR (bar.utoronto.ca) -- Sullivan et al. 2024, PMID: 39441075")
        lines.append("")

    except ImportError:
        pass  # public_expression_client not available
    except Exception:
        pass  # graceful degradation


def _add_motif_conservation_context(
    lines: list[str], result, tf_family: str,
) -> None:
    """Add curated motif conservation context for cross-species TF inference."""
    try:
        from phyto_reason.knowledge.motif_conservation import get_family_conservation

        rule = get_family_conservation(tf_family)
        if not rule:
            return

        if rule.conserved_in_plants and rule.motif_divergence == "low":
            lines.append("")
            lines.append(
                f"  🔬 Motif note: {tf_family} binding motifs are highly conserved "
                f"across plant species ({rule.motif_divergence} divergence). "
                f"This supports functional conservation of the ortholog's regulatory role."
            )
        elif rule.conserved_in_plants:
            lines.append("")
            lines.append(
                f"  🔬 Motif note: {tf_family} binding motifs show "
                f"{rule.motif_divergence} divergence across species. "
                f"Regulatory function may be partially conserved."
            )

    except ImportError:
        pass
    except Exception:
        pass


def _format_ortholog_result(
    result,
    metabolite: str = "",
    tf_family: str = "",
) -> str:
    """Format OrthologResult into LLM-consumable text."""
    lines = [
        f"Ortholog analysis for {result.query_gene_id} "
        f"from {result.source_species} -> {result.target_species}:",
        "",
    ]

    # Method note
    method_labels = {
        "ncbi": "NCBI E-utilities (gene-level ortholog search)",
        "mygene_domain": "MyGene.info (protein domain conservation -- functional orthologs)",
        "fallback_curated": "Curated knowledge (NCBI unavailable)",
    }
    lines.append(f"Method: {method_labels.get(result.method_used, result.method_used)}")
    lines.append("")

    # Hits
    if result.hits:
        lines.append(f"Orthologs found: {len(result.hits)}")
        lines.append("")
        for i, hit in enumerate(result.hits, 1):
            conf_icon = {"high": "🟢", "medium": "🟡", "low": "🔴"}.get(
                hit.confidence, "⚪"
            )
            lines.append(
                f"  {i}. {hit.target_gene_symbol} (NCBI Gene ID: {hit.target_gene_id})"
            )
            lines.append(f"     Source gene: {hit.source_gene_symbol} ({hit.source_gene_id})")
            lines.append(f"     Identity: {hit.percent_identity:.1f}% "
                         f"(aligned: {hit.alignment_length} aa)")
            lines.append(f"     Method: {hit.method.replace('_', ' ')}")
            lines.append(f"     Confidence: {conf_icon} {hit.confidence.upper()}")
            if hit.description:
                lines.append(f"     Description: {hit.description[:150]}")
            if hit.source_protein_id:
                lines.append(f"     Source protein: {hit.source_protein_id}")
            if hit.target_protein_id:
                lines.append(f"     Target protein: {hit.target_protein_id}")
            lines.append("")

        # ── Enrich with public expression data (BAR) ────────
        _enrich_with_public_expression(lines, result)
    else:
        lines.append("No orthologs found in the target species.")
        lines.append("")

    # Regulatory inference context
    if result.hits and (metabolite or tf_family):
        lines.append("Regulatory inference:")
        top_hit = result.hits[0]
        if tf_family and metabolite:
            lines.append(
                f"  The top ortholog ({top_hit.target_gene_symbol}, "
                f"{top_hit.percent_identity:.1f}% identity) may have conserved "
                f"{tf_family}-mediated regulation of {metabolite} biosynthesis. "
            )
        elif metabolite:
            lines.append(
                f"  The ortholog may be involved in {metabolite} regulation "
                f"in {result.target_species}."
            )
        lines.append("")

    # Caveats
    lines.append("Caveats:")
    lines.append("  - Ortholog evidence supports sequence conservation, but functional")
    lines.append("    conservation requires experimental validation.")
    lines.append("  - Check tissue expression patterns and promoter motif conservation.")
    lines.append("  - Cross-species inferences should be treated as hypotheses.")

    # ── Motif conservation context (curated) ──────────────
    if tf_family and result.hits:
        _add_motif_conservation_context(lines, result, tf_family)

    lines.append("")

    # Note
    if result.note:
        lines.append(f"Note: {result.note}")

    return "\n".join(lines)


def _format_species_knowledge_result(result, target_species: str, level: str | None = None) -> str:
    """Format SpeciesKnowledgeResult into LLM-consumable text.

    Args:
        level: 降级链层级 "profile" / "ortholog" / "gap"（来自 SpeciesRegistry）。
            带层级时输出置信度标注，驱动 LLM 措辞降级。
    """
    lines = [f"Cross-species analysis for '{target_species}':", ""]

    if level:
        badge = {
            "profile": "本地知识 profile（置信度正常）",
            "ortholog": "同源迁移推断（置信度降级）",
            "gap": "显式知识缺口",
        }
        lines.insert(0, f"知识层级: {badge.get(level, level)}")

    if result.known_metabolites:
        lines.append(f"Known metabolites: {', '.join(result.known_metabolites[:8])}")
    if result.known_pathways:
        lines.append(f"Known pathways: {', '.join(result.known_pathways)}")
    lines.append("")

    if result.ortholog_inferences:
        lines.append(f"Regulatory inferences ({len(result.ortholog_inferences)}):")
        for inf in result.ortholog_inferences:
            conf_emoji = {
                "high": "🟢", "medium": "🟡", "low": "🔴", "speculative": "⚪",
            }
            lines.append(
                f"  {conf_emoji.get(inf.confidence, '⚪')} "
                f"{inf.tf_family} -> {inf.metabolite}: {inf.confidence} confidence"
            )
            lines.append(f"     {inf.reasoning[:150]}")
            if inf.caveats:
                for c in inf.caveats[:2]:
                    lines.append(f"     ⚠ {c}")
    else:
        lines.append("No specific regulatory inferences available.")

    lines.append("")
    lines.append(result.note)

    return "\n".join(lines)


def handle_check_quality(arguments: dict, data: dict | None = None) -> str:
    """Check data quality and return report."""
    session_id = arguments.get("session_id", "")

    if not data or not data.get("expression"):
        return "No data uploaded. Please upload expression (CSV) and metabolite data first."

    try:
        from phyto_reason.ingestion.validators.matrix_integrity_validator import (
            validate_matrix_integrity,
        )
        from phyto_reason.ingestion.validators.missing_value_validator import (
            validate_missing_values,
        )

        expr = data.get("expression") or {}
        meta = data.get("metabolite") or {}

        qc_expr = validate_matrix_integrity(expr, name="expression")
        mv = validate_missing_values(expr)

        lines = ["Data Quality Report:", ""]
        lines.append(f"Genes: {len(expr)}")
        lines.append(f"Metabolites: {len(meta)}")

        if expr:
            first = next(iter(expr.values()), {})
            lines.append(f"Samples: {len(first)}")

        lines.append(f"Integrity: {'PASS' if qc_expr.all_passed else 'WARNING'}")
        lines.append(f"Missing values: {'PASS' if mv.passed else f'WARNING -- {mv.warnings[:2]}'}")

        return "\n".join(lines)

    except Exception as e:
        logger.error(f"Quality check failed: {e}")
        return f"Quality check error: {e}"


def handle_search_rag(arguments: dict) -> str:
    """Search the local RAG document store and return formatted results."""
    query = arguments.get("query", "")
    top_k = min(arguments.get("top_k", 5), 10)

    if not query:
        return "Error: query is required for RAG search."

    try:
        from phyto_reason.tools.rag import get_rag_engine

        engine = get_rag_engine()
        results = engine.search(query, top_k=top_k)

        if not results:
            return (
                f"No relevant documents found in the local store for: '{query}'.\n"
                f"The document index may be empty. Run a scientific pipeline or "
                f"upload data to populate it with analysis results."
            )

        context = engine.format_context(results, max_chars=4000)
        lines = [
            f"RAG search results for '{query}' ({len(results)} matches):",
            "",
            context,
        ]
        return "\n".join(lines)

    except Exception as e:
        logger.error(f"RAG search failed: {e}")
        return f"RAG search error: {e}"


def handle_web_search(arguments: dict) -> str:
    """Search the web using DuckDuckGo and return formatted results."""
    query = arguments.get("query", "")
    max_results = min(arguments.get("max_results", 8), 10)

    if not query:
        return "Error: query is required for web search."

    try:
        from phyto_reason.tools.web_search import web_search
        return web_search(query, max_results=max_results)
    except Exception as e:
        logger.error(f"Web search failed: {e}")
        return f"Web search error: {e}"


def handle_recall_memory(arguments: dict) -> str:
    """Search the agent memory stream and return formatted results."""
    query = arguments.get("query", "")
    top_k = min(arguments.get("top_k", 5), 10)

    if not query:
        return "Error: query is required for memory recall."

    try:
        from phyto_reason.agent.memory_store import memory_store

        results = memory_store.search(query, top_k=top_k)
        total = memory_store.count()

        if not results:
            return (
                f"No relevant memories found for: '{query}'.\n"
                f"The memory stream contains {total} total memories. "
                f"As we continue working together, I'll remember key findings and preferences."
            )

        lines = [f"Memory recall for '{query}' ({len(results)} of {total} total memories):", ""]

        type_emoji = {
            "observation": "📝", "reflection": "💡",
            "preference": "⭐", "result": "📊",
        }

        for i, r in enumerate(results, 1):
            emoji = type_emoji.get(r.get("memory_type", "observation"), "📝")
            lines.append(
                f"{i}. {emoji} [{r['memory_type']}] (score={r['score']:.3f}, "
                f"importance={r['importance']}/10)"
            )
            lines.append(f"   {r['content'][:300]}")
            if r.get("species"):
                lines.append(f"   Species: {r['species']}")
            if r.get("metabolite"):
                lines.append(f"   Metabolite: {r['metabolite']}")
            if r.get("tags"):
                lines.append(f"   Tags: {', '.join(r['tags'])}")
            lines.append(f"   Created: {r['created_at'][:19]}")
            lines.append("")

        return "\n".join(lines)

    except Exception as e:
        logger.error(f"Memory recall failed: {e}")
        return f"Memory recall error: {e}"


def handle_query_public_expression(arguments: dict) -> str:
    """Query public expression databases (BAR/eFP) for gene expression context.

    Use when: you need tissue-specific expression data, expressologs,
    or expression context for ortholog candidates.
    """
    gene_id = arguments.get("gene_id", "")
    if not gene_id:
        return "Error: gene_id is required (e.g., AT1G56650 for Arabidopsis)."

    try:
        from phyto_reason.knowledge.public_expression_client import get_bar_client

        bar = get_bar_client()
        ctx = bar.get_gene_context(gene_id)

        lines = [
            f"Public expression context for {gene_id}:",
            "",
        ]

        if ctx.get("symbol"):
            lines.append(f"Symbol: {ctx['symbol']}")
        if ctx.get("annotation"):
            lines.append(f"Annotation: {ctx['annotation'][:200]}")
        lines.append("")

        if ctx.get("top_tissues"):
            lines.append(f"Top expressed tissues: {', '.join(ctx['top_tissues'][:10])}")
            lines.append(f"Total tissues with data: {ctx['tissue_count']}")
        else:
            lines.append("No tissue expression data available from BAR.")

        lines.append("")

        if ctx.get("expressolog_count", 0) > 0:
            lines.append(f"Expressologs found: {ctx['expressolog_count']}")
            for el in ctx.get("expressologs", [])[:5]:
                lines.append(
                    f"  - {el['gene']} ({el['species']}): "
                    f"correlation={el['correlation']:.3f}"
                )
        else:
            lines.append("No expressologs found.")

        lines.append("")

        if ctx.get("ppi_count", 0) > 0:
            lines.append(f"Known protein interactions: {ctx['ppi_count']}")

        lines.append("")
        lines.append("Source: BAR (bar.utoronto.ca)")
        lines.append("Reference: Sullivan et al. (2024), PMID: 39441075")

        return "\n".join(lines)

    except ImportError:
        return "Public expression database client not available (requires bar_client module)."
    except Exception as e:
        logger.warning("Public expression query failed: %s", e)
        return f"Public expression query failed: {e}"


def _format_pipeline_result(state, language: str = "cn") -> str:
    """Format pipeline results for LLM consumption.

    Returns a compact but complete summary of hypotheses, evidence, and caveats.
    Language-adaptive: "cn" for Chinese, "en" for English.
    """
    hypotheses = getattr(state, "competing_hypotheses", None) or []
    mechanism_type = getattr(state, "mechanism_type", "") or "unknown"
    pathway = getattr(state, "pathway_name", "") or "unknown"
    n_samples = state.planner_state.sample_count if state.planner_state else 0
    quality_passed = getattr(state, "data_quality_passed", False)
    is_cn = language == "cn"

    # Degradation notes
    degraded = getattr(state, "degraded_tools", []) or []
    deg_notes = getattr(state, "degradation_notes", []) or []

    # Bilingual labels
    if is_cn:
        L = {
            "quality_fail": "数据质量检查未通过",
            "degraded_tools": "降级工具",
            "target": "目标代谢物",
            "pathway": "通路",
            "mechanism": "机制类型",
            "tf_found": "发现候选 TF 数",
            "score": "得分",
            "confidence": "置信度",
            "no_results": "管线已完成但未找到假说或候选基因。数据可能不足。",
            "primary_mechanism": "主要机制",
            "samples": "样本数",
            "hypotheses": "生成假说数",
            "type": "类型",
            "uncertainty": "不确定性",
            "candidate_regulators": "候选调控因子",
            "supporting": "支持证据",
            "contradictory": "矛盾证据",
            "missing": "缺失证据",
            "why_wrong": "为什么可能错误",
            "falsification": "证伪测试",
            "test": "测试",
            "if_correct": "若成立",
            "if_wrong": "若不成立",
            "full_report": "完整报告",
            "exec_progress": "执行进度",
            "degradation": "降级说明",
            "degraded_tools_list": "失败或降级的工具",
            "degraded_note": "因部分工具降级，分析置信度可能降低。",
            "unknown": "未知",
            "mechanism_labels": {
                "transcriptional_regulation": "转录调控",
                "transport_redistribution": "转运重分布",
                "stress_induced_redistribution": "胁迫诱导重分布",
                "enzymatic_regulation": "酶促调控",
                "unknown": "未知",
            },
        }
    else:
        L = {
            "quality_fail": "Data quality check failed",
            "degraded_tools": "Degraded tools",
            "target": "Target",
            "pathway": "Pathway",
            "mechanism": "Mechanism type",
            "tf_found": "TF candidates found",
            "score": "score",
            "confidence": "confidence",
            "no_results": "Pipeline completed but no hypotheses or candidates were found. Data may be insufficient.",
            "primary_mechanism": "Primary mechanism",
            "samples": "Samples",
            "hypotheses": "Hypotheses generated",
            "type": "Type",
            "uncertainty": "Uncertainty",
            "candidate_regulators": "Candidate regulators",
            "supporting": "Supporting",
            "contradictory": "Contradictory",
            "missing": "Missing evidence",
            "why_wrong": "Why might be wrong",
            "falsification": "Falsification tests",
            "test": "Test",
            "if_correct": "If correct",
            "if_wrong": "If wrong",
            "full_report": "Full Report",
            "exec_progress": "Execution Progress",
            "degradation": "Degradation Notes",
            "degraded_tools_list": "Tools that failed or were degraded",
            "degraded_note": "Confidence may be reduced due to missing evidence from degraded tools.",
            "unknown": "unknown",
            "mechanism_labels": {},
        }

    mechanism_labels = L["mechanism_labels"]
    mechanism_label = mechanism_labels.get(mechanism_type, mechanism_type) if is_cn else mechanism_type

    if not quality_passed:
        issues = getattr(state, "data_quality_issues", [])
        result = f"{L['quality_fail']}: {'; '.join(issues[:3])}"
        if degraded:
            result += f"\n{L['degraded_tools']}: {', '.join(degraded)}"
        return result

    if not hypotheses:
        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}

        # QC figures section
        lines = []
        qc_report = getattr(state, "qc_report", None) or {}
        qc_figs = qc_report.get("figure_markdown", "")
        if qc_figs and qc_report.get("n_samples", 0) > 0:
            label_qc = "数据质控" if is_cn else "Data Quality Control"
            lines.append(f"## {label_qc}")
            lines.append(f"  样本数: {qc_report['n_samples']} | 基因数: {qc_report.get('n_genes', 0)} | 组数: {len(qc_report.get('groups', []))}")
            lines.append("")
            lines.append(qc_figs)
            lines.append("")

        label_target = "目标代谢物" if is_cn else "Target"
        label_pathway = "通路" if is_cn else "Pathway"
        label_mechanism = "机制类型" if is_cn else "Mechanism"
        label_deg = "差异表达基因 (DEG)" if is_cn else "DEG"
        label_dam = "差异代谢物 (DAM)" if is_cn else "DAM"
        label_tested = "检验基因数" if is_cn else "Genes tested"
        label_fdr = "FDR-显著" if is_cn else "FDR-sig"
        label_effect = "效应量排名输出" if is_cn else "Effect-ranked"
        label_groups = "组数" if is_cn else "Groups"
        label_samples = "样本总数" if is_cn else "Total samples"
        label_no_deg = ("无表达矩阵数据" if is_cn else "No expression matrix data")
        label_empty_deg = ("所有基因均未通过预筛选" if is_cn else "No genes passed pre-filter")
        label_warn = "⚠️ 统计说明" if is_cn else "⚠️ Statistical notes"
        label_small_n = ("样本量小 -> FDR 校正严格 -> 关注效应量排名而非 p 值" if is_cn
                         else "Small n -> strict FDR -> prioritize effect size over p-value")
        label_eff_guide = ("效应量参考: η²>0.3 强效应, >0.1 中等效应, >0.05 弱效应" if is_cn
                           else "Effect size guide: η²>0.3 strong, >0.1 medium, >0.05 weak")

        lines.append(f"## {label_deg} & {label_dam}")
        lines.append("")

        if deg_report and deg_report.get("n_genes_tested", 0) > 0:
            lines.append(f"**{label_deg}**:")
            lines.append(f"  {label_tested}: {deg_report['n_genes_tested']}")
            lines.append(f"  {label_fdr}: {deg_report.get('n_fdr_significant', 0)}")
            lines.append(f"  {label_effect}: {deg_report.get('n_effect_ranked', 0)}")
            lines.append(f"  {label_groups}: {deg_report.get('n_groups', '?')}")
            lines.append(f"  {label_samples}: {deg_report.get('total_samples', '?')}")
            top_genes = deg_report.get("top_genes", [])[:10]
            if top_genes:
                lines.append(f"  Top DEG (效应量排名):")
                for g in top_genes[:10]:
                    sig = "★" if g.get("is_fdr_sig") else " "
                    lines.append(
                        f"    {sig} {g['gene_id']}: η²={g['effect_size']:.3f}, "
                        f"log2FC_max={g['max_log2fc']:.2f}, p={g['p_value']:.2e}"
                    )
            deg_warnings = deg_report.get("warnings", [])
            if deg_warnings:
                lines.append(f"  {label_warn}:")
                for w in deg_warnings:
                    lines.append(f"    - {w}")
            lines.append(f"  {label_eff_guide}")
            # DEG volcano figure
            deg_fig = deg_report.get("figure_markdown", "")
            if deg_fig:
                lines.append(f"\n{deg_fig}")
        elif deg_report:
            reason = deg_report.get("reason", label_no_deg)
            lines.append(f"**{label_deg}**: {reason}")
        else:
            lines.append(f"**{label_deg}**: {label_empty_deg}")

        lines.append("")

        if dam_report and dam_report.get("n_metabolites_tested", 0) > 0:
            lines.append(f"**{label_dam}**:")
            lines.append(f"  检验代谢物数: {dam_report['n_metabolites_tested']}")
            lines.append(f"  {label_fdr}: {dam_report.get('n_fdr_significant', 0)}")
            top_metas = dam_report.get("top_metabolites", [])[:10]
            if top_metas:
                lines.append(f"  Top DAM (效应量排名):")
                for m in top_metas[:10]:
                    sig = "★" if m.get("is_fdr_sig") else " "
                    lines.append(
                        f"    {sig} {m['metabolite']}: η²={m['effect_size']:.3f}, "
                        f"log2FC_max={m['max_log2fc']:.2f}"
                    )
            dam_warnings = dam_report.get("warnings", [])
            if dam_warnings:
                for w in dam_warnings:
                    lines.append(f"  - {w}")
            # DAM volcano figure
            dam_fig = dam_report.get("figure_markdown", "")
            if dam_fig:
                lines.append(f"\n{dam_fig}")
        elif dam_report:
            reason = dam_report.get("reason", "无代谢物数据")
            lines.append(f"**{label_dam}**: {reason}")

        # ── Multi-omics integration section ───────────────
        multiomics_report = getattr(state, "multiomics_report", None) or {}
        if multiomics_report and multiomics_report.get("n_correlation_pairs_tested", 0) > 0:
            label_multi = "多组学联合分析 (DEG-DAM 相关性)" if is_cn else "Multi-omics Integration"
            lines.append("")
            lines.append(f"**{label_multi}**:")
            lines.append(f"  检验相关对数: {multiomics_report['n_correlation_pairs_tested']}")
            lines.append(f"  FDR-显著相关对: {multiomics_report.get('n_significant_pairs', 0)}")
            lines.append(f"  共同样本数: {multiomics_report.get('n_common_samples', '?')}")
            top_pairs = multiomics_report.get("top_pairs", [])[:10]
            if top_pairs:
                lines.append(f"  Top 基因-代谢物相关对:")
                for pair in top_pairs[:10]:
                    sig = "★" if pair.get("q_value", 1) < 0.05 else " "
                    lines.append(
                        f"    {sig} {pair['gene_id']} ↔ {pair['metabolite']}: "
                        f"r={pair['correlation']:.3f}, q={pair.get('q_value', 1):.2e}"
                    )
            modules = multiomics_report.get("modules", {})
            if modules:
                lines.append(f"  共表达模块:")
                for mod_name, genes in modules.items():
                    lines.append(f"    {mod_name}: {', '.join(genes[:5])}{'...' if len(genes) > 5 else ''}")
            multi_warnings = multiomics_report.get("warnings", [])
            for w in multi_warnings:
                lines.append(f"  - {w}")

        # ── v4.7 Multi-omics extension reports ──────────────
        # Joint KEGG enrichment
        enrichment = getattr(state, "joint_enrichment_report", None) or {}
        if enrichment.get("n_pathways_tested", 0) > 0:
            label_enrich = "联合 KEGG 富集分析" if is_cn else "Joint KEGG Enrichment"
            lines.append("")
            lines.append(f"**{label_enrich}**:")
            lines.append(f"  检验通路数: {enrichment['n_pathways_tested']}")
            lines.append(f"  FDR-显著通路: {enrichment.get('n_significant_pathways', 0)}")
            sig_list = enrichment.get("significant", enrichment.get("results", []))[:5]
            for r in sig_list:
                if r.get("q_value", 1) < 0.05:
                    lines.append(f"  ★ {r.get('pathway_name', '?')}: q={r['q_value']:.2e}, FE={r.get('fold_enrichment', 0):.1f}x")
            fig_md = enrichment.get("figure_markdown", "")
            if fig_md:
                lines.append(f"\n{fig_md}\n")

        # Quadrant plot
        quadrant = getattr(state, "quadrant_plot_report", None) or {}
        if quadrant.get("n_pairs", 0) > 0:
            label_quad = "四象限图 (DEGxDAM)" if is_cn else "Quadrant Plot"
            lines.append("")
            lines.append(f"**{label_quad}**:")
            lines.append(f"  分类对数: {quadrant['n_pairs']}")
            lines.append(f"  协同上调 (Q1): {quadrant.get('Q1_coordinated_up', 0)}")
            lines.append(f"  协同下调 (Q3): {quadrant.get('Q3_coordinated_down', 0)}")
            lines.append(f"  同步率: {quadrant.get('synchronicity_ratio', 0):.1%}")
            interp = quadrant.get("interpretation", "")
            if interp:
                lines.append(f"  解读: {interp}")
            fig_md = quadrant.get("figure_markdown", "")
            if fig_md:
                lines.append(f"\n{fig_md}\n")

        # Correlation network
        corr_net = getattr(state, "correlation_network_report", None) or {}
        if corr_net.get("n_edges", 0) > 0:
            label_net = "Spearman 相关网络" if is_cn else "Spearman Correlation Network"
            lines.append("")
            lines.append(f"**{label_net}**:")
            lines.append(f"  边数: {corr_net['n_edges']}")
            lines.append(f"  节点数: {corr_net.get('n_nodes', 0)}")
            hubs = corr_net.get("hubs", [])
            if hubs:
                lines.append(f"  Hub 节点: {', '.join(h['node'] for h in hubs[:6])}")
            modules = corr_net.get("modules", [])
            if modules:
                lines.append(f"  网络模块: {len(modules)} 个")
            fig_md = corr_net.get("figure_markdown", "")
            if fig_md:
                lines.append(f"\n{fig_md}\n")

        # WGCNA
        wgcna_rpt = getattr(state, "wgcna_report", None) or {}
        if wgcna_rpt.get("n_modules", 0) > 0:
            label_wgcna = "WGCNA 共表达网络" if is_cn else "WGCNA"
            lines.append("")
            lines.append(f"**{label_wgcna}**:")
            lines.append(f"  分析基因数: {wgcna_rpt.get('n_genes_analyzed', 0)}")
            lines.append(f"  模块数: {wgcna_rpt['n_modules']}")
            lines.append(f"  软阈值 power: {wgcna_rpt.get('soft_power', '?')}")
            lines.append(f"  Scale-free R²: {wgcna_rpt.get('scale_free_r2', 0):.3f}")
            trait_corrs = wgcna_rpt.get("trait_correlations", [])
            sig_traits = [t for t in trait_corrs if t.get("q_value", 1) < 0.05]
            if sig_traits:
                lines.append(f"  显著模块-性状关联: {len(sig_traits)}")
                for t in sig_traits[:5]:
                    lines.append(f"    M{t['module_id']} ↔ {t['trait']}: r={t['correlation']:.3f}, q={t['q_value']:.2e}")
            fig_md = wgcna_rpt.get("figure_markdown", "")
            if fig_md:
                lines.append(f"\n{fig_md}\n")

        # O2PLS
        o2pls_rpt = getattr(state, "o2pls_report", None) or {}
        if o2pls_rpt.get("n_joint_components", 0) > 0:
            label_o2pls = "O2PLS 多变量整合" if is_cn else "O2PLS Integration"
            lines.append("")
            lines.append(f"**{label_o2pls}**:")
            lines.append(f"  联合成分数: {o2pls_rpt['n_joint_components']}")
            lines.append(f"  X 方差解释: {o2pls_rpt.get('x_joint_variance_explained', 0):.1%}")
            lines.append(f"  Y 方差解释: {o2pls_rpt.get('y_joint_variance_explained', 0):.1%}")
            top_genes = o2pls_rpt.get("top_x_variables", [])[:5]
            if top_genes:
                lines.append(f"  Top contributing genes (VIP):")
                for g in top_genes:
                    lines.append(f"    {g['name']}: loading={g.get('loading', 0):.3f}, VIP={g.get('vip', 0):.2f}")
            fig_md = o2pls_rpt.get("figure_markdown", "")
            if fig_md:
                lines.append(f"\n{fig_md}\n")

        lines.append("")
        lines.append(f"  {label_small_n}")
        lines.append("")

        if not deg_report and not dam_report:
            return L["no_results"]
        return "\n".join(lines)

    lines = []
    lines.append(f"{L['pathway']}: {pathway}")
    lines.append(f"{L['primary_mechanism']}: {mechanism_label}")
    lines.append(f"{L['samples']}: {n_samples}")
    lines.append(f"{L['hypotheses']}: {len(hypotheses)}")
    lines.append("")

    for h in hypotheses:
        lines.append(f"### {h.id}: {h.title}")
        lines.append(f"  {L['type']}: {h.mechanism_type.value}")
        lines.append(f"  {L['uncertainty']}: {h.uncertainty_level}")
        if h.proposed_regulators:
            lines.append(f"  {L['candidate_regulators']}: {', '.join(h.proposed_regulators[:5])}")

        if h.supporting_evidence:
            lines.append(f"  {L['supporting']} ({len(h.supporting_evidence)}):")
            for ev in h.supporting_evidence[:3]:
                lines.append(f"    + [{ev.source}] {ev.description[:100]} ({L['score']}={ev.score:.2f})")

        if h.contradictory_evidence:
            lines.append(f"  {L['contradictory']} ({len(h.contradictory_evidence)}):")
            for ev in h.contradictory_evidence[:2]:
                sev = f"[{ev.contradiction_severity}]" if ev.contradiction_severity else ""
                lines.append(f"    - {sev} {ev.description[:100]}")

        if h.missing_evidence:
            lines.append(f"  {L['missing']}:")
            for m in h.missing_evidence[:3]:
                lines.append(f"    ? {m}")

        if h.why_wrong:
            lines.append(f"  {L['why_wrong']}: {h.why_wrong[:200]}")

        if h.falsification_tests:
            lines.append(f"  {L['falsification']} ({len(h.falsification_tests)}):")
            for ft in h.falsification_tests:
                lines.append(f"    {L['test']}: {ft.test_name}")
                lines.append(f"    {L['if_correct']}: {ft.prediction[:80]}")
                if ft.contradictory_prediction:
                    lines.append(f"    {L['if_wrong']}: {ft.contradictory_prediction[:80]}")

        lines.append("")

    if hasattr(state, "hypothesis") and state.hypothesis:
        lines.append("---")
        lines.append(f"### {L['full_report']}")
        lines.append(state.hypothesis[:1500])

    progress = getattr(state, "progress_messages", []) or []
    if progress:
        lines.append("---")
        lines.append(f"### {L['exec_progress']}")
        for msg in progress:
            lines.append(f"  {msg}")
        lines.append("")

    degraded = getattr(state, "degraded_tools", []) or []
    deg_notes = getattr(state, "degradation_notes", []) or []
    if degraded:
        lines.append("---")
        lines.append(f"### {L['degradation']}")
        lines.append(f"{L['degraded_tools_list']}: {', '.join(degraded)}")
        for note in deg_notes:
            lines.append(f"  - {note}")
        lines.append(L["degraded_note"])

    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════
# Query expansion helper -- resolves gene IDs to synonyms
# ═══════════════════════════════════════════════════════════════

_GENE_ID_PATTERNS = [
    _re.compile(r'Os\d{2}g\d{7}', _re.IGNORECASE),      # RAP-DB
    _re.compile(r'AT[1-5MC]G\d{5}', _re.IGNORECASE),     # TAIR
    _re.compile(r'LOC_Os\d{2}g\d{5}', _re.IGNORECASE),   # MSU/RGAP
]


def _expand_gene_ids_in_query(query: str) -> str:
    """Detect gene IDs in a query string and expand with known synonyms.

    Example: "Os01g0884300 regulation" -> "(Os01g0884300 OR OsNAC6 OR SNAC1) regulation"
    """
    try:
        from phyto_reason.knowledge.gene_resolver import expand_gene_query

        words = query.split()
        expanded_parts = []
        modified = False

        for word in words:
            clean = word.strip('(),;"\'')
            if any(p.search(clean) for p in _GENE_ID_PATTERNS):
                try:
                    expanded = expand_gene_query(clean)
                    if expanded != clean:
                        expanded_parts.append(expanded)
                        modified = True
                        continue
                except Exception:
                    pass
            expanded_parts.append(word)

        return " ".join(expanded_parts) if modified else query

    except Exception:
        return query  # non-critical; return original query on failure


# ═══════════════════════════════════════════════════════════════
# v5.0 Fine-grained pipeline node handlers (for ReAct agent)
# ═══════════════════════════════════════════════════════════════


def _make_runtime_state(session) -> "RuntimeState":
    """Create a minimal RuntimeState from session data for standalone node execution."""
    from phyto_reason.workflows.runtime_state import RuntimeState
    from phyto_reason.models.planner_state import PlannerState

    expr = getattr(session, "expression_matrix", None)
    meta = getattr(session, "metabolite_matrix", None)

    ps = PlannerState(
        species=getattr(session, "species", ""),
        target_metabolite=getattr(session, "target_metabolite", ""),
        expression_matrix=expr,
        metabolite_matrix=meta,
        sample_metadata=getattr(session, "sample_metadata", None),
        promoter_sequences=getattr(session, "promoter_sequences", None),
        sample_count=len(next(iter(expr.values()))) if expr else 0,
        has_expression=bool(expr),
        has_metabolite=bool(meta),
    )

    # TF 预测注释（上传 FASTA 生成）注入 TF 分析链
    annotation_index = None
    ann_csv = getattr(session, "tf_annotation_csv", "")
    if ann_csv and Path(ann_csv).exists():
        try:
            from phyto_reason.ingestion.mappers.annotation_mapper import load_annotation
            annotation_index = load_annotation(ann_csv)
        except Exception as exc:
            logging.getLogger("tool_handlers").warning(
                "TF annotation csv load failed (%s): %s", ann_csv, exc)

    state = RuntimeState(
        session_id=getattr(session, "session_id", ""),
        planner_state=ps,
        annotation_index=annotation_index,
        current_node="data_qc",
        sample_alignment_report=getattr(session, "sample_alignment_report", {}) or {},
    )
    # Pre-populate cached results from earlier tool calls
    for field_name in (
        "qc_report", "deg_report", "dam_report",
        "multiomics_report", "wgcna_report",
    ):
        val = getattr(session, field_name, None)
        if val:
            setattr(state, field_name, val)
    return state


def handle_check_data_quality(arguments: dict, session=None) -> str:
    """会话级数据质控：复用 data_qc_node，输出统计摘要 + 质控图表。"""
    if not session or not (
        getattr(session, "expression_matrix", None)
        or getattr(session, "metabolite_matrix", None)
    ):
        return "No data available. Upload data first, then re-run."

    from phyto_reason.workflows.execution_nodes import data_qc_node

    state = _make_runtime_state(session)
    state.current_node = "data_qc"
    try:
        data_qc_node(state)
    except Exception as exc:
        logger.warning("QC node failed, falling back to text QC: %s", exc)
        return handle_check_quality(
            arguments,
            data={
                "expression": getattr(session, "expression_matrix", None),
                "metabolite": getattr(session, "metabolite_matrix", None),
            },
        )

    report = getattr(state, "qc_report", None) or {}
    lines = ["## Data Quality Control"]
    lines.append(
        f"Genes: {len(getattr(session, 'expression_matrix', None) or {})} | "
        f"Metabolites: {len(getattr(session, 'metabolite_matrix', None) or {})} | "
        f"QC: {'PASS' if getattr(state, 'data_quality_passed', True) else 'ISSUES FOUND'}"
    )
    issues = list(getattr(state, "data_quality_issues", []) or [])
    if issues:
        lines.append("Issues: " + "; ".join(str(i)[:120] for i in issues[:5]))
    for key in ("figure_correlation", "figure_pca", "figure_heatmap"):
        url = report.get(key)
        if url:
            lines.append(f"![{key}]({url})")
    cached = getattr(session, "figure_captions", None)
    if cached is not None:
        try:
            from phyto_reason.reports.analysis_report import FigureNumberer
            numberer = FigureNumberer(cached)
            for key in ("figure_correlation", "figure_pca", "figure_heatmap"):
                url = report.get(key)
                if url:
                    numberer.assign(url)
        except Exception:
            pass
    return "\n".join(lines)


def handle_run_deg_analysis(arguments: dict, session=None) -> str:
    """Differential expression analysis as a standalone callable tool."""
    if not session or not getattr(session, "expression_matrix", None):
        return "No expression data available. Upload data first, then re-run."

    from phyto_reason.workflows.execution_nodes import deg_analysis_node

    state = _make_runtime_state(session)
    state.current_node = "deg_analysis"

    try:
        result = deg_analysis_node(state)
        report = state.deg_report or result.get("deg_report", {})

        # Cache for downstream tools
        if report:
            try:
                session.deg_report = report
            except Exception:
                pass

        if not report:
            return "DEG analysis completed but produced no results."

        lines = ["## DEG Analysis Results"]
        lines.append(
            f"Genes tested: {report.get('n_genes_tested', 0)} | "
            f"FDR-significant: {report.get('n_fdr_significant', 0)} | "
            f"Groups: {report.get('n_groups', 0)}"
        )
        top = report.get("top_genes", [])[:10]
        if top:
            lines.append("Top differentially expressed genes (by effect size):")
            for g in top:
                lines.append(
                    f"  {g.get('gene_id', '?')}: "
                    f"eta_sq={g.get('effect_size', 0):.3f}, "
                    f"max_log2FC={g.get('max_log2fc', 0):.2f}, "
                    f"q={g.get('q_value', 0):.2e}"
                )
        for w in report.get("warnings", []):
            lines.append(f"Warning: {w}")
        figure_md = report.get("figure_markdown", "")
        if figure_md:
            lines.append(f"\n{figure_md}")
        return "\n".join(lines)

    except Exception as e:
        logger.exception("DEG analysis handler failed")
        return f"DEG analysis error: {e}"


def handle_run_dam_analysis(arguments: dict, session=None) -> str:
    """Differential metabolite analysis as a standalone callable tool."""
    if not session or not getattr(session, "metabolite_matrix", None):
        return "No metabolite data available. Upload metabolite data first."

    from phyto_reason.workflows.execution_nodes import dam_analysis_node

    state = _make_runtime_state(session)
    state.current_node = "dam_analysis"

    try:
        result = dam_analysis_node(state)
        report = state.dam_report or result.get("dam_report", {})

        if report:
            try:
                session.dam_report = report
            except Exception:
                pass

        if not report:
            return "DAM analysis completed but produced no results."

        lines = ["## DAM Analysis Results"]
        lines.append(
            f"Metabolites tested: {report.get('n_metabolites_tested', 0)} | "
            f"FDR-significant: {report.get('n_fdr_significant', 0)} | "
            f"n={report.get('n_total', '?')} | "
            f"feasibility={report.get('feasibility', 'unknown')} | "
            f"groups={report.get('group_sizes', {})} | "
            f"group_source={report.get('group_source', 'unknown')}"
        )
        if report.get("n_pairwise_tests", 0):
            lines.append(
                f"Pairwise tissue comparisons: {report.get('n_pairwise_tests', 0)} tests | "
                f"BH-significant: {report.get('n_pairwise_significant', 0)} | "
                f"method={report.get('pairwise_method', 'moderated t + BH')}"
            )
        top = report.get("top_metabolites", [])[:10]
        if top:
            lines.append("Top differentially accumulated metabolites:")
            for m in top:
                lines.append(
                    f"  {m.get('metabolite', m.get('metabolite_id', '?'))}: "
                    f"effect_size={m.get('effect_size', 0):.3f}, "
                    f"q={m.get('q_value', 0):.2e}"
                )
        for warning in report.get("warnings", []):
            lines.append(f"Warning: {warning}")
        figure_md = report.get("figure_markdown", "")
        if figure_md:
            lines.append(f"\n{figure_md}")
        pairwise_figure_md = report.get("pairwise_heatmap_markdown", "")
        if pairwise_figure_md:
            lines.append(f"\n{pairwise_figure_md}")
        return "\n".join(lines)

    except Exception as e:
        logger.exception("DAM analysis handler failed")
        return f"DAM analysis error: {e}"


def handle_run_multiomics(arguments: dict, session=None) -> str:
    """Multi-omics integration as a standalone callable tool."""
    if not session or not getattr(session, "expression_matrix", None):
        return "No expression data available for multi-omics integration."

    from phyto_reason.workflows.execution_nodes import multiomics_integration_node

    state = _make_runtime_state(session)
    state.current_node = "multiomics_integration"

    # Check prerequisites — need DEG or expression data with enough samples
    if not state.deg_report and not state.dam_report:
        return (
            "Multi-omics integration needs DEG and/or DAM results. "
            "Please run run_deg_analysis and/or run_dam_analysis first, "
            "then re-run multi-omics integration."
        )

    try:
        result = multiomics_integration_node(state)
        report = state.multiomics_report or result.get("multiomics_report", {})

        if report:
            try:
                session.multiomics_report = report
            except Exception:
                pass

        if not report:
            return "Multi-omics integration completed but produced no results."

        lines = ["## Multi-Omics Integration Results"]
        lines.append(f"Correlation pairs tested: {report.get('n_correlation_pairs_tested', 0)}")
        lines.append(f"Significant pairs (FDR<0.05): {report.get('n_significant_pairs', 0)}")
        top = report.get("top_pairs", [])[:8]
        if top:
            lines.append("Top gene-metabolite correlations:")
            for p in top:
                lines.append(
                    f"  {p.get('gene', '?')} <-> {p.get('metabolite', '?')}: "
                    f"r={p.get('correlation', 0):.3f}"
                )
        modules = report.get("modules", [])
        if modules:
            lines.append(f"Correlation modules: {len(modules)}")
        return "\n".join(lines)

    except Exception as e:
        logger.exception("Multi-omics handler failed")
        return f"Multi-omics integration error: {e}"


def handle_run_wgcna(arguments: dict, session=None) -> str:
    """WGCNA as a standalone callable tool."""
    if not session or not getattr(session, "expression_matrix", None):
        return "No expression data available for WGCNA."

    from phyto_reason.workflows.execution_nodes import wgcna_node

    state = _make_runtime_state(session)
    state.current_node = "wgcna"

    # WGCNA needs 8+ samples
    expr = session.expression_matrix
    if expr:
        n_samples = len(next(iter(expr.values())))
        if n_samples < 8:
            return (
                f"WGCNA requires 8+ samples for reliable module detection. "
                f"Current data has {n_samples} samples. "
                f"Try DEG analysis or multi-omics integration instead."
            )

    try:
        result = wgcna_node(state)
        report = state.wgcna_report or result.get("wgcna_report", {})

        if report:
            try:
                session.wgcna_report = report
            except Exception:
                pass

        if not report:
            return "WGCNA completed but produced no results."

        lines = ["## WGCNA Results"]
        lines.append(f"Modules detected: {report.get('n_modules', 0)}")
        lines.append(f"Soft power: {report.get('soft_power', '?')}")
        return "\n".join(lines)

    except Exception as e:
        logger.exception("WGCNA handler failed")
        return f"WGCNA error: {e}"


def handle_run_tf_analysis(arguments: dict, session=None) -> str:
    """TF candidate analysis as a standalone callable tool."""
    if not session:
        return "No session available."

    from phyto_reason.workflows.execution_nodes import tf_narrowing_node

    state = _make_runtime_state(session)
    state.current_node = "tf_narrowing"
    state.target_metabolite = (
        arguments.get("target_metabolite", "")
        or getattr(session, "target_metabolite", "")
    )

    # TF analysis needs DEG results
    if not state.deg_report:
        return (
            "TF analysis requires DEG results. "
            "Run DEG analysis (run_deg_analysis) first to identify "
            "differentially expressed genes, then re-run TF analysis."
        )

    try:
        result = tf_narrowing_node(state)

        if state.tf_candidates:
            try:
                session.tf_candidates = state.tf_candidates
            except Exception:
                pass

        tf_candidates = state.tf_candidates or {}
        if not tf_candidates:
            return "TF analysis completed but no candidates were identified."

        lines = ["## TF Candidate Analysis"]
        lines.append(f"Candidates identified: {len(tf_candidates)}")
        for name, candidate in list(tf_candidates.items())[:10]:
            if hasattr(candidate, "calibrated_score"):
                lines.append(
                    f"  {name}: score={candidate.calibrated_score:.3f}"
                )
            elif isinstance(candidate, dict):
                lines.append(
                    f"  {name}: score={candidate.get('calibrated_score', '?')}"
                )
            else:
                lines.append(f"  {name}")
        return "\n".join(lines)

    except Exception as e:
        logger.exception("TF analysis handler failed")
        return f"TF analysis error: {e}"


def handle_run_hypothesis_synthesis(arguments: dict, session=None) -> str:
    """Hypothesis synthesis as a standalone callable tool."""
    if not session:
        return "No session available."

    from phyto_reason.workflows.execution_nodes import hypothesis_synthesis_node

    state = _make_runtime_state(session)
    state.current_node = "hypothesis_synthesis"

    try:
        result = hypothesis_synthesis_node(state)
        hypothesis_text = state.mechanism_hypothesis or ""

        if not hypothesis_text:
            return (
                "Hypothesis synthesis completed but no hypotheses were generated. "
                "Run DEG analysis and TF analysis first to accumulate evidence, "
                "then re-run hypothesis synthesis."
            )

        return hypothesis_text

    except Exception as e:
        logger.exception("Hypothesis synthesis handler failed")
        return f"Hypothesis synthesis error: {e}"
