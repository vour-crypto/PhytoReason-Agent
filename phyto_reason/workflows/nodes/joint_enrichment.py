"""
joint_enrichment.py — 联合 KEGG 通路富集分析节点。

对 DEG + DAM 结果进行联合 KEGG 通路富集分析，使用超几何检验
(Hypergeometric test) 判断差异基因/代谢物是否在特定通路中显著富集。

对标真实植物多组学文章的标准富集分析流程。
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("workflow_nodes.joint_enrichment")


# ── 超几何检验 ──────────────────────────────────────────────

def _hypergeometric_test(k: int, M: int, n: int, N: int) -> float:
    """单尾超几何检验 (enrichment test)。

    H0: 基因在通路中随机分布
    H1: 基因在通路中富集（观察到的重叠大于期望）

    Args:
        k: 观察到的重叠数 (DEG ∩ pathway genes)
        M: 通路基因数 (pathway size)
        n: 差异基因数 (DEG count)
        N: 背景基因总数 (total genes in genome/expression matrix)

    Returns:
        p-value from hypergeometric survival function
    """
    if k <= 0 or M <= 0 or n <= 0 or N <= 0:
        return 1.0
    # Use hypergeometric SF: P(X >= k) = 1 - CDF(k-1, N, n, M)
    from scipy.stats import hypergeom
    # scipy hypergeom: M = total, n = type I, N = draws
    # P(X >= k) = hypergeom.sf(k - 1, M=N, n=M, N=n)
    # Actually: M=total balls(N), n=white balls(M), N=draws(n)
    p = hypergeom.sf(k - 1, M=N, n=M, N=n)
    return float(p)


def _fisher_combined_pvalue(pvalues: list[float]) -> float:
    """Fisher's method for combining independent p-values.

    χ² = -2 Σ ln(p_i), df = 2k
    """
    from scipy.stats import chi2
    if not pvalues:
        return 1.0
    valid = [p for p in pvalues if p > 0 and p <= 1.0]
    if not valid:
        return 1.0
    chi_sq = -2.0 * sum(np.log(p) for p in valid)
    df = 2 * len(valid)
    return float(chi2.sf(chi_sq, df))


# ── 已知植物 KEGG 通路 (fallback gene sets) ──────────────────

# 当 KEGG API 不可用时使用内置通路基因集
# 基于 KEGG Orthology → 拟南芥基因映射
PLANT_KEGG_PATHWAYS: dict[str, dict] = {
    "flavonoid_biosynthesis": {
        "map_id": "map00941",
        "name": "Flavonoid biosynthesis",
        "genes": ["CHS", "CHI", "F3H", "F3'H", "F3'5'H", "FLS", "DFR", "ANS", "LAR", "ANR",
                   "UFGT", "OMT", "C4H", "4CL", "PAL", "CYP75B1", "CYP75B2"],
    },
    "phenylpropanoid_biosynthesis": {
        "map_id": "map00940",
        "name": "Phenylpropanoid biosynthesis",
        "genes": ["PAL", "C4H", "4CL", "CCR", "CAD", "COMT", "CCoAOMT", "HCT", "C3H",
                   "F5H", "POD", "LAC"],
    },
    "terpenoid_backbone_biosynthesis": {
        "map_id": "map00900",
        "name": "Terpenoid backbone biosynthesis",
        "genes": ["DXS", "DXR", "MCT", "CMK", "MDS", "HDS", "HDR", "IDI", "GPPS",
                   "FPPS", "GGPPS"],
    },
    "carotenoid_biosynthesis": {
        "map_id": "map00906",
        "name": "Carotenoid biosynthesis",
        "genes": ["PSY", "PDS", "ZDS", "CRTISO", "LCYB", "LCYE", "BCH", "ZEP", "VDE",
                   "NSY", "NCED"],
    },
    "alkaloid_biosynthesis": {
        "map_id": "map00960",
        "name": "Tropane, piperidine and pyridine alkaloid biosynthesis",
        "genes": ["TDC", "STR", "SGD", "T6ODM", "NMT", "PMT", "MPO", "TR",
                   "CYP82E", "BBE"],
    },
    "glucosinolate_biosynthesis": {
        "map_id": "map00966",
        "name": "Glucosinolate biosynthesis",
        "genes": ["CYP79F1", "CYP79F2", "CYP79B2", "CYP79B3", "CYP83A1", "CYP83B1",
                   "SUR1", "SUR2", "UGT74B1", "UGT74C1", "SOT16", "SOT17", "SOT18"],
    },
    "anthocyanin_biosynthesis": {
        "map_id": "map00942",
        "name": "Anthocyanin biosynthesis",
        "genes": ["CHS", "CHI", "F3H", "F3'H", "DFR", "ANS", "UFGT", "OMT",
                   "MYB75", "MYB90", "MYB113", "MYB114", "TT8", "TTG1"],
    },
}


def _detect_pathway_for_metabolite(metabolite: str) -> str | None:
    """根据代谢物名称推断相关 KEGG 通路。"""
    m = metabolite.lower()
    if any(kw in m for kw in ["flavonoid", "flavonol", "anthocyanin", "proanthocyanidin"]):
        return "flavonoid_biosynthesis"
    if any(kw in m for kw in ["phenylpropanoid", "lignin", "coumarin", "cinnamate"]):
        return "phenylpropanoid_biosynthesis"
    if any(kw in m for kw in ["terpen", "isoprenoid", "gibberellin", "abscisic"]):
        return "terpenoid_backbone_biosynthesis"
    if any(kw in m for kw in ["carotenoid", "carotene", "lycopene", "xanthophyll"]):
        return "carotenoid_biosynthesis"
    if any(kw in m for kw in ["alkaloid", "berberine", "nicotine", "tropane"]):
        return "alkaloid_biosynthesis"
    if any(kw in m for kw in ["nitidine", "chelerythrine", "sanguinarine", "benzylisoquinoline"]):
        return "benzylisoquinoline_alkaloid_biosynthesis"
    if any(kw in m for kw in ["glucosinolate", "sulforaphane"]):
        return "glucosinolate_biosynthesis"
    return None


def _get_pathway_data(pathway_key: str) -> dict | None:
    """获取通路数据（优先从 KEGG API，fallback 到内置知识库）。"""
    # Try KEGG API first
    try:
        from phyto_reason.tools.pathway.kegg_client import KEGGClient
        client = KEGGClient()
        pw_info = PLANT_KEGG_PATHWAYS.get(pathway_key, {})
        map_id = pw_info.get("map_id", "")
        if map_id:
            genes = client.get_pathway_genes(map_id, organism="ath")
            if genes:
                gene_names = [g.get("name", g.get("gene_id", "")) for g in genes]
                return {
                    "name": pw_info.get("name", pathway_key),
                    "map_id": map_id,
                    "genes": gene_names,
                }
    except Exception as e:
        logger.debug("KEGG API unavailable for enrichment, using built-in: %s", e)

    # Fallback to built-in
    return PLANT_KEGG_PATHWAYS.get(pathway_key)


def _match_genes(gene_list: list[str], expression_genes: set[str]) -> set[str]:
    """将通路基因名匹配到表达矩阵中的实际基因 ID。

    使用模糊匹配：如果通路基因名是表达基因 ID 的子串，即视为匹配。
    """
    matched: set[str] = set()
    for pw_gene in gene_list:
        pw_upper = pw_gene.upper()
        for expr_gene in expression_genes:
            if pw_upper in expr_gene.upper() or expr_gene.upper() in pw_upper:
                matched.add(expr_gene)
                break
    return matched


def joint_enrichment_node(state: RuntimeState) -> dict:
    """联合 KEGG 通路富集分析节点。

    执行步骤:
      1. 从 DEG 结果获取显著差异基因列表
      2. 从 DAM 结果获取显著差异代谢物列表
      3. 识别与目标代谢物相关的 KEGG 通路
      4. 对每个通路执行超几何检验：
         - DEG 富集：差异基因是否在该通路中富集
         - 联合富集：Fisher's method 合并 DEG + DAM p-values
      5. FDR 校正多重检验
      6. 返回富集表 + 解释

    Returns:
        dict with joint_enrichment_report
    """
    try:
        # ── WorkflowPlan skip guard ────────────────────────
        from phyto_reason.workflows.execution_nodes import _should_skip_node
        if _should_skip_node(state, "joint_enrichment"):
            state.add_trace("joint_enrichment", "skipped",
                            summary={"reason": "user skipped this step"})
            return {}

        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}
        expr = None
        if state.planner_state:
            expr = state.planner_state.expression_matrix
        target = (state.planner_state.target_metabolite or "").lower().strip() if state.planner_state else ""

        if not deg_report and not dam_report:
            state.add_trace("joint_enrichment", "completed",
                            summary={"reason": "no DEG or DAM data available"})
            return {"joint_enrichment_report": {"reason": "no data"}}

        # Extract significant genes and metabolites
        sig_genes: set[str] = set()
        if deg_report:
            for g in deg_report.get("top_genes", []):
                if g.get("is_fdr_sig") or g.get("effect_size", 0) >= 0.3:
                    sig_genes.add(g["gene_id"])
            # If no FDR-sig, take top effect-size genes
            if not sig_genes and deg_report.get("n_genes_tested", 0) > 0:
                sig_genes = {g["gene_id"] for g in deg_report.get("top_genes", [])[:20]}

        sig_metas: set[str] = set()
        if dam_report:
            for m in dam_report.get("top_metabolites", []):
                if m.get("is_fdr_sig") or m.get("effect_size", 0) >= 0.3:
                    sig_metas.add(m["metabolite"])
            if not sig_metas and dam_report.get("n_metabolites_tested", 0) > 0:
                sig_metas = {m["metabolite"] for m in dam_report.get("top_metabolites", [])[:15]}

        # Background gene set: all genes in expression matrix
        background_genes: set[str] = set(expr.keys()) if expr else set()
        if not background_genes and deg_report.get("top_genes"):
            # Fallback: use all genes that had DEG testing
            background_genes = sig_genes.copy()
        N = max(len(background_genes), len(sig_genes), 1)

        if not sig_genes:
            state.add_trace("joint_enrichment", "completed",
                            summary={"reason": "no significant DEGs to test"})
            return {"joint_enrichment_report": {"reason": "no significant DEGs", "n_pathways_tested": 0}}

        # Identify pathways to test
        pathway_keys: list[str] = []
        if target:
            detected = _detect_pathway_for_metabolite(target)
            if detected:
                pathway_keys.append(detected)
            # Also add related pathways
            for key in PLANT_KEGG_PATHWAYS:
                if key not in pathway_keys:
                    pathway_keys.append(key)
        else:
            pathway_keys = list(PLANT_KEGG_PATHWAYS.keys())

        # Species profiles may provide a pathway prior for non-model plants.
        # Resolve its enzyme anchors through the uploaded annotation index so
        # the fallback remains evidence-backed rather than inventing matches.
        profile = None
        try:
            from phyto_reason.knowledge.species_registry import get_species_registry
            profile = get_species_registry().get(
                state.planner_state.species if state.planner_state else ""
            )
            if profile:
                for profile_key in profile.pathway_prior:
                    if profile_key not in pathway_keys:
                        pathway_keys.append(profile_key)
        except Exception as exc:
            logger.debug("Species pathway prior unavailable: %s", exc)

        # Limit to reasonable number
        pathway_keys = pathway_keys[:10]

        # Run hypergeometric test for each pathway
        enrichment_results: list[dict] = []
        all_pvalues: list[float] = []
        all_deg_pvalues: list[float] = []

        for pw_key in pathway_keys:
            pw_data = _get_pathway_data(pw_key)
            if profile and pw_key in profile.pathway_prior:
                info = profile.pathway_prior.get(pw_key, {}) or {}
                anchor_genes: set[str] = set()
                ann = getattr(state, "annotation_index", None)
                for enzyme in info.get("enzymes", []) or []:
                    if ann is not None:
                        terms = [str(enzyme)]
                        terms.extend(str(enzyme).split("(", 1)[:1])
                        terms.extend(info.get("search_terms", []) or [])
                        for term in terms:
                            for gene_id in ann.search_substring(term):
                                base_id = str(gene_id).split(".", 1)[0]
                                for background_id in background_genes:
                                    if str(background_id).split(".", 1)[0] == base_id:
                                        anchor_genes.add(background_id)
                if anchor_genes:
                    pw_data = {
                        "name": pw_key,
                        "map_id": f"profile:{pw_key}",
                        "genes": sorted(anchor_genes),
                    }
            if not pw_data:
                continue

            pw_genes = pw_data.get("genes", [])
            if not pw_genes:
                continue

            # Match pathway genes to background
            matched_pw_genes = _match_genes(pw_genes, background_genes)
            M_pw = len(matched_pw_genes)  # pathway gene count in background
            k_overlap = len(sig_genes & matched_pw_genes)  # DEG ∩ pathway
            n_deg = len(sig_genes)

            if M_pw == 0:
                continue

            # Hypergeometric test for DEG enrichment
            p_deg = _hypergeometric_test(k=k_overlap, M=M_pw, n=n_deg, N=N)

            # Expected overlap under null
            expected = n_deg * M_pw / N if N > 0 else 0

            # Fold enrichment
            fold_enrichment = k_overlap / expected if expected > 0 else 0.0

            # Combine independent DEG and DAM pathway evidence with Fisher's
            # method.  DAM p-values are assigned to a pathway by the same
            # metabolite-name/pathway mapping used for target detection.
            dam_pvalues: list[float] = []
            if dam_report:
                for dam_item in dam_report.get("top_metabolites", []):
                    metabolite_name = str(dam_item.get("metabolite", dam_item.get("metabolite_id", "")))
                    detected_dam_pathway = _detect_pathway_for_metabolite(metabolite_name.lower())
                    if detected_dam_pathway == pw_key:
                        raw_p = dam_item.get("p_value_raw")
                        if raw_p is None:
                            raw_p = dam_item.get("p_value")
                        if raw_p is None:
                            continue
                        try:
                            parsed_p = float(raw_p)
                        except (TypeError, ValueError):
                            continue
                        if not np.isfinite(parsed_p):
                            continue
                        # Zero can be a legitimate rounded/underflowed p-value;
                        # clamp it for log() rather than treating it as missing.
                        dam_pvalues.append(min(1.0, max(1e-300, parsed_p)))
            p_dam = _fisher_combined_pvalue(dam_pvalues) if dam_pvalues else 1.0
            p_combined = _fisher_combined_pvalue([p_deg, p_dam])

            result = {
                "pathway_id": pw_data.get("map_id", pw_key),
                "pathway_name": pw_data.get("name", pw_key),
                "pathway_key": pw_key,
                "pathway_genes_in_background": M_pw,
                "deg_overlap": k_overlap,
                "n_deg_tested": n_deg,
                "n_background": N,
                "expected_overlap": round(expected, 1),
                "fold_enrichment": round(fold_enrichment, 2),
                "p_value": round(p_combined, 6),
                "p_combined": round(p_combined, 6),
                "p_combined_raw": p_combined,
                "p_deg": round(p_deg, 6),
                "p_deg_raw": p_deg,
                "p_dam": round(p_dam, 6),
                "p_dam_raw": p_dam,
                "n_dam_overlap": len(dam_pvalues),
                "combined_test": "Fisher: X2=-2*sum(log(p_deg, p_dam)); df=4",
                "overlap_genes": sorted(sig_genes & matched_pw_genes),
            }
            enrichment_results.append(result)
            all_pvalues.append(p_combined)
            all_deg_pvalues.append(p_deg)

        if not enrichment_results:
            state.add_trace("joint_enrichment", "completed",
                            summary={"reason": "no pathways had gene matches"})
            return {"joint_enrichment_report": {"n_pathways_tested": 0, "reason": "no pathway gene matches"}}

        # FDR correction
        from phyto_reason.utils.stats_utils import compute_fdr
        q_combined = compute_fdr(all_pvalues)
        q_deg = compute_fdr(all_deg_pvalues)
        for r, q_all, q_gene in zip(enrichment_results, q_combined, q_deg):
            r["q_combined"] = round(q_all, 6)
            r["q_combined_raw"] = q_all
            r["q_value"] = round(q_all, 6)
            r["q_deg"] = round(q_gene, 6)

        # Sort by p-value
        enrichment_results.sort(key=lambda x: x["p_value"])

        # Identify significant enrichments
        sig_enrichments = [r for r in enrichment_results if r["q_value"] < 0.05]
        top_hits = enrichment_results[:5]

        joint_enrichment_report = {
            "n_pathways_tested": len(enrichment_results),
            "n_significant_pathways": len(sig_enrichments),
            "n_deg_tested": len(sig_genes),
            "n_background_genes": N,
            "target_metabolite": target,
            "method": "hypergeometric_test + Fisher_combined_DEG_DAM",
            "combined_pvalue_formula": "p_combined = chi2.sf(-2 * (log(p_deg) + log(p_dam)), df=4); p_dam=1 when no pathway-mapped DAM",
            "fdr_method": "Benjamini-Hochberg",
            "results": enrichment_results,
            "significant": sig_enrichments,
            "warnings": [],
        }

        if enrichment_results and all(r.get("n_dam_overlap", 0) == 0 for r in enrichment_results):
            joint_enrichment_report["warnings"].append(
                "All tested pathways had n_dam_overlap=0; DAM names did not map to the tested KEGG pathways."
            )
        elif any(
            r.get("n_dam_overlap", 0) > 0
            and r.get("p_deg_raw", 1.0) >= 0.999
            and r.get("p_dam_raw", 1.0) < 0.05
            for r in enrichment_results
        ):
            joint_enrichment_report["warnings"].append(
                "Enrichment is DAM-driven: DEG overlap with the pathway anchors is zero "
                "(p_deg=1.0), while DAM evidence contributes the combined result."
            )

        if len(sig_enrichments) == 0:
            joint_enrichment_report["warnings"].append(
                "FDR 校正后无显著富集通路 (q<0.05)。"
                f"可能是背景基因集过大或差异基因过少 (n_DEG={len(sig_genes)})。"
                "建议：放宽 DEG 筛选阈值或检查通路基因匹配。"
            )

        # Generate enrichment bar chart
        figure_url = None
        try:
            from phyto_reason.visualization.multiomics_figures import plot_enrichment_bars
            target_title = (state.planner_state.target_metabolite or "metabolite").strip() if state.planner_state else "metabolite"
            figure_url = plot_enrichment_bars(
                enrichment_results,
                title=f"Joint KEGG Enrichment — {target_title}",
            )
        except Exception as e:
            logger.debug("Enrichment figure generation skipped: %s", e)
        if figure_url:
            joint_enrichment_report["figure_url"] = figure_url
            joint_enrichment_report["figure_markdown"] = f"![Joint KEGG Enrichment]({figure_url})"

        state.joint_enrichment_report = joint_enrichment_report

        state.add_trace("joint_enrichment", "completed", summary={
            "n_pathways": len(enrichment_results),
            "n_sig": len(sig_enrichments),
            "top_pathway": top_hits[0]["pathway_name"] if top_hits else "",
            "top_pvalue": round(top_hits[0]["q_value"], 6) if top_hits else 1.0,
        })

        logger.info(
            "Joint enrichment: %d pathways, %d significant, top: %s",
            len(enrichment_results), len(sig_enrichments),
            top_hits[0]["pathway_name"] if top_hits else "none",
        )

        return {"joint_enrichment_report": joint_enrichment_report}

    except Exception as e:
        logger.error("joint_enrichment_node failed: %s", e, exc_info=True)
        state.add_trace("joint_enrichment", "failed", error=str(e))
        return {"joint_enrichment_report": {"error": str(e)}}
