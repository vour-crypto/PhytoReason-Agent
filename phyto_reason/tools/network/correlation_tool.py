"""
correlation_tool.py — TF-代谢物相关性分析工具。

使用 CorrelationEngine 执行 bicor / spearman / pearson，
所有筛选经过 FDR correction。
"""

from __future__ import annotations

import logging
from datetime import datetime

import numpy as np

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.tools.network.correlation_engine import CorrelationEngine
from phyto_reason.tools.network.confounder_engine import ConfounderEngine

logger = logging.getLogger("correlation_tool")


@register_tool
class CorrelationTool(BaseTool):
    tool_name = "correlation"
    description = "计算 TF 基因与代谢物之间的相关性（bicor 优先），经 FDR 校正后筛选显著对。"
    version = "2.0.0"
    parameters = [
        ToolParameter(name="tf_ids", type="array", description="TF 基因 ID 列表", required=True),
        ToolParameter(name="metabolite_names", type="array", description="代谢物名称列表", required=True),
        ToolParameter(name="expression_matrix", type="object", description="表达矩阵 {gene: {sample: val}}", required=True),
        ToolParameter(name="metabolite_matrix", type="object", description="代谢物矩阵 {metabolite: {sample: val}}", required=True),
        ToolParameter(name="method", type="string", description="相关方法: bicor|pearson|spearman", default="bicor"),
        ToolParameter(name="corr_threshold", type="number", description="最小 |r| 阈值；默认 0.8，低于 0.8 视为 exploratory", default=0.8),
        ToolParameter(name="fdr_threshold", type="number", description="FDR 阈值", default=0.05),
        ToolParameter(name="tissue_metadata", type="object", description="组织信息 {tissue: [sample, ...]}", default=None),
    ]
    supported_data_types = ["expression", "metabolite"]

    def __init__(self) -> None:
        self.engine = CorrelationEngine()
        self.confounder = ConfounderEngine()

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        for key in ["tf_ids", "metabolite_names", "expression_matrix", "metabolite_matrix"]:
            if not kwargs.get(key):
                errors.append(f"{key} 不能为空")
        return errors

    def run(self, **kwargs) -> ToolResult:
        tf_ids: list[str] = kwargs.get("tf_ids", [])
        metabolite_names: list[str] = kwargs.get("metabolite_names", [])
        expr_mat: dict = kwargs.get("expression_matrix", {})
        meta_mat: dict = kwargs.get("metabolite_matrix", {})
        method: str = kwargs.get("method", "bicor")
        corr_threshold: float = kwargs.get("corr_threshold", 0.8)
        fdr_threshold: float = kwargs.get("fdr_threshold", 0.05)
        tissue_meta: dict | None = kwargs.get("tissue_metadata")

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []
        exploratory = corr_threshold < 0.8
        if exploratory:
            warnings.append(
                f"探索性相关：显式阈值 |r|≥{corr_threshold:.2f}（正式默认 |r|≥0.80）。"
            )

        self.engine = CorrelationEngine(
            method=method,
            corr_threshold=corr_threshold,
            fdr_threshold=fdr_threshold,
        )

        tf_available = [t for t in tf_ids if t in expr_mat]
        meta_available = [m for m in metabolite_names if m in meta_mat]
        if not tf_available or not meta_available:
            return ToolResult(candidates=[], evidence_list=[],
                              warnings=["TF 列表或代谢物列表在矩阵中未找到"],
                              metadata={"n_tf_tested": 0, "n_meta_tested": 0, "n_significant": 0})

        sample_set: set[str] | None = None
        for s in expr_mat.values():
            sample_set = set(s.keys())
            break
        for s in meta_mat.values():
            if sample_set is not None:
                sample_set &= set(s.keys())
        if sample_set is None or len(sample_set) < 3:
            return ToolResult(candidates=[], evidence_list=[], warnings=["共同样本不足"],
                              metadata={"n_significant": 0})

        common_samples = sorted(sample_set)
        n = len(common_samples)

        cov_matrix = None
        if tissue_meta:
            cov_matrix = self.confounder.build_covariate_matrix(tissue_meta, common_samples)

        all_pairs: list[tuple[str, str, float, float]] = []
        for tf_id in tf_available:
            tf_vec = np.array([expr_mat[tf_id].get(s, np.nan) for s in common_samples], dtype=float)
            if np.isnan(tf_vec).any():
                continue
            for meta_name in meta_available:
                meta_vec = np.array([meta_mat[meta_name].get(s, np.nan) for s in common_samples], dtype=float)
                if np.isnan(meta_vec).any():
                    continue
                if cov_matrix is not None and cov_matrix.shape[1] > 0:
                    r, p = self.confounder.partial_corr(tf_vec, meta_vec, cov_matrix)
                else:
                    r, p = self.engine.calculate(tf_vec, meta_vec)
                if abs(r) >= corr_threshold:
                    all_pairs.append((tf_id, meta_name, r, p))

        filtered = self.engine.filter_pairs(all_pairs)
        significant_pairs = len(filtered)

        seen_tfs: dict[str, float] = {}
        for tf_id, meta_name, r, p, q in filtered:
            if tf_id not in seen_tfs or abs(r) > seen_tfs[tf_id]:
                seen_tfs[tf_id] = abs(r)

            candidate = CandidateGene(
                gene_id=tf_id, is_tf=True,
                correlation_score=round(abs(r), 4),
            )
            evidence = Evidence(
                evidence_type=EvidenceType.CORRELATION,
                source=f"correlation_tool_{self.version}",
                score=round(abs(r), 4),
                confidence=round(1.0 - min(q * 10, 0.95), 4),
                description=f"{method} r={r:.4f}, p={p:.2e}, q={q:.2e} (TF={tf_id}, metabolite={meta_name}, n={n})",
                metadata={
                    "method": method, "tf_id": tf_id,
                    "metabolite": meta_name, "r": round(r, 4),
                    "p_value": round(p, 6), "q_value": round(q, 6),
                    "n_samples": n, "fdr_passed": True,
                },
            )
            candidates.append(candidate)
            evidence_list.append(evidence)

        candidates.sort(key=lambda g: g.correlation_score, reverse=True)
        return ToolResult(
            candidates=candidates[:100],
            evidence_list=evidence_list[:100],
            warnings=warnings,
            metadata={
                "method": method, "n_tf_tested": len(tf_available),
                "n_meta_tested": len(meta_available),
                "n_samples": n,
                "n_significant": significant_pairs,
                "corr_threshold": corr_threshold,
                "fdr_threshold": fdr_threshold,
                "exploratory": exploratory,
                "fdr_status": "BH applied to tested pairs",
                "summary": (
                    f"n={n}; method={method}; |r|>={corr_threshold:.2f}; "
                    f"FDR=Benjamini-Hochberg q<{fdr_threshold:.2f}"
                    + ("; exploratory" if exploratory else "")
                ),
                "confounder_corrected": cov_matrix is not None and cov_matrix.shape[1] > 0,
            },
        )
