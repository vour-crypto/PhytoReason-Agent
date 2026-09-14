"""
wgcna_r.py — 真实 WGCNA 的 Python 实现。

使用 scipy 实现:
  - 双权重中位相关 (bicor) 邻接矩阵
  - TOM (Topological Overlap Matrix)
  - 层次聚类 + 动态模块切割
  - 模块特征基因 (eigengene)
  - 模块隶属度 (kME)

不依赖 rpy2/R，纯 Python/numpy/scipy 实现。
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool

logger = logging.getLogger("wgcna_r")


def _bicor(x: np.ndarray, y: np.ndarray) -> float:
    """双权重中位相关 (biweight midcorrelation)。"""
    n = len(x)
    if n < 4:
        return float(pearsonr(x, y)[0])
    med_x, med_y = np.median(x), np.median(y)
    mad_x = np.median(np.abs(x - med_x))
    mad_y = np.median(np.abs(y - med_y))
    mad_x = max(mad_x, 1e-6)
    mad_y = max(mad_y, 1e-6)
    u_x = (x - med_x) / (9.0 * mad_x)
    u_y = (y - med_y) / (9.0 * mad_y)
    w_x = (1 - u_x ** 2) ** 2 * (np.abs(u_x) < 1)
    w_y = (1 - u_y ** 2) ** 2 * (np.abs(u_y) < 1)
    x_bi = (x - med_x) * w_x
    y_bi = (y - med_y) * w_y
    num = np.sum(x_bi * y_bi)
    den = np.sqrt(np.sum(x_bi ** 2) * np.sum(y_bi ** 2))
    return float(num / den) if den > 1e-10 else 0.0


def _compute_adjacency(matrix: np.ndarray, power: int = 6, method: str = "bicor") -> np.ndarray:
    """计算加权邻接矩阵。"""
    n_genes = matrix.shape[0]
    adj = np.zeros((n_genes, n_genes), dtype=np.float64)

    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if method == "bicor":
                r = _bicor(matrix[i], matrix[j])
            else:
                r, _ = pearsonr(matrix[i], matrix[j])
            r = max(-1.0, min(1.0, r))
            adj[i, j] = abs(r) ** power
            adj[j, i] = adj[i, j]
    return adj


def _compute_tom(adjacency: np.ndarray) -> np.ndarray:
    """计算 TOM (Topological Overlap Matrix)。"""
    n = adjacency.shape[0]
    adj_bin = (adjacency > 0).astype(np.float64)
    k = np.sum(adjacency, axis=1)
    tom = np.zeros_like(adjacency)

    for i in range(n):
        for j in range(i, n):
            if i == j:
                tom[i, j] = 1.0
                continue
            lij = np.sum(adjacency[i] * adjacency[j])
            num = lij + adjacency[i, j]
            den = min(k[i], k[j]) + 1 - adjacency[i, j]
            tom[i, j] = num / den if den > 0 else 0.0
            tom[j, i] = tom[i, j]
    return tom


def _compute_eigengene(module_matrix: np.ndarray) -> np.ndarray:
    """通过 SVD 计算模块特征基因 (第一右奇异向量 = 跨样本表达模式)。"""
    u, s, vt = np.linalg.svd(module_matrix, full_matrices=False)
    return vt[0]


@register_tool
class RealWGCNATool(BaseTool):
    """真实 WGCNA 工具 (Python/scipy 实现)。"""
    tool_name = "wgcna_real"
    description = "真实 WGCNA 共表达网络分析 (TOM + 层次聚类 + bicor)。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="expression_matrix", type="object",
                      description="表达矩阵 {gene_id: {sample: value}}", required=True),
        ToolParameter(name="gene_ids", type="array",
                      description="要分析的基因 ID 列表", required=True),
        ToolParameter(name="soft_power", type="number", description="软阈值 power", default=6),
        ToolParameter(name="min_module_size", type="integer", description="最小模块大小", default=10),
        ToolParameter(name="method", type="string", description="相关性方法: bicor|pearson", default="bicor"),
    ]
    supported_data_types = ["expression"]

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        if not kwargs.get("expression_matrix"):
            errors.append("expression_matrix 不能为空")
        if not kwargs.get("gene_ids"):
            errors.append("gene_ids 不能为空")
        return errors

    def run(self, **kwargs) -> ToolResult:
        expr_mat: dict = kwargs["expression_matrix"]
        gene_ids: list[str] = kwargs["gene_ids"]
        soft_power: int = kwargs.get("soft_power", 6)
        min_module_size: int = kwargs.get("min_module_size", 10)
        method: str = kwargs.get("method", "bicor")

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        available = [g for g in gene_ids if g in expr_mat]
        if len(available) < min_module_size:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=[f"可用基因 {len(available)} < 最小模块 {min_module_size}"],
                metadata={"n_genes": len(available), "method": method},
            )

        sample_keys = sorted(next(iter(expr_mat.values())).keys())
        n_samples = len(sample_keys)
        gene_vecs = []
        valid_genes = []
        for g in available:
            vec = np.array([expr_mat[g].get(s, 0.0) for s in sample_keys], dtype=np.float64)
            std = np.std(vec)
            if std > 1e-10:
                gene_vecs.append((vec - np.mean(vec)) / std)
                valid_genes.append(g)

        n_genes = len(valid_genes)
        if n_genes < min_module_size:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=[f"标准化后可用基因 {n_genes} < {min_module_size}"],
                metadata={"n_genes": n_genes, "method": method},
            )

        matrix = np.array(gene_vecs)
        adj = _compute_adjacency(matrix, power=soft_power, method=method)
        tom = _compute_tom(adj)
        diss_tom = np.clip(1 - tom, 0, 1)

        # 层次聚类
        condensed = squareform(diss_tom, checks=False)
        link = linkage(condensed, method="average")
        clusters = fcluster(link, t=1 - 0.2, criterion="distance")

        # 分配模块
        module_map: dict[int, list[int]] = {}
        for idx, cid in enumerate(clusters):
            module_map.setdefault(int(cid), []).append(idx)

        merged_modules: dict[str, list[int]] = {}
        mod_counter = 1
        for cid, indices in module_map.items():
            if len(indices) >= min_module_size:
                name = f"M{mod_counter:02d}"
                merged_modules[name] = indices
                mod_counter += 1

        if not merged_modules:
            warnings.append("未检测到有效模块")

        for mod_name, indices in merged_modules.items():
            mod_genes = [valid_genes[i] for i in indices]
            mod_matrix = matrix[indices]
            eigengene = _compute_eigengene(mod_matrix)

            for i, gene in zip(indices, mod_genes):
                vec = matrix[i]
                membership = float(np.dot(vec, eigengene) / (n_samples - 1))
                membership = max(0.0, min(1.0, membership))

                candidate = CandidateGene(
                    gene_id=gene, module_id=mod_name,
                    module_membership=round(membership, 4),
                )
                evidence = Evidence(
                    evidence_type=EvidenceType.WGCNA_MODULE_MEMBERSHIP,
                    source=f"wgcna_real_{self.version}",
                    score=round(membership, 4),
                    description=f"Module={mod_name}, kME={membership:.3f}, TOM-based, method={method}",
                    metadata={
                        "module": mod_name, "module_size": len(indices),
                        "kME": round(membership, 4), "n_samples": n_samples,
                        "method": method, "soft_power": soft_power,
                    },
                )
                candidates.append(candidate)
                evidence_list.append(evidence)

        candidates.sort(key=lambda g: g.module_membership, reverse=True)
        module_sizes = {k: len(v) for k, v in merged_modules.items()}

        return ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "n_genes_input": len(gene_ids),
                "n_genes_analyzed": n_genes,
                "n_modules": len(merged_modules),
                "module_sizes": module_sizes,
                "n_samples": n_samples,
                "method": method,
                "soft_power": soft_power,
                "tool": "wgcna_real",
            },
        )
