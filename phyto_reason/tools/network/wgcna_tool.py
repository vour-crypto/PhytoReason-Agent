"""
wgcna_tool.py — WGCNA 共表达网络分析工具。

简化版 WGCNA: 基于相关性的模块检测。
实际生产环境建议使用 R/WGCNA 包。
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

logger = logging.getLogger("wgcna_tool")


@register_tool
class WGCNATool(BaseTool):
    tool_name = "wgcna"
    description = "构建加权基因共表达网络，检测共表达模块并计算模块-代谢物相关性。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="expression_matrix", type="object",
                      description="表达矩阵 (dict of {gene_id: {sample: value}})", required=True),
        ToolParameter(name="gene_ids", type="array", description="要分析的基因 ID 列表", required=True),
        ToolParameter(name="metabolite_matrix", type="object",
                      description="代谢物矩阵 (可选, 用于计算模块-代谢物相关性)", default=None),
        ToolParameter(name="soft_power", type="number", description="软阈值 power", default=6),
        ToolParameter(name="min_module_size", type="integer", description="最小模块大小", default=10),
        ToolParameter(name="merge_cut_height", type="number", description="模块合并阈值", default=0.25),
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
        expr_mat: dict = kwargs.get("expression_matrix", {})
        gene_ids: list[str] = kwargs.get("gene_ids", [])
        meta_mat: dict | None = kwargs.get("metabolite_matrix")
        soft_power: int = kwargs.get("soft_power", 6)
        min_module_size: int = kwargs.get("min_module_size", 10)

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        available_genes = [g for g in gene_ids if g in expr_mat]
        if len(available_genes) < min_module_size:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=[f"可用基因数 {len(available_genes)} 少于最小模块大小 {min_module_size}"],
                metadata={"n_genes": len(available_genes), "n_modules": 0},
            )

        # 构建表达矩阵
        sample_set: set[str] | None = None
        for samples in expr_mat.values():
            sample_set = set(samples.keys())
            break
        if sample_set is None:
            return ToolResult(
                candidates=[], evidence_list=[],
                warnings=["表达矩阵格式错误"],
                metadata={"n_genes": 0, "n_modules": 0},
            )

        common_samples = sorted(sample_set)
        n_samples = len(common_samples)
        n_genes = len(available_genes)

        gene_vecs: dict[str, np.ndarray] = {}
        for g in available_genes:
            vec = np.array([expr_mat[g].get(s, 0.0) for s in common_samples], dtype=float)
            gene_vecs[g] = (vec - np.mean(vec)) / (np.std(vec) + 1e-10)

        # 简化模块检测: 基于相关性的层次聚类
        gene_list = available_genes
        n_genes = len(gene_list)

        corr_matrix = np.zeros((n_genes, n_genes))
        for i in range(n_genes):
            for j in range(i + 1, n_genes):
                v1 = gene_vecs[gene_list[i]]
                v2 = gene_vecs[gene_list[j]]
                r = np.dot(v1, v2) / (n_samples - 1)
                r = max(-1.0, min(1.0, r))
                # 加权邻接矩阵: |corr|^soft_power
                adj = abs(r) ** soft_power
                corr_matrix[i, j] = adj
                corr_matrix[j, i] = adj

        # 简单模块分配: 基于邻接矩阵的行和聚类
        row_sums = np.sum(corr_matrix, axis=1)
        # 按连通性排序, 划分模块
        order = np.argsort(row_sums)[::-1]

        assigned = set()
        modules: dict[str, list[str]] = {}
        module_id = 1

        for idx in order:
            g = gene_list[idx]
            if g in assigned:
                continue
            # 找到与该基因高度连接的基因
            module_genes = [g]
            assigned.add(g)
            for j in range(n_genes):
                neighbor = gene_list[j]
                if neighbor in assigned:
                    continue
                if corr_matrix[idx, j] >= 0.5 ** soft_power:
                    module_genes.append(neighbor)
                    assigned.add(neighbor)

            if len(module_genes) >= min_module_size:
                mod_name = f"M{module_id:02d}"
                modules[mod_name] = module_genes
                module_id += 1

        if not modules:
            warnings.append("未检测到满足最小模块大小的共表达模块")
            return ToolResult(
                candidates=[], evidence_list=[], warnings=warnings,
                metadata={"n_genes": n_genes, "n_modules": 0},
            )

        # 计算每个模块的 eigenspace (第一主成分)
        for mod_name, mod_genes in modules.items():
            mod_mat = np.array([gene_vecs[g] for g in mod_genes])
            # 第一主成分作为模块特征基因
            u, s, vt = np.linalg.svd(mod_mat, full_matrices=False)
            eigengene = u[0] if u.shape[0] > 0 else np.zeros(n_samples)

            method = "svd" if len(mod_genes) > 1 else "single"

            for gene in mod_genes:
                vec = gene_vecs[gene]
                membership = abs(np.dot(vec, eigengene) / (n_samples - 1))
                membership = min(1.0, max(0.0, membership))

                candidate = CandidateGene(
                    gene_id=gene,
                    module_id=mod_name,
                    module_membership=round(membership, 4),
                )

                evidence = Evidence(
                    evidence_type=EvidenceType.WGCNA_MODULE_MEMBERSHIP,
                    source=f"wgcna_tool_{self.version}",
                    score=round(membership, 4),
                    confidence=round(membership * 0.9, 4),
                    description=f"Module={mod_name}, membership={membership:.3f}, n_genes={len(mod_genes)}",
                    metadata={
                        "module": mod_name,
                        "module_size": len(mod_genes),
                        "membership": round(membership, 4),
                        "n_samples": n_samples,
                        "method": method,
                    },
                )

                candidates.append(candidate)
                evidence_list.append(evidence)

        # 如果提供了代谢物矩阵, 计算模块-代谢物相关性
        module_meta_corrs = {}
        if meta_mat:
            for mod_name, mod_genes in modules.items():
                mod_vec = np.mean([gene_vecs[g] for g in mod_genes], axis=0)
                for meta_name, meta_samples in meta_mat.items():
                    meta_vec = np.array([meta_samples.get(s, 0.0) for s in common_samples], dtype=float)
                    if np.std(meta_vec) < 1e-10:
                        continue
                    meta_vec = (meta_vec - np.mean(meta_vec)) / (np.std(meta_vec) + 1e-10)
                    r = np.dot(mod_vec, meta_vec) / (n_samples - 1)
                    key = f"{mod_name}↔{meta_name}"
                    module_meta_corrs[key] = round(r, 4)

        candidates.sort(key=lambda g: g.module_membership, reverse=True)

        result = ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "n_genes_input": len(gene_ids),
                "n_genes_analyzed": n_genes,
                "n_modules": len(modules),
                "module_sizes": {k: len(v) for k, v in modules.items()},
                "n_samples": n_samples,
                "soft_power": soft_power,
                "module_metabolite_correlations": module_meta_corrs,
            },
        )
        return result
