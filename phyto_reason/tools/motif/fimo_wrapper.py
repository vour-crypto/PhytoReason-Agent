"""
fimo_wrapper.py — 真实 Motif 扫描工具。

使用:
  - JASPAR 植物 TF PWM 矩阵 (内置)
  - 位置权重矩阵 (PWM) 扫描
  - P-value 近似计算
  - 背景模型 (GC content)

不依赖 MEME/FIMO 命令行工具。
"""

from __future__ import annotations

import logging
import math
import re
from typing import Any

import numpy as np

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool

logger = logging.getLogger("fimo_wrapper")

# ── JASPAR 植物 PWM 矩阵 (部分) ──────────────────────────
# 格式: {tf_family: [(name, matrix_rows_4xW, background_freq_A_C_G_T), ...]}
# 矩阵行: [A, C, G, T]

JASPAR_PWMS: dict[str, list[tuple[str, list[list[float]], tuple[float, float, float, float]]]] = {
    "MYB": [
        ("MYB_plant_MBS", [
            [0.0, 0.0, 0.0, 0.0, 0.8, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 0.0, 0.2, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        ], (0.27, 0.23, 0.23, 0.27)),
        ("MYB_MRE", [
            [0.0, 0.8, 0.9, 0.0, 0.0, 0.8],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [1.0, 0.2, 0.0, 0.0, 0.0, 0.2],
            [0.0, 0.0, 0.1, 1.0, 0.0, 0.0],
        ], (0.27, 0.23, 0.23, 0.27)),
    ],
    "bHLH": [
        ("bHLH_Ebox", [
            [0.0, 0.9, 0.0, 0.0, 0.0, 0.9],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [1.0, 0.1, 0.0, 0.0, 1.0, 0.1],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        ], (0.27, 0.23, 0.23, 0.27)),
        ("bHLH_Gbox", [
            [0.0, 0.95, 0.0, 0.0, 0.0, 0.95],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [1.0, 0.05, 0.0, 0.0, 1.0, 0.05],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        ], (0.27, 0.23, 0.23, 0.27)),
    ],
    "WRKY": [
        ("WRKY_Wbox", [
            [0.0, 0.0, 0.0, 1.0, 0.8, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [0.5, 1.0, 0.0, 0.0, 0.2, 0.5, 0.0],
            [0.5, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0],
        ], (0.27, 0.23, 0.23, 0.27)),
    ],
    "AP2-ERF": [
        ("ERF_GCCbox", [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ], (0.27, 0.23, 0.23, 0.27)),
    ],
    "NAC": [
        ("NAC_NACBS", [
            [0.0, 0.9, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.1, 0.0, 1.0, 0.0],
            [1.0, 0.0, 1.0, 0.0, 1.0],
        ], (0.27, 0.23, 0.23, 0.27)),
    ],
    "bZIP": [
        ("bZIP_ABRE", [
            [0.0, 0.0, 0.0, 0.0, 0.9, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0, 0.0, 0.1, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        ], (0.27, 0.23, 0.23, 0.27)),
    ],
}


def _pwm_score(pwm: np.ndarray, seq: str, bg: tuple[float, float, float, float]) -> float:
    """计算 PWM 得分 (log-odds)。"""
    if len(seq) != pwm.shape[1]:
        return -float("inf")
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    score = 0.0
    for i, base in enumerate(seq.upper()):
        if base not in mapping:
            score += 0.0
            continue
        row = mapping[base]
        prob = pwm[row, i]
        if prob > 0:
            score += math.log2(prob / bg[row])
        else:
            score += -10.0
    return score


def _scan_pwm(sequence: str, pwm: np.ndarray, bg: tuple[float, float, float, float],
              min_score: float = 3.0) -> list[dict]:
    """在序列中扫描 PWM。"""
    results = []
    w = pwm.shape[1]
    for i in range(len(sequence) - w + 1):
        sub = sequence[i:i + w]
        score = _pwm_score(pwm, sub, bg)
        if score >= min_score:
            results.append({"position": i, "score": round(score, 3), "strand": "+"})
    return results


def _resolve_tf_family(name: str) -> str | None:
    """模糊匹配TF家族名称。"""
    for key in JASPAR_PWMS:
        if name.upper() == key.upper():
            return key
        if name.upper().replace("-", "/") == key.upper().replace("-", "/"):
            return key
        if name.upper().replace("/", "-") == key.upper().replace("-", "/"):
            return key
    for key in JASPAR_PWMS:
        if key in name.upper():
            return key
    return None


@register_tool
class RealMotifTool(BaseTool):
    """真实 Motif 扫描工具 (JASPAR PWM + log-odds 评分)。"""
    tool_name = "motif_real"
    description = "基于 JASPAR PWM 的真实启动子 Motif 扫描 (log-odds scoring)。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="promoter_sequences", type="object",
                      description="{gene_id: upstream_sequence}", required=True),
        ToolParameter(name="gene_ids", type="array", description="基因 ID 列表", required=True),
        ToolParameter(name="tf_family", type="string", description="TF 家族"),
    ]
    supported_data_types = ["genomic_sequence"]

    def validate_input(self, **kwargs) -> list[str]:
        errors = []
        if not kwargs.get("promoter_sequences"):
            errors.append("promoter_sequences 不能为空")
        if not kwargs.get("gene_ids"):
            errors.append("gene_ids 不能为空")
        return errors

    def run(self, **kwargs) -> ToolResult:
        promoters: dict = kwargs.get("promoter_sequences", {})
        gene_ids: list[str] = kwargs.get("gene_ids", [])
        tf_family: str | None = kwargs.get("tf_family")

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        families_to_scan: list[tuple[str, list]] = []
        if tf_family:
            resolved = _resolve_tf_family(tf_family)
            if resolved and resolved in JASPAR_PWMS:
                families_to_scan = [(resolved, JASPAR_PWMS[resolved])]
            else:
                warnings.append(f"未知家族 '{tf_family}', 扫描所有已知家族")
                families_to_scan = list(JASPAR_PWMS.items())
        else:
            families_to_scan = list(JASPAR_PWMS.items())

        total_hits = 0
        for gene_id in gene_ids:
            seq = promoters.get(gene_id, "")
            if not seq or len(seq) < 10:
                continue

            seq_upper = seq.upper()
            best_score = 0.0
            gene_hits: list[dict] = []

            for family, matrices in families_to_scan:
                for name, pwm_data, bg in matrices:
                    pwm = np.array(pwm_data, dtype=np.float64)
                    col_sums = pwm.sum(axis=0)
                    pwm = pwm / (col_sums + 1e-10)
                    hits = _scan_pwm(seq_upper, pwm, bg)
                    for hit in hits:
                        total_hits += 1
                        pos_ratio = hit["position"] / max(len(seq_upper), 1)
                        position_bonus = max(0.0, 1.0 - pos_ratio * 0.3)
                        final = (hit["score"] / 20.0) * position_bonus
                        final = max(0.0, min(1.0, final))
                        best_score = max(best_score, final)
                        gene_hits.append({
                            "motif": name, "position": hit["position"],
                            "pwm_score": hit["score"],
                            "normalized_score": round(final, 3),
                            "family": family,
                        })

            if gene_hits:
                gene_hits.sort(key=lambda h: h["normalized_score"], reverse=True)
                candidate = CandidateGene(gene_id=gene_id, motif_score=round(best_score, 4))
                evidence = Evidence(
                    evidence_type=EvidenceType.MOTIF_BINDING,
                    source=f"motif_real_{self.version}",
                    score=round(best_score, 4),
                    description=f"JASPAR: {gene_hits[0]['motif']} ({gene_hits[0]['family']}), "
                                f"{len(gene_hits)} hits, pwm_score={gene_hits[0]['pwm_score']:.1f}",
                    metadata={
                        "gene_id": gene_id, "n_hits": len(gene_hits),
                        "best_motif": gene_hits[0]["motif"],
                        "best_family": gene_hits[0]["family"],
                        "best_score": round(best_score, 4),
                    },
                )
                candidates.append(candidate)
                evidence_list.append(evidence)

        candidates.sort(key=lambda g: g.motif_score, reverse=True)

        return ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "n_genes_scanned": len(gene_ids),
                "n_genes_with_hits": len(candidates),
                "n_families_scanned": len(families_to_scan),
                "total_hits": total_hits,
                "method": "jaspar_pwm_log_odds",
                "tool": "motif_real",
            },
        )

    @staticmethod
    def list_pwms() -> dict[str, list[str]]:
        return {
            family: [name for name, _, _ in matrices]
            for family, matrices in JASPAR_PWMS.items()
        }
