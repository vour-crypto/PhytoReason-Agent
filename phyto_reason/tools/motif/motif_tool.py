"""
motif_tool.py — 启动子 Motif 扫描工具。

扫描基因上游启动子区域中的已知 TF 结合 motif。
"""

from __future__ import annotations

import logging
from datetime import datetime

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool

logger = logging.getLogger("motif_tool")


# ── 已知的 TF 家族结合 motif (IUPAC 简并碱基) ──────────────
# 来源: JASPAR / PlantCARE / cisBP

KNOWN_MOTIFS: dict[str, list[dict]] = {
    "MYB": [
        {"motif": "C[AC]NG[GT][TG]A", "name": "MYB core (MBS)", "score": 0.85},
        {"motif": "TAAC[GT]A", "name": "MYB recognition element (MRE)", "score": 0.80},
    ],
    "bHLH": [
        {"motif": "CANNTG", "name": "E-box (bHLH binding)", "score": 0.85},
        {"motif": "CACGTG", "name": "G-box (bHLH binding)", "score": 0.80},
    ],
    "WRKY": [
        {"motif": "[TG]TGAC[C/T]", "name": "W-box (WRKY binding)", "score": 0.90},
    ],
    "AP2-ERF": [
        {"motif": "GCCGCC", "name": "GCC-box (ERF binding)", "score": 0.85},
        {"motif": "AGCCGCC", "name": "GCC-box core", "score": 0.80},
    ],
    "NAC": [
        {"motif": "C[AC]GT[GA]", "name": "NAC binding site (NACBS)", "score": 0.75},
    ],
    "bZIP": [
        {"motif": "ACGT[GC]A", "name": "ABRE (bZIP binding)", "score": 0.80},
        {"motif": "CACGTG", "name": "G-box (bZIP)", "score": 0.80},  # overlaps with bHLH
    ],
    "C2H2": [
        {"motif": "AC[AG]T[GC]A", "name": "Zinc finger binding", "score": 0.65},
    ],
}


# 简单的 IUPAC 模式匹配
def _iupac_to_regex(iupac: str) -> str:
    """将 IUPAC 简并碱基转为正则表达式。"""
    mapping = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
    }
    result = ""
    for ch in iupac:
        result += mapping.get(ch.upper(), ch)
    return result


@register_tool
class MotifTool(BaseTool):
    tool_name = "motif_scan"
    description = "扫描基因上游启动子区域的已知 TF 结合 motif, 返回匹配结果。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="promoter_sequences", type="object",
                      description="启动子序列 (dict of {gene_id: upstream_sequence})", required=True),
        ToolParameter(name="tf_family", type="string", description="要扫描的 TF 家族 (如 'MYB', 'WRKY')"),
        ToolParameter(name="gene_ids", type="array", description="要扫描的基因 ID 列表", required=True),
        ToolParameter(name="upstream_length", type="integer", description="上游区域长度", default=1500),
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
        promoter_sequences: dict = kwargs.get("promoter_sequences", {})
        tf_family: str | None = kwargs.get("tf_family")
        gene_ids: list[str] = kwargs.get("gene_ids", [])

        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        motifs_to_check: list[tuple[str, list[dict]]] = []
        if tf_family:
            family_motifs = KNOWN_MOTIFS.get(tf_family)
            if not family_motifs:
                resolved = self._resolve_tf_family(tf_family)
                if resolved:
                    motifs_to_check = [(resolved, KNOWN_MOTIFS[resolved])]
                else:
                    warnings.append(f"未知 TF 家族: '{tf_family}', 使用所有已知 motif")
                    motifs_to_check = list(KNOWN_MOTIFS.items())
            else:
                motifs_to_check = [(tf_family, family_motifs)]
        else:
            motifs_to_check = list(KNOWN_MOTIFS.items())

        total_hits = 0
        for gene_id in gene_ids:
            seq = promoter_sequences.get(gene_id, "")
            if not seq or len(seq) < 10:
                warnings.append(f"基因 {gene_id}: 无启动子序列")
                continue

            seq_upper = seq.upper()
            gene_motif_score = 0.0
            gene_hits: list[dict] = []

            for family, motifs in motifs_to_check:
                for m in motifs:
                    pattern = m["motif"]
                    regex = _iupac_to_regex(pattern)
                    import re
                    matches = list(re.finditer(regex, seq_upper))
                    if matches:
                        for match in matches:
                            base_score = m["score"]
                            # 在启动子中位置越靠近转录起始位点 → 分数越高
                            pos_ratio = match.start() / max(len(seq_upper), 1)
                            position_bonus = max(0.0, 1.0 - pos_ratio * 0.5)
                            final_score = base_score * position_bonus
                            gene_motif_score = max(gene_motif_score, final_score)
                            total_hits += 1
                            gene_hits.append({
                                "motif": m["name"],
                                "pattern": pattern,
                                "position": match.start(),
                                "score": round(final_score, 3),
                            })

            if gene_hits:
                candidate = CandidateGene(
                    gene_id=gene_id,
                    motif_score=round(gene_motif_score, 4),
                )

                evidence = Evidence(
                    evidence_type=EvidenceType.MOTIF_BINDING,
                    source=f"motif_tool_{self.version}",
                    score=round(gene_motif_score, 4),
                    confidence=round(gene_motif_score * 0.85, 4),
                    description=f"Found {len(gene_hits)} motif hits in promoter (top: {gene_hits[0]['motif']})",
                    metadata={
                        "gene_id": gene_id,
                        "n_motif_hits": len(gene_hits),
                        "best_motif": gene_hits[0]["motif"],
                        "best_score": round(gene_motif_score, 4),
                        "all_hits": gene_hits[:10],
                    },
                )

                candidates.append(candidate)
                evidence_list.append(evidence)

        candidates.sort(key=lambda g: g.motif_score, reverse=True)

        result = ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "n_genes_scanned": len(gene_ids),
                "n_genes_with_hits": len(candidates),
                "n_motif_families": len(motifs_to_check),
                "n_motifs_tested": sum(len(m) for _, m in motifs_to_check),
                "total_hits": total_hits,
                "tf_family_filter": tf_family,
            },
        )
        return result

    @staticmethod
    def _resolve_tf_family(name: str) -> str | None:
        """模糊匹配TF家族名称。"""
        for key in KNOWN_MOTIFS:
            if name.upper() == key.upper():
                return key
            if name.upper().replace("-", "/") == key.upper().replace("-", "/"):
                return key
            if name.upper().replace("/", "-") == key.upper().replace("-", "/"):
                return key
        return None

    @staticmethod
    def list_known_motifs() -> dict[str, list[str]]:
        """列出所有已知的 TF 家族及其 motif 名称。"""
        return {
            family: [m["name"] for m in motifs]
            for family, motifs in KNOWN_MOTIFS.items()
        }
