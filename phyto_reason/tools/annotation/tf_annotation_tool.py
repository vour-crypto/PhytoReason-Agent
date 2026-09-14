"""
tf_annotation_tool.py — TF 基因注释工具。

识别流程:
  1. 检查 data/itak/tf_classification.txt（iTAK 输出）
  2. 检查 data/itak/tf_alignment.txt（Pfam 结构域命中）
  3. 通过内置的 TF 家族关键词进行 fallback 匹配
  4. 基因名中包含已知 TF 家族名（如 WRKY, MYB, bHLH）→ 标记为 TF

仅被标记为 TF 的基因才允许进入 tf_candidates。
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool

logger = logging.getLogger("tf_annotation_tool")

KNOWN_TF_FAMILIES: set[str] = {
    "AP2", "ARF", "ARR", "B3", "BBR", "BES1", "BZIP", "BZR",
    "C2H2", "C3H", "CAMTA", "CO-like", "CPP", "DBB", "Dof",
    "E2F-DP", "EIL", "ERF", "FAR1", "G2-like", "GATA", "GeBP",
    "GRAS", "GRF", "HB-other", "HB-PHD", "HD-ZIP", "HRT",
    "HSF", "LBD", "LFY", "LSD", "MADS", "M-type", "MYB",
    "MYB_related", "NAC", "NF-X1", "NF-YA", "NF-YB", "NF-YC",
    "NIN-like", "NLP", "NZZ", "RAV", "S1Fa-like", "SAP",
    "SBP", "SRS", "STAT", "TALE", "TCP", "Trihelix",
    "VOZ", "Whirly", "WOX", "WRKY", "YABBY", "ZF-HD",
    "bHLH", "bZIP",
}

TF_NAME_PATTERNS: list[re.Pattern] = [
    re.compile(rf"^{fam}(?:\d|_|\-|\.)", re.IGNORECASE)
    for fam in sorted(KNOWN_TF_FAMILIES, key=len, reverse=True)
] + [
    re.compile(r"(transcription.factor|TF_|_TF)"),
    re.compile(r"^(WRKY|MYB|bHLH|NAC|ERF|TCP|GRAS|ARF|AUX|IAA)", re.IGNORECASE),
]


@register_tool
class TFAnnotationTool(BaseTool):
    tool_name = "tf_annotation"
    description = "识别基因列表中的转录因子（基于 iTAK + PFAM + 名称模式）。"
    version = "1.0.0"
    parameters = [
        ToolParameter(name="gene_ids", type="array", description="待识别的基因 ID 列表", required=True),
    ]
    supported_data_types = ["annotation"]

    def validate_input(self, **kwargs) -> list[str]:
        if not kwargs.get("gene_ids"):
            return ["gene_ids 不能为空"]
        return []

    def run(self, **kwargs) -> ToolResult:
        gene_ids: list[str] = kwargs.get("gene_ids", [])
        candidates: list[CandidateGene] = []
        evidence_list: list[Evidence] = []
        warnings: list[str] = []

        itak_genes: set[str] = set()
        pfam_genes: dict[str, int] = {}

        itak_file = Path("data/itak/tf_classification.txt")
        if itak_file.exists():
            try:
                import pandas as pd
                df = pd.read_csv(itak_file, sep="\t", header=None)
                for raw in df[0].dropna().astype(str).str.strip():
                    itak_genes.add(raw)
                    clean = re.sub(r"\.\d+$", "", raw)
                    itak_genes.add(clean)
            except Exception as e:
                warnings.append(f"iTAK classification 加载失败: {e}")

        aln_file = Path("data/itak/tf_alignment.txt")
        if aln_file.exists():
            try:
                import pandas as pd
                df = pd.read_csv(aln_file, sep="\t", header=None, low_memory=False)
                pfam_mask = df[1].astype(str).str.match(r"^PF\d+")
                for raw_id in df[pfam_mask][0].dropna().astype(str).str.strip():
                    pfam_genes[raw_id] = pfam_genes.get(raw_id, 0) + 1
                    clean = re.sub(r"\.\d+$", "", raw_id)
                    pfam_genes[clean] = pfam_genes.get(clean, 0) + 1
            except Exception as e:
                warnings.append(f"iTAK alignment 加载失败: {e}")

        for gid in gene_ids:
            is_tf = False
            tf_family: str | None = None
            score = 0.0
            evidence_parts: list[str] = []

            if gid in itak_genes:
                is_tf = True
                score = max(score, 0.8)
                evidence_parts.append("iTAK_classification")

            if gid in pfam_genes:
                is_tf = True
                score = max(score, 0.7)
                evidence_parts.append(f"Pfam_domain({pfam_genes[gid]}hits)")

            if not is_tf:
                for pat in TF_NAME_PATTERNS:
                    if pat.search(gid):
                        is_tf = True
                        score = max(score, 0.5)
                        evidence_parts.append(f"name_pattern({pat.pattern[:20]})")
                        break

                gid_upper = gid.upper()
                for fam in sorted(KNOWN_TF_FAMILIES, key=len, reverse=True):
                    if fam.upper() in gid_upper:
                        tf_family = fam.capitalize()
                        break

            if is_tf:
                candidate = CandidateGene(
                    gene_id=gid, is_tf=True, tf_family=tf_family,
                    literature_score=round(score, 4),
                )
                ev = Evidence(
                    evidence_type=EvidenceType.ANNOTATION,
                    source="tf_annotation_tool",
                    score=round(score, 4),
                    description="TF identified: " + ", ".join(evidence_parts),
                    metadata={"evidence": evidence_parts},
                )
                candidate.add_evidence(ev)
                evidence_list.append(ev)
                candidates.append(candidate)

        return ToolResult(
            candidates=candidates,
            evidence_list=evidence_list,
            warnings=warnings,
            metadata={
                "n_genes_input": len(gene_ids),
                "n_tf_identified": len(candidates),
                "n_itak": len(itak_genes & set(gene_ids)),
                "n_pfam": len(pfam_genes),
            },
        )
