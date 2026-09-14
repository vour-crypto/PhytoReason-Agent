"""
markdown_report.py — Markdown 格式报告生成。
"""

from __future__ import annotations

from datetime import datetime


class MarkdownReport:
    """Markdown 格式科研报告。"""

    @staticmethod
    def header(title: str, species: str = "", metabolite: str = "",
               session_id: str = "") -> list[str]:
        lines = [
            f"# {title}",
            "",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        ]
        if session_id:
            lines.append(f"**Session:** {session_id}")
        if species:
            lines.append(f"**Species:** {species}")
        if metabolite:
            lines.append(f"**Target Metabolite:** {metabolite}")
        lines.append("")
        lines.append("---")
        lines.append("")
        return lines

    @staticmethod
    def analysis_summary(n_genes: int = 0, n_candidates: int = 0,
                          n_tools: int = 0, confidence: str = "") -> list[str]:
        return [
            "## Analysis Summary",
            "",
            f"| Metric | Value |",
            f"|---|---|",
            f"| Genes Analyzed | {n_genes} |",
            f"| Candidate TFs | {n_candidates} |",
            f"| Tools Executed | {n_tools} |",
            f"| Overall Confidence | {confidence} |",
            "",
        ]

    @staticmethod
    def candidate_table(candidates: list[dict]) -> list[str]:
        lines = [
            "## Top Candidate TFs",
            "",
            "| Rank | Gene ID | Family | Confidence | Score | Evidences |",
            "|---|---|---|---|---|---|",
        ]
        for i, c in enumerate(candidates[:20], 1):
            lines.append(
                f"| {i} | {c.get('gene_id', '')} | "
                f"{c.get('tf_family', '-')} | "
                f"{c.get('confidence_level', '')} | "
                f"{c.get('confidence_score', 0):.3f} | "
                f"{c.get('n_evidences', 0)} |"
            )
        lines.append("")
        return lines

    @staticmethod
    def evidence_section(gene_id: str, evidence_list: list[dict]) -> list[str]:
        lines = [f"### {gene_id} Evidence", ""]
        for ev in evidence_list[:5]:
            lines.append(
                f"- **{ev.get('source', '?')}**: "
                f"score={ev.get('normalized_score', 0):.2f}"
            )
        lines.append("")
        return lines

    @staticmethod
    def hypothesis(hypothesis_text: str) -> list[str]:
        return ["## Biological Hypothesis", "", hypothesis_text, ""]

    @staticmethod
    def validation(suggestions: list[str]) -> list[str]:
        lines = ["## Suggested Wet-Lab Validation", ""]
        for s in suggestions[:5]:
            lines.append(f"- {s}")
        lines.append("")
        return lines

    @staticmethod
    def dataset_table(datasets: list[dict]) -> list[str]:
        lines = [
            "## Datasets Used",
            "",
            "| File | Type | Size |",
            "|---|---|---|",
        ]
        for d in datasets:
            lines.append(
                f"| {d.get('filename', '')} | "
                f"{d.get('data_type', '')} | "
                f"{d.get('size', '-')} |"
            )
        lines.append("")
        return lines

    @staticmethod
    def tools_section(tools: list[str]) -> list[str]:
        return [
            "## Tools Executed",
            "",
            *[f"- {t}" for t in tools],
            "",
        ]
