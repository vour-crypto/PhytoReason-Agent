"""
report_builder.py — 报告 Builder。

从 RuntimeState / fusion results 构建完整科研报告。
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

from phyto_reason.reports.markdown_report import MarkdownReport


class ReportBuilder:
    """科研报告 Builder。"""

    def __init__(self) -> None:
        self._sections: list[str] = []

    def add_standard_header(self, title: str = "PhytoReason-Agent Analysis Report",
                            species: str = "", metabolite: str = "",
                            session_id: str = "") -> "ReportBuilder":
        self._sections.extend(MarkdownReport.header(title, species, metabolite, session_id))
        return self

    def add_analysis_summary(self, n_genes: int = 0, n_candidates: int = 0,
                              n_tools: int = 0, confidence: str = "") -> "ReportBuilder":
        self._sections.extend(MarkdownReport.analysis_summary(n_genes, n_candidates, n_tools, confidence))
        return self

    def add_candidates(self, candidates: list[dict]) -> "ReportBuilder":
        self._sections.extend(MarkdownReport.candidate_table(candidates))
        return self

    def add_evidence(self, gene_id: str, evidence_list: list[dict]) -> "ReportBuilder":
        self._sections.extend(MarkdownReport.evidence_section(gene_id, evidence_list))
        return self

    def add_hypothesis(self, hypothesis: str) -> "ReportBuilder":
        self._sections.extend(MarkdownReport.hypothesis(hypothesis))
        return self

    def add_validation(self, suggestions: list[str]) -> "ReportBuilder":
        self._sections.extend(MarkdownReport.validation(suggestions))
        return self

    def add_tools(self, tools: list[str]) -> "ReportBuilder":
        self._sections.extend(MarkdownReport.tools_section(tools))
        return self

    def add_section(self, title: str, content: str) -> "ReportBuilder":
        self._sections.append(f"## {title}")
        self._sections.append("")
        self._sections.append(content)
        self._sections.append("")
        return self

    def build(self) -> str:
        return "\n".join(self._sections)

    def save(self, path: str) -> str:
        content = self.build()
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        return path
