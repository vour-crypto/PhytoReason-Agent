"""
supplementary_tables.py — 补充数据表生成。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Border, Side

from phyto_reason.reports.excel_exporter import ExcelExporter


class SupplementaryTables:
    """补充数据表生成器。"""

    @staticmethod
    def export_all(
        candidates: list[dict],
        evidence_list: list[dict],
        output_dir: str,
    ) -> list[str]:
        """生成所有补充数据表。"""
        saved = []
        base = Path(output_dir)
        base.mkdir(parents=True, exist_ok=True)

        cand_path = str(base / "candidate_TFs.xlsx")
        ExcelExporter.export_candidates(candidates, cand_path)
        saved.append(cand_path)

        ev_path = str(base / "evidence_details.xlsx")
        ExcelExporter.export_evidence(evidence_list, ev_path)
        saved.append(ev_path)

        return saved
