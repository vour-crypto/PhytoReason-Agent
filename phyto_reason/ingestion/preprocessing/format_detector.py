"""
format_detector.py — 文件格式检测器。

自动检测:
  - CSV / TSV / XLSX / TXT
  - 分隔符 (comma, tab, semicolon)
  - 编码 (UTF-8, GBK, Latin-1)
  - 文件头
  - 行列数
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class FileFormatInfo(BaseModel):
    path: str = ""
    detected_format: str = "unknown"
    delimiter: str = ","
    encoding: str = "utf-8"
    has_header: bool = True
    n_rows: int = 0
    n_cols: int = 0
    columns: list[str] = Field(default_factory=list)
    is_matrix: bool = False
    warnings: list[str] = Field(default_factory=list)


class FormatDetector:
    """文件格式检测器。"""

    DELIMITERS = [",", "\t", ";", "|"]

    @staticmethod
    def detect(file_path: str | Path, sample_lines: int = 10) -> FileFormatInfo:
        """检测文件格式。"""
        path = Path(file_path)
        fmt = FileFormatInfo(path=str(path))

        ext = path.suffix.lower()
        if ext in (".xlsx", ".xls"):
            fmt.detected_format = "excel"
            return FormatDetector._detect_excel(fmt, path)

        if ext == ".fasta" or ext == ".fa":
            fmt.detected_format = "fasta"
            return fmt

        return FormatDetector._detect_text(fmt, path, sample_lines)

    @staticmethod
    def _detect_text(fmt: FileFormatInfo, path: Path, n: int) -> FileFormatInfo:
        """检测文本文件格式。"""
        # 检测编码
        for enc in ["utf-8", "gbk", "latin-1"]:
            try:
                with open(path, "r", encoding=enc) as f:
                    f.read(1000)
                fmt.encoding = enc
                break
            except (UnicodeDecodeError, UnicodeError):
                continue

        # 检测分隔符
        try:
            with open(path, "r", encoding=fmt.encoding) as f:
                head = [f.readline() for _ in range(n)]
        except Exception:
            fmt.warnings.append("Cannot read file")
            return fmt

        if not head or not head[0].strip():
            fmt.warnings.append("Empty file")
            return fmt

        best_delim = ","
        best_cols = 0
        for delim in FormatDetector.DELIMITERS:
            cols = len(head[0].split(delim))
            if cols > best_cols:
                best_cols = cols
                best_delim = delim

        fmt.delimiter = best_delim
        fmt.n_rows = sum(1 for line in head if line.strip())
        fmt.n_cols = best_cols

        try:
            reader = csv.reader(head, delimiter=best_delim)
            first_row = next(reader)
            fmt.columns = first_row
            fmt.has_header = any(c.isalpha() for c in first_row[:3])
        except Exception:
            pass

        if fmt.n_rows > 1 and fmt.n_cols > 2 and fmt.has_header:
            fmt.is_matrix = True

        fmt.detected_format = f"text/{best_delim}"
        return fmt

    @staticmethod
    def _detect_excel(fmt: FileFormatInfo, path: Path) -> FileFormatInfo:
        """检测 Excel 文件格式。"""
        try:
            from openpyxl import load_workbook
            wb = load_workbook(path, read_only=True, data_only=True)
            ws = wb.active
            if ws:
                fmt.n_rows = ws.max_row or 0
                fmt.n_cols = ws.max_column or 0
                if ws.max_row and ws.max_column:
                    first_row = [cell.value for cell in ws[1]]
                    fmt.columns = [str(c) for c in first_row if c is not None]
                    fmt.has_header = True
                    fmt.is_matrix = True
            wb.close()
        except Exception as e:
            fmt.warnings.append(f"Excel read error: {e}")

        return fmt
