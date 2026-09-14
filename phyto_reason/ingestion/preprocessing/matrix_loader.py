"""
matrix_loader.py — 矩阵加载器。

从 CSV / TSV / XLSX 加载表达矩阵和代谢物矩阵。
保留原始基因 ID。
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any

import pandas as pd


class MatrixLoader:
    """矩阵加载器。"""

    @staticmethod
    def _read_csv_any_encoding(path: Path, sep: str) -> pd.DataFrame:
        """CSV 读取：UTF-8 优先，GBK 回退（代谢物中文名常见 GBK 编码）。"""
        for enc in ("utf-8", "gbk", "gb18030"):
            try:
                return pd.read_csv(path, sep=sep, index_col=0, encoding=enc)
            except UnicodeDecodeError:
                continue
        raise RuntimeError(f"无法解码 {path}（utf-8/gbk/gb18030 均失败）")

    @staticmethod
    def load_expression_matrix(path: str | Path) -> dict[str, dict[str, float]]:
        """加载表达矩阵。首列=基因ID，首行=样本名。"""
        path = Path(path)
        ext = path.suffix.lower()

        if ext == ".xlsx":
            df = pd.read_excel(path, index_col=0)
        else:
            sep = MatrixLoader._detect_delimiter(path)
            df = MatrixLoader._read_csv_any_encoding(path, sep)

        matrix: dict[str, dict[str, float]] = {}
        for gene_id in df.index:
            row = df.loc[gene_id]
            matrix[str(gene_id)] = {
                str(col): float(v) if pd.notna(v) else float("nan")
                for col, v in row.items()
            }

        return matrix

    @staticmethod
    def load_metabolite_matrix(path: str | Path) -> dict[str, dict[str, float]]:
        """加载代谢物矩阵。首列=代谢物名，须含 Class 列。"""
        path = Path(path)
        ext = path.suffix.lower()

        if ext == ".xlsx":
            df = pd.read_excel(path, index_col=0)
        else:
            sep = MatrixLoader._detect_delimiter(path)
            df = MatrixLoader._read_csv_any_encoding(path, sep)

        # 分离 Class 列
        class_col = None
        for col in df.columns:
            if col.lower() in ("class", "type", "category", "group"):
                class_col = col
                break

        class_map: dict[str, str] = {}
        if class_col:
            class_map = df[class_col].to_dict()
            df = df.drop(columns=[class_col])

        matrix: dict[str, dict[str, float]] = {}
        for meta_id in df.index:
            row = df.loc[meta_id]
            matrix[str(meta_id)] = {
                str(col): float(v) if pd.notna(v) else float("nan")
                for col, v in row.items()
            }

        return matrix

    @staticmethod
    def load_promoter_fasta(path: str | Path) -> dict[str, str]:
        """加载 FASTA 格式启动子序列。"""
        path = Path(path)
        sequences: dict[str, str] = {}
        current_id = ""

        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_id = line[1:].split()[0]
                    sequences[current_id] = ""
                else:
                    if current_id:
                        sequences[current_id] += line.upper()

        return sequences

    @staticmethod
    def _detect_delimiter(path: Path) -> str:
        """自动检测 CSV/TSV 分隔符。

        分隔符均为 ASCII，编码无关；但文件头必须先能解码——
        UTF-8 失败按 GB18030/latin-1 回退（代谢物中文名常见 GBK 编码）。
        """
        with path.open("rb") as f:
            raw = f.read(4096)
        head = ""
        for enc in ("utf-8-sig", "gb18030", "latin-1"):
            try:
                head = raw.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        if "\t" in head:
            return "\t"
        if ";" in head:
            return ";"
        return ","
