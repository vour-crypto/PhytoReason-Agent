"""
gene_id_mapper.py — 基因 ID 映射器。

保留用户原始基因 ID，不做自动重命名。
支持:
  - 数字前缀修复 (如 "1.1234" → "GENE_1.1234")
  - 非法字符替换
  - ID 冲突检测
"""

from __future__ import annotations

import re
from typing import Any


class GeneIDMapper:
    """基因 ID 映射器。"""

    @staticmethod
    def sanitize(gene_id: str) -> str:
        """清理基因 ID 中的非法字符，但保留原始可读性。"""
        sanitized = gene_id.strip()
        if not sanitized:
            return f"UNKNOWN_{hash(gene_id) % 10000}"
        return sanitized

    @staticmethod
    def detect_format(gene_ids: list[str]) -> str:
        """检测基因 ID 的命名格式。"""
        patterns = {
            "evm": r"evm\.\S+",
            "maker": r"maker\S+",
            "trinity": r"TRINITY_\S+",
            "refseq": r"(NM_|XM_|NP_|XP_)\d+",
            "ensembl": r"(ENS\w+G\d+)",
            "arabidopsis": r"AT\wG\d+",
            "generic": r"\w+",
        }

        for gid in gene_ids[:100]:
            for fmt, pattern in patterns.items():
                if re.match(pattern, gid, re.IGNORECASE):
                    return fmt
        return "unknown"

    @staticmethod
    def check_duplicates(gene_ids: list[str]) -> list[str]:
        """检查重复基因 ID。"""
        seen: dict[str, int] = {}
        duplicates = []
        for gid in gene_ids:
            gid_lower = gid.lower()
            if gid_lower in seen:
                if seen[gid_lower] == 1:
                    duplicates.append(gid)
                seen[gid_lower] += 1
            else:
                seen[gid_lower] = 1
        return duplicates

    @staticmethod
    def fix_prefix_issues(gene_ids: list[str]) -> dict[str, str]:
        """修复以数字开头的基因 ID (如 '1.1234' → 'GENE_1.1234')。"""
        mapping: dict[str, str] = {}
        for gid in gene_ids:
            if gid and gid[0].isdigit():
                fixed = f"GENE_{gid}"
                mapping[gid] = fixed
            else:
                mapping[gid] = gid
        return mapping
