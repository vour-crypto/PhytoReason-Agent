"""
config/settings.py — 项目集中配置。

v5.0 重建：旧 config/ 包在迁移中丢失，导致 11 个模块
`from config.settings import Settings` 运行时崩溃。
本模块按现存消费方（tools/rag.py、tools/literature/pubmed_client.py、CLI）
提供所需属性，全部带默认值；敏感项优先从环境变量读取。
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from phyto_reason.platform_paths import cache_dir, output_dir, report_dir


@dataclass
class Settings:
    # ── 路径 ─────────────────────────────────────────────
    project_root: Path = field(
        default_factory=lambda: Path(__file__).resolve().parent.parent.parent
    )
    # 用户产物（分析输出等，会随备份导出）；缓存另见 cache_dir
    output_dir: Path = field(default_factory=output_dir)
    cache_dir: Path = field(default_factory=cache_dir)
    report_dir: Path = field(default_factory=report_dir)

    # ── NCBI E-utilities / PubMed ────────────────────────
    pubmed_max_retries: int = 3
    pubmed_timeout: int = 20
    pubmed_request_interval: float = 0.35  # 秒（无 API key 限速 3 次/秒）
    ncbi_tool: str = "PhytoReason-Agent"
    ncbi_email: str = os.getenv("NCBI_EMAIL", "phyto@example.com")
    ncbi_api_key: str = os.getenv("NCBI_API_KEY", "")

    # ── 物种相关（旧 species_rules 结构兼容）─────────────
    species_rules: dict = field(default_factory=dict)

    # ── 组织别名（用于样本分组映射）──────────────────────
    tissue_aliases: dict[str, list[str]] = field(default_factory=dict)

    # ── 统计阈值（CLI 单组学快捷命令用）─────────────────
    min_abs_correlation: float = 0.7
