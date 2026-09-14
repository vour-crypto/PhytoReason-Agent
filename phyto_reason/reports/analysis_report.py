"""analysis_report.py — 会话级分析摘要与 SCI 式报告渲染。

设计要点：
- 分节只覆盖**实际完成的节点**（completed_nodes 且报告非空）——
  只做差异分析就只写差异分析，不强行输出 TF/假设章节；
- 图表按生成顺序编号（Figure 1..N），编号持久在会话上，
  报告与结果面板的编号一致；
- 每图附双语题注与"如何解读"提示，面向非生信读者。
"""

from __future__ import annotations

import math
from datetime import datetime
from pathlib import Path

# 前缀 → (英文名, 中文名, 解读)
FIGURE_CAPTIONS: dict[str, tuple[str, str, str]] = {
    "qc_heatmap": ("QC clustered heatmap of top-variable genes",
                   "高变基因聚类热图",
                   "看基因在不同样本间的表达聚集成块：同组样本聚在一起说明组间差异清晰。"),
    "qc_features": ("Detected features per group", "各组检出特征数",
                    "各组柱高接近说明检测量均一；差异过大多为上样/测序量不均。"),
    "sample_correlation": ("Sample correlation heatmap", "样本相关性热图",
                           "对角线亮、组内亮组间暗为正常；组间也全亮说明分组效应弱。"),
    "pca": ("PCA score plot", "PCA 得分图",
            "同组点抱团、组间分开说明数据可靠；第一主轴常代表主要处理/组织效应。"),
    "volcano_deg": ("DEG volcano plot", "差异表达火山图",
                    "右上=显著上调、左上=显著下调基因；横轴越大变化倍数越大，纵轴越高越显著。"),
    "volcano_dam": ("DAM volcano plot", "差异代谢火山图",
                    "同火山图：右上显著上调代谢物、左上显著下调。"),
    "dam_pairwise_heatmap": ("DAM pairwise log2FC heatmap", "差异代谢物组织对热图",
                             "每格是一对组织的平均 log2FC；全部接近 0 说明组织两两间无显著差异。"),
    "correlation_network": ("Gene–metabolite correlation network", "基因-代谢物相关网络",
                            "点=基因/代谢物，线=强相关；线的数量与枢纽点提示核心调控分子。"),
    "quadrant": ("DEG × DAM quadrant plot", "差异基因×差异代谢物象限图",
                 "第一象限=基因与代谢物同向上调（协同）；Q1+Q3 占比高说明转录与代谢协同。"),
    "pathway_enrichment": ("Joint KEGG enrichment", "联合 KEGG 富集",
                           "条越长该通路被差异分子富集得越显著（-log10 q）。"),
    "wgcna_modules": ("WGCNA module sizes", "WGCNA 模块大小",
                      "彩色柱为共表达模块，grey 是未能归组的基因；关注与性状相关高的彩色模块。"),
    "wgcna_scale_free": ("WGCNA scale-free topology fit", "WGCNA 无标度拓扑拟合",
                         "软阈值 power 的选择曲线：R² 达到 ~0.8 以上的最小 power 最合适。"),
}

# 完成节点 → 报告分节（按报告出现顺序）
NODE_SECTIONS: list[tuple[str, str, str, str]] = [
    ("data_qc", "数据质控", "Quality control", "qc_report"),
    ("deg_analysis", "差异表达分析（DEG）", "Differential expression analysis", "deg_report"),
    ("dam_analysis", "差异代谢物分析（DAM）", "Differential metabolite analysis", "dam_report"),
    ("multiomics_integration", "基因-代谢物多组学关联", "Multi-omics correlation", "multiomics_report"),
    ("quadrant_plot", "DEG × DAM 象限分析", "Quadrant analysis", "quadrant_plot_report"),
    ("correlation_network", "基因-代谢物相关网络", "Correlation network", "correlation_network_report"),
    ("joint_enrichment", "联合 KEGG 富集", "Joint KEGG enrichment", "joint_enrichment_report"),
    ("wgcna", "WGCNA 共表达模块", "WGCNA co-expression modules", "wgcna_report"),
    ("o2pls", "O2PLS 联合成分分析", "O2PLS joint components", "o2pls_report"),
]

_FIGURE_KEYS = ("figure_url", "pairwise_heatmap_url")


def _slim_stats(report: dict) -> dict:
    out = {}
    for key, value in (report or {}).items():
        if isinstance(value, bool) or isinstance(value, int) or isinstance(value, float):
            out[key] = value
        elif isinstance(value, str) and 0 < len(value) <= 220:
            out[key] = value
    return out


def _report_figures(report: dict) -> list[str]:
    urls: list[str] = []
    for key in _FIGURE_KEYS:
        value = (report or {}).get(key)
        if isinstance(value, str) and value.startswith("/static/figures/"):
            urls.append(value)
    for value in (report or {}).get("figure_urls", []) or []:
        if isinstance(value, str) and value.startswith("/static/figures/"):
            urls.append(value)
    return urls


def _figure_prefix(url: str) -> str:
    stem = Path(url.replace("\\", "/")).stem  # volcano_deg_a27b08
    for prefix in sorted(FIGURE_CAPTIONS, key=len, reverse=True):
        if stem == prefix or stem.startswith(prefix + "_"):
            return prefix
    return stem.rsplit("_", 1)[0] if "_" in stem else stem


class FigureNumberer:
    """会话级图表编号器：同名文件编号稳定，新图递增。"""

    def __init__(self, registry: dict | None = None) -> None:
        self.registry = registry if registry is not None else {}
        self.next_no = 1 + max((int(v.get("no", 0)) for v in self.registry.values()), default=0)

    def assign(self, url: str) -> dict:
        name = Path(url.replace("\\", "/")).name
        entry = self.registry.get(name)
        if entry:
            return entry
        prefix = _figure_prefix(url)
        en, zh, howto = FIGURE_CAPTIONS.get(
            prefix, (prefix, prefix, "图表由对应分析生成，具体解读见该节文字。"))
        entry = {"no": self.next_no, "file": name, "url": url,
                 "caption_en": en, "caption_zh": zh, "how_to_read": howto,
                 "prefix": prefix}
        self.next_no += 1
        self.registry[name] = entry
        return entry


def build_analysis_digest(state, species: str = "", target_metabolite: str = "",
                          figure_registry: dict | None = None) -> dict:
    """从 RuntimeState（或同形对象）构建报告摘要；只收录实际完成的分节。"""
    completed = set(getattr(state, "completed_nodes", []) or [])
    numberer = FigureNumberer(figure_registry)
    sections: list[dict] = []
    for node, title_zh, title_en, report_key in NODE_SECTIONS:
        if node not in completed:
            continue
        report = getattr(state, report_key, None)
        if not report:
            continue
        urls = _report_figures(report)
        figures = [numberer.assign(url) for url in urls]
        sections.append({
            "node": node, "title_zh": title_zh, "title_en": title_en,
            "stats": _slim_stats(report), "figures": figures,
        })

    figure_urls: list[str] = []
    for section in sections:
        for figure in section["figures"]:
            if figure["url"] not in figure_urls:
                figure_urls.append(figure["url"])

    return {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "species": species,
        "target_metabolite": target_metabolite,
        "mechanism_type": str(getattr(state, "mechanism_type", "") or ""),
        "sections": sections,
        "figure_urls": figure_urls,
        "n_competing_hypotheses": len(getattr(state, "competing_hypotheses", []) or []),
        "n_tf_candidates": len(getattr(state, "tf_candidates", {}) or {}),
    }


def _fmt(value) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "—"
        return f"{value:g}"
    return str(value)


def render_markdown(digest: dict) -> str:
    """摘要 → SCI Results 风格中文 Markdown（图号引用 + 图注 + 解读）。"""
    lines: list[str] = []
    species = digest.get("species") or "未设置物种"
    target = digest.get("target_metabolite") or "未指定"
    lines.append(f"# 分析报告：{species} · 目标代谢物 {target}")
    lines.append("")
    lines.append(f"- 生成时间：{digest.get('created_at', '')}")
    lines.append(f"- 本次完成的分析：{digest.get('completed_summary', '') or '、'.join(s['title_zh'] for s in digest['sections'])}")
    lines.append("")

    if not digest["sections"]:
        lines.append("本次运行没有产生可报告的分析节。")
        return "\n".join(lines)

    for index, section in enumerate(digest["sections"], start=1):
        lines.append(f"## {index}. {section['title_zh']}（{section['title_en']}）")
        lines.append("")
        stats = section["stats"]
        if stats.get("summary"):
            lines.append(f"方法与阈值：{stats['summary']}")
        if stats.get("normalization"):
            lines.append(f"预处理：{stats['normalization']}")
        highlight_keys = [
            "n_genes_tested", "n_fdr_significant", "n_metabolites_tested",
            "n_pairwise_tests", "n_pairwise_significant", "n_significant_pairs",
            "n_correlation_pairs_tested", "n_pairs", "synchronicity_ratio",
            "n_modules", "soft_power", "scale_free_r2",
            "n_joint_components", "x_joint_variance_explained",
            "y_joint_variance_explained", "n_significant_pathways",
        ]
        bullets = [f"- {key}: {_fmt(stats[key])}" for key in highlight_keys if key in stats]
        if bullets:
            lines.extend(bullets)
        if stats.get("interpretation"):
            lines.append(f"解读：{stats['interpretation']}")
        lines.append("")
        for figure in section["figures"]:
            label = f"Figure {figure['no']}"
            lines.append(f"**{label}** {figure['caption_en']}（{figure['caption_zh']}）")
            lines.append("")
            lines.append(f"![{label} {figure['caption_en']}]({figure['url']})")
            lines.append("")
            lines.append(f"> 如何解读：{figure['how_to_read']}")
            lines.append("")
    return "\n".join(lines)
