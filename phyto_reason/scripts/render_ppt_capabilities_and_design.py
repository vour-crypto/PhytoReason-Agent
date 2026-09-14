"""Render presentation diagrams for core capabilities and design principles."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
FONT_REGULAR = Path(r"C:\Windows\Fonts\msyh.ttc")
FONT_BOLD = Path(r"C:\Windows\Fonts\msyhbd.ttc")


def box(ax, x, y, w, h, *, face, edge, linewidth=1.45, radius=0.12, zorder=3):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=linewidth, zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, *, color="#748079", linewidth=1.5, rad=0.0, zorder=2):
    patch = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=13,
        color=color, linewidth=linewidth,
        connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def title_block(ax, title, subtitle, *, bold, regular, ink, muted):
    ax.text(
        0.55, 8.55, title,
        ha="left", va="center", color=ink,
        fontsize=24, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        0.57, 8.14, subtitle,
        ha="left", va="center", color=muted,
        fontsize=11.5, fontproperties=regular,
    )


def capability_card(ax, x, y, number, title, lines, *, face, edge, regular, bold):
    box(ax, x, y, 4.20, 1.70, face=face, edge=edge, linewidth=1.55)
    ax.text(
        x + 0.28, y + 1.35, number,
        ha="left", va="center", color=edge,
        fontsize=17, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        x + 0.83, y + 1.35, title,
        ha="left", va="center", color="#25312C",
        fontsize=12.1, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        x + 0.30, y + 0.72, lines,
        ha="left", va="center", color="#52605A",
        fontsize=9.8, fontproperties=regular, linespacing=1.45,
    )


def render_capabilities(regular, bold, colors):
    fig, ax = plt.subplots(figsize=(16, 9), facecolor=colors["bg"])
    ax.set_facecolor(colors["bg"])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis("off")

    title_block(
        ax,
        "PhytoReason-Agent 六大核心能力",
        "覆盖从公共知识检索、组学分析到机制假设与实验建议的完整研究链条",
        bold=bold, regular=regular, ink=colors["ink"], muted=colors["muted"],
    )

    cards = [
        (0.75, 5.85, "01", "多源科研知识检索",
         "PubMed 等文献源 · KEGG / PlantCyc\nPlantTFDB / JASPAR · 化合物数据库",
         colors["blue_face"], colors["blue_edge"]),
        (5.90, 5.85, "02", "跨物种知识迁移",
         "物种 Profile → 同源推断 → 显式缺口\n支持模式植物与非模式药用植物",
         colors["green_face"], colors["green_edge"]),
        (11.05, 5.85, "03", "多组学统计分析",
         "QC · DEG / DAM · 相关网络\nWGCNA · O2PLS · 联合通路富集",
         colors["teal_face"], colors["teal_edge"]),
        (0.75, 1.60, "04", "代谢物与 MS/MS 注释",
         "分子式与诊断碎片 · MSI 置信等级\nMassBank → MoNA → PubChem → 本地规则",
         colors["purple_face"], colors["purple_edge"]),
        (5.90, 1.60, "05", "转录因子与 Motif 分析",
         "TF 注释 · JASPAR Motif 扫描\nPfam HMM 预测 29 个主要植物 TF 家族",
         colors["orange_face"], colors["orange_edge"]),
        (11.05, 1.60, "06", "机制假设推理",
         "证据融合 · 矛盾检测 · 证伪分析\n竞争假设排序 · 知识缺口 · 验证实验",
         colors["red_face"], colors["red_edge"]),
    ]
    for card in cards:
        capability_card(
            ax, *card[:5], face=card[5], edge=card[6],
            regular=regular, bold=bold,
        )

    # Central platform connects the six abilities.
    box(
        ax, 3.10, 4.00, 9.80, 1.20,
        face=colors["center_face"], edge=colors["center_edge"],
        linewidth=1.9, radius=0.16, zorder=4,
    )
    ax.text(
        8.00, 4.78, "PhytoReason-Agent 科研推理工作台",
        ha="center", va="center", color=colors["ink"],
        fontsize=16.5, fontproperties=bold, fontweight="bold", zorder=5,
    )
    ax.text(
        8.00, 4.35,
        "Layer 0 零数据问答  ·  Layer 1 跨物种推断  ·  Layer 2 单组学分析  ·  Layer 3 配对多组学分析",
        ha="center", va="center", color=colors["muted"],
        fontsize=10.2, fontproperties=regular, zorder=5,
    )

    top_centers = [2.85, 8.00, 13.15]
    bottom_centers = [2.85, 8.00, 13.15]
    for center_x in top_centers:
        arrow(ax, (8.00, 5.24), (center_x, 5.80), color=colors["line"], linewidth=1.25)
    for center_x in bottom_centers:
        arrow(ax, (8.00, 3.96), (center_x, 3.35), color=colors["line"], linewidth=1.25)

    ax.text(
        8.00, 0.70,
        "核心价值：把分散的数据和知识转化为可追踪、可比较、可验证的机制假设",
        ha="center", va="center", color=colors["center_edge"],
        fontsize=11.3, fontproperties=bold, fontweight="bold",
    )

    output = ROOT / "phytoreason_core_capabilities"
    fig.savefig(output.with_suffix(".png"), dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(output.with_suffix(".svg"), facecolor=fig.get_facecolor())
    plt.close(fig)


def design_card(ax, x, number, title, lines, *, face, edge, regular, bold):
    box(ax, x, 2.38, 2.75, 4.20, face=face, edge=edge, linewidth=1.55, radius=0.14)
    box(ax, x + 0.23, 5.72, 0.60, 0.60, face=edge, edge=edge, linewidth=0, radius=0.30, zorder=4)
    ax.text(
        x + 0.53, 6.02, number,
        ha="center", va="center", color="#F7F8F4",
        fontsize=12, fontproperties=bold, fontweight="bold", zorder=5,
    )
    ax.text(
        x + 0.25, 5.28, title,
        ha="left", va="center", color="#25312C",
        fontsize=12.0, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        x + 0.25, 4.03, lines,
        ha="left", va="center", color="#52605A",
        fontsize=9.7, fontproperties=regular, linespacing=1.55,
    )


def render_design(regular, bold, colors):
    fig, ax = plt.subplots(figsize=(16, 9), facecolor=colors["bg"])
    ax.set_facecolor(colors["bg"])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis("off")

    title_block(
        ax,
        "关键设计特点：面向可信科研推理",
        "系统目标不是生成更确定的答案，而是让证据、计算边界和不确定性更透明",
        bold=bold, regular=regular, ink=colors["ink"], muted=colors["muted"],
    )

    design_cards = [
        (0.55, "01", "计算与推理解耦",
         "Workflow 负责\n确定性数值计算\n\nAgent 负责\n规划、调用和解释",
         colors["teal_face"], colors["teal_edge"]),
        (3.60, "02", "全程证据溯源",
         "记录数据库与文献来源\n保存工具参数和执行结果\n\n结论能够回到\n原始证据进行检查",
         colors["blue_face"], colors["blue_edge"]),
        (6.65, "03", "诚实降级机制",
         "Profile → 同源推断 → 缺口\n远程失败返回 partial / warning\n\n不把超时、无结果\n伪装成成功",
         colors["orange_face"], colors["orange_edge"]),
        (9.70, "04", "保守假设推理",
         "比较竞争假设与矛盾证据\n主动进行证伪分析\n\n只输出候选关系\n不宣称已证明因果",
         colors["red_face"], colors["red_edge"]),
        (12.75, "05", "配置化与可复现",
         "物种、代谢物知识配置化\nTool Registry 模块化扩展\n\n统一数据模型、缓存\n测试与可重复执行",
         colors["purple_face"], colors["purple_edge"]),
    ]
    for card in design_cards:
        design_card(
            ax, *card[:4], face=card[4], edge=card[5],
            regular=regular, bold=bold,
        )

    for left_x in (3.34, 6.39, 9.44, 12.49):
        arrow(ax, (left_x, 4.48), (left_x + 0.22, 4.48), color=colors["line"], linewidth=1.3)

    box(
        ax, 2.35, 0.84, 11.30, 0.82,
        face=colors["center_face"], edge=colors["center_edge"],
        linewidth=1.7, radius=0.16,
    )
    ax.text(
        8.00, 1.25,
        "形成结果：可复现  ·  可解释  ·  不越界  ·  可扩展  ·  可验证",
        ha="center", va="center", color=colors["ink"],
        fontsize=13.3, fontproperties=bold, fontweight="bold",
    )

    output = ROOT / "phytoreason_key_design"
    fig.savefig(output.with_suffix(".png"), dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(output.with_suffix(".svg"), facecolor=fig.get_facecolor())
    plt.close(fig)


def render():
    for font_path in (FONT_REGULAR, FONT_BOLD):
        if font_path.exists():
            fontManager.addfont(font_path)
    regular = FontProperties(fname=str(FONT_REGULAR))
    bold = FontProperties(fname=str(FONT_BOLD if FONT_BOLD.exists() else FONT_REGULAR))
    mpl.rcParams["svg.fonttype"] = "none"

    colors = {
        "bg": "#F7F8F4",
        "ink": "#25312C",
        "muted": "#67736D",
        "line": "#77827C",
        "center_face": "#F3EAD8",
        "center_edge": "#8A672F",
        "blue_face": "#E1EBF3",
        "blue_edge": "#4E7899",
        "green_face": "#E5EFE2",
        "green_edge": "#5B8059",
        "teal_face": "#E0ECE8",
        "teal_edge": "#397767",
        "purple_face": "#E9E7F0",
        "purple_edge": "#6B647F",
        "orange_face": "#F5E5CF",
        "orange_edge": "#B86C2E",
        "red_face": "#F2E2DD",
        "red_edge": "#A75C49",
    }
    render_capabilities(regular, bold, colors)
    render_design(regular, bold, colors)


if __name__ == "__main__":
    render()
