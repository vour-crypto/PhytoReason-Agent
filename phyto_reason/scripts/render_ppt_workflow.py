"""Render the PhytoReason-Agent workflow diagram for presentation slides."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_STEM = ROOT / "phytoreason_agent_workflow"
FONT_REGULAR = Path(r"C:\Windows\Fonts\msyh.ttc")
FONT_BOLD = Path(r"C:\Windows\Fonts\msyhbd.ttc")


def add_box(ax, x, y, w, h, text, *, face, edge, font, size=12, weight="normal",
            radius=0.12, linewidth=1.4, color="#24312B", zorder=3):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        color=color,
        fontsize=size,
        fontproperties=font,
        fontweight=weight,
        linespacing=1.35,
        zorder=zorder + 1,
    )
    return patch


def add_arrow(ax, start, end, *, color="#66736D", rad=0.0, linewidth=1.8):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=linewidth,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=2,
        shrinkB=2,
        zorder=2,
    )
    ax.add_patch(arrow)
    return arrow


def add_group_label(ax, x, y, text, *, color, font):
    ax.text(
        x,
        y,
        text,
        ha="left",
        va="bottom",
        color=color,
        fontsize=12,
        fontproperties=font,
        fontweight="bold",
        zorder=5,
    )


def render():
    for font_path in (FONT_REGULAR, FONT_BOLD):
        if font_path.exists():
            fontManager.addfont(font_path)

    regular = FontProperties(fname=str(FONT_REGULAR))
    bold = FontProperties(fname=str(FONT_BOLD if FONT_BOLD.exists() else FONT_REGULAR))
    mpl.rcParams["svg.fonttype"] = "none"

    palette = {
        "background": "#F7F8F4",
        "ink": "#24312B",
        "muted": "#68736E",
        "line": "#6C7771",
        "input_face": "#E3EFE2",
        "input_edge": "#56805C",
        "agent_face": "#F6E5CC",
        "agent_edge": "#B66A2D",
        "knowledge_face": "#E1EBF3",
        "knowledge_edge": "#4D7898",
        "pipeline_face": "#E1ECE8",
        "pipeline_edge": "#3D7768",
        "reason_face": "#F1E2DD",
        "reason_edge": "#A65C49",
        "output_face": "#E8E6EF",
        "output_edge": "#69627E",
    }

    fig, ax = plt.subplots(figsize=(16, 9), facecolor=palette["background"])
    ax.set_facecolor(palette["background"])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis("off")

    ax.text(
        0.5, 8.55,
        "PhytoReason-Agent：从研究输入到可验证机制假设",
        ha="left", va="center", color=palette["ink"],
        fontsize=23, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        0.52, 8.13,
        "确定性统计分析提供可复现结果，Agent 负责规划、证据组织与假设推理",
        ha="left", va="center", color=palette["muted"],
        fontsize=11.5, fontproperties=regular,
    )

    # Research inputs
    add_group_label(ax, 0.45, 7.55, "研究输入", color=palette["input_edge"], font=bold)
    input_boxes = [
        (0.45, 6.30, "自然语言问题"),
        (0.45, 5.05, "转录组\n表达矩阵"),
        (0.45, 3.80, "代谢组\n代谢矩阵"),
        (0.45, 2.55, "MS/MS 谱图"),
        (0.45, 1.30, "基因 / 蛋白序列"),
    ]
    for x, y, label in input_boxes:
        add_box(
            ax, x, y, 2.25, 0.88, label,
            face=palette["input_face"], edge=palette["input_edge"],
            font=regular, size=11.5,
        )

    # Agent planning
    add_group_label(ax, 3.15, 6.75, "Agent 编排", color=palette["agent_edge"], font=bold)
    add_box(
        ax, 3.15, 4.55, 2.15, 1.80,
        "任务理解与规划\n\n识别数据类型\n选择工具与分析路径",
        face=palette["agent_face"], edge=palette["agent_edge"],
        font=regular, size=11.5, linewidth=1.7,
    )

    # Knowledge tools
    add_group_label(ax, 5.85, 7.55, "知识与工具调用", color=palette["knowledge_edge"], font=bold)
    knowledge_boxes = [
        (5.85, 6.35, "PubMed 等文献数据库"),
        (8.10, 6.35, "KEGG · PlantCyc 通路"),
        (5.85, 5.15, "PlantTFDB · JASPAR"),
        (8.10, 5.15, "同源基因 · 化合物 · 质谱库"),
    ]
    for x, y, label in knowledge_boxes:
        add_box(
            ax, x, y, 2.00, 0.84, label,
            face=palette["knowledge_face"], edge=palette["knowledge_edge"],
            font=regular, size=10.4,
        )

    # Deterministic pipeline
    add_group_label(ax, 5.85, 4.35, "确定性分析管线", color=palette["pipeline_edge"], font=bold)
    pipeline_labels = ["质量控制 QC", "DEG / DAM", "相关网络", "WGCNA / O2PLS", "联合富集"]
    pipeline_x = [5.85, 6.73, 7.61, 8.49, 9.37]
    pipeline_widths = [0.73, 0.73, 0.73, 0.73, 0.92]
    pipeline_patches = []
    for x, w, label in zip(pipeline_x, pipeline_widths, pipeline_labels):
        patch = add_box(
            ax, x, 2.25, w, 1.62, label.replace(" ", "\n", 1),
            face=palette["pipeline_face"], edge=palette["pipeline_edge"],
            font=regular, size=9.8, radius=0.09,
        )
        pipeline_patches.append((x, w))
    for (x1, w1), (x2, _w2) in zip(pipeline_patches[:-1], pipeline_patches[1:]):
        add_arrow(ax, (x1 + w1 + 0.02, 3.06), (x2 - 0.02, 3.06), color=palette["pipeline_edge"], linewidth=1.3)

    # Reasoning and output
    add_group_label(ax, 10.75, 7.00, "证据推理", color=palette["reason_edge"], font=bold)
    add_box(
        ax, 10.75, 5.10, 2.15, 1.48,
        "证据融合\n来源追踪 · 矛盾检测",
        face=palette["reason_face"], edge=palette["reason_edge"],
        font=regular, size=11.2, linewidth=1.7,
    )
    add_box(
        ax, 10.75, 2.88, 2.15, 1.48,
        "竞争假设\n构建 · 反驳 · 排序",
        face=palette["reason_face"], edge=palette["reason_edge"],
        font=regular, size=11.2, linewidth=1.7,
    )

    add_group_label(ax, 13.40, 7.55, "科研输出", color=palette["output_edge"], font=bold)
    output_boxes = [
        (13.40, 6.30, "候选调控机制"),
        (13.40, 5.05, "证据来源与置信等级"),
        (13.40, 3.80, "知识缺口"),
        (13.40, 2.55, "推荐验证实验"),
    ]
    for x, y, label in output_boxes:
        add_box(
            ax, x, y, 2.15, 0.88, label,
            face=palette["output_face"], edge=palette["output_edge"],
            font=regular, size=10.8,
        )

    # Main flow arrows
    add_arrow(ax, (2.74, 4.92), (3.10, 5.05), color=palette["line"])
    add_arrow(ax, (5.34, 5.55), (5.80, 6.15), color=palette["knowledge_edge"], rad=-0.12)
    add_arrow(ax, (5.34, 4.95), (5.80, 3.55), color=palette["pipeline_edge"], rad=0.12)
    add_arrow(ax, (10.16, 5.58), (10.70, 5.78), color=palette["knowledge_edge"])
    add_arrow(ax, (10.33, 3.07), (10.70, 5.28), color=palette["pipeline_edge"], rad=-0.18)
    add_arrow(ax, (11.83, 5.05), (11.83, 4.41), color=palette["reason_edge"])
    add_arrow(ax, (12.95, 3.62), (13.35, 4.95), color=palette["output_edge"], rad=-0.18)

    # Category cues
    ax.plot([3.20, 3.55], [0.72, 0.72], color=palette["agent_edge"], linewidth=5, solid_capstyle="round")
    ax.text(3.68, 0.72, "Agent 规划与推理", ha="left", va="center", color=palette["muted"], fontsize=10.2, fontproperties=regular)
    ax.plot([6.25, 6.60], [0.72, 0.72], color=palette["pipeline_edge"], linewidth=5, solid_capstyle="round")
    ax.text(6.73, 0.72, "确定性统计计算", ha="left", va="center", color=palette["muted"], fontsize=10.2, fontproperties=regular)
    ax.plot([9.25, 9.60], [0.72, 0.72], color=palette["knowledge_edge"], linewidth=5, solid_capstyle="round")
    ax.text(9.73, 0.72, "外部知识与数据库", ha="left", va="center", color=palette["muted"], fontsize=10.2, fontproperties=regular)

    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(OUTPUT_STEM.with_suffix(".svg"), facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    render()
