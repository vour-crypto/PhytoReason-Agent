"""Render the PhytoReason-Agent layered architecture for presentation slides."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_STEM = ROOT / "phytoreason_system_architecture"
FONT_REGULAR = Path(r"C:\Windows\Fonts\msyh.ttc")
FONT_BOLD = Path(r"C:\Windows\Fonts\msyhbd.ttc")


def rounded_box(ax, x, y, w, h, text, *, face, edge, font, size=11,
                weight="normal", linewidth=1.35, radius=0.10, zorder=3):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=linewidth, zorder=zorder,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center", color="#25312C",
        fontsize=size, fontproperties=font, fontweight=weight,
        linespacing=1.38, zorder=zorder + 1,
    )
    return patch


def layer_label(ax, y, index, title, subtitle, *, accent, bold, regular):
    ax.text(
        0.55, y + 0.34, index,
        ha="left", va="center", color=accent,
        fontsize=20, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        1.10, y + 0.48, title,
        ha="left", va="center", color="#25312C",
        fontsize=12.5, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        1.10, y + 0.18, subtitle,
        ha="left", va="center", color="#6A756F",
        fontsize=8.8, fontproperties=regular,
    )


def arrow(ax, start, end, *, color="#737E78", linewidth=1.6,
          linestyle="solid", rad=0.0, mutation=13, zorder=2):
    patch = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=mutation,
        color=color, linewidth=linewidth, linestyle=linestyle,
        connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


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
        "muted": "#6A756F",
        "interaction_face": "#E5ECE7",
        "interaction_edge": "#557162",
        "agent_face": "#F5E4CB",
        "agent_edge": "#B66A2D",
        "knowledge_face": "#E2ECF4",
        "knowledge_edge": "#4F7898",
        "tool_face": "#E8E6EF",
        "tool_edge": "#6A637F",
        "pipeline_face": "#E0ECE8",
        "pipeline_edge": "#397767",
        "reason_face": "#F1E1DC",
        "reason_edge": "#A65C49",
        "foundation_face": "#ECEDE8",
        "foundation_edge": "#777D74",
        "line": "#78827D",
    }

    fig, ax = plt.subplots(figsize=(16, 9), facecolor=colors["bg"])
    ax.set_facecolor(colors["bg"])
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis("off")

    ax.text(
        0.55, 8.55, "PhytoReason-Agent 系统架构",
        ha="left", va="center", color=colors["ink"],
        fontsize=24, fontproperties=bold, fontweight="bold",
    )
    ax.text(
        0.57, 8.14,
        "分层解耦交互、Agent 编排、科学计算、证据推理与数据基础设施",
        ha="left", va="center", color=colors["muted"],
        fontsize=11.5, fontproperties=regular,
    )

    # Layer 1: interaction
    layer_label(
        ax, 7.05, "01", "交互层", "统一研究入口",
        accent=colors["interaction_edge"], bold=bold, regular=regular,
    )
    interaction = [
        (2.35, "Web 科研工作台"),
        (5.65, "PySide6 桌面端"),
        (8.95, "CLI 命令行"),
        (12.25, "FastAPI 服务接口"),
    ]
    for x, text in interaction:
        rounded_box(
            ax, x, 6.88, 2.70, 0.78, text,
            face=colors["interaction_face"], edge=colors["interaction_edge"],
            font=regular, size=11.2,
        )

    # Layer 2: agent orchestration
    layer_label(
        ax, 5.85, "02", "Agent 编排层", "意图、工具与会话协调",
        accent=colors["agent_edge"], bold=bold, regular=regular,
    )
    rounded_box(
        ax, 2.35, 5.50, 4.05, 0.92,
        "AgentOrchestrator\n意图识别 · 任务规划 · 工具选择 · 结果整合",
        face=colors["agent_face"], edge=colors["agent_edge"],
        font=regular, size=10.8, linewidth=1.7,
    )
    rounded_box(
        ax, 6.70, 5.50, 2.45, 0.92,
        "项目与会话上下文\n物种 · 目标物 · 数据状态",
        face=colors["agent_face"], edge=colors["agent_edge"],
        font=regular, size=10.1,
    )
    rounded_box(
        ax, 9.45, 5.50, 2.45, 0.92,
        "LLM 适配层\nOpenAI / 自定义服务",
        face=colors["agent_face"], edge=colors["agent_edge"],
        font=regular, size=10.1,
    )
    rounded_box(
        ax, 12.20, 5.50, 2.75, 0.92,
        "执行追踪与安全控制\n参数校验 · 缓存 · 降级状态",
        face=colors["agent_face"], edge=colors["agent_edge"],
        font=regular, size=10.1,
    )

    # Layer 3: scientific core
    layer_label(
        ax, 3.55, "03", "科学核心层", "知识、计算与推理分工",
        accent=colors["pipeline_edge"], bold=bold, regular=regular,
    )
    rounded_box(
        ax, 2.35, 2.88, 2.85, 2.05,
        "知识层\n\nPubMed · KEGG · PlantCyc\nPlantTFDB · JASPAR\n同源基因 · 质谱数据库",
        face=colors["knowledge_face"], edge=colors["knowledge_edge"],
        font=regular, size=10.3, linewidth=1.55,
    )
    rounded_box(
        ax, 5.45, 2.88, 2.35, 2.05,
        "工具注册与调用\n\nTool Registry\n统一参数模型\n执行结果与来源记录",
        face=colors["tool_face"], edge=colors["tool_edge"],
        font=regular, size=10.3, linewidth=1.55,
    )
    rounded_box(
        ax, 8.05, 2.88, 3.10, 2.05,
        "确定性 Workflow DAG\n\n数据质控 → DEG / DAM\n相关网络 → WGCNA / O2PLS\n联合富集与科学可视化",
        face=colors["pipeline_face"], edge=colors["pipeline_edge"],
        font=regular, size=10.3, linewidth=1.75,
    )
    rounded_box(
        ax, 11.40, 2.88, 3.55, 2.05,
        "Reasoning + Fusion\n\n证据融合 · 矛盾检测\n机制一致性 · 证伪分析\n竞争假设与置信度排序",
        face=colors["reason_face"], edge=colors["reason_edge"],
        font=regular, size=10.3, linewidth=1.75,
    )

    # Layer 4: data foundation
    layer_label(
        ax, 1.55, "04", "数据基础层", "可复现、可追踪、可扩展",
        accent=colors["foundation_edge"], bold=bold, regular=regular,
    )
    foundation = [
        (2.35, 2.85, "数据摄入与校验\n矩阵解析 · 样本对齐 · 长表透视"),
        (5.45, 2.85, "核心数据模型\nEvidence · Candidate · Hypothesis"),
        (8.55, 2.85, "Ontology 与 Profiles\n物种 · 代谢物 · 通路 · TF"),
        (11.65, 3.30, "本地持久化\nWorkbenchDB · 配置 · 缓存 · 图表"),
    ]
    for x, w, text in foundation:
        rounded_box(
            ax, x, 1.28, w, 0.98, text,
            face=colors["foundation_face"], edge=colors["foundation_edge"],
            font=regular, size=9.9,
        )

    # Downstream orchestration flow
    for x in (3.70, 7.00, 10.30, 13.60):
        arrow(ax, (x, 6.84), (x, 6.47), color=colors["interaction_edge"], linewidth=1.25)

    arrow(ax, (4.38, 5.46), (3.78, 4.98), color=colors["knowledge_edge"], rad=0.08)
    arrow(ax, (6.50, 5.46), (6.63, 4.98), color=colors["tool_edge"])
    arrow(ax, (8.05, 5.46), (9.60, 4.98), color=colors["pipeline_edge"], rad=-0.06)
    arrow(ax, (10.65, 5.46), (13.15, 4.98), color=colors["reason_edge"], rad=-0.08)

    # Core collaboration flow
    arrow(ax, (5.24, 3.90), (5.40, 3.90), color=colors["knowledge_edge"], linewidth=1.35)
    arrow(ax, (7.84, 3.90), (8.00, 3.90), color=colors["tool_edge"], linewidth=1.35)
    arrow(ax, (11.19, 3.90), (11.35, 3.90), color=colors["pipeline_edge"], linewidth=1.35)

    # Foundation supports the scientific core.
    for x, target_x, target_color in (
        (3.78, 3.78, colors["knowledge_edge"]),
        (6.88, 6.63, colors["tool_edge"]),
        (9.98, 9.60, colors["pipeline_edge"]),
        (13.30, 13.15, colors["reason_edge"]),
    ):
        arrow(
            ax, (x, 2.31), (target_x, 2.83),
            color=target_color, linewidth=1.05, linestyle="dashed", mutation=11,
        )

    # Architecture principle
    ax.plot([2.35, 2.74], [0.66, 0.66], color=colors["pipeline_edge"], linewidth=5, solid_capstyle="round")
    ax.text(
        2.90, 0.66, "架构原则：确定性数值计算归 Workflow",
        ha="left", va="center", color=colors["muted"], fontsize=10.4,
        fontproperties=regular,
    )
    ax.plot([7.25, 7.64], [0.66, 0.66], color=colors["agent_edge"], linewidth=5, solid_capstyle="round")
    ax.text(
        7.80, 0.66, "不确定性规划与解释归 Agent / Reasoning",
        ha="left", va="center", color=colors["muted"], fontsize=10.4,
        fontproperties=regular,
    )
    ax.text(
        14.95, 0.66, "两者不混用",
        ha="right", va="center", color=colors["reason_edge"], fontsize=10.4,
        fontproperties=bold, fontweight="bold",
    )

    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(OUTPUT_STEM.with_suffix(".svg"), facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    render()
