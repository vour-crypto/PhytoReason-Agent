from __future__ import annotations

import uuid
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# ── 论文风格（对齐用户 D:/liang、D:/liang1 论文 R 脚本的配色与主题）──
# 主题: theme_classic 风格（无网格、黑色轴线）；字体: Times New Roman（SCI 标准）
# 分类色板: 以论文主图的青绿系为主 + NPG 红蓝对比

PALETTE = {
    "blue": "#1F78B4",    # 显著下调 / 主数据色（论文火山图显著色）
    "teal": "#10B48E",    # 主青绿（论文主图色）
    "orange": "#ed7a68",  # 珊瑚红（强调）
    "red": "#E64B35",     # 显著上调（NPG 红）
    "muted": "#999999",   # 非显著（grey60）
    "grid": "#4D4D4D",    # 阈值参考线（深灰，主题已无网格）
    "dark_green": "#1F3F3A",
    "light_teal": "#70d2bb",
    "sky": "#41C5F9",
    "yellow": "#FECE30",
}

# 分类色板（多组柱图/折线/PCA 分组按序取色）
CATEGORICAL = [
    "#1F3F3A", "#10B48E", "#70d2bb", "#41C5F9",
    "#ed7a68", "#FF2D51", "#FECE30", "#8C4356",
]

# 渐变（气泡/热图标注用）：低=steelblue 高=firebrick（论文 GO 气泡图配色）
GRADIENT_LOW = "#4682B4"
GRADIENT_HIGH = "#B22222"

_DIVERGING_CMAP = "RdBu_r"  # 相关/热图发散色，与 pheatmap 默认一致


def apply_style() -> None:
    """全局应用论文风格（幂等，可在任意时机调用）。"""
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Liberation Serif", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "svg.fonttype": "none",   # SVG 文字保持可编辑（投稿修图友好）
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.grid": False,       # theme_classic：无网格
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "black",
        "axes.linewidth": 0.9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.frameon": False,
    })


def style_axes(ax) -> None:
    """theme_classic 风格：白底、无网格、左侧与底部黑色轴线。"""
    import matplotlib.pyplot as plt

    apply_style()
    ax.set_facecolor("white")
    ax.grid(False)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("black")
        ax.spines[spine].set_linewidth(0.9)
    ax.tick_params(colors="black", labelcolor="black")


def new_figure(title: str, *, width: float = 7.0, height: float = 4.5):
    import matplotlib.pyplot as plt

    apply_style()
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    fig.patch.set_facecolor("white")
    # 论文样式：标题居中加粗
    ax.set_title(title, loc="center", fontsize=14, fontweight="bold", pad=12,
                 fontfamily="serif")
    style_axes(ax)
    return fig, ax


def diverging_cmap() -> str:
    return _DIVERGING_CMAP


def gradient_cmap() -> matplotlib.colors.LinearSegmentedColormap:
    """steelblue→firebrick 渐变（气泡图 -log10(p) 着色，论文 GO 气泡图配色）。"""
    import matplotlib.colors as mcolors

    return mcolors.LinearSegmentedColormap.from_list(
        "paper_gradient", [GRADIENT_LOW, GRADIENT_HIGH])


def save_figure(fig, prefix: str) -> str:
    """保存图表并返回 Web 可访问的 PNG URL。

    同时写出同名 SVG 矢量副本（SCI 投稿用，文字保持可编辑）；
    PNG 用于聊天/结果面板内联展示。SVG 失败不影响 PNG 主产物。
    """
    from phyto_reason.platform_paths import figures_dir
    output_dir = figures_dir()
    name = f"{prefix}_{uuid.uuid4().hex[:8]}.png"
    path = output_dir / name
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
    try:
        fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    except Exception:
        pass
    import matplotlib.pyplot as plt
    plt.close(fig)
    return f"/static/figures/{name}"
