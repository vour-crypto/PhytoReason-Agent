"""TF 预测子包测试：ORF 翻译、家族规则、端到端（夹具 HMM，离线）。"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pyhmmer")

FIXTURE_HMM = Path(__file__).parent / "fixtures" / "pfam_test_subset.hmm"


# ── 序列处理 ───────────────────────────────────────────────

def test_is_nucleotide_detection():
    from phyto_reason.tf_prediction.sequences import is_nucleotide

    assert is_nucleotide("ATGGCGTAA")
    assert not is_nucleotide("MKVLWAALVT")


def test_longest_orf_translation_forward_and_reverse():
    from phyto_reason.tf_prediction.sequences import longest_orf_protein, reverse_complement

    # 正链：ATG AAA TAA → MK（ATG 起始完整 ORF）
    forward = "ACGTATGAAATAA" + "GCG" * 40 + "TAA"
    assert longest_orf_protein(forward, min_aa=2) == "MK"

    # 反链更长的 ORF：把长 ORF 放到反义链
    long_orf_nt = "ATG" + "AAG" * 20 + "TAA"          # M K20 → 21 aa
    rev_long = reverse_complement(long_orf_nt)
    plus_short = "ATG" + "AAA" + "TAA"                # 仅 MK
    seq = plus_short + "GGG" * 5 + rev_long
    protein = longest_orf_protein(seq, min_aa=5)
    assert protein is not None and protein.startswith("M")
    assert len(protein) == 21


# ── 家族规则 ───────────────────────────────────────────────

def test_apply_family_rules_precedence_and_threshold():
    from phyto_reason.tf_prediction.predictor import apply_family_rules

    families = {
        "B3": {"rule": "any", "pfam": ["PF02362"]},
        "ARF": {"rule": "all", "pfam": ["PF02362", "PF06507"]},
        "WRKY": {"rule": "any", "pfam": ["PF03106"]},
    }
    hits = {
        "gene_arf": {"PF02362": (100.0, 1e-30), "PF06507": (80.0, 1e-20)},
        "gene_b3": {"PF02362": (90.0, 1e-25)},
        "gene_wrky": {"PF03106": (5.0, 1e-2)},   # 低于 min_bitscore
    }
    rows = apply_family_rules(hits, families, min_bitscore=20.0)
    by_gene = {row["gene_id"]: row for row in rows}
    # all 规则（ARF）优先于 any（B3）
    assert by_gene["gene_arf"]["family"] == "ARF"
    assert by_gene["gene_arf"]["domains"] == "PF02362,PF06507"
    assert by_gene["gene_b3"]["family"] == "B3"
    assert "gene_wrky" not in by_gene  # 低于阈值被拒


def test_load_family_rules_from_project_yaml():
    from phyto_reason.tf_prediction.predictor import load_family_rules

    families = load_family_rules()
    assert "WRKY" in families and "MYB" in families and "ARF" in families
    assert families["ARF"]["rule"] == "all"
    assert set(families["ARF"]["pfam"]) == {"PF02362", "PF06507"}
    for name, rule in families.items():
        assert rule["rule"] in ("any", "all")
        assert rule["pfam"], f"{name} 缺少 pfam"


# ── 端到端（夹具 HMM，离线；夹具由 scripts/make_tf_test_fixture.py 生成）──

@pytest.mark.skipif(not FIXTURE_HMM.exists(), reason="需要 pfam_test_subset.hmm 夹具")
def test_end_to_end_with_fixture_hmm(tmp_path):
    from phyto_reason.tf_prediction.predictor import predict_from_fasta
    from phyto_reason.tf_prediction.sequences import read_fasta

    # 夹具 HMM 的共识序列本身应对应 HMM 高分命中
    hmm_text = FIXTURE_HMM.read_text(encoding="ascii")
    consensus = ""
    for line in hmm_text.splitlines():
        if line.startswith("HMMER3"):
            consensus = ""
            continue
        parts = line.split()
        if len(parts) >= 22 and parts[1] and not parts[0].startswith("#"):
            # HMM 主体行：位置 + 插入状态 + 20 个氨基酸概率 + 共识字母
            letters = [ch for ch in parts if len(ch) == 1 and ch.isalpha() and ch != "-"]
            if letters:
                consensus += letters[-1]
    assert consensus, "未能从夹具 HMM 中解析共识序列"

    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(
        ">consensus_tf\n" + consensus + "\n"
        ">decoy\n" + "A" * len(consensus) + "\n",
        encoding="utf-8",
    )
    summary = predict_from_fasta(
        fasta, out_path=tmp_path / "tf_prediction.annotation.csv",
        cache_dir=tmp_path / "cache", pfam_source=FIXTURE_HMM, cpus=2,
    )
    assert summary["mode"] == "aa"
    assert summary["n_sequences"] == 2
    assert Path(summary["annotation_csv"]).exists()

    proteins = read_fasta(summary["annotation_csv"])
    # 共识序列应被预测为某 TF 家族；poly-A 诱饵不应命中
    annotated = {row.split(",")[0]: row.split(",")[1] for row in
                 Path(summary["annotation_csv"]).read_text(encoding="utf-8").splitlines()[1:]}
    assert "consensus_tf" in annotated
    assert "decoy" not in annotated
    assert annotated["consensus_tf"].startswith("Predicted TF family")


def test_paper_plot_style_applied():
    """绘图风格回归：Times New Roman 主题 + 论文配色（对齐用户 R 脚本）。"""
    import matplotlib
    matplotlib.use("Agg")
    from phyto_reason.visualization._common import CATEGORICAL, PALETTE, apply_style, gradient_cmap, new_figure

    apply_style()
    assert "serif" in matplotlib.rcParams["font.family"]
    assert "Times New Roman" in matplotlib.rcParams["font.serif"]
    assert matplotlib.rcParams["svg.fonttype"] == "none"   # SVG 文字可编辑
    assert matplotlib.rcParams["axes.grid"] is False       # theme_classic
    assert matplotlib.rcParams["axes.spines.top"] is False
    assert PALETTE["blue"] == "#1F78B4"                    # 论文火山显著蓝
    assert PALETTE["red"] == "#E64B35"                     # NPG 红
    assert PALETTE["teal"] == "#10B48E"
    assert len(CATEGORICAL) >= 8

    fig, ax = new_figure("style probe")
    assert ax.get_title(loc="center") == "style probe"
    cmap = gradient_cmap()
    assert cmap(0.0)[:3] != cmap(1.0)[:3]
    matplotlib.pyplot.close(fig)
