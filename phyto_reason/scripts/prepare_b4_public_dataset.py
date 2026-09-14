"""
prepare_b4_public_dataset.py — B4 第二数据集准备：Georgii et al. 2017 拟南芥
转录组-代谢组共分析（ArrayExpress E-MTAB-4867 + MetaboLights MTBLS355）。

下载:
  表达  https://www.ebi.ac.uk/biostudies/files/E-MTAB-4867/data.txt
  代谢  https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS355/m_MTBLS355_GC_MS_v2_maf.tsv

输出（b4_public_data/，不入 git）:
  expr_paired.csv   基因 × 配对样本（AGI locus，normalized 微阵列强度）
  metab_paired.csv  代谢物 × 配对样本（GC-MS 丰度）
  metadata.csv      sample_id, condition/tissue=处理, genotype, batch=区组

用法: python phyto_reason/scripts/prepare_b4_public_dataset.py [--download]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parents[2] / "b4_public_data"
EXPR_URL = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-4867/data.txt"
# MAF 文件只有注释无强度；GC-MS 强度矩阵在 FILES/processedGCMSData.txt
GCMS_URL = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS355/FILES/processedGCMSData.txt"
CANONICAL_RE = re.compile(r"^\d+\.(?:WT|TM|DM)\.[A-Za-z_]+\.[123]$")

# 处理标签压缩：H_HrH/H_LrH → H（湿度并入条件说明，不参与分组）
def collapse_treatment(raw: str) -> str:
    return raw.split("_")[0]


def download(url: str, target: Path) -> None:
    import requests
    response = requests.get(url, timeout=600)
    response.raise_for_status()
    target.write_bytes(response.content)
    print(f"  downloaded {target.name}: {len(response.content) // 1024}KB")


def main() -> int:
    if "--download" in sys.argv:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        download(EXPR_URL, OUT_DIR / "data.txt")
        download(GCMS_URL, OUT_DIR / "processedGCMSData.txt")

    import pandas as pd

    expr = pd.read_csv(OUT_DIR / "data.txt", sep="\t", index_col=0)
    expr.columns = [c.replace(" normalised", "").strip() for c in expr.columns]

    # GC-MS 处理数据：index=代谢物名，columns=样本名（与表达侧命名同源）
    gc = pd.read_csv(OUT_DIR / "processedGCMSData.txt", sep="\t", index_col=0)
    gc.index = [str(name).strip() for name in gc.index]
    gc.columns = [str(name).strip() for name in gc.columns]
    gc = gc.apply(pd.to_numeric, errors="coerce")

    sample_cols = [c for c in gc.columns if CANONICAL_RE.match(c)]
    metab = gc[sample_cols]
    metab = metab.T.groupby(level=0).mean().T  # 同名重复列（重复进样）取均值
    metab = metab.groupby(level=0).mean()      # 同名代谢物（同分异构体行）合并
    print(f"GC-MS canonical sample columns: {len(sample_cols)}")

    common = [s for s in expr.columns if s in metab.columns]
    print(f"paired samples: {len(common)} (expr {expr.shape[1]}, metab {metab.shape[1]})")
    if len(common) < 8:
        print("配对样本不足 8，检查命名规则", file=sys.stderr)
        return 1

    expr_paired = expr[common].dropna(axis=1, how="any")
    # 代谢物保留缺失 ≤20% 的行（管线按 NaN 处理缺失）；全空行剔除
    keep = metab[common].isna().mean(axis=1) <= 0.2
    metab_paired = metab.loc[keep, common]
    print(f"after filtering: expr {expr_paired.shape}, metab {metab_paired.shape}")

    expr_paired.to_csv(OUT_DIR / "expr_paired.csv", encoding="utf-8")
    metab_paired.to_csv(OUT_DIR / "metab_paired.csv", encoding="utf-8")

    rows = []
    for sample in common:
        block, genotype, treatment, rep = sample.split(".", 3)
        treatment = collapse_treatment(treatment)
        rows.append({
            "sample_id": sample,
            "condition": treatment,
            "tissue": treatment,
            "genotype": genotype,
            "batch": block,
            "replicate": rep,
        })
    pd.DataFrame(rows).set_index("sample_id").to_csv(OUT_DIR / "metadata.csv", encoding="utf-8")
    print(f"metadata groups: {pd.DataFrame(rows)['condition'].value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
