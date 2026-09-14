"""从 Pfam-A.hmm(.gz) 提取测试夹具 HMM 子集（PF03106 WRKY + PF00010 bHLH）。

用法: python phyto_reason/scripts/make_tf_test_fixture.py [Pfam-A.hmm.gz 路径]
输出: phyto_reason/tests/fixtures/pfam_test_subset.hmm
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

WANTED = {"PF03106", "PF00010"}


def main() -> int:
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path(__file__).resolve().parents[2] / "b4_public_data" / "Pfam-A.hmm.gz"
    out = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "pfam_test_subset.hmm"
    if not source.exists():
        print(f"Pfam 源文件不存在: {source}", file=sys.stderr)
        return 1
    opener = gzip.open if source.suffix == ".gz" else open
    kept: list[str] = []
    acc, buf = "", []
    with opener(source, "rt", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("HMMER3"):
                if acc in WANTED and buf:
                    kept.append("".join(buf))
                acc, buf = "", [line]
            elif buf:
                if line.startswith("ACC "):
                    acc = line.split()[1].split(".")[0]
                buf.append(line)
    if acc in WANTED and buf:
        kept.append("".join(buf))
    if not kept:
        print("未找到目标 HMM", file=sys.stderr)
        return 1
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(kept), encoding="ascii")
    print(f"fixture: {out} ({out.stat().st_size // 1024}KB, {len(kept)} HMMs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
