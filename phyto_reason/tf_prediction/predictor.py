"""tf_prediction — 无注释文件时的转录因子预测（Pfam HMM + 家族规则）。

方法学与 iTAK（Zheng et al. 2020, NAR）一致：
  转录本/蛋白 FASTA → (六框翻译) → pyhmmer HMM 扫描 → 家族规则判定。
家族规则见 knowledge/tf_family_rules.yaml（accession 经 EBI 逐条核验）。
"""

from __future__ import annotations

from pathlib import Path

_PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"


def _rules_path() -> Path:
    return Path(__file__).resolve().parents[1] / "knowledge" / "tf_family_rules.yaml"


def load_family_rules(path: str | Path | None = None) -> dict[str, dict]:
    import yaml
    source = Path(path) if path else _rules_path()
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    families = payload.get("families", {})
    for name, rule in families.items():
        rule["pfam"] = [str(acc) for acc in rule.get("pfam", [])]
        rule["rule"] = str(rule.get("rule", "any")).lower()
    return families


def _default_pfam_source() -> Path | None:
    """环境变量 PHYTOREASON_PFAM_HMM 可指定本地 Pfam-A.hmm(.gz)。"""
    import os
    env = os.getenv("PHYTOREASON_PFAM_HMM")
    return Path(env) if env else None


def _default_cache_dir() -> Path:
    import os
    env = os.getenv("PHYTOREASON_TF_CACHE")
    return Path(env) if env else Path.home() / ".phytoreason" / "tf_prediction"


def ensure_pfam_hmms(accessions: list[str], cache_dir: str | Path,
                     pfam_source: str | Path | None = None,
                     progress=None) -> Path:
    """从 Pfam-A.hmm（.gz）提取规则所需 HMM 子集，缓存后返回子集文件路径。

    pfam_source 可指向本地 Pfam-A.hmm / .gz；缺省自动从 EBI FTP 下载（一次性）。
    """
    import logging
    import re

    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    subset_path = cache / "pfam_tf_subset.hmm"
    wanted = {f"ACC   {acc}" for acc in accessions}
    wanted_check = set(accessions)

    if subset_path.exists() and subset_path.stat().st_size > 0:
        return subset_path

    import gzip
    import urllib.request

    if pfam_source is None:
        pfam_source = _default_pfam_source()
    if pfam_source:
        source = Path(pfam_source)
    else:
        source = cache / "Pfam-A.hmm.gz"
        if not source.exists():
            logging.getLogger("tf_prediction").info("Downloading Pfam-A.hmm.gz (~520MB, one-time) ...")
            urllib.request.urlretrieve(_PFAM_URL, source)  # noqa: S310 - 固定 EBI 官方地址
    opener = gzip.open if source.suffix == ".gz" else open
    current_lines: list[str] | None = None
    current_acc = ""
    kept: list[str] = []
    with opener(source, "rt", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("HMMER3"):
                if current_acc in wanted_check and current_lines:
                    kept.append("".join(current_lines))
                current_lines = [line]
                current_acc = ""
            elif current_lines is not None:
                if line.startswith("ACC  "):
                    current_acc = line.split()[1].split(".")[0]
                current_lines.append(line)
    if current_acc in wanted_check and current_lines:
        kept.append("".join(current_lines))

    if not kept:
        raise ValueError("Pfam 库中未找到规则所需的任何 HMM；请检查 Pfam-A 版本")
    subset_path.write_text("".join(kept), encoding="ascii")
    if progress:
        progress(len(kept))
    return subset_path


def _as_text(value) -> str:
    """pyhmmer 0.10- 返回 bytes、0.11+ 返回 str，这里统一成 str。"""
    return value.decode() if isinstance(value, bytes) else str(value)


def scan_hmms(proteins: dict[str, str], hmm_path: str | Path,
              evalue: float = 1e-5, cpus: int = 4) -> dict[str, dict[str, tuple[float, float]]]:
    """pyhmmer 扫描：返回 gene -> {pfam_acc: (bitscore, evalue)}。"""
    try:
        import pyhmmer
    except ImportError as exc:  # pragma: no cover - 依赖缺失路径
        raise RuntimeError(
            "TF 预测需要 pyhmmer：pip install pyhmmer"
        ) from exc

    alphabet = pyhmmer.easel.Alphabet.amino()
    records = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
        pyhmmer.easel.TextSequence(name=gene_id.encode(), sequence=protein).digitize(alphabet)
        for gene_id, protein in proteins.items()
    ])

    per_gene: dict[str, dict[str, tuple[float, float]]] = {}
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise ValueError(f"HMM 子集为空: {hmm_path}")
    try:  # 旧版 Pipeline 支持 cpus，新版移除
        pipeline = pyhmmer.plan7.Pipeline(alphabet, E=evalue, cpus=cpus)
    except TypeError:
        pipeline = pyhmmer.plan7.Pipeline(alphabet, E=evalue)

    for hmm in hmms:
        query_acc = _as_text(hmm.accession if hmm.accession else hmm.name).split(".")[0]
        top_hits = pipeline.search_hmm(hmm, records)
        for hit in top_hits:
            if not hit.included or hit.evalue > evalue:
                continue
            best = min(hit.domains, key=lambda d: d.c_evalue) if hit.domains else None
            score = best.score if best else hit.score
            ev = best.c_evalue if best else hit.evalue
            per_gene.setdefault(_as_text(hit.name), {})[query_acc] = (float(score), float(ev))
    return per_gene


def apply_family_rules(per_gene_hits: dict[str, dict[str, tuple[float, float]]],
                       families: dict[str, dict], min_bitscore: float = 20.0) -> list[dict]:
    """按家族规则（any/all）判定；all 优先于 any，同级取比特分最高。"""
    rows: list[dict] = []
    for gene_id, hits in per_gene_hits.items():
        candidates: list[tuple[int, float, dict]] = []
        for family, rule in families.items():
            needed = rule["pfam"]
            present = {acc: hits[acc] for acc in needed if acc in hits}
            if rule["rule"] == "all":
                if len(present) != len(needed):
                    continue
                score = min(v[0] for v in present.values())
                ev = max(v[1] for v in present.values())
                specificity = 2
            else:
                if not present:
                    continue
                best_acc = max(present, key=lambda a: present[a][0])
                score = present[best_acc][0]
                ev = present[best_acc][1]
                specificity = 1
            if score < min_bitscore:
                continue
            candidates.append((specificity, score, {
                "family": family, "score": round(score, 1), "evalue": ev,
                "domains": ",".join(sorted(present)),
            }))
        if not candidates:
            continue
        # all 规则（更特异）优先，同级取比特分最高
        best = max(candidates, key=lambda c: (c[0], c[1]))[2]
        rows.append({"gene_id": gene_id, **best})
    rows.sort(key=lambda r: r["gene_id"])
    return rows


def write_annotation_csv(rows: list[dict], out_path: str | Path) -> Path:
    """输出与 annotation_mapper.load_annotation 兼容的注释 CSV。

    列名 geneid / TF_annotation 均在 load_annotation 的自动识别范围内；
    文件名含 'annotation'，可被 validate_real_data.py 的 glob 自动拾取。
    """
    import csv

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["geneid", "TF_annotation"])
        for row in rows:
            writer.writerow([
                row["gene_id"],
                (f"Predicted TF family {row['family']} "
                 f"(HMM {row['domains']}; bitscore {row['score']}; E={row['evalue']:.1e})"),
            ])
    return out


def predict_from_fasta(fasta_path: str | Path, out_path: str | Path | None = None,
                       cache_dir: str | Path | None = None,
                       pfam_source: str | Path | None = None,
                       min_evalue: float = 1e-5, min_aa: int = 30,
                       cpus: int = 4) -> dict:
    """端到端：FASTA →（六框翻译）→ HMM 扫描 → 家族判定 → 注释 CSV。

    返回摘要 dict（不打印矩阵级细节）。
    """
    from phyto_reason.tf_prediction.sequences import protein_sequences_from_fasta

    cache = Path(cache_dir) if cache_dir else _default_cache_dir()
    proteins, mode = protein_sequences_from_fasta(fasta_path, min_aa=min_aa)
    families = load_family_rules()
    hmm_subset = ensure_pfam_hmms(sorted({acc for rule in families.values() for acc in rule["pfam"]}),
                                  cache, pfam_source=pfam_source)
    per_gene = scan_hmms(proteins, hmm_subset, evalue=min_evalue, cpus=cpus)
    rows = apply_family_rules(per_gene, families)

    if out_path is None:
        source = Path(fasta_path)
        out_path = source.parent / "tf_prediction.annotation.csv"
    csv_path = write_annotation_csv(rows, out_path)
    family_counts: dict[str, int] = {}
    for row in rows:
        family_counts[row["family"]] = family_counts.get(row["family"], 0) + 1
    return {
        "input_fasta": str(fasta_path),
        "mode": mode,
        "n_sequences": len(proteins),
        "n_predicted_tf": len(rows),
        "family_counts": dict(sorted(family_counts.items(), key=lambda kv: -kv[1])),
        "annotation_csv": str(csv_path),
    }
