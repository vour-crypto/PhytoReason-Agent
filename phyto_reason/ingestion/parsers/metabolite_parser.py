r"""
metabolite_parser.py — 代谢物矩阵解析器（纯函数）。

必须正确排除 metadata 列。
只将匹配 ^[A-Za-z]+-[0-9]+$ 的列作为 sample_id。
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from phyto_reason.ingestion.models.parsed_metabolite_matrix import ParsedMetaboliteMatrix

# ── Metadata column detection ────────────────────────────────
# Known metadata column names (lowercased, stripped)
_METADATA_EXACT: set[str] = {
    # Compound info
    "name", "tags", "tag", "id", "compound", "metabolite",
    "formula", "annotation", "annot", "annot.",
    # Database IDs
    "kegg", "hmdb", "pubchem", "chebi", "metlin",
    # Classification
    "class", "subclass", "superclass", "kingdom",
    # MS parameters
    "mz", "m/z", "rt", "rt [min]", "rt (min)", "retention time",
    "deltamass", "deltamass [ppm]", "delta mass", "mass error",
    "reference ion", "reference", "mzcloud", "mzcloud best match",
    "adduct", "charge", "polarity", "mass", "area", "intensity", "height",
    "peak area", "peak_area", "peak height", "peak_height",
    # Internal
    "_meta_class", "_meta_subclass",
}
_METADATA_EXACT_LOWER = {s.lower().strip() for s in _METADATA_EXACT}

# Sample ID pattern: matches F-1, L-2, TR-3, FR-1, S-2, Average-F, etc.
_SAMPLE_RE = re.compile(r"^[A-Za-z]+[0-9]*(-[A-Za-z0-9]+)?$", re.IGNORECASE)

# Metadata heuristic: if column name matches any of these, it's metadata
_METADATA_HEURISTIC = re.compile(
    r"|".join([
        r"\[", r"\]", r"\(", r"\)", r"/",    # brackets, slash
        r"\s{2,}",                              # multiple spaces
        r"^\s*$",                               # empty / whitespace
    ])
)


def _looks_like_metadata(col: str) -> bool:
    """Heuristic: does this column name look like metadata (not a sample ID)?"""
    c = col.strip()
    # Exact match in known metadata set
    if c.lower() in _METADATA_EXACT_LOWER:
        return True
    # Contains brackets, slashes, multiple spaces → metadata
    if _METADATA_HEURISTIC.search(c):
        return True
    # Very long names (>30 chars) → likely metadata
    if len(c) > 30:
        return True
    # Starts with "Average" or "Avg" → it's an aggregated sample column, NOT metadata
    if re.match(r"^(average|avg)[\-_\s]", c, re.IGNORECASE):
        return True
    return False


def _looks_like_sample(col: str) -> bool:
    """Heuristic: does this column name look like a sample/replicate ID?"""
    c = col.strip()
    # Direct SAMPLE_PATTERN match
    if _SAMPLE_RE.match(c):
        return True
    return False


# ── Long-table (m/z-RT-intensity) schema ─────────────────────
# Phase 6.4 Step 2: 标准长表 = sample + m/z + RT + intensity 各一列。
_LONG_SAMPLE_NAMES = ("sample", "sample_id", "sampleid", "sample id", "sample_name")
_LONG_MZ_NAMES = ("mz", "m/z")
_LONG_RT_NAMES = ("rt", "retention time", "retentiontime", "rt [min]", "rt (min)")
_LONG_INTENSITY_NAMES = ("intensity", "area", "peak area", "peak_area", "height", "peak height", "peak_height")


def _match_long_columns(columns: list[str]) -> dict[str, str] | None:
    """识别长表列名；命中返回列映射，未命中返回 None。

    保守判定：只认列名精确匹配的标准长表。若 m/z + RT + intensity 三列
    齐全但既无 sample 列、也没有任何样本样式的列，说明样本归属无法
    确定——直接拒绝，绝不把三列错当三个样本。
    """
    lowered: dict[str, str] = {}
    for col in columns:
        key = str(col).strip().lower()
        lowered.setdefault(key, str(col))

    def find(names: tuple[str, ...]) -> str | None:
        for name in names:
            if name in lowered:
                return lowered[name]
        return None

    sample_col = find(_LONG_SAMPLE_NAMES)
    mz_col = find(_LONG_MZ_NAMES)
    rt_col = find(_LONG_RT_NAMES)
    intensity_col = find(_LONG_INTENSITY_NAMES)
    if sample_col and mz_col and rt_col and intensity_col:
        return {"sample": sample_col, "mz": mz_col, "rt": rt_col, "intensity": intensity_col}
    if mz_col and rt_col and intensity_col:
        sample_like = [c for c in columns
                       if c not in {sample_col, mz_col, rt_col, intensity_col}
                       and _looks_like_sample(str(c))]
        if not sample_like:
            raise ValueError(
                "检测到 m/z + RT + intensity 列，但没有 sample 列，无法确定样本归属。"
                "请提供包含 sample 列的标准长表（sample, m/z, RT, intensity），"
                "或改用宽表（行=代谢物，列=样本）。"
            )
    return None


def _pivot_long_table(df: pd.DataFrame, cols: dict[str, str], warnings: list[str]) -> pd.DataFrame:
    """长表 → 宽矩阵（index=特征 ID `mz_<mz>_rt_<rt>`，列=样本，值=intensity）。"""
    work = df.rename(columns={
        cols["mz"]: "__mz", cols["rt"]: "__rt",
        cols["intensity"]: "__int", cols["sample"]: "__sample",
    })
    work["__mz"] = pd.to_numeric(work["__mz"], errors="coerce")
    work["__rt"] = pd.to_numeric(work["__rt"], errors="coerce")
    work["__int"] = pd.to_numeric(work["__int"], errors="coerce")
    before = len(work)
    work = work.dropna(subset=["__mz", "__int"])
    dropped = before - len(work)
    if dropped:
        warnings.append(f"Long table: dropped {dropped} rows without m/z or intensity")

    def feature_id(mz: float, rt: float) -> str:
        if pd.isna(rt):
            return f"mz_{mz:.4f}"
        return f"mz_{mz:.4f}_rt_{rt:.3f}"

    work = work.assign(__feature=[feature_id(m, r) for m, r in zip(work["__mz"], work["__rt"])])
    work["__sample"] = work["__sample"].astype(str)
    sample_order = list(dict.fromkeys(work["__sample"]))
    pivot = work.pivot_table(index="__feature", columns="__sample", values="__int", aggfunc="first")
    pivot = pivot.reindex(columns=sample_order)
    pivot.index = pivot.index.astype(str)
    return pivot


def parse_metabolite(path: str | Path, encoding: str = "utf-8") -> ParsedMetaboliteMatrix:
    """解析代谢物矩阵文件（v2 — handles Compound Discoverer / mzVault exports)."""
    import logging
    _log = logging.getLogger("ingestion.metabolite")

    path = Path(path)
    ext = path.suffix.lower()
    from phyto_reason.ingestion.encoding import resolve_text_encoding

    provenance = {"source_file": str(path)}
    encoding_warnings: list[str] = []
    if ext not in (".xlsx", ".xls"):
        encoding, detected_provenance, encoding_warnings = resolve_text_encoding(path, encoding)
        provenance.update(detected_provenance)

    # ── Auto-detect delimiter ─────────────────────────────
    try:
        if ext == ".xlsx":
            df_raw = pd.read_excel(path, header=0, index_col=None)
        elif ext == ".tsv":
            df_raw = pd.read_csv(path, sep="\t", index_col=None, encoding=encoding)
        else:
            with open(path, encoding=encoding, errors="replace") as fh:
                first_line = fh.readline()
            sep = "\t" if "\t" in first_line else ","
            df_raw = pd.read_csv(path, sep=sep, index_col=None, encoding=encoding)
    except Exception as e:
        raise ValueError(f"Failed to read metabolite file with encoding={encoding}: {e}") from e

    _log.info("Metabolite file loaded: %d rows × %d cols", df_raw.shape[0], df_raw.shape[1])

    warnings: list[str] = list(encoding_warnings)

    # ── 长表检测（Phase 6.4 Step 2）：标准 m/z-RT-intensity 长表先透视 ──
    long_columns = _match_long_columns([str(c) for c in df_raw.columns])
    long_format = long_columns is not None
    if long_format:
        df = _pivot_long_table(df_raw, long_columns, warnings)
        warnings.append(
            f"Standard long table (m/z-RT-intensity) detected: pivoted to "
            f"{df.shape[0]} features × {df.shape[1]} samples"
        )
        _log.info("Long table pivoted: %d features × %d samples", df.shape[0], df.shape[1])
    else:
        df = df_raw.copy()
        if df.shape[1] > 0:
            df.index = df_raw.iloc[:, 0]
            df = df.iloc[:, 1:]
        df = df.apply(pd.to_numeric, errors="coerce")

        # 转置检测：若行名像样本 ID，转置（长表特征名以 mz_ 开头，不会误触发）
        first_idx = str(df.index[0]) if len(df.index) > 0 else ""
        if bool(re.match(r"^(S\d|GSM|SRR|ERR|sample|replicate)", first_idx, re.IGNORECASE)):
            old_shape = df.shape
            df = df.T
            warnings.append(f"Detected samples-as-rows format: transposed {old_shape} -> {df.shape}")

    original_names = list(df.index.astype(str))

    # ── Column classification: 3-tier pipeline ─────────────
    # Tier 1: KNOWN metadata names → ALWAYS exclude (even if numeric, e.g. m/z, Pubchem)
    # Tier 2: Unknown names → content-based (numeric ratio is ground truth)
    # Tier 3: Pure name heuristic (fallback for degenerate cases)
    all_columns = [str(c) for c in df.columns]
    sample_columns: list[str] = []
    excluded_metadata: list[str] = []
    content_resolved: list[str] = []
    aggregated_columns: list[str] = []

    # Pre-compute numeric ratio for each column (Tier 2 data)
    _numeric_ratio: dict[str, float] = {}
    for col in all_columns:
        if col in df.columns:
            col_data = df[col]
            total = len(col_data)
            n_numeric = int(col_data.notna().sum()) if total > 0 else 0
            _numeric_ratio[col] = n_numeric / max(total, 1)

    for col in all_columns:
        # Tier 1: KNOWN metadata → always exclude
        # These are curated names from Compound Discoverer, mzVault, XCMS, etc.
        # Even if they contain numbers (m/z, RT, PubChem ID), they are NOT samples.
        if re.match(r"^(average|avg)[\-_\s]", col.strip(), re.IGNORECASE):
            excluded_metadata.append(col)
            aggregated_columns.append(col)
            continue
        if _looks_like_metadata(col):
            excluded_metadata.append(col)
            continue

        # Tier 2: Unknown column → content-based classification
        ratio = _numeric_ratio.get(col, 0)
        if ratio > 0.5:
            # Mostly numbers → SAMPLE data
            sample_columns.append(col)
        elif ratio > 0:
            # Has SOME numbers (sparse/mixed) → name heuristic
            if _looks_like_sample(col):
                sample_columns.append(col)
                content_resolved.append(col)
            else:
                excluded_metadata.append(col)
                content_resolved.append(col)
        else:
            # 0% numeric → pure strings (CAS, SMILES, InChI, etc.) → METADATA
            excluded_metadata.append(col)
            if _looks_like_sample(col):
                content_resolved.append(col)
                _log.info("Content override: '%s' looks like sample but 0%% numeric -> metadata", col)

    _log.info("Column classification: %d samples, %d metadata, %d content-resolved",
              len(sample_columns), len(excluded_metadata), len(content_resolved))
    _log.info("Sample columns: %s", sample_columns)
    if excluded_metadata:
        _log.info("Metadata columns (first 15): %s", excluded_metadata[:15])

    # ── Fallback: if no sample columns found, use all non-metadata columns ──
    if not sample_columns:
        sample_columns = [c for c in all_columns if c not in excluded_metadata]
        warnings.append("No sample columns matched — using all non-metadata columns as samples")
    if not sample_columns:
        warnings.append("No sample columns detected — metabolite matrix is empty")

    # 仅保留样本列
    df_samples = df[sample_columns] if sample_columns else pd.DataFrame()

    # ── 构建矩阵 ───────────────────────────────────────
    matrix: dict[str, dict[str, float]] = {}
    for met_id in df_samples.index:
        row = df_samples.loc[met_id]
        matrix[str(met_id)] = {
            str(col): float(v) if pd.notna(v) else float("nan")
            for col, v in row.items()
        }

    n_total = df_samples.shape[0] * df_samples.shape[1] if df_samples.shape[1] > 0 else 0
    n_missing = int(df_samples.isna().sum().sum()) if df_samples.shape[1] > 0 else 0
    missing_rate = n_missing / max(n_total, 1) if n_total > 0 else 0.0

    duplicates = df_samples.index[df_samples.index.duplicated()].tolist()
    if duplicates:
        warnings.append(f"Duplicated metabolite names: {duplicates[:10]}")

    # ── Final validation ──────────────────────────────────
    if not sample_columns:
        # Last resort: take ALL columns as samples
        sample_columns = [c for c in all_columns if c not in excluded_metadata]
        if not sample_columns:
            sample_columns = all_columns
        warnings.append(
            f"⚠️ No sample columns identified by name or content. "
            f"Using all {len(sample_columns)} columns as samples. "
            f"Please verify your file format: rows=metabolites, columns=samples."
        )
    elif excluded_metadata:
        warnings.append(f"Excluded {len(excluded_metadata)} metadata columns")
    if aggregated_columns:
        warnings.append("aggregated column excluded: " + ", ".join(aggregated_columns[:10]))

    # ── Classification summary for user feedback ──────────
    classification = {
        "sample_columns": sample_columns,
        "excluded_metadata": excluded_metadata,
        "content_resolved": content_resolved,  # name vs content had conflict
        "aggregated_columns": aggregated_columns,
    }

    return ParsedMetaboliteMatrix(
        matrix=matrix,
        feature_ids=[str(m) for m in df_samples.index],
        sample_ids=sample_columns,
        original_names=original_names,
        n_features=df_samples.shape[0],
        n_samples=len(sample_columns),
        missing_rate=round(missing_rate, 4),
        normalized=False,
        source_file=str(path),
        warnings=warnings,
        classification=classification,
        provenance=provenance,
    )
