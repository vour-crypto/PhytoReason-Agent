"""
expression_parser.py — 表达矩阵解析器（纯函数）。

支持 csv/tsv/xlsx。
自动检测: 基因列、样本列、重复基因ID、非数值值。
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phyto_reason.ingestion.models.parsed_expression_matrix import ParsedExpressionMatrix


def detect_delimiter(path: Path, encoding: str = "utf-8") -> str:
    """Auto-detect file delimiter by scanning first few lines."""
    ext = path.suffix.lower()
    if ext == ".tsv" or ext == ".txt":
        return "\t"
    if ext == ".xlsx":
        return "xlsx"
    # Read first 5 lines to detect delimiter
    with open(path, encoding=encoding, errors="replace") as fh:
        lines = [fh.readline() for _ in range(5)]
    # Count tabs vs commas across all sampled lines
    tabs = sum(line.count("\t") for line in lines)
    commas = sum(line.count(",") for line in lines)
    if tabs > commas:
        return "\t"
    return ","


def parse_expression(path: str | Path, encoding: str = "utf-8") -> ParsedExpressionMatrix:
    """解析表达矩阵文件。"""
    path = Path(path)
    ext = path.suffix.lower()
    from phyto_reason.ingestion.encoding import resolve_text_encoding

    provenance = {"source_file": str(path)}
    encoding_warnings: list[str] = []
    if ext not in (".xlsx", ".xls"):
        encoding, detected_provenance, encoding_warnings = resolve_text_encoding(path, encoding)
        provenance.update(detected_provenance)
    sep = detect_delimiter(path, encoding=encoding)

    if ext == ".xlsx":
        df = pd.read_excel(path, index_col=0)
    else:
        try:
            df = pd.read_csv(path, sep=sep, index_col=0, encoding=encoding)
        except (UnicodeDecodeError, UnicodeError) as exc:
            raise ValueError(f"Failed to read expression file with encoding={encoding}: {exc}") from exc
        except Exception:
            # If the detected delimiter fails, try the opposite one
            fallback_sep = "," if sep == "\t" else "\t"
            try:
                df = pd.read_csv(path, sep=fallback_sep, index_col=0, encoding=encoding)
            except (UnicodeDecodeError, UnicodeError) as exc:
                raise ValueError(f"Failed to read expression file with encoding={encoding}: {exc}") from exc

    df = df.apply(pd.to_numeric, errors="coerce")
    warnings: list[str] = []

    # Auto-detect orientation: if columns > rows*2, likely rows=features, cols=samples.
    # If rows > columns*2, likely rows=samples, cols=features — transpose.
    # Heuristic: if rows look like sample IDs (S01, GSM, SRR) and columns look like gene IDs, transpose.
    import re
    first_col = str(df.columns[0]) if len(df.columns) > 0 else ""
    first_idx = str(df.index[0]) if len(df.index) > 0 else ""
    row_is_sample = bool(re.match(r'^(S\d|GSM|SRR|ERR|sample|replicate)', first_idx, re.IGNORECASE))
    col_is_gene = bool(re.match(r'^(AT|LOC_|Solyc|GRMZM|Ntab|Os|Zm|TF_|G\d)', first_col))
    has_few_samples = df.shape[0] < 5 and df.shape[1] > df.shape[0] * 2

    if (row_is_sample and col_is_gene) or has_few_samples:
        old_shape = df.shape
        df = df.T
        warnings.append(f"Detected samples-as-rows format: transposed {old_shape} → {df.shape}")

    matrix: dict[str, dict[str, float]] = {}

    feature_ids = [str(g) for g in df.index]
    sample_ids = [str(c) for c in df.columns]

    n_total = df.shape[0] * df.shape[1]
    n_missing = int(df.isna().sum().sum())
    missing_rate = n_missing / max(n_total, 1)

    for gene_id in df.index:
        row = df.loc[gene_id]
        matrix[str(gene_id)] = {
            str(col): float(v) if pd.notna(v) else float("nan")
            for col, v in row.items()
        }

    if df.empty or df.shape[1] < 1:
        warnings.append("Empty or single-column matrix — no sample data found")
        return ParsedExpressionMatrix(
            matrix={}, feature_ids=[], sample_ids=[],
            n_features=0, n_samples=0, missing_rate=1.0,
            source_file=str(path), warnings=encoding_warnings + warnings,
            provenance=provenance,
        )

    min_tpm = 1.0
    min_detected = 2
    low_expr = []
    low_det = []
    for gid in df.index:
        vals = [v for v in matrix[str(gid)].values() if isinstance(v,(int,float)) and not pd.isna(v)]
        if vals and pd.Series(vals).mean() < min_tpm:
            low_expr.append(str(gid))
        nd = sum(1 for v in matrix[str(gid)].values() if isinstance(v,(int,float)) and not pd.isna(v) and v>0)
        if nd < min_detected:
            low_det.append(str(gid))
    if low_expr:
        warnings.append("Low expression: "+str(len(low_expr)))
    if low_det:
        warnings.append("Low detection: "+str(len(low_det)))
    duplicates = df.index[df.index.duplicated()].tolist()
    if duplicates:
        warnings.append(f"Duplicated gene IDs found: {duplicates[:10]}")

    all_nan_cols = df.columns[df.isna().all()].tolist()
    if all_nan_cols:
        warnings.append(f"Columns with all NaN: {all_nan_cols}")

    all_nan_rows = df.index[df.isna().all(axis=1)].tolist()
    if all_nan_rows:
        warnings.append(f"Rows with all NaN (skipped in analysis): {all_nan_rows[:5]}")

    return ParsedExpressionMatrix(
        matrix=matrix,
        feature_ids=feature_ids,
        sample_ids=sample_ids,
        n_features=len(feature_ids),
        n_samples=len(sample_ids),
        missing_rate=round(missing_rate, 4),
        normalization_status="raw",
        data_type="unknown",
        source_file=str(path),
        warnings=encoding_warnings + warnings,
        provenance=provenance,
    )
