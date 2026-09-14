"""
metadata_parser.py — 样本元数据解析器。
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phyto_reason.ingestion.models.sample_metadata_model import (
    SampleMetadata, SampleMetadataTable,
)

RECOGNIZED_COLUMNS = {
    "condition", "treatment", "tissue", "genotype", "replicate", "batch",
}


def parse_metadata(path: str | Path) -> SampleMetadataTable:
    """解析样本元数据文件。"""
    path = Path(path)
    df = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",", index_col=0)

    warnings: list[str] = []

    if df.index.duplicated().any():
        dups = df.index[df.index.duplicated()].tolist()
        warnings.append(f"Duplicated sample IDs in metadata: {dups}")

    samples: list[SampleMetadata] = []
    present_cols: list[str] = []
    for col in df.columns:
        col_lower = col.lower().strip()
        if col_lower in RECOGNIZED_COLUMNS:
            present_cols.append(col_lower)

    for sample_id in df.index:
        row = df.loc[sample_id]
        samples.append(SampleMetadata(
            sample_id=str(sample_id),
            condition=str(row.get("condition", row.get("Condition", ""))),
            treatment=str(row.get("treatment", row.get("Treatment", ""))),
            tissue=str(row.get("tissue", row.get("Tissue", ""))),
            genotype=str(row.get("genotype", row.get("Genotype", ""))),
            replicate=str(row.get("replicate", row.get("Replicate", ""))),
            batch=str(row.get("batch", row.get("Batch", ""))),
        ))

    missing_cols = [c for c in ["condition", "tissue"] if c not in present_cols]
    if missing_cols:
        warnings.append(f"Missing recommended columns: {missing_cols}")

    return SampleMetadataTable(
        samples=samples,
        sample_ids=[str(s) for s in df.index],
        columns_present=present_cols,
        columns_missing=missing_cols,
        source_file=str(path),
        warnings=warnings,
    )
