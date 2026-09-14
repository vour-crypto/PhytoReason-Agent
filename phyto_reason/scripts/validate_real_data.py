"""
validate_real_data.py — 真实组学数据端到端验证（L3 管线）。

用法:
    python phyto_reason/scripts/validate_real_data.py gene-expression.csv meta.csv \
        --species "Zanthoxylum nitidum" --metabolite nitidine --top-genes 3000

数据留在磁盘（脚本内部 pandas 读取），只输出紧凑摘要，不打印原始矩阵。
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd


def main() -> int:
    expr_file, meta_file = sys.argv[1], sys.argv[2]
    args = [a for a in sys.argv[3:]]
    species = ""
    metabolite = ""
    top_genes = 3000
    json_out = ""
    metadata_file = ""
    i = 0
    while i < len(args):
        if args[i] == "--species" and i + 1 < len(args):
            species = args[i + 1]; i += 2
        elif args[i] == "--metabolite" and i + 1 < len(args):
            metabolite = args[i + 1]; i += 2
        elif args[i] == "--top-genes" and i + 1 < len(args):
            top_genes = int(args[i + 1]); i += 2
        elif args[i] == "--json-out" and i + 1 < len(args):
            json_out = args[i + 1]; i += 2
        elif args[i] == "--metadata" and i + 1 < len(args):
            metadata_file = args[i + 1]; i += 2
        else:
            i += 1

    def read_csv_any_encoding(path: str) -> pd.DataFrame:
        """CSV 读取：UTF-8 优先，GBK 回退（代谢物中文名常见 GBK 编码）。"""
        for enc in ("utf-8", "gbk", "gb18030"):
            try:
                return pd.read_csv(path, index_col=0, encoding=enc)
            except UnicodeDecodeError:
                continue
        raise RuntimeError(f"无法解码 {path}（utf-8/gbk/gb18030 均失败）")

    print(f"[1/4] 加载数据… {expr_file} + {meta_file}")
    expr_df = read_csv_any_encoding(expr_file)
    meta_df = read_csv_any_encoding(meta_file)
    print(f"      表达矩阵: {expr_df.shape[0]} 基因 × {expr_df.shape[1]} 样本")

    print(f"[2/4] 方差筛选 top {top_genes} 基因…")
    top_idx = expr_df.var(axis=1).sort_values(ascending=False).head(top_genes).index
    expr_sub = expr_df.loc[top_idx]

    expr_dict = {g: {c: float(v) for c, v in row.items()} for g, row in expr_sub.iterrows()}
    meta_dict = {m: {c: float(v) for c, v in row.items()} for m, row in meta_df.iterrows()}
    print(f"      转入管线: {len(expr_dict)} 基因, {len(meta_dict)} 代谢物")

    print(f"[3/4] 运行 L3 管线（species={species or '?'}, target={metabolite or 'auto'}）…")
    # 注释文件（自动 glob 表达文件同目录的 *annotation* 文件，如
    # Integrated_Function.annotation.xlsx —— 非模式物种 TF 识别的关键输入）
    annotation_path = None
    for p in sorted(Path(expr_file).resolve().parent.glob("*annotation*")):
        if p.suffix.lower() in (".xlsx", ".csv", ".tsv"):
            annotation_path = p
            break
    ann_idx = None
    if annotation_path:
        print(f"      注释文件: {annotation_path.name}")
        from phyto_reason.ingestion.mappers.annotation_mapper import load_annotation
        ann_idx = load_annotation(annotation_path)
        print(f"      注释索引: {ann_idx.n_genes} 基因")
    else:
        print("      注释文件: 未找到（TF 识别将退化为全部基因）")

    sample_metadata = None
    if metadata_file:
        from phyto_reason.ingestion.parsers.metadata_parser import parse_metadata
        parsed_md = parse_metadata(metadata_file)
        sample_metadata = {
            "source_file": Path(metadata_file).name,
            "sample_ids": parsed_md.sample_ids,
            "columns_present": parsed_md.columns_present,
            "columns_missing": parsed_md.columns_missing,
            "samples": [sample.model_dump() for sample in parsed_md.samples],
        }
        print(f"      分组表: {Path(metadata_file).name} "
              f"({len(parsed_md.sample_ids)} 样本, 列 {parsed_md.columns_present})")

    t0 = time.time()
    from phyto_reason.workflows.workflow_runner import WorkflowRunner
    from phyto_reason.workflows.runtime_state import RuntimeState
    from phyto_reason.models.planner_state import PlannerState

    ps = PlannerState(
        species=species or "Zanthoxylum nitidum",
        target_metabolite=metabolite or "",
        sample_count=expr_df.shape[1],
        has_expression=True, has_metabolite=True, has_annotation=ann_idx is not None,
        expression_matrix=expr_dict, metabolite_matrix=meta_dict,
        sample_metadata=sample_metadata,
    )
    initial = RuntimeState(
        planner_state=ps, annotation_index=ann_idx,
        current_node="data_quality_gate",
    )
    runner = WorkflowRunner()
    result = runner.app.invoke(initial, {"recursion_limit": 50})
    state = result if isinstance(result, RuntimeState) else RuntimeState(**result)
    from phyto_reason.workflows.workflow_runner import _merge_inplace_state
    _merge_inplace_state(initial, state)  # 节点原地写入字段补回
    duration = time.time() - t0

    print(f"[4/4] 摘要（耗时 {duration:.1f}s）")
    print(f"  数据质控: {'通过' if state.data_quality_passed else '未通过'} "
          f"({len(getattr(state, 'data_quality_issues', []))} 项问题)")
    nodes = list(getattr(state, "completed_nodes", []) or [])
    print(f"  完成节点 ({len(nodes)}): {', '.join(nodes[:12])}{'…' if len(nodes) > 12 else ''}")
    print(f"  机制分类: {getattr(state, 'mechanism_type', '?')}")
    deg = getattr(state, "deg_report", {}) or {}
    dam = getattr(state, "dam_report", {}) or {}
    print(f"  DEG: {deg.get('n_fdr_significant', 0)}/{deg.get('n_genes_tested', 0)} 显著 | "
          f"DAM: {dam.get('n_fdr_significant', 0)} 显著")
    tf_cands = getattr(state, "tf_candidates", {}) or {}
    print(f"  TF 候选: {len(tf_cands)} 个")
    hyps = getattr(state, "competing_hypotheses", []) or []
    print(f"  竞争假设: {len(hyps)} 条")
    for h in hyps[:3]:
        print(f"    - {h.title[:90]} [{h.mechanism_type.value}] {h.uncertainty_level}")
    wgcna = getattr(state, "wgcna_report", {}) or {}
    if wgcna.get("n_modules", 0):
        print(f"  WGCNA: {wgcna['n_modules']} 模块 (power={wgcna.get('soft_power', '?')})")

    if json_out:
        # 节点详情导出：每个 report 仅保留标量键与列表长度（剔除 plot_base64 等大字段）
        def slim(report: dict) -> dict:
            out = {}
            for key, value in (report or {}).items():
                if isinstance(value, (int, float, str, bool)) or value is None:
                    if isinstance(value, str) and len(value) > 200:
                        continue
                    out[key] = value
                elif isinstance(value, (list, tuple)):
                    out[f"n_{key}"] = len(value)
                elif isinstance(value, dict):
                    out[f"n_{key}"] = len(value)
            return out

        report = {
            "species": species, "metabolite": metabolite,
            "n_samples": expr_df.shape[1], "duration_s": round(duration, 1),
            "data_quality_passed": bool(getattr(state, "data_quality_passed", False)),
            "mechanism_type": str(getattr(state, "mechanism_type", "")),
            "completed_nodes": nodes,
            "deg_report": slim(deg), "dam_report": slim(dam),
            "multiomics_report": slim(getattr(state, "multiomics_report", {}) or {}),
            "correlation_network_report": slim(getattr(state, "correlation_network_report", {}) or {}),
            "joint_enrichment_report": slim(getattr(state, "joint_enrichment_report", {}) or {}),
            "quadrant_plot_report": slim(getattr(state, "quadrant_plot_report", {}) or {}),
            "wgcna_report": slim(wgcna),
            "o2pls_report": slim(getattr(state, "o2pls_report", {}) or {}),
            "metabolite_validation_note": getattr(state, "metabolite_validation_note", ""),
            "n_tf_candidates": len(tf_cands),
            "n_competing_hypotheses": len(hyps),
            "figure_urls": list(getattr(state, "figure_urls", []) or []),
        }
        Path(json_out).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"  节点详情已导出: {json_out}")

    if state.error:
        print(f"  [错误] {state.error[:200]}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
