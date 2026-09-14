"""PhytoReason 转录因子预测（Pfam HMM + 家族规则，方法学对齐 iTAK）。"""

from phyto_reason.tf_prediction.predictor import (
    apply_family_rules,
    load_family_rules,
    predict_from_fasta,
    scan_hmms,
    write_annotation_csv,
)

__all__ = [
    "apply_family_rules", "load_family_rules", "predict_from_fasta",
    "scan_hmms", "write_annotation_csv",
]
