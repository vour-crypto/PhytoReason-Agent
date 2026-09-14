"""
metabolomics — 代谢物鉴定子包（v5.1, 代谢组方向）。

未知物鉴定流水线（纯 Python 确定性算法，无 R 依赖）:

    输入 MS² 谱图 (MGF/CSV)
      → 加合物推断 → 中性质量
      → 分子式枚举（CHNOPS + 七条 golden rules）
      → 诊断碎片匹配 compound_profiles（31 类规则）
      → 候选代谢物（类内 known_metabolites + 质量对齐）
      → MSI 置信度分级（Level 2/3/4，Level 1 需标准品，如实声明）

模块:
    formula_engine: 精确质量 → 候选分子式
    spectrum_io:    MGF / CSV 谱图解析（mzML 需可选依赖 pyteomics）
    annotator:      SpectrumAnnotator 主流程（物种感知：profile 标志代谢物/通路优先）
"""

from phyto_reason.metabolomics.formula_engine import (
    enumerate_formulas,
    neutral_mass_from_adduct,
    ADDUCTS,
)
from phyto_reason.metabolomics.spectrum_io import (
    parse_mgf_file,
    parse_msp_file,
    parse_fragment_csv,
    load_spectrum,
)
from phyto_reason.metabolomics.annotator import (
    SpectrumAnnotator,
    AnnotationCandidate,
    AnnotationResult,
    annotate_spectrum,
)

__all__ = [
    "enumerate_formulas",
    "neutral_mass_from_adduct",
    "ADDUCTS",
    "parse_mgf_file",
    "parse_msp_file",
    "parse_fragment_csv",
    "load_spectrum",
    "SpectrumAnnotator",
    "AnnotationCandidate",
    "AnnotationResult",
    "annotate_spectrum",
]
