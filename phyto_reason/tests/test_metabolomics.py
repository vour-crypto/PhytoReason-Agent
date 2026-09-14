"""
test_metabolomics.py — 代谢物鉴定子包测试（全本地，无网络）。

覆盖:
  1. 分子式引擎：已知化合物质量 → 正确公式召回；非法质量 → 无候选
  2. 加合物质量换算
  3. MGF / CSV 谱图解析
  4. 注释器：合成谱图 → 结构类 + MSI 分级；未知碎片 → 类级缺失（保守）
  5. 物种 profile 锚定
  6. 工具 handler 集成
"""

from __future__ import annotations

from pathlib import Path

import pytest

from phyto_reason.metabolomics.formula_engine import (
    enumerate_formulas,
    neutral_mass_from_adduct,
)
from phyto_reason.metabolomics.spectrum_io import (
    Spectrum,
    load_spectrum,
    parse_fragment_csv,
    parse_mgf_file,
    parse_msp_file,
)
from phyto_reason.metabolomics.annotator import annotate_spectrum


# ── 1. 分子式引擎 ─────────────────────────────────────────

class TestFormulaEngine:
    def test_berberine_neutral(self):
        """berberine [M]+ 336.1230 → 中性 C20H18NO4 (336.1236)。

        注意：该公式带 Lewis/Senior 软违规（离子型），排序靠后——
        检查全部合法候选而非 top-N 切片。
        """
        neutral = neutral_mass_from_adduct(336.1230, "[M]+")
        cands = enumerate_formulas(neutral, ppm=5)
        best = [f for f in cands if not f.violations]
        assert any(f.formula == "C20H18N1O4" for f in best)

    def test_quercetin(self):
        cands = enumerate_formulas(302.0427, ppm=5)
        valid = [f for f in cands if not f.violations]
        assert any(f.formula == "C15H10O7" for f in valid)

    def test_impossible_mass(self):
        """无 CHNOPS 组合能解释的质量 → 无合法候选（保守，不硬凑）。"""
        cands = enumerate_formulas(123.456, ppm=2)
        valid = [f for f in cands if not f.violations]
        assert not valid

    def test_adduct_math(self):
        assert neutral_mass_from_adduct(303.0499, "[M+H]+") == pytest.approx(302.0426, abs=1e-3)
        assert neutral_mass_from_adduct(336.1230, "[M]+") == pytest.approx(336.1236, abs=1e-3)


# ── 2. 谱图解析 ───────────────────────────────────────────

class TestSpectrumIO:
    def test_mgf_parse(self, tmp_path: Path):
        mgf = tmp_path / "test.mgf"
        mgf.write_text(
            "BEGIN IONS\n"
            "TITLE=scan_42\n"
            "PEPMASS=336.1230\n"
            "CHARGE=1+\n"
            "321.0996 80.0\n"
            "306.0761 50.0\n"
            "278.0812 30.0\n"
            "END IONS\n",
            encoding="utf-8",
        )
        spectra = parse_mgf_file(mgf)
        assert len(spectra) == 1
        s = spectra[0]
        assert s.precursor_mz == pytest.approx(336.1230)
        assert s.scan_id == "scan_42"
        assert len(s.fragments) == 3

    def test_csv_parse(self, tmp_path: Path):
        csv = tmp_path / "frag.csv"
        csv.write_text("mz,intensity\n321.1,80\n306.1,50\n", encoding="utf-8")
        spectra = parse_fragment_csv(csv)
        assert len(spectra) == 1
        assert len(spectra[0].fragments) == 2

    def test_msp_parse(self, tmp_path: Path):
        """MSP 谱库格式：Precursor_type → adduct 标注，Ion_mode → 极性。"""
        msp = tmp_path / "lib.msp"
        msp.write_text(
            "Name: Berberine\n"
            "PrecursorMZ: 336.1230\n"
            "Precursor_type: [M]+\n"
            "Ion_mode: P\n"
            "Num Peaks: 3\n"
            "321.0996 80.0\n"
            "306.0761 50.0\n"
            "278.0812 30.0\n"
            "\n"
            "Name: Salvianolic acid B\n"
            "PrecursorMZ: 717.1456\n"
            "Precursor_type: [M-H]-\n"
            "Ion_mode: N\n"
            "Num Peaks: 2\n"
            "519.0928 60.0\n"
            "321.0400 40.0\n",
            encoding="utf-8",
        )
        spectra = parse_msp_file(msp)
        assert len(spectra) == 2
        pos, neg = spectra
        assert pos.scan_id == "Berberine"
        assert pos.precursor_mz == pytest.approx(336.1230)
        assert pos.adduct == "[M]+"
        assert pos.polarity == "positive"
        assert len(pos.fragments) == 3
        assert neg.adduct == "[M-H]-"
        assert neg.polarity == "negative"

    def test_msp_adduct_used_in_annotation(self, tmp_path: Path):
        """MSP 标注的加合物直接用于注释（不重新推断）。"""
        msp = tmp_path / "neg.msp"
        msp.write_text(
            "Name: SAB\nPrecursorMZ: 717.1456\nPrecursor_type: [M-H]-\n"
            "Ion_mode: N\nNum Peaks: 2\n519.0928 60.0\n321.0400 40.0\n",
            encoding="utf-8",
        )
        spectra = load_spectrum(msp)
        r = annotate_spectrum(spectra[0])
        assert r.adduct == "[M-H]-"
        assert "文件标注" in r.note

    def test_load_unsupported(self, tmp_path: Path):
        with pytest.raises(ValueError):
            load_spectrum(tmp_path / "x.pdf")


# ── 3. 注释器 ─────────────────────────────────────────────

_BERBERINE = Spectrum(precursor_mz=336.1230, fragments=[
    (336.1230, 100), (321.0996, 80), (306.0761, 50), (278.0812, 30),
])
_QUERCETIN = Spectrum(precursor_mz=303.0499, fragments=[
    (303.0499, 100), (153.0182, 60), (137.0233, 40), (229.0495, 25),
])


class TestAnnotator:
    def test_berberine_class_and_msi(self):
        r = annotate_spectrum(_BERBERINE, species="黄连", target_metabolite="berberine")
        assert r.candidates, "应有类候选"
        top = r.candidates[0]
        assert top.class_name == "berberine_protoberberine"
        assert top.msi_level == 2
        assert top.score >= 0.5
        assert any(h["name"] == "berberine" for h in r.profile_hits)
        assert any(not f["violations"] for f in r.formula_candidates)

    def test_quercetin_flavonol(self):
        r = annotate_spectrum(_QUERCETIN)
        assert r.candidates[0].class_name == "flavonol"
        assert r.candidates[0].msi_level == 2

    def test_unknown_fragments_conservative(self):
        """随机碎片 → 无类匹配 → 候选为空（MSI 4，不硬猜）。"""
        r = annotate_spectrum(Spectrum(precursor_mz=400.1234, fragments=[
            (250.0, 100), (180.5, 60), (120.3, 40),
        ]))
        assert not r.candidates
        assert "Level 4" in r.note

    def test_negative_mode_salvianolic_acid_b(self):
        """负离子模式：丹酚酸 B [M-H]- 717.1456 → depside_oligomer。"""
        spec = Spectrum(precursor_mz=717.1456, polarity="negative", fragments=[
            (717.1456, 100), (519.0928, 60), (321.0400, 40),
        ])
        r = annotate_spectrum(spec)
        assert r.adduct == "[M-H]-"
        assert r.candidates[0].class_name == "depside_oligomer"

    def test_negative_spectrum_picks_negative_adduct(self):
        """极性为 negative 时不应推断成正离子加合物。"""
        spec = Spectrum(precursor_mz=717.1456, polarity="negative", fragments=[
            (519.0928, 100), (321.0400, 50),
        ])
        r = annotate_spectrum(spec)
        assert r.adduct == "[M-H]-"

    def test_positive_spectrum_picks_positive_adduct(self):
        spec = Spectrum(precursor_mz=336.1230, polarity="positive", fragments=[
            (321.0996, 80), (306.0761, 50),
        ])
        r = annotate_spectrum(spec)
        assert r.adduct in ("[M+H]+", "[M]+")

    def test_doubly_charged_no_false_profile_hit(self):
        """双电荷谱中性质量与 1+ 标志物不同 → 不得误报 profile 命中。

        真实数据案例：m/z 163.1228 [M+2H]2+ 中性 324.23，
        nicotine 中性 162.12 —— 之前按 precursor m/z 直接比较误报命中。
        """
        spec = Spectrum(precursor_mz=163.1228, adduct="[M+2H]2+", polarity="positive",
                        fragments=[(163.1228, 100), (120.08, 60)])
        r = annotate_spectrum(spec, species="烟草")
        assert not any(h["name"] == "nicotine" for h in r.profile_hits)
        # 双电荷中性质量正确
        assert r.neutral_mass == pytest.approx(324.2310, abs=1e-2)

    def test_species_anchor_only_when_given(self):
        r_no = annotate_spectrum(_BERBERINE)
        r_yes = annotate_spectrum(_BERBERINE, species="黄连")
        assert not r_no.profile_hits
        assert r_yes.profile_hits

    def test_no_fragments_raises(self):
        with pytest.raises(ValueError):
            annotate_spectrum(Spectrum(precursor_mz=336.1230, fragments=[]))


# ── 4. 工具 handler 集成 ──────────────────────────────────

class TestHandler:
    def test_inline_annotation(self):
        from phyto_reason.agent.tool_handlers import handle_annotate_ms2_spectrum
        out = handle_annotate_ms2_spectrum({
            "precursor_mz": 336.1230,
            "fragments": [321.0996, 306.0761, 278.0812],
            "species": "黄连",
        })
        assert "berberine_protoberberine" in out
        assert "MSI Level" in out

    def test_missing_inputs(self):
        from phyto_reason.agent.tool_handlers import handle_annotate_ms2_spectrum
        out = handle_annotate_ms2_spectrum({})
        assert "Error" in out

    def test_tool_schema_registered(self):
        from phyto_reason.agent.tool_definitions import TOOL_DEFINITIONS
        names = [t["function"]["name"] for t in TOOL_DEFINITIONS]
        assert "annotate_ms2_spectrum" in names
        assert len(names) == len(set(names))  # 保持唯一
