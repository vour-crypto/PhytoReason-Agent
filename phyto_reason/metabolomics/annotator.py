"""
annotator.py — MS² 谱图注释主流程（未知物鉴定）。

流程（纯 Python 确定性算法）:

    谱图 → 加合物推断（[M+H]+ 优先）→ 中性质量
         → 分子式枚举（CHNOPS + golden rules）
         → 诊断碎片匹配 compound_profiles 全 5 类规则（ppm 容差 + 中性丢失加分）
         → 候选代谢物（类公式/质量对齐 + 物种 profile 标志物锚定）
         → MSI 置信度分级（Level 2/3/4；Level 1 需标准品，如实声明）

保守性原则:
    - 只输出证据支持的结论；碎片不足只给类级（Level 3），不给具体化合物
    - 分子式多候选并存时如实列出，不臆断唯一式
    - 物种 profile 命中时显式标注（来源: profile），未命中不编造
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from phyto_reason.knowledge.compound_profiles import list_compound_classes, load_compound_profile
from phyto_reason.knowledge.species_registry import get_species_registry
from phyto_reason.metabolomics.formula_engine import (
    ADDUCTS,
    enumerate_formulas,
    neutral_mass_from_adduct,
)
from phyto_reason.metabolomics.spectrum_io import Spectrum

# 置信度权重（对应 YAML 的 confidence 字段）
_CONF_WEIGHT = {"high": 1.0, "medium": 0.7, "low": 0.4}


@dataclass
class AnnotationCandidate:
    """单个类/化合物注释候选。"""

    rank: int
    class_name: str
    class_file: str
    msi_level: int                      # 2=碎片+公式, 3=类级, 4=仅质量
    score: float
    matched_fragments: list[dict] = field(default_factory=list)
    neutral_losses: list[dict] = field(default_factory=list)
    name: str | None = None             # 质量对齐命中的具体代谢物
    formula: str | None = None
    mass: float | None = None
    adduct: str = ""
    pathway: str = ""
    rule_confidence: str = ""
    species_anchor: bool = False        # 是否命中物种 profile 标志物


@dataclass
class AnnotationResult:
    """注释结果。"""

    precursor_mz: float
    adduct: str
    neutral_mass: float
    formula_candidates: list[dict] = field(default_factory=list)
    candidates: list[AnnotationCandidate] = field(default_factory=list)
    profile_hits: list[dict] = field(default_factory=list)   # 物种 profile 标志物命中
    note: str = ""


def _match_ppm(obs: float, ref: float, tol_ppm: float) -> float | None:
    """ppm 容差匹配，返回误差（ppm），不匹配返回 None。"""
    if ref <= 0:
        return None
    err = abs(obs - ref) / ref * 1e6
    return err if err <= tol_ppm else None


def _top_fragments(fragments: list[tuple[float, float]], precursor_mz: float,
                   n: int = 20) -> list[float]:
    """按强度取 top-N 子离子（仅保留 < precursor 的真实碎片）。"""
    valid = [(m, i) for m, i in fragments if 0 < m < precursor_mz * 0.98]
    valid.sort(key=lambda x: -x[1])
    return [m for m, _ in valid[:n]]


class SpectrumAnnotator:
    """MS² 谱图注释器（未知物鉴定主入口）。"""

    def __init__(self, tolerance_ppm: float = 10.0, min_matched: int = 1) -> None:
        self.tol_ppm = tolerance_ppm
        self.min_matched = min_matched
        self._classes: list[dict] | None = None

    def _load_classes(self) -> list[dict]:
        if self._classes is None:
            classes: list[dict] = []
            for class_file in list_compound_classes():
                rules = load_compound_profile(class_file)
                for entry in rules.get(class_file, []):
                    entry = dict(entry)
                    entry["class_file"] = class_file
                    entry["neutral_losses"] = rules.get("neutral_losses", [])
                    classes.append(entry)
            self._classes = classes
        return self._classes

    # ── 主流程 ──────────────────────────────────────────

    def annotate(self, spectrum: Spectrum, species: str = "",
                 target_metabolite: str = "", polarity: str = "") -> AnnotationResult:
        """注释一张 MS² 谱图。

        Args:
            spectrum: 解析后的谱图（precursor_mz + fragments）。
            species: 可选物种（profile 标志物锚定 + 通路类优先）。
            target_metabolite: 可选目标代谢物（优先展示其类）。
            polarity: "positive"/"negative"/""（空 = 谱图自带或自动推断）。
        """
        if spectrum.precursor_mz <= 0:
            raise ValueError("谱图缺少 precursor m/z（PEPMASS 或 precursor 列）")
        frags = _top_fragments(spectrum.fragments, spectrum.precursor_mz)
        if not frags:
            raise ValueError("谱图没有有效子离子（MS² 碎片为空）")

        # ── 1. 加合物：文件标注（MSP Precursor_type）优先，否则极性感知推断 ──
        polarity = polarity or spectrum.polarity
        if spectrum.adduct in ADDUCTS:
            adduct = spectrum.adduct
        else:
            adduct = self._infer_adduct(spectrum.precursor_mz, polarity)
        neutral_mass = neutral_mass_from_adduct(spectrum.precursor_mz, adduct)
        formula_cands = enumerate_formulas(neutral_mass, ppm=self.tol_ppm)
        formula_ok = any(not f.violations for f in formula_cands)

        # ── 2. 物种 profile 锚定 ────────────────────────
        profile_hits: list[dict] = []
        registry = None
        if species:
            registry = get_species_registry()
            profile = registry.get(species)
            if profile:
                # 命中规则：
                #   1) 中性质量比较（谱图可为 2+ 离子，precursor 不能直接比 1+ 标志物）
                #   2) 1+ 谱附加：离子 m/z 直接一致也算命中（标志物自身就是该离子，
                #      如 berberine [M]+ 336.1230）；2+ 谱禁用此路径（防误报）
                doubly = adduct.endswith("2+")
                for m in profile.marker_metabolites:
                    marker_adduct = (m.adducts or ["[M+H]+"])[0]
                    try:
                        marker_neutral = neutral_mass_from_adduct(m.mz or 0, marker_adduct)
                    except ValueError:
                        marker_neutral = (m.mz or 0) - ADDUCTS["[M+H]+"]
                    err = None
                    if marker_neutral > 0:
                        err = _match_ppm(neutral_mass, marker_neutral, self.tol_ppm)
                    if err is None and not doubly and m.mz:
                        err = _match_ppm(spectrum.precursor_mz, m.mz, self.tol_ppm)
                    if err is not None:
                        profile_hits.append({
                            "name": m.name, "formula": m.formula,
                            "mz": m.mz, "ppm": round(err, 2), "source": m.source,
                            "fragments_expected": [f.get("mz") for f in m.fragments],
                            "target_match": (
                                target_metabolite.lower() in m.name.lower()
                                if target_metabolite else False
                            ),
                        })

        # ── 3. 诊断碎片类匹配 ───────────────────────────
        candidates: list[AnnotationCandidate] = []
        profile_pathways = set()
        if registry is not None and species:
            profile = registry.get(species)
            if profile:
                profile_pathways = set(profile.pathway_prior.keys())

        for cls in self._load_classes():
            ref_frags = cls.get("fragments", [])
            if not ref_frags:
                continue
            matched: list[dict] = []
            for obs in frags:
                best, best_err = None, None
                for rf in ref_frags:
                    err = _match_ppm(obs, rf.get("mz", 0), self.tol_ppm)
                    if err is not None and (best_err is None or err < best_err):
                        best, best_err = rf, err
                if best is not None:
                    matched.append({
                        "observed_mz": obs,
                        "reference_mz": best["mz"],
                        "ppm_error": round(best_err, 2),
                        "annotation": best.get("annotation", ""),
                    })
            if len(matched) < self.min_matched:
                continue

            # 中性丢失加分：precursor - fragment ≈ neutral loss
            losses: list[dict] = []
            for rf in ref_frags:
                delta = spectrum.precursor_mz - rf.get("mz", 0)
                for nl in cls.get("neutral_losses", []):
                    if _match_ppm(delta, nl.get("delta", 0), self.tol_ppm) is not None:
                        losses.append({"loss": nl.get("annotation", ""),
                                       "delta": nl.get("delta")})
                        break

            total_ref = len(ref_frags)
            conf = _CONF_WEIGHT.get(cls.get("confidence", "medium"), 0.7)
            score = (len(matched) / max(total_ref, 1)) * conf
            if losses:
                score = min(score + 0.15 * len(losses), 1.0)
            if cls.get("mass") and _match_ppm(neutral_mass, cls["mass"], self.tol_ppm) is not None:
                score = min(score + 0.25, 1.0)   # 类母核质量对齐，强加分
            cls_key = cls.get("class", "")
            if cls_key in profile_pathways or cls_key in profile_hits:
                score = min(score + 0.1, 1.0)    # 物种通路先验加分

            # MSI 分级
            if len(matched) >= 2 and formula_ok:
                msi = 2
            else:
                msi = 3

            # 候选代谢物：类公式/质量对齐 → 具体名；否则类级
            # 代表名取 known_metabolites 第一个具体化合物（类公式对应它，
            # 如 pyridine_alkaloid 公式 C10H14N2 = nicotine）；
            # 仅排除类级泛称（anthocyanin/flavonoid/alkaloid/terpenoid）
            _GENERIC = {"anthocyanin", "flavonoid", "alkaloid", "terpenoid"}
            name = formula = mass = None
            if cls.get("formula") and cls.get("mass"):
                err = _match_ppm(neutral_mass, cls["mass"], self.tol_ppm)
                if err is not None:
                    name = next(
                        (km for km in cls.get("known_metabolites", [])
                         if km not in _GENERIC and len(km) > 3),
                        cls.get("class"),
                    )
                    formula, mass = cls["formula"], cls["mass"]
            # 物种 profile 命中优先给出名字
            for ph in profile_hits:
                if ph["name"].lower() in (cls.get("class", "").lower()
                                          or " ".join(cls.get("known_metabolites", [])).lower()):
                    name = ph["name"]
                    formula = ph.get("formula") or formula
                    mass = ph.get("mz") or mass

            candidates.append(AnnotationCandidate(
                rank=0,
                class_name=cls_key,
                class_file=cls.get("class_file", ""),
                msi_level=msi,
                score=round(score, 3),
                matched_fragments=matched[:5],
                neutral_losses=losses[:3],
                name=name,
                formula=formula,
                mass=mass,
                adduct=adduct,
                pathway=cls.get("pathway", ""),
                rule_confidence=cls.get("confidence", "medium"),
                species_anchor=bool(profile_hits),
            ))

        candidates.sort(key=lambda c: -c.score)
        for i, c in enumerate(candidates, 1):
            c.rank = i
        candidates = candidates[:8]

        # ── 4. 备注（保守性说明）──────────────────────────
        note = self._build_note(adduct, neutral_mass, formula_cands, candidates,
                                profile_hits, species,
                                adduct_from_file=bool(spectrum.adduct))

        return AnnotationResult(
            precursor_mz=spectrum.precursor_mz,
            adduct=adduct,
            neutral_mass=neutral_mass,
            formula_candidates=[
                {"formula": f.formula, "ppm": round(f.ppm_error, 2),
                 "mass": round(f.mass, 4), "rdbe": round(f.rdbe, 2),
                 "violations": f.violations}
                for f in formula_cands[:8]
            ],
            candidates=candidates,
            profile_hits=profile_hits,
            note=note,
        )

    # ── 内部 ────────────────────────────────────────────

    def _infer_adduct(self, mz: float, polarity: str = "") -> str:
        """加合物推断：按极性优先试探，取第一个有合法分子式候选的。

        positive: [M+H]+ > [M]+ > [M+Na]+ > [M+K]+ > [M+H-H2O]+
        negative: [M-H]- > [M+Cl]- > [M+CH3COO]-
        auto:     先正后负（保守：不排除其他加合物，在 note 中说明）
        """
        if polarity == "negative":
            order = ("[M-H]-", "[M+Cl]-", "[M+CH3COO]-")
        elif polarity == "positive":
            order = ("[M+H]+", "[M]+", "[M+Na]+", "[M+K]+", "[M+H-H2O]+")
        else:
            order = ("[M+H]+", "[M]+", "[M+Na]+", "[M-H]-", "[M+K]+",
                     "[M+H-H2O]+", "[M+Cl]-")
        for adduct in order:
            try:
                neutral = neutral_mass_from_adduct(mz, adduct)
            except ValueError:
                continue
            if neutral <= 0:
                continue
            cands = enumerate_formulas(neutral, ppm=self.tol_ppm, top_k=5)
            if any(not f.violations for f in cands):
                return adduct
        return "[M+H]+" if polarity != "negative" else "[M-H]-"

    def _build_note(self, adduct, neutral_mass, formula_cands, candidates,
                    profile_hits, species, adduct_from_file: bool = False) -> str:
        parts: list[str] = []
        src = "（文件标注）" if adduct_from_file else ""
        parts.append(f"加合物 {adduct}{src}，中性质量 {neutral_mass:.4f} Da"
                     f"（±{self.tol_ppm} ppm）")
        valid = [f for f in formula_cands if not f.violations]
        if valid:
            parts.append(f"分子式候选 {len(valid)} 个: "
                         + ", ".join(f"{f.formula} ({f.ppm_error:.2f} ppm)" for f in valid[:5])
                         + " —— 非唯一，需诊断碎片/文献收敛")
        else:
            parts.append("无通过 golden rules 的分子式候选（质量可能对应多电荷/未知元素）")
        if not candidates:
            parts.append("诊断碎片未匹配到 compound_profiles 任何类 → MSI Level 4"
                         "（仅质量信息），建议查谱库（MoNA/GNPS）或文献")
        else:
            parts.append(f"类匹配 {len(candidates)} 个候选；"
                         "Level 2 = 诊断碎片+公式支持，Level 3 = 仅类级证据。"
                         "Level 1 确认必须与标准品对照 RT + MS²")
        if profile_hits and species:
            names = ", ".join(f"{h['name']}（{h['ppm']} ppm）" for h in profile_hits)
            parts.append(f"物种 profile 标志物质量命中: {names} —— 优先靶向确认")
        return " ".join(parts)


def annotate_spectrum(spectrum: Spectrum, species: str = "",
                      target_metabolite: str = "",
                      tolerance_ppm: float = 10.0,
                      polarity: str = "") -> AnnotationResult:
    """便捷入口。"""
    return SpectrumAnnotator(tolerance_ppm=tolerance_ppm).annotate(
        spectrum, species=species, target_metabolite=target_metabolite,
        polarity=polarity,
    )
