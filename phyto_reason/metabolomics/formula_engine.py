"""
formula_engine.py — 精确质量 → 候选分子式（CHNOPS + golden rules）。

确定性算法，无外部依赖。输入中性单同位素质量，枚举 CHNOPS 元素组合，
用七条 golden rules（Kind & Fiehn 2007）过滤后按 ppm 误差排序。

保守原则：只做数学上自洽的候选，不做唯一性断言；
多个候选共存时如实列出，由诊断碎片/文献进一步收敛。
"""

from __future__ import annotations

from dataclasses import dataclass

# ── 单同位素质量（Da）─────────────────────────────────────
_MASSES: dict[str, float] = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "S": 31.9720711744,
    "P": 30.97376199842,
}
_PROTON = 1.007276466621      # [M+H]+ 的质子质量
_ELECTRON = 0.000548579909    # 电子质量（[M]+/[M]- 需校正）

# 常见加合物 → (Δm, 带电离子说明)
# neutral = mz - Δm；负离子模式 neutral = mz + Δm（对 [M-H]- 而言 mz 比中性小一个质子）
ADDUCTS: dict[str, float] = {
    "[M+H]+": _PROTON,
    "[M+Na]+": 22.989218 - _ELECTRON,     # Na+ - e-
    "[M+K]+": 38.963158 - _ELECTRON,
    "[M+H-H2O]+": _PROTON + 18.010564684,  # 脱水加质子
    "[M]+": -_ELECTRON,                    # 季铵盐阳离子（如 [M]+ berberine）
    "[M-H]-": -_PROTON,                    # 负离子：neutral = mz + proton
    "[M+Cl]-": 34.96885268 - _ELECTRON,    # Cl− 离子（负模式）
    "[M+CH3COO]-": 59.013851 - _ELECTRON,  # 乙酸根负离子（负模式）
    "[M+2H]2+": 2 * _PROTON,               # 双电荷：neutral = 2*mz - 2*proton
    "[M+2H]": 2 * _PROTON,                 # MSP 常见无电荷后缀写法
}


def neutral_mass_from_adduct(mz: float, adduct: str) -> float:
    """由观测 m/z 与加合物计算中性单同位素质量。

    双电荷（2+）离子：中性 = 2×m/z − 2×质子（离子带两个质子）。
    """
    delta = ADDUCTS.get(adduct)
    if delta is None:
        raise ValueError(f"未知加合物: {adduct}（可选: {sorted(ADDUCTS)}）")
    if adduct.endswith("2+"):
        return 2 * mz - delta
    return mz - delta


@dataclass
class FormulaCandidate:
    formula: str
    mass: float
    ppm_error: float
    rdbe: float
    violations: list[str] = None  # type: ignore[assignment]
    soft_violations: list[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.violations is None:
            self.violations = []
        if self.soft_violations is None:
            self.soft_violations = []


def _rdbe(c: int, h: int, n: int, o: int, s: int, p: int) -> float:
    """不饱和度（Ring + Double Bond Equivalents）。"""
    return c - (h + s) / 2 + n / 2 + 1 - p


def _golden_rules(c: int, h: int, n: int, o: int, s: int, p: int,
                  mass: float, violations: list[str],
                  soft_violations: list[str]) -> None:
    """七条 golden rules 过滤（Kind & Fiehn 2007，校准版）。

    校准（真实化合物验证过，避免误杀）:
      - RDBE 允许半整数：even-electron 离子（[M+H]+ 等）RDBE 可为 n+0.5
        （如 berberine 阳离子 RDBE 12.5）
      - 氮规则、Lewis/Senior 价电子奇偶对离子型（季铵盐 [M]+ 等，berberine
        是教科书级例外）天然失配，均降级为软违规，不硬性排除
    """
    rdbe = _rdbe(c, h, n, o, s, p)
    if rdbe < 0 or abs(2 * rdbe - round(2 * rdbe)) > 1e-6:
        violations.append(f"RDBE 非法: {rdbe:.2f}")
    if rdbe > 40:
        violations.append(f"RDBE 过大: {rdbe:.1f}")

    # 氮规则（even-electron 离子）：质量数奇偶——离子型有已知例外 → 软违规
    if mass > 0 and h > 0:
        if (round(mass) + n) % 2 != 0:
            soft_violations.append("氮规则奇偶（可能为离子型）")

    # 元素比例合理性
    if c > 0 and h > 0 and h / c > 3.5:
        violations.append(f"H/C 过高: {h / c:.2f}")
    if c > 0 and n / c > 3:
        violations.append(f"N/C 过高: {n / c:.2f}")
    if c > 0 and o / c > 4:
        violations.append(f"O/C 过高: {o / c:.2f}")
    if s / max(c, 1) > 2:
        violations.append("S/C 过高")
    if c == 0 and (o > 3 or n > 3 or p > 1):
        violations.append("无碳却多杂原子")

    # Lewis/Senior 规则：价电子数应为偶数——离子型天然为奇 → 软违规
    if (4 * c + h + 3 * n + 2 * o + 2 * s + 3 * p) % 2 != 0:
        soft_violations.append("Lewis/Senior 价电子奇偶（可能为离子型）")

    # 杂原子比例规则
    if p > c + 1:
        violations.append("P 原子过多")
    if s > c + 2:
        violations.append("S 原子过多")
    if n > c + 3:
        violations.append("N 原子过多")


def enumerate_formulas(
    neutral_mass: float,
    ppm: float = 5.0,
    elements: str = "CHNOPS",
    max_c: int = 60,
    max_h: int = 120,
    max_n: int = 10,
    max_o: int = 30,
    max_s: int = 5,
    max_p: int = 3,
    top_k: int = 20,
) -> list[FormulaCandidate]:
    """枚举中性质量的 CHNOPS 候选分子式，golden rules 过滤后按 ppm 排序。

    Args:
        neutral_mass: 中性单同位素质量（Da）。
        ppm: 质量容差（ppm）。
        elements: 允许的元素子集，如 "CHNO"。
        top_k: 返回数量上限。

    Returns:
        按 ppm 误差升序的候选列表（已通过 golden rules 或违规数最少者优先）。
    """
    tol = neutral_mass * ppm / 1e6
    candidates: list[FormulaCandidate] = []
    has_s = "S" in elements
    has_p = "P" in elements
    has_n = "N" in elements
    has_o = "O" in elements

    # 递归/迭代枚举 C, N, O, S, P；H 由剩余质量反推
    for c in range(max_c + 1):
        m_c = c * _MASSES["C"]
        if m_c > neutral_mass + tol:
            break
        for n in range(max_n + 1) if has_n else [0]:
            m_cn = m_c + n * _MASSES["N"]
            if m_cn > neutral_mass + tol:
                break
            for o in range(max_o + 1) if has_o else [0]:
                m_cno = m_cn + o * _MASSES["O"]
                if m_cno > neutral_mass + tol:
                    break
                for s in range(max_s + 1) if has_s else [0]:
                    m_cnos = m_cno + s * _MASSES["S"]
                    if m_cnos > neutral_mass + tol:
                        break
                    for p in range(max_p + 1) if has_p else [0]:
                        m_cnosp = m_cnos + p * _MASSES["P"]
                        if m_cnosp > neutral_mass + tol:
                            continue
                        m_h = neutral_mass - m_cnosp
                        if m_h < 0:
                            continue
                        h = round(m_h / _MASSES["H"])
                        if h < 0 or h > max_h:
                            continue
                        calc_mass = m_cnosp + h * _MASSES["H"]
                        ppm_err = abs(calc_mass - neutral_mass) / neutral_mass * 1e6
                        if ppm_err > ppm:
                            continue
                        if c == 0 and h == 0 and n == 0 and o == 0 and s == 0 and p == 0:
                            continue

                        violations: list[str] = []
                        soft_violations: list[str] = []
                        _golden_rules(c, h, n, o, s, p, neutral_mass,
                                      violations, soft_violations)
                        parts = []
                        if c:
                            parts.append(f"C{c}")
                        if h:
                            parts.append(f"H{h}")
                        if n:
                            parts.append(f"N{n}")
                        if o:
                            parts.append(f"O{o}")
                        if s:
                            parts.append(f"S{s}")
                        if p:
                            parts.append(f"P{p}")
                        formula = "".join(parts)
                        candidates.append(FormulaCandidate(
                            formula=formula,
                            mass=calc_mass,
                            ppm_error=ppm_err,
                            rdbe=_rdbe(c, h, n, o, s, p),
                            violations=violations,
                            soft_violations=soft_violations,
                        ))

    # 排序：硬违规数最少优先，其次软违规数，再 ppm
    candidates.sort(key=lambda x: (len(x.violations), len(x.soft_violations), x.ppm_error))
    return candidates[:top_k]
