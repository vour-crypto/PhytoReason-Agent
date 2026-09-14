"""
spectrum_io.py — MS² 谱图解析（纯 Python，无外部依赖）。

支持:
    MGF:  通用交换格式（每谱块 BEGIN IONS ... END IONS）
    MSP:  谱库格式（NIST/MoNA；含 Precursor_type 加合物标注 + Ion_mode 极性）
    CSV:  两列/多列碎片表（mz, intensity[, annotation]）
    mzML: 可选（需要 pyteomics，pyproject 的 [project.optional-dependencies].ms2）
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

# mzML 仅在明确请求时导入（可选依赖）
try:  # pragma: no cover
    import pyteomics.mzml  # type: ignore
    _HAS_PYTEOMICS = True
except ImportError:
    _HAS_PYTEOMICS = False


@dataclass
class Spectrum:
    """一张 MS² 谱图。"""

    precursor_mz: float
    fragments: list[tuple[float, float]] = field(default_factory=list)  # (mz, intensity)
    annotations: dict[float, str] = field(default_factory=dict)        # mz → 注释
    precursor_charge: int = 1
    polarity: str = ""          # "positive" / "negative" / ""（未知，自动推断）
    adduct: str = ""            # 文件标注的加合物（如 MSP Precursor_type: "[M+H]+"）
    scan_id: str = ""
    note: str = ""


def parse_mgf_file(path: str | Path) -> list[Spectrum]:
    """解析 MGF 文件，返回全部 MS² 谱图。"""
    spectra: list[Spectrum] = []
    current: Spectrum | None = None
    title = ""

    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.upper() == "BEGIN IONS":
                current = Spectrum(precursor_mz=0.0)
                title = ""
            elif s.upper() == "END IONS":
                if current is not None and current.precursor_mz > 0:
                    if title:
                        current.scan_id = title
                    spectra.append(current)
                current = None
            elif current is not None:
                if s.upper().startswith("PEPMASS"):
                    try:
                        current.precursor_mz = float(s.split("=", 1)[1].split()[0])
                    except (ValueError, IndexError):
                        pass
                elif s.upper().startswith("CHARGE"):
                    charge_str = s.split("=", 1)[1].strip() if "=" in s else ""
                    try:
                        current.precursor_charge = int(charge_str.lstrip("+-"))
                    except (ValueError, IndexError):
                        pass
                    if charge_str.endswith("-"):
                        current.polarity = "negative"
                    elif charge_str.endswith("+"):
                        current.polarity = "positive"
                elif s.upper().startswith("TITLE"):
                    title = s.split("=", 1)[1].strip() if "=" in s else s
                elif s and s[0].isdigit() or (s and s[0] in "+-."):
                    parts = s.split()
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1]) if len(parts) > 1 else 1.0
                        current.fragments.append((mz, intensity))
                    except (ValueError, IndexError):
                        continue

    return spectra


def parse_msp_file(path: str | Path) -> list[Spectrum]:
    """解析 MSP 谱库格式（NIST/MoNA）。

    MSP 额外提供:
      - Precursor_type: 加合物标注（如 "[M+H]+" / "[M-H]-"）→ spectrum.adduct
      - Ion_mode: P（positive）/ N（negative）→ spectrum.polarity
      - Name: 化合物名 → spectrum.scan_id
    """
    spectra: list[Spectrum] = []
    current: Spectrum | None = None
    n_peaks = -1

    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.rstrip("\r\n")
            stripped = s.strip()
            if not stripped:
                if current is not None and current.precursor_mz > 0 and current.fragments:
                    spectra.append(current)
                current = None
                continue

            if current is None:
                current = Spectrum(precursor_mz=0.0)
                n_peaks = -1

            if ":" in s and not stripped[0].isdigit():
                key, _, value = s.partition(":")
                key_l = key.strip().lower()
                value = value.strip()
                if key_l == "name":
                    current.scan_id = value
                elif key_l == "precursormz":
                    try:
                        current.precursor_mz = float(value.split()[0])
                    except (ValueError, IndexError):
                        pass
                elif key_l in ("precursor_type", "precursortype"):
                    current.adduct = value.strip(" '\"")
                elif key_l in ("ion_mode", "ionmode"):
                    if value.upper().startswith("P"):
                        current.polarity = "positive"
                    elif value.upper().startswith("N"):
                        current.polarity = "negative"
                elif key_l == "num peaks":
                    try:
                        n_peaks = int(value)
                    except ValueError:
                        pass
                elif key_l in ("precursorcharge", "charge"):
                    try:
                        current.precursor_charge = int(value.lstrip("+-"))
                    except ValueError:
                        pass
            else:
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        current.fragments.append((mz, intensity))
                    except ValueError:
                        continue

    if current is not None and current.precursor_mz > 0 and current.fragments:
        spectra.append(current)
    return spectra


def parse_fragment_csv(path: str | Path) -> list[Spectrum]:
    """解析 CSV 碎片表。

    支持两种布局:
      1) 单谱：首列 mz（含/不含 precursor 行），列头或首行含 precursor m/z
      2) 多谱：列头 mz,intensity,scan_id / mz,intensity,precursor_mz
    """
    spectra: list[Spectrum] = []
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return spectra

    header = [h.strip().lower() for h in rows[0]]
    data_rows = rows[1:]

    # 多谱布局：mz, intensity, precursor_mz 或 mz,intensity,scan_id
    if len(header) >= 3 and header[1] in ("intensity", "int", "abundance"):
        pol_idx = header.index("polarity") if "polarity" in header else (
            header.index("mode") if "mode" in header else -1)
        by_scan: dict[str, Spectrum] = {}
        for r in data_rows:
            if len(r) < 2:
                continue
            try:
                mz = float(r[0])
                intensity = float(r[1])
            except ValueError:
                continue
            key = r[2] if len(r) > 2 and r[2] else "scan_1"
            if key not in by_scan:
                pm = 0.0
                if len(header) >= 4 and header[3] == "precursor_mz" and len(r) > 3:
                    try:
                        pm = float(r[3])
                    except ValueError:
                        pm = 0.0
                pol = ""
                if pol_idx >= 0 and len(r) > pol_idx:
                    pol = "positive" if r[pol_idx].strip().lower().startswith("pos") \
                        else "negative" if r[pol_idx].strip().lower().startswith("neg") else ""
                by_scan[key] = Spectrum(precursor_mz=pm, scan_id=key, polarity=pol)
            by_scan[key].fragments.append((mz, intensity))
        return list(by_scan.values())

    # 单谱布局：mz[, intensity[, annotation]]
    frags: list[tuple[float, float]] = []
    annotations: dict[float, str] = {}
    precursor_mz = 0.0
    for r in data_rows:
        if not r or not r[0].strip():
            continue
        try:
            mz = float(r[0])
        except ValueError:
            continue
        if mz > 5000:
            continue
        intensity = float(r[1]) if len(r) > 1 and r[1].strip() else 1.0
        if len(r) > 2 and r[2].strip():
            annotations[mz] = r[2].strip()
        if len(r) > 3 and r[3].strip():
            try:
                precursor_mz = float(r[3])
            except ValueError:
                pass
        frags.append((mz, intensity))

    if frags:
        spectra.append(Spectrum(
            precursor_mz=precursor_mz,
            fragments=frags,
            annotations=annotations,
            scan_id="csv",
        ))
    return spectra


def load_spectrum(path: str | Path, polarity: str = "") -> list[Spectrum]:
    """按扩展名自动解析谱图文件（mgf/csv/tsv；mzML 需 pyteomics）。

    Args:
        polarity: 强制指定 "positive"/"negative"；空则按文件名启发式
            （含 neg/negative → negative，含 pos/positive → positive），
            再回落到谱图自带信息（MGF CHARGE / CSV polarity 列）。
    """
    p = Path(path)
    ext = p.suffix.lower()
    if ext == ".mgf":
        spectra = parse_mgf_file(p)
    elif ext == ".msp":
        spectra = parse_msp_file(p)
    elif ext in (".csv", ".tsv", ".txt"):
        spectra = parse_fragment_csv(p)
    elif ext == ".mzml":
        if not _HAS_PYTEOMICS:
            raise RuntimeError(
                "mzML 解析需要可选依赖 pyteomics: pip install pyteomics"
            )
        spectra = _parse_mzml(p)
    else:
        raise ValueError(f"不支持的谱图格式: {ext}（支持 .mgf/.msp/.csv/.tsv/.mzml）")

    # 极性：显式参数 > 文件名启发式 > 谱图自带
    if not polarity:
        name = p.name.lower()
        if "neg" in name:
            polarity = "negative"
        elif "pos" in name:
            polarity = "positive"
    for s in spectra:
        if polarity:
            s.polarity = polarity
    return spectra


def _parse_mzml(path: Path) -> list[Spectrum]:  # pragma: no cover
    """mzML 解析（pyteomics，仅可选依赖存在时使用）。"""
    spectra: list[Spectrum] = []
    with pyteomics.mzml.MzML(str(path)) as reader:
        for spec in reader:
            if spec.get("ms level") != 2:
                continue
            try:
                mz_arr = spec["m/z array"]
                int_arr = spec["intensity array"]
            except KeyError:
                continue
            precursors = spec.get("precursorList", {}).get("precursor", [])
            pm = 0.0
            if precursors:
                try:
                    pm = float(precursors[0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"])
                except (KeyError, IndexError, TypeError):
                    pm = 0.0
            spectra.append(Spectrum(
                precursor_mz=pm,
                fragments=[(float(m), float(i)) for m, i in zip(mz_arr, int_arr)],
                scan_id=str(spec.get("id", "")),
            ))
    return spectra
