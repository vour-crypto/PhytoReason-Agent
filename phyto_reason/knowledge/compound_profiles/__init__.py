"""
compound_profiles — 代谢物类 MS/MS 诊断规则知识库加载器。

每个代谢物类一个 YAML（如 alkaloids.yaml / terpenoids.yaml / flavonoids.yaml），
schema 与物种无关：类内条目定义母核、诊断碎片离子、加合物与生源通路。

v5.0 起为 MS/MS 诊断离子规则的唯一规范入口
（旧 tools/chemical_expert.py 及其硬编码路径已移除）。

用法:
    from phyto_reason.knowledge.compound_profiles import (
        load_compound_profile, list_compound_classes)

    rules = load_compound_profile("alkaloids")
    classes = list_compound_classes()
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

_PROFILES_DIR = Path(__file__).resolve().parent


@lru_cache(maxsize=16)
def _load_yaml(class_name: str) -> dict[str, Any]:
    """Load a single compound-class YAML (cached)."""
    path = _PROFILES_DIR / f"{class_name}.yaml"
    if not path.exists():
        raise FileNotFoundError(
            f"compound profile '{class_name}' not found in {_PROFILES_DIR} "
            f"(available: {', '.join(list_compound_classes()) or 'none'})"
        )
    with path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def load_compound_profile(class_name: str) -> dict[str, Any]:
    """Load a compound-class rule set, e.g. ``load_compound_profile("alkaloids")``.

    Returns the full dict from the YAML (class entries + shared neutral_losses).
    Unknown class raises FileNotFoundError with the available list.
    """
    return _load_yaml(class_name)


def list_compound_classes() -> list[str]:
    """List available compound classes (YAML files in this directory)."""
    return sorted(p.stem for p in _PROFILES_DIR.glob("*.yaml"))


def _normalize(name: str) -> str:
    """代谢物名归一化：小写、空格/连字符 → 下划线。"""
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def _all_entries() -> list[tuple[str, dict]]:
    """(class_file, entry) 扁平列表：全部代谢物类的全部类条目。"""
    result: list[tuple[str, dict]] = []
    for class_file in list_compound_classes():
        rules = _load_yaml(class_file)
        for entry in rules.get(class_file, []):
            result.append((class_file, entry))
    return result


def get_class_info(metabolite_or_class: str) -> dict | None:
    """按类名或代谢物名检索一个类条目。

    先精确匹配类名（如 "flavonol"），再按 known_metabolites 匹配
    （如 "quercetin" → flavonol 条目）。
    """
    target = _normalize(metabolite_or_class)
    for class_file, entry in _all_entries():
        if _normalize(entry.get("class", "")) == target:
            return {"class_file": class_file, **entry}
    return find_diagnostic_rules(metabolite_or_class)


def find_diagnostic_rules(metabolite: str) -> dict | None:
    """按代谢物名称检索诊断规则（known_metabolites 双向子串匹配）。

    例如 ``find_diagnostic_rules("quercetin")`` 返回 flavonol 类条目
    （含母核/诊断碎片/加合物），未收录返回 None。
    """
    if not metabolite:
        return None
    target = _normalize(metabolite)
    for class_file, entry in _all_entries():
        for km in entry.get("known_metabolites", []):
            km_norm = _normalize(km)
            if target in km_norm or km_norm in target:
                return {"class_file": class_file, **entry}
    return None
