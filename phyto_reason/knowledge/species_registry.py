"""
species_registry.py — 物种配置化注册表（v5.0 核心）。

物种知识从代码解绑为数据：每个物种一个 YAML（knowledge/species_profiles/）。
新增物种 = 新增一个 YAML，无需改代码。

查询统一走三级降级链（非模式物种优先的核心承诺）:

    1. 有本地 profile  → 直接答（置信度正常）
    2. 无 profile      → SpeciesKnowledge 同源迁移（置信度降级 + caveats）
    3. 两者都没有      → 显式知识缺口消息（明说"不知道"，不编造）

降级链把"保守性哲学"从措辞规范扩展到物种维度：
模式物种直接答、非模式物种迁移着答、两者都缺就明确说不知道。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

from phyto_reason.knowledge.species_knowledge import (
    SpeciesKnowledge,
    SpeciesKnowledgeResult,
)

_PROFILES_DIR = Path(__file__).resolve().parent / "species_profiles"

# 物种知识状态（knowledge_for 的降级层级）
PROFILE_LEVEL = "profile"          # 有本地 profile
ORTHOLOG_LEVEL = "ortholog"        # 同源迁移
GAP_LEVEL = "gap"                  # 显式知识缺口


@dataclass
class MarkerMetabolite:
    """物种标志性代谢物（targeted confirmation / L0-L2 的输入锚点）。"""

    name: str
    formula: str | None = None
    mz: float | None = None                 # 单同位素 [M+H]+（或注释指定加合物）
    adducts: list[str] = field(default_factory=list)
    fragments: list[dict] = field(default_factory=list)   # [{mz, formula, annotation}]
    reference_rt: float | None = None
    source: str = "文献"                     # 药典 / 文献 / 推测


@dataclass
class SpeciesProfile:
    """单个物种的知识配置（对应 species_profiles/<slug>.yaml）。"""

    slug: str
    scientific_name: str
    common_name: str = ""
    taxonomy_id: int | None = None
    genome_status: str = "unavailable"      # available | partial | unavailable
    marker_metabolites: list[MarkerMetabolite] = field(default_factory=list)
    pathway_prior: dict[str, dict] = field(default_factory=dict)
    refs: list[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, slug: str, data: dict[str, Any]) -> "SpeciesProfile":
        # 物种元信息在 YAML 的 `species:` 包装键下（见 species_profiles/README.md）
        meta = data.get("species") or data
        return cls(
            slug=slug,
            scientific_name=meta.get("scientific_name", slug),
            common_name=meta.get("common_name", ""),
            taxonomy_id=meta.get("taxonomy_id"),
            genome_status=meta.get("genome_status", "unavailable"),
            marker_metabolites=[
                MarkerMetabolite(**m) for m in data.get("marker_metabolites", [])
            ],
            pathway_prior=data.get("pathway_prior", {}),
            refs=data.get("refs", []),
        )

    def to_summary(self) -> str:
        """供 prompt / 报告使用的紧凑摘要。"""
        mets = ", ".join(m.name for m in self.marker_metabolites) or "无"
        paths = ", ".join(self.pathway_prior.keys()) or "无"
        return (
            f"{self.scientific_name}"
            f"{('（' + self.common_name + '）') if self.common_name else ''}: "
            f"标志代谢物 [{mets}]; 通路先验 [{paths}]; "
            f"基因组 {self.genome_status}"
        )


class SpeciesRegistry:
    """物种 profile 注册表：加载、模糊匹配、降级链查询。"""

    def __init__(self, profiles_dir: Path | None = None) -> None:
        self._dir = profiles_dir or _PROFILES_DIR
        self._profiles: dict[str, SpeciesProfile] = {}
        self._knowledge = SpeciesKnowledge()

    # ── 加载与匹配 ──────────────────────────────────────

    def load_all(self) -> dict[str, SpeciesProfile]:
        """加载 profiles 目录下所有 YAML（slug → profile）。"""
        self._profiles = {}
        for path in sorted(self._dir.glob("*.yaml")):
            with path.open(encoding="utf-8") as fh:
                data = yaml.safe_load(fh) or {}
            slug = path.stem
            profile = SpeciesProfile.from_dict(slug, data)
            profile.slug = slug  # 以文件名 slug 为准
            self._profiles[slug] = profile
        return self._profiles

    def get(self, species: str) -> SpeciesProfile | None:
        """模糊匹配：slug / 学名 / 俗名（子串双向匹配，沿用 SpeciesKnowledge.lookup 风格）。"""
        if not self._profiles:
            self.load_all()
        query = species.lower().strip()
        if not query:
            return None

        for slug, p in self._profiles.items():
            if query in slug.lower() or slug.lower() in query:
                return p
        for p in self._profiles.values():
            sci = p.scientific_name.lower()
            common = p.common_name.lower()
            if query in sci or sci in query or (common and (query in common or common in query)):
                return p
        return None

    def list_species(self) -> list[str]:
        if not self._profiles:
            self.load_all()
        return sorted(self._profiles.keys())

    def has_profile(self, species: str) -> bool:
        return self.get(species) is not None

    # ── 降级链查询 ──────────────────────────────────────

    def knowledge_for(
        self,
        species: str,
        metabolite: str = "",
        tf_family: str = "",
        *,
        use_ortholog_fallback: bool = True,
    ) -> tuple[SpeciesKnowledgeResult, str]:
        """物种感知查询，返回 (结果, 知识层级)。

        层级: "profile"（直接答）→ "ortholog"（同源迁移，降级）→ "gap"（显式缺口）。
        """
        profile = self.get(species)

        if profile is not None:
            return self._from_profile(profile, metabolite=metabolite, tf_family=tf_family), PROFILE_LEVEL

        if use_ortholog_fallback:
            result = self._knowledge.infer_regulation(
                target_species=species, metabolite=metabolite, tf_family=tf_family
            )
            result.note = (
                f"[降级] 该物种无本地 profile（species_profiles/），"
                f"推断基于同源迁移与一般保守性规则。{result.note}"
            )
            return result, ORTHOLOG_LEVEL

        return self._gap_result(species, metabolite=metabolite), GAP_LEVEL

    def describe_gap(self, species: str) -> str:
        """显式知识缺口消息（供 prompt 注入 / 报告输出）。"""
        profile = self.get(species)
        if profile is not None:
            return f"该物种有本地 profile（{profile.to_summary()}）。"
        return (
            f"物种 '{species}' 无本地 profile 且公共知识库无条目。"
            f"推断只能基于同源迁移，置信度降级；"
            f"若需要直接知识，请为该物种添加 species_profiles/ 配置文件。"
        )

    # ── 内部 ────────────────────────────────────────────

    def _from_profile(
        self, profile: SpeciesProfile, *, metabolite: str = "", tf_family: str = ""
    ) -> SpeciesKnowledgeResult:
        """profile 直答：标志代谢物 + 通路先验 + 同源推断补充。"""
        known_mets = [m.name for m in profile.marker_metabolites]
        known_paths = list(profile.pathway_prior.keys())

        # 用 profile 的通路/代谢物信息驱动同源推断（tf_family → metabolite 已知关系）
        result = self._knowledge.infer_regulation(
            target_species=profile.scientific_name,
            metabolite=metabolite or (known_mets[0] if known_mets else ""),
            tf_family=tf_family,
        )
        result.known_metabolites = known_mets
        result.known_pathways = known_paths
        result.note = (
            f"[profile] {profile.to_summary()}. "
            f"{result.note if result.ortholog_inferences else '无额外同源推断。'}"
        )
        return result

    def _gap_result(self, species: str, *, metabolite: str = "") -> SpeciesKnowledgeResult:
        return SpeciesKnowledgeResult(
            query_species=species,
            query_metabolite=metabolite,
            note=(
                f"[缺口] 物种 '{species}' 无本地 profile，且未启用同源迁移。"
                f"无法提供基于知识的推断——该问题需要物种配置或文献检索。"
            ),
        )


# ── 模块单例 ─────────────────────────────────────────────

@lru_cache(maxsize=1)
def get_species_registry() -> SpeciesRegistry:
    """全局 SpeciesRegistry 单例（缓存）。"""
    return SpeciesRegistry()
