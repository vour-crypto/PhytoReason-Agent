"""
test_species_registry.py — SpeciesRegistry 三级降级链测试。

覆盖:
  1. profile YAML 加载与反序列化（MarkerMetabolite / pathway_prior）
  2. 模糊匹配（slug / 学名 / 俗名）
  3. knowledge_for 三级降级链: profile → ortholog → gap
  4. describe_gap 显式缺口消息
  5. compound_profiles 加载器
"""

from __future__ import annotations

import pytest

from phyto_reason.knowledge.species_registry import (
    PROFILE_LEVEL,
    ORTHOLOG_LEVEL,
    GAP_LEVEL,
    SpeciesRegistry,
    get_species_registry,
)
from phyto_reason.knowledge.compound_profiles import (
    load_compound_profile,
    list_compound_classes,
)


@pytest.fixture(scope="module")
def registry() -> SpeciesRegistry:
    return SpeciesRegistry()


# ── 1. 加载与反序列化 ─────────────────────────────────────

class TestLoad:
    def test_load_all_contains_starter_profile(self, registry: SpeciesRegistry):
        profiles = registry.load_all()
        assert "zanthoxylum_nitidum" in profiles

    def test_profile_fields_parsed(self, registry: SpeciesRegistry):
        p = registry.get("zanthoxylum_nitidum")
        assert p is not None
        assert p.scientific_name == "Zanthoxylum nitidum"
        assert p.common_name == "两面针"
        assert p.taxonomy_id == 210110
        assert p.genome_status == "unavailable"
        assert any(m.name == "nitidine" for m in p.marker_metabolites)
        nitidine = next(m for m in p.marker_metabolites if m.name == "nitidine")
        assert nitidine.formula == "C21H18NO4"
        assert nitidine.mz == pytest.approx(348.1230)
        assert nitidine.fragments, "nitidine 应有诊断碎片"
        assert "benzylisoquinoline_alkaloid_biosynthesis" in p.pathway_prior
        # known_regulators 显式为空（不编造）
        assert p.pathway_prior["benzylisoquinoline_alkaloid_biosynthesis"]["known_regulators"] == []


# ── 2. 模糊匹配 ───────────────────────────────────────────

class TestLookup:
    @pytest.mark.parametrize("query", [
        "zanthoxylum_nitidum",        # slug
        "Zanthoxylum nitidum",        # 学名
        "两面针",                      # 俗名
        "zanthoxylum",                # 属名子串
        "nitidum",                    # 种加词
    ])
    def test_fuzzy_match(self, registry: SpeciesRegistry, query: str):
        assert registry.get(query) is not None

    def test_unknown_species(self, registry: SpeciesRegistry):
        assert registry.get("Epimedium pubescens") is None
        assert registry.has_profile("淫羊藿") is False


# ── 3. 降级链 ─────────────────────────────────────────────

class TestDegradeChain:
    def test_level_profile(self, registry: SpeciesRegistry):
        result, level = registry.knowledge_for(
            "Zanthoxylum nitidum", metabolite="nitidine", tf_family="MYB"
        )
        assert level == PROFILE_LEVEL
        assert "nitidine" in result.known_metabolites
        assert "benzylisoquinoline" in result.known_pathways[0]
        assert "[profile]" in result.note

    def test_level_ortholog_fallback(self, registry: SpeciesRegistry):
        result, level = registry.knowledge_for(
            "Epimedium pubescens", metabolite="icariin", tf_family="MYB"
        )
        assert level == ORTHOLOG_LEVEL
        assert "[降级]" in result.note

    def test_level_gap_without_fallback(self, registry: SpeciesRegistry):
        result, level = registry.knowledge_for(
            "Epimedium pubescens", use_ortholog_fallback=False
        )
        assert level == GAP_LEVEL
        assert "[缺口]" in result.note

    def test_describe_gap_known(self, registry: SpeciesRegistry):
        assert "有本地 profile" in registry.describe_gap("两面针")

    def test_describe_gap_unknown(self, registry: SpeciesRegistry):
        msg = registry.describe_gap("Herba epimedii")
        assert "无本地 profile" in msg
        assert "置信度降级" in msg


# ── 3b. 多物种 profile（Phase 1 新增）─────────────────────

class TestPhase1Profiles:
    @pytest.mark.parametrize("slug,common,sci,met", [
        ("zanthoxylum_nitidum", "两面针", "Zanthoxylum nitidum", "nitidine"),
        ("coptis_chinensis", "黄连", "Coptis chinensis", "berberine"),
        ("scutellaria_baicalensis", "黄芩", "Scutellaria baicalensis", "baicalein"),
        ("salvia_miltiorrhiza", "丹参", "Salvia miltiorrhiza", "tanshinone_IIA"),
        ("glycyrrhiza_uralensis", "甘草", "Glycyrrhiza uralensis", "glycyrrhizin"),
        ("nicotiana_tabacum", "烟草", "Nicotiana tabacum", "nicotine"),
    ])
    def test_all_profiles_load(self, registry, slug, common, sci, met):
        p = registry.get(slug)
        assert p is not None
        assert p.common_name == common
        assert p.scientific_name == sci
        assert any(m.name == met for m in p.marker_metabolites)
        assert p.pathway_prior, "每个 profile 必须有通路先验"

    def test_common_name_matching(self, registry):
        for name in ["两面针", "黄连", "黄芩", "丹参", "甘草"]:
            assert registry.get(name) is not None, name

    def test_marker_metabolite_has_source(self, registry):
        for slug in registry.list_species():
            for m in registry.get(slug).marker_metabolites:
                assert m.source in {"药典", "文献", "推测"}

    def test_known_regulators_evidence_backed(self, registry):
        """调控因子诚实性：为空 = 显式未知；非空 = 必须带 PMID 证据（不编造）。"""
        for slug in registry.list_species():
            for pw, info in registry.get(slug).pathway_prior.items():
                for reg in info.get("known_regulators") or []:
                    assert "PMID" in reg, f"{slug}/{pw} 调控因子 {reg} 缺少 PMID 证据"


# ── 3c. 工具层降级链集成（无网络，纯本地）──────────────────

class TestHandlerIntegration:
    def test_search_knowledge_base_with_species(self, registry):
        from phyto_reason.agent.tool_handlers import handle_search_knowledge_base
        out = handle_search_knowledge_base({
            "metabolite": "berberine",
            "species": "黄连",
        })
        assert "物种知识: 黄连" in out
        assert "本地 profile" in out
        assert "berberine" in out

    def test_search_knowledge_base_unknown_species(self):
        from phyto_reason.agent.tool_handlers import handle_search_knowledge_base
        out = handle_search_knowledge_base({
            "metabolite": "berberine",
            "species": "Epimedium pubescens",
        })
        assert "同源迁移" in out or "知识缺口" in out

    def test_cross_species_infer_profile_level(self, registry):
        from phyto_reason.agent.tool_handlers import handle_cross_species_infer
        out = handle_cross_species_infer({
            "target_species": "丹参",
            "metabolite": "tanshinone_IIA",
            "tf_family": "MYB",
        })
        assert "知识层级: 本地知识 profile" in out
        assert "tanshinone" in out

    def test_cross_species_infer_gap_level(self):
        from phyto_reason.agent.tool_handlers import handle_cross_species_infer
        out = handle_cross_species_infer({
            "target_species": "Herba epimedii",
            "metabolite": "icariin",
        })
        assert "知识层级" in out


# ── 4. 全局单例 ───────────────────────────────────────────

def test_global_registry_singleton():
    assert get_species_registry() is get_species_registry()


# ── 5. compound_profiles 加载器 ───────────────────────────

class TestCompoundProfiles:
    def test_list_classes(self):
        classes = list_compound_classes()
        assert "alkaloids" in classes

    def test_load_alkaloids(self):
        rules = load_compound_profile("alkaloids")
        alkaloids = rules["alkaloids"]
        assert any(a["class"] == "berberine_protoberberine" for a in alkaloids)
        assert "neutral_losses" in rules

    def test_unknown_class_raises(self):
        with pytest.raises(FileNotFoundError):
            load_compound_profile("terpenoids_not_yet")
