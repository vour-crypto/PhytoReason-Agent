"""
test_onboarding.py — Auto-onboarding 流程测试（全 mock，无网络）。

覆盖:
  1. taxid 解析（PLANT_TAXIDS 查表 / 未知物种返回 None 不崩溃）
  2. LLM 起草 → YAML 渲染 → SpeciesRegistry 可回读
  3. 引用清洗：检索结果外的 PMID 被丢弃（防编造）
  4. write_profile 拒绝覆盖已存在 profile
  5. LLM 不可用时抛 ValueError
  6. genome_status 非法值校验
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml

from phyto_reason.knowledge.species_onboarding import (
    draft_species_profile,
    write_profile,
    OnboardingDraft,
)
from phyto_reason.knowledge.species_registry import SpeciesRegistry


def _mock_literature(pmids=("12345678", "87654321")):
    """构造可用的 LiteratureResult mock 列表。"""
    lit = []
    for i, pmid in enumerate(pmids):
        r = MagicMock()
        r.pmid = pmid
        r.title = f"Biosynthesis of icariin in Epimedium ({i})"
        r.snippet = "Marker metabolites icariin, epimedin C."
        lit.append(r)
    return lit


def _canned_draft() -> dict:
    """LLM 返回的草稿 dict（含一个编造 PMID 用于测试清洗）。"""
    return {
        "scientific_name": "Epimedium pubescens",
        "marker_metabolites": [
            {"name": "Icariin", "formula": "C33H40O15", "mz": 677.2440,
             "adducts": ["[M+H]+"],
             "fragments": [{"mz": 531.1862, "annotation": "[M+H-Rha]+"}],
             "source": "文献"},
            {"name": "epimedin C", "formula": None, "mz": None,
             "adducts": [], "fragments": [], "source": "推测"},
        ],
        "pathway_prior": {
            "flavonol_glycoside_biosynthesis": {
                "genes": [], "enzymes": ["CHS", "CHI", "FNS", "UGT"],
                "known_regulators": [],
            }
        },
        "refs": ["PMID: 12345678", "PMID: 99999999"],   # 后者是编造的
        "knowledge_gaps": ["无染色体级基因组"],
    }


def _mock_llm(draft: dict):
    client = MagicMock()
    client.is_available.return_value = True
    client.chat_structured.return_value = draft
    return client


# ── 1. taxid 解析 ─────────────────────────────────────────

class TestTaxid:
    def test_resolve_from_table(self):
        # 走 PLANT_TAXIDS 查表（丹参已在表内），不触发 NCBI
        from phyto_reason.knowledge.species_onboarding import _resolve_taxid
        assert _resolve_taxid("Salvia miltiorrhiza", None) == 22663
        # 显式参数优先
        assert _resolve_taxid("Salvia miltiorrhiza", 999) == 999

    def test_unknown_species_no_crash(self):
        with patch(
            "phyto_reason.knowledge.ortholog_finder.OrthologFinder._resolve_taxid",
            return_value=None,
        ):
            draft = draft_species_profile("Herba ignota", llm=_mock_llm(_canned_draft()))
        assert any("taxonomy_id" in w for w in draft.warnings)


# ── 2. 起草 + 渲染 + 回读 ─────────────────────────────────

@patch("phyto_reason.knowledge.species_onboarding.PubMedRetriever.smart_search",
       return_value=MagicMock(
           results=_mock_literature(), total_hits=2, query="x", query_note=""))
class TestDraft:
    def test_yaml_roundtrip(self, _smart, tmp_path: Path):
        draft = draft_species_profile(
            "Epimedium pubescens", common_name="柔毛淫羊藿",
            llm=_mock_llm(_canned_draft()),
        )
        assert draft.slug == "epimedium_pubescens"
        assert "Icariin".lower() in draft.yaml_text.lower()

        # 渲染结果能被 SpeciesRegistry 正常加载
        profile_file = tmp_path / "epimedium_pubescens.yaml"
        profile_file.write_text(draft.yaml_text, encoding="utf-8")
        reg = SpeciesRegistry(profiles_dir=tmp_path)
        p = reg.get("epimedium_pubescens")
        assert p is not None
        assert p.common_name == "柔毛淫羊藿"
        assert any(m.name == "icariin" for m in p.marker_metabolites)

    def test_fake_pmid_sanitized(self, _smart):
        draft = draft_species_profile("Epimedium pubescens", llm=_mock_llm(_canned_draft()))
        assert "PMID: 12345678" in draft.yaml_text
        assert "99999999" not in draft.yaml_text
        assert any("防编造" in w for w in draft.warnings)

    def test_no_literature_warns(self, _smart, tmp_path: Path):
        # 检索结果为空 → 警告 + 自动药典条目
        smart = MagicMock(results=[], total_hits=0, query="x", query_note="")
        with patch("phyto_reason.knowledge.species_onboarding.PubMedRetriever.smart_search",
                   return_value=smart):
            draft = draft_species_profile(
                "Epimedium pubescens", common_name="柔毛淫羊藿",
                llm=_mock_llm(_canned_draft()),
            )
        assert any("无文献支撑" in w for w in draft.warnings)
        assert "中国药典" in draft.yaml_text

    def test_llm_unavailable_raises(self, _smart):
        client = MagicMock()
        client.is_available.return_value = False
        with pytest.raises(ValueError, match="LLM 不可用"):
            draft_species_profile("Epimedium pubescens", llm=client)

    def test_invalid_genome_status(self, _smart):
        with pytest.raises(ValueError, match="genome_status"):
            draft_species_profile("Epimedium pubescens", genome_status="nope",
                                  llm=_mock_llm(_canned_draft()))


# ── 3. write_profile ──────────────────────────────────────

class TestWrite:
    def test_refuse_overwrite(self, tmp_path: Path):
        target = tmp_path / "x.yaml"
        target.write_text("old", encoding="utf-8")
        draft = OnboardingDraft(slug="x", yaml_text="new", warnings=[])
        with pytest.raises(FileExistsError):
            write_profile(draft, profiles_dir=tmp_path)

    def test_force_overwrite(self, tmp_path: Path):
        target = tmp_path / "x.yaml"
        target.write_text("old", encoding="utf-8")
        draft = OnboardingDraft(slug="x", yaml_text="new", warnings=[])
        path = write_profile(draft, profiles_dir=tmp_path, force=True)
        assert path.read_text(encoding="utf-8") == "new"

    def test_write_new(self, tmp_path: Path):
        draft = OnboardingDraft(slug="y", yaml_text="content", warnings=[])
        path = write_profile(draft, profiles_dir=tmp_path)
        assert path.exists()
        assert path.read_text(encoding="utf-8") == "content"


# ── 4. 防御性渲染（真实 LLM 暴露的问题）───────────────────

class TestDefensiveRender:
    def test_pathway_prior_as_list(self, tmp_path: Path):
        """json_object 模式下 LLM 可能把 pathway_prior 输出成列表——不能崩。"""
        from phyto_reason.knowledge.species_onboarding import _render_yaml
        data = {
            "scientific_name": "Epimedium pubescens",
            "marker_metabolites": [{"name": "Icariin", "source": "药典"}],
            "pathway_prior": [
                {"pathway": "Flavonoid biosynthesis", "enzymes": ["CHS"], "known_regulators": []},
                {"name": "Phenylpropanoid", "enzymes": ["PAL"]},
                "garbage",  # 非 dict 条目
            ],
            "refs": [],
        }
        text, dropped = _render_yaml("Epimedium pubescens", "淫羊藿", "unavailable", None, data, [])
        assert dropped >= 1  # garbage 被丢弃
        assert "flavonoid_biosynthesis" in text
        assert "phenylpropanoid" in text

    def test_single_char_refs_silently_dropped(self):
        from phyto_reason.knowledge.species_onboarding import _sanitize_refs
        clean, warnings = _sanitize_refs(["待", "人", "工", "确", "认", "PMID: 12345678"], ["12345678"], "")
        assert clean == ["PMID: 12345678"]
        assert all("单字" not in w for w in warnings)  # 单字静默丢弃
