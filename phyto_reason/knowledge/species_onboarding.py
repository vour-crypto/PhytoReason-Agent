"""
species_onboarding.py — Auto-onboarding：LLM 辅助生成物种 profile。

论文方法学卖点："零人工配置接入新物种"。流程:

    1. 物种名 → taxid（ortholog_finder.PLANT_TAXIDS 查表 → NCBI 兜底）
    2. PubMed 检索该物种次生代谢文献（pubmed_retriever，失败不阻塞）
    3. LLM 起草（chat_structured + JSON schema）：marker 代谢物（含 formula/mz/碎片/来源）、
       通路先验（酶）、参考文献（PMID）
    4. SafetyGuard 校验 + 引用清洗（只保留检索结果中真实出现的 PMID，防编造）
    5. 渲染为 species_profiles/<slug>.yaml，--write 落盘（拒绝覆盖已有 profile）

用法:
    from phyto_reason.knowledge.species_onboarding import draft_species_profile, write_profile

    draft = draft_species_profile("Epimedium pubescens", common_name="柔毛淫羊藿")
    print(draft.yaml_text)
    path = write_profile(draft)          # 落盘（已存在则拒绝）
    path = write_profile(draft, force=True)

CLI:
    plantomics onboarding "Epimedium pubescens" --common 柔毛淫羊藿 --write
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from phyto_reason.knowledge.ortholog_finder import PLANT_TAXIDS
from phyto_reason.knowledge.pubmed_retriever import PubMedRetriever
from phyto_reason.llm.llm_client import LLMClient
from phyto_reason.llm.safety_guard import SafetyGuard

_PROFILES_DIR = Path(__file__).resolve().parent / "species_profiles"

# LLM 结构化输出的 JSON schema（DeepSeek json_object 模式）
DRAFT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "scientific_name": {"type": "string"},
        "marker_metabolites": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "formula": {"type": ["string", "null"]},
                    "mz": {"type": ["number", "null"]},
                    "adducts": {"type": "array", "items": {"type": "string"}},
                    "fragments": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "mz": {"type": "number"},
                                "annotation": {"type": "string"},
                            },
                        },
                    },
                    "source": {"type": "string", "enum": ["药典", "文献", "推测"]},
                },
                "required": ["name", "source"],
            },
        },
        "pathway_prior": {
            "type": "object",
            "additionalProperties": {
                "type": "object",
                "properties": {
                    "genes": {"type": "array", "items": {"type": "string"}},
                    "enzymes": {"type": "array", "items": {"type": "string"}},
                    "known_regulators": {"type": "array", "items": {"type": "string"}},
                },
            },
        },
        "refs": {"type": "array", "items": {"type": "string"}},
        "knowledge_gaps": {"type": "array", "items": {"type": "string"}},
    },
    "required": ["scientific_name", "marker_metabolites", "pathway_prior", "refs"],
}

_DRAFT_SYSTEM = """你是药用植物次生代谢领域的知识起草助手。你的任务是为一个物种生成
物种知识 profile 的 JSON 草稿（供 PhytoReason-Agent 使用）。

严格规则（保守性原则，违反即失败）:
1. marker_metabolites: 该物种的标志性次生代谢物 3-8 个；name 用英文小写（下划线分隔）；
   formula/mz 只填你有把握的，不确定填 null；fragments 只列文献明确的诊断碎片，最多 3 个；
   source 三选一: 药典（中国药典收录）> 文献（有文献支撑）> 推测（其余）。
2. pathway_prior: 1-3 条已知生源通路；key 用 snake_case 英文；enzymes 只列该通路通识明确的酶；
   known_regulators 除非文献明确给出已知转录因子，一律留空数组（不编造调控因子）。
3. refs: 只引用【用户提供的文献列表】中真实出现的 PMID（格式 "PMID: 数字"），
   绝不编造 PMID；如无文献可用，写 "待人工确认"。
4. knowledge_gaps: 明确列出该物种缺失的知识（如无基因组、无已知调控因子等）。
5. 只输出 JSON 对象本身，不要任何额外文字、markdown 代码块或注释。"""


@dataclass
class OnboardingDraft:
    """Auto-onboarding 产物：可直接写入 species_profiles/ 的草稿。"""

    slug: str
    yaml_text: str
    warnings: list[str] = field(default_factory=list)
    source_pmids: list[str] = field(default_factory=list)
    source_titles: list[str] = field(default_factory=list)


# ── taxid 解析 ────────────────────────────────────────────

def _resolve_taxid(scientific_name: str, taxonomy_id: int | None) -> int | None:
    """taxid：显式参数 > PLANT_TAXIDS 查表 > NCBI（网络，失败返回 None）。"""
    if taxonomy_id:
        return int(taxonomy_id)
    key = scientific_name.lower().strip()
    if key in PLANT_TAXIDS:
        return int(PLANT_TAXIDS[key])
    try:
        from phyto_reason.knowledge.ortholog_finder import OrthologFinder
        return OrthologFinder()._resolve_taxid(scientific_name) or None
    except Exception:
        return None


# ── 文献检索（失败不阻塞）────────────────────────────────

def _search_literature(scientific_name: str, n_abstracts: int) -> tuple[list[Any], list[str]]:
    """返回 (检索结果, 警告)。网络失败时警告并返回空。"""
    try:
        result = PubMedRetriever().smart_search(
            f"{scientific_name} secondary metabolites biosynthesis", max_results=n_abstracts
        )
        return result.results, []
    except Exception as e:
        return [], [f"PubMed 检索失败（{e}）——草稿将无文献支撑，请人工确认所有条目"]


# ── LLM 起草 ──────────────────────────────────────────────

def _build_user_prompt(
    scientific_name: str,
    common_name: str,
    genome_status: str,
    taxonomy_id: int | None,
    literature: list[Any],
) -> str:
    lit_block = "\n".join(
        f"- PMID {r.pmid}: {r.title}"
        + (f" — {r.snippet[:120]}" if getattr(r, "snippet", "") else "")
        for r in literature[:8]
    ) or "（无可用文献）"
    return f"""物种: {scientific_name}
俗名: {common_name or "未知"}
NCBI taxonomy_id: {taxonomy_id if taxonomy_id else "未知"}
基因组状态: {genome_status}

可用文献（refs 只能引用这里的 PMID）:
{lit_block}

请生成该物种的 profile JSON 草稿。"""


def _sanitize_refs(refs: list[str], retrieved_pmids: list[str], common_name: str) -> tuple[list[str], list[str]]:
    """引用清洗：只保留检索结果真实 PMID / 药典类条目，丢弃编造 PMID。

    返回 (清洗后 refs, 警告)。
    """
    warnings: list[str] = []
    clean: list[str] = []
    retrieved = set(retrieved_pmids)

    for ref in refs:
        if not isinstance(ref, str) or not ref.strip():
            continue
        if len(ref.strip()) < 2:
            continue  # LLM 偶尔把占位文本拆成单字，静默丢弃
        m = re.match(r"^PMID[:\s]*(\d{4,9})$", ref.strip())
        if m:
            pmid = m.group(1)
            if pmid in retrieved:
                clean.append(f"PMID: {pmid}")
            else:
                warnings.append(f"丢弃未出现在检索结果中的 PMID: {pmid}（防编造）")
        elif "药典" in ref:
            clean.append(ref)
        elif "NCBI" in ref or "Taxonomy" in ref:
            clean.append(ref)
        else:
            warnings.append(f"丢弃无法核实的引用: {ref[:60]}")

    if common_name and not any("药典" in r for r in clean):
        clean.append(f"中国药典 2020 版（{common_name}）")
        warnings.append("自动添加药典条目（标志代谢物需人工对照药典确认）")
    return clean, warnings


def _render_yaml(
    scientific_name: str,
    common_name: str,
    genome_status: str,
    taxonomy_id: int | None,
    data: dict[str, Any],
    refs: list[str],
) -> tuple[str, int]:
    """按 species_profiles/README.md 的 schema 渲染 YAML。

    json_object 模式下 LLM 可能输出结构不规范的字段（如 pathway_prior 为列表、
    mz 为字符串等），此处全部防御性处理：坏条目丢弃并计数。
    返回 (yaml_text, 丢弃条目数)。
    """
    meta: dict[str, Any] = {"scientific_name": scientific_name}
    if common_name:
        meta["common_name"] = common_name
    if taxonomy_id:
        meta["taxonomy_id"] = int(taxonomy_id)
    meta["genome_status"] = genome_status

    dropped = 0

    def _to_float(v: Any) -> float | None:
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    # ── marker_metabolites（可能为 dict 或 list）──────────
    raw_mets = data.get("marker_metabolites") or []
    if isinstance(raw_mets, dict):
        raw_mets = list(raw_mets.values())
    metabolites = []
    for m in raw_mets:
        if not isinstance(m, dict) or not m.get("name"):
            dropped += 1
            continue
        entry: dict[str, Any] = {
            "name": str(m["name"]).strip().lower().replace(" ", "_"),
            "formula": m.get("formula") or None,
            "mz": _to_float(m.get("mz")),
            "adducts": m.get("adducts") or [],
            "fragments": [
                {"mz": _to_float(f.get("mz")), "annotation": f.get("annotation", "")}
                for f in (m.get("fragments") or [])[:3]
                if isinstance(f, dict) and _to_float(f.get("mz"))
            ],
            "reference_rt": None,
            "source": m.get("source", "推测"),
        }
        if not entry["mz"]:
            entry["mz"] = None
        metabolites.append(entry)

    # ── pathway_prior（模型可能输出列表 → 尝试转换）──────
    raw_pw = data.get("pathway_prior") or {}
    if isinstance(raw_pw, list):
        converted: dict[str, Any] = {}
        for item in raw_pw:
            if not isinstance(item, dict):
                dropped += 1
                continue
            name = (item.get("pathway") or item.get("name")
                    or item.get("pathway_name") or "")
            if isinstance(name, str) and name.strip():
                converted[name] = item
        raw_pw = converted
    elif not isinstance(raw_pw, dict):
        raw_pw = {}

    pathway_prior: dict[str, Any] = {}
    for pw, info in raw_pw.items():
        if not pw or not isinstance(info, dict):
            dropped += 1
            continue
        pathway_prior[str(pw).strip().lower().replace(" ", "_")] = {
            "genes": info.get("genes") or [],
            "enzymes": info.get("enzymes") or [],
            "known_regulators": info.get("known_regulators") or [],  # 显式空 = 未知
        }

    text = yaml.safe_dump(
        {
            "species": meta,
            "marker_metabolites": metabolites,
            "pathway_prior": pathway_prior,
            "refs": refs,
        },
        allow_unicode=True,
        sort_keys=False,
        default_flow_style=False,
    )
    return "# 由 Auto-onboarding 生成（LLM 起草 + SafetyGuard 清洗），请人工确认后使用。\n" + text, dropped


# ── 主流程 ────────────────────────────────────────────────

def draft_species_profile(
    scientific_name: str,
    common_name: str = "",
    genome_status: str = "unavailable",
    taxonomy_id: int | None = None,
    llm: LLMClient | None = None,
    n_abstracts: int = 5,
) -> OnboardingDraft:
    """生成物种 profile 草稿（LLM 起草，不落盘）。

    Raises:
        ValueError: LLM 未配置或不可用（需要 .env 的 OPENAI_API_KEY / LLM_MODEL）。
    """
    if genome_status not in {"available", "partial", "unavailable"}:
        raise ValueError(f"genome_status 必须为 available|partial|unavailable，得到 {genome_status!r}")

    client = llm or LLMClient()
    if not client.is_available():
        raise ValueError(
            "LLM 不可用：请在项目根 .env 配置 OPENAI_API_KEY（DeepSeek 兼容接口），"
            "或传入已配置的 llm 实例"
        )

    warnings: list[str] = []
    slug = scientific_name.strip().lower().replace(" ", "_")

    # 1. taxid
    taxid = _resolve_taxid(scientific_name, taxonomy_id)
    if not taxid:
        warnings.append("未解析到 taxonomy_id（PLANT_TAXIDS 未收录且 NCBI 不可用）——可后续手动补")

    # 2. 文献
    literature, lit_warnings = _search_literature(scientific_name, n_abstracts)
    warnings.extend(lit_warnings)
    retrieved_pmids = [str(r.pmid) for r in literature if getattr(r, "pmid", None)]

    # 3. LLM 起草
    messages = [
        {"role": "system", "content": _DRAFT_SYSTEM},
        {"role": "user", "content": _build_user_prompt(
            scientific_name, common_name, genome_status, taxid, literature
        )},
    ]
    data = client.chat_structured(messages, response_format=DRAFT_SCHEMA, temperature=0.1)
    if data.get("parse_error"):
        raise ValueError(f"LLM 结构化输出失败: {data.get('error') or data.get('raw', '')[:200]}")

    # 4. SafetyGuard 校验
    verdict = SafetyGuard.check(
        yaml.safe_dump(data, allow_unicode=True),
        context={"source_pmids": retrieved_pmids},
    )
    if not verdict.passed:
        warnings.extend(f"[SafetyGuard] {w}" for w in verdict.warnings[:5])

    # 5. 引用清洗 + 渲染
    refs, ref_warnings = _sanitize_refs(data.get("refs") or [], retrieved_pmids, common_name)
    warnings.extend(ref_warnings)
    if not literature:
        warnings.append("无文献支撑：请人工核对所有代谢物与通路条目")

    yaml_text, dropped = _render_yaml(
        scientific_name, common_name, genome_status, taxid, data, refs
    )
    if dropped:
        warnings.append(f"LLM 输出中 {dropped} 个结构不规范条目已丢弃（请人工核对）")
    return OnboardingDraft(
        slug=slug,
        yaml_text=yaml_text,
        warnings=warnings,
        source_pmids=retrieved_pmids,
        source_titles=[r.title for r in literature[:8]],
    )


def write_profile(
    draft: OnboardingDraft,
    profiles_dir: Path | None = None,
    force: bool = False,
) -> Path:
    """将草稿写入 species_profiles/<slug>.yaml。

    Raises:
        FileExistsError: 目标已存在且未传 force=True。
    """
    target_dir = profiles_dir or _PROFILES_DIR
    target = target_dir / f"{draft.slug}.yaml"
    if target.exists() and not force:
        raise FileExistsError(
            f"{target} 已存在。用 force=True 覆盖，或改用其他 slug。"
        )
    target.write_text(draft.yaml_text, encoding="utf-8")
    return target
