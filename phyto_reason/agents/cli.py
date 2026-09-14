"""
PhytoReason-Agent CLI — AgentOrchestrator 的薄包装（唯一 CLI 入口）。

用法:
  plantomics --help
  plantomics ask "黄芩素在黄芩里的已知调控因子有哪些？" --species "黄芩" --metabolite baicalein
  plantomics species list
  plantomics species show zanthoxylum_nitidum
"""

from __future__ import annotations

import sys

import click

from phyto_reason.agent.orchestrator import AgentOrchestrator
from phyto_reason.knowledge.species_registry import get_species_registry


@click.group()
def main() -> None:
    """PhytoReason-Agent v5.0 — 多物种药用植物次生代谢调控科研推理 Agent."""


@main.command()
@click.argument("message")
@click.option("--species", default="", help="会话物种（学名/俗名/slug，注入 profile 或降级声明）")
@click.option("--metabolite", default="", help="目标代谢物")
@click.option("--session", "session_id", default="", help="会话 ID（多轮对话续接）")
@click.option("--compact/--no-compact", default=False, help="长会话压缩模式")
def ask(message: str, species: str, metabolite: str, session_id: str, compact: bool) -> None:
    """向 Agent 提问（L0-L3 自动路由）。"""
    orchestrator = AgentOrchestrator()
    sid = session_id or ""

    # 会话物种/目标代谢物注入（进入 system prompt 的物种感知模板）
    if species or metabolite:
        sid = orchestrator.upload_data(
            sid, species=species or "", target_metabolite=metabolite or ""
        )

    result = orchestrator.handle(message, session_id=sid or None)
    click.echo(result.get("response", "(无输出)"))
    click.echo(f"\n[会话 {result.get('session_id', '')}] "
               f"工具调用: {', '.join(result.get('tool_calls_made', [])) or '无'}",
               err=True)
    if result.get("error"):
        click.echo(f"[错误] {result['error']}", err=True)


@main.command()
@click.argument("scientific_name")
@click.option("--common", "common_name", default="", help="俗名/中文名")
@click.option("--genome", "genome_status", default="unavailable",
              type=click.Choice(["available", "partial", "unavailable"]),
              help="基因组状态")
@click.option("--taxid", type=int, default=None, help="NCBI taxonomy ID（不填则自动解析）")
@click.option("--write", "do_write", is_flag=True, help="确认后写入 species_profiles/")
@click.option("--force", is_flag=True, help="覆盖已存在的 profile")
def onboarding(scientific_name: str, common_name: str, genome_status: str,
               taxid: int | None, do_write: bool, force: bool) -> None:
    """Auto-onboarding: LLM 起草新物种 profile（需配置 .env 的 OPENAI_API_KEY）。

    流程: taxid 解析 → PubMed 检索 → LLM 起草 → SafetyGuard 引用清洗 → YAML。
    默认只打印草稿；--write 落盘 species_profiles/<slug>.yaml。
    """
    from phyto_reason.knowledge.species_onboarding import draft_species_profile, write_profile

    try:
        draft = draft_species_profile(
            scientific_name, common_name=common_name,
            genome_status=genome_status, taxonomy_id=taxid,
        )
    except ValueError as e:
        click.echo(f"[错误] {e}", err=True)
        raise SystemExit(1)

    click.echo(draft.yaml_text)
    if draft.warnings:
        click.echo("\n[警告]", err=True)
        for w in draft.warnings:
            click.echo(f"  ⚠ {w}", err=True)
    if draft.source_pmids:
        click.echo(f"\n[引用文献] {len(draft.source_pmids)} 篇: "
                   f"{', '.join(draft.source_pmids[:5])}", err=True)

    if do_write:
        try:
            path = write_profile(draft, force=force)
            click.echo(f"\n✓ 已写入 {path}", err=True)
        except FileExistsError as e:
            click.echo(f"[错误] {e}", err=True)
            raise SystemExit(1)


@main.command("predict-tf")
@click.argument("fasta")
@click.option("--out", default="",
              help="输出注释 CSV 路径；缺省 <fasta目录>/tf_prediction.annotation.csv")
@click.option("--pfam", "pfam_source", default="",
              help="本地 Pfam-A.hmm(.gz) 路径；缺省自动从 EBI 下载（一次性，~520MB）")
@click.option("--cache", "cache_dir", default="", help="HMM 子集缓存目录")
@click.option("--min-evalue", default=1e-5, type=float, help="HMM 扫描 E 值阈值")
@click.option("--min-aa", default=30, type=int, help="核酸输入最短 ORF 氨基酸数")
@click.option("--cpus", default=4, type=int, help="并行线程数")
def predict_tf(fasta: str, out: str, pfam_source: str, cache_dir: str,
               min_evalue: float, min_aa: int, cpus: int) -> None:
    """无注释文件时预测转录因子：FASTA（转录本或蛋白）→ Pfam HMM → 家族判定。

    输出的 tf_prediction.annotation.csv 与表达矩阵放在同一目录时，
    L3 验证脚本会自动拾取为 TF 注释（无需改管线）。
    """
    from phyto_reason.tf_prediction.predictor import predict_from_fasta

    try:
        summary = predict_from_fasta(
            fasta,
            out_path=out or None,
            cache_dir=cache_dir or None,
            pfam_source=pfam_source or None,
            min_evalue=min_evalue,
            min_aa=min_aa,
            cpus=cpus,
        )
    except RuntimeError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"输入: {summary['input_fasta']}（{summary['mode']}，{summary['n_sequences']} 序列）")
    click.echo(f"预测 TF: {summary['n_predicted_tf']} 条")
    for family, count in list(summary["family_counts"].items())[:15]:
        click.echo(f"  {family}: {count}")
    if len(summary["family_counts"]) > 15:
        click.echo(f"  ... 共 {len(summary['family_counts'])} 个家族")
    click.echo(f"注释输出: {summary['annotation_csv']}")


@main.command()
@click.argument("spectrum_file")
@click.option("--species", default="", help="物种（profile 标志物锚定）")
@click.option("--target", "target_metabolite", default="", help="目标代谢物")
@click.option("--ppm", "tolerance_ppm", default=10.0, type=float, help="质量容差")
@click.option("--scan", type=int, default=1, help="注释第 N 张谱图（默认 1，0=全部）")
@click.option("--mode", "polarity", default="",
              type=click.Choice(["", "positive", "negative"]),
              help="电离模式（默认自动：文件名 pos/neg 或 MGF CHARGE）")
@click.option("--summary", is_flag=True,
              help="紧凑模式：每张谱一行摘要（谱库批量筛选用，不烧 token）")
def annotate(spectrum_file: str, species: str, target_metabolite: str,
             tolerance_ppm: float, scan: int, polarity: str, summary: bool) -> None:
    """未知物鉴定：MS² 谱图 → 分子式/结构类/候选代谢物 + MSI 置信度。"""
    from phyto_reason.metabolomics.spectrum_io import load_spectrum
    from phyto_reason.metabolomics.annotator import annotate_spectrum
    from phyto_reason.agent.tool_handlers import _format_annotation

    try:
        spectra = load_spectrum(spectrum_file, polarity=polarity)
    except (ValueError, RuntimeError) as e:
        click.echo(f"[错误] {e}", err=True)
        raise SystemExit(1)
    if not spectra:
        click.echo("[错误] 未解析到谱图", err=True)
        raise SystemExit(1)

    # 跳过无碎片谱图（MSP 的 Num Peaks: 0 / 只有前体信息的条目）
    valid = [s for s in spectra if s.fragments]
    if len(valid) != len(spectra):
        click.echo(f"[跳过 {len(spectra) - len(valid)} 张无碎片谱图（Num Peaks: 0）]", err=True)
    spectra = valid
    if not spectra:
        click.echo("[错误] 全部谱图均无碎片（无法做 MS² 鉴定）", err=True)
        raise SystemExit(1)

    if summary:
        for i, s in enumerate(spectra, 1):
            r = annotate_spectrum(s, species=species,
                                  target_metabolite=target_metabolite,
                                  tolerance_ppm=tolerance_ppm,
                                  polarity=polarity)
            top = r.candidates[0] if r.candidates else None
            cls = f"{top.class_name} ({top.name})" if top and top.name else (
                top.class_name if top else "无类匹配")
            msi = f"L{top.msi_level}" if top else "L4"
            hits = ",".join(h["name"] for h in r.profile_hits)
            pol = s.polarity or polarity
            click.echo(
                f"#{i} m/z {s.precursor_mz:9.4f} {r.adduct} {pol:8s} | "
                f"{cls:40s} | MSI {msi} | {top.score if top else 0.0:.2f}"
                + (f" | profile: {hits}" if hits else "")
            )
        click.echo(f"[共 {len(spectra)} 张谱图]", err=True)
        return

    if scan == 0:
        blocks = [
            _format_annotation(annotate_spectrum(
                s, species=species, target_metabolite=target_metabolite,
                tolerance_ppm=tolerance_ppm, polarity=polarity),
                label=f"谱图 {i + 1}（{s.scan_id or f'{s.precursor_mz:.4f}'}）")
            for i, s in enumerate(spectra)
        ]
        click.echo("\n\n".join(blocks))
    else:
        if scan > len(spectra):
            click.echo(f"[错误] 文件只有 {len(spectra)} 张谱图", err=True)
            raise SystemExit(1)
        spec = spectra[scan - 1]
        result = annotate_spectrum(
            spec, species=species, target_metabolite=target_metabolite,
            tolerance_ppm=tolerance_ppm, polarity=polarity)
        click.echo(_format_annotation(result, label=spec.scan_id or f"m/z {spec.precursor_mz:.4f}"))


@main.group()
def species() -> None:
    """物种 profile 查询。"""


@species.command("list")
def species_list() -> None:
    """列出已收录物种。"""
    registry = get_species_registry()
    for slug in registry.list_species():
        profile = registry.get(slug)
        click.echo(f"{slug} — {profile.scientific_name} ({profile.common_name})")


@species.command("show")
@click.argument("name")
def species_show(name: str) -> None:
    """显示物种 profile 摘要与降级链层级。"""
    registry = get_species_registry()
    profile = registry.get(name)
    if profile is None:
        click.echo(registry.describe_gap(name))
        sys.exit(1)
    click.echo(profile.to_summary())
    result, level = registry.knowledge_for(name)
    click.echo(f"知识层级: {level}")
    click.echo(f"备注: {result.note}")


if __name__ == "__main__":
    main()
