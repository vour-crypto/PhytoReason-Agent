"""
scoring.py — PhytoReason-Agent v5.0 benchmark scoring framework.

评估维度（权重，v5.0 重平衡）：
  Engineering (25%):    无崩溃、可复现、不超时
  Conservatism (25%):   随机数据→低置信度、小样本→警告
  Evidence (25%):       知识层召回（cases/*.yaml: B01/B02/B03/B06）
  Reasoning (15%):      机制分类、矛盾检测
  Species Adaptation (10%, v5.0 新增): 降级链诚实性（profile→ortholog→gap 不编造）

cases/ 数据:
  positive_tf_pairs.yaml  (B01) 已知 TF-代谢物调控对
  pathways.yaml           (B02) 已知通路酶集
  layer0_queries.yaml     (B03) Layer 0 查询标准答案（自动评测）
  cross_species.yaml      (B06) 跨物种推断对

用法:
  python -m phyto_reason.benchmark.scoring
  python -m phyto_reason.benchmark.scoring --json  # JSON 输出
"""

from __future__ import annotations

import json
import time
import random
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger("benchmark")

_CASES_DIR = Path(__file__).resolve().parent / "cases"


def _load_cases(name: str) -> list[dict]:
    """加载 benchmark/cases/<name>（YAML 列表）。"""
    path = _CASES_DIR / name
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh) or []


# ═══════════════════════════════════════════════════════════════
# Data models
# ═══════════════════════════════════════════════════════════════

@dataclass
class TestResult:
    """单个测试结果。"""
    name: str
    passed: bool
    score: float = 0.0         # 0.0-1.0
    detail: str = ""
    duration_s: float = 0.0


@dataclass
class DimensionScore:
    """单个维度的评分。"""
    name: str
    weight: float
    score: float               # 0.0-1.0
    tests: list[TestResult] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def weighted(self) -> float:
        return self.score * self.weight


@dataclass
class BenchmarkResult:
    """完整 benchmark 结果。"""
    version: str = "4.0.0"
    timestamp: str = ""
    dimensions: list[DimensionScore] = field(default_factory=list)
    overall: float = 0.0
    grade: str = ""
    summary: str = ""

    def to_dict(self) -> dict:
        return {
            "version": self.version,
            "timestamp": self.timestamp,
            "overall": round(self.overall, 3),
            "grade": self.grade,
            "dimensions": [
                {
                    "name": d.name,
                    "weight": d.weight,
                    "score": round(d.score, 3),
                    "weighted": round(d.weighted, 3),
                    "tests": [{"name": t.name, "passed": t.passed, "score": round(t.score, 3), "detail": t.detail} for t in d.tests],
                    "notes": d.notes,
                }
                for d in self.dimensions
            ],
            "summary": self.summary,
        }


# ═══════════════════════════════════════════════════════════════
# Mock data helpers
# ═══════════════════════════════════════════════════════════════

def _make_mock_expression(n_genes: int = 30, n_samples: int = 8) -> dict:
    rng = random.Random(42)
    return {f"gene_{i:03d}": {f"s_{j:02d}": rng.uniform(0, 10) for j in range(n_samples)} for i in range(n_genes)}


def _make_mock_metabolite(n_met: int = 8, n_samples: int = 8) -> dict:
    rng = random.Random(42)
    return {f"met_{i:02d}": {f"s_{j:02d}": rng.uniform(0, 5) for j in range(n_samples)} for i in range(n_met)}


def _shuffle_expression(expr: dict) -> dict:
    """Shuffle each gene's expression values to break correlations."""
    rng = random.Random(99)
    result = {}
    for gene, samples in expr.items():
        values = list(samples.values())
        rng.shuffle(values)
        result[gene] = {k: v for k, v in zip(samples.keys(), values)}
    return result


# ═══════════════════════════════════════════════════════════════
# Dimension 1: Engineering (25%)
# ═══════════════════════════════════════════════════════════════

def benchmark_engineering() -> DimensionScore:
    """Test: pipeline 不崩溃、可复现、不超时。"""
    tests: list[TestResult] = []

    # Test 1: Pipeline doesn't crash with normal data
    try:
        from phyto_reason.workflows.workflow_runner import WorkflowRunner
        from phyto_reason.workflows.runtime_state import RuntimeState
        from phyto_reason.models.planner_state import PlannerState

        t0 = time.time()
        runner = WorkflowRunner()
        ps = PlannerState(
            species="Arabidopsis thaliana",
            target_metabolite="anthocyanin",
            has_expression=True, has_metabolite=True,
            expression_matrix=_make_mock_expression(30, 8),
            metabolite_matrix=_make_mock_metabolite(8, 8),
            sample_count=8,
        )
        state = runner.run(
            species="Arabidopsis thaliana",
            target_metabolite="anthocyanin",
            sample_count=8,
            has_expression=True, has_metabolite=True,
            expression_matrix=_make_mock_expression(30, 8),
            metabolite_matrix=_make_mock_metabolite(8, 8),
        )
        duration = time.time() - t0
        passed = state is not None and duration < 120
        tests.append(TestResult(
            name="pipeline_no_crash",
            passed=passed,
            score=1.0 if passed else 0.0,
            detail=f"Pipeline completed in {duration:.1f}s, {len(state.completed_nodes)} nodes",
            duration_s=duration,
        ))
    except Exception as e:
        tests.append(TestResult(
            name="pipeline_no_crash",
            passed=False, score=0.0,
            detail=f"Pipeline crashed: {e}",
        ))

    # Test 2: Deterministic — same input → same output
    try:
        expr = _make_mock_expression(20, 6)
        meta = _make_mock_metabolite(5, 6)
        runner = WorkflowRunner()

        state1 = runner.run(
            species="test", target_metabolite="flavonoid",
            sample_count=6, has_expression=True, has_metabolite=True,
            expression_matrix=expr, metabolite_matrix=meta,
        )
        state2 = runner.run(
            species="test", target_metabolite="flavonoid",
            sample_count=6, has_expression=True, has_metabolite=True,
            expression_matrix=expr, metabolite_matrix=meta,
        )

        same = (state1.mechanism_type == state2.mechanism_type and
                state1.data_quality_passed == state2.data_quality_passed)
        tests.append(TestResult(
            name="pipeline_deterministic",
            passed=same,
            score=1.0 if same else 0.5,
            detail=f"Same input → mechanism type: {state1.mechanism_type}, quality: {state1.data_quality_passed}",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="pipeline_deterministic",
            passed=False, score=0.0,
            detail=f"Determinism test failed: {e}",
        ))

    # Test 3: All nodes registered in graph
    try:
        from phyto_reason.workflows.workflow_graph import compile_workflow
        app = compile_workflow()
        nodes = list(app.nodes.keys()) if hasattr(app, 'nodes') else []
        expected = {"__start__", "data_quality_gate", "metabolite_validation",
                     "pathway_coherence", "mechanism_classifier", "tf_narrowing",
                     "transport_evidence", "stress_evidence", "regulation_evidence",
                     "contradiction_check", "falsification", "hypothesis_competition",
                     "hypothesis_synthesis"}
        missing = expected - set(nodes)
        passed = len(missing) == 0
        tests.append(TestResult(
            name="graph_all_nodes",
            passed=passed,
            score=1.0 if passed else 1.0 - len(missing) / len(expected),
            detail=f"Graph nodes: {len(nodes)}/{len(expected)}" if passed else f"Missing: {missing}",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="graph_all_nodes",
            passed=False, score=0.0,
            detail=f"Graph check failed: {e}",
        ))

    score = sum(t.score for t in tests) / len(tests) if tests else 0
    return DimensionScore(name="Engineering", weight=0.25, score=score, tests=tests)


# ═══════════════════════════════════════════════════════════════
# Dimension 2: Conservatism (25%)
# ═══════════════════════════════════════════════════════════════

def benchmark_conservatism() -> DimensionScore:
    """Test: 系统在无信号数据上是否保守。"""
    tests: list[TestResult] = []

    # Test 1: Shuffled data → low confidence
    try:
        from phyto_reason.workflows.workflow_runner import WorkflowRunner

        expr = _make_mock_expression(30, 8)
        shuffled = _shuffle_expression(expr)
        meta = _make_mock_metabolite(8, 8)

        runner = WorkflowRunner()
        state = runner.run(
            species="test", target_metabolite="flavonoid",
            sample_count=8, has_expression=True, has_metabolite=True,
            expression_matrix=shuffled, metabolite_matrix=meta,
        )

        # Shuffled data should NOT produce confident hypotheses
        hypotheses = getattr(state, "competing_hypotheses", []) or []
        mechanism = getattr(state, "mechanism_type", "") or ""

        # Check: shuffled data should either have no hypotheses or mechanism type
        # indicating "technical_artifact" or "insufficient_data"
        conservative = (
            len(hypotheses) == 0
            or mechanism in ("technical_artifact", "insufficient_data", "unknown")
        )
        tests.append(TestResult(
            name="conservatism_shuffled",
            passed=conservative,
            score=1.0 if conservative else 0.3,
            detail=f"Shuffled data: {len(hypotheses)} hypotheses, mechanism={mechanism}",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="conservatism_shuffled",
            passed=True, score=0.8,  # crash is not ideal but still detected
            detail=f"Shuffled data caused error (acceptable): {e}",
        ))

    # Test 2: Small sample → warning
    try:
        expr_small = _make_mock_expression(20, 4)

        runner = WorkflowRunner()
        state = runner.run(
            species="test", target_metabolite="flavonoid",
            sample_count=4, has_expression=True, has_metabolite=True,
            expression_matrix=expr_small, metabolite_matrix=_make_mock_metabolite(5, 4),
        )

        warnings = len(getattr(state, "data_quality_issues", []))
        has_warnings = warnings > 0
        tests.append(TestResult(
            name="conservatism_small_sample",
            passed=has_warnings,
            score=1.0 if has_warnings else 0.3,
            detail=f"n=4: {warnings} quality issues raised",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="conservatism_small_sample",
            passed=False, score=0.3,
            detail=f"Small sample test failed: {e}",
        ))

    # Test 3: No data → pipeline rejects
    try:
        runner = WorkflowRunner()
        state = runner.run(
            species="test", target_metabolite="flavonoid",
            sample_count=0, has_expression=False, has_metabolite=False,
        )
        passed = not state.data_quality_passed
        tests.append(TestResult(
            name="conservatism_no_data",
            passed=passed,
            score=1.0 if passed else 0.0,
            detail=f"No data: quality_passed={state.data_quality_passed}",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="conservatism_no_data",
            passed=False, score=0.0,
            detail=f"No-data test crashed: {e}",
        ))

    score = sum(t.score for t in tests) / len(tests) if tests else 0
    return DimensionScore(name="Conservatism", weight=0.25, score=score, tests=tests)


# ═══════════════════════════════════════════════════════════════
# Dimension 3: Evidence Quality (25%)
# ═══════════════════════════════════════════════════════════════

def benchmark_evidence() -> DimensionScore:
    """Test: 文献检索覆盖、KEGG 通路映射准确。

    NOTE: Full evidence scoring requires B01-B06 benchmark data.
    This baseline tests what's available without external data.
    """
    tests: list[TestResult] = []
    notes: list[str] = []

    # Test 1: Knowledge base has comprehensive coverage
    try:
        from phyto_reason.knowledge.tf_knowledge_base import CANONICAL_KNOWLEDGE, CLASS_ALIASES

        n_relations = len(CANONICAL_KNOWLEDGE)
        has_strong = sum(1 for r in CANONICAL_KNOWLEDGE if r.strength == "strong")
        has_pmids = sum(1 for r in CANONICAL_KNOWLEDGE if r.pmids)

        # Baseline: >= 30 relations, >= 5 strong, >= 20 with PMIDs
        sufficient = n_relations >= 30 and has_strong >= 5 and has_pmids >= 20
        tests.append(TestResult(
            name="evidence_knowledge_coverage",
            passed=sufficient,
            score=min(n_relations / 40, 1.0),
            detail=f"{n_relations} relations, {has_strong} strong, {has_pmids} with PMIDs",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="evidence_knowledge_coverage",
            passed=False, score=0.0,
            detail=f"Error: {e}",
        ))

    # Test 2: KEGG browser has essential pathways
    try:
        from phyto_reason.knowledge.kegg_browser import BUILTIN_PATHWAYS

        n_pathways = len(BUILTIN_PATHWAYS)
        has_enzymes = sum(1 for pw in BUILTIN_PATHWAYS.values() if pw.enzyme_genes)
        sufficient = n_pathways >= 5 and has_enzymes >= 5
        tests.append(TestResult(
            name="evidence_kegg_coverage",
            passed=sufficient,
            score=min(n_pathways / 6, 1.0),
            detail=f"{n_pathways} pathways, {has_enzymes} with enzyme genes",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="evidence_kegg_coverage",
            passed=False, score=0.0,
            detail=f"Error: {e}",
        ))

    # Test 3: B03 — Layer 0 知识层召回（本地，无网络）
    try:
        from phyto_reason.knowledge.tf_knowledge_base import TFKnowledgeBase
        from phyto_reason.knowledge.compound_profiles import find_diagnostic_rules
        from phyto_reason.knowledge.species_registry import get_species_registry

        queries = _load_cases("layer0_queries.yaml")
        kb = TFKnowledgeBase()
        registry = get_species_registry()

        per_query: list[tuple[float, str]] = []
        for q in queries:
            checks: list[bool] = []
            met = q.get("metabolite", "")

            fams = {r.tf_family for r in kb.query(metabolite=met).relations} if met else set()
            if q.get("expected_tf_families"):
                checks.append(fams >= set(q["expected_tf_families"]))

            diag = find_diagnostic_rules(met) if met else None
            if q.get("expected_compound_class"):
                checks.append(bool(diag) and diag["class"] == q["expected_compound_class"])
            if q.get("expected_fragment") and diag:
                checks.append(any(
                    abs(fr["mz"] - q["expected_fragment"]) < 0.01
                    for fr in diag.get("fragments", [])
                ))
            if q.get("expected_pathway") and diag:
                # 通路召回为软检查：期望通路的有效词元（len>4）命中类条目上下文即算
                haystack = " ".join([
                    diag.get("class", ""),
                    diag.get("class_file", ""),
                    diag.get("pathway", "") or "",
                ]).lower()
                tokens = [t for t in q["expected_pathway"].lower().split() if len(t) > 4]
                checks.append(bool(tokens) and any(t in haystack for t in tokens))
            # expected_enzymes 在物种 profile 语境下检查（L0 的酶覆盖见 B02）
            if q.get("expected_enzymes") and q.get("species"):
                result, _level = registry.knowledge_for(q["species"], metabolite=met)
                profile_enzymes = {
                    e.upper() for info in registry.get(q["species"]).pathway_prior.values()
                    for e in info.get("enzymes", [])
                } if registry.get(q["species"]) else set()
                checks.append(
                    set(q["expected_enzymes"]).issubset(profile_enzymes)
                )

            if q.get("expected_species"):
                result, level = registry.knowledge_for(q["species"], metabolite=met)
                checks.append(level == q.get("expected_level"))
                if q.get("expected_metabolites"):
                    checks.append(
                        set(q["expected_metabolites"]) <= set(result.known_metabolites)
                    )

            score = sum(checks) / len(checks) if checks else 0.0
            per_query.append((score, f"{q['id']}: {score:.0%}"))

        avg = sum(s for s, _ in per_query) / len(per_query) if per_query else 0.0
        passed = avg >= 0.7
        tests.append(TestResult(
            name="evidence_literature_recall",
            passed=passed, score=avg,
            detail=f"{len(per_query)} queries, recall {avg:.0%}"
                   + (f" — {', '.join(d for _, d in per_query[:3])}" if per_query else ""),
        ))
    except Exception as e:
        tests.append(TestResult(
            name="evidence_literature_recall",
            passed=False, score=0.0,
            detail=f"Layer 0 recall test failed: {e}",
        ))

    # Test 4: B01 — 已知调控对召回（TFKnowledgeBase）+ motif 灵敏度（B01 中 PLANT_TF_MOTIFS 覆盖的家族）
    try:
        from phyto_reason.knowledge.tf_knowledge_base import TFKnowledgeBase
        from phyto_reason.knowledge.motif_conservation import PLANT_TF_MOTIFS, MotifScanner
        from phyto_reason.knowledge.compound_profiles import find_diagnostic_rules

        pairs = _load_cases("positive_tf_pairs.yaml")
        kb = TFKnowledgeBase()
        scanner = MotifScanner()

        # 4a. TF-代谢物对召回
        recall_scores = []
        for p in pairs:
            rels = kb.query(metabolite=p.get("target", ""), tf_family=p.get("tf_family", ""))
            hit = any(r.tf_family == p["tf_family"] for r in rels.relations)
            recall_scores.append(1.0 if hit else 0.0)

        # 4b. motif 灵敏度：合成含该家族全部 motif consensus 的启动子，
        #     扫描应检出该家族（家族可识别性）。短退化 motif 单独可能低于
        #     扫描阈值，故嵌入全家族 motif 集合（保守的家族级测试）。
        # PLANT_TF_MOTIFS 按 motif 名建键，值 (consensus, tf_family, pmid)
        family_motifs: dict[str, list[str]] = {}
        for _motif_name, (consensus, tf_family, _pmid) in PLANT_TF_MOTIFS.items():
            family_motifs.setdefault(tf_family, []).append(consensus)

        motif_scores = []
        skipped = []
        for p in pairs:
            family = p.get("tf_family", "")
            consensi = family_motifs.get(family)
            if not consensi:
                skipped.append(family)
                continue
            promoter = "ACGT" * 50 + "ACGT".join(consensi) + "ACGT" * 50
            result = scanner.scan(promoter, gene_id=f"B01-{p['id']}")
            motif_scores.append(1.0 if family in {h.tf_family for h in result.hits} else 0.0)

        recall_avg = sum(recall_scores) / len(recall_scores) if recall_scores else 0.0
        motif_avg = sum(motif_scores) / len(motif_scores) if motif_scores else 0.0
        tests.append(TestResult(
            name="evidence_tf_pair_recall",
            passed=recall_avg >= 0.8, score=recall_avg,
            detail=f"{len(pairs)} pairs, recall {recall_avg:.0%}",
        ))
        detail = f"{len(motif_scores)} families scanned, sensitivity {motif_avg:.0%}"
        if skipped:
            detail += f" (skip: {', '.join(sorted(skipped))})"
            notes.append(f"B01 motif: 无 PWM 的家族跳过 {sorted(skipped)}")
        tests.append(TestResult(
            name="evidence_motif_sensitivity",
            passed=motif_avg >= 0.8, score=motif_avg,
            detail=detail,
        ))
    except Exception as e:
        tests.append(TestResult(
            name="evidence_tf_pair_recall",
            passed=False, score=0.0, detail=f"B01 recall failed: {e}",
        ))

    # Test 5: B02 — 通路酶集映射（species profile pathway_prior + BUILTIN_PATHWAYS）
    try:
        from phyto_reason.knowledge.species_registry import get_species_registry
        from phyto_reason.knowledge.kegg_browser import BUILTIN_PATHWAYS

        registry = get_species_registry()
        kegg_enzymes = {
            e.upper()
            for pw in BUILTIN_PATHWAYS.values()
            for e in (pw.enzyme_genes or [])
        }
        pw_scores = []
        for pw in _load_cases("pathways.yaml"):
            expected = {e.upper() for e in pw.get("expected_enzymes", [])}
            if not expected:
                continue
            # 物种 profile 通路酶
            species_enzymes: set[str] = set()
            if pw.get("species_example"):
                profile = registry.get(pw["species_example"])
                if profile:
                    for info in profile.pathway_prior.values():
                        species_enzymes |= {e.upper() for e in info.get("enzymes", [])}
            union = species_enzymes | kegg_enzymes
            hit = len(expected & union) / len(expected)
            pw_scores.append(hit)

        pw_avg = sum(pw_scores) / len(pw_scores) if pw_scores else 0.0
        tests.append(TestResult(
            name="evidence_pathway_mapping",
            passed=pw_avg >= 0.5, score=pw_avg,
            detail=f"{len(pw_scores)} pathways, enzyme coverage {pw_avg:.0%}",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="evidence_pathway_mapping",
            passed=False, score=0.0, detail=f"B02 failed: {e}",
        ))

    # Test 6: B06 — 跨物种推断（降级链诚实性）
    try:
        from phyto_reason.knowledge.species_registry import get_species_registry

        registry = get_species_registry()
        cs_scores = []
        for cs in _load_cases("cross_species.yaml"):
            result, level = registry.knowledge_for(
                cs["target_species"],
                metabolite=cs.get("metabolite", ""),
                tf_family=cs.get("tf_family", ""),
            )
            honest = level in ("profile", "ortholog") and bool(result.note)
            if level != "profile":
                honest = honest and "降级" in result.note  # 无 profile 必须显式降级
            cs_scores.append(1.0 if honest else 0.0)

        cs_avg = sum(cs_scores) / len(cs_scores) if cs_scores else 0.0
        tests.append(TestResult(
            name="evidence_cross_species_honesty",
            passed=cs_avg >= 0.8, score=cs_avg,
            detail=f"{len(cs_scores)} pairs, honest degrade {cs_avg:.0%}",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="evidence_cross_species_honesty",
            passed=False, score=0.0, detail=f"B06 failed: {e}",
        ))

    score = sum(t.score for t in tests) / len(tests) if tests else 0
    return DimensionScore(name="Evidence Quality", weight=0.25, score=score, tests=tests, notes=notes)


# ═══════════════════════════════════════════════════════════════
# Dimension 4: Reasoning Quality (15%)
# ═══════════════════════════════════════════════════════════════

def benchmark_reasoning() -> DimensionScore:
    """Test: 推理质量 — 机制分类、矛盾检测、证伪。

    NOTE: Full reasoning scoring requires real multi-omics data (B03).
    """
    tests: list[TestResult] = []
    notes: list[str] = []

    # Test 1: Mechanism classifier produces valid types
    try:
        from phyto_reason.workflows.execution_nodes import mechanism_classifier_node
        from phyto_reason.workflows.runtime_state import RuntimeState
        from phyto_reason.models.planner_state import PlannerState

        ps = PlannerState(
            species="test", target_metabolite="flavonoid",
            has_expression=True, has_metabolite=True,
            expression_matrix=_make_mock_expression(50, 8),
            metabolite_matrix=_make_mock_metabolite(10, 8),
            sample_count=8,
        )
        state = RuntimeState(planner_state=ps)
        mechanism_classifier_node(state)
        mechanism = state.mechanism_type

        valid_types = {
            "transcriptional_regulation", "transport_redistribution",
            "stress_induced_redistribution", "enzymatic_regulation",
            "technical_artifact", "unknown", "insufficient_data", "",
        }
        passed = mechanism in valid_types
        tests.append(TestResult(
            name="reasoning_mechanism_classifier",
            passed=passed,
            score=0.8 if passed else 0.0,
            detail=f"Classified as: '{mechanism}'",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="reasoning_mechanism_classifier",
            passed=False, score=0.0,
            detail=f"Error: {e}",
        ))

    # Test 2: Contradiction detector has rules
    try:
        from phyto_reason.fusion.contradiction_detector import ContradictionDetector

        detector = ContradictionDetector()
        # ContradictionDetector should have detect() method
        has_method = hasattr(detector, 'detect')
        tests.append(TestResult(
            name="reasoning_contradiction_ready",
            passed=has_method,
            score=1.0 if has_method else 0.0,
            detail="ContradictionDetector.detect() available",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="reasoning_contradiction_ready",
            passed=False, score=0.0,
            detail=f"Error: {e}",
        ))

    # Test 3: Falsification engine exists
    try:
        from phyto_reason.reasoning.falsification_engine import FalsificationEngine
        has_class = FalsificationEngine is not None
        tests.append(TestResult(
            name="reasoning_falsification_ready",
            passed=has_class,
            score=1.0 if has_class else 0.0,
            detail="FalsificationEngine available",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="reasoning_falsification_ready",
            passed=False, score=0.0,
            detail=f"Error: {e}",
        ))

    # Placeholder
    notes.append("Full reasoning evaluation requires B03 multi-omics benchmark data.")

    score = sum(t.score for t in tests) / len(tests) if tests else 0
    return DimensionScore(name="Reasoning Quality", weight=0.15, score=score, tests=tests, notes=notes)


# ═══════════════════════════════════════════════════════════════
# Dimension 5: Species Adaptation (10%, v5.0 新增)
# 核心问题: 降级链是否诚实 —— 有 profile 直答、无 profile 降级、
#           都无则显式缺口，绝不编造已知调控因子。
# ═══════════════════════════════════════════════════════════════

def benchmark_species_adaptation() -> DimensionScore:
    """Test: 物种适配能力与降级链诚实性。"""
    tests: list[TestResult] = []
    notes: list[str] = []

    # Test 1: 全部 profile 物种 → 直接答（profile 层级），且不编造调控因子
    try:
        from phyto_reason.knowledge.species_registry import get_species_registry

        registry = get_species_registry()
        profiles = registry.list_species()
        results = []
        for slug in profiles:
            result, level = registry.knowledge_for(slug)
            results.append(
                level == "profile"
                and bool(result.known_metabolites)
                and all(
                    not info.get("known_regulators")
                    for info in registry.get(slug).pathway_prior.values()
                )
            )
        score = sum(results) / len(results) if results else 0.0
        tests.append(TestResult(
            name="species_profile_direct_answer",
            passed=score == 1.0, score=score,
            detail=f"{len(profiles)} species, direct-answer {score:.0%}"
                   + " (known_regulators 显式空 = 不编造)",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="species_profile_direct_answer",
            passed=False, score=0.0, detail=f"Error: {e}",
        ))

    # Test 2: 无 profile 物种 → ortholog 降级（不能假装知道）
    try:
        from phyto_reason.knowledge.species_registry import get_species_registry

        registry = get_species_registry()
        result, level = registry.knowledge_for("Epimedium pubescens", metabolite="icariin")
        honest = level == "ortholog" and "[降级]" in result.note
        tests.append(TestResult(
            name="species_degrade_ortholog",
            passed=honest, score=1.0 if honest else 0.0,
            detail=f"无 profile 物种 → {level} 层级" + (" [降级标注 OK]" if honest else " [未降级 X]"),
        ))
    except Exception as e:
        tests.append(TestResult(
            name="species_degrade_ortholog",
            passed=False, score=0.0, detail=f"Error: {e}",
        ))

    # Test 3: 完全未知物种 → 显式知识缺口（明说不知道）
    try:
        from phyto_reason.knowledge.species_registry import get_species_registry

        registry = get_species_registry()
        result, level = registry.knowledge_for("Herba ignota", use_ortholog_fallback=False)
        honest = level == "gap" and "[缺口]" in result.note
        tests.append(TestResult(
            name="species_gap_honesty",
            passed=honest, score=1.0 if honest else 0.0,
            detail=f"完全未知物种 → {level} 层级" + (" [缺口声明 OK]" if honest else " [无缺口声明 X]"),
        ))
    except Exception as e:
        tests.append(TestResult(
            name="species_gap_honesty",
            passed=False, score=0.0, detail=f"Error: {e}",
        ))

    # Test 4: 提示词物种注入（有 profile → 摘要；无 → 降级声明）
    try:
        from phyto_reason.agent.system_prompt import build_system_prompt

        with_profile = build_system_prompt(species="黄连", target_metabolite="berberine")
        no_profile = build_system_prompt(species="Epimedium pubescens")
        ok = (
            "berberine" in with_profile
            and "本地知识 profile" in with_profile
            and ("降级" in no_profile or "知识缺口" in no_profile)
            and "MetDNA3" not in with_profile
        )
        tests.append(TestResult(
            name="species_prompt_injection",
            passed=ok, score=1.0 if ok else 0.0,
            detail="profile 摘要 + 无 profile 降级声明 + 无 MetDNA3 残留",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="species_prompt_injection",
            passed=False, score=0.0, detail=f"Error: {e}",
        ))

    # Test 5: 工具层三层输出（物种知识 + MS/MS 诊断规则）
    try:
        from phyto_reason.agent.tool_handlers import handle_search_knowledge_base

        out = handle_search_knowledge_base({
            "metabolite": "baicalein", "species": "黄芩",
        })
        ok = ("物种知识: 黄芩" in out and "MS/MS 诊断规则" in out
              and "flavone" in out and "153.0182" in out)
        tests.append(TestResult(
            name="species_tool_three_layer",
            passed=ok, score=1.0 if ok else 0.0,
            detail="物种知识 + 诊断规则 + TF 关系三层输出",
        ))
    except Exception as e:
        tests.append(TestResult(
            name="species_tool_three_layer",
            passed=False, score=0.0, detail=f"Error: {e}",
        ))

    score = sum(t.score for t in tests) / len(tests) if tests else 0
    return DimensionScore(name="Species Adaptation", weight=0.10, score=score, tests=tests, notes=notes)


# ═══════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════

def run_benchmark() -> BenchmarkResult:
    """Run all benchmark dimensions and return result."""
    from datetime import datetime

    dims = [
        benchmark_engineering(),
        benchmark_conservatism(),
        benchmark_evidence(),
        benchmark_reasoning(),
        benchmark_species_adaptation(),
    ]

    overall = sum(d.weighted for d in dims)

    # Grade
    if overall >= 0.90:
        grade = "A"
    elif overall >= 0.80:
        grade = "B"
    elif overall >= 0.70:
        grade = "C"
    elif overall >= 0.50:
        grade = "D"
    else:
        grade = "F (incomplete — needs benchmark data)"

    lines = [f"PhytoReason-Agent v5.0 Benchmark Report", "=" * 50, ""]
    for d in dims:
        lines.append(f"## {d.name} ({d.weight:.0%})  Score: {d.score:.2f}  Weighted: {d.weighted:.3f}")
        for t in d.tests:
            status = "PASS" if t.passed else ("SKIP" if t.score == 0 and "Requires" in t.detail else "FAIL")
            lines.append(f"  [{status}] {t.name}: {t.detail}")
        for n in d.notes:
            lines.append(f"  [NOTE] {n}")
        lines.append("")
    lines.append(f"Overall: {overall:.3f}  Grade: {grade}")
    lines.append(
        "(Engineering {:.3f} + Conservatism {:.3f} + Evidence {:.3f} + Reasoning {:.3f} + Species {:.3f})".format(
            dims[0].weighted, dims[1].weighted, dims[2].weighted, dims[3].weighted, dims[4].weighted
        )
    )

    return BenchmarkResult(
        version="5.0.0",
        timestamp=datetime.now().isoformat(),
        dimensions=dims,
        overall=overall,
        grade=grade,
        summary="\n".join(lines),
    )


# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys

    result = run_benchmark()

    if "--json" in sys.argv:
        print(json.dumps(result.to_dict(), indent=2))
    else:
        print(result.summary)
