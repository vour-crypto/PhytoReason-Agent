# PhytoReason-Agent v5.0 Benchmark 设计方案

> **用途**：本文件是 benchmark 数据收集的操作指南。
> 目标：不测"准确率"（调控关系是概率性的），而是测**保守性、可追溯性、工程健全性**。
> 运行：`python -m phyto_reason.benchmark.scoring`（已实现，cases/ 已就绪）

---

## 1. 评估维度与权重

| 维度 | 权重 | 核心问题 |
|------|------|----------|
| **工程健全性** | 25% | pipeline 不崩溃、可复现、不超时？ |
| **保守性** | 25% | 数据不足时系统能否拒绝输出高置信度结论？ |
| **证据质量** | 25% | 知识层召回（B01/B02/B03/B06）？KEGG 映射是否准确？ |
| **推理质量** | 15% | 机制分类是否合理？矛盾是否被检测？ |
| **物种适配**（v5.0 新增） | 10% | 降级链是否诚实（profile 直答 / ortholog 降级 / gap 显式缺口，不编造）？ |

**打分原则：保守性 > 证据质量 > 召回率。** 宁可保守的系统比"什么都敢说"的更有科学价值。

---

## 2. 数据收集清单

### B01: 已知 TF-靶基因调控对（正向案例）

| 属性 | 值 |
|------|-----|
| 用途 | 验证证据采集能力：已知关系是否出现在候选集中 |
| 数量 | **50+ 对** |
| 来源 | TRANSFAC, JASPAR, PlantTFDB, 已发表 ChIP-seq 论文 |
| 优先级 | 🔴 最高 |

**格式：**
```yaml
- id: POS-001
  tf_gene: "AtMYB75"
  tf_family: "MYB"
  target_gene: "DFR"
  species: "Arabidopsis thaliana"
  metabolite: "anthocyanin"
  evidence_type: "ChIP-seq + mutant phenotype"
  pmid: "25228336"
  confidence: "confirmed"  # experimentally validated
```

**建议数据源：**
1. [PlantTFDB](http://planttfdb.gao-lab.org/) — TF 家族注释
2. [JASPAR](https://jaspar.genereg.net/) — TF binding profiles
3. 搜 `"ChIP-seq" AND "secondary metabolism" AND "Arabidopsis"` 在 PubMed
4. 综述文章：DOI 10.1146/annurev-arplant-050213-035642


### B02: 已知 Pathway 基因集（通路映射验证）

| 属性 | 值 |
|------|-----|
| 用途 | 验证 KEGG/内置通路映射的正确性 |
| 数量 | **20 条通路** |
| 来源 | KEGG, PlantCyc |
| 优先级 | 🟡 高 |

**格式：**
```yaml
- id: PWY-001
  pathway_id: "map00950"
  pathway_name: "Isoquinoline alkaloid biosynthesis"
  expected_enzymes: ["BBE", "CYP80G2", "SMT", "CNMT", "6OMT", "4'OMT"]
  species_example: "Coptis japonica"
  metabolite: "berberine"
```


### B03: 文献查询标准答案（Layer 0 评估）

| 属性 | 值 |
|------|-----|
| 用途 | 评估文献检索和知识库查询的覆盖度 |
| 数量 | **20 条查询** |
| 来源 | 人工标注（你知道答案的问题） |
| 优先级 | 🔴 最高 |

**格式：**
```yaml
- id: L0-001
  query: "What TFs are known to regulate berberine biosynthesis?"
  query_zh: "已知调控小檗碱生物合成的转录因子有哪些？"
  expected_tf_families: ["bHLH", "WRKY", "ERF"]
  expected_pathway: "isoquinoline alkaloid biosynthesis"
  min_expected_pmids: 3
  evaluation: "manual"  # 需要人工判断回答质量
  acceptable_terms: ["bHLH", "WRKY", "ERF", "CjbHLH1", "isoquinoline"]

- id: L0-002
  query: "Does MYB regulate anthocyanin in plants?"
  query_zh: "MYB 转录因子是否调控植物花青素合成？"
  expected_tf_families: ["MYB", "bHLH", "WD40"]
  expected_mechanism: "MBW complex"
  min_expected_pmids: 5

- id: L0-003
  query: "What is the biosynthetic pathway of nicotine?"
  query_zh: "尼古丁的生物合成途径是什么？"
  expected_pathway: "pyridine alkaloid biosynthesis"
  expected_enzymes: ["PMT", "MPO", "QPT", "A622"]
  min_expected_pmids: 3
```

**建议的 20 条查询覆盖：**
- 5 条：特定代谢物调控（berberine, nicotine, anthocyanin, artemisinin, taxol）
- 5 条：TF 家族功能（MYB, bHLH, WRKY, ERF, NAC 在次生代谢中的作用）
- 5 条：跨物种推断（Arabidopsis→tobacco, rice→maize, etc.）
- 3 条：通路查询
- 2 条：实验设计建议


### B04: 阴性对照数据（保守性测试）

| 属性 | 值 |
|------|-----|
| 用途 | 验证系统在无信号数据上的保守性 |
| 数量 | **10 套** |
| 来源 | **程序生成** |
| 优先级 | 🟡 高 |

**生成方法：**
```python
# 随机打乱表达矩阵
import numpy as np
import pandas as pd

# 加载真实数据
df = pd.read_csv("real_expression.csv", index_col=0)

# 生成 10 套阴性对照
for i in range(10):
    shuffled = df.apply(lambda col: np.random.permutation(col.values), axis=0)
    shuffled.to_csv(f"negative_control_{i:02d}.csv")
```

**期望行为：**
- 系统应拒绝给出高置信度结论
- `data_quality_gate` 可能通过（数据格式正确），但后续节点应降级
- 最终输出应包含"insufficient evidence"或"weak"


### B05: 小样本边界案例

| 属性 | 值 |
|------|-----|
| 用途 | 验证系统在边界条件下的行为 |
| 数量 | **5 套**（从正常数据中取子集） |
| 来源 | 程序生成 |
| 优先级 | 🟢 中 |

**生成方法：**
- 从完整数据（≥12 样本）中随机取 3、4、6、8、10 个样本
- 期望：样本 < 6 时系统应警告"样本不足，相关性不稳定"


### B06: 物种-代谢物映射对（跨物种推理验证）

| 属性 | 值 |
|------|-----|
| 用途 | 验证跨物种 TF→代谢物的推断正确性 |
| 数量 | **30 对** |
| 来源 | NCBI Taxonomy + 文献 |
| 优先级 | 🟢 中 |

**格式：**
```yaml
- id: CS-001
  source_species: "Arabidopsis thaliana"
  target_species: "Coptis chinensis"
  tf_family: "bHLH"
  metabolite: "berberine"
  known_conservation: "plausible"  # 有证据表明可能保守
  note: "CjbHLH1 identified in Coptis, but mechanism different from Arabidopsis"
  pmid: "25296257"
```

---

## 3. 目录结构

```
benchmark/
├── README.md                    # 本文件
├── scoring.py                   # 评分脚本（待实现）
├── cases/
│   ├── positive_tf_pairs.yaml   # B01: 已知调控对
│   ├── pathways.yaml            # B02: 通路基因集
│   ├── layer0_queries.yaml      # B03: 文献查询标准答案
│   ├── cross_species.yaml       # B06: 跨物种映射
│   └── boundary_cases.yaml      # B05: 边界案例定义
├── data/
│   ├── negative_controls/       # B04: 程序生成的阴性对照
│   └── small_samples/           # B05: 小样本子集
└── results/
    └── (评分输出存放)
```

---

## 4. 评分脚本设计

```python
# benchmark/scoring.py (伪代码)

@dataclass
class BenchmarkResult:
    engineering_score: float     # 0-1
    conservatism_score: float    # 0-1
    evidence_score: float        # 0-1
    reasoning_score: float       # 0-1
    overall: float               # 加权总分
    details: list[TestResult]
    failures: list[str]


class BenchmarkRunner:
    def run_all(self) -> BenchmarkResult:
        eng = self._test_engineering()    # 30%
        con = self._test_conservatism()   # 30%
        evd = self._test_evidence()       # 25%
        rea = self._test_reasoning()      # 15%
        return BenchmarkResult(...)

    def _test_engineering(self) -> float:
        # 1. 用 B01 数据跑 pipeline，验证不崩溃
        # 2. 跑两次相同输入，验证输出一致（deterministic）
        # 3. 计时：单次分析 < 120s

    def _test_conservatism(self) -> float:
        # 1. 用 B04 阴性对照数据跑 pipeline
        # 2. 期望：所有结果置信度 ≤ "weak"
        # 3. 期望：系统能识别"无信号"
        # 4. 用 B05 小样本数据测试警告机制

    def _test_evidence(self) -> float:
        # 1. 用 B01 已知 TF-靶基因对验证召回
        # 2. 用 B03 查询验证文献检索覆盖
        # 3. 用 B02 验证 KEGG 通路映射

    def _test_reasoning(self) -> float:
        # 1. 对已知通路，验证 mechanism_classifier 分类
        # 2. 验证 contradiction_check 能检测已知 artifact
        # 3. 验证 hypothesis 中包含"为什么可能错"和"缺失证据"
```

---

## 5. 快速开始：最小可用的 benchmark

如果你时间有限，**优先收集以下最小集合**：

| 编号 | 内容 | 数量 | 预计时间 |
|------|------|------|----------|
| B01 | 已知 TF-靶基因调控对 | 20 对 | 2-4 小时（搜文献） |
| B03 | Layer 0 查询 + 标准答案 | 10 条 | 1-2 小时（你自己想 10 个问题） |
| B04 | 阴性对照数据 | 3 套 | 10 分钟（程序生成） |

这 3 样足够跑第一轮 benchmark，得到有意义的评分。

---

## 6. 数据收集后如何跑 benchmark

```bash
# 1. 把 B01 数据放到 benchmark/data/positive/
# 2. 把 B03 查询写成 benchmark/cases/layer0_queries.yaml
# 3. 生成阴性对照
python benchmark/generate_negatives.py

# 4. 跑 benchmark
python benchmark/run.py --all

# 5. 查看报告
cat benchmark/results/report.md
```

---

## 7. 期望的 benchmark 基线

| 维度 | 最低可接受 | 良好 | 目标 |
|------|-----------|------|------|
| 工程健全性 | 0.80 | 0.90 | 0.95 |
| 保守性 | 0.70 | 0.85 | 0.90 |
| 证据质量 | 0.50 | 0.70 | 0.80 |
| 推理质量 | 0.50 | 0.65 | 0.75 |

保守性和工程健全性是最重要的 —— 宁可系统说 "我不知道" 而不是胡说。
