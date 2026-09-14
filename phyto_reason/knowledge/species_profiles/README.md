# species_profiles/ — 物种知识配置

**新增物种 = 新增一个 YAML**（文件名即 slug，如 `coptis_chinensis.yaml`），
无需改任何代码。由 `phyto_reason.knowledge.species_registry.SpeciesRegistry`
加载，查询走三级降级链（有 profile → 直接答；无 → 同源迁移降级；都无 → 显式知识缺口）。

## Schema

```yaml
# <slug>.yaml
species:
  scientific_name: Zanthoxylum nitidum   # 学名（拉丁双名）
  common_name: 两面针                    # 俗名/中文名（可选）
  taxonomy_id: 210110                    # NCBI Taxonomy ID（可选，用于 ortholog_finder）
  genome_status: unavailable             # available | partial | unavailable
                                         #   available → 可做 motif/TF 注释等基因组依赖分析

marker_metabolites:                      # 标志代谢物（targeted confirmation 的输入锚点）
  - name: nitidine                       # 代谢物名（小写英文，匹配 KEGG/HMDB 命名）
    formula: C21H18NO4                   # 中性分子式（可选）
    mz: 348.1230                         # 单同位素 m/z（注释中写明加合物）
    adducts: ["[M]+"]                    # 常见电离加合物
    fragments:                           # 诊断碎片离子（可选）
      - { mz: 332.0923, formula: "C20H14NO4", annotation: "[M-CH3]+" }
    reference_rt: null                   # 参考保留时间（可选）
    source: 药典                         # 数据来源：药典 / 文献 / 推测

pathway_prior:                           # 该物种已知通路先验
  benzylisoquinoline_alkaloid_biosynthesis:
    genes: []                            # 已知通路基因（可选）
    enzymes: ["CYP719", "BBE"]           # 已知酶
    known_regulators: []                 # 已知调控因子——留空 = 显式"未知"（不编造）

refs:                                    # 知识来源（药典/文献/数据库）
  - "中国药典 2020 版"
```

## 降级链（SpeciesRegistry.knowledge_for）

| 层级 | 条件 | 行为 | 置信度 |
|---|---|---|---|
| `profile` | 有本地 YAML | 直接答：标志代谢物 + 通路先验 + 同源推断补充 | 正常 |
| `ortholog` | 无 YAML | `SpeciesKnowledge.infer_regulation` 同源迁移 | 降级 + caveats |
| `gap` | 都无 | 显式缺口消息，明说"不知道" | — |

## 示例

```python
from phyto_reason.knowledge.species_registry import get_species_registry

reg = get_species_registry()
reg.list_species()                    # ['zanthoxylum_nitidum']
reg.get("两面针")                      # SpeciesProfile（模糊匹配）
result, level = reg.knowledge_for("Zanthoxylum nitidum", metabolite="nitidine")
# level == "profile"
```

## 已收录物种

| slug | 物种 | 标志代谢物 | 通路 | 基因组 |
|---|---|---|---|---|
| `zanthoxylum_nitidum` | 两面针 | nitidine / chelerythrine / sanguinarine | BIA → 苯并菲啶 | unavailable |
| `coptis_chinensis` | 黄连 | berberine / coptisine / palmatine / jatrorrhizine | BIA → 原小檗碱 | unavailable |
| `scutellaria_baicalensis` | 黄芩 | baicalein / baicalin / wogonin / wogonoside | 黄酮（flavone） | available |
| `salvia_miltiorrhiza` | 丹参 | tanshinone IIA / cryptotanshinone / tanshinone I / 丹酚酸 B | 二萜醌 + 酚酸 | available |
| `glycyrrhiza_uralensis` | 甘草 | glycyrrhizin / 甘草次酸 / liquiritigenin / liquiritin | 三萜皂苷 + 黄酮 | partial |
