# compound_profiles/ — 代谢物类 MS/MS 诊断规则

每个代谢物类一个 YAML（`alkaloids` / `flavonoids` / `terpenoids` / `saponins` /
`phenolic_acids`），由 `phyto_reason.knowledge.compound_profiles` 加载
（v5.0 起为 MS/MS 诊断离子规则的唯一规范入口，旧 chemical_expert 已移除）。

## Schema

```yaml
# <class_file>.yaml
flavonoids:                      # 顶层键 = 文件名（类库名）
  - class: flavonol              # 结构类名（get_class_info 精确匹配用）
    formula: C15H10O7            # 代表性母核分子式（可选）
    mass: 302.0427               # 代表性单同位素质量 Da（可选）
    fragments:                   # 诊断碎片离子
      - { mz: 303.0499, formula: "C15H11O7", annotation: "quercetin [M+H]+" }
    adducts: ["[M+H]+", "[M-H]-"]
    known_metabolites: [...]     # 该类典型代谢物（名称→规则检索用）
    pathway: 生源通路
    confidence: high             # high | medium | low

neutral_losses:                  # 类通用中性丢失（可选）
  - { delta: 162.0528, annotation: "C6H10O5 (hexose)" }
```

## 检索 API

```python
from phyto_reason.knowledge.compound_profiles import (
    load_compound_profile,       # 加载整个类库 dict
    list_compound_classes,       # ['alkaloids', 'flavonoids', ...]
    find_diagnostic_rules,       # 代谢物名 → 类条目（known_metabolites 双向子串匹配）
    get_class_info,              # 类名或代谢物名 → 类条目
)

load_compound_profile("flavonoids")            # 全部黄酮类规则
find_diagnostic_rules("quercetin")             # → flavonol 条目（含诊断碎片/加合物）
find_diagnostic_rules("xanthohumol")           # → None（未收录）
```

## 已收录类库

| 类库 | 类条目 | 典型代谢物 |
|---|---|---|
| alkaloids | 11（BIA/原小檗碱/单萜吲哚/托烷/喹诺里西啶/甾体糖苷碱…） | berberine, caffeine, atropine, solanine… |
| flavonoids | 6（黄酮/黄酮醇/二氢黄酮/异黄酮/黄烷-3-醇/花青素） | quercetin, baicalein, naringenin, genistein… |
| terpenoids | 6（环烯醚萜/倍半萜内酯/半日花烷/二萜醌/五环三萜/达玛烷） | artemisinin, tanshinone IIA, oleanolic acid… |
| saponins | 3（齐墩果烷/达玛烷/甾体皂苷） | glycyrrhizin, ginsenoside Rg1, dioscin… |
| phenolic_acids | 4（羟基肉桂酸/羟基苯甲酸/绿原酸/缩酚酸寡聚体） | caffeic acid, chlorogenic acid, salvianolic acid B… |
