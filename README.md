# PhytoReason-Agent

药用植物次生代谢的科研推理 Agent。把组学分析流水线和 LLM 推理拆成两层:确定性计算(差异分析、聚类、网络)交给代码,不确定性的事情(证据解释、候选调控因子排序)交给 LLM,两边不越界。


## 为什么有这个项目

多数药用植物没有基因组,公共数据稀疏。拿这类物种直接问 LLM"某某成分受什么调控",它会给你编一个看起来很合理的答案。另一方面,实验室里多组学分析的流程(质控、差异分析、WGCNA、网络构建)高度重复,每次都是同一套代码换个物种再跑一遍。

这个项目把两件事合成一个系统:分析流程沉淀为可复用的工具节点,LLM 负责调度工具、解释证据、组织假设——并且**在被问到自己不知道的物种时,明说不知道**,而不是编。

## 现在能做什么

- **带降级的知识问答**。内置 6 个物种知识档案(黄芩、黄连、丹参、甘草、烟草、两面针)。问档案内物种,直接回答;问没有档案的物种,走同源迁移并压低置信度;两边都没有,输出显式知识缺口。
- **多组学分析流水线**。上传表达矩阵或代谢物丰度表,依次执行质控、差异分析(DEG/DAM)、WGCNA、转录-代谢关联网络,产出候选调控因子及证据链。
- **MS/MS 谱图注释**。MGF/CSV/MSP 上传后给出分子式、结构类推断和 MSI 分级,支持物种锚定;对苄基异喹啉类生物碱和类黄酮内置了诊断离子规则库。
- **转录因子预测**。上传蛋白序列 FASTA,pyhmmer 跑 HMM,注释自动注入后续会话分析。
- **三个入口**:命令行(`plantomics`)、HTTP 服务(`uvicorn`)、PySide6 桌面客户端。

## 快速开始

Python 3.10+。

```bash
git clone https://github.com/vour-crypto/PhytoReason-Agent.git
cd PhytoReason-Agent
pip install -e ".[desktop,ms2,tfpred]"
```

配置 LLM 后端(项目根目录建 `.env`,兼容 OpenAI 接口的服务都可以):

```
OPENAI_API_KEY=sk-xxx
OPENAI_BASE_URL=https://xxx/v1
LLM_MODEL=xxx
```

命令行试一下:

```bash
plantomics onboarding                 # 首次引导,检查环境
plantomics ask "黄芩素在黄芩里的已知调控因子有哪些？" --species 黄芩 --metabolite baicalein
plantomics species list               # 看看现在有哪些物种档案
plantomics predict_tf protein.fasta   # 转录因子预测
plantomics annotate spectra.msp       # MS/MS 注释
```

起本地服务:

```bash
uvicorn phyto_reason.api.app:app --reload
```

桌面客户端:

```bash
plantomics-desktop
```

## 几个设计决定

**确定性计算和 LLM 推理分层。** 早期版本试过把统计分析也交给 LLM,结论是它连基础统计都会出错,而且错了还理直气壮。v5 起所有数值计算收归流水线节点,LLM 只看到计算结果,负责解释和推理。这条边界是整个系统里最重要的约束。

**"不知道"是一等公民。** 物种知识查询走三级降级:有档案直答 → 无档案走直系同源迁移(置信度降级、附加 caveats)→ 迁移也没有就返回显式知识缺口。系统输出里的每条结论都带证据来源,Plausible / Weak / Insufficient 三档,不给"X 调控 Y"这种断言式表述。

**物种知识是数据不是代码。** 每个物种一个 YAML(标记代谢物、通路先验、已知调控因子,允许留空)。加物种不改代码。

## 目录结构

```
phyto_reason/
├── agent/            # AgentOrchestrator,会话与工具调度
├── workflows/        # 20 节点分析 DAG(质控→差异→网络→假设)
├── knowledge/
│   ├── species_profiles/    # 物种档案 YAML
│   └── compound_profiles/   # 代谢物类 MS/MS 诊断规则
├── tools/            # 22 个工具:检索/统计/谱图/预测
├── reasoning/ fusion/  # 证据融合、假设证伪与竞争排序
├── api/  desktop/  agents/   # HTTP 服务 / PySide6 客户端 / CLI
└── tests/            # 30 个测试模块
```

## 现状与局限

说实话的部分:

- 物种档案只有 6 个,化合物规则库只覆盖生物碱和类黄酮两类,大部分物种问过来还是走降级或缺口。
- benchmark 框架是自建的,题目量还小,评分口径主观成分不低,数字仅供参考。
- 桌面端是单机应用,API 没做鉴权,**不要**把服务直接暴露到公网。
- MS/MS 注释依赖规则库覆盖,超出覆盖范围的化合物仍需人工解谱,这个工具只是初筛。
- 所有输出都是"候选假设",不是生物学结论。验证要回实验室。

## Roadmap

- **2026 Q4**:PyInstaller 打包 Windows 安装包;Docker 镜像;API 加最简单的鉴权。
- **2026 Q4**:benchmark 扩容并把评分口径写进文档。
- **视立项**:育种方向(水稻)的知识接入,和合作方谈。
- **按需**:新物种档案、新化合物类规则库。

## 数据说明

`phyto_reason/tests/fixtures/` 下的矩阵是合成数据,仅用于测试。`massbank_sample.json` 是 MassBank 公开数据的**截断流样例**——保留截断是故意的,`parse_massbank_records` 需要处理真实下载中流被掐断的情况。仓库不含任何未发表的课题数据。

本项目的分析代码用于一篇撰写中的研究论文(两面针多组学,投稿准备中)。

## License

MIT
