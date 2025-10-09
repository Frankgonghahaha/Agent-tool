# 📖 智能体简介
  本智能体主要作用是根据工程微生物识别智能体的识别结果和用户输入的水质参数（温度、pH、盐度和好氧/缺氧/厌氧）进行菌剂的设计，最后输出出最佳的设计结果。
---
# 📃 准备文件
  - Sheet1_Complementarity.xlsx : 候选名单里物种之前的代谢互补指数和代谢竞争指数
  - Sheet2_Species_environment.xlsx：Tool_GenomeSPOT 计算所得微生物生境条件
  - Sheet3_Function_enzyme_kact.xlsx: Tool_DLkcat 计算所得酶降解速率
  - Sheet4_species_enzyme.xlsx：降解功能微生物包含的降解酶的种类
  - Sheet5 Sheet5_PhyloMInt.csv： 微生物数据库总的PhyloMInt表
---
# 🚀 计算过程
**Step 1 运行Design_pipeline_1.py**
---
  ## 工作内容
  `Design_pipline_1.py` 主要完成 3 件事：
  1. **功能菌候选筛选**：按目标工况（温度、pH、盐度、氧环境）筛选满足条件的功能菌，并汇总每个菌的 `kcat_max / kcat_mean`；  
  2. **单菌综合打分（S_microbe）**：对功能菌与通过环境筛选的互补菌计算单菌得分，输出到 `Result2_candidate_scores.csv`；  
  3. **两两互作矩阵**：融合 Sheet1（专家/实验知识）与 PhyloMint（基因组距离推断）生成 `(competition, complementarity, delta)` 指标，输出到 `Result3_pair_Com_index.csv`。
---
  ## 🧭 运行方式
  脚本以**交互式**方式获取目标工况：
  ```text
  目标温度 (°C)：
  目标 pH：
  目标盐度（% NaCl）：
  氧环境（好氧 / 厌氧 / 缺氧）：
  ```
  - 氧环境会自动映射为：好氧 → `tolerant`；厌氧/缺氧 → `not tolerant`。
---
## **功能菌候选筛选**
根据输入的目标工况，匹配适应的功能微生物和互补微生物，并汇总匹配物种、含有酶种类、酶降解速率等信息，输入表Result1_candidate_function_species.csv
---
## **单菌综合打分（S_microbe）**
根据计算公式对单菌进行环境适应度和降解功能能力方面的综合打分，打分计算公式为：
```
S_microbe
= 0.5 * Norm01(kcat_max)
+ 0.4 * env_soft_score
+ 0.1 * Norm01(enzyme_diversity)
```
- `Norm01` 为列内 min–max 归一化（整列常数→1）；kcat_max 是所含降解酶最大降解速率。
- `env_soft_score` 是该物种的环境评分具体下下一小节。
- `enzyme_diversity` 来自 `Sheet4` 中该物种的酶名列表长度。

输出列（节选）：
- `species, kcat_max, kcat_mean, enzyme_diversity, environment_match(=env_soft_score), source, S_microbe`
---

## 🧠 软环境匹配评分（0~1）
- **温度/ pH**：区间内按“**三角隶属**”随偏离最适点线性下降；超出区间按“**指数衰减**”降分；
- **盐度**：不超过上限记 1；超出上限按“**指数衰减**”；
- **氧气**：匹配=1.0，未知=0.5，不匹配=0.2；
- 自动忽略缺失项并**重归一化权重**。

默认权重（可在函数 `soft_env_score_row` 中改）：
- 温度 `wT=0.35`、pH `wPH=0.35`、盐度 `wS=0.10`、氧气 `wO2=0.20`。

**结果列：**
- `env_match_all`：布尔型，是否同时满足四项硬条件；
- `env_soft_score`：0~1 的软匹配分（用于后续 `S_microbe` 的环境分量）。

---

## 🧪 功能菌候选与 Kcat 汇总
- 从 `Sheet4_species_enzyme.xlsx` 获取物种-酶映射；
- 在 `Sheet3_Function_enzyme_kact.xlsx` 中按酶名聚合其 `Kcat value (1/s)` 的**中位数**；
- 每个物种：`kcat_max` = 其酶的 Kcat 中位数最大值；`kcat_mean` = 其酶的 Kcat 中位数的平均；
- 使用目标工况筛掉 `env_match_all=False` 的功能菌；保留的候选表按 `kcat_max` 降序输出：  
  **`Result1_candidate_function_species.csv`**（脚本里常量 `OUTPUT_PATH`）。

候选表附带互补菌通过/失败统计：
- `complement_total / complement_pass / complement_fail`；
- `complement_pass_names / complement_fail_names`（以 `;` 分隔）。

---

## 🔢 单菌综合打分 S_microbe
对 **功能菌 + 通过软评分的互补菌** 输出 `Result2_candidate_scores.csv`。评分构成：
```
S_microbe
= 0.5 * Norm01(kcat_max)
+ 0.4 * env_soft_score
+ 0.1 * Norm01(enzyme_diversity)
```
- `Norm01` 为列内 min–max 归一化（整列常数→1）；
- `enzyme_diversity` 来自 `Sheet4` 中该物种的酶名列表长度。

输出列（节选）：
- `species, kcat_max, kcat_mean, enzyme_diversity, environment_match(=env_soft_score), source, S_microbe`

---

## 🔁 两两互作矩阵（Sheet1 + PhyloMint）
- 若 `(功能菌, 对应互补菌)` 在 **Sheet1** 中找到，则**优先**使用 Sheet1 的 `competition, complementarity, delta`；
- 否则尝试在 **PhyloMint** 中用基因组配对的平均值；
- 结果对称补齐（便于任意方向检索）；输出到：  
  **`Result3_pair_Com_index.csv`**。

主要列：
- `functional_species, complement_species, competition_index, complementarity_index, delta_index, source`

---

## 📤 脚本输出一览
- `Result1_candidate_function_species.csv`：满足目标工况的功能菌候选（含 `kcat_*` 与互补菌统计）；
- `Result2_candidate_scores.csv`：功能菌+互补菌的单菌综合打分；
- `Result3_pair_Com_index.csv`：两两互作矩阵（融合 Sheet1 与 PhyloMint）。

---

## 🧰 常见问题
- **列名不匹配**：脚本对四张表的列名/工作表名有假设，若与你的数据不同，请在“固定列映射”部分调整。  
- **环境字段缺失**：软评分会自动略过缺失维度并重归一化；但 `env_match_all` 的硬条件需要字段完整才会严格判定。  
- **PhyloMint 缺失**：脚本会提示并仅用 Sheet1 的配对关系生成矩阵。
---

# 🚀 计算过程
**Step 2 运行Design_pipeline_2.py**
---
## 📌 工具简介
- 目标：在给定候选物种集合与两两互作指标的基础上，搜索得到**综合得分最高**的菌群组合。
- 评分考虑：单菌打分、两两互作（增益与竞争）、成员平均 `kcat_max`、组合规模惩罚。

---
## 📥 输入文件格式
以pipeline1输出的结果表Result2_candidate_scores.csv和Result3_pair_Com_index.csv为输入。

## 🧮 目标函数（组合得分）

对成员集合 $M$（规模 $|M|=N$）的最终得分：

$
S_{\text{consort}}
= \alpha\, \overline{S_{\text{microbe}}}
+ \beta\, \overline{\Delta^+}
- \gamma\, \overline{\mathrm{Comp}^+}
+ \lambda\, \overline{\mathrm{kcat}}
- \mu\, N


- $\overline{S_{\text{microbe}}}$：成员在 `scores` 表中 `S_microbe` 的**算术平均**。
- $\overline{\Delta^+}$：成员两两 `delta_index` 中 **>0** 的平均。
- $\overline{\mathrm{Comp}^+}$：成员两两 `competition_index` 中 **>0** 的平均（惩罚项）。
- $\overline{\mathrm{kcat}}$（即 `avg_kcat`）：成员在 `scores` 表中 `kcat_max` 的**算术平均**（缺失的成员将被忽略于平均）。
- $N$：组合规模（成员数）。



---

## ⚙️ 命令行参数与默认值

| 参数 | 说明 | 默认值 | 必填 |
|---|---|---|:---:|
| `--scores` | `Result2_candidate_scores.csv` 路径 | 无 | ✅ |
| `--pairs` | `Result3_pair_Com_index.csv` 路径 | 无 | ✅ |
| `--out` | 最优组合输出 | `Result4_optimal_consortia.csv` |  |
| `--members_out` | 成员贡献排名输出 | `Result4_members_rank.csv` |  |
| `--topN` | 先按 `S_microbe` 取前 N 个入搜索池 | `50` |  |
| `--kmin` / `--kmax` | 组合最小/最大规模 | `2` / `5` |  |
| `--mode` | `greedy` 或 `exhaustive` | `greedy` |  |
| `--topK` | 输出前 K 个最优组合 | `20` |  |
| `--require_functional` | 是否要求含至少 1 个功能菌（1/0） | `1` |  |
| `--alpha` | 权重：平均 `S_microbe` | `0.2` |  |
| `--beta` | 权重：平均正 `delta_index` | `0.2` |  |
| `--gamma` | 权重：平均正 `competition_index`（惩罚） | `0.1` |  |
| `--lambda_` | 权重：平均 `kcat_max`（加分） | `0.35` |  |
| `--mu` | 权重：规模惩罚 | `-0.05` |  |

---

## 🚀 运行示例

### 1) 贪心搜索（推荐快速试跑）
```bash
python Design_pipeline_2.py   --scores ./Result2_candidate_scores.csv   --pairs  ./Result3_pair_Com_index.csv   --kmin 2 --kmax 6 --mode greedy --topK 30
```

### 2) 穷举搜索（小规模/严格）
```bash
python Design_pipeline_2.py   --scores ./Result2_candidate_scores.csv   --pairs  ./Result3_pair_Com_index.csv   --topN 30 --kmin 2 --kmax 4 --mode exhaustive --topK 50
```
---

### 📊 输出示例
- 脚本输出不同组合的得分结果：`Result4_optimal_consortia.csv`
以下为示例数据（前 5 行）：

| consortium_id | members | size | avg_S_microbe | avg_delta_pos | avg_comp_pos | avg_kcat | S_consort | used_pairs | source_count |
|----------------|----------|------|----------------|----------------|---------------|------------|-------------|---------------|----------------|
| C001 | Bacillus subtilis; Pseudomonas putida; Rhodococcus erythropolis | 3 | 0.84 | 0.21 | 0.08 | 35.2 | 0.563 | 3 | {"functional":1,"complement":2} |
| C002 | Pseudomonas putida; E.coli K12; Acinetobacter baumannii | 3 | 0.78 | 0.19 | 0.10 | 28.7 | 0.524 | 3 | {"functional":1,"complement":2} |
| C003 | Bacillus subtilis; Acinetobacter baumannii | 2 | 0.81 | 0.23 | 0.05 | 31.5 | 0.518 | 1 | {"functional":1,"complement":1} |
| C004 | Rhodococcus erythropolis; Pseudomonas putida | 2 | 0.79 | 0.18 | 0.07 | 29.4 | 0.503 | 1 | {"functional":1,"complement":1} |
| C005 | Bacillus subtilis; E.coli K12; Streptomyces coelicolor | 3 | 0.82 | 0.16 | 0.09 | 33.8 | 0.499 | 3 | {"functional":1,"complement":2} |

> 注：
> - `avg_S_microbe`：组合内成员平均的单菌打分；
> - `avg_delta_pos` / `avg_comp_pos`：正向互补性与竞争性平均；
> - `avg_kcat`：组合成员 `kcat_max` 的平均；
> - `S_consort`：综合评分结果；
> - `used_pairs`：该组合在互作矩阵中命中的成对数量；
> - `source_count`：功能菌/互补菌来源计数（JSON 结构）。



