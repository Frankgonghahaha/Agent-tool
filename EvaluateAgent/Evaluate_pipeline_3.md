# 🧪 Evaluate_pipeline_3.py — 社区在推荐培养基上的两步优化（DBP 摄入，MICOM medium 耦合）

本脚本用于在**推荐培养基**下，对给定目录中的所有单菌模型构建**群落（MICOM Community）**，并进行**两阶段优化**：

1) **阶段一：**最大化社区增长率（Biomass_max）  
2) **阶段二：**在约束 `growth ≥ α · Biomass_max` 下，**最小化 `EX_dbp_m` 通量**（负值越小代表摄入越强），得到社区对 DBP 的最大总摄入能力。

脚本使用 **MICOM 的 `community.medium`** 设置培养基，从而在**社区层与成员层**建立交换通量的**守恒耦合**。

---

## ✨ 主要特性

- 递归读取 `--model-dir` 下的全部 `.xml/.sbml` 并构建社区  
- 从 `--medium-csv` 读取培养基（支持列 `reaction + flux` 或 `reaction + suggested_upper_bound`）  
- **自动加入 `EX_dbp_m = 20`** 到培养基上限（可在代码中改为其他值）  
- 使用 `community.medium = {...}` 应用培养基并建立耦合  
- 两阶段优化；可选 **α 扫描**；可选 **鲁棒性分析（Robust）**：逐一移除群落成员  
- 输出 Excel：  
  - `Simulate result`：总览（成员数、群落组成、阶段二 Biomass 与 DBP flux）  
  - `Microbial growth`：阶段二下各成员的增长率  
  - `Robust`（可选）：移除单个物种后的 Biomass 与 DBP flux

---

## 📥 输入

### 1) 模型目录（`--model-dir`）
- 递归搜索 `.xml/.sbml`，每个文件作为一个成员  
- 构建社区时会**统一外液舱室**标记（将 `ext/external/C_e` 等同义词并到 `e`）  
- 脚本会尝试自动识别并设置各成员的**biomass 反应为 objective**

### 2) 培养基 CSV（`--medium-csv`）
- 必须包含列 `reaction`，另一个列如下二选一：
  - `flux`（表示上界）
  - `suggested_upper_bound`（表示上界）
- 上界会被写入 `community.medium` 中（正值=最大摄入）。脚本内还会**强制**加入 `EX_dbp_m = 20`。

---

## ▶️ 用法

```bash
python Evaluate_pipeline_3.py \
  --model-dir "/Users/frank_gong/文档/生物智能体/硬盘备份/Candidate_models" \
  --medium-csv "/Users/frank_gong/文档/生物智能体/硬盘备份/信息表/Result5_Medium_components.csv"
```

可选参数：
- `--alpha 0.7`：阶段二的增长率下界系数（默认 0.7）
- `--out <输出文件>`：输出路径（默认写到 `model-dir/Result6_community.xlsx`）
- `--alpha-scan`：在一组 α 上批量执行阶段二（默认扫描 0.50, 0.60, 0.70, 0.80, 0.90）
- `--alphas "0.55,0.65,0.75"`：自定义 α 列表（配合 `--alpha-scan`）
- `--Robust`：启用鲁棒性分析（依次移除一个成员重跑两阶段）

示例（含 α 扫描与鲁棒性）：

```bash
python Evaluate_pipeline_3.py \
  --model-dir "/Users/frank_gong/文档/生物智能体/硬盘备份/Candidate_models" \
  --medium-csv "/Users/frank_gong/文档/生物智能体/硬盘备份/信息表/Result5_Medium_components.csv" \
  --alpha 0.7 \
  --alpha-scan \
  --Robust \
  --out "/Users/frank_gong/文档/生物智能体/硬盘备份/Candidate_models/Result6_community.xlsx"
```

---

## 📦 输出

生成一个 Excel（默认 `Result6_community.xlsx`），包含：

1. **Simulate result**  
   - `model_count`：成员数量  
   - `Communitity`：群落组成（成员名称，逗号分隔）  
   - `Biomass`：阶段二的群落增长率（钳制 ≤ Biomass_max，避免数值抖动）  
   - `DBP flux`：阶段二搜索到的 `EX_dbp_m` 固定通量（负值=摄入）

2. **Microbial growth**  
   - `member`：成员 ID  
   - `growth_rate`：阶段二条件下的成员增长率  
   > 注：脚本内会过滤掉可能出现的 `medium` 汇总行。

3. **Robust**（当使用 `--Robust` 时出现）  
   - `Removal species`：被移除的成员  
   - `DBP flux`：该移除情形下阶段二的 `EX_dbp_m`  
   - `Biomass`：该移除情形下阶段二的群落增长率

如启用 `--alpha-scan`，还会导出 `Result6_community_alpha_scan.csv`（每个 α 的 `(growth_at_f, dbp_flux_f)`）。

---

## 🛠️ 常见问题

- **阶段一 Biomass_max = 0**  
  - 培养基过于严格或关键 EX 缺失；确认 `community.medium` 中存在必要的 `EX_*` 并且上界为正  
  - 单菌模型目标/质量问题（biomass 反应未识别或不可行）

- **阶段二输出的 growth 略高于 Biomass_max**  
  - 数值抖动：已在输出做“钳制”，仅展示 `min(growth(f), Biomass_max + 1e-6)`

- **`EX_dbp_m` 不存在**  
  - 社区里没有名为 `EX_dbp_m` 的反应；脚本会提示并返回 `nan`，请在构建模型时确保该交换反应存在。

---

## 📘 依赖

- Python ≥ 3.9  
- `cobra` ≥ 0.27  
- `micom`（与 `optlang`）  
- `pandas`

安装示例：
```bash
pip install cobra micom pandas
```

---

## 📜 License
MIT
