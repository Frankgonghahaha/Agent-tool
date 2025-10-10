
# Evalaute_pipeline_1  
该代码主要是确定 **pipeline 2** 群落模拟的培养基成分及设置通量上限值。

---

## ✨ 功能

- 递归遍历某目录下所有 `.xml/.sbml` 模型；
- 使用 [MICOM](https://github.com/micom-dev/micom) 的 `complete_*` 工作流（优先 `workflows.media.complete_community_medium`；若不可用则回退到 `media.complete_medium`）；
- 对每个模型打印**非零** `EX_*_m` 进口（全部）；
- 汇总所有模型，按 `p75_flux` 分档生成**推荐培养基**：
  - `<0.01 → 0.01`，`[0.01,0.1) → 0.1`，`[0.1,1) → 1`，`(1,10) → 10`，其余 `clip` 到 `--max-import`
- 导出为 CSV（两列：`reaction, flux`）。
---

## 📦 依赖

- `micom`, `cobra`
-
- 安装示例：
```bash
pip install micom cobra pandas
```

> 某些求解器（如 `glpk`, `cplex`, `gurobi`）可能需要单独安装与许可证。`micom` 默认使用可用的后端。

---

## 🚀 使用

```bash
python pipeline1_cli.py \
  --model-dir "/path/to/models" \
  --out "/path/to/Result5_Medium.csv" \
  --community-growth 0.1 \
  --min-growth 0.1 \
  --max-import 20 \
  --candidate-ex EX_glc__D_m=10 EX_o2_m=20
```

### 参数说明

- `--model-dir`：**必填**。包含 SBML 模型的目录（递归查找 `.xml/.sbml`）。
- `--out`：**必填**。输出 CSV 路径（两列：`reaction, flux`）。
- `--community-growth`：社区生长率参数（默认 `0.1`）。
- `--min-growth`：最小生长率参数（默认 `0.1`）。
- `--max-import`：单一底物最大进口上界（默认 `20`）。
- `--candidate-ex`：候选 EX 列表，形如 `EX_id=value`，多项用空格分隔。默认：`EX_glc__D_m=10 EX_o2_m=20`。

---

## 🧠 推荐培养基的生成逻辑

1. 对每个模型调用 MICOM 的 `complete_*`，得到该模型在指定参数下**需要进口**的 `EX_*_m` 通量；
2. 过滤掉通量≤0 的项；
3. 按反应名聚合，统计 `mean/p50/p75/max`；
4. 用 `p75` 做分档映射（`0.01/0.1/1/10`），并 `clip` 到 `--max-import`；
5. 导出为 `reaction, flux`。

> 注意：此“推荐培养基”仅包含**净摄入**项；如需在验证阶段允许排出（例如 CO₂ / H₂O / H⁺ 等平衡出口），请在后续步骤中自行扩展白名单。

---

## 📄 输出

- `--out` 指定路径的 CSV，列：
  - `reaction`：社区层交换反应（如 `EX_glc__D_m`）；
  - `flux`：推荐上界（已分档）。

---

## 🔎 日志示例

- 每个模型的“非零 `EX_*_m`（全部）”
- “基于单菌扫描得到的社区培养基建议”预览（前 50 行）
- “已导出推荐培养基：/path/to/Result5_Medium.csv”





