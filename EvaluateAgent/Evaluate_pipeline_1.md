# 🧬 Evaluate_pipeline_1.py — 通用代谢反应填充脚本

一个用于 **批量为 SBML 模型添加代谢反应** 的通用脚本。

本工具支持从 CSV 文件读取反应方程，为指定目录下的所有 SBML 模型自动补充反应，并可选自动补齐相关的 **EX/TRANS** 反应与代谢物。

---

## 🚀 功能简介

- 从 CSV 文件批量导入任意数量的代谢反应；
- 自动创建缺失的代谢物（按 `_c`, `_e`, `_m`, `_p` 等后缀推断舱室，默认 `c`）；
- 自动补齐 EX/TRANS 反应：
  - `TRANS_<base>` （胞外 ↔ 胞质）
  - `EX_<base>_e` （与外界交换）
- 自动补充 `o2` 的 `TRANS_o2` 与 `EX_o2_e`；
- 递归处理目录下所有 `.xml` / `.sbml` 模型；
- 输出修改后的 SBML 到指定目录。

---

## 🧩 输入文件格式

### 🔹 reaction.csv

| id | Reaction equation | lb | ub | name |
|----|--------------------|----|----|------|
| RXN_001 | `dbp_c + h2o_c -> phthalate_c + 2 btoh_c` | 0 | 1000 | DBP hydrolysis |
| RXN_002 | `o2_c + glc__D_c <-> 2 pyr_c + 2 h2o_c` | -1000 | 1000 | Glycolysis |

**列说明：**

| 列名 | 是否必需 | 说明 |
|------|-----------|------|
| `Reaction equation` / `equation` | ✅ | 反应方程，支持 `->`（不可逆）或 `<->`（可逆） |
| `id` / `Reaction id` | ❌ | 反应 ID（缺省自动生成 RXN_00001…） |
| `lb`, `ub` | ❌ | 下界与上界；若缺省，自动根据箭头类型设置：`->` → (0, 1000)，`<->` → (-1000, 1000) |
| `name` | ❌ | 反应的可读名称 |

---

## ⚙️ 命令行参数

```bash
python Evaluate_pipeline_1.py \
  --model-dir <模型目录> \
  --reactions <反应表 CSV> \
  --out <输出目录> \
  [--only <部分文件名,逗号分隔>] \
  [--auto-extrans]
```

### 参数说明

| 参数 | 说明 |
|------|------|
| `--model-dir` | 输入模型目录，递归搜索 `.xml` 或 `.sbml` 文件 |
| `--reactions` | 包含代谢反应定义的 CSV 文件 |
| `--out` | 输出目录，写入修改后的 SBML |
| `--only` | （可选）仅处理指定文件，逗号分隔 |
| `--auto-extrans` | （可选）自动为反应中出现的代谢物生成 EX 与 TRANS 反应 |

---

## 🧪 示例

### 基本用法（不自动补 EX/TRANS）

```bash
python Evaluate_pipeline_1.py \
  --model-dir "/Users/frank_gong/文档/生物智能体/硬盘备份/添加路径_功能微生物" \
  --reactions "/Users/frank_gong/文档/生物智能体/硬盘备份/信息表/reaction.csv" \
  --out "/Users/frank_gong/文档/生物智能体/硬盘备份/已添加路径"
```

### 自动补 EX/TRANS 反应(推荐✨)

```bash
python Evaluate_pipeline_1.py \
  --model-dir "/Users/frank_gong/文档/生物智能体/硬盘备份/添加路径_功能微生物" \
  --reactions "/Users/frank_gong/文档/生物智能体/硬盘备份/信息表/reaction.csv" \
  --out "/Users/frank_gong/文档/生物智能体/硬盘备份/已添加路径" \
  --auto-extrans
```

### 仅处理部分模型

```bash
python Evaluate_pipeline_1.py \
  --model-dir "./models" \
  --reactions "./reaction.csv" \
  --out "./output" \
  --only "A.xml,B.xml"
```

---

## 📦 输出结果

- 输出目录结构与原模型目录一致；
- 文件名保持不变；
- 每个模型会打印日志，列出新增的反应；
- 若反应已存在，不会重复添加。

---

## 🔍 工作逻辑简述

1. 读取 CSV 并解析每个方程；
2. 根据箭头符号判断方向性与默认上下界；
3. 为每个涉及的代谢物创建对象；
4. 构建 Reaction 并加入模型；
5. 若开启 `--auto-extrans`：
   - 自动补 `EX_<base>_e` 和 `TRANS_<base>`；
   - 若存在 `o2`，同时补 `EX_o2_e` 与 `TRANS_o2`；
6. 输出新的 SBML 模型文件。

---

## 📘 依赖

- Python ≥ 3.9  
- [cobra](https://pypi.org/project/cobra/) ≥ 0.27  
- pandas ≥ 1.5  

安装依赖：
```bash
pip install cobra pandas
```

---

## 🧠 提示与建议

- 若要统一添加多步代谢途径，可在 CSV 中按顺序列出多条方程；
- 若部分代谢物缺 compartment 后缀（如 `_c/_e`），默认归入胞质；
- 对模型规模较大的情况建议逐步测试单个反应以确认无误。

---

## 🧩 作者与许可

**Author:** Frank Gong  
**License:** MIT  

欢迎自由使用与修改。如在科研中使用，请注明引用。
