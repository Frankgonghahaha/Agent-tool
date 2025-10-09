
# 工具简介
- 该工具主要用于预测降解酶对于特定底物的降解速率，所使用工具为 **DLKcat**，github 链接为：https://github.com/SysBioChalmers/DLKcat

# DLKcat 预测批处理脚本（参数化版）

该脚本用于：
1) 从 Excel 中读取蛋白序列与底物信息并生成 **`input.tsv`**；
2) 调用 **`prediction_for_input.py`** 完成 Kcat 预测；
3) 将 **`Kcat value (1/s)`** 写回 Excel，默认输出为 **`Sheet3_Function_enzyme_kact.xlsx`**。

---
## 安装
  ```
  git clone https://github.com/SysBioChalmers/DLKcat.git
  ```
## ✅ 运行环境
- Python 3.8+
- 已安装 PyTorch（脚本会预检），如需 CPU 版：
  ```bash
  python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
  ```
- 依赖包：
  ```bash
  pip install pandas openpyxl
  ```

---

## 📦 输入文件格式
输入的 Excel 须包含以下列（区分大小写）：
- `Sequence` ：降解酶的蛋白序列
- `substrate_name`：底物的名称，需要是 PubChem 的标准化名称
- `substrate_smiles`：底物的 SMILES 结构，需要是 PubChem 的结果

### 🧮 输入文件示例
| Enzyme   | Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | substrate_name    | substrate_smiles                |
|:---------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------|:--------------------------------|
| patE     | MSALTAAAEEYQRLRTEFREKGLGGRIGFGVRPAVVVVDLITGFTDRRSPLAGDLDTQIDATKILLALARKAQVPIIFSTVAYDAELQEAGAWIGKIPSNKYLVEGSQWVEIDERLEQQPGETTLVKKYASCFFGTDLAARLISRRIDTVIIVGCTTSGCVRATAVDACSYGFHTIVVEDAVGDRAALPHTASLFDIDAKYGDVVGLDEASAYLESVPSSS                                                                                                                                                                                                                                                                                                                              | Dibutyl Phthalate | O=C(C1=CC=CC=C1C(=O)OCCCC)OCCCC |
| pehA     | MEIVLVHGGWVGGWVWDGVADELRRMGHEVIAPTLRGLEDGDVDRSGVTMSMMARDLIDQVRELTQLDIVLVGHSGGGPLIQLVAEAMPERIGRVVFVDAWVLRDGETINDVLPDPLVAATKALASQSDDNTIVMPPELWAASMQDMSPFEQQQLAALEPRLVPSPAGWSDEPIRLDRFWASSIPSSYVFLAQDQAVPAEIYQAAAGRLDSPRTIEIDGSHLVMLTHPERLARALDAVIA                                                                                                                                                                                                                                                                                                           | Dibutyl Phthalate | O=C(C1=CC=CC=C1C(=O)OCCCC)OCCCC |
| GoEst15  | MGVPNTAVADTRVSTAEGIVSGVRGRRSRRGTVAWRGIPYAAPPVGGGRFDAPQPVAPWPGVRRCENFGDAAVQDKLLTATGLGRFQPVSEDCLTLNVVAPDTVASAPRPVMVFIHGGAYILGTAATPLHDGTHLARAQDVVVVTIQYRFGPFGMLDFSGYSTPERHFDENPGLKDHIAALQWVQNIAAFGGDPDNVTIFGESAGGTSVLCLLAAPGARGLFHRAIAESPATYLAISRESAALFADEFLRLLADPSRRSTDATSPREPIDPHEAARRLDAATPEELHRAGGRLMGFARHADSVEPAPFGPAYGVPSLPQSPYAAARDGNTAPVPLIIGTNRDEGKLFDKFWNLLPDTQRTLLAIQDETARDEIMNQYAGGPDDQLRLTADSIFWAPVTAFADAHREVAPTYVYRFDFHTRALARAGLGATHGTELFAVFGGYRMPVFAGLATADWRAAGAMVDEMQSRWGDFARGRPPGPDWPAYDAAHRPVMVLDRTSHVESDPDASRRQAWDRGRQLLVASGGSPEAPSDVAALDV | Dibutyl Phthalate | O=C(C1=CC=CC=C1C(=O)OCCCC)OCCCC |

---

## 📊 输出文件示例
| Enzyme   | Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | substrate_name    | substrate_smiles                |   Kcat value (1/s) |
|:---------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------|:--------------------------------|-------------------:|
| patE     | MSALTAAAEEYQRLRTEFREKGLGGRIGFGVRPAVVVVDLITGFTDRRSPLAGDLDTQIDATKILLALARKAQVPIIFSTVAYDAELQEAGAWIGKIPSNKYLVEGSQWVEIDERLEQQPGETTLVKKYASCFFGTDLAARLISRRIDTVIIVGCTTSGCVRATAVDACSYGFHTIVVEDAVGDRAALPHTASLFDIDAKYGDVVGLDEASAYLESVPSSS                                                                                                                                                                                                                                                                                                                              | Dibutyl Phthalate | O=C(C1=CC=CC=C1C(=O)OCCCC)OCCCC |             0.0029 |
| pehA     | MEIVLVHGGWVGGWVWDGVADELRRMGHEVIAPTLRGLEDGDVDRSGVTMSMMARDLIDQVRELTQLDIVLVGHSGGGPLIQLVAEAMPERIGRVVFVDAWVLRDGETINDVLPDPLVAATKALASQSDDNTIVMPPELWAASMQDMSPFEQQQLAALEPRLVPSPAGWSDEPIRLDRFWASSIPSSYVFLAQDQAVPAEIYQAAAGRLDSPRTIEIDGSHLVMLTHPERLARALDAVIA                                                                                                                                                                                                                                                                                                           | Dibutyl Phthalate | O=C(C1=CC=CC=C1C(=O)OCCCC)OCCCC |             0.0004 |
| GoEst15  | MGVPNTAVADTRVSTAEGIVSGVRGRRSRRGTVAWRGIPYAAPPVGGGRFDAPQPVAPWPGVRRCENFGDAAVQDKLLTATGLGRFQPVSEDCLTLNVVAPDTVASAPRPVMVFIHGGAYILGTAATPLHDGTHLARAQDVVVVTIQYRFGPFGMLDFSGYSTPERHFDENPGLKDHIAALQWVQNIAAFGGDPDNVTIFGESAGGTSVLCLLAAPGARGLFHRAIAESPATYLAISRESAALFADEFLRLLADPSRRSTDATSPREPIDPHEAARRLDAATPEELHRAGGRLMGFARHADSVEPAPFGPAYGVPSLPQSPYAAARDGNTAPVPLIIGTNRDEGKLFDKFWNLLPDTQRTLLAIQDETARDEIMNQYAGGPDDQLRLTADSIFWAPVTAFADAHREVAPTYVYRFDFHTRALARAGLGATHGTELFAVFGGYRMPVFAGLATADWRAAGAMVDEMQSRWGDFARGRPPGPDWPAYDAAHRPVMVLDRTSHVESDPDASRRQAWDRGRQLLVASGGSPEAPSDVAALDV | Dibutyl Phthalate | O=C(C1=CC=CC=C1C(=O)OCCCC)OCCCC |             0.016  |

> 输出文件相比输入多出一列 **`Kcat value (1/s)`**，为预测的酶学常数结果。

---

## 🧭 命令行参数

| 参数 | 是否必填 | 说明 |
|---|---|---|
| `--script-path` | ✅ | `prediction_for_input.py` 的完整路径 |
| `--file-path` | ✅ | 输入 Excel 文件路径 |
| `--output-path` | ❌ | 输出 Excel 路径（默认与输入同目录，文件名为 `Sheet3_Function_enzyme_kact.xlsx`） |

> 运行时脚本会在 **输出 Excel 同目录** 生成 `input.tsv`，并在常见位置自动探测 `output.tsv`。

---

## 🚀 使用示例

### 1. 使用默认输出文件名（与输入同目录）
```bash
python DLKcat.py   --script-path "/path/to/prediction_for_input.py"   --file-path "/path/to/降解功能酶序列.xlsx"
```

### 2. 指定输出 Excel 完整路径
```bash
python DLKcat.py   --script-path "/path/to/prediction_for_input.py"   --file-path "/path/to/降解功能酶序列.xlsx"   --output-path "/path/to/results/Sheet3_Function_enzyme_kact.xlsx"
```

---

## 📂 生成的文件
- `input.tsv`：给 `prediction_for_input.py` 的输入文件（位于输出 Excel 同目录）
- `output.tsv`：预测脚本生成（脚本会在常见位置自动查找）
- `Sheet3_Function_enzyme_kact.xlsx`：包含追加的 `Kcat value (1/s)` 列

---

## ❗ 常见问题
- **报错：未检测到 PyTorch**：请先安装与本机环境匹配的 `torch`。CPU 版安装命令见上。
- **找不到 `output.tsv`**：请检查 `prediction_for_input.py` 实际写入路径，或把其输出目录调整为与 `input.tsv` 同目录。
