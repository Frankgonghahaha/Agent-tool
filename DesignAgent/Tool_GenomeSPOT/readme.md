# 工具简介
本工具的功能是计算微生物的生境条件（温度、pH、盐度和是否耐氧）。本工具以微生物蛋白序列为输入，适应生境条件为输出。所使用工具为GenomeSPOT，github链接为：https://github.com/cultivarium/GenomeSPOT/blob/main/README.md
# GenomeSPOT 批量运行脚本

本脚本提供一个 **最小化命令行接口**，用于批量运行
[GenomeSPOT](https://github.com/genome-spot) 的预测任务，对多个 `.faa`
文件进行处理，并汇总所有结果生成一个汇总表格。

------------------------------------------------------------------------

## 🧬 功能概述

`run_genomespot_batch.py` 可以自动化执行 GenomeSPOT
对蛋白序列文件（`.faa`）的预测，自动匹配对应的 contig
文件，并将每个样本的预测结果汇总为一份 CSV 表格。

------------------------------------------------------------------------

## ⚙️ 功能特点

-   自动检测 GenomeSPOT 的 `models/` 目录（可从包内或本地副本获取）\
-   支持输入单个 `.faa` 文件或包含多个 `.faa` 的目录\
-   可自动匹配同名的 contigs 文件（后缀 `.fna`, `.fa`, `.fasta` 等）\
-   自动将每个样本的 `*.predictions.tsv` 输出文件汇总成一个 CSV\
-   自动标准化环境参数（如温度、pH、盐度、氧耐受性等）字段名称

------------------------------------------------------------------------

## 🧩 运行依赖

-    Requirements:

  Python version >=3.8.16 and <3.12
  ```
  git clone https://github.com/cultivarium/GenomeSPOT.git
  cd GenomeSPOT
  pip install .
  pip install -r requirements.txt
  ```
------------------------------------------------------------------------

## 🚀 使用方法

### 示例 1 --- `.faa` 与 contigs 在同一文件夹内

``` bash
python run_genomespot_batch.py     --input "/path/to/faa_dir"     --workdir outputs_genomespot
```

### 示例 2 --- `.faa` 与 contigs 分别放在不同目录中

``` bash
python run_genomespot_batch.py     --input "/path/to/faa_dir"     --contigs-dir "/path/to/contigs_dir"     --workdir outputs_genomespot
```

### 示例 3 --- 指定模型目录与汇总结果路径

``` bash
python run_genomespot_batch.py     --input "./faa_files"     --workdir "./outputs"     --models "./models"     --summary "./outputs/summary_predictions.csv"
```

------------------------------------------------------------------------

## 📁 输出文件结构

运行后会生成以下文件：

    outputs_genomespot/
    ├── sample1.predictions.tsv
    ├── sample2.predictions.tsv
    └── summary_predictions.csv

其中 `summary_predictions.csv` 是所有样本预测结果的汇总文件，示例如下：

  样本名    最适温度(°C)   最适pH   最适盐度   氧耐受性
  --------- -------------- -------- ---------- ----------
  Sample1   37.5           7.2      2.5        tolerant

------------------------------------------------------------------------

## ⚠️ 注意事项

-   请确保 GenomeSPOT 已正确安装，或其 `models/` 文件夹可用。\
-   若找不到对应的 contigs 文件（与 `.faa` 同名），该样本会被跳过。\
-   如果某个 TSV 文件列名不符合预期（如缺少
    `target`、`value`），程序会输出警告。

------------------------------------------------------------------------

## 🧠 汇总逻辑说明

脚本会从每个样本的 `*.predictions.tsv` 文件中提取并标准化参数名称：

  原始字段名     规范化字段名
  -------------- -----------------------
  Temp_optimum   temperature_optimum_C
  pH_optimum     ph_optimum
  NaCl_optimum   salinity_optimum

------------------------------------------------------------------------

## 📜 输出汇总字段

最终的汇总 CSV 文件包括以下字段：

-   strain（样本名）\
-   temperature_optimum_C / temperature_minimum / temperature_maximum\
-   ph_optimum / ph_minimum / ph_maximum\
-   salinity_optimum / salinity_minimum / salinity_maximum\
-   oxygen_tolerance

------------------------------------------------------------------------

## 🧾 许可证

MIT License © 2025 GenomeSPOT 项目贡献者

  
