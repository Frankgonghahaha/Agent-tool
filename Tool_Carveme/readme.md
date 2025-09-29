本工具的作用是将下载下来的微生物全基因组格式转换为代谢模型。所经历步骤为先通过prodigal将全基因组.fa文件编译出可转录的氨基酸文件.aa格式。之后通过Carveme将.aa构建为基因组规模代谢模型.xml格式。
本工具提供两种工作模式分别对应于输入.aa格式文件和.fasta/fa/fna全基因组格式文件。
# build_gsmm_from_aa.py
批量化将微生物 **氨基酸序列文件 (.aa/.faa)** 转换为 **基因组规模代谢模型 (GSMM, SBML .xml)** 的 Python 脚本。  
支持两种模式：
1. 直接从已有 `.aa/.faa` 文件构建模型  
2. 从基因组 `.fa/.fna` 文件先运行 Prodigal，再进入 CarveMe 流程

---

## 功能
- **Prodigal**：从基因组 `.fa/.fna` 文件预测编码区并导出蛋白 `.aa`  
- **CarveMe**：从 `.aa/.faa` 文件构建代谢模型（.xml）  
- **批处理支持**：可多线程同时处理多个输入文件  
- **参数透传**：可将 CarveMe 的参数（如 `--fbc2`, `--init M9`）透传  

---

## 使用说明：
1）输入.aa格式文件

python build_gsmm_from_aa.py \
  --input_path /path/to/aa_folder \
  --output_path /path/to/models_out \
  --carve_extra --fbc2 --init M9
  
2）输入.fasta/fa/fna全基因组文件

python build_gsmm_from_aa.py \
  --genomes_path /path/to/genomes \
  --input_path  /path/to/aa_folder \
  --output_path /path/to/models_out \
  --carve_extra --fbc2 --init M9
  
---
## 参数说明
--input_path
必选，包含 .aa/.faa 文件的目录
--output_path
必选，输出 .xml 模型的目录
--genomes_path
可选，包含 .fa/.fna 基因组文件的目录（启用时会先跑 Prodigal）
--threads
并行线程数（默认 4）
--overwrite
若目标 .xml 已存在，是否覆盖
--carve_extra
透传给 CarveMe 的额外参数
--prodigal_mode
Prodigal 的模式，meta（宏基因组）或 single（单菌基因组），默认 meta
