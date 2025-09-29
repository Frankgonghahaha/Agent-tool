# GSMM Builder (Prodigal + CarveMe)

本项目提供了一个 Python 脚本，用于 **批量将氨基酸序列文件 (.aa/.faa) 或基因组文件 (.fa/.fna/fasta) 构建为基因组规模代谢模型 (GSMM, SBML .xml 格式)**。  
支持直接输入氨基酸文件，或先用 **Prodigal** 从基因组预测基因再构建。

---

## 📌 功能特性

- 批量处理 `.aa` / `.faa` → `.xml`  
- 可选从基因组 `.fa` / `.fna` / `.fasta`先调用 **Prodigal** 生成 `.aa` 文件  
- 支持并行计算（默认 4 线程，可配置）  
- 支持覆盖/跳过已有模型  
- 可透传参数给 **CarveMe**（如 `--fbc2`、`--universal-model` 等）  

---

## ⚙️ 依赖环境
- 必须安装的外部工具：
  - [CarveMe](https://github.com/cdanielmachado/carveme) （提供 `carve` 命令）  
  - [Prodigal](https://github.com/hyattpd/Prodigal) （提供 `prodigal` 命令）

---

## 🚀 使用方法

### 场景 1：已有Prodigal .aa / .faa 文件
```bash
python build_gsmm.py   --input_path /path/to/aa_files   --output_path /path/to/output_models   --threads 8   --overwrite   --carve_extra --fbc2 --universal-model bacteria
```

### 场景 2：输入基因组文件（先跑 Prodigal，再 CarveMe）
```bash
python build_gsmm.py   --input_path /path/to/aa_out   --output_path /path/to/output_models   --genomes_path /path/to/genomes   --threads 8   --prodigal_mode meta   --carve_extra --fbc2 --universal-model bacteria
```

### 常用参数说明
- `--input_path`：包含 `.aa/.faa` 的目录  
- `--output_path`：保存 `.xml` 的目录  
- `--threads`：并行线程数（默认 4）  
- `--overwrite`：覆盖已有同名 `.xml` 文件  
- `--carve_cmd`：CarveMe 命令名（默认 `carve`）  
- `--carve_extra`：透传参数给 CarveMe，例如 `--fbc2`  
- `--genomes_path`：可选，包含基因组文件的目录（.fa/.fna），会先用 Prodigal 生成 `.aa`  
- `--prodigal_cmd`：Prodigal 命令名（默认 `prodigal`）  
- `--prodigal_mode`：Prodigal 模式（`meta` 或 `single`，默认 meta）  

---

## 📂 输入输出示例

### 输入目录
```
input_path/
├── sample1.aa
├── sample2.faa
genomes_path/
├── genome1.fa
├── genome2.fna
```

### 输出目录
```
output_path/
├── sample1.xml
├── sample2.xml
├── genome1.xml
├── genome2.xml
```

---

## ⚠️ 注意事项

- 如果提供了 `--genomes_path`，会将 Prodigal 输出的 `.aa` 写入 `--input_path` 指定的目录  
- CarveMe 的执行参数可通过 `--carve_extra` 灵活传入  
- 确保 `carve` 和 `prodigal` 已经正确安装并在 PATH 中可执行  
