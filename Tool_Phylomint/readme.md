# PhyloMInt Pipeline Wrapper

本项目提供了一个 Python 脚本 `run_phylomint.py`，用于简化 **PhyloMInt** 的调用流程，并在运行结束后自动处理结果表格。  
支持自动运行 **Step1–4（运行 PhyloMInt、结果转换、表头清理）** 和 **Step5（物种互补分析，结果输出为 Excel）**。  

---

## 📌 功能特性

### ✅ Step 1–4: PhyloMInt 执行与结果处理
1. 调用 PhyloMInt 运行代谢模型分析  
2. 将输出 TSV 文件转换为 CSV  
3. 修改 CSV 表头：  
   - 第一列 → `A`  
   - 第二列 → `B`  
   - 第三列 → `Competition`  
   - 第四列 → `Complementarity`  
   - 删除第 5 列  
4. 去掉 A、B 列中的 `_CDS` 后缀  

### ✅ Step 5: 物种互补分析（可选）
- 从 `--function-species-csv` 指定的物种列表文件中读取目标物种  
- 在 PhyloMInt 结果表中查询目标物种的互补信息  
- 计算 **Metabolic Complementarity Strength = Complementarity - Competition**  
- 仅保留 >0 的记录  
- 输出到 **Excel 文件**：  
  - `Summary` sheet（总表，包含所有物种结果，按物种和互补强度排序）  
  - 每个物种一个独立 sheet  

---

## ⚙️ 依赖环境

- Python >= 3.8  
- 必须安装的库：
  ```bash
  pip install pandas xlsxwriter
  ```
- 需提前部署并可执行的 [PhyloMInt](https://github.com/)  

---

## 🚀 使用方法

### 场景 1：完整流程（Step1–5）
```bash
python run_phylomint.py   --phylo /work/.../PhyloMint/PhyloMInt   --models /work/.../代谢模型   --output /work/.../PhyloMInt_output.csv   --function-species-csv /work/.../species_list.csv
```

### 场景 2：已有处理好的 CSV，直接进入 Step5
```bash
python run_phylomint.py   --output /work/.../processed_Phylomint.csv   --function-species-csv /work/.../species_list.csv   --skip-preprocess
```

---

## 📂 输入文件说明

### 1. PhyloMInt 输出结果
脚本会自动调用 PhyloMInt 生成 TSV，并转换为清理后的 CSV。

### 2. species_list.csv（可选）
- 格式：一列物种名（有无列名都可）  
- 示例：
  ```csv
  species_name
  Escherichia_coli
  Bacillus_subtilis
  Pseudomonas_aeruginosa
  ```

---

## 📊 输出文件说明

### 1. 清理后的 CSV
文件名由 `--output` 参数指定，例如：
```
PhyloMInt_output.csv
```

### 2. 互补微生物识别结果（Step5）
如果启用了 `--function-species-csv`，将生成：
```
互补微生物识别结果.xlsx
```

#### 📑 Summary sheet
- 汇总所有物种的互补信息  
- 列结构：
  - `Target Species`
  - `Complementarity Species`
  - `Competition Index`
  - `Complementarity Index`
  - `Metabolic Complementarity Strength`

#### 📑 物种独立 sheet
- 每个物种一个工作表  
- 列结构：
  - `Complementarity Species`
  - `Competition Index`
  - `Complementarity Index`
  - `Metabolic Complementarity Strength`

---

## ⚠️ 注意事项

- Excel sheet 名称限制为 **31 字符**，超过会自动截断  
- 请确保 `species_list.csv` 中的物种名与结果表 **B 列**一致（已去掉 `_CDS` 后缀）  
- 结果仅保留 **Metabolic Complementarity Strength > 0** 的记录  
- Summary 表中结果 **按 Target Species 升序、互补强度降序** 排列  
