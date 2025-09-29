本工具的作用是将下载下来的微生物全基因组格式转换为代谢模型。所经历步骤为先通过prodigal将全基因组.fa文件编译出可转录的氨基酸文件.aa格式。之后通过Carveme将.aa构建为基因组规模代谢模型.xml格式。



使用方式：
python build_gsmm_from_aa.py \
  --input_path /path/to/aa_folder \
  --output_path /path/to/models_out \
  --threads 6 \
  --carve_extra --fbc2
