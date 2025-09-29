#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd


def rename_csv_columns(csv_path: Path):
    """
    Step 3:
    - 修改 CSV 文件的前四列列名
    Step 4:
    - 去掉 A、B 列的 "_CDS" 后缀
    """
    if not csv_path.exists():
        print(f"[ERROR] 找不到 CSV 文件: {csv_path}", file=sys.stderr)
        sys.exit(4)

    df = pd.read_csv(csv_path)
    if len(df.columns) < 5:
        print(f"[ERROR] CSV 列数不足 5 列", file=sys.stderr)
        sys.exit(6)

    # Step 3: 重命名前四列，删除第五列
    new_columns = ['A', 'B', 'Competition', 'Complementarity', 'TO_DROP']
    df = df.iloc[:, :5]
    df.columns = new_columns
    df = df.drop(columns=['TO_DROP'])

    # Step 4: 去掉 _CDS 后缀
    df["A"] = df["A"].astype(str).str.replace("_CDS", "", regex=False)
    df["B"] = df["B"].astype(str).str.replace("_CDS", "", regex=False)

    df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    print(f"[OK] 已修改表头并清理 _CDS 后缀: {csv_path}")
    return df


def process_species(input_csv: Path, species_csv: Path):
    """
    Step 5: 物种互补分析
    - 从 species_csv 读取物种名
    - 在结果表 B 列查询
    - 计算 Metabolic Complementarity Strength
    - 过滤 >0 的结果
    - 保存到 Excel 文件：Summary 总表（按 Target Species 排序）+ 每个物种一个 sheet
    """
    df = pd.read_csv(input_csv)

    try:
        species_list = pd.read_csv(species_csv).iloc[:, 0].dropna().astype(str).tolist()
    except Exception as e:
        print(f"[ERROR] 读取物种列表失败: {e}", file=sys.stderr)
        sys.exit(9)

    out_path = input_csv.parent / "互补微生物识别结果.xlsx"
    if out_path.exists():
        out_path.unlink()

    summary_records = []

    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        for sp in species_list:
            sub_df = df[df["B"] == sp].copy()
            if sub_df.empty:
                print(f"[WARN] 物种 {sp} 未找到记录")
                continue

            sub_df["Metabolic Complementarity Strength"] = sub_df["Complementarity"] - sub_df["Competition"]
            sub_df = sub_df[sub_df["Metabolic Complementarity Strength"] > 0]

            if sub_df.empty:
                print(f"[INFO] 物种 {sp} 没有正向互补记录")
                continue

            result_df = pd.DataFrame({
                "Complementarity Species": sub_df["A"],
                "Competition Index": sub_df["Competition"],
                "Complementarity Index": sub_df["Complementarity"],
                "Metabolic Complementarity Strength": sub_df["Metabolic Complementarity Strength"]
            })

            # 添加到 Summary 总表
            temp = result_df.copy()
            temp.insert(0, "Target Species", sp)
            summary_records.append(temp)

            # 保存单独物种表
            result_df.to_excel(writer, sheet_name=sp[:31], index=False)

            print(f"[{sp}] Complementarity Species: {', '.join(result_df['Complementarity Species'].astype(str))}")

        # 写入 Summary 总表
        if summary_records:
            summary_df = pd.concat(summary_records, ignore_index=True)
            summary_df = summary_df.sort_values(
                by=["Target Species", "Metabolic Complementarity Strength"],
                ascending=[True, False]
            )
            summary_df.to_excel(writer, sheet_name="Summary", index=False)

    print(f"[OK] 识别结果已保存到 {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="调用 PhyloMInt 并自动输出 CSV（含表头修改、后缀清理、物种分析）"
    )
    parser.add_argument("--phylo", help="PhyloMInt 可执行文件路径（仅在执行 Step1–4 时必需）")
    parser.add_argument("--models", help="模型目录路径（仅在执行 Step1–4 时必需）")
    parser.add_argument("--output", required=True, help="输出或输入的 CSV 文件完整路径（含文件名）")
    parser.add_argument("--function-species-csv", help="可选：包含物种名的 CSV 文件，启用则执行 Step5 物种分析")
    parser.add_argument("--skip-preprocess", action="store_true",
                        help="跳过 Step1–4，直接对已处理好的 CSV 运行 Step5")
    args = parser.parse_args()

    output_csv = Path(args.output).expanduser().resolve()
    outdir = output_csv.parent
    outdir.mkdir(parents=True, exist_ok=True)

    if not args.skip_preprocess:
        # Step 1–4 正常执行
        if not args.phylo or not args.models:
            print("[ERROR] 执行 Step1–4 需要同时提供 --phylo 和 --models 参数", file=sys.stderr)
            sys.exit(10)

        phylo = Path(args.phylo).expanduser().resolve()
        models = Path(args.models).expanduser().resolve()

        # 临时 TSV 文件
        temp_tsv = outdir / (output_csv.stem + ".tsv")

        # Step 1: 运行 PhyloMInt
        cmd = [str(phylo), "-d", str(models), "--outdir", str(outdir), "-o", temp_tsv.name]
        print("运行命令:", " ".join(cmd))
        subprocess.run(cmd, check=True)

        if not temp_tsv.exists():
            print(f"[ERROR] 未找到 TSV: {temp_tsv}", file=sys.stderr)
            sys.exit(2)

        # Step 2: TSV 转 CSV
        df = pd.read_csv(temp_tsv, sep="\t", header=0)
        df.to_csv(output_csv, sep=",", index=False)
        temp_tsv.unlink()
        print(f"[OK] 已生成 CSV 文件: {output_csv}")

        # Step 3 & 4: 修改表头并清理 _CDS 后缀
        rename_csv_columns(output_csv)

    else:
        print("[INFO] 跳过 Step1–4，直接进入 Step5")

    # Step 5: 如果提供了物种查询表，执行物种互补分析
    if args.function_species_csv:
        process_species(output_csv, Path(args.function_species_csv))


if __name__ == "__main__":
    main()
