import argparse
import os
import sys
import subprocess
import pandas as pd
from pathlib import Path

# 预检：确保当前解释器已安装 PyTorch（子进程将复用该解释器）
try:
    import torch  # noqa: F401
except Exception:
    raise SystemExit(
        "未检测到 PyTorch。请在当前解释器中安装：\n"
        f"  {sys.executable} -m pip install torch --index-url https://download.pytorch.org/whl/cpu\n"
        "或使用对应 CUDA 版本的安装命令。"
    )

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "批量生成 DLKcat 所需 input.tsv，调用 prediction_for_input.py 进行预测，"
            "并把 Kcat 结果写回到 Excel。"
        )
    )
    parser.add_argument(
        "--script-path",
        required=True,
        help="prediction_for_input.py 的完整路径（将以当前解释器调用）",
    )
    parser.add_argument(
        "--file-path",
        required=True,
        help="输入 Excel 文件路径，需包含列：Sequence、substrate_name、substrate_smiles",
    )
    parser.add_argument(
        "--output-path",
        default=None,
        help=(
            "输出 Excel 路径（默认为输入文件所在目录下的 'Sheet3_Function_enzyme_kact.xlsx'）。"
        ),
    )
    return parser.parse_args()

def main():
    args = parse_args()

    script_path = Path(args.script_path).expanduser().resolve()
    input_excel = Path(args.file_path).expanduser().resolve()

    if args.output_path:
        output_excel = Path(args.output_path).expanduser().resolve()
    else:
        output_excel = input_excel.with_name("Sheet3_Function_enzyme_kact.xlsx")

    if not script_path.exists():
        raise FileNotFoundError(f"找不到预测脚本: {script_path}")
    if not input_excel.exists():
        raise FileNotFoundError(f"找不到输入 Excel: {input_excel}")

    # 读取 Excel
    data = pd.read_excel(input_excel)

    # 校验必要列
    required_cols = {"Sequence", "substrate_name", "substrate_smiles"}
    missing = required_cols - set(map(str, data.columns))
    if missing:
        raise KeyError(
            "输入 Excel 缺少必要列: " + ", ".join(sorted(missing))
        )

    sequences = data["Sequence"]
    names = data["substrate_name"]
    smiles = data["substrate_smiles"]

    # 生成 input.tsv：默认写到输出 Excel 同目录，避免权限/路径混乱
    input_tsv = output_excel.with_name("input.tsv")
    input_df = pd.DataFrame({
        "Substrate Name": names,
        "Substrate SMILES": smiles,
        "Protein Sequence": sequences,
    })
    input_df.to_csv(input_tsv, sep="\t", index=False)
    print(f"[INFO] 写入预测输入: {input_tsv}")

    # 在预测脚本所在目录执行，保证其相对导入（如 import model）可用
    script_dir = script_path.parent
    predict_cmd = [sys.executable, str(script_path), str(input_tsv)]
    print(f"[INFO] 运行预测脚本: {' '.join(predict_cmd)}\n      cwd={script_dir}")

    subprocess.run(predict_cmd, cwd=str(script_dir), check=True)
    print("[INFO] 预测完成。")

    # 探测 output.tsv 的位置
    candidate_outputs = [
        # 与输入同目录
        input_tsv.with_name("output.tsv"),
        # 预测脚本目录
        script_dir / "output.tsv",
        # 常见的用户自定义路径（与原项目类似）
        input_excel.with_name("output.tsv"),
    ]
    output_tsv = next((p for p in candidate_outputs if p.exists()), None)

    if output_tsv is None:
        raise FileNotFoundError(
            "未找到 output.tsv。请检查 prediction_for_input.py 实际输出路径，"
            "或在脚本中补充 candidate_outputs。"
        )

    print(f"[INFO] 读取预测结果: {output_tsv}")
    pred_df = pd.read_csv(output_tsv, sep="\t")

    # 将 Kcat 列追加回原表
    if "Kcat value (1/s)" not in pred_df.columns:
        raise KeyError("预测结果缺少列: 'Kcat value (1/s)'")

    data["Kcat value (1/s)"] = pred_df["Kcat value (1/s)"]

    # 写出 Excel，文件名按需已更改为 Sheet3_Function_enzyme_kact.xlsx
    output_excel.parent.mkdir(parents=True, exist_ok=True)
    data.to_excel(output_excel, index=False)
    print(f"[INFO] 已写出结果 Excel: {output_excel}")

if __name__ == "__main__":
    main()
