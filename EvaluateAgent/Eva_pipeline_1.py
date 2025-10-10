
# -*- coding: utf-8 -*-
"""
Pipeline 1 (CLI) — 生成推荐培养基组分

功能概述：
- 遍历给定目录下的全部 SBML 模型（.xml/.sbml）
- 使用 MICOM 的 complete_*（优先 workflows，缺省回退到 media.complete_medium）
  计算每个模型在给定 growth 参数下的最小进口需求（EX_*_m 通量）
- 汇总所有模型的非零进口，按 p75_flux 分档生成“推荐培养基”并导出到指定路径
- 控制台打印：每个模型的非零 EX_*_m（全部）、推荐培养基表

注意：本脚本不做“单菌生长验证”与“缺失 EX 的补齐迭代”。
"""

import os
import argparse
import tempfile
import shutil
import warnings
from typing import List
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
from micom import Community

# ---------- 首选 workflows，回退到单模型 API ----------
HAVE_WORKFLOW = False
try:
    from micom.workflows.media import complete_community_medium as wf_complete
    HAVE_WORKFLOW = True
except Exception:
    from micom.media import complete_medium as single_complete
    HAVE_WORKFLOW = False

warnings.filterwarnings("ignore", category=FutureWarning)

EXTERNAL_COMPARTMENT_SYNONYMS = {"C_e", "ext", "external", "extracellular"}
TARGET_EXTERNAL = "e"


def discover_models(models_dir: str) -> List[str]:
    """递归发现 .xml/.sbml"""
    paths: List[str] = []
    for root, _dirs, files in os.walk(models_dir):
        for fn in files:
            if fn.lower().endswith((".xml", ".sbml")):
                paths.append(os.path.join(root, fn))
    paths.sort()
    return paths


def normalize_external_compartment(model):
    """将明确的外液同义舱室并到 'e'。"""
    changed = 0
    for met in model.metabolites:
        comp = (met.compartment or "").strip()
        if comp in EXTERNAL_COMPARTMENT_SYNONYMS:
            met.compartment = TARGET_EXTERNAL
            changed += 1
    if changed:
        print(f" · 统一外液舱室：修正 {changed} 个代谢物 compartment → 'e'")
    return model


def build_singleton_community(sbml_path: str, name: str) -> Community:
    tax = (
        pd.DataFrame({
            "id": [name],
            "file": [sbml_path],
            "abundance": [1.0],
            "biomass": [None],
        })
        .set_index("id", drop=False)
    )
    return Community(tax, name=f"COMM_{name}")


def run_with_workflow(com: Community, name: str, community_growth: float, min_growth: float, max_import: float,
                      candidate_ex: pd.DataFrame) -> pd.DataFrame:
    tmpd = tempfile.mkdtemp(prefix=f"cm_wf_{name}_")
    try:
        pickle_name = f"{name}.pickle"
        pickle_path = os.path.join(tmpd, pickle_name)
        com.to_pickle(pickle_path)
        manifest = pd.DataFrame({
            "sample_id": [name],
            "file": [pickle_name],
        })
        fixed = wf_complete(
            manifest=manifest,
            model_folder=tmpd,
            medium=candidate_ex,
            community_growth=community_growth,
            min_growth=min_growth,
            max_import=max_import,
            summarize=True,
            threads=1,
        )
        out = fixed[["reaction", "flux"]].copy()
        out.insert(0, "model", name)
        return out
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def run_with_single_api(com: Community, name: str, community_growth: float, min_growth: float, max_import: float,
                        candidate_ex: pd.DataFrame) -> pd.DataFrame:
    med_series = pd.Series({r: float(v) for r, v in candidate_ex.values})
    fixed = single_complete(
        model=com,
        medium=med_series,
        growth=community_growth,
        min_growth=min_growth,
        max_import=max_import,
    )
    df = pd.DataFrame({"reaction": fixed.index, "flux": fixed.values})
    df.insert(0, "model", name)
    return df


def run_for_model(sbml_path: str, community_growth: float, min_growth: float, max_import: float,
                  candidate_ex: pd.DataFrame) -> pd.DataFrame:
    """统一外液 → 写入临时 SBML → 构建 Community → 调用相应 API。"""
    name = os.path.splitext(os.path.basename(sbml_path))[0]
    tmpd = tempfile.mkdtemp(prefix=f"cm_sbml_{name}_")
    try:
        model = read_sbml_model(sbml_path)
        model = normalize_external_compartment(model)
        tmp_sbml = os.path.join(tmpd, f"{name}.xml")
        write_sbml_model(model, tmp_sbml)
        com = build_singleton_community(tmp_sbml, name)

        if HAVE_WORKFLOW:
            return run_with_workflow(com, name, community_growth, min_growth, max_import, candidate_ex)
        else:
            return run_with_single_api(com, name, community_growth, min_growth, max_import, candidate_ex)
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def recommend_medium(all_rows: pd.DataFrame, max_import: float) -> pd.DataFrame:
    """基于全部单菌结果生成推荐培养基（reaction, suggested_upper_bound）。"""
    df_pos = all_rows[all_rows["flux"] > 1e-12].copy()
    if df_pos.empty:
        return pd.DataFrame(columns=["reaction", "suggested_upper_bound"])

    grp = df_pos.groupby("reaction")["flux"]
    summary = pd.DataFrame({
        "models_with_need": grp.size(),
        "mean_flux": grp.mean(),
        "p50_flux": grp.quantile(0.5),
        "p75_flux": grp.quantile(0.75),
        "max_flux": grp.max(),
    }).sort_values(["models_with_need", "p75_flux", "mean_flux"], ascending=False)

    # 分档：<0.01→0.01；[0.01,0.1)→0.1；[0.1,1)→1；(1,10)→10；>=10→clip 为 max_import
    v = summary["p75_flux"].astype(float)
    sug = v.copy()
    sug = sug.where(v >= 0.01, 0.01)
    sug = sug.mask((v >= 0.01) & (v < 0.1), 0.1)
    sug = sug.mask((v >= 0.1) & (v < 1.0), 1.0)
    sug = sug.mask((v > 1.0) & (v < 10.0), 10.0)
    summary["suggested_upper_bound"] = sug.clip(upper=max_import)

    medium_df = summary.reset_index()[["reaction", "suggested_upper_bound"]].copy()
    return medium_df


def main():
    parser = argparse.ArgumentParser(description="Pipeline1 CLI：扫描单菌 EX 需求并生成推荐培养基")
    parser.add_argument("--model-dir", required=True, help="代谢模型目录（递归搜索 .xml/.sbml）")
    parser.add_argument("--out", required=True, help="推荐培养基输出 CSV 路径（reaction,flux）")
    parser.add_argument("--community-growth", type=float, default=0.1, help="community growth（默认 0.1）")
    parser.add_argument("--min-growth", type=float, default=0.1, help="min growth（默认 0.1）")
    parser.add_argument("--max-import", type=float, default=20.0, help="单一底物最大进口（默认 20）")
    parser.add_argument("--candidate-ex", nargs="*", default=["EX_glc__D_m=10", "EX_o2_m=20"],
                        help="候选 EX 列表，形如 EX_id=value，多项用空格分隔")
    args = parser.parse_args()

    # 解析候选 EX
    cand_pairs = []
    for item in args.candidate_ex:
        if "=" in item:
            rid, val = item.split("=", 1)
            try:
                cand_pairs.append((rid.strip(), float(val)))
            except Exception:
                pass
    if not cand_pairs:
        cand_pairs = [("EX_glc__D_m", 10.0), ("EX_o2_m", 20.0)]
    candidate_ex = pd.DataFrame(cand_pairs, columns=["reaction", "flux"])

    print("使用 API：", "workflows.complete_community_medium" if HAVE_WORKFLOW else "media.complete_medium  (回退模式)")
    print(f"参数：community_growth={args.community_growth}, min_growth={args.min_growth}, max_import={args.max_import}")
    print("候选培养基：", ", ".join([f"{r}={v}" for r, v in cand_pairs]))

    # 遍历模型
    model_paths = discover_models(args.model_dir)
    if not model_paths:
        print(f"[错误] 未在目录下发现 SBML：{args.model_dir}")
        return
    print(f"共发现 {len(model_paths)} 个模型，将逐一计算…")

    all_rows = []
    for fpath in model_paths:
        print("\n" + "=" * 80)
        print(f"▶ 计算：{os.path.relpath(fpath, args.model_dir)}")
        try:
            df = run_for_model(fpath, args.community_growth, args.min_growth, args.max_import, candidate_ex)
            df_nz = df[df["flux"] > 1e-12].sort_values("flux", ascending=False)
            if not df_nz.empty:
                print("非零 EX_*_m（全部）：")
                print(df_nz.to_string(index=False))
            else:
                print("无非零 EX_*_m 项")
            all_rows.append(df)
        except Exception as e:
            print(f"[异常] {type(e).__name__}: {e}")

    if not all_rows:
        print("\n没有可汇总的结果。")
        return

    out = pd.concat(all_rows, ignore_index=True)

    # 生成推荐培养基
    medium_df = recommend_medium(out, args.max_import)
    if medium_df.empty:
        print("\n[培养基建议] 没有正需求通量，无法给出建议。")
        return

    print("\n=== 基于单菌扫描得到的社区培养基建议 ===")
    print("(按需要该 EX 的模型数量降序，给出建议上界 suggested_upper_bound)\n")
    # 打印前 50 行概览
    preview_cols = ["reaction", "suggested_upper_bound"]
    print(medium_df.head(50)[preview_cols].to_string(index=False))

    # 导出推荐培养基
    out_path = os.path.abspath(args.out)
    try:
        medium_df.rename(columns={"suggested_upper_bound": "flux"}).to_csv(out_path, index=False)
        print("\n✅ 已导出推荐培养基：", out_path)
    except Exception as e:
        print("[错误] 推荐培养基导出失败：", e)


if __name__ == "__main__":
    main()
