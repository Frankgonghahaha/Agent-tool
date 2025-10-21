#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import json
import pandas as pd
import numpy as np
from itertools import combinations
# 并行相关
from joblib import Parallel, delayed
import multiprocessing

# ======== 1) 基本配置：你的数据目录（可按需修改） ========
BASE_DIR = "/Users/frank_gong/文档/生物智能体/20251015/信息表"
PATH_SHEET1 = os.path.join(BASE_DIR, "Sheet1_Complementarity.xlsx")
PATH_SHEET2 = os.path.join(BASE_DIR, "Sheet2_Species_environment.xlsx")
PATH_SHEET3 = os.path.join(BASE_DIR, "Sheet3_Function_enzyme_kact.xlsx")
PATH_SHEET4 = os.path.join(BASE_DIR, "Sheet4_species_enzyme.xlsx")
OUTPUT_PATH = os.path.join(BASE_DIR, "Result1_candidate_function_species.csv")

# PhyloMint 文件路径与列名（如有不同，可在此调整）
PATH_PHYLOMINT = os.path.join(BASE_DIR, "Sheet5_PhyloMint.csv")
PHYLO_COL_A = "A"
PHYLO_COL_B = "B"
PHYLO_COL_COMPETITION = "Competition"
PHYLO_COL_COMPLEMENTARITY = "Complementarity"
PHYLO_SUFFIX = "_CDS"  # functional/complement species → PhyloMint 的 ID 后缀

# ======== 2) 读取&清洗工具函数 ========
def read_first_sheet(path):
    """读 Excel 第一张表；若不存在报错更清晰。"""
    if not os.path.exists(path):
        raise FileNotFoundError(f"文件不存在: {path}")
    x = pd.ExcelFile(path)
    return pd.read_excel(path, sheet_name=x.sheet_names[0])

def normalize_name(x):
    """统一名称：去空白、全小写、半角化（只做简单处理）。"""
    if pd.isna(x):
        return None
    s = str(x).strip()
    s = s.replace("（", "(").replace("）", ")").replace("，", ",").replace("、", ",")
    s = " ".join(s.split())
    return s

def split_enzyme_list(s):
    """把 Sheet4 里的酶列表拆分；支持‘,’或‘、’等分隔。"""
    if pd.isna(s) or str(s).strip() == "":
        return []
    s = str(s).replace("、", ",").replace("；", ",").replace(";", ",")
    parts = [p.strip() for p in s.split(",") if p.strip()]
    return parts


def to_float_safe(x):
    try:
        return float(x)
    except Exception:
        return np.nan

# ======== 2.1) 软环境评分工具函数（三角隶属 + 指数尾部） ========
import math

def _tri_or_tail(x, xmin, xopt, xmax, k=math.log(10)):
    """三角隶属 + 区间外指数尾部；允许 xopt 缺失。返回 [0,1] 或 None（无法评分）。"""
    try:
        import numpy as _np
    except Exception:  # 兜底，避免未导入 numpy 的问题
        _np = None
    def _isnan(v):
        try:
            return _np.isnan(v) if _np is not None else pd.isna(v)
        except Exception:
            return pd.isna(v)
    if xmin is None or xmax is None or _isnan(xmin) or _isnan(xmax):
        return None
    if x is None or _isnan(x):
        return None
    rng = xmax - xmin
    if rng <= 0:
        return None
    # 区间内：三角/梯形
    if xmin <= x <= xmax:
        if (xopt is None) or _isnan(xopt) or not (xmin <= xopt <= xmax):
            center = (xmin + xmax) / 2.0
            half = max(rng / 2.0, 1e-9)
            return max(0.0, 1.0 - abs(x - center) / half)
        left = max(xopt - xmin, 1e-9)
        right = max(xmax - xopt, 1e-9)
        return max(0.0, 1.0 - abs(x - xopt) / max(left, right))
    # 区间外：指数尾部
    dist = xmin - x if x < xmin else x - xmax
    return math.exp(-k * (dist / max(rng, 1e-9)))


def _salt_soft(salt_star, salt_max, k=math.log(10)):
    """盐度软评分：未超过上限得 1；超过用指数衰减。缺值返回 None。"""
    if salt_max is None or pd.isna(salt_max) or salt_star is None or pd.isna(salt_star):
        return None
    if salt_star <= salt_max:
        return 1.0
    scale = max(1.0, salt_max)
    return math.exp(-k * ((salt_star - salt_max) / scale))


def _o2_soft(o2_value, o2_target):
    """氧环境软评分：匹配=1.0；未知/模糊=0.5；不匹配=0.2。"""
    val = str(o2_value).strip().lower() if pd.notna(o2_value) else ""
    tgt = str(o2_target or "").strip().lower()
    if val == tgt and val != "":
        return 1.0
    if val in {"", "unknown", "nan", "none"}:
        return 0.5
    return 0.2


def soft_env_score_row(row, T_star, pH_star, salt_star, o2_target,
                        wT=0.35, wPH=0.35, wS=0.10, wO2=0.20):
    """对单行环境数据计算软环境分，自动跳过缺失维度并重归一化。"""
    Tmin = row.get("temperature_minimum")
    Tmax = row.get("temperature_maximum")
    Top  = row.get("temperature_optimum_C")
    pHmin= row.get("ph_minimum")
    pHmax= row.get("ph_maximum")
    pHopt= row.get("ph_optimum")
    Smax = row.get("salinity_maximum")
    o2   = row.get("oxygen_tolerance")

    sT  = _tri_or_tail(T_star, Tmin, Top, Tmax)
    sPH = _tri_or_tail(pH_star, pHmin, pHopt, pHmax)
    sS  = _salt_soft(salt_star, Smax)
    sO2 = _o2_soft(o2, o2_target)

    parts, weights = [], []
    for sc, w in [(sT, wT), (sPH, wPH), (sS, wS), (sO2, wO2)]:
        if sc is not None:
            parts.append(sc * w)
            weights.append(w)
    if not parts:
        return 0.0
    return float(sum(parts) / sum(weights))

# ======== 3) 打印提示并收集目标工况 ========
print("请输入你的目标工况（按提示输入数值/选项）")
T_target = to_float_safe(input("目标温度 (°C)：").strip())
pH_target = to_float_safe(input("目标 pH：").strip())
sal_target = to_float_safe(input("目标盐度（% NaCl）：").strip())
O2_text = normalize_name(input("氧环境（好氧 / 厌氧 / 缺氧）：").strip())

if any(np.isnan(v) for v in [T_target, pH_target, sal_target]):
    print("❌ 温度 / pH / 盐度 请输入数字。")
    sys.exit(1)

# 氧环境映射：好氧 -> tolerant；厌氧/缺氧 -> not tolerant
if O2_text in ["好氧", "有氧", "aerobic", "氧气", "氧化"]:
    O2_target = "tolerant"
elif O2_text in ["厌氧", "缺氧", "微氧", "anaerobic", "anoxic", "microaerobic", "低氧"]:
    O2_target = "not tolerant"
else:
    print("⚠️ 未识别的氧环境输入，默认使用‘好氧’ → tolerant")
    O2_target = "tolerant"

print(f"\n✅ 目标工况：T={T_target}°C, pH={pH_target}, 盐度={sal_target}% NaCl, 氧环境={O2_text} → 映射为 {O2_target}\n")

# ======== 4) 读取四张表 ========
print("读取数据表...")
# 固定工作表名：
# Sheet1_Complementarity.xlsx -> Sheet1
# Sheet2_Species_environment.xlsx -> prediction
# Sheet3_Function_enzyme_kact.xlsx -> Sheet1
# Sheet4_species_enzyme.xlsx -> Sheet1

df_comp = pd.read_excel(PATH_SHEET1, sheet_name="Sheet1")
df_env  = pd.read_excel(PATH_SHEET2, sheet_name="prediction")
# —— 兼容不同环境表的物种标识列名 —— 
# 期望存在 “strain” 列；若没有，自动从常见命名中选择并重命名为 “strain”
_strain_aliases = {
    "strain", "species", "Species"
}
if "strain" not in df_env.columns:
    candidates = [c for c in df_env.columns if str(c).strip() in _strain_aliases]
    if candidates:
        chosen = candidates[0]
        df_env = df_env.rename(columns={chosen: "strain"})
        print(f"⚠️ 环境表未找到列 'strain'，已使用 '{chosen}' 列作为物种标识并重命名为 'strain'。")
    else:
        raise KeyError(f"环境表缺少用于物种标识的列（期望 'strain' 或常见别名）。当前列: {list(df_env.columns)}")
df_kcat = pd.read_excel(PATH_SHEET3, sheet_name="Sheet1")
df_map  = pd.read_excel(PATH_SHEET4, sheet_name="Sheet1")

# ======== 5) 固定列映射（数据接入层：确定版） ========
# 互补关系（Sheet1_Complementarity.xlsx / Sheet1）
df_comp = df_comp.rename(columns={
    "Species": "functional_species",
    "Complementarity Species": "complement_species",
    "Competition": "competition_index",
    "Complementarity": "complementarity_index",
    "Delta": "delta_index",
})[[
    "functional_species",
    "complement_species",
    "competition_index",
    "complementarity_index",
    "delta_index",
]]

# 规范名称，保证与 df_map["species"] 一致
df_comp["functional_species"] = df_comp["functional_species"].apply(normalize_name)
df_comp["complement_species"] = df_comp["complement_species"].apply(normalize_name)

# 环境（Sheet2_Species_environment.xlsx / prediction）
# 直接使用文件中的标准列名；若缺少 salinity_optimum 列则不强制
env_cols_expected = [
    "strain",
    "temperature_optimum_C",
    "temperature_minimum",
    "temperature_maximum",
    "ph_optimum",
    "ph_minimum",
    "ph_maximum",
    "salinity_optimum",
    "salinity_minimum",
    "salinity_maximum",
    "oxygen_tolerance",
]
existing_env_cols = [c for c in env_cols_expected if c in df_env.columns]
df_env_use = df_env[existing_env_cols].copy()

# 类型标准化
df_env_use["oxygen_tolerance"] = (
    df_env_use["oxygen_tolerance"].astype(str).str.strip().str.lower()
)
for c in [
    "temperature_optimum_C",
    "temperature_minimum",
    "temperature_maximum",
    "ph_optimum",
    "ph_minimum",
    "ph_maximum",
    "salinity_optimum",
    "salinity_minimum",
    "salinity_maximum",
]:
    if c in df_env_use.columns:
        df_env_use[c] = pd.to_numeric(df_env_use[c], errors="coerce")

# 物种-酶映射（Sheet4_species_enzyme.xlsx / Sheet1）
df_map = df_map.rename(columns={
    "Functional Species": "species",
    "Function Enzyme": "enzymes",
})[["species", "enzymes"]]

def split_enzyme_list_fixed(s):
    if pd.isna(s) or str(s).strip() == "":
        return []
    s = str(s).replace("、", ",").replace("；", ",").replace(";", ",")
    return [p.strip() for p in s.split(",") if p.strip()]

df_map["species"] = df_map["species"].apply(normalize_name)
df_map["enzymes_list"] = df_map["enzymes"].apply(split_enzyme_list_fixed)
df_map = df_map[["species", "enzymes_list"]].dropna(subset=["species"]).reset_index(drop=True)

# 酶-Kcat（Sheet3_Function_enzyme_kact.xlsx / Sheet1）
df_kcat = df_kcat.rename(columns={
    "Enzyme": "enzyme",
    "Kcat value (1/s)": "kcat",
})[["enzyme", "kcat"]]

df_kcat["_enzyme_norm"] = df_kcat["enzyme"].astype(str).str.strip().str.lower()
df_kcat["_kcat"] = pd.to_numeric(df_kcat["kcat"], errors="coerce")
kcat_map = df_kcat.groupby("_enzyme_norm")["_kcat"].median().to_dict()

# ======== 6) 计算功能菌的 Kcat_max / Kcat_mean ========
def enzyme_to_kcat(enzyme_name):
    if enzyme_name is None:
        return np.nan
    key = str(enzyme_name).strip().lower()
    return kcat_map.get(key, np.nan)

def species_kcat_stats(enz_list):
    vals = [enzyme_to_kcat(e) for e in (enz_list or [])]
    vals = [v for v in vals if not pd.isna(v)]
    if not vals:
        return np.nan, np.nan
    return float(np.max(vals)), float(np.mean(vals))

species_rec = []
for _, row in df_map.iterrows():
    s = row["species"]
    enz = row["enzymes_list"]
    kmax, kmean = species_kcat_stats(enz)
    species_rec.append({
        "species": s,
        "enzymes": ",".join(enz),
        "kcat_max": kmax,
        "kcat_mean": kmean,
    })
df_func = pd.DataFrame(species_rec)

# ======== 7) 环境匹配：按目标工况筛选 ========
df_env_use["match_temp"] = True
if {"temperature_minimum", "temperature_maximum"}.issubset(df_env_use.columns):
    df_env_use["match_temp"] = (
        (df_env_use["temperature_minimum"] <= T_target)
        & (T_target <= df_env_use["temperature_maximum"])
    )

df_env_use["match_ph"] = True
if {"ph_minimum", "ph_maximum"}.issubset(df_env_use.columns):
    df_env_use["match_ph"] = (
        (df_env_use["ph_minimum"] <= pH_target)
        & (pH_target <= df_env_use["ph_maximum"])
    )

df_env_use["match_salt"] = True
if "salinity_maximum" in df_env_use.columns:
    df_env_use["match_salt"] = (
        (sal_target <= df_env_use["salinity_maximum"]) | df_env_use["salinity_maximum"].isna()
    )

# 氧匹配（好氧 → tolerant；厌氧/缺氧 → not tolerant）
df_env_use["match_o2"] = (df_env_use["oxygen_tolerance"].str.lower() == O2_target)

# 总匹配
df_env_use["env_match_all"] = df_env_use[["match_temp", "match_ph", "match_salt", "match_o2"]].all(axis=1)

# 软环境匹配分（0~1）：基于目标工况的连续评分
df_env_use["env_soft_score"] = df_env_use.apply(
    lambda r: soft_env_score_row(r, T_target, pH_target, sal_target, O2_target),
    axis=1
)

# 为互补菌检测准备：按 strain 建立环境匹配的查找表
_env_cols_needed = [
    "strain","env_match_all","env_soft_score","fail_reasons","match_o2","match_temp","match_ph","match_salt",
    "temperature_minimum","temperature_maximum","ph_minimum","ph_maximum","salinity_maximum","oxygen_tolerance"
]
cols_present = [c for c in _env_cols_needed if c in df_env_use.columns]
df_env_lookup = df_env_use[cols_present].copy()

# 构建字典：strain -> 环境是否通过 及 明细
def _reason_row(r):
    reasons = []
    if not r.get("match_o2", True):
        reasons.append("oxygen_mismatch")
    if not r.get("match_temp", True):
        reasons.append("temperature_out_of_range")
    if not r.get("match_ph", True):
        reasons.append("pH_out_of_range")
    if not r.get("match_salt", True):
        reasons.append("salinity_exceeds_max")
    return ";".join(reasons) if reasons else "ok"

df_env_lookup["fail_reasons"] = df_env_lookup.apply(_reason_row, axis=1)
_env_by_strain = {
    str(r["strain"]): {k: r.get(k) for k in ["env_match_all","fail_reasons","env_soft_score"]}
    for _, r in df_env_lookup.iterrows()
}

# ======== 8) 合并功能信息与环境表 ========
merged1 = df_func.merge(df_env_use, how="left", left_on="species", right_on="strain")

# ======== 9) 输出候选：先看环境完全匹配的功能菌 ========
candidates = merged1[ merged1["env_match_all"] == True ].copy()

# 排序：优先看 kcat_max（降序），再看 kcat_mean
candidates = candidates.sort_values(by=["kcat_max","kcat_mean"], ascending=[False, False])

# 选择输出列
out_cols = [
    "species", "enzymes", "kcat_max", "kcat_mean",
    "temperature_optimum_C","temperature_minimum","temperature_maximum",
    "ph_optimum","ph_minimum","ph_maximum",
    "salinity_optimum","salinity_minimum","salinity_maximum",
    "oxygen_tolerance"
]
out_cols = [c for c in out_cols if c in candidates.columns]
# === 对每个候选功能菌，统计其互补微生物在目标工况下的通过情况 ===
comp_stats = {"complement_total": [], "complement_pass": [], "complement_fail": [],
              "complement_pass_names": [], "complement_fail_names": []}

for _, crow in candidates.iterrows():
    func_name = crow.get("species")
    # 找到该功能菌的所有互补菌（去重）
    comps = df_comp.loc[df_comp["functional_species"] == func_name, "complement_species"].dropna().astype(str).unique().tolist()
    total = len(comps)
    pass_list, fail_list = [], []
    for c in comps:
        info = _env_by_strain.get(c)
        if info is None:
            fail_list.append(f"{c} (no_env_record)")
        else:
            if bool(info.get("env_match_all", False)):
                pass_list.append(c)
            else:
                reason = info.get("fail_reasons", "fail")
                fail_list.append(f"{c} ({reason})")
    comp_stats["complement_total"].append(total)
    comp_stats["complement_pass"].append(len(pass_list))
    comp_stats["complement_fail"].append(len(fail_list))
    comp_stats["complement_pass_names"].append(";".join(pass_list))
    comp_stats["complement_fail_names"].append(";".join(fail_list))

# 追加到输出表
candidates_out = candidates[out_cols].reset_index(drop=True)
for k, v in comp_stats.items():
    candidates_out[k] = v

# ======== 10) 打印预览 & 保存 ========
print("\n===== 满足目标工况的‘功能菌’候选（按 kcat_max 排序，前 10 条）=====")
print(candidates_out.head(10).to_string(index=False))

candidates_out.to_csv(OUTPUT_PATH, index=False, encoding="utf-8")
print(f"\n✅ 已保存候选菌清单到：{OUTPUT_PATH}")

print("\n[核对] comp列:", list(df_comp.columns))
print("[核对] env 列:", list(df_env_use.columns))
print("[核对] kcat列:", list(df_kcat.columns))
print("[核对] map 列:", list(df_map.columns))

print("\n提示：")
print("1) 这一步仅完成‘功能菌’在目标工况下的筛选与功能强度汇总（Kcat）。")
print("2) 下一步可把互补微生物（Sheet1）引入，做组合打分与群落优化。")
print("3) 若环境表的 ‘strain’ 名称与物种表的 ‘species’ 不一致，需要提供一个映射表做更强的对齐。")
# ======== 11) 功能菌 + 互补菌 打分模块 (S_microbe) ========

print("\n===== 对功能菌和互补菌进行综合打分 (S_microbe) =====")

# 1️⃣ 先准备功能菌表（已经在 candidates_out 里）
func_records = []
for _, r in candidates_out.iterrows():
    # 使用软环境分：候选功能菌通常已在区间内，但软分可区分远近
    _envrow = df_env_use.loc[df_env_use["strain"] == r["species"]]
    _envsoft = float(_envrow.iloc[0]["env_soft_score"]) if (not _envrow.empty and pd.notna(_envrow.iloc[0]["env_soft_score"])) else 1.0
    func_records.append({
        "species": r["species"],
        "kcat_max": r.get("kcat_max", np.nan),
        "kcat_mean": r.get("kcat_mean", np.nan),
        "enzyme_diversity": len(str(r.get("enzymes", "")).split(",")) if pd.notna(r.get("enzymes")) else 0,
        "environment_match": _envsoft,
        "source": "functional"
    })

# 2️⃣ 收集所有通过环境匹配的互补菌名单
comp_pass_all = []
for names in candidates_out["complement_pass_names"]:
    if isinstance(names, str) and names.strip():
        comp_pass_all.extend([n.strip() for n in names.split(";") if n.strip()])
comp_pass_all = sorted(set(comp_pass_all))

print(f"共检测到 {len(comp_pass_all)} 个通过目标工况的互补菌。")

# 3️⃣ 互补菌打分数据收集
comp_records = []
for cname in comp_pass_all:
    # 环境匹配信息
    env_row = df_env_use.loc[df_env_use["strain"] == cname]
    if not env_row.empty and pd.notna(env_row.iloc[0].get("env_soft_score", np.nan)):
        env_score = float(env_row.iloc[0]["env_soft_score"])  # 0~1 连续分
    else:
        env_score = 0.0  # 缺失则视为不匹配

    # 酶与 Kcat 信息
    func_row = df_map.loc[df_map["species"] == cname]
    if func_row.empty:
        kmax, kmean, enz_list = np.nan, np.nan, []
    else:
        enz_list = func_row.iloc[0]["enzymes_list"]
        kmax, kmean = species_kcat_stats(enz_list)
    enz_div = len(enz_list)

    comp_records.append({
        "species": cname,
        "kcat_max": kmax,
        "kcat_mean": kmean,
        "enzyme_diversity": enz_div,
        "environment_match": env_score,
        "source": "complement"
    })

# 4️⃣ 合并功能菌与互补菌
df_all = pd.DataFrame(func_records + comp_records)

# 5️⃣ 归一化并计算 S_microbe
def normalize_01(series):
    if series.min() == series.max():
        return np.ones_like(series)
    return (series - series.min()) / (series.max() - series.min())

df_all["f_Kcat"] = normalize_01(df_all["kcat_max"].fillna(0))
df_all["f_EnvMatch"] = df_all["environment_match"].fillna(0)
df_all["f_EnzymeDiversity"] = normalize_01(df_all["enzyme_diversity"].fillna(0))

# 权重
wK, wE, wD = 0.5, 0.4, 0.1
df_all["S_microbe"] = (
    wK * df_all["f_Kcat"]
    + wE * df_all["f_EnvMatch"]
    + wD * df_all["f_EnzymeDiversity"]
)

# 6️⃣ 排序与输出
df_all = df_all.sort_values(by="S_microbe", ascending=False)
output_path = os.path.join(BASE_DIR, "Result2_candidate_scores.csv")
df_all.to_csv(output_path, index=False, encoding="utf-8")

print(f"\n✅ 已保存功能菌 + 互补菌打分结果到: {output_path}")
print(df_all.head(10).to_string(index=False))

# ======== 12) 生成 species 间互补/竞争指数表（融合 Sheet1 与 PhyloMint） ========
print("\n===== 基于 Sheet1 与 PhyloMint 生成互作矩阵 =====")

# 12.1 构建 species 全集（功能菌 + 通过环境匹配的互补菌）
species_all = sorted(pd.Series(df_all["species"].astype(str).tolist()).dropna().unique().tolist())
print(f"将生成 {len(species_all)} 个物种的两两互作。")

# 12.2 建立“功能菌-其对应互补菌”的集合（来自 Sheet1_Complementarity.xlsx）
corresponding_pairs_set = set(zip(
    df_comp["functional_species"].astype(str),
    df_comp["complement_species"].astype(str)
))

# 12.3 读取 PhyloMint
if not os.path.exists(PATH_PHYLOMINT):
    print(f"⚠️ 未找到 PhyloMint 文件：{PATH_PHYLOMINT}，仅输出 Sheet1 内已有配对。")
    df_phy = pd.DataFrame(columns=[PHYLO_COL_A, PHYLO_COL_B, PHYLO_COL_COMPETITION, PHYLO_COL_COMPLEMENTARITY])
else:
    df_phy = pd.read_csv(PATH_PHYLOMINT)
    # 校验并标准化
    for c in [PHYLO_COL_A, PHYLO_COL_B]:
        if c not in df_phy.columns:
            raise ValueError(f"PhyloMint 缺少列：{c}")
        df_phy[c] = df_phy[c].astype(str).str.strip()
    for c in [PHYLO_COL_COMPETITION, PHYLO_COL_COMPLEMENTARITY]:
        if c not in df_phy.columns:
            raise ValueError(f"PhyloMint 缺少列：{c}")
        df_phy[c] = pd.to_numeric(df_phy[c], errors="coerce")

# —— 将 PhyloMint 预聚合为无序键字典，便于 O(1) 查找 ——
if 'df_phy' in locals() and not df_phy.empty:
    df_phy["_A"] = df_phy[PHYLO_COL_A].astype(str).str.strip()
    df_phy["_B"] = df_phy[PHYLO_COL_B].astype(str).str.strip()
    _key_arr = np.where(df_phy["_A"] <= df_phy["_B"],
                        df_phy["_A"] + "||" + df_phy["_B"],
                        df_phy["_B"] + "||" + df_phy["_A"])
    df_phy["_key"] = _key_arr
    _phy_mean = (df_phy.groupby("_key", as_index=True)[[PHYLO_COL_COMPETITION, PHYLO_COL_COMPLEMENTARITY]].mean())
    phy_dict = _phy_mean.to_dict(orient="index")
else:
    phy_dict = {}

def _phylomint_pair_avg(a, b):
    """A=a_CDS,B=b_CDS 与 A=b_CDS,B=a_CDS 的 Competition/Complementarity 合并求均值"""
    a_id = f"{a}{PHYLO_SUFFIX}"
    b_id = f"{b}{PHYLO_SUFFIX}"
    key = a_id + "||" + b_id if a_id <= b_id else b_id + "||" + a_id
    rec = phy_dict.get(key)
    if rec is None:
        return np.nan, np.nan, np.nan
    comp_mean = float(rec.get(PHYLO_COL_COMPETITION, np.nan))
    compl_mean = float(rec.get(PHYLO_COL_COMPLEMENTARITY, np.nan))
    delta = float(compl_mean - comp_mean) if (not pd.isna(comp_mean) and not pd.isna(compl_mean)) else np.nan
    return comp_mean, compl_mean, delta

# —— Sheet1 的功能↔互补关系预构建成字典（有序键）——
sheet1_dict = {
    (str(r["functional_species"]), str(r["complement_species"])): (
        float(r["competition_index"]),
        float(r["complementarity_index"]),
        float(r["delta_index"])
    )
    for _, r in df_comp.iterrows()
}

def _choose_sheet1_if_corresponding(a, b):
    """
    若 (a,b) 是“功能菌-其对应互补菌”的关系（出现在 Sheet1 中），
    则直接返回 Sheet1 数值（注意 df_comp 已按均值聚合/或按单行使用）。
    """
    rec = sheet1_dict.get((a, b))
    if rec is not None:
        c, m, d = rec
        return c, m, d, True
    rec2 = sheet1_dict.get((b, a))
    if rec2 is not None:
        c, m, d = rec2
        return c, m, d, True
    return np.nan, np.nan, np.nan, False

# 并行计算互作矩阵
cpu_total = multiprocessing.cpu_count()
n_jobs = max(1, cpu_total - 1)
print(f"💡 并行计算互作矩阵，使用 {n_jobs}/{cpu_total} 个CPU核心")

def _compute_pair_record(a, b):
    comp_mean, compl_mean, delta, used_sheet1 = _choose_sheet1_if_corresponding(a, b)
    source = "from_sheet1" if used_sheet1 else "from_phylomint"
    if not used_sheet1:
        comp_mean, compl_mean, delta = _phylomint_pair_avg(a, b)
    return {
        "functional_species": a,
        "complement_species": b,
        "competition_index": comp_mean,
        "complementarity_index": compl_mean,
        "delta_index": delta,
        "source": source
    }

interaction_records = Parallel(n_jobs=n_jobs, backend="loky", prefer="processes")(
    delayed(_compute_pair_record)(a, b) for a, b in combinations(species_all, 2)
)

df_interact = pd.DataFrame(interaction_records)

# 12.5 对称补齐（便于任意方向检索）
df_rev = df_interact.rename(columns={
    "functional_species": "complement_species",
    "complement_species": "functional_species"
})
df_all_pairs = pd.concat([df_interact, df_rev], ignore_index=True)
df_all_pairs = df_all_pairs.sort_values(by=["functional_species", "complement_species"]).reset_index(drop=True)

# 12.6 保存输出
merged_out_path = os.path.join(BASE_DIR, "Result3_pair_Com_index.csv")
df_all_pairs.to_csv(merged_out_path, index=False, encoding="utf-8")
print(f"✅ 已生成融合互作矩阵：{merged_out_path}")

# 12.7 预览
print(df_all_pairs.head(12).to_string(index=False)) 