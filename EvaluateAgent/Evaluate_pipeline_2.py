# -*- coding: utf-8 -*-
from __future__ import annotations
"""
Eva_pipeline 2 — 社区在推荐培养基上的两步优化（DBP 摄入，应用 MICOM community.medium）

流程：
1) 从 --model-dir 递归读取全部 SBML（.xml/.sbml），构建 MICOM Community。
2) 读取 --medium-csv（reaction, flux），并强制加入 EX_dbp_m=20。
3) 使用 MICOM community.medium 应用培养基并建立耦合：通过设置 community.medium，保证社区层与成员层 EX 通量守恒。
4) 第一步：最大化群落生长 -> 得到 Biomass_max。
5) 第二步：加入约束 growth ≥ alpha * Biomass_max，目标改为“最小化 EX_dbp_m”，求得 DBP 摄入最大值（负值越小越强）。

注意：
- 若社区中不存在 EX_dbp_m，则第二步会提示并返回摄入=nan。
- 若第一步不可行（Biomass_max=0），第二步将直接跳过。
"""

print("RUNNING FILE:", __file__)
import os
import argparse
import tempfile
import shutil
import warnings
from typing import List, Tuple, Optional, Any
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
from micom import Community
def apply_medium_via_micom(comm: Community, medium_df: pd.DataFrame) -> Tuple[int, int]:
    """使用 MICOM 的 community.medium 设置培养基上限（正值=最大摄入）。
    返回 (applied_count, missing_count)。"""
    # 只保留正上限；将 NaN/负值置 0
    medium_df = medium_df.copy()
    medium_df["upper"] = pd.to_numeric(medium_df["upper"], errors="coerce").fillna(0.0)
    medium_df.loc[medium_df["upper"] < 0, "upper"] = 0.0
    pairs = [(str(rid).strip(), float(up)) for rid, up in zip(medium_df["reaction"], medium_df["upper"])]
    # 统计存在/缺失
    applied, missing = 0, 0
    med_dict = {}
    for rid, up in pairs:
        if rid in comm.reactions:
            med_dict[rid] = up
            applied += 1
        else:
            missing += 1
    # 赋值给 MICOM，建立 medium 守恒耦合
    try:
        comm.medium = med_dict
    except Exception as e:
        print("[错误] 设置 community.medium 失败：", e)
        raise
    return applied, missing

warnings.filterwarnings("ignore", category=FutureWarning)

EXTERNAL_COMPARTMENT_SYNONYMS = {"C_e", "ext", "external", "extracellular"}
TARGET_EXTERNAL = "e"
MAX_IMPORT = 1000.0  # 社区层 EX 上界

DBP_EX_ID = "EX_dbp_m"

# 二分与阈值默认参数（不暴露 CLI）
EPSILON = 0.1      # 最小要求的摄入强度 ε（f ≤ -ε）
BIS_TOL = 1e-3     # 二分终止阈值
MAX_ITER = 30      # 二分最大迭代
# ---------- 固定反应上下界的辅助上下文管理器 ----------
from contextlib import contextmanager

@contextmanager
def fixed_bound(comm: Community, rxn_id: str, value: float):
    """将反应 rxn_id 的下上界同时固定为 value，with 退出时恢复。"""
    rxn = comm.reactions.get_by_id(rxn_id)
    old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
    try:
        rxn.lower_bound = value
        rxn.upper_bound = value
        yield rxn
    finally:
        rxn.lower_bound = old_lb
        rxn.upper_bound = old_ub

# ---------- CT 最大增长（固定 EX 通量） ----------
def ct_max_growth_under_ex(comm: Community, ex_id: str, f_value: float) -> float:
    """固定某交换通量为 f_value，返回 CT(f=1.0) 下的最大社区增长率。"""
    if ex_id not in comm.reactions:
        return float("nan")
    with fixed_bound(comm, ex_id, f_value):
        try:
            sol = comm.cooperative_tradeoff(fraction=1.0, fluxes=False)
            return float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
        except Exception:
            return 0.0

# ---------- 不加增长约束时的最小 EX ----------
def unconstrained_min_ex(comm: Community, ex_id: str) -> float:
    """在不加增长约束的情况下最小化 ex_id（返回最小值，负值越小代表摄入越强）。"""
    if ex_id not in comm.reactions:
        return float("nan")
    rxn = comm.reactions.get_by_id(ex_id)
    old_obj, old_dir = comm.objective, comm.objective_direction
    try:
        comm.objective = rxn
        comm.objective_direction = "min"
        sol = comm.optimize()
        val = float(getattr(sol, "objective_value", 0.0) or 0.0)
        return val
    finally:
        comm.objective = old_obj
        comm.objective_direction = old_dir


# ---------- 辅助函数 ----------
def find_biomass_rxn_id_model(model) -> Optional[str]:
    """稳健识别单菌模型的 biomass 反应。
    1) 从 objective 表达式映射 optlang 变量 → Reaction（forward/reverse_variable）
    2) 兜底用 id/name 含 'biomass' 或 'growth' 的反应
    返回反应 id 或 None。
    """
    try:
        expr = model.objective.expression
        coeffs = expr.as_coefficients_dict()
        for var, coef in coeffs.items():
            try:
                if float(coef) == 0.0:
                    continue
            except Exception:
                continue
            for r in model.reactions:
                try:
                    if getattr(r, "forward_variable", None) is var or getattr(r, "reverse_variable", None) is var:
                        return r.id
                except Exception:
                    continue
        if len(coeffs) == 1:
            var = next(iter(coeffs.keys()))
            vname = (getattr(var, "name", "") or "").lower()
            for r in model.reactions:
                if r.id.lower() in vname:
                    return r.id
    except Exception:
        pass
    for r in model.reactions:
        rid = r.id.lower()
        rname = (r.name or "").lower()
        if ("biomass" in rid) or ("growth" in rid) or ("biomass" in rname) or ("growth" in rname):
            return r.id
    return None


# ---------- 工具函数 ----------
def discover_models(models_dir: str) -> List[str]:
    """递归发现 .xml/.sbml 模型"""
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


def build_community_from_dir(models_dir: str) -> Tuple[Community, str]:
    """把目录下所有 SBML 作为成员构建一个社区。返回 (community, tmp_dir)。"""
    model_paths = discover_models(models_dir)
    if not model_paths:
        raise FileNotFoundError(f"未在目录下发现 SBML：{models_dir}")

    tmpd = tempfile.mkdtemp(prefix="community_from_dir_")
    try:
        rows = []
        for f in model_paths:
            name = os.path.splitext(os.path.basename(f))[0]
            model = read_sbml_model(f)
            model = normalize_external_compartment(model)
            # 检测 biomass 反应 ID 并尽量设为该模型目标
            biomass_id = find_biomass_rxn_id_model(model)
            print(f"  成员 {name} 的 biomass 反应: {biomass_id if biomass_id else '未检测到'}")
            if biomass_id:
                try:
                    model.objective = biomass_id
                except Exception:
                    pass
            tmp_sbml = os.path.join(tmpd, f"{name}.xml")
            write_sbml_model(model, tmp_sbml)
            rows.append({"id": name, "file": tmp_sbml, "abundance": 1.0, "biomass": biomass_id})
        tax = pd.DataFrame(rows).set_index("id", drop=False)
        com = Community(tax, name="COMM_DBP")
        return com, tmpd
    except Exception:
        shutil.rmtree(tmpd, ignore_errors=True)
        raise


def read_medium_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # 兼容列名 flux 或 suggested_upper_bound
    if "reaction" not in df.columns:
        raise ValueError("培养基 CSV 必须包含 'reaction' 列")
    if "flux" in df.columns:
        df = df[["reaction", "flux"]].rename(columns={"flux": "upper"})
    elif "suggested_upper_bound" in df.columns:
        df = df[["reaction", "suggested_upper_bound"]].rename(columns={"suggested_upper_bound": "upper"})
    else:
        raise ValueError("培养基 CSV 必须包含 'flux' 或 'suggested_upper_bound' 列")
    # 合法化
    df["reaction"] = df["reaction"].astype(str)
    df["upper"] = pd.to_numeric(df["upper"], errors="coerce").fillna(0.0)
    return df


def medium_plus_dbp(med: pd.DataFrame, dbp_upper: float = 20.0) -> pd.DataFrame:
    """在培养基清单上追加/提升 EX_dbp_m 的上限"""
    med = med.copy()
    idx = med.index[med["reaction"] == DBP_EX_ID]
    if len(idx) > 0:
        med.loc[idx, "upper"] = med.loc[idx, "upper"].clip(lower=dbp_upper)
    else:
        med = pd.concat([med, pd.DataFrame({"reaction": [DBP_EX_ID], "upper": [dbp_upper]})], ignore_index=True)
    return med


def soft_apply_medium(comm: Community, medium_df: pd.DataFrame, ub_default: float = MAX_IMPORT) -> Tuple[int, int]:
    """
    软应用培养基：仅对 medium_df 的 EX 设置 lb/ub（lb=-upper, ub=ub_default），不关闭其它 EX。
    返回 (applied_count, skipped_count)。
    """
    applied, skipped = 0, 0
    for _, row in medium_df.iterrows():
        rid = str(row["reaction"]).strip()
        up = float(row["upper"])
        if rid in comm.reactions:
            rxn = comm.reactions.get_by_id(rid)
            rxn.lower_bound = -abs(up)
            rxn.upper_bound = ub_default
            applied += 1
        else:
            skipped += 1
    return applied, skipped


def find_community_biomass_rxn(comm: Community) -> Optional[str]:
    """更稳健地识别社区生长反应：
    ① 从 objective 表达式映射变量 → Reaction；② 常见固定 ID；③ 关键词兜底。
    """
    # ① 从 objective 表达式解析
    try:
        expr = comm.objective.expression
        coeffs = expr.as_coefficients_dict()
        for var, coef in coeffs.items():
            try:
                if float(coef) == 0.0:
                    continue
            except Exception:
                continue
            for r in comm.reactions:
                try:
                    if getattr(r, "forward_variable", None) is var or getattr(r, "reverse_variable", None) is var:
                        return r.id
                except Exception:
                    continue
        if len(coeffs) == 1:
            var = next(iter(coeffs.keys()))
            vname = (getattr(var, "name", "") or "").lower()
            for r in comm.reactions:
                if r.id.lower() in vname:
                    return r.id
    except Exception:
        pass
    # ② 常见固定 ID
    for cid in ["community_biomass", "community_growth", "Community_biomass"]:
        if cid in comm.reactions:
            return cid
    # ③ 关键词兜底
    for r in comm.reactions:
        low = r.id.lower()
        if "community" in low and ("biomass" in low or "growth" in low):
            return r.id
    return None


def step1_max_growth(comm: Community) -> float:
    """第一步：最大化群落生长，返回 Biomass_max。使用 CT 优化（如不可行则回退普通优化）。"""
    try:
        sol = comm.cooperative_tradeoff(fraction=1.0, fluxes=False)
        mu = float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
        return mu
    except Exception:
        sol = comm.optimize()
        mu = float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
        return mu



def step2_max_dbp_uptake(comm: Community, alpha: float, biomass_max: float) -> Tuple[float, float]:
    """
    二阶段（循环参数模拟风格）：
      - 目标：在 μ_comm ≥ α·μ* 的前提下，强制 f ≤ -ε，并最小化 EX_dbp_m（越负越好）
      - 实现：对 f 做二分查找。每一步将 EX_dbp_m 固定为 f，调用 CT(f=1.0) 求最大增长，判定是否 ≥ α·μ*。
    返回 (growth_at_f*, f_star)。
    """
    if DBP_EX_ID not in comm.reactions:
        print(f"[警告] 社区中不存在 {DBP_EX_ID}，无法度量 DBP 摄入。")
        return float("nan"), float("nan")

    # 先算不加增长约束时的最小摄入 f_min
    f_min = unconstrained_min_ex(comm, DBP_EX_ID)
    print(f"[阶段二初始化] 不加增长约束时的最小 {DBP_EX_ID} = {f_min:.6f}")

    target = alpha * biomass_max
    print(f"[阶段二目标] 需要满足增长 μ_comm ≥ α·μ* = {alpha:.3f} × {biomass_max:.6f} = {target:.6f}")

    # 右端点：至少要达到 -ε 的摄入
    f_hi = -float(EPSILON)
    mu_hi = ct_max_growth_under_ex(comm, DBP_EX_ID, f_hi)
    print(f"[可行性判定] f_hi={f_hi:.6f} → μ*(f_hi)={mu_hi:.6f}")

    if mu_hi < target - 1e-12:
        # 连 -ε 都不可行，则在 [0, -ε] 里找最小可行的负值
        left, right = 0.0, f_hi
        mu_left = ct_max_growth_under_ex(comm, DBP_EX_ID, left)
        print(f"[右端点不可行] 尝试在 [0, {f_hi:.6f}] 内二分；μ*(0)={mu_left:.6f}")
        if mu_left < target - 1e-12:
            print("[失败] 连 f=0 都不满足增长约束，无法给出负摄入。")
            return mu_left, 0.0
        best = left
        for _ in range(MAX_ITER):
            mid = 0.5 * (left + right)
            mu_mid = ct_max_growth_under_ex(comm, DBP_EX_ID, mid)
            if mu_mid >= target:
                best = mid
                right = mid
            else:
                left = mid
            if abs(right - left) <= BIS_TOL:
                break
        return ct_max_growth_under_ex(comm, DBP_EX_ID, best), best

    # 右端点可行：在 [f_min, f_hi] 内找“最负且可行”的 f
    # 注意 f_min 可能比 f_hi 更负（更小）；确保区间有序：f_left ≤ f_right
    f_left = min(f_min, f_hi)
    f_right = max(f_min, f_hi)
    # 保证 f_left ≤ f_right 且 f_right=f_hi（可行），若 f_left 不可行则向右收缩
    mu_left = ct_max_growth_under_ex(comm, DBP_EX_ID, f_left)
    print(f"[初始化区间] f_left={f_left:.6f} → μ*(f_left)={mu_left:.6f}")
    if mu_left < target - 1e-12:
        # 向右收缩直到可行
        left, right = f_left, f_hi
        best = right
        for _ in range(MAX_ITER):
            mid = 0.5 * (left + right)
            mu_mid = ct_max_growth_under_ex(comm, DBP_EX_ID, mid)
            if mu_mid >= target:
                best = mid
                left = mid  # 注意这里向更负的一侧推进需小心，但由于 mid 介于左负更大与右较少负之间，这样递推可收敛
            else:
                right = mid
            if abs(right - left) <= BIS_TOL:
                break
        return ct_max_growth_under_ex(comm, DBP_EX_ID, best), best

    # 若 f_left 本身也可行，说明整个区间都可行，则取最负的 f_left
    return ct_max_growth_under_ex(comm, DBP_EX_ID, f_left), f_left


# ---------- 主程序 ----------
def main():
    parser = argparse.ArgumentParser(description="Pipeline 2：推荐培养基 + 两步优化（最大生长 → 在α下最大化 DBP 摄入）")
    parser.add_argument("--model-dir", required=True, help="代谢模型目录（递归搜索 .xml/.sbml）")
    parser.add_argument("--medium-csv", required=True, help="培养基 CSV（两列：reaction, flux 或 suggested_upper_bound）")
    parser.add_argument("--alpha", type=float, default=0.7, help="第二步生长下界系数 alpha（默认 0.7）")
    parser.add_argument("--out", help="导出结果 CSV 路径（默认写入 model-dir/Result6_community.csv）")
    args = parser.parse_args()

    print("参数：")
    print("  model_dir  =", args.model_dir)
    print("  medium_csv =", args.medium_csv)
    print("  alpha      =", args.alpha)

    # 读取与修正培养基
    medium = read_medium_csv(args.medium_csv)
    medium = medium_plus_dbp(medium, dbp_upper=20.0)
    print("\n=== 推荐培养基（含 EX_dbp_m=20） ===")
    print(medium.to_string(index=False))

    # 构建社区
    com, tmpd = build_community_from_dir(args.model_dir)
    try:
        provided_biomass = sum(1 for x in getattr(com.taxa, 'biomass', []) if x)
    except Exception:
        provided_biomass = 'N/A'
    print(f"\n[社区信息] 成员数={len(com.taxa)}, 已提供成员biomass={provided_biomass}")
    try:
        # 使用 MICOM 的 community.medium 应用培养基（建立 community↔members 守恒耦合）
        applied, missing = apply_medium_via_micom(com, medium)
        print(f"\n[培养基应用] 通过 MICOM 耦合设置 {applied} 个 EX；输入中社区不存在的 {missing} 个（已忽略）。")
        # 可选：确认 EX_dbp_m 是否在当前 medium 中
        if DBP_EX_ID not in com.medium:
            print(f"[警告] 当前 medium 中未包含 {DBP_EX_ID}；将无法限制社区对 DBP 的总摄入。")

        # 第一步：最大生长
        biomass_max = step1_max_growth(com)
        print(f"\n阶段一：最大群落生长 Biomass_max = {biomass_max:.6f}")

        # 第二步：在 alpha*Biomass_max 下最小化 EX_dbp_m
        if biomass_max <= 1e-12:
            print("\n[警告] Biomass_max 为 0，跳过阶段二优化。")
            stage2_growth, dbp_flux = float("nan"), float("nan")
        else:
            stage2_growth, dbp_flux = step2_max_dbp_uptake(com, args.alpha, biomass_max)
            print(f"\n阶段二：growth ≥ {args.alpha:.2f} * {biomass_max:.6f} = {args.alpha*biomass_max:.6f}")
            print(f"  解得 growth(f*) = {stage2_growth:.6f}")
            print(f"  f* (EX_dbp_m 固定通量) = {dbp_flux:.6f}  （负值=摄入，越负代表摄入越强）")

        # 结果输出路径：若未显式提供 --out，则默认写入 model-dir/Result6_community.csv
        default_out = os.path.join(os.path.abspath(args.model_dir), "Result6_community.csv")
        out_path = os.path.abspath(args.out) if args.out else default_out
        pd.DataFrame([{
            "model_count": len(com.taxa),
            "biomass_max": biomass_max,
            "alpha": args.alpha,
            "dbp_flux_stage2": dbp_flux,
            "notes": ""
        }]).to_csv(out_path, index=False)
        print("\n✅ 已导出结果：", out_path)

    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


if __name__ == "__main__":
    main()
