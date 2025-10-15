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
from typing import List, Tuple, Optional, Any, Iterable
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
EPSILON = 1     # 最小要求的摄入强度 ε（f ≤ -ε）
BIS_TOL = 1e-3     # 二分终止阈值
MAX_ITER = 30      # 二分最大迭代
MU_TOL = 1e-4      # 目标增长与阈值的接近度容差
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


def build_community_from_dir(models_dir: str) -> Tuple[Community, str, pd.DataFrame, List[str]]:
    """把目录下所有 SBML 作为成员构建一个社区。返回 (community, tmp_dir, taxa DataFrame, member_names)。"""
    model_paths = discover_models(models_dir)
    if not model_paths:
        raise FileNotFoundError(f"未在目录下发现 SBML：{models_dir}")

    tmpd = tempfile.mkdtemp(prefix="community_from_dir_")
    try:
        rows = []
        member_names: List[str] = []
        for f in model_paths:
            name = os.path.splitext(os.path.basename(f))[0]
            member_names.append(name)
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
        return com, tmpd, tax, member_names
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




# ---------- 汇总 DBP_HYDRO_BTOH 通量的辅助函数 ----------

def _extract_flux_df(sol: Any) -> Optional[pd.DataFrame]:
    """尽量从 MICOM 解对象中拿到 fluxes DataFrame。不同版本属性名可能不同。"""
    if sol is None:
        return None
    fx = getattr(sol, "fluxes", None)
    if isinstance(fx, pd.DataFrame):
        return fx
    # 一些版本可能挂在 solution 或 results 下
    for attr in ("solution", "results"):
        obj = getattr(sol, attr, None)
        if obj is not None:
            fx2 = getattr(obj, "fluxes", None)
            if isinstance(fx2, pd.DataFrame):
                return fx2
    return None

def _sum_flux_pattern(fluxes: Optional[pd.DataFrame], pattern: str) -> float:
    """在 MICOM fluxes 矩阵中，汇总列名包含 pattern 的通量总和（排除 'medium' 行）。
    若拿不到 fluxes，返回 0.0。"""
    if fluxes is None or not isinstance(fluxes, pd.DataFrame) or fluxes.empty:
        return 0.0
    cols = [c for c in fluxes.columns if pattern.lower() in str(c).lower()]
    if not cols:
        return 0.0
    rows = [r for r in fluxes.index if str(r) != "medium"]
    try:
        sub = fluxes.loc[rows, cols]
        sub = pd.DataFrame(sub)  # 确保是 DataFrame
        vals = pd.to_numeric(sub.values.ravel(), errors="coerce")
        vals = pd.Series(vals).fillna(0.0)
        return float(vals.sum())
    except Exception:
        return 0.0



# ---- 新增成员生长表辅助 ----
def _extract_members_df(sol: Any) -> Optional[pd.DataFrame]:
    for attr in ("members",):
        df = getattr(sol, attr, None)
        if isinstance(df, pd.DataFrame):
            return df.copy()
    for parent in ("solution", "results"):
        obj = getattr(sol, parent, None)
        if obj is not None:
            df = getattr(obj, "members", None)
            if isinstance(df, pd.DataFrame):
                return df.copy()
    return None

def _members_growth_table(sol: Any) -> Optional[pd.DataFrame]:
    df = _extract_members_df(sol)
    if df is None or df.empty:
        return None
    # 去掉合成的 'medium' 行（有些 MICOM 版本会在 members 里包含 medium 聚合行）
    try:
        mask_idx = ~df.index.astype(str).str.lower().eq("medium")
    except Exception:
        mask_idx = pd.Series([True] * len(df), index=df.index)
    if "id" in df.columns:
        try:
            mask_id = ~df["id"].astype(str).str.lower().eq("medium")
        except Exception:
            mask_id = True
        df = df[mask_idx & (mask_id if isinstance(mask_id, pd.Series) else mask_id)]
    else:
        df = df[mask_idx]
    if df.empty:
        return None
    # 兼容不同列名，优先使用包含 'growth' 的列
    growth_col = None
    for c in df.columns:
        if "growth" in str(c).lower():
            growth_col = c
            break
    if growth_col is None:
        return None
    out = pd.DataFrame({
        "member": df.index.astype(str) if df.index.name else df.get("id", df.index).astype(str),
        "growth_rate": pd.to_numeric(df[growth_col], errors="coerce")
    })
    out = out.fillna(0.0)
    # 终极清洗：去掉 member 为 'medium' 的行（忽略大小写与首尾空白）
    try:
        out = out[~out["member"].astype(str).str.strip().str.lower().eq("medium")]
    except Exception:
        pass
    return out

def step2_max_dbp_uptake(comm: Community, alpha: float, biomass_max: float) -> Tuple[float, float, Optional[pd.DataFrame]]:
    """
    二阶段（改进版）：
      - 目标：使增长 μ*(f) 尽量贴近阈值 target=α·μ*（且 μ*(f)≥target），同时 f 尽可能负（最大化 DBP 摄入）。
      - 实现：对 f 做二分，每次固定 EX_dbp_m=f 计算 μ*(f)；
        * 若 μ*(f) ≥ target，则记录为候选并收缩右端（逼近最负仍可行的边界→使 μ*(f) 向 target 贴近）；
        * 若 μ*(f) < target，则提升左端（减少负载）。
      - 返回： (growth_at_f*, f_star, members_growth_df)
    """
    if DBP_EX_ID not in comm.reactions:
        print(f"[警告] 社区中不存在 {DBP_EX_ID}，无法度量 DBP 摄入。")
        return float("nan"), float("nan"), None

    # 先算不加增长约束时的最小摄入 f_min
    f_min = unconstrained_min_ex(comm, DBP_EX_ID)
    print(f"[阶段二初始化] 不加增长约束时的最小 {DBP_EX_ID} = {f_min:.6f}")

    target = alpha * biomass_max
    print(f"[阶段二目标] 需要满足增长 μ_comm ≥ α·μ* = {alpha:.3f} × {biomass_max:.6f} = {target:.6f}")

    # 右端点：至少要达到 -ε 的摄入
    f_hi = -float(EPSILON)
    mu_hi = ct_max_growth_under_ex(comm, DBP_EX_ID, f_hi)
    print(f"[可行性判定] f_hi={f_hi:.6f} → μ*(f_hi)={mu_hi:.6f}")

    # 情况 A：右端点不可行 → 在 [0, f_hi] 内找“最轻负载”的可行负 f
    if mu_hi < target - 1e-12:
        left, right = 0.0, f_hi
        mu_left = ct_max_growth_under_ex(comm, DBP_EX_ID, left)
        print(f"[右端点不可行] 尝试在 [0, {f_hi:.6f}] 内二分；μ*(0)={mu_left:.6f}")
        if mu_left < target - 1e-12:
            print("[失败] 连 f=0 都不满足增长约束，无法给出负摄入。")
            # 返回 f=0 处解
            with fixed_bound(comm, DBP_EX_ID, 0.0):
                sol = comm.cooperative_tradeoff(fraction=1.0, fluxes=True)
                g = float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
                mg = _members_growth_table(sol)
                return g, 0.0, mg
        best_f, best_mu = left, mu_left
        for _ in range(MAX_ITER):
            mid = 0.5 * (left + right)
            mu_mid = ct_max_growth_under_ex(comm, DBP_EX_ID, mid)
            if mu_mid >= target:
                best_f, best_mu = mid, mu_mid
                right = mid
            else:
                left = mid
            if abs(right - left) <= BIS_TOL or abs(mu_mid - target) <= MU_TOL:
                break
        with fixed_bound(comm, DBP_EX_ID, best_f):
            sol = comm.cooperative_tradeoff(fraction=1.0, fluxes=True)
            g = float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
            mg = _members_growth_table(sol)
            return g, best_f, mg

    # 情况 B：右端点可行 → 在 [f_min, f_hi] 搜索“最负仍可行”，并尽量贴近 target
    f_left = min(f_min, f_hi)
    f_right = max(f_min, f_hi)
    mu_left = ct_max_growth_under_ex(comm, DBP_EX_ID, f_left)
    print(f"[初始化区间] f_left={f_left:.6f} → μ*(f_left)={mu_left:.6f}")

    if mu_left < target - 1e-12:
        left, right = f_left, f_hi
        best_f, best_mu = f_hi, mu_hi
        for _ in range(MAX_ITER):
            mid = 0.5 * (left + right)
            mu_mid = ct_max_growth_under_ex(comm, DBP_EX_ID, mid)
            if mu_mid >= target:
                best_f, best_mu = mid, mu_mid
                right = mid  # 可行→收缩右端，逼近阈值
            else:
                left = mid
            if abs(right - left) <= BIS_TOL or abs(mu_mid - target) <= MU_TOL:
                break
        with fixed_bound(comm, DBP_EX_ID, best_f):
            sol = comm.cooperative_tradeoff(fraction=1.0, fluxes=True)
            g = float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
            mg = _members_growth_table(sol)
            return g, best_f, mg

    # 若 f_left 本身也可行，整个区间都可行 → 取最负 f_left
    with fixed_bound(comm, DBP_EX_ID, f_left):
        sol = comm.cooperative_tradeoff(fraction=1.0, fluxes=True)
        g = float(getattr(sol, "growth_rate", getattr(sol, "objective_value", 0.0)) or 0.0)
        mg = _members_growth_table(sol)
        return g, f_left, mg


# ---------- 主程序 ----------
# ---------- α扫描相关辅助 ----------
def _parse_alphas(s: Optional[str]) -> List[float]:
    if not s:
        return [0.5, 0.6, 0.7, 0.8, 0.9]
    vals: List[float] = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            v = float(tok)
            if v <= 0 or v > 1.0:
                continue
            vals.append(v)
        except Exception:
            continue
    return sorted(set(vals)) or [0.5, 0.6, 0.7, 0.8, 0.9]


def run_alpha_scan(comm: Community, biomass_max: float, alphas: Iterable[float]) -> pd.DataFrame:
    rows = []
    for a in alphas:
        g, f, _ = step2_max_dbp_uptake(comm, float(a), biomass_max)
        rows.append({
            "alpha": float(a),
            "biomass_max": biomass_max,
            "growth_at_f": g,
            "dbp_flux_f": f
        })
        print(f"  [α-scan] alpha={a:.3f} → growth(f)={g:.6f}, f*(EX_dbp_m)={f:.6f}")
    return pd.DataFrame(rows)

# ---------- 群落成员名称提取辅助 ----------
def _community_member_names(comm: Community) -> List[str]:
    """尽量稳健地拿到群落成员名称：优先 'id' 列；其次 index；再退而求其次用 'file' 列的基名去掉扩展名。"""
    try:
        taxa = comm.taxa
        if isinstance(taxa, pd.DataFrame) and not taxa.empty:
            # 1) 优先用 'id' 列
            if 'id' in taxa.columns:
                vals = [str(x).strip() for x in taxa['id'].tolist()]
                vals = [v for v in vals if v]
                if vals:
                    return vals
            # 2) 再用 index
            try:
                idx_vals = [str(x).strip() for x in taxa.index.tolist()]
                idx_vals = [v for v in idx_vals if v]
                # 排除纯数字索引的情况（如 0,1,2...）
                if idx_vals and not all(v.isdigit() for v in idx_vals):
                    return idx_vals
            except Exception:
                pass
            # 3) 最后用 'file' 列推断
            if 'file' in taxa.columns:
                base = []
                for p in taxa['file'].tolist():
                    try:
                        b = os.path.splitext(os.path.basename(str(p)))[0]
                        if b:
                            base.append(b)
                    except Exception:
                        continue
                if base:
                    return base
    except Exception:
        pass
    return []

def main():
    parser = argparse.ArgumentParser(description="Pipeline 2：推荐培养基 + 两步优化（最大生长 → 在α下最大化 DBP 摄入）")
    parser.add_argument("--model-dir", required=True, help="代谢模型目录（递归搜索 .xml/.sbml）")
    parser.add_argument("--medium-csv", required=True, help="培养基 CSV（两列：reaction, flux 或 suggested_upper_bound）")
    parser.add_argument("--alpha", type=float, default=0.7, help="第二步生长下界系数 alpha（默认 0.7）")
    parser.add_argument("--out", help="导出结果 CSV 路径（默认写入 model-dir/Result6_community.csv）")
    parser.add_argument("--alpha-scan", action="store_true", help="启用 α 扫描（默认扫描 0.5,0.6,0.7,0.8,0.9）")
    parser.add_argument("--alphas", type=str, help="自定义 α 列表，逗号分隔，如: 0.5,0.6,0.7")
    parser.add_argument("--Robust", action="store_true", help="启用稳健性分析：顺序去除每个物种并重新优化")
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
    com, tmpd, tax_base, member_names = build_community_from_dir(args.model_dir)
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

        # 可选：α 扫描
        if biomass_max > 1e-12 and args.alpha_scan:
            scan_list = _parse_alphas(args.alphas)
            print("\n[α 扫描] 将在以下 α 上运行二阶段：", ", ".join(f"{x:.2f}" for x in scan_list))
            df_scan = run_alpha_scan(com, biomass_max, scan_list)
            scan_out = os.path.join(os.path.abspath(args.model_dir), "Result6_community_alpha_scan.csv")
            df_scan.to_csv(scan_out, index=False)
            print("[α 扫描] 已导出：", scan_out)

        # 第二步：在 alpha*Biomass_max 下最小化 EX_dbp_m
        if biomass_max <= 1e-12:
            print("\n[警告] Biomass_max 为 0，跳过阶段二优化。")
            stage2_growth, dbp_flux, stage2_members_df = float("nan"), float("nan"), None
        else:
            stage2_growth, dbp_flux, stage2_members_df = step2_max_dbp_uptake(com, args.alpha, biomass_max)
            # 钳制输出，防止数值抖动导致增长率看似超过最大值
            biomass_stage2 = min(stage2_growth, biomass_max + 1e-6)
            # 钳制输出，防止数值抖动导致增长率看似超过最大值
            growth_report = biomass_stage2
            print(f"\n阶段二：growth ≥ {args.alpha:.2f} * {biomass_max:.6f} = {args.alpha*biomass_max:.6f}")
            print(f"  解得 growth(f*) = {growth_report:.6f}")
            print(f"  f* (EX_dbp_m 固定通量) = {dbp_flux:.6f}  （负值=摄入，越负代表摄入越强）")

        # 结果输出路径：写入 Excel，多 sheet
        default_out_xlsx = os.path.join(os.path.abspath(args.model_dir), "Result6_community.xlsx")
        out_path = os.path.abspath(args.out) if args.out else default_out_xlsx
        if not out_path.lower().endswith(".xlsx"):
            root, _ext = os.path.splitext(out_path)
            out_path = root + ".xlsx"

        # 优先从 MICOM 的 taxa 表中读取 id 作为群落组成；若失败再回退到文件名列表
        try:
            ids = com.taxa["id"].astype(str).tolist()
        except Exception:
            try:
                ids = [str(x) for x in com.taxa.index.tolist()]
            except Exception:
                ids = list(member_names)
        community_str = ",".join(ids)

        # 若阶段二未运行，Biomass 使用 nan
        try:
            biomass_stage2
        except NameError:
            biomass_stage2 = float('nan')

        summary_df = pd.DataFrame([{
            "model_count": len(com.taxa),
            "Communitity": community_str,
            "Biomass": biomass_stage2,
            "DBP flux": dbp_flux,
        }])

        # 第二张表：成员生长（若可用）
        mg_df = None
        try:
            mg_df = stage2_members_df if isinstance(stage2_members_df, pd.DataFrame) else None
        except NameError:
            mg_df = None

        # ---------- Robustness analysis (rebuild from tax_base each time) ----------
        robust_results = []
        if args.Robust:
            print("\n[Robustness] 启用顺序去除每个物种并重新优化：")
            if not isinstance(tax_base, pd.DataFrame) or tax_base.empty or "id" not in tax_base.columns:
                print("[错误] 无法执行 Robust：tax_base 不是有效的 DataFrame 或缺少 'id' 列。")
            else:
                taxa_ids = tax_base["id"].astype(str).tolist()
                for remove_id in taxa_ids:
                    print(f"\n[Robustness] 移除物种: {remove_id}")
                    # 从原始 taxa 基础表移除该成员，确保始终是 DataFrame
                    taxa_new = tax_base[tax_base["id"].astype(str) != str(remove_id)].copy()
                    if taxa_new.empty:
                        print("  [Robustness] 移除后社区为空，跳过。")
                        robust_results.append({
                            "Removal species": remove_id,
                            "DBP flux": float("nan"),
                            "Biomass": float("nan")
                        })
                        continue
                    try:
                        com_new = Community(taxa_new, name="COMM_DBP_Robust")
                        apply_medium_via_micom(com_new, medium)
                        biomass_max_rb = step1_max_growth(com_new)
                        if biomass_max_rb <= 1e-12:
                            robust_results.append({
                                "Removal species": remove_id,
                                "DBP flux": float("nan"),
                                "Biomass": float("nan")
                            })
                            print("  [Robustness] Biomass_max=0, 跳过阶段二。")
                            continue
                        stage2_growth_rb, dbp_flux_rb, _ = step2_max_dbp_uptake(com_new, args.alpha, biomass_max_rb)
                        biomass_stage2_rb = min(stage2_growth_rb, biomass_max_rb + 1e-6)
                        robust_results.append({
                            "Removal species": remove_id,
                            "DBP flux": dbp_flux_rb,
                            "Biomass": biomass_stage2_rb
                        })
                        print(f"  [Robustness] Biomass={biomass_stage2_rb:.6f}, DBP flux={dbp_flux_rb:.6f}")
                    except Exception as e:
                        robust_results.append({
                            "Removal species": remove_id,
                            "DBP flux": float("nan"),
                            "Biomass": float("nan")
                        })
                        print(f"  [Robustness] 发生错误: {e}")

        with pd.ExcelWriter(out_path) as writer:
            summary_df.to_excel(writer, index=False, sheet_name="Simulate result")
            if mg_df is not None and not mg_df.empty:
                try:
                    mg_df = mg_df[~mg_df["member"].astype(str).str.strip().str.lower().eq("medium")]
                except Exception:
                    pass
                mg_df.to_excel(writer, index=False, sheet_name="Microbial growth")
            else:
                # 若没有成员表，写一个空结构占位
                pd.DataFrame(columns=["member", "growth_rate"]).to_excel(writer, index=False, sheet_name="Microbial growth")
            # 写入 Robust sheet
            if args.Robust:
                robust_df = pd.DataFrame(robust_results, columns=["Removal species", "DBP flux", "Biomass"])
                robust_df.to_excel(writer, index=False, sheet_name="Robust")

        print("\n✅ 已导出结果：", out_path)

    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


if __name__ == "__main__":
    main()
