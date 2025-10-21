#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import itertools
import math
import numpy as np
import pandas as pd
from collections import defaultdict

# ========== 工具函数 ==========
def normalize_name(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    s = s.replace("（","(").replace("）",")").replace("，",",").replace("、",",")
    s = " ".join(s.split())
    return s

def safe_mean(vals):
    """对非空可迭代求均值；若为空返回0."""
    vals = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return float(np.mean(vals)) if vals else 0.0

def read_scores(path_scores):
    df = pd.read_csv(path_scores)
    # 必要列: species, S_microbe
    need = ["species", "S_microbe"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"[scores] 缺少必要列: {missing}, 现有: {list(df.columns)}")

    # 可选列默认处理
    if "environment_match" not in df.columns:
        df["environment_match"] = 1.0
    if "source" not in df.columns:
        df["source"] = "unknown"

    # 类型清洗
    df["species"] = df["species"].map(normalize_name)
    df["S_microbe"] = pd.to_numeric(df["S_microbe"], errors="coerce").fillna(0)
    if "kcat_max" in df.columns:
        df["kcat_max"] = pd.to_numeric(df["kcat_max"], errors="coerce")

    df = df.dropna(subset=["species"]).reset_index(drop=True)
    return df

def read_pairs(path_pairs):
    df = pd.read_csv(path_pairs)
    # 需要列: functional_species, complement_species, competition_index, complementarity_index, delta_index
    need = ["functional_species", "complement_species", "competition_index", "complementarity_index", "delta_index"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"[pairs] 缺少必要列: {missing}, 现有: {list(df.columns)}")
    df["functional_species"] = df["functional_species"].map(normalize_name)
    df["complement_species"] = df["complement_species"].map(normalize_name)
    for c in ["competition_index","complementarity_index","delta_index"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def build_pair_lookup(df_pairs):
    """构建 pair -> (comp, compl, delta) 的查找。对称存储。"""
    pair2vals = {}
    for _, r in df_pairs.iterrows():
        a = r["functional_species"]; b = r["complement_species"]
        if a is None or b is None:
            continue
        comp = r["competition_index"]; compl = r["complementarity_index"]; delta = r["delta_index"]
        pair2vals[(a,b)] = (comp, compl, delta)
        pair2vals[(b,a)] = (comp, compl, delta)
    return pair2vals

def aggregate_pair_metrics(members, pair2vals):
    """对一个组合，聚合 pairwise 指标，返回 (avg_delta_pos, avg_comp_pos, used_pairs)"""
    deltas_pos = []
    comps_pos  = []
    used = 0
    for i in range(len(members)):
        for j in range(i+1, len(members)):
            a, b = members[i], members[j]
            vals = pair2vals.get((a,b))
            if not vals:
                continue
            comp, compl, delta = vals
            if delta is not None and not (isinstance(delta, float) and math.isnan(delta)) and delta > 0:
                deltas_pos.append(float(delta))
            if comp is not None and not (isinstance(comp, float) and math.isnan(comp)) and comp > 0:
                comps_pos.append(float(comp))
            used += 1
    return safe_mean(deltas_pos), safe_mean(comps_pos), used

def average_S_microbe(members, score_df):
    vals = []
    for s in members:
        row = score_df.loc[score_df["species"] == s]
        v = float(row.iloc[0]["S_microbe"]) if not row.empty else 0.0
        vals.append(v)
    return safe_mean(vals)

def average_kcat(members, score_df):
    """从 scores 表中读取成员的 kcat_max 求平均；缺失或非数值会自动跳过。"""
    vals = []
    for s in members:
        row = score_df.loc[score_df["species"] == s]
        if row.empty:
            continue
        v = row.iloc[0].get("kcat_max", np.nan)
        try:
            v = float(v)
        except Exception:
            v = np.nan
        if not (isinstance(v, float) and math.isnan(v)):
            vals.append(v)
    return safe_mean(vals)

def count_source(members, score_df):
    cnt = defaultdict(int)
    for s in members:
        row = score_df.loc[score_df["species"] == s]
        src = row.iloc[0].get("source", "unknown") if not row.empty else "unknown"
        cnt[str(src)] += 1
    return dict(cnt)

def calc_S_consort(members, score_df, pair2vals, alpha, beta, gamma, lambd, mu):
    avg_micro = average_S_microbe(members, score_df)
    avg_delta_pos, avg_comp_pos, used_pairs = aggregate_pair_metrics(members, pair2vals)
    avg_kcat = average_kcat(members, score_df)
    N = len(members)
    score = (
        alpha * avg_micro
        + beta  * avg_delta_pos
        - gamma * avg_comp_pos
        + lambd * avg_kcat
        - mu    * N
    )
    detail = {
        "avg_S_microbe": avg_micro,
        "avg_delta_pos": avg_delta_pos,
        "avg_comp_pos":  avg_comp_pos,
        "avg_kcat": avg_kcat,
        "size": N,
        "used_pairs": used_pairs
    }
    return float(score), detail

# ========== 搜索策略 ==========
def greedy_search(species_list, score_df, pair2vals, kmin, kmax, alpha, beta, gamma, lambd, mu, require_functional=True, topK=20):
    """从最高 S_microbe 开始逐步加点，记录每个规模的最优组合。"""
    # 以 S_microbe 排序
    s2score = {r["species"]: float(r["S_microbe"]) for _, r in score_df.iterrows()}
    species_sorted = sorted([s for s in species_list if s in s2score], key=lambda s: s2score[s], reverse=True)

    best_by_k = {}
    current = []

    # 至少包含一个 functional 的约束检查
    def ok_combo(members):
        if not require_functional:
            return True
        rows = score_df[score_df["species"].isin(members)]
        return (rows["source"] == "functional").any()

    for s in species_sorted:
        # 从单个种子开始，往上扩展到 kmax
        current = [s]
        if kmin <= 1 and ok_combo(current):
            score, detail = calc_S_consort(current, score_df, pair2vals, alpha, beta, gamma, lambd, mu)
            best_by_k[1] = (score, current[:], detail)

        for _ in range(1, kmax):
            # 在剩余物种里找使得增益最大的那个
            candidates = [x for x in species_sorted if x not in current]
            best_gain = None
            best_next = None
            best_detail_tmp = None
            for c in candidates:
                trial = current + [c]
                if not ok_combo(trial):
                    continue
                sc, det = calc_S_consort(trial, score_df, pair2vals, alpha, beta, gamma, lambd, mu)
                if best_gain is None or sc > best_gain:
                    best_gain = sc
                    best_next = c
                    best_detail_tmp = det
            if best_next is None:
                break
            current.append(best_next)
            k = len(current)
            if k >= kmin:
                best_prev = best_by_k.get(k)
                if best_prev is None or best_gain > best_prev[0]:
                    best_by_k[k] = (best_gain, current[:], best_detail_tmp)

    # 整理输出 TopK
    rows = []
    for k in sorted(best_by_k.keys()):
        sc, members, det = best_by_k[k]
        src_cnt = count_source(members, score_df)
        rows.append({
            "consortium_id": f"greedy_k{k}",
            "members": ";".join(members),
            "size": k,
            "avg_S_microbe": det["avg_S_microbe"],
            "avg_delta_pos": det["avg_delta_pos"],
            "avg_comp_pos": det["avg_comp_pos"],
            "avg_kcat": det["avg_kcat"],
            "S_consort": sc,
            "used_pairs": det["used_pairs"],
            "source_count": src_cnt
        })
    rows = sorted(rows, key=lambda r: r["S_consort"], reverse=True)[:topK]
    return pd.DataFrame(rows)

def exhaustive_search(species_list, score_df, pair2vals, kmin, kmax, alpha, beta, gamma, lambd, mu, require_functional=True, topK=50, hard_cap=500000):
    """全组合搜索（注意组合数可能很大，提供 hard_cap 保护）"""
    rows = []
    total_checked = 0

    # 功能约束
    def ok_combo(members):
        if not require_functional:
            return True
        rows_ = score_df[score_df["species"].isin(members)]
        return (rows_["source"] == "functional").any()

    for k in range(kmin, kmax+1):
        for combo in itertools.combinations(species_list, k):
            if total_checked > hard_cap:
                break
            total_checked += 1
            combo = list(combo)
            if not ok_combo(combo):
                continue
            sc, det = calc_S_consort(combo, score_df, pair2vals, alpha, beta, gamma, lambd, mu)
            src_cnt = count_source(combo, score_df)
            rows.append({
                "consortium_id": f"exhaustive_k{k}_{total_checked}",
                "members": ";".join(combo),
                "size": k,
                "avg_S_microbe": det["avg_S_microbe"],
                "avg_delta_pos": det["avg_delta_pos"],
                "avg_comp_pos": det["avg_comp_pos"],
                "avg_kcat": det["avg_kcat"],
                "S_consort": sc,
                "used_pairs": det["used_pairs"],
                "source_count": src_cnt
            })
    if not rows:
        return pd.DataFrame(columns=[
            "consortium_id","members","size","avg_S_microbe",
            "avg_delta_pos","avg_comp_pos","avg_kcat",
            "S_consort","used_pairs","source_count"
        ])
    rows = sorted(rows, key=lambda r: r["S_consort"], reverse=True)[:topK]
    return pd.DataFrame(rows)

# ========== 主流程 ==========
def main():
    ap = argparse.ArgumentParser(description="Design_pipeline_2: 菌群组合优化（基于 S_microbe + pairwise 互作 + 平均kcat）")
    ap.add_argument("--scores", required=True, help="Result2_candidate_scores.csv 路径")
    ap.add_argument("--pairs", required=True, help="Result3_pair_Com_index.csv 路径")
    ap.add_argument("--out", default=None, help="输出 CSV（默认 Result4_optimal_consortia.csv）")
    ap.add_argument("--members_out", default=None, help="成员贡献排名输出 CSV（默认 Result4_members_rank.csv）")

    # 物种池控制
    ap.add_argument("--topN", type=int, default=50, help="按 S_microbe 取前 N 个物种作为搜索池")
    ap.add_argument("--kmin", type=int, default=2, help="组合最小大小")
    ap.add_argument("--kmax", type=int, default=5, help="组合最大大小")
    ap.add_argument("--mode", choices=["greedy","exhaustive"], default="greedy", help="搜索模式")
    ap.add_argument("--topK", type=int, default=20, help="输出前 K 个最优组合")
    ap.add_argument("--require_functional", type=int, default=1, help="是否强制组合包含至少一个功能菌（1/0）")

    # 目标函数权重
    ap.add_argument("--alpha", type=float, default=0.2, help="权重: 平均 S_microbe")
    ap.add_argument("--beta",  type=float, default=0.2, help="权重: 平均正 Delta")
    ap.add_argument("--gamma", type=float, default=0.1, help="权重: 平均正 Competition（惩罚）")
    ap.add_argument("--lambda_",type=float, default=0.35, help="权重: 平均 kcat_max（加分）")
    ap.add_argument("--mu",    type=float, default=-0.05, help="权重: 规模 N（惩罚）")

    args = ap.parse_args()

    out_csv = args.out or "Result4_optimal_consortia.csv"
    out_members = args.members_out or "Result5_members_rank.csv"

    # 读取
    scores = read_scores(args.scores)
    pairs  = read_pairs(args.pairs)

    # 物种池：按 S_microbe 取 TopN
    scores = scores.sort_values(by="S_microbe", ascending=False).reset_index(drop=True)
    species_pool = scores["species"].tolist()[:args.topN]

    # pair lookup
    pair2vals = build_pair_lookup(pairs)

    # 搜索
    if args.mode == "greedy":
        df_best = greedy_search(
            species_pool, scores, pair2vals,
            args.kmin, args.kmax,
            args.alpha, args.beta, args.gamma, args.lambda_, args.mu,
            require_functional=bool(args.require_functional),
            topK=args.topK
        )
    else:
        df_best = exhaustive_search(
            species_pool, scores, pair2vals,
            args.kmin, args.kmax,
            args.alpha, args.beta, args.gamma, args.lambda_, args.mu,
            require_functional=bool(args.require_functional),
            topK=args.topK
        )

    if df_best.empty:
        print("⚠️ 未找到满足条件的组合。请检查参数或扩大物种池/topN。")
        return

    # 保存最优组合
    df_best.to_csv(out_csv, index=False, encoding="utf-8")
    print(f"✅ 已保存最优组合到: {out_csv}")
    print(df_best.head(10).to_string(index=False))

    # === 统计成员在最优解中的“出现频率 & 贡献” ===
    # 用简单启发：成员出现频次 × 成员的 S_microbe
    best_members = []
    for _, r in df_best.iterrows():
        mems = [m.strip() for m in str(r["members"]).split(";") if m.strip()]
        best_members.extend(mems)
    freq = pd.Series(best_members).value_counts().rename("freq").to_frame()
    s_map = scores.set_index("species")["S_microbe"]
    src_map = scores.set_index("species")["source"]
    env_map = scores.set_index("species")["environment_match"] if "environment_match" in scores.columns else None
    rank_df = freq.join(s_map, how="left").join(src_map, how="left")
    if env_map is not None:
        rank_df = rank_df.join(env_map, how="left")
    rank_df["freq_weighted_score"] = rank_df["freq"] * rank_df["S_microbe"].fillna(0)
    rank_df = rank_df.reset_index().rename(columns={"index":"species"})
    rank_df = rank_df.sort_values(by=["freq_weighted_score","freq","S_microbe"], ascending=[False, False, False])

    rank_df.to_csv(out_members, index=False, encoding="utf-8")
    print(f"✅ 已保存成员贡献排序到: {out_members}")
    print(rank_df.head(15).to_string(index=False))


if __name__ == "__main__":
    main()