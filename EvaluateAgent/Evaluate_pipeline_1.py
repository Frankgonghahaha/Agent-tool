# -*- coding: utf-8 -*-
"""
Evaluate_pipeline_1.py — 通用“按表填充代谢反应”的脚本（参数化，支持任意反应）

功能概述
- 从 CSV 读取任意数量的反应（行），为目录下所有 SBML 模型批量添加这些反应；
- 自动创建缺失的代谢物（按 _c/_e/_m/_p 后缀推断舱室；缺省为 c）；
- 可选：自动为方程中出现的物质补充 EX_*_e 与 TRANS_*（e<->c），并为 o2 自动补全 TRANS_o2 与 EX_o2_e；
- 反应边界（LB/UB）可从 CSV 指定；若未指定则按箭头推断：
    - "A + B -> C"（不可逆）默认 lb=0, ub=1000
    - "A + B <-> C"（可逆）默认 lb=-1000, ub=1000

CSV 格式（列名，大小写不敏感）
必需：
- Reaction equation 或 equation： 反应方程，形如：
    "2 glc__D_c + o2_c -> 2 pyr_c + 2 h2o_c"
可选：
- id 或 Reaction id：反应 ID（如缺省，自动生成 RXN_00001, RXN_00002...）
- lb / lower_bound：下界（数值）
- ub / upper_bound：上界（数值）
- name：反应可读名称

示例用法：
python add_pathway.py \
  --model-dir "/path/to/models" \
  --reactions "/path/to/reaction.csv" \
  --out "/path/to/output" \
  --auto-extrans
"""

from __future__ import annotations
import os
import re
import argparse
from typing import Dict, List, Tuple, Optional

import pandas as pd
from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model, write_sbml_model

# ---------------------------
# 工具 & 解析
# ---------------------------

COMP_SUFFIX = {"c", "e", "p", "m"}

ARROW_REV = "<->"
ARROW_IRR = "->"

def discover_models(model_dir: str, only_files: List[str] | None = None) -> List[str]:
    if only_files:
        return [os.path.join(model_dir, f) for f in only_files]
    hits: List[str] = []
    for root, _dirs, files in os.walk(model_dir):
        for fn in files:
            if fn.lower().endswith((".xml", ".sbml")):
                hits.append(os.path.join(root, fn))
    return sorted(hits)


def _norm_colname(c: str) -> str:
    return c.strip().lower().replace(" ", "_")


def load_reaction_table(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    cols = {_norm_colname(c): c for c in df.columns}
    # 找到 equation 列
    eq_col = None
    for cand in ["reaction_equation", "equation"]:
        if cand in cols:
            eq_col = cols[cand]
            break
    if eq_col is None:
        raise ValueError("CSV 需包含列 'Reaction equation' 或 'equation'")
    # 兼容其它列
    id_col = None
    for cand in ["id", "reaction_id"]:
        if cand in cols:
            id_col = cols[cand]
            break
    lb_col = None
    for cand in ["lb", "lower_bound"]:
        if cand in cols:
            lb_col = cols[cand]
            break
    ub_col = None
    for cand in ["ub", "upper_bound"]:
        if cand in cols:
            ub_col = cols[cand]
            break
    name_col = cols.get("name")

    # 标准化输出列
    out = []
    for i, row in df.iterrows():
        eq = str(row[eq_col]).strip()
        if not eq:
            continue
        rid = str(row[id_col]).strip() if id_col and pd.notna(row[id_col]) else ""
        lb = float(row[lb_col]) if lb_col and pd.notna(row[lb_col]) else None
        ub = float(row[ub_col]) if ub_col and pd.notna(row[ub_col]) else None
        rname = str(row[name_col]).strip() if name_col and pd.notna(row[name_col]) else ""
        out.append({"id": rid, "equation": eq, "lb": lb, "ub": ub, "name": rname})
    return pd.DataFrame(out)


def _parse_term(term: str) -> Tuple[float, str]:
    """解析一个项，如 "2 glc__D_c" 或 "o2_c" → (coef, met_id)。"""
    term = term.strip()
    if not term:
        raise ValueError("empty term")
    m = re.match(r"^([+-]?\d*(?:\.\d+)?)\s*([A-Za-z0-9_\.]+)$", term)
    if m:
        coef_s, met = m.groups()
        coef = float(coef_s) if coef_s not in (None, "", "+", "-") else ( -1.0 if coef_s == "-" else 1.0 )
    else:
        # 无系数时，默认1
        coef, met = 1.0, term
    return coef, met


def parse_equation_to_stoich(equation: str) -> Tuple[Dict[str, float], bool]:
    """把字符串方程解析为 {met: coef}，底物系数为负，产物为正；返回 (stoich, reversible)。"""
    eq = equation.replace("=", "").strip()
    reversible = ARROW_REV in eq
    arrow = ARROW_REV if reversible else ARROW_IRR
    if arrow not in eq:
        raise ValueError(f"方程缺少箭头 '->' 或 '<->'：{equation}")
    left, right = [s.strip() for s in eq.split(arrow, 1)]
    stoich: Dict[str, float] = {}
    if left:
        for t in left.split("+"):
            coef, met = _parse_term(t)
            stoich[met] = stoich.get(met, 0.0) - coef
    if right:
        for t in right.split("+"):
            coef, met = _parse_term(t)
            stoich[met] = stoich.get(met, 0.0) + coef
    return stoich, reversible


# ---------------------------
# COBRA 构建
# ---------------------------

def _get_or_make_met(model: Model, met_id: str) -> Metabolite:
    met = model.metabolites.get_by_id(met_id) if met_id in model.metabolites else None
    if met is not None:
        return met
    # 推断 compartment
    comp = "c"
    parts = met_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1] in COMP_SUFFIX:
        comp = parts[1]
    met = Metabolite(met_id, name=met_id, compartment=comp)
    model.add_metabolites([met])
    return met


def _ensure_reaction(model: Model, rxn_id: str, stoich: Dict[str, float], lb: float, ub: float, name: str = "") -> Reaction:
    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        return rxn
    rxn = Reaction(rxn_id)
    rxn.name = name or rxn_id
    rxn.lower_bound = lb
    rxn.upper_bound = ub
    mets: Dict[Metabolite, float] = {}
    for met_id, coef in stoich.items():
        mets[_get_or_make_met(model, met_id)] = float(coef)
    rxn.add_metabolites(mets)
    model.add_reactions([rxn])
    return rxn


def _ensure_exchange(model: Model, base: str, lb: float = -1000.0, ub: float = 1000.0) -> Reaction:
    met_e = f"{base}_e"
    stoich = {met_e: -1.0}
    rid = f"EX_{base}_e"
    return _ensure_reaction(model, rid, stoich, lb, ub, name=rid)


def _ensure_transport_ce(model: Model, base: str, lb: float = -1000.0, ub: float = 1000.0) -> Reaction:
    met_c, met_e = f"{base}_c", f"{base}_e"
    stoich = {met_e: -1.0, met_c: +1.0}  # e -> c 方向为负
    rid = f"TRANS_{base}"
    return _ensure_reaction(model, rid, stoich, lb, ub, name=rid)


def _base_of(met_id: str) -> str:
    parts = met_id.rsplit("_", 1)
    return parts[0] if len(parts) == 2 and parts[1] in COMP_SUFFIX else met_id


# ---------------------------
# 主流程
# ---------------------------

def add_reactions_to_model(model: Model, rows: List[Dict[str, object]], auto_extrans: bool = False) -> List[str]:
    created: List[str] = []
    for i, row in enumerate(rows, start=1):
        rid = str(row.get("id") or "").strip()
        eq = str(row.get("equation") or "").strip()
        lb = row.get("lb")
        ub = row.get("ub")
        rname = str(row.get("name") or "").strip()
        if not eq:
            continue
        stoich, reversible = parse_equation_to_stoich(eq)
        # 默认边界
        if lb is None or ub is None:
            if reversible:
                lb = -1000.0 if lb is None else float(lb)
                ub =  1000.0 if ub is None else float(ub)
            else:
                lb = 0.0 if lb is None else float(lb)
                ub = 1000.0 if ub is None else float(ub)
        else:
            lb = float(lb); ub = float(ub)
        if not rid:
            rid = f"RXN_{i:05d}"
        # 添加反应
        _ensure_reaction(model, rid, stoich, lb, ub, name=rname)
        created.append(rid)
        # 可选：自动补 EX/TRANS
        if auto_extrans:
            bases = { _base_of(mid) for mid in stoich.keys() }
            for b in bases:
                _get_or_make_met(model, f"{b}_c")
                _get_or_make_met(model, f"{b}_e")
                _ensure_transport_ce(model, b, -1000.0, 1000.0)
                _ensure_exchange(model, b, -1000.0, 1000.0)
            # 氧气特别处理
            if any(_base_of(m) == "o2" for m in stoich.keys()):
                _get_or_make_met(model, "o2_c"); _get_or_make_met(model, "o2_e")
                _ensure_transport_ce(model, "o2", -1000.0, 1000.0)
                _ensure_exchange(model, "o2", -1000.0, 1000.0)
    return created


def process_one(in_path: str, out_dir: str, rxn_df: pd.DataFrame, auto_extrans: bool) -> None:
    print(f"\n>>> 处理模型：{in_path}")
    model = read_sbml_model(in_path)
    rows = rxn_df.to_dict(orient="records")
    created = add_reactions_to_model(model, rows, auto_extrans=auto_extrans)
    os.makedirs(out_dir, exist_ok=True)
    name = os.path.splitext(os.path.basename(in_path))[0]
    out_path = os.path.join(out_dir, f"{name}.xml")
    write_sbml_model(model, out_path)
    print(f"  ✅ 已保存模型：{out_path}")
    if created:
        print("  + 本次新增/确保存在的反应：", ", ".join(created))
    else:
        print("  * 未新增任何反应（可能已存在）。")


def main():
    parser = argparse.ArgumentParser(description="通用：按 CSV 批量为 SBML 添加反应，可选自动补 EX/TRANS")
    parser.add_argument("--model-dir", required=True, help="输入模型目录（递归搜索 .xml/.sbml）")
    parser.add_argument("--reactions", required=True, help="反应 CSV（列：Reaction equation / equation；可选 id/lb/ub/name）")
    parser.add_argument("--out", required=True, help="输出目录（写入修改后的 SBML）")
    parser.add_argument("--only", default="", help="仅处理这些文件名（逗号分隔，相对于 --model-dir），可选")
    parser.add_argument("--auto-extrans", action="store_true", help="为方程中出现的物质自动补 EX_*_e 与 TRANS_* 以及 o2 的 EX/TRANS")
    args = parser.parse_args()

    only_files = [s.strip() for s in args.only.split(",") if s.strip()] if args.only else []
    rxn_df = load_reaction_table(args.reactions)

    targets = discover_models(args.model_dir, only_files)
    if not targets:
        print(f"⚠ 未在目录中找到 SBML 模型：{args.model_dir}")
        return

    print(f"将处理 {len(targets)} 个模型；输出目录：{args.out}")
    for fp in targets:
        try:
            process_one(fp, args.out, rxn_df, auto_extrans=args.auto_extrans)
        except Exception as e:
            print(f"  [跳过] {fp} 发生错误：{e}")


if __name__ == "__main__":
    main()
