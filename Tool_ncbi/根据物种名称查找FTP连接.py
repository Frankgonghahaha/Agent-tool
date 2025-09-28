#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hard-coded, verbose version for PyCharm:
- 扫描指定目录递归所有 *.xlsx
- 只读取 sheet = 'Bacteria'
- 读取列 'Strain ID/Microorganism'（兼容若干别名）
- 调 NCBI Assembly 搜索各菌株的 WGS/Assembly
- 评分挑选最佳装配，生成结果 CSV + 一键下载命令
- 打印详细日志，便于观察进度

依赖：
    pip install biopython pandas openpyxl

可选（自动下载命中装配）：
    安装 NCBI datasets CLI 并将 DO_DOWNLOAD=True
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
"""

# ---------------------- CONFIG ----------------------
ROOT_DIR   = r"/Users/frank_gong/文档/生物智能体/智能体数据库/MIBPOP-main/Organism"
OUT_CSV    = r"/Users/frank_gong/文档/生物智能体/wgs_lookup_results.csv"
EMAIL      = "710244143@qq.com"   # <<< 修改为你的真实邮箱
API_KEY    = None                       # 可选：你的 NCBI API key（提高速率限额）
SHEET_NAME = "Bacteria"                 # 你的需求：sheet 改为 Bacteria
DO_DOWNLOAD   = False                   # True 将尝试用 datasets CLI 实时下载
ENTREZ_SLEEP  = 0.34                    # 每次 E-utilities 调用后的休眠（秒）
MAX_RET_PER_Q = 10                      # 每次检索最多拿多少条候选
# ----------------------------------------------------

import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import Entrez
from xml.etree import ElementTree as ET

ASSEMBLY_LEVEL_ORDER = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1
}

def log(msg: str) -> None:
    print(msg, flush=True)

def clean_str(s: Optional[str]) -> str:
    return (s or "").strip()

def species_from_strain_name(name: str) -> str:
    """
    从菌株字符串中提取物种层级关键词：
    - 'Rhodococcus ruber YC-YT1' -> 'Rhodococcus ruber'
    - 'Dehalobacter sp. E2' -> 'Dehalobacter'
    """
    name = clean_str(name)
    parts = name.split()
    if not parts:
        return ""
    if len(parts) >= 2 and parts[1].lower() not in {"sp.", "spp.", "cf.", "aff."}:
        return " ".join(parts[:2])
    return parts[0]

def score_match(query: str, docsum: Dict) -> int:
    """
    给候选装配打分：是否包含原查询字样、是否匹配物种层级、是否 RefSeq、装配等级等
    """
    q = query.lower()
    score = 0
    fields = []
    for key in ["Organism", "SpeciesName", "AssemblyName", "AssemblyAccession", "Title", "SubmitterOrganization"]:
        val = clean_str(docsum.get(key, "")).lower()
        if val:
            fields.append(val)
    subtype = docsum.get("SubType", [])
    subval  = docsum.get("SubName", [])
    for k, v in zip(subtype, subval):
        fields.append(f"{k}:{(v or '').lower()}")
        fields.append((v or "").lower())
    hay = " | ".join(fields)

    if q in hay:
        score += 10
    sp = species_from_strain_name(query).lower()
    if sp and sp in hay:
        score += 5
    if clean_str(docsum.get("FtpPath_RefSeq")):
        score += 2
    level = clean_str(docsum.get("AssemblyStatus"))
    score += ASSEMBLY_LEVEL_ORDER.get(level, 0)
    return score

def choose_best(doc_summaries: List[Dict], query: str) -> Optional[Dict]:
    if not doc_summaries:
        return None
    ranked = sorted(doc_summaries, key=lambda d: score_match(query, d), reverse=True)
    return ranked[0]

def entrez_esearch_assembly(term: str, retmax: int = MAX_RET_PER_Q, email: str = "", api_key: Optional[str] = None) -> List[str]:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.esearch(db="assembly", term=term, sort="relevance", retmax=retmax, rettype="uilist")
    rec = Entrez.read(handle)
    handle.close()
    return rec.get("IdList", [])

def entrez_esummary_assembly(uids: List[str], email: str = "", api_key: Optional[str] = None) -> List[Dict]:
    if not uids:
        return []
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    handle = Entrez.esummary(db="assembly", id=",".join(uids), retmode="xml")
    data = handle.read()
    handle.close()
    root = ET.fromstring(data)
    out = []
    for doc in root.findall(".//DocumentSummary"):
        d = {child.tag: (child.text or "") for child in doc}
        # 把 FTP 字段也带上
        for tag in ["FtpPath_RefSeq", "FtpPath_GenBank"]:
            elem = doc.find(tag)
            d[tag] = elem.text if elem is not None else ""
        # 采集 SubType/SubName 用于识别 strain 等
        subtypes = [e.text for e in doc.findall("SubType/Item")]
        subnames = [e.text for e in doc.findall("SubName/Item")]
        d["SubType"] = subtypes
        d["SubName"] = subnames
        out.append(d)
    return out

def build_search_terms(name: str) -> List[str]:
    """
    从严格到宽松构建多个检索式，优先返回最新装配（latest）
    """
    q = re.sub(r'[\[\]"]+', " ", name).strip()
    terms = [
        f'"{q}"[All Fields] AND latest[filter]',
        f'"{q}"[Organism] AND latest[filter]',
        f'("{q}"[All Fields]) AND (genome repres*[All Fields]) AND latest[filter]',
    ]
    sp = species_from_strain_name(q)
    if sp and sp.lower() != q.lower():
        terms.extend([
            f'"{sp}"[Organism] AND latest[filter]',
            f'"{sp}"[All Fields] AND latest[filter]',
        ])
    return terms

def datasets_cmd_for_accession(acc: str, out_prefix: str = None) -> str:
    safe_prefix = re.sub(r'[^A-Za-z0-9._-]+', "_", out_prefix or acc)
    return f"datasets download genome accession {acc} --include genome,gff3,protein --filename {safe_prefix}.zip"

def try_download_with_datasets(acc: str, out_dir: Path) -> Optional[int]:
    """
    可选：若系统安装了 datasets CLI，则尝试直接下载。
    """
    import shutil, subprocess
    if shutil.which("datasets") is None:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "datasets", "download", "genome", "accession", acc,
        "--include", "genome,gff3,protein",
        "--filename", str(out_dir / f"{acc}.zip")
    ]
    try:
        return subprocess.call(cmd)
    except Exception:
        return None

def best_assembly_for_name(name: str, email: str, api_key: Optional[str], sleep: float = 0.34):
    """
    返回：(status, best_docsum_or_None, all_docs)
    status ∈ {"exact_or_best", "species_fallback", "not_found"}
    """
    terms = build_search_terms(name)
    all_docs: List[Dict] = []
    for i, term in enumerate(terms):
        log(f"  - 尝试检索式 {i+1}/{len(terms)}: {term}")
        try:
            uids = entrez_esearch_assembly(term, retmax=MAX_RET_PER_Q, email=email, api_key=api_key)
            log(f"    · 命中 UID 数：{len(uids)}")
            time.sleep(sleep)
            if not uids:
                continue
            docs = entrez_esummary_assembly(uids, email=email, api_key=api_key)
            log(f"    · 解析候选记录数：{len(docs)}")
            time.sleep(sleep)
            if not docs:
                continue
            all_docs.extend(docs)
            best = choose_best(docs, name)
            if best:
                level = clean_str(best.get("AssemblyStatus"))
                acc   = clean_str(best.get("AssemblyAccession"))
                org   = clean_str(best.get("Organism"))
                log(f"    · 暂选最佳：{org} | {acc} | {level}")
                status = "exact_or_best" if i < 2 else ("species_fallback" if i >= 3 else "exact_or_best")
                return status, best, all_docs
        except Exception as e:
            log(f"    ! 检索异常：{e}")
            time.sleep(1.0)
            continue
    return "not_found", None, all_docs

def main():
    start_ts = time.time()
    log("=== WGS/Assembly 批量检索开始 ===")
    log(f"[配置] ROOT_DIR   = {ROOT_DIR}")
    log(f"[配置] SHEET_NAME = {SHEET_NAME}")
    log(f"[配置] OUT_CSV    = {OUT_CSV}")
    log(f"[配置] EMAIL      = {EMAIL}")
    log(f"[配置] ENTREZ_SLEEP = {ENTREZ_SLEEP}s, MAX_RET_PER_Q = {MAX_RET_PER_Q}")
    if DO_DOWNLOAD:
        log("[配置] DO_DOWNLOAD = True（将尝试调用 'datasets' CLI 下载）")
    else:
        log("[配置] DO_DOWNLOAD = False（仅生成下载命令）")

    root = Path(ROOT_DIR).expanduser()
    if not root.exists():
        log(f"[错误] ROOT_DIR 不存在：{root}")
        sys.exit(1)

    # 1) 枚举 Excel
    xlsx_files = sorted(root.rglob("*.xlsx"))
    log(f"[1/4] 共发现 .xlsx 文件：{len(xlsx_files)} 个")
    if not xlsx_files:
        log("没有找到任何 Excel，结束。")
        sys.exit(0)

    # 2) 读取菌名
    col_candidates = [
        "Strain ID/Microorganism",
        "Strain ID / Microorganism",
        "Strain",
        "Microorganism"
    ]
    all_names: List[str] = []
    files_scanned = 0
    for xf in xlsx_files:
        files_scanned += 1
        try:
            xls = pd.ExcelFile(xf)
            if SHEET_NAME not in xls.sheet_names:
                log(f"  [跳过] {xf.name}: 不含 sheet '{SHEET_NAME}'")
                continue
            df = pd.read_excel(xf, sheet_name=SHEET_NAME)
            col_name = None
            for c in df.columns:
                if c in col_candidates:
                    col_name = c; break
            if col_name is None:
                # case-insensitive fallback
                lower_set = {s.lower() for s in col_candidates}
                for c in df.columns:
                    if c.strip().lower() in lower_set:
                        col_name = c; break
            if col_name is None:
                log(f"  [跳过] {xf.name}: 找不到菌株列（{col_candidates}）")
                continue
            names = df[col_name].dropna().astype(str).map(str.strip).tolist()
            names = [n for n in names if n and n.lower() not in {"nan", "none"}]
            all_names.extend(names)
            log(f"  [读取] {xf.name}: 从列 '{col_name}' 提取到 {len(names)} 个名称")
        except Exception as e:
            log(f"  [警告] 读取失败 {xf.name}: {e}")

    if not all_names:
        log("[2/4] 没有在任何表中找到菌株名称，结束。")
        sys.exit(0)

    unique_names = sorted(set(all_names))
    log(f"[2/4] 去重后共有菌株/物种名称：{len(unique_names)} 个")

    # 3) 检索每个名称
    rows = []
    found, species_fallback, not_found = 0, 0, 0
    for idx, name in enumerate(unique_names, 1):
        log(f"\n[3/4] ({idx}/{len(unique_names)}) 查询：{name}")
        status, best, _ = best_assembly_for_name(name, EMAIL, API_KEY, sleep=ENTREZ_SLEEP)
        if best:
            acc   = clean_str(best.get("AssemblyAccession"))
            org   = clean_str(best.get("Organism"))
            title = clean_str(best.get("Title"))
            level = clean_str(best.get("AssemblyStatus"))
            refseq_ftp  = clean_str(best.get("FtpPath_RefSeq"))
            genbank_ftp = clean_str(best.get("FtpPath_GenBank"))
            bioproject  = clean_str(best.get("BioProjectAccn"))
            biosample   = clean_str(best.get("BioSampleAccn"))
            cmd         = datasets_cmd_for_accession(acc, out_prefix=f"{org}_{acc}")

            log(f"    · 结果：{status} | {org} | {acc} | {level}")
            if refseq_ftp:
                log(f"    · RefSeq FTP：{refseq_ftp}")
            elif genbank_ftp:
                log(f"    · GenBank FTP：{genbank_ftp}")
            else:
                log("    · 无 FTP 路径（可能下线或未归档）")

            dl_code = ""
            if DO_DOWNLOAD and acc:
                code = try_download_with_datasets(acc, Path("datasets_downloads"))
                dl_code = f"datasets_exit={code}"
                log(f"    · 尝试下载：exit={code}")

            rows.append({
                "query_name": name,
                "match_status": status,
                "assembly_accession": acc,
                "organism_name": org,
                "title": title,
                "assembly_level": level,
                "bioproject": bioproject,
                "biosample": biosample,
                "ftp_refseq": refseq_ftp,
                "ftp_genbank": genbank_ftp,
                "datasets_cmd": cmd,
                "notes": dl_code
            })

            if status == "exact_or_best":
                found += 1
            elif status == "species_fallback":
                species_fallback += 1
        else:
            log("    · 未找到装配记录")
            rows.append({
                "query_name": name,
                "match_status": "not_found",
                "assembly_accession": "",
                "organism_name": "",
                "title": "",
                "assembly_level": "",
                "bioproject": "",
                "biosample": "",
                "ftp_refseq": "",
                "ftp_genbank": "",
                "datasets_cmd": "",
                "notes": ""
            })
            not_found += 1

        time.sleep(ENTREZ_SLEEP)

    # 4) 写出结果
    out_path = Path(OUT_CSV).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False, encoding="utf-8-sig")

    elapsed = time.time() - start_ts
    log("\n[4/4] 写出结果完成")
    log(f"      路径：{out_path}")
    log(f"      合计：{len(unique_names)} 个查询 | 命中（exact_or_best）：{found} | 物种回退：{species_fallback} | 未找到：{not_found}")
    log(f"      用时：{elapsed:.1f} s")
    log("=== 全部完成 ===")

if __name__ == "__main__":
    main()