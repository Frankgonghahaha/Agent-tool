#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量下载 NCBI 全基因组 FASTA（仅使用 Excel 提供的 ftp_refseq 列）：
- 输入 Excel: /Users/frank_gong/文档/生物智能体/智能体数据库/降解功能微生物全基因组/NCBI查询到有全基因组.xlsx
  需要列: 'assembly_accession', 'query_name', 'ftp_refseq'
- 输出目录: /Users/frank_gong/文档/生物智能体/智能体数据库/降解功能微生物全基因组
  输出命名: {query_name}.fasta
"""

import io
import os
import sys
import time
import gzip
import logging
import pathlib
from typing import Optional

from urllib.parse import urlsplit, urlunsplit
import pandas as pd
import requests

# -------- 你的路径配置 --------
EXCEL_PATH = "/Users/frank_gong/文档/生物智能体/智能体数据库/降解功能微生物全基因组/NCBI查询到有全基因组.xlsx"
OUT_DIR    = "/Users/frank_gong/文档/生物智能体/智能体数据库/降解功能微生物全基因组"

# 列名（如你的表头不同，可在这里改）
COL_ACC  = "assembly_accession"
COL_NAME = "query_name"
COL_FTP  = "ftp_refseq"

# 行为配置
SKIP_IF_EXISTS = True    # 输出的 {query_name}.fasta 已存在时是否跳过
SLEEP_BETWEEN  = 0.3     # 每条下载之间的暂停秒数（礼貌速率）
TIMEOUT        = 120     # 每个请求的超时秒
RETRY_TIMES    = 2       # 失败重试次数（针对网络波动；总尝试=RETRY_TIMES+1）

# HTTP 头（建议写上自己的邮箱以友好标识）
UA = "genome-downloader/1.2 (contact: your_email@example.com)"


def to_https(url: str) -> str:
    """把 ftp:// 或 http:// 统一转成 https://，保持路径不变"""
    if not isinstance(url, str):
        return url
    url = url.strip()
    if not url:
        return url
    parts = urlsplit(url)
    if parts.scheme in ("ftp", "http"):
        parts = parts._replace(scheme="https")
        return urlunsplit(parts)
    return url


def ensure_outdir(path: str) -> None:
    pathlib.Path(path).expanduser().mkdir(parents=True, exist_ok=True)


def read_mapping(excel_path: str, col_acc: str, col_name: str, col_ftp: str) -> pd.DataFrame:
    df = pd.read_excel(excel_path, dtype=str)
    for c in (col_acc, col_name, col_ftp):
        if c not in df.columns:
            raise KeyError(f"在 Excel 中找不到列: {c}；现有列: {list(df.columns)}")
    # 清洗
    df[col_acc]  = df[col_acc].astype(str).str.strip()
    df[col_name] = df[col_name].astype(str).str.strip()
    df[col_ftp]  = df[col_ftp].astype(str).str.strip()
    df = df.dropna(subset=[col_acc, col_name, col_ftp])
    df = df[(df[col_acc] != "") & (df[col_name] != "") & (df[col_ftp] != "") & (df[col_ftp].str.lower() != "na")]
    return df


def http_get(url: str) -> Optional[bytes]:
    """带超时与简单重试的 GET。200 返回内容；4xx 直接放弃；5xx/网络异常走重试。"""
    for attempt in range(1, RETRY_TIMES + 2):
        try:
            r = requests.get(url, headers={"User-Agent": UA}, timeout=TIMEOUT, stream=True)
            if r.status_code == 200:
                return r.content
            logging.warning(f"GET {url} → HTTP {r.status_code}")
            if 400 <= r.status_code < 500:
                return None  # 客户端错误不重试
        except requests.RequestException as e:
            logging.warning(f"GET 异常({attempt}): {e}")
        time.sleep(0.8 * attempt)
    return None


def ftp_from_excel(acc: str, ftp_url: str) -> Optional[bytes]:
    """
    从 Excel 提供的 ftp_refseq 地址下载 *_genomic.fna.gz 或 *_genomic.fna
    返回 FASTA bytes（若是 .gz 会解压）
    """
    if not isinstance(ftp_url, str) or not ftp_url or ftp_url.lower() == "na":
        return None

    ftp_url = to_https(ftp_url).rstrip("/")  # ★ 把 ftp:// 转成 https://
    base = ftp_url.split("/")[-1]            # 例: GCF_000001405.40_GRCh38.p14

    # 1) 尝试 gz
    gz_url = f"{ftp_url}/{base}_genomic.fna.gz"
    raw = http_get(gz_url)
    if raw:
        try:
            with io.BytesIO(raw) as bio:
                with gzip.GzipFile(fileobj=bio) as gzf:
                    data = gzf.read()
                    if data:
                        return data
        except OSError as e:
            logging.warning(f"解压失败(可能文件损坏) {acc}: {gz_url} -> {e}")

    # 2) 回退未压缩 .fna
    fna_url = f"{ftp_url}/{base}_genomic.fna"
    raw = http_get(fna_url)
    if raw:
        return raw

    # 3) 再回退 cds_from_genomic（某些条目只有 CDS）
    cds_gz = f"{ftp_url}/{base}_cds_from_genomic.fna.gz"
    raw = http_get(cds_gz)
    if raw:
        try:
            with io.BytesIO(raw) as bio:
                with gzip.GzipFile(fileobj=bio) as gzf:
                    data = gzf.read()
                    if data:
                        return data
        except OSError:
            pass

    cds_fna = f"{ftp_url}/{base}_cds_from_genomic.fna"
    raw = http_get(cds_fna)
    if raw:
        return raw

    return None


def save_fasta_bytes(fasta_bytes: bytes, out_path: str) -> None:
    pathlib.Path(os.path.dirname(out_path)).mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        f.write(fasta_bytes)


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
    )

    ensure_outdir(OUT_DIR)
    df = read_mapping(EXCEL_PATH, COL_ACC, COL_NAME, COL_FTP)
    logging.info(f"读取到 {len(df)} 条记录（含有效 ftp_refseq）")

    ok, fail = 0, 0
    for i, row in df.iterrows():
        acc = str(row[COL_ACC]).strip()
        name = str(row[COL_NAME]).strip()
        ftp_url = str(row[COL_FTP]).strip()

        if not acc or not name or not ftp_url:
            logging.warning(f"第{i}行无效（空 acc/name/ftp），跳过")
            continue

        out_path = os.path.join(OUT_DIR, f"{name}.fasta")
        if SKIP_IF_EXISTS and os.path.exists(out_path):
            logging.info(f"[{i+1}/{len(df)}] 已存在，跳过: {out_path}")
            continue

        logging.info(f"[{i+1}/{len(df)}] 下载 {acc} → {name}.fasta")
        fasta = ftp_from_excel(acc, ftp_url)

        if fasta:
            save_fasta_bytes(fasta, out_path)
            logging.info(f"✅ 保存成功: {out_path} ({len(fasta):,} bytes)")
            ok += 1
        else:
            logging.error(f"❌ 失败: {acc} → {name}.fasta  (ftp_refseq={ftp_url})")
            fail += 1

        time.sleep(SLEEP_BETWEEN)

    logging.info(f"完成。成功 {ok} 条，失败 {fail} 条。输出目录: {OUT_DIR}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n中断退出。")
        sys.exit(1)