#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal GenomeSPOT runner (per-sample TSV only)
- Input : a directory (or a single .faa file)
- Effect: For each .faa, run GenomeSPOT and write <workdir>/<stem>.predictions.tsv

Usage examples:
  # contigs and .faa in the same folder (same stem)
  python run_genomespot_batch.py \
      --input "/Users/frank_gong/文档/生物智能体/硬盘备份/faa文件" \
      --workdir outputs_genomespot

  # contigs kept separately (same stem, different dir)
  python run_genomespot_batch.py \
      --input "/Users/frank_gong/文档/生物智能体/硬盘备份/faa文件" \
      --contigs-dir "/Users/frank_gong/文档/生物智能体/硬盘备份/fa文件" \
      --workdir outputs_genomespot
"""

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import warnings
from pandas.errors import PerformanceWarning
import re


CONTIG_SUFFIXES = (".fna", ".fa", ".fasta", ".ffn", ".fsa")


def find_models_dir() -> str:
    """Locate GenomeSPOT models/ directory from installed package or local repo."""
    try:
        import genome_spot  # type: ignore
        pkg_dir = Path(genome_spot.__file__).parent
        models = pkg_dir / "models"
        if models.exists():
            return str(models)
    except Exception:
        pass
    local = Path(__file__).resolve().parent / "models"
    if local.exists():
        return str(local)
    raise RuntimeError(
        "Could not locate GenomeSPOT models/ directory. Install genome_spot and ensure models/ exist, or pass --models explicitly."
    )


def collect_faa_paths(inp: Path):
    if inp.is_dir():
        return sorted([p for p in inp.glob("*.faa") if p.is_file()])
    elif inp.is_file() and inp.suffix.lower() == ".faa":
        return [inp]
    else:
        raise ValueError(f"Input must be a .faa file or a directory containing .faa files: {inp}")


def find_contigs_for(faa_path: Path, contigs_dir: Path | None = None) -> Path | None:
    stem = faa_path.stem
    # 1) alongside the .faa
    for suf in CONTIG_SUFFIXES:
        p = faa_path.with_suffix(suf)
        if p.exists():
            return p
    # 2) in the provided contigs_dir
    if contigs_dir:
        for suf in CONTIG_SUFFIXES:
            p = contigs_dir / f"{stem}{suf}"
            if p.exists():
                return p
    return None


def run_genomespot(models_dir: str, contigs: Path, faa: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = faa.stem
    out_prefix = out_dir / stem
    cmd = [
        sys.executable, "-m", "genome_spot.genome_spot",
        "-m", models_dir,
        "-c", str(contigs),
        "-p", str(faa),
        "-o", str(out_prefix),
    ]
    print(f"[RUN] {stem}")
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"[ERROR] GenomeSPOT failed for {stem}:\n{res.stderr}", file=sys.stderr)
    else:
        print(f"[OK] {stem} → {out_prefix.with_suffix('.predictions.tsv')}")
    return out_prefix.with_suffix(".predictions.tsv")


def main():
    ap = argparse.ArgumentParser(description="Minimal per-sample GenomeSPOT runner")
    ap.add_argument("--input", required=True, help="Path to a .faa file or a directory containing .faa files")
    ap.add_argument("--workdir", default="outputs_genomespot", help="Directory to store per-sample outputs")
    ap.add_argument("--models", default=None, help="Path to GenomeSPOT models/ (optional; auto-detect if omitted)")
    ap.add_argument("--contigs-dir", default=None, help="Optional directory with contigs (.fna/.fa/.fasta) matching .faa stems")
    ap.add_argument("--summary", default=None, help="Optional path to write aggregated CSV (defaults to <workdir>/summary_predictions.csv)")
    args = ap.parse_args()

    inp = Path(args.input).expanduser().resolve()
    workdir = Path(args.workdir).expanduser().resolve()
    contigs_dir = Path(args.contigs_dir).expanduser().resolve() if args.contigs_dir else None

    models_dir = args.models or find_models_dir()
    faa_list = collect_faa_paths(inp)
    if not faa_list:
        print(f"No .faa files found in {inp}", file=sys.stderr)
        sys.exit(2)

    for faa in faa_list:
        contigs = find_contigs_for(faa, contigs_dir)
        if contigs is None:
            print(f"[WARN] No contigs (.fna/.fa/.fasta) found for {faa.name}. Skipping.", file=sys.stderr)
            continue
        run_genomespot(models_dir, contigs, faa, workdir)

    # === Aggregate all per-sample TSVs into one CSV ===
    tsv_paths = sorted(workdir.glob("*.predictions.tsv"))
    rows = []
    for tsv in tsv_paths:
        try:
            df = pd.read_csv(tsv, sep="\t", header=0)
            # Expect rows like: target, value, error, units, is_novel, warning
            if "target" not in df.columns or "value" not in df.columns:
                print(f"[WARN] {tsv.name} has unexpected columns: {list(df.columns)}", file=sys.stderr)
                continue
            # Normalize target names for robust lookup
            df["_target_norm"] = (
                df["target"].astype(str)
                .str.strip()
                .str.lower()
                .str.replace(" ", "_", regex=False)
                .str.replace("-", "_", regex=False)
            )
            lookup = dict(zip(df["_target_norm"], df["value"]))

            def get_by_targets(cands):
                for k in cands:
                    if k in lookup:
                        return lookup[k]
                return np.nan

            row = {
                "strain": tsv.stem.replace(".predictions", ""),
                "temperature_optimum_C": get_by_targets(["temperature_optimum", "temp_optimum", "topt"]),
                "temperature_minimum":   get_by_targets(["temperature_minimum", "temperature_min", "temp_min", "tmin"]),
                "temperature_maximum":   get_by_targets(["temperature_maximum", "temperature_max", "temp_max", "tmax"]),
                "ph_optimum":            get_by_targets(["ph_optimum", "pH_optimum".lower()]),
                "ph_minimum":            get_by_targets(["ph_minimum", "ph_min"]),
                "ph_maximum":            get_by_targets(["ph_maximum", "ph_max"]),
                "salinity_optimum":      get_by_targets(["salinity_optimum", "salt_optimum", "nacl_optimum"]),
                "salinity_minimum":      get_by_targets(["salinity_minimum", "salinity_min", "salt_min", "nacl_min"]),
                "salinity_maximum":      get_by_targets(["salinity_maximum", "salinity_max", "salt_max", "nacl_max"]),
                "oxygen_tolerance":      get_by_targets(["oxygen", "o2", "oxygen_tolerance"]),
            }

            # If all metrics are NaN, dump debug info once
            if all(pd.isna(v) for k, v in row.items() if k != "strain"):
                print(f"[DEBUG] Could not map targets in {tsv.name}. Available targets:", file=sys.stderr)
                print(df[["target", "value"]].to_string(index=False), file=sys.stderr)

            rows.append(row)
        except Exception as e:
            print(f"[WARN] Failed to read {tsv.name}: {e}", file=sys.stderr)
            continue

    if rows:
        columns_order = [
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
        agg = pd.DataFrame(rows, columns=columns_order)
        summary_path = Path(args.summary).expanduser().resolve() if args.summary else (workdir / "summary_predictions.csv")
        agg.to_csv(summary_path, index=False)
        print(f"[OK] Wrote aggregated CSV: {summary_path} (n={len(agg)})")
    else:
        print("[WARN] No per-sample TSVs found to aggregate.", file=sys.stderr)

    print("[DONE] Finished running GenomeSPOT for all available samples.")


if __name__ == "__main__":
    main()