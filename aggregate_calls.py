#!/usr/bin/env python3
"""Aggregate per-gene calls.tsv files into a single summary table.

Reads <ASM_ROOT>/<SAMPLE>/<gene_lc>/<TAG>/calls.tsv for each gene and writes
a single tab-separated file with columns:

    sample  gene  R1  R2  D1  D2  source

`source` is `em-refined` if a `calls.baseline.tsv` sibling exists (meaning the
EM stage overrode the baseline), otherwise `baseline`.

Usage:
    aggregate_calls.py --asm-root asm_v2 --sample mySample \\
        [--genes HLA-A HLA-B ...] [--out final_calls.tsv]

Defaults to the 6 typed genes and writes
<asm-root>/<sample>/<sample>.final_calls.tsv.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]


def read_calls(path: Path):
    """Return list of (global_hap, assignment, allele) sorted by global_hap."""
    rows = []
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_h = header.index("global_hap")
            i_a = header.index("assignment")
            i_l = header.index("allele")
        except ValueError:
            sys.stderr.write(f"[warn] {path}: unexpected header {header}\n")
            return []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_h, i_a, i_l):
                continue
            rows.append((parts[i_h], parts[i_a], parts[i_l]))
    rows.sort(key=lambda r: int(r[0]) if r[0].isdigit() else r[0])
    return rows


def collect(asm_root: Path, sample: str, genes):
    out_rows = []
    for gene in genes:
        gene_lc = gene.lower()
        d = asm_root / sample / gene_lc / gene
        calls = d / "calls.tsv"
        if not calls.exists():
            out_rows.append((sample, gene, "NA", "NA", "NA", "NA", "missing"))
            continue
        rows = read_calls(calls)
        # bucket by R / D
        rs = [a for (_, t, a) in rows if t == "R"]
        ds = [a for (_, t, a) in rows if t == "D"]
        # pad / truncate to two each
        rs = (rs + ["NA", "NA"])[:2]
        ds = (ds + ["NA", "NA"])[:2]
        source = "em-refined" if (d / "calls.baseline.tsv").exists() else "baseline"
        out_rows.append((sample, gene, rs[0], rs[1], ds[0], ds[1], source))
    return out_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    rows = collect(args.asm_root, args.sample, args.genes)
    out_path = args.out or (args.asm_root / args.sample / f"{args.sample}.final_calls.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("sample\tgene\tR1\tR2\tD1\tD2\tsource\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    sys.stderr.write(f"[aggregate] wrote {out_path} ({len(rows)} genes)\n")


if __name__ == "__main__":
    main()
