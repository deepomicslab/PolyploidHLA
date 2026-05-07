#!/usr/bin/env python3
"""Aggregate per-gene calls.tsv files into a single summary table.

Reads <ASM_ROOT>/<SAMPLE>/<gene_lc>/<TAG>/calls.tsv for each gene and writes
a single tab-separated file. The full allele calls are preserved, and extra
2-field / report columns are emitted for low-resolution truth comparison and
low-coverage genes.

    sample  gene  R1_full  R2_full  D1_full  D2_full  R1_2field ...

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
import sys
from pathlib import Path
from typing import Optional

DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]


def allele_2field(allele: str) -> str:
    if not allele or allele == "NA" or "*" not in allele:
        return allele or "NA"
    gene, fields = allele.replace("HLA-", "").split("*", 1)
    parts = fields.replace("G", "").replace("N", "").split(":")
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def fasta_n_fraction(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    seq = []
    with path.open() as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    s = "".join(seq)
    return None if not s else s.count("N") / len(s)


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


def collect(asm_root: Path, sample: str, genes, mask_warn: float):
    out_rows = []
    for gene in genes:
        gene_lc = gene.lower()
        d = asm_root / sample / gene_lc / gene
        calls = d / "calls.tsv"
        mask_values = [fasta_n_fraction(d / f"hap{i}.fa") for i in range(1, 5)]
        mask_values = [x for x in mask_values if x is not None]
        mean_mask = sum(mask_values) / len(mask_values) if mask_values else None
        if not calls.exists():
            out_rows.append({
                "sample": sample, "gene": gene,
                "R1_full": "NA", "R2_full": "NA", "D1_full": "NA", "D2_full": "NA",
                "source": "missing", "mean_mask_fraction": "NA",
                "report_level": "missing", "warning": "missing_calls_tsv",
            })
            continue
        rows = read_calls(calls)
        rs = [a for (_, t, a) in rows if t == "R"]
        ds = [a for (_, t, a) in rows if t == "D"]
        rs = (rs + ["NA", "NA"])[:2]
        ds = (ds + ["NA", "NA"])[:2]
        source = "em-refined" if (d / "calls.baseline.tsv").exists() else "baseline"
        high_mask = mean_mask is not None and mean_mask >= mask_warn
        report_level = "2-field" if high_mask else "full"
        warning = "high_mask_report_2field" if high_mask else ""
        out_rows.append({
            "sample": sample, "gene": gene,
            "R1_full": rs[0], "R2_full": rs[1], "D1_full": ds[0], "D2_full": ds[1],
            "R1_2field": allele_2field(rs[0]), "R2_2field": allele_2field(rs[1]),
            "D1_2field": allele_2field(ds[0]), "D2_2field": allele_2field(ds[1]),
            "R1_report": allele_2field(rs[0]) if high_mask else rs[0],
            "R2_report": allele_2field(rs[1]) if high_mask else rs[1],
            "D1_report": allele_2field(ds[0]) if high_mask else ds[0],
            "D2_report": allele_2field(ds[1]) if high_mask else ds[1],
            "source": source,
            "mean_mask_fraction": "NA" if mean_mask is None else f"{mean_mask:.4f}",
            "report_level": report_level,
            "warning": warning,
        })
    return out_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--mask-warn", type=float, default=0.15,
                    help="mean hap FASTA N fraction above which report columns "
                    "are downgraded to 2-field")
    args = ap.parse_args()

    rows = collect(args.asm_root, args.sample, args.genes, args.mask_warn)
    out_path = args.out or (args.asm_root / args.sample / f"{args.sample}.final_calls.tsv")
    cols = [
        "sample", "gene",
        "R1_full", "R2_full", "D1_full", "D2_full",
        "R1_2field", "R2_2field", "D1_2field", "D2_2field",
        "R1_report", "R2_report", "D1_report", "D2_report",
        "source", "mean_mask_fraction", "report_level", "warning",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    sys.stderr.write(f"[aggregate] wrote {out_path} ({len(rows)} genes)\n")


if __name__ == "__main__":
    main()
