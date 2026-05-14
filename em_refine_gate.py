#!/usr/bin/env python3
"""Decide whether EM-refined calls should override baseline calls.

The gate is intentionally truth-free. It uses only quality signals produced by
the pipeline: EM fit residual, EM top fraction, selected allele support, and
masked assembly fraction. The goal is to keep the strong EM gains for ordinary
genes while preventing high-mask, high-ambiguity loci from replacing a stable
baseline quartet with a poorly supported EM solution.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Optional


def parse_gene_float_map(spec: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in spec.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise SystemExit(f"invalid per-gene parameter entry: {item!r}")
        gene, value = item.split("=", 1)
        try:
            out[gene.strip()] = float(value)
        except ValueError:
            raise SystemExit(f"invalid float in per-gene parameter entry: {item!r}")
    return out


def gene_float(gene: str, default: float, overrides: dict[str, float]) -> float:
    return overrides.get(gene, default)


def fasta_n_fraction(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    seq = []
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    s = "".join(seq)
    return None if not s else s.count("N") / len(s)


def mean_mask_fraction(gene_dir: Path) -> float:
    vals = [fasta_n_fraction(gene_dir / f"hap{i}.fa") for i in range(1, 5)]
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else 0.0


def read_summary(path: Path) -> dict[str, str]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else {}


def two_field(allele: str) -> str:
    allele = allele.strip().replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    rest = rest.rstrip("GP")
    fields = rest.split(":")
    return f"{gene}*{fields[0]}:{fields[1]}" if len(fields) >= 2 else f"{gene}*{fields[0]}"


def read_selected(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for side in ("R", "D"):
        out.extend(two_field(r["allele"]) for r in rows if r.get("assignment") == side)
    return out[:4]


def read_tf_counts(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[row["allele_2field"]] = float(row["fraction"])
            except (KeyError, ValueError):
                continue
    return out


def quartet_residual(quartet: list[str], fractions: dict[str, float], chi_r: float) -> Optional[float]:
    if len(quartet) < 4:
        return None
    expected = {}
    for allele in quartet[:2]:
        expected[allele] = expected.get(allele, 0.0) + chi_r / 2.0
    for allele in quartet[2:4]:
        expected[allele] = expected.get(allele, 0.0) + (1.0 - chi_r) / 2.0
    keys = set(fractions) | set(expected)
    return sum(abs(fractions.get(allele, 0.0) - expected.get(allele, 0.0)) for allele in keys)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--gene-dir", required=True, type=Path)
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--em-calls", required=True, type=Path)
    ap.add_argument("--baseline-calls", required=True, type=Path)
    ap.add_argument("--tf-counts", required=True, type=Path)
    ap.add_argument("--max-diff", type=float, required=True)
    ap.add_argument("--high-mask", type=float, default=0.40)
    ap.add_argument("--ambiguous-diff", type=float, default=0.35)
    ap.add_argument("--ambiguous-top-frac", type=float, default=0.42)
    ap.add_argument("--selected-min-frac", type=float, default=0.0005)
    ap.add_argument("--baseline-min-frac", type=float, default=0.005)
    ap.add_argument("--baseline-near-tie", type=float, default=0.001,
                    help="keep a fully supported high-mask baseline quartet when its residual is within this margin")
    ap.add_argument("--class-i-baseline-genes", default="HLA-A",
                    help="comma-separated class-I genes where a high EM residual can keep a fully supported baseline")
    ap.add_argument("--class-i-baseline-diff", type=float, default=0.20,
                    help="minimum EM sum_abs_diff for class-I baseline protection")
    ap.add_argument("--class-i-baseline-top-frac", type=float, default=0.50,
                    help="maximum EM top_frac for class-I baseline protection")
    ap.add_argument("--gene-max-diff", default="",
                    help="comma-separated per-gene overrides, e.g. HLA-DRB1=0.35,HLA-DPB1=0.30")
    ap.add_argument("--gene-high-mask", default="",
                    help="comma-separated per-gene overrides for --high-mask")
    ap.add_argument("--gene-ambiguous-diff", default="",
                    help="comma-separated per-gene overrides for --ambiguous-diff")
    ap.add_argument("--gene-ambiguous-top-frac", default="",
                    help="comma-separated per-gene overrides for --ambiguous-top-frac")
    ap.add_argument("--gene-selected-min-frac", default="",
                    help="comma-separated per-gene overrides for --selected-min-frac")
    ap.add_argument("--gene-baseline-min-frac", default="",
                    help="comma-separated per-gene overrides for --baseline-min-frac")
    ap.add_argument("--gene-baseline-near-tie", default="",
                    help="comma-separated per-gene overrides for --baseline-near-tie")
    args = ap.parse_args()

    max_diff = gene_float(args.gene, args.max_diff, parse_gene_float_map(args.gene_max_diff))
    high_mask = gene_float(args.gene, args.high_mask, parse_gene_float_map(args.gene_high_mask))
    ambiguous_diff = gene_float(args.gene, args.ambiguous_diff, parse_gene_float_map(args.gene_ambiguous_diff))
    ambiguous_top_frac = gene_float(args.gene, args.ambiguous_top_frac, parse_gene_float_map(args.gene_ambiguous_top_frac))
    selected_min_frac = gene_float(args.gene, args.selected_min_frac, parse_gene_float_map(args.gene_selected_min_frac))
    baseline_min_frac = gene_float(args.gene, args.baseline_min_frac, parse_gene_float_map(args.gene_baseline_min_frac))
    baseline_near_tie = gene_float(args.gene, args.baseline_near_tie, parse_gene_float_map(args.gene_baseline_near_tie))

    summary = read_summary(args.summary)
    diff = float(summary.get("sum_abs_diff", "inf"))
    top_frac = float(summary.get("top_frac", "0"))
    chi_r = float(summary.get("chi_r_fit", "0"))
    mask = mean_mask_fraction(args.gene_dir)
    selected = read_selected(args.em_calls)
    baseline_selected = read_selected(args.baseline_calls)
    tf_counts = read_tf_counts(args.tf_counts)
    min_selected = min((tf_counts.get(a, 0.0) for a in selected), default=0.0)
    baseline_min = min((tf_counts.get(a, 0.0) for a in baseline_selected), default=0.0)
    baseline_supported = sum(
        1 for allele in baseline_selected
        if tf_counts.get(allele, 0.0) >= baseline_min_frac
    )
    baseline_residual = quartet_residual(baseline_selected, tf_counts, chi_r)
    class_i_baseline_genes = {
        g.strip() for g in args.class_i_baseline_genes.split(",") if g.strip()
    }

    keep = False
    reason = []
    override_reason = []
    if diff >= max_diff:
        if baseline_supported == 0 and selected and min_selected >= selected_min_frac:
            override_reason.append(
                f"sumAbsDiff={diff:.4f}>={max_diff:.4f} but baseline_unsupported"
            )
        else:
            keep = True
            reason.append(f"sumAbsDiff={diff:.4f}>={max_diff:.4f}")
    if (mask >= high_mask and diff >= ambiguous_diff
            and top_frac <= ambiguous_top_frac
            and baseline_supported >= 3):
        keep = True
        reason.append(
            f"ambiguous_high_mask mask={mask:.4f} diff={diff:.4f} "
            f"top_frac={top_frac:.4f} baseline_supported={baseline_supported}/4"
        )
    if (baseline_residual is not None and mask >= high_mask
            and baseline_supported >= 4
            and baseline_residual - diff <= baseline_near_tie):
        keep = True
        reason.append(
            f"high_mask_baseline_near_tie mask={mask:.4f} "
            f"baseline_residual={baseline_residual:.4f} diff={diff:.4f} "
            f"delta={baseline_residual - diff:.4f}"
        )
    if (args.gene in class_i_baseline_genes
            and diff >= args.class_i_baseline_diff
            and top_frac <= args.class_i_baseline_top_frac
            and baseline_supported >= 4):
        keep = True
        reason.append(
            f"class_i_high_em_residual gene={args.gene} diff={diff:.4f} "
            f"top_frac={top_frac:.4f} baseline_supported={baseline_supported}/4"
        )
    if selected and min_selected < selected_min_frac and diff >= ambiguous_diff:
        keep = True
        reason.append(
            f"low_selected_support min_selected={min_selected:.6f} diff={diff:.4f}"
        )

    if keep:
        print("KEEP_BASELINE\t" + ";".join(reason))
    else:
        print(
            "OVERRIDE\t"
            f"{';'.join(override_reason)};"
            f"sumAbsDiff={diff:.4f};top_frac={top_frac:.4f};mask={mask:.4f};"
            f"max_diff={max_diff:.4f};high_mask={high_mask:.4f};"
            f"min_selected={min_selected:.6f};baseline_min={baseline_min:.6f};"
            f"baseline_supported={baseline_supported}/4;"
            f"baseline_residual={baseline_residual if baseline_residual is not None else 'NA'}"
        )


if __name__ == "__main__":
    main()
