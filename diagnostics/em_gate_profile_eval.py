#!/usr/bin/env python3
"""Replay EM gate profiles on existing baseline and EM outputs.

This is diagnostic-only: it does not modify per-gene calls.tsv files. It lets us
test per-gene gate thresholds against a truth set after EM has already produced
<gene>.calls.tsv and <gene>.summary.tsv files.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from em_refine_gate import mean_mask_fraction, quartet_residual, read_selected, read_summary, read_tf_counts
from evaluate_calls import (
    load_g_group,
    load_truth,
    norm_allele,
    normalize_for_display,
    overlap,
    overlap_g_group_truth_resolution,
)


DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]
GENE_CLASS = {
    "HLA-A": "I", "HLA-B": "I", "HLA-C": "I",
    "HLA-DRB1": "II", "HLA-DPB1": "II", "HLA-DQB1": "II",
}


def parse_gene_float_map(spec: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in spec.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        gene, value = item.split("=", 1)
        out[gene.strip()] = float(value)
    return out


def builtin_profiles() -> dict[str, dict[str, float]]:
    return {
        "global_max_070": {},
        "pipeline_default": {"HLA-DRB1": 0.35, "HLA-DPB1": 0.20, "HLA-DQB1": 0.50},
        "class2_max_050": {"HLA-DRB1": 0.50, "HLA-DPB1": 0.50, "HLA-DQB1": 0.50},
        "class2_max_035": {"HLA-DRB1": 0.35, "HLA-DPB1": 0.35, "HLA-DQB1": 0.50},
        "class2_max_025": {"HLA-DRB1": 0.25, "HLA-DPB1": 0.25, "HLA-DQB1": 0.50},
        "dpb1_strict_020": {"HLA-DPB1": 0.20},
        "drb1_strict_020": {"HLA-DRB1": 0.20},
        "dpb1_noem_drb1_035": {"HLA-DRB1": 0.35, "HLA-DPB1": 0.00, "HLA-DQB1": 0.50},
        "drb1_dpb1_noem": {"HLA-DRB1": 0.00, "HLA-DPB1": 0.00, "HLA-DQB1": 0.50},
    }


def read_summary_rows(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["sample"], row["set"], row["gene"])
            if key not in seen:
                rows.append(key)
                seen.add(key)
    return rows


def read_call(path: Path) -> tuple[list[str], list[str]]:
    recipient: list[str] = []
    donor: list[str] = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("assignment") == "R":
                recipient.append(row.get("allele", "NA"))
            elif row.get("assignment") == "D":
                donor.append(row.get("allele", "NA"))
    return (recipient + ["NA", "NA"])[:2], (donor + ["NA", "NA"])[:2]


def baseline_path(gene_dir: Path) -> Path:
    saved_baseline = gene_dir / "calls.baseline.tsv"
    return saved_baseline if saved_baseline.exists() else gene_dir / "calls.tsv"


def choose_call_path(args, sample: str, gene: str, gene_max_diff: dict[str, float]) -> tuple[str, Path]:
    gene_dir = args.asm_root / sample / gene.lower() / gene
    baseline_calls = baseline_path(gene_dir)
    em_dir = args.em_root / sample / "em_refine"
    em_calls = em_dir / f"{gene}.calls.tsv"
    summary_path = em_dir / f"{gene}.summary.tsv"
    tf_counts_path = em_dir / f"{gene}.tf_counts.tsv"
    if not (baseline_calls.exists() and em_calls.exists() and summary_path.exists() and tf_counts_path.exists()):
        return "KEEP_BASELINE", baseline_calls

    summary = read_summary(summary_path)
    diff = float(summary.get("sum_abs_diff", "inf"))
    top_frac = float(summary.get("top_frac", "0"))
    chi_r = float(summary.get("chi_r_fit", "0"))
    max_diff = gene_max_diff.get(gene, args.max_diff)
    mask = mean_mask_fraction(gene_dir)
    selected = read_selected(em_calls)
    baseline_selected = read_selected(baseline_calls)
    tf_counts = read_tf_counts(tf_counts_path)
    min_selected = min((tf_counts.get(allele, 0.0) for allele in selected), default=0.0)
    baseline_supported = sum(
        1 for allele in baseline_selected
        if tf_counts.get(allele, 0.0) >= args.baseline_min_frac
    )
    baseline_residual = quartet_residual(baseline_selected, tf_counts, chi_r)

    keep = False
    if diff >= max_diff:
        if not (baseline_supported == 0 and selected and min_selected >= args.selected_min_frac):
            keep = True
    if (mask >= args.high_mask and diff >= args.ambiguous_diff
            and top_frac <= args.ambiguous_top_frac and baseline_supported >= 3):
        keep = True
    if (baseline_residual is not None and mask >= args.high_mask
            and baseline_supported >= 4 and baseline_residual - diff <= args.baseline_near_tie):
        keep = True
    if (gene in args.class_i_baseline_genes and diff >= args.class_i_baseline_diff
            and top_frac <= args.class_i_baseline_top_frac and baseline_supported >= 4):
        keep = True
    if selected and min_selected < args.selected_min_frac and diff >= args.ambiguous_diff:
        keep = True
    return ("KEEP_BASELINE", baseline_calls) if keep else ("OVERRIDE", em_calls)


def score_gene(truth, gene: str, recipient: list[str], donor: list[str], gmap) -> tuple[int, int, int]:
    ok_2field = 0
    ok_g_group = 0
    total = 0
    for side, pred in (("PATIENT", recipient), ("DONOR", donor)):
        truth_vals = truth[side][gene]
        truth_2field = normalize_for_display(truth_vals, "2field", gmap)
        pred_2field = sorted(norm_allele(value, "2field", gmap) for value in pred)
        ok_2field += overlap(truth_2field, pred_2field)
        ok_g_group += overlap_g_group_truth_resolution(truth_vals, pred, gmap)
        total += len(truth_vals)
    return ok_2field, ok_g_group, total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--em-root", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--truth-dir", required=True, type=Path)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--profile", action="append", default=[],
                        help="extra profile as name:HLA-DRB1=0.35,HLA-DPB1=0.20")
    parser.add_argument("--max-diff", type=float, default=0.70)
    parser.add_argument("--high-mask", type=float, default=0.40)
    parser.add_argument("--ambiguous-diff", type=float, default=0.35)
    parser.add_argument("--ambiguous-top-frac", type=float, default=0.42)
    parser.add_argument("--selected-min-frac", type=float, default=0.0005)
    parser.add_argument("--baseline-min-frac", type=float, default=0.005)
    parser.add_argument("--baseline-near-tie", type=float, default=0.001)
    parser.add_argument("--class-i-baseline-genes", default="HLA-A")
    parser.add_argument("--class-i-baseline-diff", type=float, default=0.20)
    parser.add_argument("--class-i-baseline-top-frac", type=float, default=0.50)
    args = parser.parse_args()

    args.class_i_baseline_genes = {
        gene.strip() for gene in args.class_i_baseline_genes.split(",") if gene.strip()
    }
    profiles = builtin_profiles()
    for profile in args.profile:
        name, spec = profile.split(":", 1)
        profiles[name] = parse_gene_float_map(spec)

    rows = read_summary_rows(args.summary)
    truth_by_set = {
        set_id: load_truth(args.truth_dir / f"truth_typing-set-{set_id}.tsv")
        for set_id in sorted({set_id for _, set_id, _ in rows})
    }
    gmap = load_g_group(args.g_group)

    print("profile\tscore2\tscoreg\tI\tII\toverrides\tHLA-A\tHLA-B\tHLA-C\tHLA-DRB1\tHLA-DPB1\tHLA-DQB1")
    for profile_name, gene_max_diff in profiles.items():
        total_2field = 0
        total_g_group = 0
        total = 0
        by_class: defaultdict[str, int] = defaultdict(int)
        by_gene: defaultdict[str, int] = defaultdict(int)
        actions: Counter[tuple[str, str]] = Counter()
        for sample, set_id, gene in rows:
            action, call_path = choose_call_path(args, sample, gene, gene_max_diff)
            actions[(gene, action)] += 1
            recipient, donor = read_call(call_path)
            ok_2field, ok_g_group, denom = score_gene(truth_by_set[set_id], gene, recipient, donor, gmap)
            total_2field += ok_2field
            total_g_group += ok_g_group
            total += denom
            by_class[GENE_CLASS[gene]] += ok_2field
            by_gene[gene] += ok_2field
        overrides = sum(count for (gene, action), count in actions.items() if action == "OVERRIDE")
        fields = [
            profile_name,
            f"{total_2field}/{total}",
            f"{total_g_group}/{total}",
            str(by_class["I"]),
            str(by_class["II"]),
            str(overrides),
        ] + [str(by_gene[gene]) for gene in DEFAULT_GENES]
        print("\t".join(fields))


if __name__ == "__main__":
    main()
