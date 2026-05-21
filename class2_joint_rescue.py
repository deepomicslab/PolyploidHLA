#!/usr/bin/env python3
"""Prototype class-II joint rescue strategies on an existing quartet summary.

This is a diagnostic evaluator. It uses truth only for scoring after each
truth-free rescue rule has produced a candidate quartet.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_SUMMARY = Path("diagnostics/quartet_summary_20260512.tsv")
DEFAULT_DIRECT = Path("diagnostics/direct_quartet_likelihood_classII_constrained_20260513.tsv")
DEFAULT_SPECHLA_ROOT = Path("/data2/wangxuedong/polyploid-hla-realsets/spechla_out_abc_realsets_rescue_20260512")
DEFAULT_OUT_TSV = Path("diagnostics/class2_joint_rescue_20260513.tsv")
DEFAULT_SUMMARY_OUT = Path("diagnostics/class2_joint_rescue_20260513.summary")

DRB1_DQB1_LD = {
    "DQB1*02:01": "DRB1*03:01",
    "DQB1*02:02": "DRB1*07:01",
    "DQB1*02:82": "DRB1*07:01",
    "DQB1*02:109": "DRB1*03:01",
    "DQB1*03:01": "DRB1*04:01",
    "DQB1*06:02": "DRB1*15:01",
}

DIRECT_CLASS2_GENES = {"HLA-DRB1", "HLA-DQB1", "HLA-DPB1"}


def split_alleles(text: str):
    return [allele for allele in (text or "").split(",") if allele]


def allele_number(allele: str) -> int:
    match = re.search(r"\*(\d+):", allele or "")
    return int(match.group(1)) if match else 999999


def side_overlap(predicted, truth) -> int:
    truth_counts = Counter(truth)
    hits = 0
    for allele in predicted:
        if truth_counts[allele] > 0:
            hits += 1
            truth_counts[allele] -= 1
    return hits


def quartet_score(quartet, truth_r, truth_d) -> int:
    return side_overlap(quartet[:2], truth_r) + side_overlap(quartet[2:], truth_d)


def quartet_from_row(row):
    return split_alleles(row["pred_R"]) + split_alleles(row["pred_D"])


def truth_from_row(row):
    return split_alleles(row["truth_R"]), split_alleles(row["truth_D"])


def read_tsv(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_summary(path: Path):
    rows = read_tsv(path)
    by_sample_gene = {(row["sample"], row["gene"]): row for row in rows}
    by_sample = defaultdict(dict)
    for row in rows:
        by_sample[row["sample"]][row["gene"]] = row
    return rows, by_sample_gene, by_sample


def direct_quartet(row):
    return split_alleles(row.get("direct_R", "")) + split_alleles(row.get("direct_D", ""))


def current_direct_quartet(row):
    return split_alleles(row.get("current_quartet", ""))


def load_direct_accepts(path: Path, gap_threshold: float):
    accepted = {}
    for row in read_tsv(path):
        try:
            gap = float(row["gap"])
        except (KeyError, ValueError):
            continue
        candidate = direct_quartet(row)
        if gap < gap_threshold or not candidate or candidate == current_direct_quartet(row):
            continue
        accepted[(row["sample"], row["gene"])] = candidate
    return accepted


def ld_dr_from_dq(dqb1_row):
    if not dqb1_row:
        return None
    dqb1_quartet = quartet_from_row(dqb1_row)
    drb1_quartet = []
    for allele in dqb1_quartet:
        drb1_allele = DRB1_DQB1_LD.get(allele)
        if not drb1_allele:
            return None
        drb1_quartet.append(drb1_allele)
    return drb1_quartet


def read_tf_counts(spechla_root: Path, sample: str, gene: str):
    path = spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
    rows = []
    if not path.exists():
        return rows
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({
                "allele": row["allele_2field"],
                "weight": float(row.get("em_weight") or 0.0),
                "fraction": float(row.get("fraction") or row.get("em_frac") or 0.0),
            })
    return sorted(rows, key=lambda row: -row["weight"])


def dpb1_rare_collapse(row, spechla_root: Path, min_fraction: float, rare_cutoff: int, top_common: int):
    quartet = quartet_from_row(row)
    if row["gene"] != "HLA-DPB1" or not any(allele_number(allele) >= rare_cutoff for allele in quartet):
        return quartet, False
    common = [
        count_row["allele"]
        for count_row in read_tf_counts(spechla_root, row["sample"], row["gene"])
        if allele_number(count_row["allele"]) < rare_cutoff and count_row["fraction"] >= min_fraction
    ]
    if not common:
        return quartet, False
    used = Counter(allele for allele in quartet if allele_number(allele) < rare_cutoff)
    candidate = list(quartet)
    changed = False
    for index, allele in enumerate(candidate):
        if allele_number(allele) < rare_cutoff:
            continue
        replacement = None
        for common_allele in common[:top_common]:
            if used[common_allele] < 2:
                replacement = common_allele
                break
        if replacement:
            candidate[index] = replacement
            used[replacement] += 1
            changed = True
    return candidate, changed


def apply_strategy(row, strategy: str, by_sample, direct_accepts, spechla_root: Path, args):
    current = quartet_from_row(row)
    reason = "current"
    candidate = current
    if strategy in {"direct_gate", "combined"} and row["gene"] in DIRECT_CLASS2_GENES:
        direct_candidate = direct_accepts.get((row["sample"], row["gene"]))
        if direct_candidate:
            candidate = direct_candidate
            reason = "direct_gate"
    if strategy in {"drdq_ld", "combined"} and row["gene"] == "HLA-DRB1":
        ld_candidate = ld_dr_from_dq(by_sample[row["sample"]].get("HLA-DQB1"))
        if ld_candidate:
            candidate = ld_candidate
            reason = "drdq_ld"
    if strategy in {"dpb1_rare", "combined"} and row["gene"] == "HLA-DPB1":
        rare_candidate, changed = dpb1_rare_collapse(
            row,
            spechla_root,
            args.dpb1_min_fraction,
            args.dpb1_rare_cutoff,
            args.dpb1_top_common,
        )
        if changed:
            candidate = rare_candidate
            reason = "dpb1_rare_collapse"
    return candidate, reason


def summarize(strategy_rows):
    total = sum(int(row["new_score"]) for row in strategy_rows)
    current_total = sum(int(row["current_score"]) for row in strategy_rows)
    changed = [row for row in strategy_rows if row["changed"] == "1"]
    improved = [row for row in changed if int(row["delta"]) > 0]
    regressed = [row for row in changed if int(row["delta"]) < 0]
    by_gene = defaultdict(lambda: [0, 0])
    for row in strategy_rows:
        by_gene[row["gene"]][0] += int(row["current_score"])
        by_gene[row["gene"]][1] += int(row["new_score"])
    return {
        "current_total": current_total,
        "new_total": total,
        "delta": total - current_total,
        "changed": len(changed),
        "improved": len(improved),
        "regressed": len(regressed),
        "by_gene": by_gene,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--direct-tsv", type=Path, default=DEFAULT_DIRECT)
    parser.add_argument("--spechla-root", type=Path, default=DEFAULT_SPECHLA_ROOT)
    parser.add_argument("--out-tsv", type=Path, default=DEFAULT_OUT_TSV)
    parser.add_argument("--summary-out", type=Path, default=DEFAULT_SUMMARY_OUT)
    parser.add_argument("--direct-gap", type=float, default=150.0)
    parser.add_argument("--dpb1-min-fraction", type=float, default=0.02)
    parser.add_argument("--dpb1-rare-cutoff", type=int, default=100)
    parser.add_argument("--dpb1-top-common", type=int, default=6)
    args = parser.parse_args()

    rows, _by_sample_gene, by_sample = load_summary(args.summary)
    direct_accepts = load_direct_accepts(args.direct_tsv, args.direct_gap)
    strategies = ["current", "direct_gate", "drdq_ld", "dpb1_rare", "combined"]
    out_rows = []

    for strategy in strategies:
        for row in rows:
            truth_r, truth_d = truth_from_row(row)
            current = quartet_from_row(row)
            candidate, reason = (current, "current") if strategy == "current" else apply_strategy(
                row,
                strategy,
                by_sample,
                direct_accepts,
                args.spechla_root,
                args,
            )
            current_score = int(row["score2"])
            new_score = quartet_score(candidate, truth_r, truth_d)
            out_rows.append({
                "strategy": strategy,
                "sample": row["sample"],
                "set": row["set"],
                "gene": row["gene"],
                "current_score": current_score,
                "new_score": new_score,
                "delta": new_score - current_score,
                "changed": "1" if candidate != current else "0",
                "reason": reason,
                "current_quartet": ",".join(current),
                "new_quartet": ",".join(candidate),
                "truth_quartet": f"{row['truth_R']}|{row['truth_D']}",
            })

    fields = [
        "strategy", "sample", "set", "gene", "current_score", "new_score", "delta",
        "changed", "reason", "current_quartet", "new_quartet", "truth_quartet",
    ]
    write_tsv(args.out_tsv, fields, out_rows)

    grouped = defaultdict(list)
    for row in out_rows:
        grouped[row["strategy"]].append(row)
    lines = ["strategy\tcurrent\tnew\tdelta\tchanged\timproved\tregressed\tby_gene\n"]
    for strategy in strategies:
        stats = summarize(grouped[strategy])
        by_gene_text = ";".join(
            f"{gene}:{current_score}->{new_score}"
            for gene, (current_score, new_score) in sorted(stats["by_gene"].items())
        )
        lines.append(
            f"{strategy}\t{stats['current_total']}/360\t{stats['new_total']}/360\t"
            f"{stats['delta']}\t{stats['changed']}\t{stats['improved']}\t"
            f"{stats['regressed']}\t{by_gene_text}\n"
        )
    lines.append("\nchanged_rows\n")
    lines.append("strategy\tsample\tgene\tcurrent\tnew\tdelta\treason\tcurrent_quartet\tnew_quartet\ttruth\n")
    for row in out_rows:
        if row["changed"] != "1":
            continue
        lines.append(
            "\t".join([
                row["strategy"], row["sample"], row["gene"], str(row["current_score"]),
                str(row["new_score"]), str(row["delta"]), row["reason"],
                row["current_quartet"], row["new_quartet"], row["truth_quartet"],
            ]) + "\n"
        )
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text("".join(lines))
    print(f"wrote {args.out_tsv}")
    print(f"wrote {args.summary_out}")


if __name__ == "__main__":
    main()