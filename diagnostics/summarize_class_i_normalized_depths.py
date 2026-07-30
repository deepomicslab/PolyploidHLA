#!/usr/bin/env python3
"""Re-gate stored class-I read evidence and summarize accuracy across depths."""

from __future__ import annotations

import argparse
import csv
import glob
import re
from collections import defaultdict
from pathlib import Path


DEPTH_PATTERN = re.compile(r"cov0*([0-9]+)x")
MIXTURE_PATTERN = re.compile(r"graft([0-9]+)")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def depth_from_condition(condition: str) -> int:
    match = DEPTH_PATTERN.search(condition)
    if not match:
        raise ValueError(f"coverage not found in condition: {condition}")
    return int(match.group(1))


def mixture_from_condition(condition: str) -> int:
    match = MIXTURE_PATTERN.search(condition)
    if not match:
        raise ValueError(f"mixture ratio not found in condition: {condition}")
    return int(match.group(1))


def normalize_row(row: dict[str, str], args: argparse.Namespace) -> dict[str, object]:
    baseline = row["baseline_quartet"]
    proposal = row["joint_quartet"]
    baseline_correct = int(row["baseline_correct"])
    proposal_correct = int(row["joint_correct"])
    if baseline == proposal:
        accept = False
        decision = "baseline_same"
        normalized_bf = 0.0
        private_ratio = 1.0
    else:
        informative_pairs = max(1, int(row["informative_pairs"]))
        normalized_bf = float(row["log_bayes_factor"]) / informative_pairs
        baseline_private = int(row["private_pairs_baseline_only"])
        proposal_private = int(row["private_pairs_proposal_only"])
        private_ratio = (proposal_private + 1.0) / (baseline_private + 1.0)
        accept = (
            normalized_bf >= args.min_log_bayes_factor_per_informative_pair
            and int(row["discriminating_pairs"]) >= args.min_discriminating_pairs
            and proposal_private >= args.min_proposal_private_pairs
            and private_ratio >= args.min_private_pair_ratio
            and proposal_private >= baseline_private
        )
        decision = "proposal" if accept else "baseline"
    selected_correct = proposal_correct if accept else baseline_correct
    return {
        **row,
        "depth": depth_from_condition(row["condition"]),
        "normalized_v1_decision": decision,
        "normalized_v1_quartet": proposal if accept else baseline,
        "normalized_v1_correct": selected_correct,
        "normalized_v1_delta": selected_correct - baseline_correct,
        "log_bayes_factor_per_informative_pair": f"{normalized_bf:.6f}",
        "private_pair_ratio": f"{private_ratio:.6f}",
    }


def summarize(rows: list[dict[str, object]], profile: str) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        depth = f"{row['depth']}x"
        mixture = f"{mixture_from_condition(str(row['condition']))}%"
        groups[("overall", "all")].append(row)
        groups[("gene", str(row["gene"]))].append(row)
        groups[("depth", depth)].append(row)
        groups[("mixture", mixture)].append(row)
        groups[("depth_gene", f"{depth}:{row['gene']}")].append(row)
        groups[("gene_mixture", f"{row['gene']}:{mixture}")].append(row)
        groups[("depth_mixture", f"{depth}:{mixture}")].append(row)
        groups[("depth_mixture_gene", f"{depth}:{mixture}:{row['gene']}")].append(row)
        groups[("condition", str(row["condition"]))].append(row)
    scope_order = {
        "overall": 0,
        "gene": 1,
        "depth": 2,
        "mixture": 3,
        "depth_gene": 4,
        "gene_mixture": 5,
        "depth_mixture": 6,
        "depth_mixture_gene": 7,
        "condition": 8,
    }
    output = []
    for (scope, name), selected in sorted(groups.items(), key=lambda item: (scope_order[item[0][0]], item[0][1])):
        baseline = sum(int(row["baseline_correct"]) for row in selected)
        normalized = sum(int(row["normalized_v1_correct"]) for row in selected)
        copies = 4 * len(selected)
        output.append({
            "profile": profile,
            "scope": scope,
            "name": name,
            "sample_loci": len(selected),
            "truth_copies": copies,
            "baseline_correct_copies": baseline,
            "normalized_correct_copies": normalized,
            "delta_correct_copies": normalized - baseline,
            "baseline_copy_recall": f"{baseline / copies:.6f}",
            "normalized_copy_recall": f"{normalized / copies:.6f}",
            "baseline_exact_quartets": sum(int(row["baseline_correct"]) == 4 for row in selected),
            "normalized_exact_quartets": sum(int(row["normalized_v1_correct"]) == 4 for row in selected),
            "normalized_exact_quartet_accuracy": f"{sum(int(row['normalized_v1_correct']) == 4 for row in selected) / len(selected):.6f}",
            "accepted_proposals": sum(row["normalized_v1_decision"] == "proposal" for row in selected),
            "improved_loci": sum(int(row["normalized_v1_delta"]) > 0 for row in selected),
            "regressed_loci": sum(int(row["normalized_v1_delta"]) < 0 for row in selected),
            "neutral_loci": sum(int(row["normalized_v1_delta"]) == 0 for row in selected),
            "errors": sum(bool(row.get("error")) for row in selected),
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", action="append", required=True)
    parser.add_argument("--detail-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    parser.add_argument("--profile", default="class_i_normalized_v1")
    parser.add_argument("--min-log-bayes-factor-per-informative-pair", type=float, default=0.10)
    parser.add_argument("--min-discriminating-pairs", type=int, default=3)
    parser.add_argument("--min-proposal-private-pairs", type=int, default=10)
    parser.add_argument("--min-private-pair-ratio", type=float, default=5.0)
    args = parser.parse_args()

    paths = [Path(path) for pattern in args.input_glob for path in sorted(glob.glob(pattern))]
    if not paths:
        raise RuntimeError("no input files matched")
    rows = [normalize_row(row, args) for path in paths for row in read_tsv(path)]
    rows.sort(key=lambda row: (int(row["depth"]), str(row["experiment"]), str(row["condition"]), str(row["sample"]), str(row["gene"])))
    keys = [(row["experiment"], row["condition"], row["sample"], row["gene"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate per-locus keys")
    write_tsv(args.detail_output, rows)
    write_tsv(args.summary_output, summarize(rows, args.profile))
    print(f"[normalized-depths] loci={len(rows)} detail={args.detail_output}")
    print(f"[normalized-depths] summary={args.summary_output}")


if __name__ == "__main__":
    main()