#!/usr/bin/env python3
"""Combine class-I read-gated and class-II joint 50x quartet evaluations."""

from __future__ import annotations

import argparse
import csv
import glob
from collections import defaultdict
from pathlib import Path


PRIMARY_GENES = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")


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


def normalize_condition(condition: str) -> str:
    return condition.replace("cov0050x", "cov50x")


def class_i_rows(paths: list[Path], log_bayes_factor_threshold: float) -> list[dict[str, object]]:
    output = []
    for path in paths:
        for row in read_tsv(path):
            if "experiment" not in row:
                raise ValueError(f"class-I detail schema required: {path}")
            output.append({
                "experiment": row["experiment"],
                "condition": normalize_condition(row["condition"]),
                "sample": row["sample"],
                "gene": row["gene"],
                "strategy_component": "class_i_read_bayes_gate",
                "gate_log_bayes_factor_threshold": f"{log_bayes_factor_threshold:.6f}",
                "decision": row["hybrid_decision"],
                "baseline_quartet": row["baseline_quartet"],
                "proposal_quartet": row["joint_quartet"],
                "selected_quartet": row["hybrid_quartet"],
                "truth_quartet": row["truth_quartet"],
                "baseline_correct": int(row["baseline_correct"]),
                "selected_correct": int(row["hybrid_correct"]),
                "delta_correct": int(row["hybrid_correct"]) - int(row["baseline_correct"]),
                "baseline_exact": int(row["baseline_correct"]) == 4,
                "selected_exact": int(row["hybrid_correct"]) == 4,
                "log_bayes_factor": row.get("log_bayes_factor", ""),
                "informative_pairs": row.get("informative_pairs", ""),
                "discriminating_pairs": row.get("discriminating_pairs", ""),
                "proposal_pair_wins": row.get("proposal_pair_wins", ""),
                "private_pairs_baseline_only": row.get("private_pairs_baseline_only", ""),
                "private_pairs_proposal_only": row.get("private_pairs_proposal_only", ""),
                "error": row.get("error", ""),
            })
    return output


def class_ii_rows(path: Path, experiment_prefix: str) -> list[dict[str, object]]:
    output = []
    for row in read_tsv(path):
        if not row["experiment"].startswith(experiment_prefix):
            continue
        output.append({
            "experiment": row["experiment"],
            "condition": normalize_condition(row["condition"]),
            "sample": row["sample"],
            "gene": row["gene"],
            "strategy_component": "class_ii_source_agnostic_joint",
            "gate_log_bayes_factor_threshold": "",
            "decision": "joint",
            "baseline_quartet": row["baseline_quartet"],
            "proposal_quartet": row["joint_quartet"],
            "selected_quartet": row["joint_quartet"],
            "truth_quartet": row["truth_quartet"],
            "baseline_correct": int(row["baseline_correct"]),
            "selected_correct": int(row["joint_correct"]),
            "delta_correct": int(row["joint_correct"]) - int(row["baseline_correct"]),
            "baseline_exact": int(row["baseline_correct"]) == 4,
            "selected_exact": int(row["joint_correct"]) == 4,
            "log_bayes_factor": "",
            "informative_pairs": "",
            "discriminating_pairs": "",
            "proposal_pair_wins": "",
            "private_pairs_baseline_only": "",
            "private_pairs_proposal_only": "",
            "error": "",
        })
    return output


def summarize(rows: list[dict[str, object]], profile: str) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[("overall", "all")].append(row)
        groups[("gene", str(row["gene"]))].append(row)
        groups[("condition", str(row["condition"]))].append(row)
        groups[("component", str(row["strategy_component"]))].append(row)
    output = []
    scope_order = {"overall": 0, "component": 1, "gene": 2, "condition": 3}
    for (scope, name), selected in sorted(groups.items(), key=lambda item: (scope_order[item[0][0]], item[0][1])):
        baseline_correct = sum(int(row["baseline_correct"]) for row in selected)
        strategy_correct = sum(int(row["selected_correct"]) for row in selected)
        truth_copies = 4 * len(selected)
        output.append({
            "profile": profile,
            "scope": scope,
            "name": name,
            "sample_loci": len(selected),
            "truth_copies": truth_copies,
            "baseline_correct_copies": baseline_correct,
            "strategy_correct_copies": strategy_correct,
            "delta_correct_copies": strategy_correct - baseline_correct,
            "baseline_copy_recall": f"{baseline_correct / truth_copies:.6f}",
            "strategy_copy_recall": f"{strategy_correct / truth_copies:.6f}",
            "baseline_exact_quartets": sum(bool(row["baseline_exact"]) for row in selected),
            "strategy_exact_quartets": sum(bool(row["selected_exact"]) for row in selected),
            "strategy_exact_quartet_accuracy": f"{sum(bool(row['selected_exact']) for row in selected) / len(selected):.6f}",
            "changed_loci": sum(row["baseline_quartet"] != row["selected_quartet"] for row in selected),
            "improved_loci": sum(int(row["delta_correct"]) > 0 for row in selected),
            "regressed_loci": sum(int(row["delta_correct"]) < 0 for row in selected),
            "neutral_loci": sum(int(row["delta_correct"]) == 0 for row in selected),
            "errors": sum(bool(row["error"]) for row in selected),
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--class-i-glob", required=True)
    parser.add_argument("--class-ii", required=True, type=Path)
    parser.add_argument("--experiment-prefix", default="accuracy_depth_v1_cov0050_")
    parser.add_argument("--log-bayes-factor-threshold", type=float, default=180.0)
    parser.add_argument("--profile", default="class_i_read_bf180+class2_joint_v2")
    parser.add_argument("--detail-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    args = parser.parse_args()

    class_i_paths = [Path(path) for path in sorted(glob.glob(args.class_i_glob))]
    if not class_i_paths:
        raise RuntimeError(f"no class-I files matched {args.class_i_glob}")
    rows = class_i_rows(class_i_paths, args.log_bayes_factor_threshold)
    rows.extend(class_ii_rows(args.class_ii, args.experiment_prefix))
    rows.sort(key=lambda row: (str(row["experiment"]), str(row["condition"]), str(row["sample"]), str(row["gene"])))
    keys = [(row["experiment"], row["condition"], row["sample"], row["gene"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate per-locus keys in combined evaluation")
    if any(row["gene"] not in PRIMARY_GENES for row in rows):
        raise RuntimeError("non-primary gene in combined evaluation")
    write_tsv(args.detail_output, rows)
    write_tsv(args.summary_output, summarize(rows, args.profile))
    print(f"[50x-summary] loci={len(rows)} detail={args.detail_output}")
    print(f"[50x-summary] summary={args.summary_output}")


if __name__ == "__main__":
    main()