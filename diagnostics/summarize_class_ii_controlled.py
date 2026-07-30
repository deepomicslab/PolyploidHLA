#!/usr/bin/env python3
"""Summarize class-II accuracy with depth and mixture ratio held fixed."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


CONDITION_PATTERN = re.compile(r"graft([0-9]+)_cov0*([0-9]+)x")
CLASS_II_GENES = {"HLA-DPB1", "HLA-DQB1", "HLA-DRB1"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, nargs="+")
    parser.add_argument("--mixtures", type=int, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    for input_path in args.input:
        with input_path.open(newline="") as handle:
            rows.extend(
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if row["gene"] in CLASS_II_GENES
            )
    keys = [(row["experiment"], row["condition"], row["sample"], row["gene"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate per-locus keys across class-II inputs")

    groups: dict[tuple[int, int, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        match = CONDITION_PATTERN.fullmatch(row["condition"])
        if not match:
            raise ValueError(f"unsupported condition: {row['condition']}")
        mixture, depth = map(int, match.groups())
        if args.mixtures and mixture not in args.mixtures:
            continue
        groups[(depth, mixture, row["gene"])].append(row)

    output = []
    for (depth, mixture, gene), selected in sorted(groups.items()):
        baseline = sum(int(row["baseline_correct"]) for row in selected)
        strategy_field = "selected_correct" if "selected_correct" in selected[0] else "joint_correct"
        strategy = sum(int(row[strategy_field]) for row in selected)
        truth_copies = 4 * len(selected)
        output.append({
            "scope": "depth_mixture_gene",
            "depth": f"{depth}x",
            "mixture": f"{mixture}%",
            "gene": gene,
            "sample_loci": len(selected),
            "truth_copies": truth_copies,
            "baseline_correct_copies": baseline,
            "strategy_correct_copies": strategy,
            "delta_correct_copies": strategy - baseline,
            "baseline_copy_recall": f"{baseline / truth_copies:.6f}",
            "strategy_copy_recall": f"{strategy / truth_copies:.6f}",
            "improved_loci": sum(int(row["delta_correct"]) > 0 for row in selected),
            "regressed_loci": sum(int(row["delta_correct"]) < 0 for row in selected),
            "neutral_loci": sum(int(row["delta_correct"]) == 0 for row in selected),
            "errors": sum(bool(row.get("error", "")) for row in selected),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(output[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output)
    print(f"[class-ii-controlled] cells={len(output)} output={args.output}")


if __name__ == "__main__":
    main()