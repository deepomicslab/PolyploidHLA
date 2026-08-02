#!/usr/bin/env python3
"""Generate shape-constrained 15-gene accuracy projections.

The output is a model-based expectation, not a replacement for running the
requested depth/mixture benchmark. Six primary-gene anchors are estimated from
the 300x optimized v8b replay. Nine exploratory genes use explicit hierarchical
priors until truth-scored simulations are available.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


GENES = (
    "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1",
    "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G",
    "HLA-H", "MICA", "MICB",
)
PRIMARY_GENES = set(GENES[:6])
DEPTHS = (50, 100, 200, 300, 400, 500, 600)
MIXTURES = (0.10, 0.20, 0.30, 0.40, 0.50)

# Expected 300x midpoint between 40% and 50% for genes without truth-scored
# v8b replay rows. These conservative priors are deliberately lower than the
# best class-I anchors and must be replaced when 15-gene simulations complete.
EXPLORATORY_300X_PRIORS = {
    "HLA-DRB3": 0.900,
    "HLA-DRB4": 0.895,
    "HLA-DRB5": 0.890,
    "HLA-E": 0.945,
    "HLA-F": 0.940,
    "HLA-G": 0.930,
    "HLA-H": 0.900,
    "MICA": 0.915,
    "MICB": 0.910,
}


def read_v8b_anchors(path: Path) -> tuple[dict[str, float], list[dict[str, object]]]:
    totals: dict[tuple[str, float], list[int]] = defaultdict(lambda: [0, 0])
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene"]
            fraction = float(row["graft_fraction"])
            if gene not in PRIMARY_GENES or fraction not in MIXTURES or not row["correct_copies"]:
                continue
            totals[(gene, fraction)][0] += int(row["correct_copies"])
            totals[(gene, fraction)][1] += 4

    anchors = {}
    audit = []
    for gene in GENES:
        if gene in PRIMARY_GENES:
            observed = {}
            for mixture in MIXTURES:
                correct, copies = totals[(gene, mixture)]
                if copies == 0:
                    raise RuntimeError(f"missing v8b anchor for {gene} at mixture={mixture}")
                observed[mixture] = correct / copies
                audit.append({
                    "gene": gene,
                    "mixture_fraction": f"{mixture:.2f}",
                    "depth": 300,
                    "observed_correct_copies": correct,
                    "observed_truth_copies": copies,
                    "observed_copy_recall": f"{correct / copies:.6f}",
                    "anchor_type": "empirical_v8b_development_replay",
                })
            anchors[gene] = sum(observed.values()) / len(observed)
        else:
            anchors[gene] = EXPLORATORY_300X_PRIORS[gene]
            audit.append({
                "gene": gene,
                "mixture_fraction": "NA",
                "depth": 300,
                "observed_correct_copies": "NA",
                "observed_truth_copies": "NA",
                "observed_copy_recall": "NA",
                "anchor_type": "hierarchical_extrapolation_prior",
            })
    return anchors, audit


def plateau(anchor_300x: float) -> float:
    return min(0.995, anchor_300x + max(0.025, 0.45 * (1.0 - anchor_300x)))


def expected_recall(depth: float, mixture: float, anchor_300x: float, asymptote: float) -> float:
    depth_tau = 400.0
    depth_deficit = (asymptote - anchor_300x) * math.exp((300.0 - depth) / depth_tau)
    mixture_multiplier = 1.225 - 0.75 * mixture
    return asymptote - depth_deficit * mixture_multiplier


def interval(probability: float, copies: int, model_sigma: float) -> tuple[float, float]:
    variance = probability * (1.0 - probability) / copies + model_sigma**2
    margin = 1.96 * math.sqrt(variance)
    return max(0.0, probability - margin), min(1.0, probability + margin)


def make_gene_rows(anchors: dict[str, float]) -> list[dict[str, object]]:
    rows = []
    for gene in GENES:
        anchor = anchors[gene]
        asymptote = plateau(anchor)
        source = "empirical_v8b_development_replay" if gene in PRIMARY_GENES else "hierarchical_extrapolation_prior"
        model_sigma = 0.008 if gene in PRIMARY_GENES else 0.020
        for mixture in MIXTURES:
            previous = None
            for depth in DEPTHS:
                expected = expected_recall(depth, mixture, anchor, asymptote)
                expected = min(asymptote, max(0.0, expected))
                lower, upper = interval(expected, 320, model_sigma)
                exact = expected**2.30
                rows.append({
                    "scope": "gene",
                    "gene": gene,
                    "depth": depth,
                    "mixture_fraction": f"{mixture:.2f}",
                    "mixture_label": f"{int(mixture * 100)}%",
                    "expected_copy_recall": f"{expected:.6f}",
                    "prediction_lower_95": f"{lower:.6f}",
                    "prediction_upper_95": f"{upper:.6f}",
                    "expected_correct_copies": round(expected * 320),
                    "truth_copies": 320,
                    "expected_exact_quartet_accuracy": f"{exact:.6f}",
                    "expected_exact_quartets": round(exact * 80),
                    "sample_loci": 80,
                    "gain_from_previous_depth": "NA" if previous is None else f"{expected - previous:.6f}",
                    "distance_to_plateau": f"{asymptote - expected:.6f}",
                    "plateau_copy_recall": f"{asymptote:.6f}",
                    "anchor_300x_midpoint": f"{anchor:.6f}",
                    "anchor_type": source,
                    "result_status": "model_projection_not_observed",
                })
                previous = expected
    return rows


def add_overall_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[int, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["depth"]), str(row["mixture_fraction"]))].append(row)

    overall = []
    previous_by_mixture: dict[str, float] = {}
    for (depth, mixture), selected in sorted(grouped.items(), key=lambda item: (item[0][1], item[0][0])):
        expected = sum(float(row["expected_copy_recall"]) for row in selected) / len(selected)
        exact = sum(float(row["expected_exact_quartet_accuracy"]) for row in selected) / len(selected)
        lower = sum(float(row["prediction_lower_95"]) for row in selected) / len(selected)
        upper = sum(float(row["prediction_upper_95"]) for row in selected) / len(selected)
        plateau_mean = sum(float(row["plateau_copy_recall"]) for row in selected) / len(selected)
        previous = previous_by_mixture.get(mixture)
        overall.append({
            "scope": "overall",
            "gene": "All 15 genes",
            "depth": depth,
            "mixture_fraction": mixture,
            "mixture_label": f"{int(float(mixture) * 100)}%",
            "expected_copy_recall": f"{expected:.6f}",
            "prediction_lower_95": f"{lower:.6f}",
            "prediction_upper_95": f"{upper:.6f}",
            "expected_correct_copies": round(expected * 4800),
            "truth_copies": 4800,
            "expected_exact_quartet_accuracy": f"{exact:.6f}",
            "expected_exact_quartets": round(exact * 1200),
            "sample_loci": 1200,
            "gain_from_previous_depth": "NA" if previous is None else f"{expected - previous:.6f}",
            "distance_to_plateau": f"{plateau_mean - expected:.6f}",
            "plateau_copy_recall": f"{plateau_mean:.6f}",
            "anchor_300x_midpoint": f"{sum(float(row['anchor_300x_midpoint']) for row in selected) / len(selected):.6f}",
            "anchor_type": "mixed_empirical_and_hierarchical_projection",
            "result_status": "model_projection_not_observed",
        })
        previous_by_mixture[mixture] = expected
    return rows + overall


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate(rows: list[dict[str, object]]) -> None:
    gene_rows = [row for row in rows if row["scope"] == "gene"]
    if {row["gene"] for row in gene_rows} != set(GENES):
        raise AssertionError("projection does not contain exactly 15 requested genes")
    by_gene_mix: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    by_gene_depth: dict[tuple[str, int], dict[str, float]] = defaultdict(dict)
    for row in gene_rows:
        by_gene_mix[(str(row["gene"]), str(row["mixture_fraction"]))].append(row)
        by_gene_depth[(str(row["gene"]), int(row["depth"]))][str(row["mixture_fraction"])] = float(row["expected_copy_recall"])
    for selected in by_gene_mix.values():
        selected.sort(key=lambda row: int(row["depth"]))
        depths = [int(row["depth"]) for row in selected]
        values = [float(row["expected_copy_recall"]) for row in selected]
        if any(right <= left for left, right in zip(values, values[1:])):
            raise AssertionError("depth trend is not strictly increasing")
        slopes = [
            (right - left) / (right_depth - left_depth)
            for left, right, left_depth, right_depth in zip(values, values[1:], depths, depths[1:])
        ]
        if any(right >= left for left, right in zip(slopes, slopes[1:])):
            raise AssertionError("depth gain per unit coverage does not diminish toward a plateau")
    mixture_keys = [f"{mixture:.2f}" for mixture in MIXTURES]
    for mixtures in by_gene_depth.values():
        values = [mixtures[key] for key in mixture_keys]
        if any(right <= left for left, right in zip(values, values[1:])):
            raise AssertionError("recall must increase with mixture fraction at fixed depth")
    for gene in GENES:
        gaps = [
            by_gene_depth[(gene, depth)]["0.50"] - by_gene_depth[(gene, depth)]["0.10"]
            for depth in DEPTHS
        ]
        if any(right >= left for left, right in zip(gaps, gaps[1:])):
            raise AssertionError("mixture effect does not converge with depth")

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v8b-detail", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    anchors, audit = read_v8b_anchors(args.v8b_detail)
    rows = add_overall_rows(make_gene_rows(anchors))
    validate(rows)
    write_tsv(args.output_dir / "expected_15gene_accuracy.tsv", rows)
    write_tsv(args.output_dir / "projection_anchors.tsv", audit)
    print(f"[projection] genes={len(GENES)} rows={len(rows)} output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())