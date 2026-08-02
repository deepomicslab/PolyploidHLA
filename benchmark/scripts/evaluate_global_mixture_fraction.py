#!/usr/bin/env python3
"""Evaluate sample-wide mixture-fraction estimates across sequencing depths."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr


METRICS = ("Pearson r", "Spearman rho", "R2", "CCC", "MAE", "RMSE")
MIXTURE_PATTERN = re.compile(r"graft(\d+)")


def read_estimates(path: Path) -> list[dict[str, object]]:
    estimates: dict[tuple[int, str, str, str], float] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            depth = int(row["depth"])
            key = (depth, row["experiment"], row["condition"], row["sample"])
            minority_fraction = 1.0 - float(row["major_fraction_prior"])
            previous = estimates.setdefault(key, minority_fraction)
            if not np.isclose(previous, minority_fraction):
                raise ValueError(f"inconsistent sample-wide estimate for {key}")

    rows = []
    for (depth, experiment, condition, sample), estimate in sorted(estimates.items()):
        match = MIXTURE_PATTERN.search(condition)
        if match is None:
            raise ValueError(f"cannot parse true mixture fraction from {condition}")
        rows.append({
            "depth": depth,
            "experiment": experiment,
            "condition": condition,
            "sample": sample,
            "true_mixture_fraction": int(match.group(1)) / 100.0,
            "estimated_mixture_fraction": estimate,
        })
    return rows


def concordance_correlation(true: np.ndarray, estimated: np.ndarray) -> float:
    covariance = np.cov(true, estimated, ddof=0)[0, 1]
    denominator = true.var() + estimated.var() + (true.mean() - estimated.mean()) ** 2
    return float(2.0 * covariance / denominator)


def metric_values(true: np.ndarray, estimated: np.ndarray) -> dict[str, float]:
    residual = estimated - true
    return {
        "Pearson r": float(pearsonr(true, estimated).statistic),
        "Spearman rho": float(spearmanr(true, estimated).statistic),
        "R2": float(1.0 - np.square(residual).sum() / np.square(true - true.mean()).sum()),
        "CCC": concordance_correlation(true, estimated),
        "MAE": float(np.abs(residual).mean()),
        "RMSE": float(np.sqrt(np.square(residual).mean())),
    }


def bootstrap_intervals(
    true: np.ndarray,
    estimated: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> dict[str, tuple[float, float]]:
    values: dict[str, list[float]] = defaultdict(list)
    completed = 0
    while completed < iterations:
        indices = rng.integers(0, len(true), len(true))
        sampled_true = true[indices]
        if np.unique(sampled_true).size < 2:
            continue
        for metric, value in metric_values(sampled_true, estimated[indices]).items():
            values[metric].append(value)
        completed += 1
    return {
        metric: tuple(float(value) for value in np.percentile(values[metric], [2.5, 97.5]))
        for metric in METRICS
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("per_locus_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260731)
    args = parser.parse_args()

    estimates = read_estimates(args.per_locus_tsv)
    metric_rows = []
    rng = np.random.default_rng(args.seed)
    for depth in sorted({int(row["depth"]) for row in estimates}):
        selected = [row for row in estimates if row["depth"] == depth]
        true = np.array([float(row["true_mixture_fraction"]) for row in selected])
        estimated = np.array([float(row["estimated_mixture_fraction"]) for row in selected])
        point = metric_values(true, estimated)
        intervals = bootstrap_intervals(true, estimated, args.bootstrap_iterations, rng)
        for metric in METRICS:
            lower, upper = intervals[metric]
            metric_rows.append({
                "depth": depth,
                "samples": len(selected),
                "metric": metric,
                "value": f"{point[metric]:.6f}",
                "lower_95": f"{lower:.6f}",
                "upper_95": f"{upper:.6f}",
                "result_status": "empirical_sample_level",
            })

    serializable_estimates = [
        {**row,
         "true_mixture_fraction": f"{float(row['true_mixture_fraction']):.6f}",
         "estimated_mixture_fraction": f"{float(row['estimated_mixture_fraction']):.6f}"}
        for row in estimates
    ]
    write_tsv(args.output_dir / "global_mixture_fraction_estimates.tsv", serializable_estimates)
    write_tsv(args.output_dir / "global_mixture_fraction_metrics.tsv", metric_rows)
    print(f"[mixture] samples={len(estimates)} depths={len({row['depth'] for row in estimates})} metrics={len(metric_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())