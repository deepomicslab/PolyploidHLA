#!/usr/bin/env python3
"""Simulate mixture-fraction prediction across depth and mixture conditions."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr


DEPTHS = (50, 100, 200, 300, 400, 500, 600)
MIXTURES = (0.10, 0.20, 0.30, 0.40, 0.50)
GLOBAL_PATTERN = re.compile(r"GLOBAL\s+chi_R=([0-9.]+)")
CONDITION_PATTERN = re.compile(r"graft(\d+)_cov(\d+)x")


def read_unlabeled_300x_anchors(log_root: Path) -> dict[float, np.ndarray]:
    grouped: dict[float, list[float]] = defaultdict(list)
    pattern = "accuracy_main_v4_shard*/graft*_cov300x/SIM*.polyphase.log"
    for path in log_root.glob(pattern):
        condition = CONDITION_PATTERN.fullmatch(path.parent.name)
        if condition is None:
            continue
        mixture = int(condition.group(1)) / 100.0
        if mixture not in MIXTURES:
            continue
        matches = GLOBAL_PATTERN.findall(path.read_text(errors="ignore"))
        if matches:
            chi = float(matches[-1])
            grouped[mixture].append(min(chi, 1.0 - chi))
    return {key: np.asarray(values) for key, values in grouped.items()}


def concordance_correlation(true: np.ndarray, predicted: np.ndarray) -> float:
    covariance = np.cov(true, predicted, ddof=0)[0, 1]
    denominator = true.var() + predicted.var() + (true.mean() - predicted.mean()) ** 2
    return float(2.0 * covariance / denominator)


def simulate(
    anchors: dict[float, np.ndarray], replicates: int, seed: int
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rng = np.random.default_rng(seed)
    condition_rows = []
    depth_values: dict[int, tuple[list[float], list[float]]] = {
        depth: ([], []) for depth in DEPTHS
    }
    sampled_residuals = {
        mixture: rng.choice(anchors[mixture] - mixture, size=replicates, replace=True)
        for mixture in MIXTURES
    }
    for depth in DEPTHS:
        for mixture in MIXTURES:
            predicted = np.clip(
                mixture + sampled_residuals[mixture] * np.sqrt(300.0 / depth), 0.0, 0.5
            )
            error = predicted - mixture
            depth_values[depth][0].extend([mixture] * replicates)
            depth_values[depth][1].extend(predicted.tolist())
            condition_rows.append({
                "depth": depth,
                "mixture_fraction": f"{mixture:.2f}",
                "replicates": replicates,
                "mean_predicted_fraction": f"{predicted.mean():.6f}",
                "prediction_lower_95": f"{np.quantile(predicted, 0.025):.6f}",
                "prediction_upper_95": f"{np.quantile(predicted, 0.975):.6f}",
                "bias": f"{error.mean():.6f}",
                "median_error": f"{np.median(error):.6f}",
                "error_q25": f"{np.quantile(error, 0.25):.6f}",
                "error_q75": f"{np.quantile(error, 0.75):.6f}",
                "error_sd": f"{error.std(ddof=1):.6f}",
                "error_lower_95": f"{np.quantile(error, 0.025):.6f}",
                "error_upper_95": f"{np.quantile(error, 0.975):.6f}",
                "mae": f"{np.abs(error).mean():.6f}",
                "rmse": f"{np.sqrt(np.square(error).mean()):.6f}",
                "within_0.02": f"{np.mean(np.abs(error) <= 0.02):.6f}",
                "within_0.05": f"{np.mean(np.abs(error) <= 0.05):.6f}",
                "anchor_source": "unlabeled_global_chi_300x_residuals",
                "result_status": "model_based_monte_carlo",
            })

    metric_rows = []
    for depth, (true_values, predicted_values) in depth_values.items():
        true = np.asarray(true_values)
        predicted = np.asarray(predicted_values)
        residual = predicted - true
        metrics = {
            "Pearson r": pearsonr(true, predicted).statistic,
            "Spearman rho": spearmanr(true, predicted).statistic,
            "R2": 1.0 - np.square(residual).sum() / np.square(true - true.mean()).sum(),
            "CCC": concordance_correlation(true, predicted),
            "MAE": np.abs(residual).mean(),
            "RMSE": np.sqrt(np.square(residual).mean()),
        }
        for metric, value in metrics.items():
            metric_rows.append({"depth": depth, "metric": metric, "value": f"{value:.6f}"})
    return condition_rows, metric_rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_results(
    condition_rows: list[dict[str, object]],
    metric_rows: list[dict[str, object]],
    anchors: dict[float, np.ndarray],
    path: Path,
) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))
    colors = plt.cm.viridis(np.linspace(0.08, 0.88, len(DEPTHS)))
    for depth, color in zip(DEPTHS, colors):
        selected = [row for row in condition_rows if row["depth"] == depth]
        truth = np.array([float(row["mixture_fraction"]) for row in selected])
        mean = np.array([float(row["mean_predicted_fraction"]) for row in selected])
        axes[0, 0].plot(truth, mean, color=color, marker="o", linewidth=1.25, markersize=3.5, label=f"{depth}x")
    axes[0, 0].plot([0.1, 0.5], [0.1, 0.5], color="#555555", linestyle="--", linewidth=1, label="Identity")
    axes[0, 0].set(xlabel="True mixture fraction", ylabel="Mean predicted fraction", title="Calibration by condition")
    axes[0, 0].set_xticks(MIXTURES)
    axes[0, 0].legend(frameon=False, ncol=2, fontsize=7)

    def heatmap(ax: plt.Axes, field: str, title: str, cmap: str, vmin: float, vmax: float, percent: bool = False) -> None:
        matrix = np.array([
            [float(next(row[field] for row in condition_rows if row["depth"] == depth and float(row["mixture_fraction"]) == mixture)) for depth in DEPTHS]
            for mixture in MIXTURES
        ])
        image = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                label = f"{matrix[row_index, column_index] * 100:.0f}%" if percent else f"{matrix[row_index, column_index]:.3f}"
                ax.text(column_index, row_index, label, ha="center", va="center", fontsize=7,
                        color="white" if matrix[row_index, column_index] > (vmin + vmax) / 2 else "#202020")
        ax.set_xticks(range(len(DEPTHS)), DEPTHS)
        ax.set_yticks(range(len(MIXTURES)), [f"{int(value * 100)}%" for value in MIXTURES])
        ax.set(xlabel="Sequencing depth", ylabel="True mixture fraction", title=title)
        fig.colorbar(image, ax=ax, fraction=0.045, pad=0.03)

    heatmap(axes[0, 1], "mae", "Mean absolute error", "magma_r", 0.0, 0.24)
    mixture_colors = ("#246A73", "#4D908E", "#90BE6D", "#F9C74F", "#E76F51")
    for mixture, color in zip(MIXTURES, mixture_colors):
        absolute_error = np.sort(np.abs(anchors[mixture] - mixture))
        cumulative_probability = np.arange(1, len(absolute_error) + 1) / len(absolute_error)
        axes[1, 0].step(
            absolute_error, cumulative_probability, where="post",
            color=color, linewidth=1.8, label=f"{int(mixture * 100)}%",
        )
    axes[1, 0].axvline(0.05, color="#555555", linestyle="--", linewidth=1)
    axes[1, 0].set(
        xlabel="Absolute prediction error", ylabel="Cumulative fraction",
        title="Error distributions by mixture at 300x",
    )
    axes[1, 0].set_xlim(0.0, 0.50)
    axes[1, 0].set_ylim(0.0, 1.02)
    axes[1, 0].legend(frameon=False, ncol=3, fontsize=7, title="Mixture fraction")

    agreement_metrics = ("Pearson r", "Spearman rho", "R2", "CCC")
    metric_colors = ("#246A73", "#E76F51", "#577590", "#7A8F3A")
    for metric, color in zip(agreement_metrics, metric_colors):
        selected = [row for row in metric_rows if row["metric"] == metric]
        axes[1, 1].plot(
            [int(row["depth"]) for row in selected], [float(row["value"]) for row in selected],
            color=color, marker="o", linewidth=1.4, markersize=3.5, label=metric,
        )
    axes[1, 1].set(xlabel="Sequencing depth", ylabel="Agreement", title="Cross-mixture agreement")
    axes[1, 1].set_xticks(DEPTHS)
    axes[1, 1].set_ylim(-0.5, 1.02)
    axes[1, 1].legend(frameon=False, fontsize=7)
    for ax in axes.flat:
        ax.grid(axis="y", color="#E8E6E1", linewidth=0.55)
    fig.suptitle("Mixture-fraction prediction performance", x=0.07, y=0.985, ha="left", fontsize=14, fontweight="bold")
    fig.subplots_adjust(top=0.91, bottom=0.08, left=0.08, right=0.97, hspace=0.38, wspace=0.30)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--replicates", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260731)
    args = parser.parse_args()

    anchors = read_unlabeled_300x_anchors(args.log_root)
    if set(anchors) != set(MIXTURES):
        raise RuntimeError(f"missing abundance anchors: {sorted(set(MIXTURES) - set(anchors))}")
    condition_rows, metric_rows = simulate(anchors, args.replicates, args.seed)
    write_tsv(args.output_dir / "abundance_prediction_conditions.tsv", condition_rows)
    write_tsv(args.output_dir / "abundance_prediction_metrics.tsv", metric_rows)
    plot_results(
        condition_rows, metric_rows, anchors,
        args.output_dir / "abundance_prediction_performance.png",
    )
    print(f"[abundance] conditions={len(condition_rows)} depths={len(DEPTHS)} mixtures={len(MIXTURES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())