#!/usr/bin/env python3
"""Evaluate current mixture-fraction prediction on existing simulations."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr

DEPTHS = (50, 100, 300)
MIXTURES = (0.10, 0.20, 0.30, 0.40, 0.50)
CONDITION_PATTERN = re.compile(r"graft(\d+)_cov0*(\d+)x")
CHI_PATTERN = re.compile(r"\bchi_R=([0-9.]+)")
GLOBAL_PATTERN = re.compile(r"^GLOBAL\s+(.*)$", re.MULTILINE)
KEY_VALUE_PATTERN = re.compile(r"(\w+)=([^\s]+)")
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def discover_cases(run_root: Path) -> list[tuple[Path, int, float, Path]]:
    patterns = {
        50: "accuracy_depth_v1_cov0050_shard*/graft*_cov50x/SIM*/spechla_out/SIM*/*.pooled_continuous.vcf.gz",
        100: "accuracy_depth_v1_cov0100_shard*/graft*_cov100x/SIM*/spechla_out/SIM*/*.pooled_continuous.vcf.gz",
        300: "accuracy_main_v4_shard*/graft*_cov300x/SIM*/spechla_out/SIM*/*.pooled_continuous.vcf.gz",
    }
    cases = []
    for depth, pattern in patterns.items():
        for vcf in run_root.glob(pattern):
            condition = next(
                (CONDITION_PATTERN.fullmatch(part) for part in vcf.parts
                 if CONDITION_PATTERN.fullmatch(part)),
                None,
            )
            if condition is None:
                continue
            mixture = int(condition.group(1)) / 100.0
            sample = vcf.name.split(".pooled_continuous", 1)[0]
            prior_path = vcf.parent / f"{sample}.chimerism.txt"
            if mixture in MIXTURES and prior_path.exists():
                cases.append((vcf, depth, mixture, prior_path))
    return sorted(cases, key=lambda item: (item[1], item[2], str(item[0])))


def parse_prior(path: Path) -> float:
    matches = CHI_PATTERN.findall(path.read_text(errors="replace"))
    if not matches:
        raise RuntimeError(f"missing chi_R prior: {path}")
    return float(matches[-1])


def evaluate_case(
    case: tuple[Path, int, float, Path], estimator: Path, bootstrap: int,
) -> dict[str, object]:
    vcf, depth, mixture, prior_path = case
    command = [
        sys.executable, str(estimator), str(vcf), "--prior-chi", str(parse_prior(prior_path)),
        "--bootstrap", str(bootstrap),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    match = GLOBAL_PATTERN.search(completed.stdout)
    if match is None:
        raise RuntimeError(f"missing GLOBAL output: {vcf}")
    fields = dict(KEY_VALUE_PATTERN.findall(match.group(1)))
    raw = None if fields.get("raw_chi_R") == "NA" else float(fields["raw_chi_R"])
    reported = None if fields.get("chi_R") == "NA" else float(fields["chi_R"])
    return {
        "depth": depth,
        "mixture_fraction": mixture,
        "sample": vcf.name.split(".pooled_continuous", 1)[0],
        "status": fields.get("status", "UNKNOWN"),
        "reasons": fields.get("reasons", "unknown"),
        "raw_chi_R": raw,
        "reported_chi_R": reported,
        "raw_absolute_error": None if raw is None else abs(raw - mixture),
        "reported_absolute_error": None if reported is None else abs(reported - mixture),
        "vcf": str(vcf),
    }


def concordance_correlation(true: np.ndarray, predicted: np.ndarray) -> float:
    covariance = np.cov(true, predicted, ddof=0)[0, 1]
    denominator = true.var() + predicted.var() + (true.mean() - predicted.mean()) ** 2
    return float(2.0 * covariance / denominator)


def summarize_conditions(cases: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[int, float], list[dict[str, object]]] = defaultdict(list)
    for case in cases:
        grouped[(int(case["depth"]), float(case["mixture_fraction"]))].append(case)
    rows = []
    for (depth, mixture), selected in sorted(grouped.items()):
        raw = np.asarray([case["raw_chi_R"] for case in selected if case["raw_chi_R"] is not None])
        passed = np.asarray([
            case["reported_chi_R"] for case in selected if case["reported_chi_R"] is not None
        ])
        raw_error = np.abs(raw - mixture)
        pass_error = np.abs(passed - mixture)
        status_counts = Counter(str(case["status"]) for case in selected)
        rows.append({
            "depth": depth,
            "mixture_fraction": f"{mixture:.2f}",
            "total": len(selected),
            "raw_callable": len(raw),
            "pass": len(passed),
            "low_confidence": status_counts["LOW_CONFIDENCE"],
            "model_mismatch": status_counts["MODEL_MISMATCH"],
            "call_rate": f"{len(passed) / len(selected):.6f}",
            "mean_predicted_fraction": f"{passed.mean():.6f}" if len(passed) else "NA",
            "prediction_lower_95": f"{np.quantile(passed, 0.025):.6f}" if len(passed) else "NA",
            "prediction_upper_95": f"{np.quantile(passed, 0.975):.6f}" if len(passed) else "NA",
            "raw_mae": f"{raw_error.mean():.6f}" if len(raw_error) else "NA",
            "pass_mae": f"{pass_error.mean():.6f}" if len(pass_error) else "NA",
            "raw_within_0.05": f"{np.mean(raw_error <= 0.05):.6f}" if len(raw_error) else "NA",
            "pass_within_0.05": f"{np.mean(pass_error <= 0.05):.6f}" if len(pass_error) else "NA",
            "raw_severe_rate": f"{np.mean(raw_error > 0.10):.6f}" if len(raw_error) else "NA",
            "pass_severe_rate": f"{np.mean(pass_error > 0.10):.6f}" if len(pass_error) else "NA",
            "anchor_source": "current_estimator_empirical_vcf_replay",
            "result_status": "development_qc_replay",
        })
    return rows


def summarize_metrics(cases: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for depth in DEPTHS:
        depth_cases = [case for case in cases if case["depth"] == depth]
        selected = [case for case in depth_cases if case["reported_chi_R"] is not None]
        true = np.asarray([case["mixture_fraction"] for case in selected], dtype=float)
        predicted = np.asarray([case["reported_chi_R"] for case in selected], dtype=float)
        residual = predicted - true
        metrics = {
            "Pearson r": pearsonr(true, predicted).statistic,
            "Spearman rho": spearmanr(true, predicted).statistic,
            "R2": 1.0 - np.square(residual).sum() / np.square(true - true.mean()).sum(),
            "CCC": concordance_correlation(true, predicted),
            "MAE": np.abs(residual).mean(),
            "RMSE": np.sqrt(np.square(residual).mean()),
            "Call rate": len(selected) / len(depth_cases),
        }
        rows.extend(
            {"depth": depth, "metric": metric, "value": f"{value:.6f}",
             "result_status": "development_qc_replay"}
            for metric, value in metrics.items()
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_results(
    cases: list[dict[str, object]], condition_rows: list[dict[str, object]], path: Path,
) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))
    colors = {50: "#355C7D", 100: "#2A9D8F", 300: "#E76F51"}

    for depth in DEPTHS:
        selected = [row for row in condition_rows if row["depth"] == depth]
        truth = np.asarray([float(row["mixture_fraction"]) for row in selected])
        mean = np.asarray([float(row["mean_predicted_fraction"]) for row in selected])
        lower = np.asarray([float(row["prediction_lower_95"]) for row in selected])
        upper = np.asarray([float(row["prediction_upper_95"]) for row in selected])
        axes[0, 0].errorbar(
            truth, mean, yerr=[mean - lower, upper - mean], color=colors[depth],
            marker="o", linewidth=1.3, markersize=4, capsize=2, label=f"{depth}x",
        )
    axes[0, 0].plot([0.1, 0.5], [0.1, 0.5], color="#555555", linestyle="--", linewidth=1)
    axes[0, 0].set(
        xlabel="True mixture fraction", ylabel="Mean PASS estimate",
        title="Empirical calibration (95% range)", xlim=(0.075, 0.525), ylim=(0.05, 0.525),
    )
    axes[0, 0].legend(frameon=False)

    matrix = np.full((len(MIXTURES), len(DEPTHS)), np.nan)
    labels = {}
    for row in condition_rows:
        row_index = MIXTURES.index(float(row["mixture_fraction"]))
        column_index = DEPTHS.index(int(row["depth"]))
        matrix[row_index, column_index] = float(row["pass_mae"])
        labels[(row_index, column_index)] = f"{float(row['pass_mae']):.3f}\n{row['pass']}/{row['total']}"
    image = axes[0, 1].imshow(matrix, cmap="YlOrRd", aspect="auto", vmin=0.0, vmax=0.06)
    for (row_index, column_index), label in labels.items():
        axes[0, 1].text(column_index, row_index, label, ha="center", va="center", fontsize=7)
    axes[0, 1].set_xticks(range(len(DEPTHS)), DEPTHS)
    axes[0, 1].set_yticks(range(len(MIXTURES)), [f"{int(value * 100)}%" for value in MIXTURES])
    axes[0, 1].set(xlabel="Sequencing depth", ylabel="True mixture fraction", title="PASS MAE and calls / total")
    fig.colorbar(image, ax=axes[0, 1], fraction=0.045, pad=0.03)

    mixture_colors = ("#264653", "#2A9D8F", "#8AB17D", "#E9C46A", "#E76F51")
    for mixture, color in zip(MIXTURES, mixture_colors):
        errors = sorted(
            float(case["reported_absolute_error"]) for case in cases
            if case["depth"] == 300 and case["mixture_fraction"] == mixture
            and case["reported_absolute_error"] is not None
        )
        cumulative = np.arange(1, len(errors) + 1) / len(errors)
        axes[1, 0].step(errors, cumulative, where="post", color=color, linewidth=1.7,
                        label=f"{int(mixture * 100)}%")
    axes[1, 0].axvline(0.05, color="#555555", linestyle="--", linewidth=1)
    axes[1, 0].set(
        xlabel="Absolute prediction error", ylabel="Cumulative fraction",
        title="PASS error distributions at 300x", xlim=(0.0, 0.25), ylim=(0.0, 1.02),
    )
    axes[1, 0].legend(frameon=False, ncol=3, fontsize=7)

    statuses = ("PASS", "LOW_CONFIDENCE", "MODEL_MISMATCH", "RAW_NA")
    status_colors = ("#2A9D8F", "#E9C46A", "#E76F51", "#8D99AE")
    bottoms = np.zeros(len(DEPTHS))
    for status, color in zip(statuses, status_colors):
        values = []
        for depth in DEPTHS:
            depth_cases = [case for case in cases if case["depth"] == depth]
            if status == "RAW_NA":
                count = sum(case["raw_chi_R"] is None for case in depth_cases)
            else:
                count = sum(case["raw_chi_R"] is not None and case["status"] == status for case in depth_cases)
            values.append(count / len(depth_cases))
        axes[1, 1].bar(range(len(DEPTHS)), values, bottom=bottoms, color=color, label=status)
        bottoms += np.asarray(values)
    axes[1, 1].set_xticks(range(len(DEPTHS)), DEPTHS)
    axes[1, 1].set(
        xlabel="Sequencing depth", ylabel="Fraction of samples",
        title="QC disposition and reporting coverage", ylim=(0.0, 1.0),
    )
    axes[1, 1].legend(frameon=False, fontsize=7, loc="lower right")

    for ax in axes.flat:
        ax.grid(axis="y", color="#E8E6E1", linewidth=0.55)
    fig.suptitle("Mixture-fraction prediction: current estimator", x=0.07, y=0.985,
                 ha="left", fontsize=14, fontweight="bold")
    fig.text(0.07, 0.012, "Empirical development replay; QC thresholds require independent validation.",
             fontsize=7.5, color="#555555")
    fig.subplots_adjust(top=0.91, bottom=0.09, left=0.08, right=0.97, hspace=0.38, wspace=0.30)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--estimator", type=Path, default=PROJECT_ROOT / "estimate_chi_pooled.py")
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    discovered = discover_cases(args.run_root)
    if not discovered:
        raise RuntimeError(f"no pooled-continuous VCFs found under {args.run_root}")
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        cases = list(executor.map(
            lambda case: evaluate_case(case, args.estimator, args.bootstrap), discovered,
        ))
    condition_rows = summarize_conditions(cases)
    metric_rows = summarize_metrics(cases)
    serializable_cases = [
        {key: "NA" if value is None else value for key, value in case.items()}
        for case in cases
    ]
    write_tsv(args.output_dir / "abundance_prediction_cases.tsv", serializable_cases)
    write_tsv(args.output_dir / "abundance_prediction_conditions.tsv", condition_rows)
    write_tsv(args.output_dir / "abundance_prediction_metrics.tsv", metric_rows)
    plot_results(cases, condition_rows, args.output_dir / "abundance_prediction_performance.png")
    print(f"[abundance] cases={len(cases)} conditions={len(condition_rows)} depths={DEPTHS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
