#!/usr/bin/env python3
"""Plot publication figures for the 15-gene performance projection."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


GENES = (
    "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1",
    "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G",
    "HLA-H", "MICA", "MICB",
)
MIXTURE_COLORS = {
    "10%": "#246A73",
    "20%": "#4D908E",
    "30%": "#90BE6D",
    "40%": "#F9C74F",
    "50%": "#E76F51",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.labelweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.title_fontsize": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def numeric(row: dict[str, str], field: str) -> float:
    return float(row[field])


def group_rows(rows: list[dict[str, str]], field: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row[field]].append(row)
    return grouped


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=400, bbox_inches="tight")


def make_depth_facet_figure(
    gene_rows: list[dict[str, str]],
    value_field: str,
    lower_field: str,
    upper_field: str,
    title: str,
    subtitle: str,
    y_label: str,
    caption: str,
) -> plt.Figure:
    fig, axes = plt.subplots(3, 5, figsize=(11.2, 7.5), sharex=True, sharey=True)
    projection_by_gene = group_rows(gene_rows, "gene")
    for ax, gene in zip(axes.flat, GENES):
        by_mixture = group_rows(projection_by_gene[gene], "mixture_label")
        for mixture, selected in by_mixture.items():
            selected.sort(key=lambda row: int(row["depth"]))
            depth = np.array([int(row["depth"]) for row in selected])
            values = np.array([numeric(row, value_field) for row in selected])
            lower = np.array([numeric(row, lower_field) for row in selected])
            upper = np.array([numeric(row, upper_field) for row in selected])
            color = MIXTURE_COLORS[mixture]
            ax.fill_between(depth, lower, upper, color=color, alpha=0.09, linewidth=0)
            ax.plot(depth, values, color=color, marker="o", markersize=2.8, linewidth=1.35)
        ax.set_title(gene, fontsize=8.5, pad=5, bbox={"facecolor": "#F2F0EA", "edgecolor": "#B8B5AC", "boxstyle": "square,pad=0.28"})
        ax.set_xlim(35, 615)
        ax.set_ylim(0.65, 1.005)
        ax.set_xticks([50, 100, 200, 300, 400, 500, 600])
        ax.set_xticklabels(["50", "100", "200", "300", "400", "500", "600"], rotation=45, ha="right", fontsize=6.5)
        ax.set_yticks([0.70, 0.80, 0.90, 1.00])
        ax.set_yticklabels(["70%", "80%", "90%", "100%"])
        ax.grid(axis="y", color="#E8E6E1", linewidth=0.55)
    fig.suptitle(title, x=0.07, y=0.985, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.07, 0.951, subtitle, fontsize=8.5, color="#404040")
    fig.supxlabel("Sequencing depth", y=0.045, fontsize=9, fontweight="bold")
    fig.supylabel(y_label, x=0.025, fontsize=9, fontweight="bold")
    handles = [
        Line2D([], [], color=color, marker="o", label=mixture, linewidth=1.4, markersize=4)
        for mixture, color in MIXTURE_COLORS.items()
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.97, 0.985), frameon=False, ncol=5, title="Mixture fraction")
    fig.text(
        0.07, 0.012,
        caption,
        fontsize=7.5, color="#555555",
    )
    fig.subplots_adjust(top=0.88, bottom=0.10, left=0.08, right=0.97, hspace=0.42, wspace=0.26)
    return fig


def make_mixture_metric_figure(metric_rows: list[dict[str, str]]) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.7))
    colors = {
        "Pearson r": "#246A73",
        "Spearman rho": "#E76F51",
        "R2": "#577590",
        "CCC": "#7A8F3A",
        "MAE": "#C8553D",
        "RMSE": "#3D5A80",
    }
    groups = (("Pearson r", "Spearman rho", "R2", "CCC"), ("MAE", "RMSE"))
    for ax, metrics in zip(axes, groups):
        for metric in metrics:
            selected = sorted(
                (row for row in metric_rows if row["metric"] == metric),
                key=lambda row: int(row["depth"]),
            )
            depth = np.array([int(row["depth"]) for row in selected])
            values = np.array([numeric(row, "value") for row in selected])
            lower = np.array([numeric(row, "lower_95") for row in selected])
            upper = np.array([numeric(row, "upper_95") for row in selected])
            ax.errorbar(
                depth, values, yerr=np.vstack((values - lower, upper - values)),
                color=colors[metric], marker="o", markersize=4.5, linewidth=1.5,
                capsize=3, label=metric,
            )
        ax.set_xticks([50, 100, 300])
        ax.set_xlabel("Sequencing depth")
        ax.grid(axis="y", color="#E8E6E1", linewidth=0.6)
        ax.legend(frameon=False)
    axes[0].set_ylabel("Agreement")
    axes[0].set_ylim(0.35, 1.01)
    axes[0].set_title("Correlation and concordance", loc="left")
    axes[1].set_ylabel("Absolute fraction error")
    axes[1].set_ylim(0.0, 0.075)
    axes[1].set_title("Estimation error", loc="left")
    fig.suptitle("Global mixture-fraction estimation", x=0.08, y=0.98, ha="left", fontsize=13, fontweight="bold")
    fig.text(0.08, 0.91, "Each depth includes 240 sample-condition estimates; error bars show bootstrap 95% confidence intervals", fontsize=8.3, color="#404040")
    fig.subplots_adjust(top=0.78, bottom=0.17, left=0.09, right=0.98, wspace=0.32)
    return fig


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--mixture-metrics", type=Path)
    args = parser.parse_args()

    configure_style()
    projection = read_tsv(args.data_tsv)
    gene_rows = [row for row in projection if row["scope"] == "gene"]
    recall_figure = make_depth_facet_figure(
        gene_rows,
        "expected_copy_recall", "prediction_lower_95", "prediction_upper_95",
        "Depth response across 15 HLA-related genes",
        "Curves show copy-level recall by mixture fraction; shaded bands indicate 95% model intervals",
        "Copy-level recall",
        "Intervals combine assumed binomial sampling uncertainty with gene-level model uncertainty; no hypothesis tests were performed.",
    )
    save_figure(recall_figure, args.output_dir, "15gene_depth_response")
    plt.close(recall_figure)
    figure_count = 1
    legacy_abundance = args.output_dir / "15gene_abundance_accuracy.png"
    legacy_abundance.unlink(missing_ok=True)
    if args.mixture_metrics is not None:
        metric_figure = make_mixture_metric_figure(read_tsv(args.mixture_metrics))
        save_figure(metric_figure, args.output_dir, "global_mixture_fraction_metrics")
        plt.close(metric_figure)
        figure_count += 1
    print(f"[plot] backend=matplotlib figures={figure_count} format=png output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())