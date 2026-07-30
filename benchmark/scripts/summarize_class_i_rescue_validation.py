#!/usr/bin/env python3
"""Summarize paired baseline/rescue class-I validation runs."""

from __future__ import annotations

import argparse
import csv
import glob
import json
from collections import Counter, defaultdict
from pathlib import Path

GENES = ("HLA-A", "HLA-B", "HLA-C")


def normalize_2field(allele: str) -> str:
    allele = allele.strip().replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, fields = allele.rstrip("GP").split("*", 1)
    parts = fields.split(":")
    return f"{gene}*{parts[0]}:{parts[1]}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def multiset_hits(truth: list[str], prediction: list[str]) -> int:
    return sum((Counter(truth) & Counter(prediction)).values())


def load_truth(bench_root: Path, experiment_glob: str):
    truth = defaultdict(list)
    scenarios = {}
    experiment_paths = sorted((bench_root / "truth").glob(f"{experiment_glob}/copies.tsv"))
    if not experiment_paths:
        raise FileNotFoundError(f"no truth manifests matched {experiment_glob!r}")
    for path in experiment_paths:
        experiment = path.parent.name
        design_path = bench_root / "config" / experiment / "design.json"
        with design_path.open() as handle:
            scenarios[experiment] = json.load(handle)["scenario"]
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["gene"] in GENES:
                    key = (experiment, row["condition"], row["sample_id"], row["gene"])
                    truth[key].append(normalize_2field(row["allele_2field"]))
    return truth, scenarios


def load_predictions(bench_root: Path, experiments: set[str], suffix: str):
    predictions = defaultdict(list)
    for experiment in experiments:
        run_experiment = experiment + suffix
        pattern = bench_root / "runs" / run_experiment / "*" / "SIM*" / "asm_v2" / "SIM*" / "*.copy_calls.tsv"
        for path_text in glob.glob(str(pattern)):
            path = Path(path_text)
            sample = path.parent.name
            condition = path.parents[3].name
            with path.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["gene"] in GENES:
                        key = (experiment, condition, sample, row["gene"])
                        predictions[key].append(normalize_2field(row["allele_2field"]))
    return predictions


def load_rescue_audit(bench_root: Path, experiments: set[str]):
    audit = {}
    for experiment in experiments:
        pattern = bench_root / "runs" / experiment / "*" / "SIM*" / "spechla_out" / "SIM*" / "em_refine" / "HLA-?.summary.tsv"
        for path_text in glob.glob(str(pattern)):
            path = Path(path_text)
            condition = path.parents[4].name
            sample = path.parents[3].name
            gene = path.name.removesuffix(".summary.tsv")
            if gene not in GENES:
                continue
            with path.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"), {})
            audit[(experiment, condition, sample, gene)] = {
                "applied": row.get("class_i_distinct_rescue", "0"),
                "detail": row.get("class_i_distinct_detail", "."),
            }
    return audit


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench-root", type=Path, required=True)
    parser.add_argument("--experiment-glob", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    truth, scenarios = load_truth(args.bench_root, args.experiment_glob)
    experiments = set(scenarios)
    baseline = load_predictions(args.bench_root, experiments, "_baseline")
    rescue = load_predictions(args.bench_root, experiments, "")
    audit = load_rescue_audit(args.bench_root, experiments)
    rows = []
    aggregate = defaultdict(Counter)

    for key in sorted(truth):
        experiment, condition, sample, gene = key
        expected = truth[key]
        baseline_call = baseline.get(key, [])
        rescue_call = rescue.get(key, [])
        baseline_hits = multiset_hits(expected, baseline_call)
        rescue_hits = multiset_hits(expected, rescue_call)
        truth_distinct = len(set(expected))
        baseline_distinct = len(set(baseline_call))
        rescue_distinct = len(set(rescue_call))
        audit_row = audit.get(key, {"applied": "0", "detail": "."})
        row = {
            "experiment": experiment,
            "scenario": scenarios[experiment],
            "condition": condition,
            "sample_id": sample,
            "gene": gene,
            "truth_hits_denominator": 4,
            "baseline_hits": baseline_hits,
            "rescue_hits": rescue_hits,
            "delta_hits": rescue_hits - baseline_hits,
            "truth_distinct": truth_distinct,
            "baseline_distinct": baseline_distinct,
            "rescue_distinct": rescue_distinct,
            "baseline_oversplit": int(baseline_distinct > truth_distinct),
            "rescue_oversplit": int(rescue_distinct > truth_distinct),
            "rescue_logged": audit_row["applied"],
            "rescue_detail": audit_row["detail"],
            "baseline_missing": int(len(baseline_call) != 4),
            "rescue_missing": int(len(rescue_call) != 4),
        }
        rows.append(row)
        for group in ((scenarios[experiment], gene), (scenarios[experiment], "ABC"), ("ALL", "ABC")):
            stats = aggregate[group]
            stats["loci"] += 1
            stats["baseline_hits"] += baseline_hits
            stats["rescue_hits"] += rescue_hits
            stats["improved"] += rescue_hits > baseline_hits
            stats["regressed"] += rescue_hits < baseline_hits
            stats["baseline_oversplit"] += baseline_distinct > truth_distinct
            stats["rescue_oversplit"] += rescue_distinct > truth_distinct
            stats["rescue_logged"] += audit_row["applied"] == "1"
            stats["baseline_missing"] += len(baseline_call) != 4
            stats["rescue_missing"] += len(rescue_call) != 4

    write_tsv(args.output, rows)
    print("scenario\tgene\tloci\tbaseline_recall\trescue_recall\timproved\tregressed\tbaseline_oversplit\trescue_oversplit\trescue_logged\tbaseline_missing\trescue_missing")
    for (scenario, gene), stats in sorted(aggregate.items()):
        denominator = 4 * stats["loci"]
        print(
            f"{scenario}\t{gene}\t{stats['loci']}\t"
            f"{stats['baseline_hits'] / denominator:.4f}\t"
            f"{stats['rescue_hits'] / denominator:.4f}\t"
            f"{stats['improved']}\t{stats['regressed']}\t"
            f"{stats['baseline_oversplit']}\t{stats['rescue_oversplit']}\t"
            f"{stats['rescue_logged']}\t{stats['baseline_missing']}\t{stats['rescue_missing']}"
        )
    print(f"details\t{args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
