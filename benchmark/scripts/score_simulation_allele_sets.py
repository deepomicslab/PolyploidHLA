#!/usr/bin/env python3
"""Score simulated four-copy HLA allele sets without source assignment."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_GENES = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")
CONDITION_PATTERN = re.compile(r"graft\d+_cov[^/]+")
SAMPLE_PATTERN = re.compile(r"SIM\d{4}")


def normalize_2field(allele: str) -> str:
    allele = (allele or "NA").replace("HLA-", "").replace("G", "")
    if allele == "NA" or "*" not in allele:
        return allele
    gene, fields = allele.split("*", 1)
    parts = fields.split(":")
    return f"{gene}*{':'.join(parts[:2])}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def path_token(path: Path, pattern: re.Pattern[str]) -> str:
    return next(part for part in path.parts if pattern.fullmatch(part))


def load_truth(bench_root: Path, experiment_glob: str, genes: tuple[str, ...]):
    truth: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    paths = sorted((bench_root / "truth").glob(f"{experiment_glob}/copies.tsv"))
    if not paths:
        raise FileNotFoundError(f"no truth manifests matched {experiment_glob!r}")
    for path in paths:
        for row in read_tsv(path):
            if row["gene"] in genes:
                truth[(row["condition"], row["sample_id"], row["gene"])].append(
                    normalize_2field(row["allele_2field"])
                )
    return truth


def load_final_calls(bench_root: Path, experiment_glob: str, genes: tuple[str, ...]):
    calls: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    pattern = f"{experiment_glob}/*/SIM*/asm_v2/SIM*/*.copy_calls.tsv"
    for path in sorted((bench_root / "runs").glob(pattern)):
        condition = path_token(path, CONDITION_PATTERN)
        sample = path_token(path, SAMPLE_PATTERN)
        for row in read_tsv(path):
            if row["gene"] in genes:
                calls[(condition, sample, row["gene"])].append(
                    normalize_2field(row["allele_2field"])
                )
    return calls


def load_pre_rescue_calls(bench_root: Path, experiment_glob: str, genes: tuple[str, ...]):
    calls: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    pattern = f"{experiment_glob}/*/SIM*/asm_v2/SIM*/hla-*/HLA-*/calls.tsv"
    for final_path in sorted((bench_root / "runs").glob(pattern)):
        gene = final_path.parent.name
        if gene not in genes:
            continue
        condition = path_token(final_path, CONDITION_PATTERN)
        sample = path_token(final_path, SAMPLE_PATTERN)
        rescue_input = final_path.with_name("calls.class2_joint_input.tsv")
        source_path = rescue_input if rescue_input.exists() else final_path
        calls[(condition, sample, gene)] = [
            normalize_2field(row["allele"]) for row in read_tsv(source_path)
        ]
    return calls


def multiset_hits(truth: list[str], prediction: list[str]) -> int:
    return sum((Counter(truth) & Counter(prediction)).values())


def score_stage(stage: str, truth, calls):
    rows = []
    for (condition, sample, gene), truth_alleles in sorted(truth.items()):
        predicted = calls.get((condition, sample, gene), [])
        hits = multiset_hits(truth_alleles, predicted)
        rows.append({
            "stage": stage,
            "condition": condition,
            "sample_id": sample,
            "gene": gene,
            "truth_alleles": ",".join(sorted(truth_alleles)),
            "predicted_alleles": ",".join(sorted(predicted)),
            "called_alleles": len(predicted),
            "correct_alleles": hits,
            "allele_set_recall": f"{hits / len(truth_alleles):.6f}",
            "exact_quartet": int(hits == len(truth_alleles) and len(predicted) == len(truth_alleles)),
        })
    return rows


def summarize(rows: list[dict[str, object]]):
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        stage = str(row["stage"])
        groups[(stage, "overall", "all")].append(row)
        groups[(stage, "condition", str(row["condition"]))].append(row)
        groups[(stage, "gene", str(row["gene"]))].append(row)

    output = []
    for (stage, scope, name), selected in sorted(groups.items()):
        distribution = Counter(int(row["correct_alleles"]) for row in selected)
        correct = sum(int(row["correct_alleles"]) for row in selected)
        truth_copies = 4 * len(selected)
        exact = sum(int(row["exact_quartet"]) for row in selected)
        output.append({
            "stage": stage,
            "scope": scope,
            "name": name,
            "sample_loci": len(selected),
            "correct_alleles": correct,
            "truth_alleles": truth_copies,
            "allele_set_recall": f"{correct / truth_copies:.6f}",
            "exact_quartets": exact,
            "exact_quartet_accuracy": f"{exact / len(selected):.6f}",
            **{f"hit_{hits}_of_4": distribution[hits] for hits in range(5)},
        })
    return output


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench-root", required=True, type=Path)
    parser.add_argument("--experiment-glob", default="accuracy_main_v4_shard*")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES))
    parser.add_argument("--include-pre-rescue", action="store_true")
    args = parser.parse_args()

    genes = tuple(args.genes)
    truth = load_truth(args.bench_root, args.experiment_glob, genes)
    rows = score_stage("final", truth, load_final_calls(args.bench_root, args.experiment_glob, genes))
    if args.include_pre_rescue:
        rows.extend(score_stage(
            "pre_rescue",
            truth,
            load_pre_rescue_calls(args.bench_root, args.experiment_glob, genes),
        ))
    write_tsv(args.out_dir / "allele_set_metrics.tsv", rows)
    write_tsv(args.out_dir / "allele_set_summary.tsv", summarize(rows))
    print(f"[score] sample_loci={len(rows)} stages={len({row['stage'] for row in rows})}")
    print(f"[score] metrics={args.out_dir / 'allele_set_metrics.tsv'}")
    print(f"[score] summary={args.out_dir / 'allele_set_summary.tsv'}")


if __name__ == "__main__":
    main()