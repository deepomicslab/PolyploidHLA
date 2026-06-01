#!/usr/bin/env python3
"""Score a GenDx validation ASM root against the existing quartet summary."""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


SLOTS = ("R1_2field", "R2_2field", "D1_2field", "D2_2field")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_alleles(text: str) -> list[str]:
    return [item for item in (text or "").split(",") if item]


def read_sample_set(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    samples = set()
    with path.open() as handle:
        for line in handle:
            sample = line.strip()
            if sample:
                samples.add(sample)
    return samples


def norm_2field(allele: str) -> str:
    allele = (allele or "").strip().replace("HLA-", "").replace("G", "").rstrip("P")
    if not allele or allele == "NA" or "*" not in allele:
        return allele or "NA"
    gene, rest = allele.split("*", 1)
    parts = rest.split(":")
    if parts and parts[-1].isalpha():
        parts[-1] = parts[-1][:-1]
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def overlap(truth_vals: list[str], pred_vals: list[str]) -> int:
    counts = Counter(norm_2field(value) for value in truth_vals)
    hits = 0
    for pred in pred_vals:
        value = norm_2field(pred)
        if counts[value] > 0:
            counts[value] -= 1
            hits += 1
    return hits


def score_quartet(pred: list[str], truth_r: list[str], truth_d: list[str]) -> int:
    return overlap(truth_r, pred[:2]) + overlap(truth_d, pred[2:4])


def load_final_calls(path: Path) -> dict[str, list[str]]:
    calls = {}
    for row in read_rows(path):
        calls[row["gene"]] = [norm_2field(row.get(slot, "NA")) for slot in SLOTS]
    return calls


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quartet-summary", required=True, type=Path)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--samples-file", type=Path, default=None)
    args = parser.parse_args()

    sample_filter = read_sample_set(args.samples_file)

    final_cache: dict[str, dict[str, list[str]]] = {}
    detail_rows: list[dict[str, object]] = []
    gene_stats = defaultdict(lambda: Counter(total_rows=0, baseline_score=0, validation_score=0, improved=0, regressed=0, same=0, missing=0))
    all_stats = Counter(total_rows=0, baseline_score=0, validation_score=0, improved=0, regressed=0, same=0, missing=0)

    for row in read_rows(args.quartet_summary):
        sample = row["sample"]
        if sample_filter is not None and sample not in sample_filter:
            continue
        gene = row["gene"]
        if sample not in final_cache:
            final_path = args.asm_root / sample / f"{sample}.final_calls.tsv"
            final_cache[sample] = load_final_calls(final_path)
        pred = final_cache[sample].get(gene)
        truth_r = split_alleles(row["truth_R"])
        truth_d = split_alleles(row["truth_D"])
        baseline_score = int(row["score2"])
        missing = pred is None
        validation_score = 0 if missing else score_quartet(pred, truth_r, truth_d)
        delta = validation_score - baseline_score
        if missing:
            status = "missing_final_call"
        elif delta > 0:
            status = "improved"
        elif delta < 0:
            status = "regressed"
        else:
            status = "same"

        for stats in (all_stats, gene_stats[gene]):
            stats["total_rows"] += 1
            stats["baseline_score"] += baseline_score
            stats["validation_score"] += validation_score
            stats[status if status in {"improved", "regressed", "same"} else "missing"] += 1

        detail_rows.append({
            "sample": sample,
            "set": row["set"],
            "gene": gene,
            "baseline_score": baseline_score,
            "validation_score": validation_score,
            "delta": delta,
            "status": status,
            "truth_R": row["truth_R"],
            "truth_D": row["truth_D"],
            "baseline_pred_R": row["pred_R"],
            "baseline_pred_D": row["pred_D"],
            "validation_pred_R": "" if missing else ",".join(pred[:2]),
            "validation_pred_D": "" if missing else ",".join(pred[2:4]),
        })

    fields = [
        "sample", "set", "gene", "baseline_score", "validation_score", "delta", "status",
        "truth_R", "truth_D", "baseline_pred_R", "baseline_pred_D", "validation_pred_R", "validation_pred_D",
    ]
    write_rows(args.out, fields, detail_rows)

    summary_rows: list[dict[str, object]] = []
    for label, stats in [("ALL", all_stats)] + sorted(gene_stats.items()):
        summary_rows.append({
            "gene": label,
            "rows": stats["total_rows"],
            "baseline_score": stats["baseline_score"],
            "validation_score": stats["validation_score"],
            "delta": stats["validation_score"] - stats["baseline_score"],
            "improved_rows": stats["improved"],
            "regressed_rows": stats["regressed"],
            "same_rows": stats["same"],
            "missing_rows": stats["missing"],
        })
    write_rows(args.summary, [
        "gene", "rows", "baseline_score", "validation_score", "delta",
        "improved_rows", "regressed_rows", "same_rows", "missing_rows",
    ], summary_rows)
    print(f"wrote {args.summary}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()