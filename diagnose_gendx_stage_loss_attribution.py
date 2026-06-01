#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

SLOTS = ("R1_2field", "R2_2field", "D1_2field", "D2_2field")
SIDES = ("R", "R", "D", "D")


def read_tsv(path: Path):
    if not path or not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_alleles(text: str):
    return [item for item in (text or "").split(",") if item]


def side_misses(truth_r, truth_d, pred_r, pred_d):
    misses = []
    for side, truth, pred in (("R", truth_r, pred_r), ("D", truth_d, pred_d)):
        pred_counts = Counter(pred)
        for allele in truth:
            if pred_counts[allele] > 0:
                pred_counts[allele] -= 1
            else:
                misses.append((side, allele))
    return misses


def load_final_by_key(asm_root: Path):
    out = {}
    for sample_dir in asm_root.iterdir():
        if not sample_dir.is_dir():
            continue
        final_path = sample_dir / f"{sample_dir.name}.final_calls.tsv"
        for row in read_tsv(final_path):
            out[(sample_dir.name, row.get("gene", ""))] = row
    return out


def load_tf_counts(spechla_root: Path, sample: str, gene: str):
    path = spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
    rows = read_tsv(path)
    out = {}
    for index, row in enumerate(rows, 1):
        allele = row.get("allele_2field", "")
        if not allele:
            continue
        try:
            frac = float(row.get("fraction") or row.get("em_frac") or 0.0)
        except ValueError:
            frac = 0.0
        try:
            weight = float(row.get("em_weight") or 0.0)
        except ValueError:
            weight = 0.0
        out[allele] = {"rank": index, "frac": frac, "weight": weight}
    return out


def load_raw_support(path: Path | None):
    out = {}
    if path is None:
        return out
    for row in read_tsv(path):
        key = (row.get("sample", ""), row.get("gene", ""), row.get("side", ""), row.get("missing_allele", ""))
        out[key] = row
    return out


def classify_miss(allele: str, side: str, pred_r, pred_d, tf_row, raw_row):
    opposite = pred_d if side == "R" else pred_r
    same = pred_r if side == "R" else pred_d
    if allele in opposite or same.count(allele) > 0:
        return "side_or_copy_assignment"
    raw_class = raw_row.get("raw_support_class", "") if raw_row else ""
    if raw_class == "binning_loss":
        return "read_binning_loss"
    if raw_class in {"low_raw_support", "no_raw_support_by_unique_kmers"}:
        return raw_class
    if tf_row and tf_row["frac"] > 0:
        return "candidate_available_not_selected"
    return "absent_from_em_candidate"


def score_from_summary(rows):
    return sum(int(row.get("score2", "0") or 0) for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quartet-summary", required=True, type=Path)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--raw-support", type=Path, default=None)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary-out", required=True, type=Path)
    args = parser.parse_args()

    summary_rows = read_tsv(args.quartet_summary)
    final_by_key = load_final_by_key(args.asm_root)
    raw_support = load_raw_support(args.raw_support)
    tf_cache = {}
    out_rows = []

    for row in summary_rows:
        sample = row["sample"]
        gene = row["gene"]
        truth_r = split_alleles(row["truth_R"])
        truth_d = split_alleles(row["truth_D"])
        pred_r = split_alleles(row["pred_R"])
        pred_d = split_alleles(row["pred_D"])
        final = final_by_key.get((sample, gene), {})
        tf_key = (sample, gene)
        if tf_key not in tf_cache:
            tf_cache[tf_key] = load_tf_counts(args.spechla_root, sample, gene)
        tf_counts = tf_cache[tf_key]
        pred_total_counts = Counter(pred_r + pred_d)
        truth_total_counts = Counter(truth_r + truth_d)
        for side, allele in side_misses(truth_r, truth_d, pred_r, pred_d):
            tf_row = tf_counts.get(allele, {})
            raw_row = raw_support.get((sample, gene, side, allele), {})
            stage = classify_miss(allele, side, pred_r, pred_d, tf_row, raw_row)
            out_rows.append({
                "set": row.get("set", ""),
                "sample": sample,
                "gene": gene,
                "side": side,
                "missing_allele": allele,
                "current_score": row.get("score2", ""),
                "stage_class": stage,
                "pred_total_count": pred_total_counts[allele],
                "truth_total_count": truth_total_counts[allele],
                "em_rank": tf_row.get("rank", ""),
                "em_frac": f"{tf_row.get('frac', 0.0):.8f}" if tf_row else "",
                "em_weight": f"{tf_row.get('weight', 0.0):.3f}" if tf_row else "",
                "raw_support_class": raw_row.get("raw_support_class", ""),
                "full_support_pairs_k31": raw_row.get("full_support_pairs_k31", ""),
                "binned_support_pairs_k31": raw_row.get("binned_support_pairs_k31", ""),
                "retained_fraction_k31": raw_row.get("retained_fraction_k31", ""),
                "copy_identifiability": final.get("copy_identifiability", ""),
                "copy_fit_error": final.get("copy_fit_error", ""),
                "mean_mask_fraction": final.get("mean_mask_fraction", ""),
                "warning": final.get("warning", ""),
                "pred_R": row.get("pred_R", ""),
                "pred_D": row.get("pred_D", ""),
                "truth_R": row.get("truth_R", ""),
                "truth_D": row.get("truth_D", ""),
            })

    fields = [
        "set", "sample", "gene", "side", "missing_allele", "current_score", "stage_class",
        "pred_total_count", "truth_total_count", "em_rank", "em_frac", "em_weight", "raw_support_class",
        "full_support_pairs_k31", "binned_support_pairs_k31", "retained_fraction_k31",
        "copy_identifiability", "copy_fit_error", "mean_mask_fraction", "warning", "pred_R", "pred_D", "truth_R", "truth_D",
    ]
    write_tsv(args.out, fields, out_rows)

    stats = defaultdict(Counter)
    for miss in out_rows:
        for scope, name in (("ALL", "ALL"), ("gene", miss["gene"]), ("stage", miss["stage_class"]), ("gene_stage", f"{miss['gene']}|{miss['stage_class']}")):
            stats[(scope, name)]["missing_copies"] += 1

    summary_out = []
    summary_out.append({
        "scope": "score",
        "name": "baseline",
        "missing_copies": len(out_rows),
        "score": f"{score_from_summary(summary_rows)}/{len(summary_rows) * 4}",
    })
    for (scope, name), counts in sorted(stats.items()):
        summary_out.append({
            "scope": scope,
            "name": name,
            "missing_copies": counts["missing_copies"],
            "score": "",
        })
    write_tsv(args.summary_out, ["scope", "name", "missing_copies", "score"], summary_out)
    print(args.summary_out)


if __name__ == "__main__":
    main()
