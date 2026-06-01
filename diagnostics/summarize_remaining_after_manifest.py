#!/usr/bin/env python3
"""Summarize remaining missed truth copies after applying a proposal manifest.

Truth is used only for offline validation. The manifest itself is assumed to be
produced by a truth-free proposal generator.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_calls import load_g_group, normalize_for_display, overlap  # noqa: E402


def read_tsv(path: Path):
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


def side_quartet(quartet, side: str):
    return quartet[:2] if side == "R" else quartet[2:]


def current_quartet(row):
    return split_alleles(row["pred_R"]) + split_alleles(row["pred_D"])


def score_quartet(row, quartet, gmap) -> int:
    truth_r = normalize_for_display(split_alleles(row["truth_R"]), "2field", gmap)
    truth_d = normalize_for_display(split_alleles(row["truth_D"]), "2field", gmap)
    pred_r = sorted(quartet[:2])
    pred_d = sorted(quartet[2:])
    return overlap(truth_r, pred_r) + overlap(truth_d, pred_d)


def proposal_map(rows):
    return {(row["sample"], row["gene"]): split_alleles(row["proposed_quartet"]) for row in rows}


def copy_key(row):
    return (row["set"], row["sample"], row["gene"], row["side"], row["missing_allele"])


def summarize_remaining(stage_rows, proposals):
    by_copy = defaultdict(list)
    for row in stage_rows:
        by_copy[copy_key(row)].append(row)

    detail_rows = []
    fixed_stage = Counter()
    remaining_stage = Counter()
    fixed_gene_stage = Counter()
    remaining_gene_stage = Counter()

    for key, rows in sorted(by_copy.items()):
        set_label, sample, gene, side, allele = key
        first = rows[0]
        truth_side = split_alleles(first["truth_R"] if side == "R" else first["truth_D"])
        pred_side = split_alleles(first["pred_R"] if side == "R" else first["pred_D"])
        truth_count = truth_side.count(allele)
        baseline_pred_count = pred_side.count(allele)
        baseline_missing = max(0, truth_count - baseline_pred_count)
        proposed = proposals.get((sample, gene), current_quartet(first))
        proposed_count = side_quartet(proposed, side).count(allele)
        remaining = max(0, truth_count - proposed_count)
        fixed = max(0, baseline_missing - remaining)
        stage = first["stage_class"]
        raw_support = first.get("raw_support_class", "")
        if fixed:
            fixed_stage[stage] += fixed
            fixed_gene_stage[(gene, stage)] += fixed
        if remaining:
            remaining_stage[stage] += remaining
            remaining_gene_stage[(gene, stage)] += remaining
        detail_rows.append({
            "set": set_label,
            "sample": sample,
            "gene": gene,
            "side": side,
            "allele": allele,
            "stage_class": stage,
            "raw_support_class": raw_support,
            "baseline_missing": baseline_missing,
            "fixed_by_manifest": fixed,
            "remaining_missing": remaining,
            "proposed_quartet": ",".join(proposed),
        })

    summary_rows = []
    for stage in sorted(set(fixed_stage) | set(remaining_stage)):
        summary_rows.append({
            "scope": "stage",
            "name": stage,
            "fixed_copies": fixed_stage[stage],
            "remaining_copies": remaining_stage[stage],
        })
    for key in sorted(set(fixed_gene_stage) | set(remaining_gene_stage)):
        gene, stage = key
        summary_rows.append({
            "scope": "gene_stage",
            "name": f"{gene}|{stage}",
            "fixed_copies": fixed_gene_stage[key],
            "remaining_copies": remaining_gene_stage[key],
        })
    return detail_rows, summary_rows


def score_summary(quartet_rows, proposals, gmap):
    baseline = 0
    optimized = 0
    rows = 0
    regress = improve = neutral = 0
    for row in quartet_rows:
        key = (row["sample"], row["gene"])
        current = current_quartet(row)
        proposed = proposals.get(key, current)
        base_score = score_quartet(row, current, gmap)
        new_score = score_quartet(row, proposed, gmap)
        baseline += base_score
        optimized += new_score
        rows += 1
        if key in proposals:
            delta = new_score - base_score
            if delta > 0:
                improve += 1
            elif delta < 0:
                regress += 1
            else:
                neutral += 1
    return [{
        "scope": "score",
        "name": "manifest",
        "rows": rows,
        "improve": improve,
        "regress": regress,
        "neutral": neutral,
        "delta": optimized - baseline,
        "baseline": f"{baseline}/{rows * 4}",
        "optimized": f"{optimized}/{rows * 4}",
    }]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-attribution", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--quartet-summary", type=Path)
    parser.add_argument("--g-group", type=Path, default=SCRIPT_ROOT / "resources" / "spechla" / "db" / "HLA" / "hla_nom_g.txt")
    parser.add_argument("--out-detail", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    args = parser.parse_args()

    stage_rows = read_tsv(args.stage_attribution)
    proposals = proposal_map(read_tsv(args.manifest))
    detail_rows, summary_rows = summarize_remaining(stage_rows, proposals)

    if args.quartet_summary:
        gmap = load_g_group(args.g_group) if args.g_group.exists() else {}
        summary_rows = score_summary(read_tsv(args.quartet_summary), proposals, gmap) + summary_rows

    write_tsv(args.out_detail, [
        "set", "sample", "gene", "side", "allele", "stage_class", "raw_support_class",
        "baseline_missing", "fixed_by_manifest", "remaining_missing", "proposed_quartet",
    ], detail_rows)
    write_tsv(args.out_summary, ["scope", "name", "fixed_copies", "remaining_copies", "rows", "improve", "regress", "neutral", "delta", "baseline", "optimized"], summary_rows)


if __name__ == "__main__":
    main()
