#!/usr/bin/env python3
"""Evaluate final_calls.tsv against a small truth_typing.tsv table.

Truth format matches this project:
    HLA Typing  A A B B C C DRB1 DRB1 DQB1 DQB1 DPB1 DPB1
    DONOR       ...
    PATIENT     ...

Outputs per-resolution overlap counts for PATIENT(R) and DONOR(D) calls.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def norm_allele(allele: str, level: str) -> str:
    if not allele or allele == "NA":
        return "NA"
    allele = allele.replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    rest = rest.replace("G", "").replace("N", "")
    fields = rest.split(":")
    if level == "exact_noG":
        return f"{gene}*{':'.join(fields)}"
    if level == "3field":
        return f"{gene}*{':'.join(fields[:3])}"
    if level == "2field":
        return f"{gene}*{':'.join(fields[:2])}"
    raise ValueError(level)


def load_truth(path: Path):
    with path.open() as f:
        rows = [r for r in csv.reader(f, delimiter="\t") if r]
    if len(rows) < 3:
        raise SystemExit(f"truth file is empty or malformed: {path}")
    header = rows[0][1:]
    truth = defaultdict(lambda: defaultdict(list))
    for row in rows[1:]:
        side = row[0]
        for gene_short, allele in zip(header, row[1:]):
            gene = f"HLA-{gene_short}"
            truth[side][gene].append(f"{gene_short}*{allele}")
    return truth


def load_calls(path: Path, call_set: str):
    out = defaultdict(dict)
    with path.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            gene = row["gene"]
            suffix = call_set
            if suffix == "report" and "R1_report" not in row:
                suffix = "full"
            out["PATIENT"][gene] = [row.get(f"R1_{suffix}", "NA"), row.get(f"R2_{suffix}", "NA")]
            out["DONOR"][gene] = [row.get(f"D1_{suffix}", "NA"), row.get(f"D2_{suffix}", "NA")]
    return out


def overlap(truth_vals, pred_vals):
    truth_counter = Counter(truth_vals)
    hits = 0
    for pred in pred_vals:
        if truth_counter[pred] > 0:
            hits += 1
            truth_counter[pred] -= 1
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--call-set", choices=("report", "full", "2field"), default="report",
                    help="which final_calls columns to evaluate")
    args = ap.parse_args()

    truth = load_truth(args.truth)
    calls = load_calls(args.calls, args.call_set)
    print(f"call_set\t{args.call_set}")
    for level in ("exact_noG", "3field", "2field"):
        total = 0
        ok = 0
        mismatches = []
        for side in ("PATIENT", "DONOR"):
            for gene, truth_vals in truth[side].items():
                pred_vals = calls[side].get(gene, [])
                if not pred_vals:
                    continue
                t = sorted(norm_allele(x, level) for x in truth_vals)
                p = sorted(norm_allele(x, level) for x in pred_vals)
                h = overlap(t, p)
                ok += h
                total += len(t)
                if h != len(t):
                    mismatches.append((side, gene, h, len(t), ",".join(t), ",".join(p)))
        print(f"{level}\t{ok}/{total}")
        for m in mismatches:
            print("MISMATCH\t" + "\t".join(map(str, m)))


if __name__ == "__main__":
    main()
