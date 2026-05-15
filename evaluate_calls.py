#!/usr/bin/env python3
"""Evaluate final_calls.tsv against a small truth_typing.tsv table.

Truth format matches this project:
    HLA Typing  A A B B C C DRB1 DRB1 DQB1 DQB1 DPB1 DPB1
    DONOR       ...
    PATIENT     ...

Outputs 2-field and G group overlap counts for PATIENT(R) and DONOR(D) calls.
"""
from __future__ import annotations

import argparse
import csv
import os
from collections import Counter, defaultdict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BUNDLED_SPECHLA = SCRIPT_DIR / "resources" / "spechla"
DEFAULT_LEGACY_SPECHLA = SCRIPT_DIR.parent / "SpecHLA"
DEFAULT_SPECHLA = Path(os.environ.get("SPECHLA", DEFAULT_BUNDLED_SPECHLA if DEFAULT_BUNDLED_SPECHLA.exists() else DEFAULT_LEGACY_SPECHLA))
DEFAULT_G_GROUP = DEFAULT_SPECHLA / "db" / "HLA" / "hla_nom_g.txt"


def load_g_group(path: Path):
    gmap = {}
    if not path.exists():
        return gmap
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(";")
            if len(parts) < 3:
                continue
            gene = parts[0].rstrip("*")
            members = parts[-2].split("/") if parts[-2] else []
            group = parts[-1] or parts[-2]
            if not group:
                continue
            group_name = f"{gene}*{group}"
            for member in members:
                gmap[f"{gene}*{member}"] = group_name
            gmap[group_name] = group_name
    return gmap


def strip_expr_suffix(field: str) -> str:
    return field[:-1] if field and field[-1].isalpha() and field[-1] != "G" else field


def clean_allele(allele: str) -> str:
    if not allele or allele == "NA":
        return "NA"
    allele = allele.replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    fields = rest.split(":")
    fields[-1] = strip_expr_suffix(fields[-1])
    return f"{gene}*{':'.join(fields)}"


def normalize_truth_gene_allele(gene_short: str, allele: str):
    gene_short = gene_short.strip()
    allele = allele.strip()
    if gene_short in {"DRB345", "DPB345"}:
        if not allele or allele in {"NA", "-", "."}:
            return "HLA-DRB345", "NA"
        allele = allele.replace("HLA-", "")
        if allele.startswith(("DRB3*", "DRB4*", "DRB5*")):
            return "HLA-DRB345", allele
        if "*" in allele:
            return "HLA-DRB345", allele
        fields = [field for field in allele.rstrip(":").split(":") if field]
        if not fields:
            return "HLA-DRB345", "NA"
        locus_code = fields[0]
        locus = {"03": "DRB3", "3": "DRB3", "04": "DRB4", "4": "DRB4", "05": "DRB5", "5": "DRB5"}.get(locus_code)
        if locus and len(fields) >= 2:
            return "HLA-DRB345", f"{locus}*{':'.join(fields[1:])}"
        return "HLA-DRB345", f"DRB345*{allele}"
    gene = f"HLA-{gene_short}"
    if "*" in allele:
        return gene, allele.replace("HLA-", "")
    return gene, f"{gene_short}*{allele}"


def norm_allele(allele: str, level: str, gmap=None) -> str:
    allele = clean_allele(allele)
    if allele == "NA" or "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    if level == "g_group":
        if allele.endswith("G"):
            return allele
        candidates = [allele]
        fields = rest.split(":")
        if len(fields) == 2:
            candidates.extend([f"{gene}*{rest}:01", f"{gene}*{rest}:01:01"])
        elif len(fields) == 3:
            candidates.append(f"{gene}*{rest}:01")
        for candidate in candidates:
            if gmap and candidate in gmap:
                return gmap[candidate]
        return allele
    rest = rest.replace("G", "")
    fields = rest.split(":")
    if level == "2field":
        return f"{gene}*{':'.join(fields[:2])}"
    raise ValueError(level)


def field_count(allele: str) -> int:
    allele = clean_allele(allele)
    if "*" not in allele:
        return 0
    return len(allele.split("*", 1)[1].replace("G", "").split(":"))


def g_groups_matching_prefix(allele: str, gmap):
    allele = clean_allele(allele)
    if not gmap or "*" not in allele:
        return set()
    prefix = allele + ":"
    return {group for member, group in gmap.items() if not member.endswith("G") and member.startswith(prefix)}


def g_group_truth_target(allele: str, gmap):
    allele = clean_allele(allele)
    if field_count(allele) == 2:
        groups = g_groups_matching_prefix(allele, gmap)
        if len(groups) != 1:
            return "2field", norm_allele(allele, "2field", gmap)
    as_group = norm_allele(allele, "g_group", gmap)
    if allele.endswith("G") or as_group != allele or field_count(allele) >= 3:
        return "g_group", as_group
    return "2field", norm_allele(allele, "2field", gmap)


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
            gene, normalized = normalize_truth_gene_allele(gene_short, allele)
            truth[side][gene].append(normalized)
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


def overlap_g_group_truth_resolution(truth_vals, pred_vals, gmap):
    truth_targets = [g_group_truth_target(x, gmap) for x in truth_vals]
    pred_by_level = {
        "2field": Counter(norm_allele(x, "2field", gmap) for x in pred_vals),
        "g_group": Counter(norm_allele(x, "g_group", gmap) for x in pred_vals),
    }
    hits = 0
    for target_level, target_value in truth_targets:
        if pred_by_level[target_level][target_value] > 0:
            hits += 1
            pred_by_level[target_level][target_value] -= 1
    return hits


def normalize_for_display(vals, level, gmap):
    if level == "g_group":
        return sorted(value for _, value in (g_group_truth_target(x, gmap) for x in vals))
    return sorted(norm_allele(x, level, gmap) for x in vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--call-set", choices=("report", "full", "2field", "g_group"), default="report",
                    help="which final_calls columns to evaluate")
    ap.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP,
                    help="WMDA hla_nom_g.txt used for G group conversion")
    args = ap.parse_args()

    truth = load_truth(args.truth)
    calls = load_calls(args.calls, args.call_set)
    gmap = load_g_group(args.g_group)
    print(f"call_set\t{args.call_set}")
    print(f"g_group_file\t{args.g_group}")
    for level in ("2field", "g_group"):
        total = 0
        ok = 0
        mismatches = []
        for side in ("PATIENT", "DONOR"):
            for gene, truth_vals in truth[side].items():
                pred_vals = calls[side].get(gene, [])
                if not pred_vals:
                    continue
                if level == "g_group":
                    t = normalize_for_display(truth_vals, level, gmap)
                    p = sorted(norm_allele(x, level, gmap) for x in pred_vals)
                    h = overlap_g_group_truth_resolution(truth_vals, pred_vals, gmap)
                else:
                    t = normalize_for_display(truth_vals, level, gmap)
                    p = sorted(norm_allele(x, level, gmap) for x in pred_vals)
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
