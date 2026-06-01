#!/usr/bin/env python3
"""Combine proposal manifests and score them against a quartet summary.

Proposal generation is assumed to be truth-free. This script uses truth only for
post-hoc validation of the emitted manifests.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path):
    if not path.exists():
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


def score_quartet(quartet, truth_r, truth_d) -> int:
    def overlap(truth_vals, pred_vals):
        counts = Counter(truth_vals)
        hits = 0
        for allele in pred_vals:
            if counts[allele] > 0:
                counts[allele] -= 1
                hits += 1
        return hits
    return overlap(truth_r, quartet[:2]) + overlap(truth_d, quartet[2:4])


def current_quartet(row):
    return split_alleles(row["pred_R"]) + split_alleles(row["pred_D"])


def proposed_quartet(row):
    for field in ("proposed_quartet", "new_2field_quartet", "direct_quartet", "proposed"):
        value = row.get(field, "")
        if value:
            return value
    direct_r = row.get("direct_R", "")
    direct_d = row.get("direct_D", "")
    if direct_r or direct_d:
        return ",".join(split_alleles(direct_r) + split_alleles(direct_d))
    return ""


def load_manifest(spec: str):
    if "=" in spec:
        label, path_text = spec.split("=", 1)
    else:
        path_text = spec
        label = Path(path_text).stem
    path = Path(path_text)
    rows = []
    for row in read_tsv(path):
        proposed = proposed_quartet(row)
        if not proposed:
            continue
        rows.append({
            "sample": row["sample"],
            "set": row.get("set", ""),
            "gene": row["gene"],
            "rule": row.get("rule") or row.get("gate") or row.get("strategy") or label,
            "source_label": label,
            "source": str(path),
            "proposed_quartet": proposed,
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--manifest", action="append", required=True,
                        help="Input manifest in priority order. Use label=path or path.")
    parser.add_argument("--out-manifest", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--out-conflicts", required=True, type=Path)
    args = parser.parse_args()

    summary_rows = read_tsv(args.summary)
    proposals = {}
    conflicts = []
    for spec in args.manifest:
        for row in load_manifest(spec):
            key = (row["sample"], row["gene"])
            if key in proposals:
                conflicts.append({
                    "sample": row["sample"],
                    "gene": row["gene"],
                    "kept_rule": proposals[key]["rule"],
                    "skipped_rule": row["rule"],
                    "kept_source": proposals[key]["source_label"],
                    "skipped_source": row["source_label"],
                    "kept_quartet": proposals[key]["proposed_quartet"],
                    "skipped_quartet": row["proposed_quartet"],
                })
                continue
            proposals[key] = row

    manifest_rows = []
    by_rule = defaultdict(Counter)
    by_gene = defaultdict(Counter)
    baseline_total = 0
    optimized_total = 0
    for row in summary_rows:
        key = (row["sample"], row["gene"])
        truth_r = split_alleles(row["truth_R"])
        truth_d = split_alleles(row["truth_D"])
        current = current_quartet(row)
        base_score = score_quartet(current, truth_r, truth_d)
        proposal = proposals.get(key)
        proposed = split_alleles(proposal["proposed_quartet"]) if proposal else current
        new_score = score_quartet(proposed, truth_r, truth_d)
        baseline_total += base_score
        optimized_total += new_score
        if not proposal:
            continue
        delta = new_score - base_score
        verdict = "improve" if delta > 0 else "regress" if delta < 0 else "neutral"
        for bucket, name in ((by_rule, proposal["rule"]), (by_gene, row["gene"])):
            bucket[name]["rows"] += 1
            bucket[name]["delta"] += delta
            bucket[name][verdict] += 1
        manifest_rows.append({
            "sample": row["sample"],
            "set": row.get("set", proposal.get("set", "")),
            "gene": row["gene"],
            "rule": proposal["rule"],
            "source_label": proposal["source_label"],
            "source": proposal["source"],
            "current_score": base_score,
            "proposed_score": new_score,
            "delta": delta,
            "verdict": verdict,
            "current_quartet": ",".join(current),
            "proposed_quartet": ",".join(proposed),
            "truth_R": row["truth_R"],
            "truth_D": row["truth_D"],
        })

    summary_out = [{
        "scope": "ALL",
        "name": "combined",
        "rows": len(manifest_rows),
        "improve": sum(1 for row in manifest_rows if row["verdict"] == "improve"),
        "regress": sum(1 for row in manifest_rows if row["verdict"] == "regress"),
        "neutral": sum(1 for row in manifest_rows if row["verdict"] == "neutral"),
        "delta": optimized_total - baseline_total,
        "baseline": f"{baseline_total}/{len(summary_rows) * 4}",
        "optimized": f"{optimized_total}/{len(summary_rows) * 4}",
    }]
    for scope, bucket in (("rule", by_rule), ("gene", by_gene)):
        for name, stats in sorted(bucket.items()):
            summary_out.append({
                "scope": scope,
                "name": name,
                "rows": stats["rows"],
                "improve": stats["improve"],
                "regress": stats["regress"],
                "neutral": stats["neutral"],
                "delta": stats["delta"],
                "baseline": "",
                "optimized": "",
            })

    write_tsv(args.out_manifest, [
        "sample", "set", "gene", "rule", "source_label", "source",
        "current_score", "proposed_score", "delta", "verdict",
        "current_quartet", "proposed_quartet", "truth_R", "truth_D",
    ], manifest_rows)
    write_tsv(args.out_summary, [
        "scope", "name", "rows", "improve", "regress", "neutral", "delta", "baseline", "optimized",
    ], summary_out)
    write_tsv(args.out_conflicts, [
        "sample", "gene", "kept_rule", "skipped_rule", "kept_source", "skipped_source",
        "kept_quartet", "skipped_quartet",
    ], conflicts)

    print(f"baseline\t{baseline_total}/{len(summary_rows) * 4}")
    print(f"optimized\t{optimized_total}/{len(summary_rows) * 4}")
    print(f"proposals\t{len(manifest_rows)}")
    print(f"conflicts\t{len(conflicts)}")


if __name__ == "__main__":
    main()
