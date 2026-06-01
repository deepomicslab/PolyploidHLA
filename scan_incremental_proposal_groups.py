#!/usr/bin/env python3
"""Scan proposal groups for incremental validation gain.

Proposal files are assumed to be generated without truth. Truth from the summary
is used only here for offline validation.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from score_proposal_manifests import current_quartet, proposed_quartet, score_quartet, split_alleles


def read_tsv(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_summary(path: Path):
    rows = read_tsv(path)
    by_key = {(row["sample"], row["gene"]): row for row in rows}
    return rows, by_key


def load_base_proposals(paths):
    proposals = {}
    for path in paths:
        for row in read_tsv(path):
            proposed = proposed_quartet(row)
            if not proposed:
                continue
            proposals.setdefault((row["sample"], row["gene"]), split_alleles(proposed))
    return proposals


def row_score(summary_row, quartet):
    return score_quartet(
        quartet,
        split_alleles(summary_row["truth_R"]),
        split_alleles(summary_row["truth_D"]),
    )


def load_candidate_groups(paths, summary_by_key, base_keys):
    groups = defaultdict(list)
    for path in paths:
        for row in read_tsv(path):
            proposed = proposed_quartet(row)
            if not proposed:
                continue
            key = (row["sample"], row["gene"])
            if key not in summary_by_key or key in base_keys:
                continue
            group_key = (str(path), row.get("strategy") or row.get("rule") or path.stem, row["gene"])
            groups[group_key].append((key, row, split_alleles(proposed)))
    return groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--base-manifest", action="append", default=[], type=Path)
    parser.add_argument("--candidate", action="append", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    summary_rows, summary_by_key = load_summary(args.summary)
    base_proposals = load_base_proposals(args.base_manifest)
    base_scores = {}
    for row in summary_rows:
        key = (row["sample"], row["gene"])
        quartet = base_proposals.get(key, current_quartet(row))
        base_scores[key] = row_score(row, quartet)

    groups = load_candidate_groups(args.candidate, summary_by_key, set(base_proposals))
    out_rows = []
    for (source, rule, gene), rows in groups.items():
        stats = Counter()
        sample_rows = []
        seen = set()
        for key, row, proposed in rows:
            if key in seen:
                continue
            seen.add(key)
            summary_row = summary_by_key[key]
            base_score = base_scores[key]
            proposed_score = row_score(summary_row, proposed)
            delta = proposed_score - base_score
            verdict = "improve" if delta > 0 else "regress" if delta < 0 else "neutral"
            stats[verdict] += 1
            stats["delta"] += delta
            sample_rows.append(f"{key[0]}:{delta}")
        out_rows.append({
            "source": source,
            "rule": rule,
            "gene": gene,
            "rows": len(seen),
            "improve": stats["improve"],
            "regress": stats["regress"],
            "neutral": stats["neutral"],
            "delta": stats["delta"],
            "samples": ";".join(sample_rows),
        })

    out_rows.sort(key=lambda row: (-int(row["delta"]), int(row["regress"]), -int(row["improve"]), row["gene"], row["rule"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        fields = ["source", "rule", "gene", "rows", "improve", "regress", "neutral", "delta", "samples"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(out_rows)

    for row in out_rows[:30]:
        print("\t".join(str(row[field]) for field in ["source", "rule", "gene", "rows", "improve", "regress", "neutral", "delta", "samples"]))


if __name__ == "__main__":
    main()