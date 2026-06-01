#!/usr/bin/env python3
"""Combine GenDx noA4 rescue proposal manifests and score them post hoc.

The proposal sources are truth-free diagnostics. Truth columns from the quartet
summary are used only here to measure overlap and conflicts.
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

from hla_ld_maps import load_drb1_dqb1_map  # noqa: E402

DEFAULT_DRB1_DQB1_LD_MAP = SCRIPT_ROOT / "resources" / "drb1_dqb1_ld.tsv"


def read_rows(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, fields, rows) -> None:
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


def load_summary(path: Path):
    rows = read_rows(path)
    by_key = {(row["sample"], row["gene"]): row for row in rows}
    return rows, by_key


def current_quartet(row):
    return split_alleles(row["pred_R"]) + split_alleles(row["pred_D"])


def truth_pair(row):
    return split_alleles(row["truth_R"]), split_alleles(row["truth_D"])


def apply_drb1_ld_side_fix(proposed, current, dqb1_current, drb1_to_dqb1):
    fixed = list(proposed)
    current_counts = Counter(current)
    proposed_counts = Counter(proposed)
    introduced = []
    for allele, count in proposed_counts.items():
        extra = count - current_counts.get(allele, 0)
        introduced.extend([allele] * max(0, extra))
    for allele in introduced:
        partner = drb1_to_dqb1.get(allele)
        if not partner:
            continue
        r_has = partner in dqb1_current[:2]
        d_has = partner in dqb1_current[2:4]
        if r_has == d_has:
            continue
        target_indices = [0, 1] if r_has else [2, 3]
        if any(fixed[index] == allele for index in target_indices):
            continue
        source_indices = [index for index, value in enumerate(fixed) if value == allele]
        if not source_indices:
            continue
        source_index = source_indices[0]
        target_index = target_indices[0]
        fixed[source_index], fixed[target_index] = fixed[target_index], fixed[source_index]
    return fixed


def add_proposal(proposals, conflicts, source_row):
    key = (source_row["sample"], source_row["gene"])
    if key in proposals:
        conflicts.append({
            "sample": source_row["sample"],
            "gene": source_row["gene"],
            "kept_rule": proposals[key]["rule"],
            "skipped_rule": source_row["rule"],
            "kept_quartet": proposals[key]["proposed_quartet"],
            "skipped_quartet": source_row["proposed_quartet"],
        })
        return
    proposals[key] = source_row


def add_rebalance(proposals, conflicts, path: Path, strategy: str) -> None:
    for row in read_rows(path):
        if row.get("strategy") != strategy or not row.get("proposed"):
            continue
        add_proposal(proposals, conflicts, {
            "sample": row["sample"],
            "set": row["set"],
            "gene": row["gene"],
            "rule": strategy,
            "source": str(path),
            "proposed_quartet": row["proposed"],
        })


def add_drb1_match_ld(proposals, conflicts, path: Path, summary_by_key, strategy: str, drb1_to_dqb1) -> None:
    for row in read_rows(path):
        if row.get("strategy") != strategy or row.get("gene") != "HLA-DRB1" or not row.get("proposed"):
            continue
        summary_row = summary_by_key[(row["sample"], row["gene"])]
        dqb1_row = summary_by_key.get((row["sample"], "HLA-DQB1"))
        proposed = split_alleles(row["proposed"])
        rule = "DRB1_match_guard"
        if dqb1_row:
            fixed = apply_drb1_ld_side_fix(proposed, current_quartet(summary_row), current_quartet(dqb1_row), drb1_to_dqb1)
            if fixed != proposed:
                proposed = fixed
                rule = "DRB1_match_guard_ld_side_fix"
        add_proposal(proposals, conflicts, {
            "sample": row["sample"],
            "set": row["set"],
            "gene": row["gene"],
            "rule": rule,
            "source": str(path),
            "proposed_quartet": ",".join(proposed),
        })


def add_manifest(proposals, conflicts, path: Path, source_rule: str) -> None:
    for row in read_rows(path):
        proposed = row.get("new_2field_quartet") or row.get("direct_quartet")
        if not proposed:
            continue
        add_proposal(proposals, conflicts, {
            "sample": row["sample"],
            "set": row.get("set", ""),
            "gene": row["gene"],
            "rule": row.get("gate") or row.get("rule") or source_rule,
            "source": str(path),
            "proposed_quartet": proposed,
        })


def add_proposed_quartet_manifest(proposals, conflicts, path: Path, source_rule: str) -> None:
    for row in read_rows(path):
        proposed = row.get("proposed_quartet") or row.get("new_2field_quartet")
        if not proposed:
            continue
        add_proposal(proposals, conflicts, {
            "sample": row["sample"],
            "set": row.get("set", ""),
            "gene": row["gene"],
            "rule": row.get("rule") or source_rule,
            "source": str(path),
            "proposed_quartet": proposed,
        })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--rebal-tsv", required=True, type=Path)
    parser.add_argument("--match-tsv", required=True, type=Path)
    parser.add_argument("--dpb1-manifest", required=True, type=Path)
    parser.add_argument("--direct-manifest", required=True, type=Path)
    parser.add_argument("--dqb1-anchor-manifest", type=Path, default=None)
    parser.add_argument("--drb1-dqb1-ld-map", type=Path, default=DEFAULT_DRB1_DQB1_LD_MAP)
    parser.add_argument("--out-manifest", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--out-conflicts", required=True, type=Path)
    args = parser.parse_args()

    summary_rows, summary_by_key = load_summary(args.summary)
    drb1_to_dqb1, _dqb1_to_drb1 = load_drb1_dqb1_map(args.drb1_dqb1_ld_map)
    proposals = {}
    conflicts = []
    add_rebalance(proposals, conflicts, args.rebal_tsv, "C0_one_EM_AF_PHASE_REBAL")
    add_drb1_match_ld(proposals, conflicts, args.match_tsv, summary_by_key, "C0_one_EM_AF_PHASE_MATCH_NODPB1", drb1_to_dqb1)
    add_manifest(proposals, conflicts, args.dpb1_manifest, "DPB1_guarded")
    add_manifest(proposals, conflicts, args.direct_manifest, "direct_guarded")
    if args.dqb1_anchor_manifest:
        add_proposed_quartet_manifest(proposals, conflicts, args.dqb1_anchor_manifest, "DQB1_DRB1_anchor_guarded")

    manifest_rows = []
    baseline_total = 0
    combined_total = 0
    by_rule = defaultdict(lambda: {"rows": 0, "delta": 0, "improve": 0, "regress": 0, "neutral": 0})
    by_gene = defaultdict(lambda: {"rows": 0, "delta": 0, "improve": 0, "regress": 0, "neutral": 0})
    for row in summary_rows:
        key = (row["sample"], row["gene"])
        truth_r, truth_d = truth_pair(row)
        current = current_quartet(row)
        current_score = score_quartet(current, truth_r, truth_d)
        baseline_total += current_score
        proposal = proposals.get(key)
        proposed = split_alleles(proposal["proposed_quartet"]) if proposal else current
        proposed_score = score_quartet(proposed, truth_r, truth_d)
        combined_total += proposed_score
        if not proposal:
            continue
        delta = proposed_score - current_score
        verdict = "improve" if delta > 0 else "regress" if delta < 0 else "neutral"
        by_rule[proposal["rule"]]["rows"] += 1
        by_rule[proposal["rule"]]["delta"] += delta
        by_rule[proposal["rule"]][verdict] += 1
        by_gene[row["gene"]]["rows"] += 1
        by_gene[row["gene"]]["delta"] += delta
        by_gene[row["gene"]][verdict] += 1
        manifest_rows.append({
            "sample": row["sample"],
            "set": row["set"],
            "gene": row["gene"],
            "rule": proposal["rule"],
            "source": proposal["source"],
            "current_score": current_score,
            "proposed_score": proposed_score,
            "delta": delta,
            "verdict": verdict,
            "current_quartet": ",".join(current),
            "proposed_quartet": ",".join(proposed),
            "truth_R": row["truth_R"],
            "truth_D": row["truth_D"],
        })

    manifest_fields = [
        "sample", "set", "gene", "rule", "source", "current_score",
        "proposed_score", "delta", "verdict", "current_quartet",
        "proposed_quartet", "truth_R", "truth_D",
    ]
    write_rows(args.out_manifest, manifest_fields, manifest_rows)

    summary_rows_out = [{
        "scope": "ALL",
        "name": "combined",
        "rows": len(manifest_rows),
        "improve": sum(1 for row in manifest_rows if row["verdict"] == "improve"),
        "regress": sum(1 for row in manifest_rows if row["verdict"] == "regress"),
        "neutral": sum(1 for row in manifest_rows if row["verdict"] == "neutral"),
        "delta": combined_total - baseline_total,
        "baseline": f"{baseline_total}/{len(summary_rows) * 4}",
        "combined": f"{combined_total}/{len(summary_rows) * 4}",
    }]
    for scope, data in (("rule", by_rule), ("gene", by_gene)):
        for name, stats in sorted(data.items()):
            summary_rows_out.append({
                "scope": scope,
                "name": name,
                "rows": stats["rows"],
                "improve": stats["improve"],
                "regress": stats["regress"],
                "neutral": stats["neutral"],
                "delta": stats["delta"],
                "baseline": "",
                "combined": "",
            })
    write_rows(args.out_summary, ["scope", "name", "rows", "improve", "regress", "neutral", "delta", "baseline", "combined"], summary_rows_out)
    write_rows(args.out_conflicts, ["sample", "gene", "kept_rule", "skipped_rule", "kept_quartet", "skipped_quartet"], conflicts)

    print(f"baseline\t{baseline_total}/{len(summary_rows) * 4}")
    print(f"combined\t{combined_total}/{len(summary_rows) * 4}")
    print(f"proposals\t{len(manifest_rows)}")
    print(f"conflicts\t{len(conflicts)}")


if __name__ == "__main__":
    main()