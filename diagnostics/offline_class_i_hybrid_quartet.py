#!/usr/bin/env python3
"""Gate source-agnostic class-I quartet proposals with direct read-pair evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPT_ROOT))

from diagnostics.direct_read_quartet_likelihood import (
    DEFAULT_IMGT,
    build_informative_kmers,
    load_candidate_sequences,
    split_alleles,
)
from diagnostics.offline_joint_quartet_posterior import unique_mixture_groupings


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def iter_fastq_pairs(fq1: Path, fq2: Path):
    with open_text(fq1) as handle1, open_text(fq2) as handle2:
        while True:
            header1 = handle1.readline()
            header2 = handle2.readline()
            if not header1 or not header2:
                return
            sequence1 = handle1.readline().strip().upper()
            sequence2 = handle2.readline().strip().upper()
            handle1.readline()
            handle2.readline()
            quality1 = handle1.readline().strip()
            quality2 = handle2.readline().strip()
            if not quality1 or not quality2:
                raise ValueError("truncated paired FASTQ")
            yield sequence1, quality1, sequence2, quality2


def quality_weight(quality: str) -> float:
    if not quality:
        return 0.0
    minimum_q = min(max(0, ord(char) - 33) for char in quality)
    return 1.0 - 10.0 ** (-minimum_q / 10.0)


def mate_evidence(
    sequence: str,
    quality: str,
    kmer_owners: dict[str, tuple[str, ...]],
    k: int,
) -> tuple[Counter[str], Counter[str]]:
    evidence: Counter[str] = Counter()
    private: Counter[str] = Counter()
    seen = set()
    for start in range(len(sequence) - k + 1):
        kmer = sequence[start:start + k]
        if "N" in kmer or kmer in seen:
            continue
        seen.add(kmer)
        owners = kmer_owners.get(kmer)
        if not owners:
            continue
        weight = quality_weight(quality[start:start + k]) / len(owners)
        for allele in owners:
            evidence[allele] += weight
        if len(owners) == 1:
            private[owners[0]] += 1
    return evidence, private


def pair_evidence(
    sequence1: str,
    quality1: str,
    sequence2: str,
    quality2: str,
    kmer_owners: dict[str, tuple[str, ...]],
    k: int,
    concordance_bonus: float,
) -> tuple[dict[str, float], set[str]]:
    mate1, private1 = mate_evidence(sequence1, quality1, kmer_owners, k)
    mate2, private2 = mate_evidence(sequence2, quality2, kmer_owners, k)
    alleles = set(mate1) | set(mate2)
    combined = {
        allele: mate1[allele] + mate2[allele]
        + concordance_bonus * min(mate1[allele], mate2[allele])
        for allele in alleles
    }
    private = {allele for allele in set(private1) | set(private2) if private1[allele] or private2[allele]}
    return combined, private


def logsumexp(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def grouping_weights(
    grouping: tuple[tuple[str, str], tuple[str, str]], major_fraction: float
) -> dict[str, float]:
    weights: dict[str, float] = defaultdict(float)
    for allele in grouping[0]:
        weights[allele] += major_fraction / 2.0
    for allele in grouping[1]:
        weights[allele] += (1.0 - major_fraction) / 2.0
    return dict(weights)


def pair_log_likelihood(
    evidence: dict[str, float],
    weights: dict[str, float],
    score_scale: float,
) -> float:
    return logsumexp([
        math.log(max(weight, 1e-12)) + score_scale * evidence.get(allele, 0.0)
        for allele, weight in weights.items()
    ])


def quartet_log_evidence(
    pair_rows: list[dict[str, float]],
    quartet: tuple[str, ...],
    major_fractions: tuple[float, ...],
    score_scale: float,
) -> tuple[float, tuple[tuple[str, str], tuple[str, str]], float]:
    states = []
    for grouping in unique_mixture_groupings(tuple(sorted(quartet))):
        for major_fraction in major_fractions:
            weights = grouping_weights(grouping, major_fraction)
            score = sum(pair_log_likelihood(row, weights, score_scale) for row in pair_rows)
            states.append((score, grouping, major_fraction))
    marginal = logsumexp([state[0] for state in states]) - math.log(len(states))
    best = max(states, key=lambda state: state[0])
    return marginal, best[1], best[2]


def compare_quartets(
    pair_rows: list[dict[str, float]],
    baseline: tuple[str, ...],
    proposal: tuple[str, ...],
    score_scale: float = 0.35,
    major_fraction_min: float = 0.50,
    major_fraction_max: float = 0.95,
    major_fraction_step: float = 0.025,
) -> dict[str, object]:
    fractions = tuple(
        round(major_fraction_min + index * major_fraction_step, 6)
        for index in range(
            int(round((major_fraction_max - major_fraction_min) / major_fraction_step)) + 1
        )
    )
    baseline_score, baseline_grouping, baseline_fraction = quartet_log_evidence(
        pair_rows, baseline, fractions, score_scale
    )
    proposal_score, proposal_grouping, proposal_fraction = quartet_log_evidence(
        pair_rows, proposal, fractions, score_scale
    )
    discriminating_pairs = 0
    proposal_pair_wins = 0
    for row in pair_rows:
        baseline_pair = max(
            pair_log_likelihood(row, grouping_weights(grouping, fraction), score_scale)
            for grouping in unique_mixture_groupings(tuple(sorted(baseline)))
            for fraction in fractions
        )
        proposal_pair = max(
            pair_log_likelihood(row, grouping_weights(grouping, fraction), score_scale)
            for grouping in unique_mixture_groupings(tuple(sorted(proposal)))
            for fraction in fractions
        )
        if abs(proposal_pair - baseline_pair) >= 0.05:
            discriminating_pairs += 1
            proposal_pair_wins += proposal_pair > baseline_pair
    return {
        "baseline_log_evidence": baseline_score,
        "proposal_log_evidence": proposal_score,
        "log_bayes_factor": proposal_score - baseline_score,
        "discriminating_pairs": discriminating_pairs,
        "proposal_pair_wins": proposal_pair_wins,
        "baseline_major_fraction": baseline_fraction,
        "proposal_major_fraction": proposal_fraction,
        "baseline_grouping": baseline_grouping,
        "proposal_grouping": proposal_grouping,
    }


def fastq_paths(fastq_bench_root: Path, row: dict[str, str]) -> tuple[Path, Path]:
    base = fastq_bench_root / "reads" / row["experiment"] / row["condition"] / row["sample"]
    return Path(f"{base}.R1.fastq.gz"), Path(f"{base}.R2.fastq.gz")


def proposal_is_supported(
    comparison: dict[str, object],
    baseline_private: int,
    proposal_private: int,
    args: argparse.Namespace,
) -> bool:
    informative_pairs = max(1, int(comparison["informative_pairs"]))
    normalized_log_bayes_factor = float(comparison["log_bayes_factor"]) / informative_pairs
    private_ratio = (proposal_private + 1.0) / (baseline_private + 1.0)
    return (
        float(comparison["log_bayes_factor"]) >= args.min_log_bayes_factor
        and normalized_log_bayes_factor >= args.min_log_bayes_factor_per_informative_pair
        and comparison["discriminating_pairs"] >= args.min_discriminating_pairs
        and proposal_private >= args.min_proposal_private_pairs
        and private_ratio >= args.min_private_pair_ratio
        and proposal_private + args.private_pair_slack >= baseline_private
    )


def evaluate_row(row: dict[str, str], args: argparse.Namespace) -> dict[str, object]:
    baseline = tuple(sorted(split_alleles(row["baseline_quartet"])))
    proposal = tuple(sorted(split_alleles(row["joint_quartet"])))
    output: dict[str, object] = dict(row)
    if baseline == proposal:
        output.update({
            "hybrid_decision": "baseline_same",
            "hybrid_quartet": ",".join(baseline),
            "log_bayes_factor": "0.000000",
            "informative_pairs": 0,
            "discriminating_pairs": 0,
            "proposal_pair_wins": 0,
            "private_pairs_baseline_only": 0,
            "private_pairs_proposal_only": 0,
            "hybrid_correct": row.get("baseline_correct", ""),
            "hybrid_delta": 0,
            "error": "",
        })
        return output
    candidates = sorted(set(baseline) | set(proposal))
    sequences = load_candidate_sequences(args.imgt, row["gene"], candidates)
    missing = [allele for allele, values in sequences.items() if not values]
    fq1, fq2 = fastq_paths(args.fastq_bench_root, row)
    if missing or not fq1.exists() or not fq2.exists():
        output.update({"hybrid_decision": "error", "error": ",".join(missing) or f"{fq1},{fq2}"})
        return output
    kmer_owners = build_informative_kmers(
        sequences, args.k, args.max_full_alleles_per_2field, args.max_owner_fraction
    )
    pair_rows = []
    private_pair_counts: Counter[str] = Counter()
    total_pairs = 0
    for sequence1, quality1, sequence2, quality2 in iter_fastq_pairs(fq1, fq2):
        total_pairs += 1
        evidence, private = pair_evidence(
            sequence1, quality1, sequence2, quality2, kmer_owners, args.k,
            args.concordance_bonus,
        )
        if sum(evidence.values()) < args.min_pair_evidence:
            continue
        pair_rows.append(evidence)
        for allele in private:
            private_pair_counts[allele] += 1
    comparison = compare_quartets(pair_rows, baseline, proposal, args.score_scale)
    comparison["informative_pairs"] = len(pair_rows)
    baseline_only = set(baseline) - set(proposal)
    proposal_only = set(proposal) - set(baseline)
    baseline_private = sum(private_pair_counts[allele] for allele in baseline_only)
    proposal_private = sum(private_pair_counts[allele] for allele in proposal_only)
    accept = proposal_is_supported(comparison, baseline_private, proposal_private, args)
    selected = proposal if accept else baseline
    baseline_correct = int(row["baseline_correct"]) if row.get("baseline_correct", "") else None
    proposal_correct = int(row["joint_correct"]) if row.get("joint_correct", "") else None
    selected_correct = proposal_correct if accept else baseline_correct
    output.update({
        "hybrid_decision": "proposal" if accept else "baseline",
        "hybrid_quartet": ",".join(selected),
        "log_bayes_factor": f"{float(comparison['log_bayes_factor']):.6f}",
        "log_bayes_factor_per_informative_pair": (
            f"{float(comparison['log_bayes_factor']) / max(1, len(pair_rows)):.6f}"
        ),
        "total_pairs": total_pairs,
        "informative_pairs": len(pair_rows),
        "discriminating_pairs": comparison["discriminating_pairs"],
        "proposal_pair_wins": comparison["proposal_pair_wins"],
        "private_pairs_baseline_only": baseline_private,
        "private_pairs_proposal_only": proposal_private,
        "private_pair_ratio": f"{(proposal_private + 1.0) / (baseline_private + 1.0):.6f}",
        "hybrid_correct": selected_correct if selected_correct is not None else "",
        "hybrid_delta": selected_correct - baseline_correct if selected_correct is not None else "",
        "error": "",
    })
    return output


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    evaluated = [row for row in rows if not row.get("error")]
    baseline = sum(int(row["baseline_correct"]) for row in evaluated)
    hybrid = sum(int(row["hybrid_correct"]) for row in evaluated)
    return [{
        "rows": len(rows),
        "evaluated_rows": len(evaluated),
        "proposal_accepted": sum(row["hybrid_decision"] == "proposal" for row in evaluated),
        "baseline_correct": baseline,
        "hybrid_correct": hybrid,
        "delta_correct": hybrid - baseline,
        "improved_loci": sum(int(row["hybrid_delta"]) > 0 for row in evaluated),
        "regressed_loci": sum(int(row["hybrid_delta"]) < 0 for row in evaluated),
        "errors": len(rows) - len(evaluated),
    }]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposal-tsv", required=True, type=Path)
    parser.add_argument("--fastq-bench-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--genes", nargs="+", default=["HLA-A", "HLA-B", "HLA-C"])
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-full-alleles-per-2field", type=int, default=25)
    parser.add_argument("--max-owner-fraction", type=float, default=0.75)
    parser.add_argument("--score-scale", type=float, default=0.35)
    parser.add_argument("--concordance-bonus", type=float, default=0.5)
    parser.add_argument("--min-pair-evidence", type=float, default=1.0)
    parser.add_argument("--min-log-bayes-factor", type=float, default=5.0)
    parser.add_argument("--min-log-bayes-factor-per-informative-pair", type=float, default=-math.inf)
    parser.add_argument("--min-discriminating-pairs", type=int, default=3)
    parser.add_argument("--min-proposal-private-pairs", type=int, default=0)
    parser.add_argument("--min-private-pair-ratio", type=float, default=0.0)
    parser.add_argument("--private-pair-slack", type=int, default=0)
    parser.add_argument("--conditions", nargs="+")
    parser.add_argument("--samples", nargs="+")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = [
        evaluate_row(row, args)
        for row in read_tsv(args.proposal_tsv)
        if row["gene"] in args.genes
        and (args.conditions is None or row["condition"] in args.conditions)
        and (args.samples is None or row["sample"] in args.samples)
    ]
    if not rows:
        raise RuntimeError("no proposal rows matched")
    write_tsv(args.out, rows)
    summary_path = args.summary or args.out.with_suffix(".summary.tsv")
    write_tsv(summary_path, summarize(rows))
    print(f"[class1-hybrid] rows={len(rows)} out={args.out}")
    print(f"[class1-hybrid] summary={summary_path}")


if __name__ == "__main__":
    main()