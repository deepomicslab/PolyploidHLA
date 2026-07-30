#!/usr/bin/env python3
"""Replay a private-read-gated class-I rescue without modifying caller output."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from diagnostics.direct_read_quartet_likelihood import (
    DEFAULT_IMGT,
    clean_allele,
    iter_fastq_pairs,
    iter_kmers,
    load_imgt_cached,
    revcomp,
)

DEFAULT_GENES = ("HLA-A", "HLA-B", "HLA-C")


def normalize_2field(allele: str) -> str:
    allele = allele.strip().replace("HLA-", "").rstrip("GP")
    if "*" not in allele:
        return allele
    gene, fields = allele.split("*", 1)
    parts = fields.split(":")
    return f"{gene}*{parts[0]}:{parts[1]}" if len(parts) >= 2 else allele


def multiset_hits(truth: list[str], prediction: list[str]) -> int:
    return sum((Counter(truth) & Counter(prediction)).values())


def load_truth(bench_root: Path, experiment_glob: str, genes=DEFAULT_GENES):
    gene_set = set(genes)
    truth = defaultdict(list)
    for path in sorted((bench_root / "truth").glob(f"{experiment_glob}/copies.tsv")):
        experiment = path.parent.name
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["gene"] in gene_set:
                    key = (experiment, row["condition"], row["sample_id"], row["gene"])
                    truth[key].append(normalize_2field(row["allele_2field"]))
    if not truth:
        raise FileNotFoundError(f"no truth rows matched {experiment_glob!r}")
    return truth


def load_calls(bench_root: Path, experiment_glob: str, run_suffix: str, genes=DEFAULT_GENES):
    gene_set = set(genes)
    calls = defaultdict(list)
    pattern = bench_root / "runs" / f"{experiment_glob}{run_suffix}" / "*" / "SIM*" / "asm_v2" / "SIM*" / "*.copy_calls.tsv"
    for path_text in glob.glob(str(pattern)):
        path = Path(path_text)
        run_experiment = path.parents[4].name
        experiment = run_experiment.removesuffix(run_suffix) if run_suffix else run_experiment
        condition = path.parents[3].name
        sample = path.parent.name
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["gene"] in gene_set:
                    calls[experiment, condition, sample, row["gene"]].append(
                        normalize_2field(row["allele_2field"])
                    )
    return calls


def read_em_counts(path: Path) -> dict[str, float]:
    counts = {}
    if not path.exists():
        return counts
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[normalize_2field(row["allele_2field"])] = float(row["em_weight"])
    return counts


@lru_cache(maxsize=4)
def load_2field_sequences(imgt: str, gene: str) -> dict[str, tuple[str, ...]]:
    sequences = defaultdict(list)
    gene_prefix = gene.replace("HLA-", "") + "*"
    for full_name, sequence in load_imgt_cached(imgt).items():
        clean_name = clean_allele(full_name)
        if not clean_name:
            continue
        allele = normalize_2field(clean_name)
        if allele.startswith(gene_prefix):
            sequences[allele].append(sequence.upper())
    return {allele: tuple(values) for allele, values in sequences.items()}


@lru_cache(maxsize=4096)
def allele_kmers(
    imgt: str,
    gene: str,
    allele: str,
    k: int,
    max_full_alleles: int,
) -> frozenset[str]:
    kmers = set()
    for sequence in load_2field_sequences(imgt, gene).get(allele, ())[:max_full_alleles]:
        for source in (sequence, revcomp(sequence)):
            kmers.update(iter_kmers(source, k))
    return frozenset(kmers)


def build_cached_informative_kmers(
    imgt: Path,
    gene: str,
    candidates: list[str],
    k: int,
    max_full_alleles: int,
) -> dict[str, tuple[str, ...]]:
    kmer_owners = defaultdict(set)
    for allele in candidates:
        for kmer in allele_kmers(str(imgt), gene, allele, k, max_full_alleles):
            kmer_owners[kmer].add(allele)
    max_owners = max(1, math.floor(len(candidates) * 0.75))
    return {
        kmer: tuple(sorted(owners))
        for kmer, owners in kmer_owners.items()
        if 0 < len(owners) < len(candidates) and len(owners) <= max_owners
    }


def private_pair_support(
    fq1: Path,
    fq2: Path,
    gene: str,
    baseline: list[str],
    incoming_candidates: list[str],
    imgt: Path,
    k: int,
    max_full_alleles: int,
) -> dict[str, Counter[str]]:
    candidates = sorted(set(baseline + incoming_candidates))
    sequences = load_2field_sequences(str(imgt), gene)
    comparisons: dict[str, dict[str, tuple[str, ...]]] = {}
    for incoming in incoming_candidates:
        comparison = sorted(set(baseline + [incoming]))
        if any(not sequences[allele] for allele in comparison):
            continue
        comparisons[incoming] = build_cached_informative_kmers(
            imgt,
            gene,
            comparison,
            k,
            max_full_alleles,
        )
    reverse_index: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for incoming, owners_by_kmer in comparisons.items():
        for kmer, owners in owners_by_kmer.items():
            if len(owners) == 1:
                reverse_index[kmer].append((incoming, owners[0]))

    support = {incoming: Counter() for incoming in comparisons}
    for seq1, seq2 in iter_fastq_pairs(fq1, fq2):
        pair_kmers = set(iter_kmers(seq1, k)) | set(iter_kmers(seq2, k))
        seen = set()
        for kmer in pair_kmers:
            seen.update(reverse_index.get(kmer, ()))
        for incoming, allele in seen:
            support[incoming][allele] += 1
    return support


def support_cache_key(
    fq1: Path,
    fq2: Path,
    gene: str,
    baseline: list[str],
    incoming_candidates: list[str],
    imgt: Path,
    k: int,
    max_full_alleles: int,
) -> str:
    def file_identity(path: Path) -> dict[str, int | str]:
        stat = path.stat()
        return {
            "path": str(path.resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }

    payload = {
        "version": 1,
        "fq1": file_identity(fq1),
        "fq2": file_identity(fq2),
        "imgt": file_identity(imgt),
        "gene": gene,
        "baseline": baseline,
        "incoming_candidates": incoming_candidates,
        "k": k,
        "max_full_alleles": max_full_alleles,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def cached_private_pair_support(
    cache_dir: Path | None,
    fq1: Path,
    fq2: Path,
    gene: str,
    baseline: list[str],
    incoming_candidates: list[str],
    imgt: Path,
    k: int,
    max_full_alleles: int,
) -> dict[str, Counter[str]]:
    if cache_dir is None:
        return private_pair_support(
            fq1, fq2, gene, baseline, incoming_candidates, imgt, k, max_full_alleles
        )

    cache_key = support_cache_key(
        fq1, fq2, gene, baseline, incoming_candidates, imgt, k, max_full_alleles
    )
    cache_path = cache_dir / cache_key[:2] / f"{cache_key}.json"
    if cache_path.exists():
        with cache_path.open() as handle:
            saved = json.load(handle)
        return {
            incoming: Counter(values)
            for incoming, values in saved["support"].items()
        }

    support = private_pair_support(
        fq1, fq2, gene, baseline, incoming_candidates, imgt, k, max_full_alleles
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_suffix(f".json.tmp")
    with temporary.open("w") as handle:
        json.dump(
            {"cache_key": cache_key, "support": support},
            handle,
            sort_keys=True,
            separators=(",", ":"),
        )
    temporary.replace(cache_path)
    return support


def propose_rescue(
    baseline: list[str],
    counts: dict[str, float],
    comparison_support: dict[str, Counter[str]],
    enable_four_distinct: bool,
    four_distinct_top_n: int,
    enable_em_gap_private_override: bool,
    em_gap_top_n: int,
    min_em_fraction: float,
    min_em_count: float,
    min_em_gap: float,
    min_candidate_private: int,
    weak_singleton_max: int,
    candidate_weak_ratio: float,
    blocked_candidates: set[str] | None = None,
    require_weak_singleton: bool = False,
):
    blocked_candidates = blocked_candidates or set()
    distinct_count = len(set(baseline))
    if (len(baseline) != 4 or not counts
            or distinct_count not in ({3, 4} if enable_four_distinct else {3})):
        return baseline, "ineligible", "", "", 0.0, []
    total = sum(counts.values()) or 1.0
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    if len(ranked) < 4:
        return baseline, "ineligible", "", "", 0.0, []
    fifth_count = ranked[4][1] if len(ranked) > 4 else 0.0
    gap = ranked[3][1] / fifth_count if fifth_count else float("inf")
    em_gap_override = (
        distinct_count == 3
        and gap < min_em_gap
        and enable_em_gap_private_override
    )
    if distinct_count == 4:
        top_n = four_distinct_top_n
    elif em_gap_override:
        top_n = em_gap_top_n
    else:
        top_n = 4
    candidate_rows = [
        row for row in ranked[:top_n]
        if row[0] not in baseline and row[0] not in blocked_candidates
    ]
    if not candidate_rows:
        return baseline, f"no_unselected_top{top_n}", "", "", 0.0, []
    if distinct_count == 3 and not em_gap_override:
        candidate_rows = [
            row for row in candidate_rows
            if row[1] >= min_em_count and row[1] / total >= min_em_fraction
        ]
        if not candidate_rows:
            return baseline, "em_gate", "", "", 0.0, []
    supported_rows = [
        row for row in candidate_rows
        if comparison_support.get(row[0], Counter())[row[0]] >= min_candidate_private
    ]
    if supported_rows:
        candidate, candidate_count = max(supported_rows, key=lambda row: (row[1], row[0]))
    else:
        candidate, candidate_count = max(
            candidate_rows,
            key=lambda row: (
                comparison_support.get(row[0], Counter())[row[0]],
                row[1],
                row[0],
            ),
        )
    support = comparison_support.get(candidate, Counter())
    candidate_fraction = candidate_count / total
    if distinct_count == 3 and gap < min_em_gap and not em_gap_override:
        return baseline, "em_gate", candidate, "", candidate_fraction, []
    if support[candidate] < min_candidate_private:
        return baseline, "private_gate", candidate, "", candidate_fraction, [row[0] for row in candidate_rows]

    copy_counts = Counter(baseline)
    singletons = [allele for allele in baseline if copy_counts[allele] == 1]
    weakest = min(singletons, key=lambda allele: (support[allele], allele))
    if (support[weakest] <= weak_singleton_max
            and support[candidate] >= candidate_weak_ratio * max(1, support[weakest])):
        replaced = weakest
    elif require_weak_singleton:
        return baseline, "weak_singleton_gate", candidate, weakest, candidate_fraction, [row[0] for row in candidate_rows]
    elif em_gap_override:
        return baseline, "em_gap_replacement_gate", candidate, weakest, candidate_fraction, [row[0] for row in candidate_rows]
    elif distinct_count == 4:
        return baseline, "four_distinct_replacement_gate", candidate, weakest, candidate_fraction, [row[0] for row in candidate_rows]
    else:
        replaced = next(allele for allele, copies in copy_counts.items() if copies == 2)
    rescued = list(baseline)
    rescued[rescued.index(replaced)] = candidate
    return rescued, "rescue", candidate, replaced, candidate_fraction, [row[0] for row in candidate_rows]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench-root", required=True, type=Path)
    parser.add_argument("--experiment-glob", default="accuracy_main_v4_shard*")
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES))
    parser.add_argument("--run-suffix", default="",
                        help="optional run experiment suffix, for example _baseline")
    parser.add_argument("--enable-four-distinct", action="store_true")
    parser.add_argument("--four-distinct-top-n", type=int, default=8)
    parser.add_argument("--enable-em-gap-private-override", action="store_true")
    parser.add_argument("--em-gap-top-n", type=int, default=8)
    parser.add_argument("--enable-second-pass", action="store_true")
    parser.add_argument("--second-min-candidate-private", type=int, default=50)
    parser.add_argument("--second-weak-singleton-max", type=int, default=5)
    parser.add_argument("--second-candidate-weak-ratio", type=float, default=5.0)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--support-cache",
        type=Path,
        help="cache directory for reusable per-locus private-read support",
    )
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-full-alleles", type=int, default=25)
    parser.add_argument("--min-em-fraction", type=float, default=0.005)
    parser.add_argument("--min-em-count", type=float, default=20.0)
    parser.add_argument("--min-em-gap", type=float, default=1.5)
    parser.add_argument("--min-candidate-private", type=int, default=30)
    parser.add_argument("--weak-singleton-max", type=int, default=10)
    parser.add_argument("--candidate-weak-ratio", type=float, default=3.0)
    args = parser.parse_args()

    truth = load_truth(args.bench_root, args.experiment_glob, args.genes)
    calls = load_calls(args.bench_root, args.experiment_glob, args.run_suffix, args.genes)
    rows = []
    aggregate = defaultdict(Counter)
    for key in sorted(truth):
        experiment, condition, sample, gene = key
        baseline = calls.get(key, [])
        run_root = args.bench_root / "runs" / f"{experiment}{args.run_suffix}" / condition / sample
        sample_root = run_root / "spechla_out" / sample
        counts = read_em_counts(sample_root / "em_refine" / f"{gene}.tf_counts.tsv")
        preliminary, reason, candidate, replaced, candidate_fraction, evidence_candidates = propose_rescue(
            baseline, counts, {}, args.enable_four_distinct,
            args.four_distinct_top_n,
            args.enable_em_gap_private_override, args.em_gap_top_n,
            args.min_em_fraction, args.min_em_count,
            args.min_em_gap, 0, args.weak_singleton_max, args.candidate_weak_ratio,
        )
        comparison_support: dict[str, Counter[str]] = {}
        if evidence_candidates:
            short = gene.replace("HLA-", "")
            comparison_support = cached_private_pair_support(
                args.support_cache,
                sample_root / f"{short}.R1.fq.gz",
                sample_root / f"{short}.R2.fq.gz",
                gene,
                baseline,
                evidence_candidates,
                args.imgt,
                args.k,
                args.max_full_alleles,
            )
        rescued, reason, candidate, replaced, candidate_fraction, evidence_candidates = propose_rescue(
            baseline, counts, comparison_support, args.enable_four_distinct,
            args.four_distinct_top_n,
            args.enable_em_gap_private_override, args.em_gap_top_n,
            args.min_em_fraction, args.min_em_count,
            args.min_em_gap, args.min_candidate_private, args.weak_singleton_max,
            args.candidate_weak_ratio,
        )
        second_reason = "not_run"
        second_candidate = ""
        second_replaced = ""
        second_candidate_fraction = 0.0
        second_support: Counter[str] = Counter()
        if args.enable_second_pass and reason == "rescue":
            (
                _second_preliminary,
                second_reason,
                second_candidate,
                second_replaced,
                second_candidate_fraction,
                second_evidence_candidates,
            ) = propose_rescue(
                rescued, counts, {}, args.enable_four_distinct,
                args.four_distinct_top_n,
                args.enable_em_gap_private_override, args.em_gap_top_n,
                args.min_em_fraction, args.min_em_count,
                args.min_em_gap, 0, args.second_weak_singleton_max,
                args.second_candidate_weak_ratio,
                blocked_candidates={replaced}, require_weak_singleton=True,
            )
            second_comparison_support: dict[str, Counter[str]] = {}
            if second_evidence_candidates:
                short = gene.replace("HLA-", "")
                second_comparison_support = cached_private_pair_support(
                    args.support_cache,
                    sample_root / f"{short}.R1.fq.gz",
                    sample_root / f"{short}.R2.fq.gz",
                    gene,
                    rescued,
                    second_evidence_candidates,
                    args.imgt,
                    args.k,
                    args.max_full_alleles,
                )
            (
                second_rescued,
                second_reason,
                second_candidate,
                second_replaced,
                second_candidate_fraction,
                _second_evidence_candidates,
            ) = propose_rescue(
                rescued, counts, second_comparison_support, args.enable_four_distinct,
                args.four_distinct_top_n,
                args.enable_em_gap_private_override, args.em_gap_top_n,
                args.min_em_fraction, args.min_em_count,
                args.min_em_gap, args.second_min_candidate_private,
                args.second_weak_singleton_max, args.second_candidate_weak_ratio,
                blocked_candidates={replaced}, require_weak_singleton=True,
            )
            second_support = second_comparison_support.get(second_candidate, Counter())
            if second_reason == "rescue":
                rescued = second_rescued
                reason = "rescue_twice"
        expected = truth[key]
        old_hits = multiset_hits(expected, baseline)
        new_hits = multiset_hits(expected, rescued)
        selected_support = comparison_support.get(candidate, Counter())
        row = {
            "experiment": experiment,
            "condition": condition,
            "sample_id": sample,
            "gene": gene,
            "baseline_hits": old_hits,
            "private_v2_hits": new_hits,
            "delta_hits": new_hits - old_hits,
            "decision": reason,
            "candidate": candidate or ".",
            "replaced": replaced or ".",
            "candidate_em_fraction": f"{candidate_fraction:.6f}",
            "candidate_private_pairs": selected_support[candidate] if candidate else 0,
            "replaced_private_pairs": selected_support[replaced] if replaced else 0,
            "second_decision": second_reason,
            "second_candidate": second_candidate or ".",
            "second_replaced": second_replaced or ".",
            "second_candidate_em_fraction": f"{second_candidate_fraction:.6f}",
            "second_candidate_private_pairs": second_support[second_candidate] if second_candidate else 0,
            "second_replaced_private_pairs": second_support[second_replaced] if second_replaced else 0,
            "baseline": ",".join(baseline),
            "private_v2": ",".join(rescued),
        }
        rows.append(row)
        for group in ((gene, condition), (gene, "ALL"), ("ABC", "ALL")):
            stats = aggregate[group]
            stats["loci"] += 1
            stats["old"] += old_hits
            stats["new"] += new_hits
            stats["rescued"] += reason in {"rescue", "rescue_twice"}
            stats["improved"] += new_hits > old_hits
            stats["regressed"] += new_hits < old_hits

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("gene\tcondition\tloci\tbaseline_recall\tprivate_v2_recall\trescued\timproved\tregressed")
    for (gene, condition), stats in sorted(aggregate.items()):
        denominator = 4 * stats["loci"]
        print(
            f"{gene}\t{condition}\t{stats['loci']}\t{stats['old']/denominator:.4f}\t"
            f"{stats['new']/denominator:.4f}\t{stats['rescued']}\t"
            f"{stats['improved']}\t{stats['regressed']}"
        )
    print(f"details\t{args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())