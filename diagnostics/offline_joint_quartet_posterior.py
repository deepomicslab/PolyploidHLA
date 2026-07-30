#!/usr/bin/env python3
"""Offline source-agnostic four-copy caller using preserved EM evidence.

The caller searches allele multisets, so repeated alleles shared by the two
mixture components are valid states. For each multiset it marginalizes every
2+2 major/minor grouping and a grid of major-component fractions. It does not
assign biological source labels. Truth is read only after calling, for
evaluation.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_GENES = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")
CHI_PATTERN = re.compile(r"chi_R=([0-9.]+)")
MAJOR_PATTERN = re.compile(r"major=([0-9.]+)")


def normalize_2field(allele: str) -> str:
    allele = (allele or "").replace("HLA-", "").replace("G", "")
    if not allele or allele == "NA" or "*" not in allele:
        return ""
    gene, fields = allele.split("*", 1)
    parts = fields.split(":")
    return f"{gene}*{':'.join(parts[:2])}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def logsumexp(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def unique_mixture_groupings(quartet: tuple[str, ...]) -> list[tuple[tuple[str, str], tuple[str, str]]]:
    """Return unique major-pair/minor-pair groupings for a multiset."""
    partitions = set()
    for major_indices in itertools.combinations(range(4), 2):
        major_index_set = set(major_indices)
        major_group = tuple(sorted(quartet[index] for index in major_indices))
        minor_group = tuple(sorted(quartet[index] for index in range(4) if index not in major_index_set))
        partitions.add((major_group, minor_group))
    return sorted(partitions)


def expected_fractions(
    grouping: tuple[tuple[str, str], tuple[str, str]], major_fraction: float
) -> dict[str, float]:
    expected: dict[str, float] = defaultdict(float)
    for allele in grouping[0]:
        expected[allele] += major_fraction / 2.0
    for allele in grouping[1]:
        expected[allele] += (1.0 - major_fraction) / 2.0
    return dict(expected)


def partition_log_likelihood(
    observed: dict[str, float],
    candidates: tuple[str, ...],
    grouping: tuple[tuple[str, str], tuple[str, str]],
    major_fraction: float,
    concentration: float,
    noise: float,
) -> float:
    expected = expected_fractions(grouping, major_fraction)
    background = noise / len(candidates)
    return concentration * sum(
        observed.get(allele, 0.0)
        * math.log((1.0 - noise) * expected.get(allele, 0.0) + background)
        for allele in candidates
    )


def quartet_log_posterior(
    observed: dict[str, float],
    candidates: tuple[str, ...],
    quartet: tuple[str, ...],
    major_fraction_values: tuple[float, ...],
    major_fraction_prior: float,
    major_fraction_prior_sd: float,
    concentration: float,
    noise: float,
) -> tuple[float, float, tuple[tuple[str, str], tuple[str, str]], float]:
    groupings = unique_mixture_groupings(quartet)
    states = []
    for major_fraction in major_fraction_values:
        prior = -0.5 * (
            (major_fraction - major_fraction_prior) / major_fraction_prior_sd
        ) ** 2
        for grouping in groupings:
            likelihood = partition_log_likelihood(
                observed, candidates, grouping, major_fraction, concentration, noise
            )
            states.append((likelihood + prior, major_fraction, grouping))
    marginal = logsumexp([state[0] for state in states]) - math.log(len(states))
    best_state = max(states, key=lambda state: state[0])
    state_values = sorted((state[0] for state in states), reverse=True)
    state_gap = state_values[0] - state_values[1] if len(state_values) > 1 else math.inf
    return marginal, best_state[1], best_state[2], state_gap


def call_quartet(
    counts: dict[str, float],
    baseline: tuple[str, ...],
    major_fraction_prior: float,
    top_n: int = 8,
    max_top_n: int = 12,
    candidate_min_fraction: float = 0.002,
    major_fraction_min: float = 0.50,
    major_fraction_max: float = 0.95,
    major_fraction_step: float = 0.025,
    major_fraction_prior_sd: float = 0.08,
    concentration: float = 80.0,
    noise: float = 0.02,
) -> dict[str, object]:
    total = sum(counts.values())
    if total <= 0:
        raise ValueError("EM counts are empty")
    ranked = sorted(counts, key=lambda allele: (-counts[allele], allele))
    selected = ranked[:top_n]
    selected.extend(
        allele
        for allele in ranked[top_n:max_top_n]
        if counts[allele] / total >= candidate_min_fraction
    )
    candidates = tuple(sorted(set(selected) | set(baseline)))
    observed = {allele: counts.get(allele, 0.0) / total for allele in candidates}
    observed_total = sum(observed.values())
    observed = {allele: value / observed_total for allele, value in observed.items()}
    major_fraction_values = tuple(
        round(major_fraction_min + index * major_fraction_step, 6)
        for index in range(
            int(round((major_fraction_max - major_fraction_min) / major_fraction_step)) + 1
        )
    )
    scored = []
    for quartet in itertools.combinations_with_replacement(candidates, 4):
        score, major_fraction, grouping, state_gap = quartet_log_posterior(
            observed,
            candidates,
            quartet,
            major_fraction_values,
            major_fraction_prior,
            major_fraction_prior_sd,
            concentration,
            noise,
        )
        scored.append((score, quartet, major_fraction, grouping, state_gap))
    scored.sort(key=lambda item: (-item[0], item[1]))
    best = scored[0]
    posterior_gap = best[0] - scored[1][0] if len(scored) > 1 else math.inf
    return {
        "quartet": best[1],
        "score": best[0],
        "major_fraction": best[2],
        "major_group": best[3][0],
        "minor_group": best[3][1],
        "mixture_state_gap": best[4],
        "posterior_gap": posterior_gap,
        "candidate_count": len(candidates),
    }


def read_baseline(path: Path) -> tuple[str, ...]:
    rows = read_tsv(path)
    alleles = [normalize_2field(row.get("allele", "")) for row in rows]
    return tuple(sorted(allele for allele in alleles if allele))


def read_counts(path: Path) -> dict[str, float]:
    counts = defaultdict(float)
    for row in read_tsv(path):
        allele = normalize_2field(row.get("allele_2field", ""))
        if allele:
            counts[allele] += float(row.get("em_weight", 0.0) or 0.0)
    return dict(counts)


def read_major_fraction_prior(path: Path) -> float:
    text = path.read_text()
    match = MAJOR_PATTERN.search(text)
    if match:
        value = float(match.group(1))
        return value if value >= 0.5 else 1.0 - value
    match = CHI_PATTERN.search(text)
    if not match:
        return 0.75
    value = float(match.group(1))
    return value if value >= 0.5 else 1.0 - value


def truth_by_experiment(bench_root: Path, experiment_glob: str):
    truth = defaultdict(list)
    for path in sorted((bench_root / "truth").glob(f"{experiment_glob}/copies.tsv")):
        experiment = path.parent.name
        for row in read_tsv(path):
            truth[(experiment, row["condition"], row["sample_id"], row["gene"])].append(
                normalize_2field(row["allele_2field"])
            )
    return truth


def multiset_hits(truth: tuple[str, ...] | list[str], prediction: tuple[str, ...]) -> int:
    return sum((Counter(truth) & Counter(prediction)).values())


def discover_rows(args: argparse.Namespace) -> list[dict[str, object]]:
    truth = truth_by_experiment(args.bench_root, args.experiment_glob)
    rows = []
    run_glob = f"{args.experiment_glob}/*/SIM*/spechla_out/SIM*/em_refine"
    for em_dir in sorted((args.bench_root / "runs").glob(run_glob)):
        sample = em_dir.parent.name
        run_sample_dir = em_dir.parents[2]
        condition = run_sample_dir.parent.name
        experiment = run_sample_dir.parent.parent.name
        asm_sample = run_sample_dir / "asm_v2" / sample
        chi_path = em_dir.parent / f"{sample}.chimerism.txt"
        if not chi_path.exists():
            continue
        major_fraction_prior = read_major_fraction_prior(chi_path)
        for gene in args.genes:
            counts_path = em_dir / f"{gene}.tf_counts.tsv"
            calls_path = asm_sample / gene.lower() / gene / "calls.tsv"
            if not counts_path.exists() or not calls_path.exists():
                continue
            baseline = read_baseline(calls_path)
            if len(baseline) != 4:
                continue
            result = call_quartet(
                read_counts(counts_path),
                baseline,
                major_fraction_prior,
                top_n=args.top_n,
                max_top_n=args.max_top_n,
                candidate_min_fraction=args.candidate_min_fraction,
                major_fraction_prior_sd=args.major_fraction_prior_sd,
                concentration=args.concentration,
                noise=args.noise,
            )
            prediction = tuple(result["quartet"])
            truth_quartet = tuple(sorted(truth.get((experiment, condition, sample, gene), [])))
            hits = multiset_hits(truth_quartet, prediction) if len(truth_quartet) == 4 else ""
            baseline_hits = multiset_hits(truth_quartet, baseline) if len(truth_quartet) == 4 else ""
            rows.append({
                "experiment": experiment,
                "condition": condition,
                "sample": sample,
                "gene": gene,
                "major_fraction_prior": f"{major_fraction_prior:.4f}",
                "fitted_major_fraction": f"{float(result['major_fraction']):.4f}",
                "baseline_quartet": ",".join(baseline),
                "joint_quartet": ",".join(prediction),
                "joint_major_group": ",".join(result["major_group"]),
                "joint_minor_group": ",".join(result["minor_group"]),
                "posterior_gap": f"{float(result['posterior_gap']):.6f}",
                "mixture_state_gap": f"{float(result['mixture_state_gap']):.6f}",
                "candidate_count": result["candidate_count"],
                "truth_quartet": ",".join(truth_quartet),
                "baseline_correct": baseline_hits,
                "joint_correct": hits,
                "delta_correct": int(hits) - int(baseline_hits) if hits != "" else "",
            })
    return rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups = defaultdict(list)
    for row in rows:
        groups[("overall", "all")].append(row)
        groups[("condition", row["condition"])].append(row)
        groups[("gene", row["gene"])].append(row)
    output = []
    for (scope, name), selected in sorted(groups.items()):
        baseline = sum(int(row["baseline_correct"]) for row in selected)
        joint = sum(int(row["joint_correct"]) for row in selected)
        total = 4 * len(selected)
        output.append({
            "scope": scope,
            "name": name,
            "sample_loci": len(selected),
            "baseline_correct": baseline,
            "joint_correct": joint,
            "truth_copies": total,
            "baseline_accuracy": f"{baseline / total:.6f}",
            "joint_accuracy": f"{joint / total:.6f}",
            "delta_correct": joint - baseline,
            "improved_loci": sum(int(row["delta_correct"]) > 0 for row in selected),
            "regressed_loci": sum(int(row["delta_correct"]) < 0 for row in selected),
        })
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench-root", required=True, type=Path)
    parser.add_argument("--experiment-glob", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES))
    parser.add_argument("--top-n", type=int, default=8)
    parser.add_argument("--max-top-n", type=int, default=12)
    parser.add_argument("--candidate-min-fraction", type=float, default=0.002)
    parser.add_argument(
        "--major-fraction-prior-sd", "--chi-prior-sd",
        dest="major_fraction_prior_sd", type=float, default=0.08,
    )
    parser.add_argument("--concentration", type=float, default=80.0)
    parser.add_argument("--noise", type=float, default=0.02)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = discover_rows(args)
    if not rows:
        raise RuntimeError("no eligible runs found")
    write_tsv(args.out, rows)
    summary_path = args.summary or args.out.with_suffix(".summary.tsv")
    write_tsv(summary_path, summarize(rows))
    print(f"[joint-quartet] rows={len(rows)} out={args.out}")
    print(f"[joint-quartet] summary={summary_path}")


if __name__ == "__main__":
    main()