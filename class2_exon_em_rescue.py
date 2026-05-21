#!/usr/bin/env python3
"""Fast class-II exon-candidate + EM rescue diagnostics.

This intentionally avoids expensive IMGT genotype panel construction. It tests:
  1. exon_inject_em: raw-exon-supported families are forced into EM quartet search.
  4. ensemble_gate: only override low-confidence current calls when exon support and
     EM residual agree.
  5. gene_specific_gate: per-gene support/residual thresholds.

Truth from quartet_summary.tsv is used only for reporting score deltas.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_calls import load_g_group, normalize_for_display, overlap  # noqa: E402
from iterative_remap_em import fit_4hap, quartet_residual_from_counts, two_field  # noqa: E402
from rerank_multi_strategy import read_calls_quartet, read_final, read_tf_counts, read_chi_r_fit  # noqa: E402

CLASS2_GENES = ["HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]

GENE_GATES = {
    "HLA-DRB1": dict(min_support=3, max_rank=3, min_gain=0.015, max_current=3, top_n=12),
    "HLA-DQB1": dict(min_support=3, max_rank=3, min_gain=0.020, max_current=3, top_n=12),
    "HLA-DPB1": dict(min_support=6, max_rank=4, min_gain=0.040, max_current=2, top_n=10),
}


def quartet_key(quartet: Sequence[str]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    return tuple(sorted(quartet[:2])), tuple(sorted(quartet[2:4]))


def score_pair(pred: Sequence[str], truth: Sequence[str]) -> int:
    pred_counts = Counter(pred)
    truth_counts = Counter(truth)
    return sum(min(pred_counts[key], truth_counts[key]) for key in pred_counts)


def truth_score(quartet: Sequence[str], truth_r: Sequence[str], truth_d: Sequence[str], gmap) -> int:
    pred_r = [two_field(value) for value in quartet[:2]]
    pred_d = [two_field(value) for value in quartet[2:4]]
    truth_r2 = normalize_for_display(truth_r, "2field", gmap)
    truth_d2 = normalize_for_display(truth_d, "2field", gmap)
    return score_pair(pred_r, truth_r2) + score_pair(pred_d, truth_d2)


def parse_top5(value: str) -> List[Tuple[str, int, int]]:
    out = []
    for rank, item in enumerate((value or "").split(";"), 1):
        if not item or ":" not in item:
            continue
        allele, support = item.rsplit(":", 1)
        try:
            support_i = int(float(support))
        except ValueError:
            continue
        out.append((two_field(allele), support_i, rank))
    return out


def load_exon_candidates(path: Path) -> Dict[Tuple[str, str], Dict[str, Tuple[int, int]]]:
    out: Dict[Tuple[str, str], Dict[str, Tuple[int, int]]] = defaultdict(dict)
    if not path.exists():
        return out
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["sample"], row["gene"])
            for allele, support, rank in parse_top5(row.get("top5_exon_support", "")):
                prev = out[key].get(allele)
                if prev is None or (support, -rank) > (prev[0], -prev[1]):
                    out[key][allele] = (support, rank)
    return out


def load_summary(path: Path, genes: Sequence[str]) -> List[dict]:
    gene_set = set(genes)
    with path.open() as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["gene"] in gene_set]


def add(out: List[str], value: str) -> None:
    value = two_field(value)
    if value and value != "NA" and value not in out:
        out.append(value)


def current_from_final(final_path: Path, gene: str) -> Optional[Tuple[str, str, str, str]]:
    finals = read_final(final_path)
    if gene not in finals:
        return None
    return tuple(two_field(value) for value in finals[gene])


def best_exon_for_quartet(quartet: Sequence[str], exon: Dict[str, Tuple[int, int]]) -> Tuple[int, int]:
    support_best = 0
    rank_best = 999
    for allele in quartet:
        hit = exon.get(two_field(allele))
        if hit is None:
            continue
        support, rank = hit
        if support > support_best or (support == support_best and rank < rank_best):
            support_best = support
            rank_best = rank
    return support_best, rank_best


def has_exon_new_allele(current: Sequence[str], quartet: Sequence[str], exon: Dict[str, Tuple[int, int]]) -> bool:
    current_set = set(two_field(value) for value in current)
    return any(two_field(value) not in current_set and two_field(value) in exon for value in quartet)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--raw-exon", required=True, type=Path)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--report-prefix", required=True, type=Path)
    parser.add_argument("--genes", nargs="+", default=CLASS2_GENES)
    parser.add_argument("--top-n", type=int, default=12)
    parser.add_argument("--min-frac", type=float, default=0.002)
    args = parser.parse_args()

    rows = load_summary(args.summary, args.genes)
    gmap = load_g_group(args.g_group)
    exon_candidates = load_exon_candidates(args.raw_exon)
    args.report_prefix.parent.mkdir(parents=True, exist_ok=True)

    strategies = ["exon_inject_em", "ensemble_gate", "gene_specific_gate"]
    totals = defaultdict(lambda: Counter())
    baseline_total = 0
    loci = 0
    out_tsv = args.report_prefix.with_suffix(".tsv")
    fields = [
        "strategy", "set", "sample", "gene", "status", "current_score", "chosen_score",
        "delta", "current", "chosen", "truth_R", "truth_D", "em_gain", "chosen_residual",
        "current_residual", "exon_support", "exon_rank", "reason", "exon_candidates",
    ]
    with out_tsv.open("w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            sample = row["sample"]
            gene = row["gene"]
            sample_out = args.out_root / sample
            sample_asm = args.asm_root / sample
            current = current_from_final(sample_asm / f"{sample}.final_calls.tsv", gene)
            if current is None:
                continue
            baseline = read_calls_quartet(sample_asm / gene.lower() / gene / "calls.baseline.tsv")
            if baseline is None:
                baseline = read_calls_quartet(sample_asm / gene.lower() / gene / "calls.tsv")
            tf_rows = read_tf_counts(sample_out / "em_refine" / f"{gene}.tf_counts.tsv")
            counts = dict(tf_rows)
            chi = read_chi_r_fit(sample_out / "em_refine" / f"{gene}.summary.tsv")
            if not counts or chi is None:
                continue
            exon = exon_candidates.get((sample, gene), {})
            truth_r = [value for value in row["truth_R"].split(",") if value]
            truth_d = [value for value in row["truth_D"].split(",") if value]
            current_score = truth_score(current, truth_r, truth_d, gmap)
            current_residual = quartet_residual_from_counts(counts, current, chi)
            baseline_total += current_score
            loci += 1

            force_names = set(current)
            if baseline:
                force_names.update(baseline)
            force_names.update(exon.keys())
            force_quartets = [current]
            if baseline:
                force_quartets.append(tuple(two_field(value) for value in baseline))
            top_n = GENE_GATES[gene]["top_n"]
            fit = fit_4hap(
                counts,
                chi,
                top_n=max(args.top_n, top_n),
                min_frac=args.min_frac,
                per_gene_chi=True,
                chi_lo=max(0.005, chi - 0.06),
                chi_hi=min(0.5, chi + 0.06),
                chi_step=0.01,
                chi_prior_lambda=0.10,
                force_names=force_names,
                force_quartets=force_quartets,
            )
            if fit is None or fit[0] is None:
                continue
            proposed, proposed_residual, proposed_chi = fit
            support, rank = best_exon_for_quartet(proposed, exon)
            em_gain = current_residual - proposed_residual
            exon_summary = ";".join(
                f"{allele}:{support_i}@{rank_i}"
                for allele, (support_i, rank_i) in sorted(exon.items(), key=lambda item: (-item[1][0], item[1][1], item[0]))[:5]
            )

            for strategy in strategies:
                chosen = current
                reason = "hold"
                gate = GENE_GATES[gene]
                if strategy == "exon_inject_em":
                    if em_gain >= 0.0:
                        chosen = proposed
                        reason = f"best_em chi={proposed_chi:.3f}"
                elif strategy == "ensemble_gate":
                    if (current_score <= 3 and support >= 2 and rank <= 5
                            and em_gain >= 0.015 and has_exon_new_allele(current, proposed, exon)):
                        chosen = proposed
                        reason = f"low_conf+exon+em chi={proposed_chi:.3f}"
                elif strategy == "gene_specific_gate":
                    if (current_score <= gate["max_current"] and support >= gate["min_support"]
                            and rank <= gate["max_rank"] and em_gain >= gate["min_gain"]
                            and has_exon_new_allele(current, proposed, exon)):
                        chosen = proposed
                        reason = f"gene_gate chi={proposed_chi:.3f}"
                chosen_score = truth_score(chosen, truth_r, truth_d, gmap)
                delta = chosen_score - current_score
                status = "override" if quartet_key(chosen) != quartet_key(current) else "hold"
                totals[strategy]["rows"] += 1
                totals[strategy]["current"] += current_score
                totals[strategy]["after"] += chosen_score
                totals[strategy]["override"] += int(status == "override")
                totals[strategy]["imp"] += int(delta > 0)
                totals[strategy]["reg"] += int(delta < 0)
                totals[strategy]["neu"] += int(delta == 0 and status == "override")
                writer.writerow({
                    "strategy": strategy,
                    "set": row["set"],
                    "sample": sample,
                    "gene": gene,
                    "status": status,
                    "current_score": current_score,
                    "chosen_score": chosen_score,
                    "delta": delta,
                    "current": ",".join(current),
                    "chosen": ",".join(chosen),
                    "truth_R": ",".join(normalize_for_display(truth_r, "2field", gmap)),
                    "truth_D": ",".join(normalize_for_display(truth_d, "2field", gmap)),
                    "em_gain": f"{em_gain:.6f}",
                    "chosen_residual": f"{proposed_residual:.6f}",
                    "current_residual": f"{current_residual:.6f}",
                    "exon_support": support,
                    "exon_rank": rank if rank != 999 else "NA",
                    "reason": reason,
                    "exon_candidates": exon_summary,
                })

    summary_path = args.report_prefix.with_suffix(".summary")
    with summary_path.open("w") as handle:
        handle.write(f"# classII_loci={loci} baseline={baseline_total}/{loci * 4}\n")
        handle.write("strategy\trows\toverride\timp\treg\tneu\tnet_loci\tdelta_score\tafter\n")
        for strategy in strategies:
            stats = totals[strategy]
            delta_score = stats["after"] - stats["current"]
            net_loci = stats["imp"] - stats["reg"]
            handle.write(
                f"{strategy}\t{stats['rows']}\t{stats['override']}\t{stats['imp']}\t"
                f"{stats['reg']}\t{stats['neu']}\t{net_loci}\t{delta_score}\t"
                f"{stats['after']}/{loci * 4}\n"
            )
    print(f"wrote {out_tsv}")
    print(f"wrote {summary_path}")
    print(f"baseline={baseline_total}/{loci * 4}")
    for strategy in strategies:
        stats = totals[strategy]
        print(
            f"{strategy}\toverride={stats['override']} imp={stats['imp']} reg={stats['reg']} "
            f"delta={stats['after'] - stats['current']} after={stats['after']}/{loci * 4}"
        )


if __name__ == "__main__":
    main()
