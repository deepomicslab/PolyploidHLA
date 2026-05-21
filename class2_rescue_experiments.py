#!/usr/bin/env python3
"""Class-II rescue strategy experiments.

This diagnostic harness tests several truth-free rescue ideas on class-II loci:

1. exon_inject_em: inject raw-exon-supported families into the candidate pool,
   then choose the best EM dose residual quartet.
4. ensemble_gate: override only when EM and an independent AF/discriminating
   evidence channel both improve.
5. gene_specific_ensemble: same as ensemble_gate with stricter per-gene gates.
6. disc_panel: use only discriminating pooled-continuous sites plus EM.

Truth from quartet_summary.tsv is used only for post-hoc scoring.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from caller_free_4hap import collect_obs_af, imgt_genotypes_at_sites, load_imgt, parse_bed  # noqa: E402
from evaluate_calls import load_g_group, normalize_for_display, overlap  # noqa: E402
from rerank_multi_strategy import (  # noqa: E402
    af_disc,
    af_l1,
    build_candidates_C1,
    collect_phase_evidence,
    em_residual,
    enum_full,
    phase_score,
    read_calls_quartet,
    read_chi_r_fit,
    read_final,
    read_tf_counts,
    two_field,
)

CLASS2_GENES = ["HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]

GENE_GATES = {
    "HLA-DRB1": {
        "exon_min_support": 3,
        "exon_max_rank": 3,
        "em_gain": 0.020,
        "af_gain": 0.002,
        "disc_gain": 0.002,
        "phase_gain": 0.0002,
        "max_current_score": 3,
    },
    "HLA-DQB1": {
        "exon_min_support": 3,
        "exon_max_rank": 3,
        "em_gain": 0.030,
        "af_gain": 0.003,
        "disc_gain": 0.003,
        "phase_gain": 0.0003,
        "max_current_score": 3,
    },
    "HLA-DPB1": {
        "exon_min_support": 6,
        "exon_max_rank": 4,
        "em_gain": 0.050,
        "af_gain": 0.005,
        "disc_gain": 0.005,
        "phase_gain": 0.0005,
        "max_current_score": 2,
    },
}


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


def quartet_key(quartet: Sequence[str]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    return tuple(sorted(quartet[:2])), tuple(sorted(quartet[2:4]))


def parse_top5_support(value: str) -> List[Tuple[str, int, int]]:
    out = []
    for rank, item in enumerate((value or "").split(";"), 1):
        if not item or ":" not in item:
            continue
        allele, count = item.rsplit(":", 1)
        try:
            support = int(float(count))
        except ValueError:
            continue
        out.append((two_field(allele), support, rank))
    return out


def load_exon_candidates(path: Path) -> Dict[Tuple[str, str], Dict[str, Tuple[int, int]]]:
    candidates: Dict[Tuple[str, str], Dict[str, Tuple[int, int]]] = defaultdict(dict)
    if not path.exists():
        return candidates
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["sample"], row["gene"])
            for allele, support, rank in parse_top5_support(row.get("top5_exon_support", "")):
                prev = candidates[key].get(allele)
                if prev is None or (support, -rank) > (prev[0], -prev[1]):
                    candidates[key][allele] = (support, rank)
    return candidates


def load_summary_rows(path: Path, genes: Iterable[str]) -> List[dict]:
    gene_set = set(genes)
    with path.open() as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["gene"] in gene_set]


def add_unique(out: List[str], allele: str) -> None:
    allele = two_field(allele)
    if allele and allele != "NA" and allele not in out:
        out.append(allele)


def make_candidate_pool(current: Sequence[str], baseline: Optional[Sequence[str]],
                        tf_rows: Sequence[Tuple[str, float]],
                        exon_rows: Dict[str, Tuple[int, int]],
                        gene: str, top_k: int, include_exon: bool) -> List[str]:
    out: List[str] = []
    for allele in current:
        add_unique(out, allele)
    if baseline:
        for allele in baseline:
            add_unique(out, allele)
    for allele, _frac in tf_rows[:top_k]:
        add_unique(out, allele)
    if include_exon:
        gate = GENE_GATES[gene]
        for allele, (support, rank) in sorted(exon_rows.items(), key=lambda item: (-item[1][0], item[1][1], item[0])):
            if support >= gate["exon_min_support"] and rank <= gate["exon_max_rank"]:
                add_unique(out, allele)
    return out


def read_locus(row: dict, args, gmap, include_exon: bool, exon_candidates):
    sample = row["sample"]
    gene = row["gene"]
    sample_out = args.out_root / sample
    sample_asm = args.asm_root / sample
    final_path = sample_asm / f"{sample}.final_calls.tsv"
    finals = read_final(final_path)
    if gene not in finals:
        return None
    current = finals[gene]
    em_dir = sample_out / "em_refine"
    summary_path = em_dir / f"{gene}.summary.tsv"
    chi_r = read_chi_r_fit(summary_path)
    if chi_r is None:
        return None
    tf_rows = read_tf_counts(em_dir / f"{gene}.tf_counts.tsv")
    fractions = {allele: frac for allele, frac in tf_rows}
    baseline_path = sample_asm / gene.lower() / gene / "calls.baseline.tsv"
    if not baseline_path.exists():
        baseline_path = sample_asm / gene.lower() / gene / "calls.tsv"
    baseline = read_calls_quartet(baseline_path)
    exon_rows = exon_candidates.get((sample, gene), {})
    candidates = make_candidate_pool(current, baseline, tf_rows, exon_rows, gene, args.top_k, include_exon)
    if len(candidates) < 2:
        return None

    contig = gene.replace("-", "_")
    start, end = parse_bed(str(args.bed), contig)
    if start is None:
        return None
    pc_vcf = sample_out / f"{sample}.pooled_continuous.vcf.gz"
    if not pc_vcf.exists():
        return None
    obs = collect_obs_af(str(pc_vcf), contig, (start, end))
    sites = sorted(obs.keys())
    if len(sites) < 30:
        return None
    obs_af = np.asarray([obs[key][0] for key in sites], dtype=float)
    weight = np.minimum(np.asarray([obs[key][1] for key in sites], dtype=float), 200.0)

    ref_seq = pysam.FastaFile(str(args.ref)).fetch(contig, start, end).upper()
    prefix = gene.replace("HLA-", "") + "*"
    imgt = load_imgt(str(args.imgt), prefix)
    keep = set(candidates)
    imgt_sub = {name: seq for name, seq in imgt.items() if two_field(name) in keep}
    if not imgt_sub:
        return None
    geno_all = imgt_genotypes_at_sites(imgt_sub, ref_seq, start, sites)
    candidates = [allele for allele in candidates if allele in geno_all]
    if any(allele not in candidates for allele in current):
        return None
    geno = {allele: geno_all[allele] for allele in candidates}
    stack = np.stack([geno[allele] for allele in candidates])
    sites_used, pair_counts = collect_phase_evidence(sample_out / f"{sample}.merge.bam", contig, sites)

    return {
        "sample": sample,
        "set": row["set"],
        "gene": gene,
        "current": tuple(current),
        "baseline": tuple(baseline) if baseline else None,
        "candidates": candidates,
        "candidates_g": build_candidates_C1(candidates, gmap),
        "fractions": fractions,
        "chi_r": chi_r,
        "obs_af": obs_af,
        "weight": weight,
        "geno": geno,
        "stack": stack,
        "sites": sites_used,
        "pair_counts": pair_counts,
        "exon_rows": exon_rows,
        "current_score": int(row["score2"]),
        "truth_r": [x for x in row["truth_R"].split(",") if x],
        "truth_d": [x for x in row["truth_D"].split(",") if x],
    }


def score_channels(quartet: Sequence[str], locus: dict) -> dict:
    em = em_residual(quartet, locus["fractions"], locus["chi_r"])
    af, af_n = af_l1(quartet, locus["geno"], locus["obs_af"], locus["weight"], locus["chi_r"])
    disc, disc_n = af_disc(quartet, locus["stack"], locus["geno"], locus["obs_af"], locus["weight"], locus["chi_r"])
    phase, phase_n = phase_score(quartet, locus["sites"], locus["pair_counts"], locus["geno"], locus["chi_r"])
    return {
        "em": em,
        "af": af,
        "disc": disc,
        "phase": phase,
        "af_n": af_n,
        "disc_n": disc_n,
        "phase_n": phase_n,
    }


def best_by_channel(locus: dict, candidates: Sequence[str], channel: str) -> Optional[Tuple[Tuple[str, str, str, str], dict]]:
    best = None
    seen = set()
    for quartet in enum_full(candidates):
        key = quartet_key(quartet)
        if key in seen:
            continue
        seen.add(key)
        if any(allele not in locus["geno"] for allele in quartet):
            continue
        scores = score_channels(quartet, locus)
        value = scores.get(channel)
        if value is None or not math.isfinite(value):
            continue
        if best is None or value < best[1][channel]:
            best = (tuple(quartet), scores)
    return best


def channel_gain(current_scores: dict, alt_scores: dict, channel: str) -> float:
    cur = current_scores.get(channel)
    alt = alt_scores.get(channel)
    if cur is None or alt is None or not math.isfinite(cur) or not math.isfinite(alt):
        return -float("inf")
    return float(cur - alt)


def exon_support_for_quartet(locus: dict, quartet: Sequence[str]) -> Tuple[int, int]:
    best_support = 0
    best_rank = 999
    for allele in quartet:
        support_rank = locus["exon_rows"].get(two_field(allele))
        if support_rank is None:
            continue
        support, rank = support_rank
        if support > best_support or (support == best_support and rank < best_rank):
            best_support, best_rank = support, rank
    return best_support, best_rank


def apply_strategy(name: str, locus: dict, gmap) -> Tuple[Tuple[str, ...], str, dict]:
    current = tuple(locus["current"])
    current_scores = score_channels(current, locus)
    gene = locus["gene"]
    gate = GENE_GATES[gene]
    candidates = locus["candidates"]
    reason = []

    if name == "exon_inject_em":
        best = best_by_channel(locus, candidates, "em")
        if not best:
            return current, "no_candidate", current_scores
        alt, scores = best
        gain = channel_gain(current_scores, scores, "em")
        support, rank = exon_support_for_quartet(locus, alt)
        if support >= gate["exon_min_support"] and rank <= gate["exon_max_rank"] and gain >= gate["em_gain"]:
            return alt, f"override em_gain={gain:.4f} exon={support}@{rank}", scores
        return current, f"hold em_gain={gain:.4f} exon={support}@{rank}", current_scores

    if name == "disc_panel":
        best = best_by_channel(locus, candidates, "disc")
        if not best:
            return current, "no_disc", current_scores
        alt, scores = best
        em_gain = channel_gain(current_scores, scores, "em")
        disc_gain = channel_gain(current_scores, scores, "disc")
        if em_gain >= gate["em_gain"] and disc_gain >= gate["disc_gain"]:
            return alt, f"override em_gain={em_gain:.4f} disc_gain={disc_gain:.4f}", scores
        return current, f"hold em_gain={em_gain:.4f} disc_gain={disc_gain:.4f}", current_scores

    if name in {"ensemble_gate", "gene_specific_ensemble"}:
        # Prefer the best EM candidate, but require independent support.
        best = best_by_channel(locus, candidates, "em")
        if not best:
            return current, "no_candidate", current_scores
        alt, scores = best
        em_gain = channel_gain(current_scores, scores, "em")
        af_gain = channel_gain(current_scores, scores, "af")
        disc_gain = channel_gain(current_scores, scores, "disc")
        phase_gain = channel_gain(current_scores, scores, "phase")
        support, rank = exon_support_for_quartet(locus, alt)
        independent = (
            af_gain >= gate["af_gain"]
            or disc_gain >= gate["disc_gain"]
            or phase_gain >= gate["phase_gain"]
        )
        if name == "gene_specific_ensemble":
            low_conf = locus["current_score"] <= gate["max_current_score"]
            exon_ok = support >= gate["exon_min_support"] and rank <= gate["exon_max_rank"]
        else:
            low_conf = locus["current_score"] <= 3
            exon_ok = support >= 2 and rank <= 5
        if em_gain >= gate["em_gain"] and independent and low_conf and exon_ok:
            return alt, (
                f"override em={em_gain:.4f} af={af_gain:.4f} disc={disc_gain:.4f} "
                f"phase={phase_gain:.4f} exon={support}@{rank}"
            ), scores
        return current, (
            f"hold em={em_gain:.4f} af={af_gain:.4f} disc={disc_gain:.4f} "
            f"phase={phase_gain:.4f} exon={support}@{rank} low_conf={int(low_conf)}"
        ), current_scores

    raise ValueError(name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--ref", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--imgt", required=True, type=Path)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--raw-exon", required=True, type=Path)
    parser.add_argument("--report-prefix", required=True, type=Path)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--genes", nargs="+", default=CLASS2_GENES)
    args = parser.parse_args()

    args.report_prefix.parent.mkdir(parents=True, exist_ok=True)
    gmap = load_g_group(args.g_group)
    exon_candidates = load_exon_candidates(args.raw_exon)
    summary_rows = load_summary_rows(args.summary, args.genes)
    strategies = ["exon_inject_em", "disc_panel", "ensemble_gate", "gene_specific_ensemble"]

    out_tsv = args.report_prefix.with_suffix(".tsv")
    fields = [
        "strategy", "set", "sample", "gene", "status", "current_score", "chosen_score",
        "delta", "chosen", "current", "truth_R", "truth_D", "reason",
        "n_candidates", "exon_candidates", "em", "af", "disc", "phase",
    ]
    totals = defaultdict(lambda: Counter())
    baseline_by_locus = {}

    with out_tsv.open("w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for index, row in enumerate(summary_rows, 1):
            locus = read_locus(row, args, gmap, include_exon=True, exon_candidates=exon_candidates)
            if locus is None:
                print(f"[{index}/{len(summary_rows)}] skip {row['sample']} {row['gene']}", flush=True)
                continue
            current_score = truth_score(locus["current"], locus["truth_r"], locus["truth_d"], gmap)
            baseline_by_locus[(locus["sample"], locus["gene"])] = current_score
            exon_summary = ";".join(
                f"{allele}:{support}@{rank}"
                for allele, (support, rank) in sorted(locus["exon_rows"].items(), key=lambda item: (-item[1][0], item[1][1], item[0]))[:5]
            )
            for strategy in strategies:
                chosen, reason, scores = apply_strategy(strategy, locus, gmap)
                chosen_score = truth_score(chosen, locus["truth_r"], locus["truth_d"], gmap)
                delta = chosen_score - current_score
                status = "override" if quartet_key(chosen) != quartet_key(locus["current"]) else "hold"
                totals[strategy]["rows"] += 1
                totals[strategy]["after"] += chosen_score
                totals[strategy]["current"] += current_score
                totals[strategy]["override"] += int(status == "override")
                totals[strategy]["imp"] += int(delta > 0)
                totals[strategy]["reg"] += int(delta < 0)
                totals[strategy]["neu"] += int(delta == 0 and status == "override")
                writer.writerow({
                    "strategy": strategy,
                    "set": locus["set"],
                    "sample": locus["sample"],
                    "gene": locus["gene"],
                    "status": status,
                    "current_score": current_score,
                    "chosen_score": chosen_score,
                    "delta": delta,
                    "chosen": ",".join(chosen),
                    "current": ",".join(locus["current"]),
                    "truth_R": ",".join(normalize_for_display(locus["truth_r"], "2field", gmap)),
                    "truth_D": ",".join(normalize_for_display(locus["truth_d"], "2field", gmap)),
                    "reason": reason,
                    "n_candidates": len(locus["candidates"]),
                    "exon_candidates": exon_summary,
                    "em": "" if scores.get("em") is None else f"{scores['em']:.6f}",
                    "af": "" if scores.get("af") is None else f"{scores['af']:.6f}",
                    "disc": "" if scores.get("disc") is None else f"{scores['disc']:.6f}",
                    "phase": "" if scores.get("phase") is None else f"{scores['phase']:.6f}",
                })
            print(f"[{index}/{len(summary_rows)}] {locus['sample']} {locus['gene']} cur={current_score}/4", flush=True)

    summary_path = args.report_prefix.with_suffix(".summary")
    baseline = sum(baseline_by_locus.values())
    with summary_path.open("w") as handle:
        handle.write(f"# loci={len(baseline_by_locus)} baseline_classII={baseline}/{len(baseline_by_locus) * 4}\n")
        handle.write("strategy\trows\toverride\timp\treg\tneu\tnet_loci\tdelta_score\tafter\n")
        for strategy in strategies:
            values = totals[strategy]
            delta_score = values["after"] - values["current"]
            net_loci = values["imp"] - values["reg"]
            handle.write(
                f"{strategy}\t{values['rows']}\t{values['override']}\t{values['imp']}\t"
                f"{values['reg']}\t{values['neu']}\t{net_loci}\t{delta_score}\t"
                f"{values['after']}/{len(baseline_by_locus) * 4}\n"
            )
    print(f"wrote {out_tsv}")
    print(f"wrote {summary_path}")
    print(f"baseline_classII={baseline}/{len(baseline_by_locus) * 4}")
    for strategy in strategies:
        values = totals[strategy]
        print(
            f"{strategy}\toverride={values['override']} imp={values['imp']} "
            f"reg={values['reg']} delta={values['after'] - values['current']} "
            f"after={values['after']}/{len(baseline_by_locus) * 4}"
        )


if __name__ == "__main__":
    main()
