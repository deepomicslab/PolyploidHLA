#!/usr/bin/env python3
"""Run multiple single-sample quartet-calling experiments.

This diagnostic harness keeps truth out of every decision and uses truth only
for the final report. It reuses the existing Locus builder, so every strategy
sees the same candidate pool and evidence channels: EM fractions, pooled AF,
IMGT genotypes, sample-wide chi_R, and fragment phase where available.
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

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from dose_aware_quartet import joint_residual, phase_residual  # noqa: E402
from evaluate_calls import load_g_group, load_truth, normalize_for_display, overlap  # noqa: E402
from rerank_multi_strategy import (  # noqa: E402
    GENES,
    build_locus,
    em_residual,
    enum_full,
    enum_one,
    enum_swap,
    quartet_key,
    read_calls_quartet,
    two_field,
)


Strategy = Dict[str, object]


STRATEGIES: List[Strategy] = [
    dict(name="joint_full_fixed", search="full", chi="fixed",
         comps=("em", "af"), margin=0.05, prior=0.00),
    dict(name="joint_full_hier", search="full", chi="hier",
         comps=("em", "af"), margin=0.05, prior=0.00),
    dict(name="joint_full_phase_hier", search="full", chi="hier",
         comps=("em", "af", "phase"), margin=0.05, prior=0.00),
    dict(name="disagreement_resolver", search="hybrid", chi="hier",
         comps=("em", "af"), margin=0.08, prior=0.12),
    dict(name="one_position_repair", search="one", chi="hier",
         comps=("em", "af"), margin=0.04, prior=0.04),
    dict(name="rd_swap_repair", search="swap", chi="hier",
         comps=("em", "af", "phase"), margin=0.04, prior=0.02),
    dict(name="private_support_repair", search="private_one", chi="hier",
         comps=("em", "af"), margin=0.03, prior=0.02),
]


def read_summary_rows(path: Path, samples: Sequence[str], genes: Sequence[str], max_rows: int):
    sample_set = set(samples or [])
    gene_set = set(genes or [])
    rows = []
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            if sample_set and row["sample"] not in sample_set:
                continue
            if gene_set and row["gene"] not in gene_set:
                continue
            rows.append(row)
            if max_rows and len(rows) >= max_rows:
                break
    return rows


def score_pair(pred: Sequence[str], truth: Sequence[str]) -> int:
    cp = Counter(pred)
    ct = Counter(truth)
    return sum(min(cp[k], ct[k]) for k in cp)


def truth_score(q: Sequence[str], truth_p: Sequence[str], truth_d: Sequence[str], gmap) -> int:
    pred_r = [two_field(x) for x in q[:2]]
    pred_d = [two_field(x) for x in q[2:4]]
    truth_r = normalize_for_display(truth_p, "2field", gmap)
    truth_d = normalize_for_display(truth_d, "2field", gmap)
    return score_pair(pred_r, truth_r) + score_pair(pred_d, truth_d)


def chi_grid(locus, mode: str, width: float, step: float) -> np.ndarray:
    if mode == "fixed":
        return np.asarray([locus.chi_r], dtype=np.float64)
    lo = max(0.005, float(locus.chi_r) - width)
    hi = min(0.995, float(locus.chi_r) + width)
    if lo > hi:
        lo, hi = hi, lo
    return np.arange(lo, hi + 1e-9, step, dtype=np.float64)


def quartet_distance(q: Sequence[str], cur: Sequence[str]) -> float:
    same_r = score_pair([two_field(x) for x in q[:2]], [two_field(x) for x in cur[:2]])
    same_d = score_pair([two_field(x) for x in q[2:4]], [two_field(x) for x in cur[2:4]])
    return float(4 - same_r - same_d) / 4.0


def support_fraction(locus, allele: str) -> float:
    return float(locus.fractions.get(two_field(allele), 0.0))


def hybrid_candidates(locus, top_n: int) -> List[str]:
    out: List[str] = []
    def add(a):
        a = two_field(a)
        if a and a not in out:
            out.append(a)
    for a in locus.cur:
        add(a)
    # c0 order is current, baseline, then EM top alleles in the existing builder.
    for a in locus.c0[:max(top_n + 4, 8)]:
        add(a)
    for a, _ in sorted(locus.fractions.items(), key=lambda kv: -kv[1])[:top_n]:
        add(a)
    return out


def enumerate_hybrid(locus, top_n: int) -> List[Tuple[str, str, str, str]]:
    cands = hybrid_candidates(locus, top_n)
    return enum_full(cands)


def enumerate_private_one(locus) -> List[Tuple[str, str, str, str]]:
    total = sum(locus.fractions.values()) or 1.0
    top = [a for a, _ in sorted(locus.fractions.items(), key=lambda kv: -kv[1])[:6]]
    cur = list(locus.cur)
    out = []
    for idx, allele in enumerate(cur):
        frac = locus.fractions.get(allele, 0.0) / total
        expected = locus.chi_r / 2.0 if idx < 2 else (1.0 - locus.chi_r) / 2.0
        if frac > max(0.03, 0.25 * expected):
            continue
        for repl in top:
            if repl == allele:
                continue
            q = list(cur)
            q[idx] = repl
            out.append(tuple(q))
    return out


def search_space(locus, kind: str, top_n: int) -> List[Tuple[str, str, str, str]]:
    if kind == "full":
        return enum_full(locus.c0)
    if kind == "one":
        return enum_one(locus.cur, locus.c0)
    if kind == "swap":
        return enum_swap(locus.cur)
    if kind == "hybrid":
        return enumerate_hybrid(locus, top_n)
    if kind == "private_one":
        return enumerate_private_one(locus)
    raise ValueError(kind)


def component_values(q, locus, grid: np.ndarray) -> Dict[str, Tuple[float, float, int]]:
    geno = locus.geno_c0
    vals: Dict[str, Tuple[float, float, int]] = {}
    vals["em"] = (em_residual(q, locus.fractions, locus.chi_r), locus.chi_r, len(locus.fractions))
    af, af_chi, af_n = joint_residual(q, geno, locus.obs_af, locus.weight, grid,
                                      scorer="binom_nll", locus=locus)
    vals["af"] = (af / max(1, af_n), af_chi, af_n)
    ph, ph_chi, ph_n = phase_residual(q, locus, grid)
    vals["phase"] = (ph, ph_chi, ph_n)
    return vals


def normalize_component(raw: Dict[Tuple[str, str, str, str], Dict[str, Tuple[float, float, int]]],
                        comp: str) -> Dict[Tuple[str, str, str, str], float]:
    values = [v[comp][0] for v in raw.values() if math.isfinite(v[comp][0])]
    if not values:
        return {q: 0.0 for q in raw}
    arr = np.asarray(values, dtype=np.float64)
    lo = float(np.min(arr))
    hi = float(np.quantile(arr, 0.75))
    scale = max(hi - lo, 1e-9)
    out = {}
    for q, vals in raw.items():
        v = vals[comp][0]
        out[q] = 0.0 if not math.isfinite(v) else (v - lo) / scale
    return out


def evaluate_strategy(locus, strategy: Strategy, args) -> Dict[str, object]:
    qs = search_space(locus, str(strategy["search"]), args.hybrid_top)
    qs.append(tuple(locus.cur))
    dedup = []
    seen = set()
    for q in qs:
        q = tuple(two_field(x) for x in q)
        key = quartet_key(q)
        if key in seen:
            continue
        if any(a not in locus.geno_c0 for a in q):
            continue
        seen.add(key)
        dedup.append(q)
    if not dedup:
        return dict(status="skip", reason="no_candidates")
    grid = chi_grid(locus, str(strategy["chi"]), args.chi_width, args.chi_step)
    raw = {q: component_values(q, locus, grid) for q in dedup}
    comps = [c for c in strategy["comps"] if any(math.isfinite(raw[q][c][0]) for q in raw)]
    if not comps:
        return dict(status="skip", reason="no_components")
    norm = {c: normalize_component(raw, c) for c in comps}
    scored = []
    for q in dedup:
        score = sum(norm[c][q] for c in comps) / len(comps)
        score += float(strategy.get("prior", 0.0)) * quartet_distance(q, locus.cur)
        scored.append((score, q))
    scored.sort(key=lambda x: x[0])
    best_score, best_q = scored[0]
    second_score = scored[1][0] if len(scored) > 1 else float("inf")
    cur_q = tuple(two_field(x) for x in locus.cur)
    cur_raw_q = next((q for q in raw if quartet_key(q) == quartet_key(cur_q)), cur_q)
    cur_score = next((s for s, q in scored if quartet_key(q) == quartet_key(cur_q)), None)
    gap = second_score - best_score if math.isfinite(second_score) else float("inf")
    margin = float(strategy["margin"])
    if quartet_key(best_q) == quartet_key(cur_q):
        chosen = cur_q
        status = "hold_current_best"
    elif gap < margin:
        chosen = cur_q
        status = "hold_low_gap"
    else:
        chosen = best_q
        status = "propose"
    conflict = any(raw[best_q][c][0] > raw[cur_raw_q][c][0] for c in comps)
    return dict(
        status=status,
        chosen=chosen,
        best=best_q,
        cur=cur_q,
        best_score=best_score,
        cur_score=cur_score if cur_score is not None else float("nan"),
        gap=gap,
        comps=",".join(comps),
        best_em=raw[best_q]["em"][0],
        best_af=raw[best_q]["af"][0],
        best_phase=raw[best_q]["phase"][0],
        best_chi=raw[best_q]["af"][1],
        conflict=conflict,
        abstain=(gap < margin or conflict),
    )


def sample_set_label(row_set: str) -> str:
    return row_set.replace("set-", "").replace("set_", "")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--out-root", required=True, type=Path)
    ap.add_argument("--truth-dir", required=True, type=Path)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--imgt", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument("--sample", action="append", default=[])
    ap.add_argument("--gene", action="append", default=[])
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--top-k", type=int, default=8)
    ap.add_argument("--hybrid-top", type=int, default=5)
    ap.add_argument("--chi-width", type=float, default=0.05)
    ap.add_argument("--chi-step", type=float, default=0.01)
    args = ap.parse_args()

    rows = read_summary_rows(args.summary, args.sample, args.gene, args.max_rows)
    truths = {lab: load_truth(args.truth_dir / f"truth_typing-set-{lab}.tsv")
              for lab in ("a", "b", "c")}
    gmap = load_g_group(SCRIPT_DIR / "resources" / "spechla" / "db" / "HLA" / "hla_nom_g.txt")
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_tsv = args.out_prefix.with_suffix(".tsv")
    fields = [
        "strategy", "set", "sample", "gene", "status", "cur_score",
        "chosen_score", "best_score_truth", "delta", "gap", "abstain",
        "conflict", "best", "chosen", "current", "truth_R", "truth_D",
        "best_model_score", "cur_model_score", "best_em", "best_af",
        "best_phase", "best_chi", "components",
    ]
    agg = defaultdict(lambda: defaultdict(int))
    cur_total_by_locus = {}
    with out_tsv.open("w") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for idx, row in enumerate(rows, 1):
            sample = row["sample"]
            gene = row["gene"]
            lab = sample_set_label(row["set"])
            try:
                locus = build_locus(sample, gene, lab, args.asm_root, args.out_root,
                                    args.ref, args.bed, args.imgt, gmap, args.top_k)
            except Exception as exc:
                print(f"skip {sample} {gene}: {exc}", file=sys.stderr, flush=True)
                continue
            if locus is None:
                continue
            truth = truths[lab]
            truth_r = truth["PATIENT"][gene]
            truth_d = truth["DONOR"][gene]
            cur_q = tuple(two_field(x) for x in locus.cur)
            cur_truth_score = truth_score(cur_q, truth_r, truth_d, gmap)
            cur_total_by_locus[(sample, gene)] = cur_truth_score
            for strat in STRATEGIES:
                result = evaluate_strategy(locus, strat, args)
                if result.get("status") == "skip":
                    continue
                chosen = result["chosen"]
                best = result["best"]
                chosen_truth_score = truth_score(chosen, truth_r, truth_d, gmap)
                best_truth_score = truth_score(best, truth_r, truth_d, gmap)
                delta = chosen_truth_score - cur_truth_score
                name = str(strat["name"])
                agg[name]["rows"] += 1
                agg[name]["after"] += chosen_truth_score
                agg[name]["cur"] += cur_truth_score
                agg[name]["prop"] += 1 if result["status"] == "propose" else 0
                agg[name]["imp"] += 1 if delta > 0 else 0
                agg[name]["reg"] += 1 if delta < 0 else 0
                agg[name]["neu"] += 1 if delta == 0 and result["status"] == "propose" else 0
                agg[name]["abstain"] += 1 if result["abstain"] else 0
                writer.writerow({
                    "strategy": name,
                    "set": lab,
                    "sample": sample,
                    "gene": gene,
                    "status": result["status"],
                    "cur_score": cur_truth_score,
                    "chosen_score": chosen_truth_score,
                    "best_score_truth": best_truth_score,
                    "delta": delta,
                    "gap": f"{result['gap']:.5f}",
                    "abstain": int(bool(result["abstain"])),
                    "conflict": int(bool(result["conflict"])),
                    "best": ",".join(best),
                    "chosen": ",".join(chosen),
                    "current": ",".join(cur_q),
                    "truth_R": ",".join(normalize_for_display(truth_r, "2field", gmap)),
                    "truth_D": ",".join(normalize_for_display(truth_d, "2field", gmap)),
                    "best_model_score": f"{result['best_score']:.5f}",
                    "cur_model_score": f"{result['cur_score']:.5f}",
                    "best_em": f"{result['best_em']:.5f}",
                    "best_af": f"{result['best_af']:.5f}",
                    "best_phase": f"{result['best_phase']:.5f}",
                    "best_chi": f"{result['best_chi']:.4f}",
                    "components": result["comps"],
                })
            print(f"[{idx}/{len(rows)}] {sample} {gene} cur={cur_truth_score}/4", flush=True)

    summary_path = args.out_prefix.with_suffix(".summary")
    baseline = sum(cur_total_by_locus.values())
    with summary_path.open("w") as fh:
        fh.write(f"# loci={len(cur_total_by_locus)} baseline={baseline}\n")
        fh.write("strategy\trows\tprop\timp\treg\tneu\tabstain\tnet_loci\tdelta_score\tafter\n")
        for strat in STRATEGIES:
            name = str(strat["name"])
            a = agg[name]
            delta_score = a["after"] - a["cur"]
            net_loci = a["imp"] - a["reg"]
            fh.write(f"{name}\t{a['rows']}\t{a['prop']}\t{a['imp']}\t{a['reg']}\t"
                     f"{a['neu']}\t{a['abstain']}\t{net_loci}\t{delta_score}\t{a['after']}\n")
    print(f"wrote {out_tsv}")
    print(f"wrote {summary_path}")
    print(f"baseline={baseline}")
    for strat in STRATEGIES:
        name = str(strat["name"])
        a = agg[name]
        print(f"{name}\tprop={a['prop']} imp={a['imp']} reg={a['reg']} "
              f"delta={a['after'] - a['cur']} after={a['after']}")


if __name__ == "__main__":
    main()