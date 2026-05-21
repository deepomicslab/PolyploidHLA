#!/usr/bin/env python3
"""Dose-aware quartet solver (truth-free, single sample).

Given the per-locus evidence already produced by the pipeline:
  - top-K candidate 2-field alleles from EM tf_counts
  - per-site pooled-continuous AF + depth (within the gene region)
  - IMGT genotype matrix at those sites for each candidate
solve

    argmin_{R-multiset, D-multiset, chi_R}  weighted L1 (obs_AF, predicted_AF)

with explicit MULTISET enumeration (R can be {a,a}, {a,b}, etc. — no
ordering constraints inside R or inside D, but R != D pairing is
preserved).  chi_R is searched on a small grid PER quartet hypothesis,
because the optimal chi_R depends on the hypothesis (this is the
identifiability fix).

The override is applied only if the gap between the best and the second
best quartet residual exceeds a margin tau.  When ambiguous, fall back
to the pipeline's current call (no-op).

This file is a diagnostic harness:
  - build a Locus (reuses rerank_multi_strategy.build_locus)
  - solve via dose-aware solver
  - compare to pipeline current and to truth (truth used for evaluation
    only; never seen by the solver)
  - report leave-set-out-style aggregates so the tau gate can be
    cross-validated.

Output:
  diagnostics/dose_aware_<tag>.tsv      one row per (sample, gene)
  diagnostics/dose_aware_<tag>.summary  per-tau aggregate, by_set
"""
from __future__ import annotations

import argparse
import csv
import itertools
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

# Reuse the locus builder and helpers.
from rerank_multi_strategy import (  # noqa: E402
    GENES, two_field, build_locus, Locus, phase_score,
)
from evaluate_calls import load_g_group, load_truth, normalize_for_display, overlap  # noqa: E402


# -------- core solver --------

SCORERS = ("af_l1", "count_l1", "binom_nll", "poisson_nll", "phase_nll",
           "phase_plus_binom")


# -------- vectorized phase precompute --------

class PhaseTable:
    """Per-locus precomputed structures for fast phase NLL evaluation.

    For each candidate allele c and each pair p=(i,j), idx[c, p] in {0,1,2,3}
    encodes (geno[c,i], geno[c,j]) as 2*g_i + g_j.  Pairs where any
    candidate has unknown geno (-1) are excluded.

    obs_counts[p, b] = fragment count for bucket b at pair p, b in {00,01,10,11}.
    """
    __slots__ = ("idx", "obs_counts", "alleles", "n_pairs", "total_obs")

    def __init__(self, idx, obs_counts, alleles, total_obs):
        self.idx = idx               # (|C|, P) uint8
        self.obs_counts = obs_counts # (P, 4) float32
        self.alleles = alleles       # list of allele names ordered as idx rows
        self.n_pairs = idx.shape[1] if idx.size else 0
        self.total_obs = total_obs   # scalar


def build_phase_table(locus: Locus) -> Optional[PhaseTable]:
    if not locus.pair_counts:
        return None
    pairs = list(locus.pair_counts.items())
    if not pairs:
        return None
    # collect candidate alleles
    alleles = [a for a in locus.c0 if a in locus.geno_c0]
    if len(alleles) < 2:
        return None
    geno_stack = np.stack([locus.geno_c0[a] for a in alleles])  # (|C|, n_sites)
    # filter to pairs where every candidate has known genotype on BOTH sites
    n_cands = geno_stack.shape[0]
    keep_pairs = []
    keep_obs = []
    keep_idx = []
    for (i, j), mat in pairs:
        if i >= geno_stack.shape[1] or j >= geno_stack.shape[1]:
            continue
        gi = geno_stack[:, i]
        gj = geno_stack[:, j]
        if (gi == -1).any() or (gj == -1).any():
            continue
        # idx encoded
        col = (2 * gi + gj).astype(np.uint8)
        keep_idx.append(col)
        # flatten 2x2 to (4,) order: (0,0),(0,1),(1,0),(1,1)
        flat = mat.reshape(-1).astype(np.float32)
        keep_obs.append(flat)
        keep_pairs.append((i, j))
    if len(keep_pairs) < 5:
        return None
    idx = np.stack(keep_idx, axis=1)             # (|C|, P)
    obs_counts = np.stack(keep_obs, axis=0)      # (P, 4)
    total_obs = float(obs_counts.sum()) + 1e-9
    return PhaseTable(idx, obs_counts, alleles, total_obs)


_PHASE_CACHE: Dict[int, Optional[PhaseTable]] = {}


def get_phase_table(locus: Locus) -> Optional[PhaseTable]:
    key = id(locus)
    if key not in _PHASE_CACHE:
        _PHASE_CACHE[key] = build_phase_table(locus)
    return _PHASE_CACHE[key]


def phase_residual(quartet, locus: Locus, chi_grid: np.ndarray,
                   smooth: float = 1e-3
                   ) -> Tuple[float, float, int]:
    """Vectorised per-fragment co-occurrence NLL with chi_R grid search.

    NLL is normalised by total fragment-pair counts so it is on the same
    'mean negative log-likelihood per fragment-observation' scale across
    loci of different depth.
    """
    pt = get_phase_table(locus)
    if pt is None:
        return float("inf"), 0.0, 0
    name_to_row = {a: r for r, a in enumerate(pt.alleles)}
    rows = [name_to_row.get(a) for a in quartet]
    if any(r is None for r in rows):
        return float("inf"), 0.0, 0
    P = pt.n_pairs
    if P == 0:
        return float("inf"), 0.0, 0
    # For each hap, per-pair bucket index in {0,1,2,3}
    iR1 = pt.idx[rows[0]]; iR2 = pt.idx[rows[1]]
    iD1 = pt.idx[rows[2]]; iD2 = pt.idx[rows[3]]
    # mass per (pair, bucket) for R = count of haps in {R1,R2} that map to bucket b
    arange = np.arange(P)
    m_R = np.zeros((P, 4), dtype=np.float32)
    m_R[arange, iR1] += 1.0; m_R[arange, iR2] += 1.0
    m_D = np.zeros((P, 4), dtype=np.float32)
    m_D[arange, iD1] += 1.0; m_D[arange, iD2] += 1.0
    chi = np.asarray(chi_grid, dtype=np.float32).reshape(-1, 1, 1)  # (C,1,1)
    pred = (chi / 2.0) * m_R[None, :, :] + ((1.0 - chi) / 2.0) * m_D[None, :, :]
    pred = pred + smooth
    pred = pred / pred.sum(axis=2, keepdims=True)
    # NLL = -Σ obs * log(pred)
    log_pred = np.log(pred)
    nll = -(pt.obs_counts[None, :, :] * log_pred).sum(axis=(1, 2))  # (C,)
    nll_norm = nll / pt.total_obs
    k = int(np.argmin(nll_norm))
    return float(nll_norm[k]), float(np.asarray(chi_grid).flatten()[k]), int(P)


def joint_residual(quartet, geno: Dict[str, np.ndarray],
                   obs_af: np.ndarray, weight: np.ndarray,
                   chi_grid: np.ndarray,
                   scorer: str = "binom_nll",
                   eps: float = 5e-3,
                   locus: Optional[Locus] = None) -> Tuple[float, float, int]:
    """Score one quartet, jointly optimising chi_r over chi_grid.

    Scorers:
      af_l1            : depth-weighted, depth-normalised |AF_obs-AF_pred|
                         (rate-space, dilutes high-depth informative sites).
      count_l1         : unnormalised |alt_count - dp*pred|.
      binom_nll        : binomial NLL on (alt, dp) per site.
      poisson_nll      : Poisson NLL on alt with rate=dp*pred.
      phase_nll        : per-fragment co-occurrence NLL across pairs of
                         heterozygous sites (uses locus.pair_counts).
      phase_plus_binom : phase_nll + alpha * binom_nll, alpha auto-scaled
                         so the two channels contribute on similar
                         magnitudes (median of best-quartet residuals
                         is not stable here; we use a fixed 1:1 mix in
                         normalised units — see implementation).
    """
    if any(a not in geno for a in quartet):
        return float("inf"), 0.0, 0
    G = [geno[a] for a in quartet]
    valid = (G[0] != -1) & (G[1] != -1) & (G[2] != -1) & (G[3] != -1)
    n = int(valid.sum())
    if n < 30:
        return float("inf"), 0.0, n
    if scorer == "phase_nll":
        return phase_residual(quartet, locus, chi_grid)
    if scorer == "phase_plus_binom":
        ph = phase_residual(quartet, locus, chi_grid)
        # also compute binom NLL at the chi minimising phase, for consistency
        if ph[0] == float("inf"):
            return ph
        # one-pass binom at chi=ph[1]
        chi_arr = np.array([ph[1]], dtype=np.float64)
        bn = joint_residual(quartet, geno, obs_af, weight, chi_arr,
                            scorer="binom_nll", eps=eps, locus=locus)
        # combine in standardized units: just sum (phase NLL is mean per
        # read, binom NLL is sum across sites; rescale binom by 1/n_sites
        # so both are roughly per-observation)
        nrm_binom = bn[0] / max(1.0, float(n))
        return float(ph[0] + nrm_binom), float(ph[1]), int(ph[2])
    R_dose = 0.5 * (G[0].astype(np.float32) + G[1].astype(np.float32))
    D_dose = 0.5 * (G[2].astype(np.float32) + G[3].astype(np.float32))
    chi = chi_grid.reshape(-1, 1)            # (C,1)
    pred = R_dose[None, :] * chi + D_dose[None, :] * (1.0 - chi)  # (C,S)
    pred = np.clip(pred, eps, 1.0 - eps)
    dp = weight.astype(np.float32)            # (S,)
    alt = obs_af.astype(np.float32) * dp      # (S,)
    mask = valid.astype(np.float32)           # (S,)
    if scorer == "af_l1":
        diff = np.abs(obs_af[None, :] - pred)
        use = (weight * valid).astype(np.float32)
        denom = float(use.sum()) + 1e-9
        res = (diff * use[None, :]).sum(axis=1) / denom
    elif scorer == "count_l1":
        diff = np.abs(alt[None, :] - pred * dp[None, :])
        res = (diff * mask[None, :]).sum(axis=1)
    elif scorer == "poisson_nll":
        rate = pred * dp[None, :]
        # -log P(alt|rate) ∝ rate - alt*log(rate)  (drop alt! constant)
        contrib = rate - alt[None, :] * np.log(rate + 1e-12)
        res = (contrib * mask[None, :]).sum(axis=1)
    else:  # binom_nll
        ref = dp[None, :] - alt[None, :]
        contrib = -(alt[None, :] * np.log(pred)
                    + ref * np.log(1.0 - pred))
        res = (contrib * mask[None, :]).sum(axis=1)
    k = int(np.argmin(res))
    return float(res[k]), float(chi_grid[k]), n


def dose_aware_quartet(locus: Locus,
                       chi_grid: np.ndarray,
                       cset_key: str = "c0",
                       scorer: str = "binom_nll"
                       ) -> List[Tuple[Tuple[str, str, str, str], float, float, int]]:
    """Enumerate multiset quartets over candidate pool and rank by
    joint-chi residual.  No ordering constraint within R or within D.
    Returns sorted list of (quartet, residual, chi_used, n_sites)."""
    cands = locus.c0 if cset_key == "c0" else locus.c1
    geno = locus.geno_c0 if cset_key == "c0" else locus.geno_c1
    cands = [a for a in cands if a in geno]
    if len(cands) < 2:
        return []
    pairs = list(itertools.combinations_with_replacement(cands, 2))
    out = []
    for R in pairs:
        for D in pairs:
            q = (R[0], R[1], D[0], D[1])
            r, c, n = joint_residual(q, geno, locus.obs_af, locus.weight,
                                     chi_grid, scorer=scorer, locus=locus)
            if r == float("inf"):
                continue
            out.append((q, r, c, n))
    out.sort(key=lambda x: x[1])
    return out


def solve_with_gate(locus: Locus, chi_grid: np.ndarray, tau: float,
                    cset_key: str = "c0"
                    ) -> Tuple[Optional[List[str]], dict]:
    """Return (proposed_quartet_or_None, info).  None means 'do not
    override'.  info contains diagnostics."""
    ranked = dose_aware_quartet(locus, chi_grid, cset_key)
    if len(ranked) < 2:
        return None, dict(reason="too_few_quartets", n=len(ranked))
    best, second = ranked[0], ranked[1]
    gap = second[1] - best[1]
    cur_q = tuple(locus.cur)
    info = dict(
        best=best[0], best_res=best[1], best_chi=best[2], best_n=best[3],
        second=second[0], second_res=second[1], gap=gap,
        cur=cur_q,
    )
    # cur residual for context
    cur_r, cur_c, cur_n = joint_residual(cur_q, locus.geno_c0, locus.obs_af,
                                         locus.weight, chi_grid)
    info["cur_res"] = cur_r
    info["cur_chi"] = cur_c
    if best[0] == cur_q:
        info["reason"] = "best_equals_current"
        return None, info
    if gap <= tau:
        info["reason"] = "ambiguous_gap_below_tau"
        return None, info
    info["reason"] = "override"
    return list(best[0]), info


# -------- truth-side scoring --------

def score_quartet_2field(q: Sequence[str], truth_p, truth_d, gmap) -> int:
    """2-field hits: greedy match R against truth_p (PATIENT) and D against truth_d."""
    p = sorted([two_field(q[0]), two_field(q[1])])
    d = sorted([two_field(q[2]), two_field(q[3])])
    tp = sorted([two_field(x) for x in truth_p])
    td = sorted([two_field(x) for x in truth_d])
    return overlap(tp, p) + overlap(td, d)


# -------- driver --------

def discover_samples(fqs_root: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for setdir, lab in [("set_a", "a"), ("set B", "b"), ("set C", "c")]:
        d = fqs_root / setdir
        if not d.exists():
            continue
        for fq in d.glob("*_R1_001.fastq.gz"):
            sname = fq.name.replace("_R1_001.fastq.gz", "")
            out[sname] = lab
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", type=Path, required=True)
    ap.add_argument("--out-root", type=Path, required=True,
                    help="spechla_out_<tag> root (per-sample tf_counts/AF/BAM live here)")
    ap.add_argument("--ref", type=Path, required=True)
    ap.add_argument("--bed", type=Path, required=True)
    ap.add_argument("--imgt", type=Path, required=True)
    ap.add_argument("--g-group", type=Path, required=True)
    ap.add_argument("--truth-dir", type=Path, required=True,
                    help="contains truth_typing-set-{a,b,c}.tsv")
    ap.add_argument("--fqs-root", type=Path, required=True,
                    help="contains set_a/, set B/, set C/ subdirs")
    ap.add_argument("--top-k", type=int, default=12,
                    help="EM top-K alleles to keep as candidate pool")
    ap.add_argument("--chi-lo", type=float, default=0.02)
    ap.add_argument("--chi-hi", type=float, default=0.50)
    ap.add_argument("--chi-step", type=float, default=0.01)
    ap.add_argument("--taus", type=str,
                    default="0.0,0.0005,0.001,0.002,0.005,0.01",
                    help="comma-separated identifiability gate values to evaluate")
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--cset", choices=("c0", "c1"), default="c0")
    ap.add_argument("--scorers", type=str,
                    default="af_l1,count_l1,binom_nll,poisson_nll",
                    help="comma-separated subset of {af_l1,count_l1,binom_nll,poisson_nll}")
    args = ap.parse_args()

    chi_grid = np.arange(args.chi_lo, args.chi_hi + 1e-9, args.chi_step,
                         dtype=np.float64)
    taus_raw = [float(x) for x in args.taus.split(",")]
    scorer_list = [s for s in args.scorers.split(",") if s]
    for s in scorer_list:
        if s not in SCORERS:
            raise SystemExit(f"unknown scorer: {s}")
    sample_set = discover_samples(args.fqs_root)
    if not sample_set:
        raise SystemExit(f"no samples found under {args.fqs_root}")
    truths = {lab: load_truth(args.truth_dir / f"truth_typing-set-{lab}.tsv")
              for lab in ("a", "b", "c")}
    gmap = load_g_group(args.g_group)

    out_tsv = args.out_prefix.with_suffix(".tsv")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fh = out_tsv.open("w")
    w = csv.writer(fh, delimiter="\t")
    header = ["set", "sample", "gene", "cur", "cur_score"]
    for s in scorer_list:
        header += [f"{s}_best", f"{s}_second", f"{s}_best_res",
                   f"{s}_second_res", f"{s}_gap", f"{s}_chi",
                   f"{s}_best_score", f"{s}_delta"]
    w.writerow(header)

    # accumulators per (scorer, tau): use absolute taus for af_l1
    # (small floats) and SAME taus interpreted as absolute residual gap
    # for the count-based scorers (which are on count scale).  Also
    # auto-build a multiplicative tau grid for count scorers based on
    # observed gap magnitudes.
    agg = {(s, tau): dict(prop=0, imp=0, reg=0, neu=0, hold=0, after=0,
                          by_set=defaultdict(lambda: dict(prop=0, imp=0, reg=0, neu=0)),
                          by_gene=defaultdict(lambda: dict(prop=0, imp=0, reg=0, neu=0)))
           for s in scorer_list for tau in taus_raw}

    # also collect gaps per scorer for tau-grid auto-tuning
    gap_collect = {s: [] for s in scorer_list}

    baseline_total = 0
    n_loci = 0

    for sample, lab in sorted(sample_set.items()):
        tr = truths[lab]
        for gene in GENES:
            try:
                locus = build_locus(sample, gene, lab, args.asm_root,
                                    args.out_root, args.ref, args.bed,
                                    args.imgt, gmap, args.top_k)
            except Exception as e:
                print(f"# build_locus FAIL {sample} {gene}: {e}",
                      file=sys.stderr)
                continue
            if locus is None:
                continue
            n_loci += 1
            cur_q = list(locus.cur)
            cur_score = score_quartet_2field(
                cur_q,
                tr.get("PATIENT", {}).get(gene, []),
                tr.get("DONOR", {}).get(gene, []),
                gmap)
            baseline_total += cur_score
            row = [lab, sample, gene, ",".join(cur_q), cur_score]

            for scorer in scorer_list:
                ranked = dose_aware_quartet(locus, chi_grid, args.cset, scorer)
                if len(ranked) < 2:
                    row += ["", "", "", "", "", "", "", ""]
                    for tau in taus_raw:
                        agg[(scorer, tau)]["after"] += cur_score
                        agg[(scorer, tau)]["hold"] += 1
                    continue
                best, second = ranked[0], ranked[1]
                gap = second[1] - best[1]
                gap_collect[scorer].append(gap)
                best_score = score_quartet_2field(
                    list(best[0]),
                    tr.get("PATIENT", {}).get(gene, []),
                    tr.get("DONOR", {}).get(gene, []),
                    gmap)
                row += [",".join(best[0]), ",".join(second[0]),
                        f"{best[1]:.4f}", f"{second[1]:.4f}",
                        f"{gap:.4f}", f"{best[2]:.3f}",
                        best_score, best_score - cur_score]
                for tau in taus_raw:
                    a = agg[(scorer, tau)]
                    if best[0] == tuple(cur_q) or gap <= tau:
                        a["after"] += cur_score
                        a["hold"] += 1
                        continue
                    a["after"] += best_score
                    a["prop"] += 1
                    d = best_score - cur_score
                    bs = a["by_set"][lab]; bg = a["by_gene"][gene]
                    bs["prop"] += 1; bg["prop"] += 1
                    if d > 0:
                        a["imp"] += 1; bs["imp"] += 1; bg["imp"] += 1
                    elif d < 0:
                        a["reg"] += 1; bs["reg"] += 1; bg["reg"] += 1
                    else:
                        a["neu"] += 1; bs["neu"] += 1; bg["neu"] += 1
            w.writerow(row)
    fh.close()

    summary = args.out_prefix.with_suffix(".summary")
    with summary.open("w") as f:
        f.write(f"# loci processed = {n_loci}\n")
        f.write(f"# baseline (sum of cur 2field scores) = {baseline_total}\n")
        for scorer in scorer_list:
            f.write(f"\n## scorer = {scorer}\n")
            gaps = sorted(gap_collect[scorer])
            if gaps:
                qs = [gaps[int(p * (len(gaps) - 1))] for p in (0.5, 0.75, 0.9, 0.95, 0.99)]
                f.write(f"# gap quantiles 50/75/90/95/99 = "
                        f"{qs[0]:.4f} {qs[1]:.4f} {qs[2]:.4f} {qs[3]:.4f} {qs[4]:.4f}\n")
            f.write("# tau\tprop\timp\treg\tneu\thold\tnet\tafter\n")
            for tau in taus_raw:
                a = agg[(scorer, tau)]
                net = a["imp"] - a["reg"]
                f.write(f"{tau:.4f}\t{a['prop']}\t{a['imp']}\t{a['reg']}\t"
                        f"{a['neu']}\t{a['hold']}\t{net}\t{a['after']}\n")
            f.write("# by_set imp/reg/neu/net at each tau\n")
            for tau in taus_raw:
                a = agg[(scorer, tau)]
                parts = [f"tau={tau:.4f}"]
                for s in sorted(a["by_set"]):
                    bs = a["by_set"][s]
                    net = bs["imp"] - bs["reg"]
                    parts.append(f"set-{s}={bs['imp']}/{bs['reg']}/{bs['neu']}/{net}")
                f.write("\t".join(parts) + "\n")
            f.write("# by_gene imp/reg/neu/net at each tau\n")
            for tau in taus_raw:
                a = agg[(scorer, tau)]
                parts = [f"tau={tau:.4f}"]
                for g in sorted(a["by_gene"]):
                    bg = a["by_gene"][g]
                    net = bg["imp"] - bg["reg"]
                    parts.append(f"{g}={bg['imp']}/{bg['reg']}/{bg['neu']}/{net}")
                f.write("\t".join(parts) + "\n")

    print(f"wrote {out_tsv}")
    print(f"wrote {summary}")
    print(f"baseline total = {baseline_total}")
    for scorer in scorer_list:
        print(f"\n[{scorer}]")
        for tau in taus_raw:
            a = agg[(scorer, tau)]
            net = a["imp"] - a["reg"]
            print(f"  tau={tau:.4f}  prop={a['prop']:3d}  "
                  f"imp={a['imp']:2d}  reg={a['reg']:2d}  "
                  f"neu={a['neu']:2d}  net={net:+d}  after={a['after']}")


if __name__ == "__main__":
    main()
