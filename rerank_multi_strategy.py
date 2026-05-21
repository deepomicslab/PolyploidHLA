#!/usr/bin/env python3
"""Multi-strategy truth-free rerank diagnostic.

For every (sample, gene) of the current rescue run this tool evaluates a
small set of independent override strategies, each one defined by a
fixed combination of:

  - candidate construction:
      C0:  final + baseline + EM top-10
      C1:  C0 collapsed by IMGT G-group equivalence (gene-agnostic)

  - search space relative to the current quartet:
      S_full : enumerate all (R1<=R2, D1<=D2) over candidates
      S_swap : only consider pairwise R<->D swaps of the current quartet
               (5 neighbours: swap R1/D1, R1/D2, R2/D1, R2/D2, R1<->R2 vs D1<->D2)
      S_one  : only consider replacing exactly one of (R1,R2,D1,D2) with
               another candidate

  - score channel:
      M_em : EM L1 fraction residual (read EM, chi_R fitted)
      M_af : AF weighted L1 residual on pooled-continuous biallelic sites
      M_disc : AF L1 restricted to "discriminating sites" (sites where
               at least one candidate genotype disagrees with another in
               the candidate pool); weighted by depth
      M_phase : read-pair linkage support, computed on merge.bam between
                pairs of biallelic sites covered by the same fragment

  - decision rule:
      D_both  : alt must beat current on TWO listed channels by fixed eps
      D_strict: alt must beat current on ALL listed channels by fixed eps

All thresholds are fixed a priori and not tuned on truth.

Truth is consulted only for post-hoc per-strategy precision/recall and
per-set/per-gene Δscore.

Outputs:

  diagnostics/rerank_multi_<tag>.tsv      one row per (sample, gene, strategy)
  diagnostics/rerank_multi_<tag>.summary  per-strategy improve/regress/net
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
from caller_free_4hap import (  # noqa: E402
    collect_obs_af, imgt_genotypes_at_sites, load_imgt, parse_bed)
from evaluate_calls import (  # noqa: E402
    load_g_group, load_truth, normalize_for_display, overlap)

import pysam  # noqa: E402


GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]


# ---------------- generic helpers ----------------

def two_field(allele: str) -> str:
    a = allele.replace("HLA-", "")
    if "*" not in a:
        return a
    g, rest = a.split("*", 1)
    rest = rest.rstrip("GP")
    parts = rest.split(":")
    return f"{g}*{parts[0]}:{parts[1]}" if len(parts) >= 2 else f"{g}*{parts[0]}"


def to_g_group(allele: str, gmap) -> str:
    """Map a 2-field allele to its G-group label using the WMDA table.

    The gmap loaded by evaluate_calls is keyed on full allele names; we
    look up the 2-field representative by trying common 4-digit forms
    (the gmap stored in evaluate_calls is `{full_allele: g_group}`).
    Returns the 2-field name itself when no mapping is found."""
    a = two_field(allele)
    # Try a few canonical resolutions.
    for sub in (a + ":01", a + ":01:01", a + "G"):
        g = gmap.get(sub)
        if g:
            return g
    return a


def read_final(path: Path) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["gene"]
            q = [r.get("R1_2field", ""), r.get("R2_2field", ""),
                 r.get("D1_2field", ""), r.get("D2_2field", "")]
            out[g] = [two_field(x) for x in q]
    return out


def read_calls_quartet(path: Path) -> Optional[List[str]]:
    if not path.exists():
        return None
    rs: List[str] = []
    ds: List[str] = []
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            allele = r.get("allele") or r.get("allele_2field") or ""
            if not allele:
                continue
            side = r.get("assignment", "")
            if side == "R":
                rs.append(two_field(allele))
            elif side == "D":
                ds.append(two_field(allele))
    if len(rs) < 2 or len(ds) < 2:
        return None
    return rs[:2] + ds[:2]


def read_tf_counts(path: Path) -> List[Tuple[str, float]]:
    out: List[Tuple[str, float]] = []
    if not path.exists():
        return out
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out.append((r["allele_2field"], float(r["fraction"])))
            except (KeyError, ValueError):
                continue
    return out


def read_chi_r_fit(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return None
    try:
        return float(rows[0].get("chi_r_fit", ""))
    except ValueError:
        return None


# ---------------- residuals ----------------

def em_residual(quartet: Sequence[str], fractions: Dict[str, float],
                chi_r: float) -> float:
    expected: Dict[str, float] = {}
    for a in quartet[:2]:
        expected[a] = expected.get(a, 0.0) + chi_r / 2.0
    for a in quartet[2:4]:
        expected[a] = expected.get(a, 0.0) + (1.0 - chi_r) / 2.0
    keys = set(fractions) | set(expected)
    return sum(abs(fractions.get(a, 0.0) - expected.get(a, 0.0)) for a in keys)


def af_l1(quartet: Sequence[str], geno: Dict[str, np.ndarray],
          obs_af: np.ndarray, weight: np.ndarray, chi_r: float
          ) -> Tuple[Optional[float], int]:
    if any(a not in geno for a in quartet):
        return None, 0
    G = [geno[a] for a in quartet]
    valid = (G[0] != -1) & (G[1] != -1) & (G[2] != -1) & (G[3] != -1)
    n = int(valid.sum())
    if n < 30:
        return None, n
    pred = 0.5 * (G[0] + G[1]) * chi_r + 0.5 * (G[2] + G[3]) * (1.0 - chi_r)
    diff = np.abs(obs_af - pred)
    denom = float((weight * valid).sum()) + 1e-9
    return float((weight * diff * valid).sum() / denom), n


def af_disc(quartet: Sequence[str], cand_geno_stack: np.ndarray,
            geno: Dict[str, np.ndarray],
            obs_af: np.ndarray, weight: np.ndarray, chi_r: float
            ) -> Tuple[Optional[float], int]:
    """AF L1 restricted to 'discriminating' sites (sites where the
    candidate pool's genotypes are not all identical or all unknown).

    cand_geno_stack: shape (n_cands, n_sites) int8 in {-1,0,1}.
    """
    if any(a not in geno for a in quartet):
        return None, 0
    G = [geno[a] for a in quartet]
    valid = (G[0] != -1) & (G[1] != -1) & (G[2] != -1) & (G[3] != -1)
    # discriminating: set of distinct values in {0,1} across candidates >= 2
    rows = cand_geno_stack
    has0 = (rows == 0).any(axis=0)
    has1 = (rows == 1).any(axis=0)
    disc = has0 & has1
    use = valid & disc
    n = int(use.sum())
    if n < 10:
        return None, n
    pred = 0.5 * (G[0] + G[1]) * chi_r + 0.5 * (G[2] + G[3]) * (1.0 - chi_r)
    diff = np.abs(obs_af - pred)
    denom = float((weight * use).sum()) + 1e-9
    return float((weight * diff * use).sum() / denom), n


# ---------------- read-phasing channel ----------------

def collect_phase_evidence(bam_path: Path, contig: str, sites_full: Sequence[Tuple[int, str, str]],
                           max_sites: int = 60, max_pairs: int = 400
                           ) -> Tuple[List[Tuple[int, str, str]], Dict[Tuple[int, int], np.ndarray]]:
    """For pairs of biallelic sites observed on the SAME read/fragment,
    return a co-occurrence count matrix per pair.

    Returns:
        sites_used: subset of input sites actually polled (capped at max_sites)
        pair_counts: {(i, j): np.array shape(2,2)} where entry [a,b] is
                     number of fragments carrying allele a (0=ref,1=alt)
                     at site i AND allele b at site j. i<j into sites_used.
    """
    if not bam_path.exists():
        return [], {}
    if not sites_full:
        return [], {}
    sites_used = sites_full[:max_sites]
    site_index = {(p, r, a): i for i, (p, r, a) in enumerate(sites_used)}
    n = len(sites_used)
    if n < 2:
        return sites_used, {}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    # per fragment: name -> dict {site_idx: 0|1}
    frag: Dict[str, Dict[int, int]] = defaultdict(dict)
    minp = min(p for p, _, _ in sites_used)
    maxp = max(p for p, _, _ in sites_used)
    try:
        for read in bam.fetch(contig, minp, maxp + 1):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if read.is_duplicate:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos is None or qpos is None:
                    continue
                if rpos < minp or rpos > maxp:
                    continue
                # Look up if this rpos matches any site
                # We allow only SNP sites (single-base ref/alt).
                # Build a small per-rpos lookup once outside loop.
                pass
    except (ValueError, OSError):
        bam.close()
        return sites_used, {}
    # Rebuild with a faster pre-indexed approach.
    pos_to_sites: Dict[int, List[Tuple[int, str, str]]] = defaultdict(list)
    for i, (p, r, a) in enumerate(sites_used):
        if len(r) == 1 and len(a) == 1:
            pos_to_sites[p].append((i, r, a))
    if not pos_to_sites:
        bam.close()
        return sites_used, {}
    try:
        for read in bam.fetch(contig, minp, maxp + 1):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            if read.is_duplicate:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            qname = read.query_name
            obs = frag.get(qname)
            if obs is None:
                obs = {}
                frag[qname] = obs
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos is None or qpos is None:
                    continue
                if rpos not in pos_to_sites:
                    continue
                base = seq[qpos].upper()
                for i, ref, alt in pos_to_sites[rpos]:
                    if base == ref:
                        prev = obs.get(i)
                        if prev is None:
                            obs[i] = 0
                        elif prev != 0:
                            obs[i] = -2  # conflict on this fragment
                    elif base == alt:
                        prev = obs.get(i)
                        if prev is None:
                            obs[i] = 1
                        elif prev != 1:
                            obs[i] = -2
    finally:
        bam.close()
    pair_counts: Dict[Tuple[int, int], np.ndarray] = {}
    for obs in frag.values():
        keys = [k for k, v in obs.items() if v in (0, 1)]
        if len(keys) < 2:
            continue
        keys.sort()
        for a, b in itertools.combinations(keys, 2):
            mat = pair_counts.get((a, b))
            if mat is None:
                mat = np.zeros((2, 2), dtype=np.int32)
                pair_counts[(a, b)] = mat
            mat[obs[a], obs[b]] += 1
    if len(pair_counts) > max_pairs:
        # keep most-populated pairs only
        items = sorted(pair_counts.items(), key=lambda kv: -int(kv[1].sum()))
        pair_counts = dict(items[:max_pairs])
    return sites_used, pair_counts


def phase_score(quartet: Sequence[str], sites_used: Sequence[Tuple[int, str, str]],
                pair_counts: Dict[Tuple[int, int], np.ndarray],
                geno: Dict[str, np.ndarray], chi_r: float
                ) -> Tuple[Optional[float], int]:
    """Compute a likelihood-style score for the quartet against observed
    fragment pair counts.

    For each pair of sites (i, j), the predicted joint distribution
    P(g_i, g_j) under the 4-hap model with chi_R is:
        sum over the 4 haplotypes h of weight(h) * 1[h_i = g_i, h_j = g_j]
    where weight = chi_R/2 for the two R haps and (1-chi_R)/2 for the two
    D haps.
    Score = - sum over pairs of sum_{a,b} c_{a,b} * log P(a,b)
    Lower is better. Pairs with any candidate genotype unknown are
    skipped.
    """
    if any(a not in geno for a in quartet):
        return None, 0
    G = [geno[a] for a in quartet]
    weights = np.array([chi_r / 2.0, chi_r / 2.0,
                        (1.0 - chi_r) / 2.0, (1.0 - chi_r) / 2.0])
    nsites = len(sites_used)
    site_idx_map = {i: i for i in range(nsites)}
    total_ll = 0.0
    n_pairs = 0
    for (i, j), mat in pair_counts.items():
        gi = np.array([G[h][i] for h in range(4)])
        gj = np.array([G[h][j] for h in range(4)])
        if (gi == -1).any() or (gj == -1).any():
            continue
        # build joint pred 2x2
        pred = np.zeros((2, 2))
        for h in range(4):
            pred[gi[h], gj[h]] += weights[h]
        # smooth to avoid log(0) — fixed prior, not tuned
        pred = pred + 1e-3
        pred = pred / pred.sum()
        ll = -float((mat * np.log(pred)).sum())
        total_ll += ll
        n_pairs += 1
    if n_pairs < 5:
        return None, n_pairs
    return total_ll / max(1.0, float(sum(int(m.sum()) for m in pair_counts.values()))), n_pairs


# ---------------- candidate / search ----------------

def build_candidates_C0(final_q, baseline_q, tf_top, top_k):
    seen: List[str] = []
    def add(a):
        if a and a not in seen:
            seen.append(a)
    for a in final_q:
        add(a)
    if baseline_q:
        for a in baseline_q:
            add(a)
    for a, _f in tf_top[:top_k]:
        add(a)
    return seen


def build_candidates_C1(c0: Sequence[str], gmap) -> List[str]:
    """Collapse by G-group: keep only one 2-field per G-group (the first
    one in order — order is final, then baseline, then EM top), so the
    representative is biased toward the final/baseline call. This is
    gene-agnostic (uses WMDA hla_nom_g)."""
    seen_g: Dict[str, str] = {}
    out: List[str] = []
    for a in c0:
        key = to_g_group(a, gmap)
        if key in seen_g:
            continue
        seen_g[key] = a
        out.append(a)
    return out


def quartet_key(q: Sequence[str]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    return (tuple(sorted(q[:2])), tuple(sorted(q[2:4])))


def enum_full(cands: Sequence[str]) -> List[Tuple[str, str, str, str]]:
    out = []
    pairs = []
    for i in range(len(cands)):
        for j in range(i, len(cands)):
            pairs.append((cands[i], cands[j]))
    for r1, r2 in pairs:
        for d1, d2 in pairs:
            out.append((r1, r2, d1, d2))
    return out


def enum_swap(cur: Sequence[str]) -> List[Tuple[str, str, str, str]]:
    """R<->D side swaps only."""
    R1, R2, D1, D2 = cur
    out = []
    out.append((D1, R2, R1, D2))
    out.append((D2, R2, D1, R1))
    out.append((R1, D1, R2, D2))
    out.append((R1, D2, D1, R2))
    out.append((D1, D2, R1, R2))  # full side swap
    return out


def enum_one(cur: Sequence[str], cands: Sequence[str]
             ) -> List[Tuple[str, str, str, str]]:
    """Replace exactly one position with each other candidate."""
    out = []
    for i in range(4):
        for c in cands:
            if c == cur[i]:
                continue
            q = list(cur)
            q[i] = c
            out.append(tuple(q))
    return out


# ---------------- per-locus engine ----------------

class Locus:
    __slots__ = ("sample", "gene", "set_label", "cur", "fractions", "chi_r",
                 "c0", "c1", "geno_c0", "geno_c1", "obs_af", "weight",
                 "sites", "pair_counts", "stack_c0", "stack_c1")

    def __init__(self, sample, gene, set_label, cur, fractions, chi_r,
                 c0, c1, geno_c0, geno_c1, obs_af, weight, sites, pair_counts,
                 stack_c0, stack_c1):
        self.sample = sample; self.gene = gene; self.set_label = set_label
        self.cur = cur; self.fractions = fractions; self.chi_r = chi_r
        self.c0 = c0; self.c1 = c1
        self.geno_c0 = geno_c0; self.geno_c1 = geno_c1
        self.obs_af = obs_af; self.weight = weight
        self.sites = sites; self.pair_counts = pair_counts
        self.stack_c0 = stack_c0; self.stack_c1 = stack_c1


def build_locus(sample, gene, set_label, asm_root, out_root, ref_path, bed_path,
                imgt_path, gmap, top_k):
    final_path = asm_root / sample / f"{sample}.final_calls.tsv"
    finals = read_final(final_path)
    if gene not in finals:
        return None
    cur = finals[gene]
    em_dir = out_root / sample / "em_refine"
    summary = em_dir / f"{gene}.summary.tsv"
    chi_r = read_chi_r_fit(summary)
    if chi_r is None:
        return None
    tf_rows = read_tf_counts(em_dir / f"{gene}.tf_counts.tsv")
    fractions = {a: f for a, f in tf_rows}

    gene_lc = gene.lower()
    base_path = asm_root / sample / gene_lc / gene / "calls.baseline.tsv"
    if not base_path.exists():
        base_path = asm_root / sample / gene_lc / gene / "calls.tsv"
    baseline_q = read_calls_quartet(base_path)

    c0 = build_candidates_C0(cur, baseline_q, tf_rows, top_k)
    if len(c0) < 2:
        return None
    c1 = build_candidates_C1(c0, gmap)

    contig = gene.replace("-", "_")
    s, e = parse_bed(str(bed_path), contig)
    if s is None:
        return None
    pc_vcf = out_root / sample / f"{sample}.pooled_continuous.vcf.gz"
    if not pc_vcf.exists():
        return None
    obs = collect_obs_af(str(pc_vcf), contig, (s, e))
    sites = sorted(obs.keys())
    if len(sites) < 30:
        return None
    obs_af = np.asarray([obs[k][0] for k in sites], dtype=float)
    obs_dp = np.asarray([obs[k][1] for k in sites], dtype=float)
    weight = np.minimum(obs_dp, 200.0)

    ref_seq = pysam.FastaFile(str(ref_path)).fetch(contig, s, e).upper()
    prefix = gene.replace("HLA-", "") + "*"
    imgt = load_imgt(str(imgt_path), prefix)
    keep = set(c0)
    imgt_sub = {n: t for n, t in imgt.items() if two_field(n) in keep}
    if not imgt_sub:
        return None
    geno_2field = imgt_genotypes_at_sites(imgt_sub, ref_seq, s, sites)
    geno_c0 = {a: geno_2field[a] for a in c0 if a in geno_2field}
    geno_c1 = {a: geno_2field[a] for a in c1 if a in geno_2field}
    if any(a not in geno_c0 for a in cur):
        return None
    stack_c0 = np.stack([geno_c0[a] for a in c0 if a in geno_c0])
    stack_c1 = np.stack([geno_c1[a] for a in c1 if a in geno_c1]) if geno_c1 else stack_c0

    bam_path = out_root / sample / f"{sample}.merge.bam"
    sites_used, pair_counts = collect_phase_evidence(bam_path, contig, sites)

    return Locus(sample, gene, set_label, cur, fractions, chi_r,
                 c0, c1, geno_c0, geno_c1, obs_af, weight,
                 sites_used, pair_counts, stack_c0, stack_c1)


# ---------------- strategies ----------------

# Each strategy defined by:
#   name, candidates_field, search_kind, channels (list), eps_per_channel
# Decision: alt must beat current on EVERY listed channel by its eps
# (intersection of margin conditions).

STRATEGIES = [
    # baseline two-channel from previous run, recorded for reference
    dict(name="C0_full_EM_AF",   cset="c0", search="full", channels=["em", "af"]),
    # restricted search spaces
    dict(name="C0_swap_EM_AF",   cset="c0", search="swap", channels=["em", "af"]),
    dict(name="C0_one_EM_AF",    cset="c0", search="one",  channels=["em", "af"]),
    # G-group collapsed candidate pool (gene-agnostic)
    dict(name="C1_full_EM_AF",   cset="c1", search="full", channels=["em", "af"]),
    dict(name="C1_swap_EM_AF",   cset="c1", search="swap", channels=["em", "af"]),
    # discriminating-site AF instead of plain AF
    dict(name="C0_full_EM_DISC", cset="c0", search="full", channels=["em", "disc"]),
    dict(name="C0_swap_EM_DISC", cset="c0", search="swap", channels=["em", "disc"]),
    # phasing-augmented (restricted to swap/one to keep runtime sane)
    dict(name="C0_swap_EM_PHASE",         cset="c0", search="swap", channels=["em", "phase"]),
    dict(name="C0_one_EM_PHASE",          cset="c0", search="one",  channels=["em", "phase"]),
    dict(name="C0_swap_AF_PHASE",         cset="c0", search="swap", channels=["af", "phase"]),
    dict(name="C0_swap_EM_AF_PHASE",      cset="c0", search="swap", channels=["em", "af", "phase"]),
    dict(name="C0_one_EM_AF_PHASE",       cset="c0", search="one",  channels=["em", "af", "phase"]),
    dict(name="C1_swap_EM_AF_PHASE",      cset="c1", search="swap", channels=["em", "af", "phase"]),
]

EPS = dict(em=0.01, af=0.005, disc=0.005, phase=0.0005)


_SCORE_CACHE: Dict[Tuple[int, str, Tuple, str], float] = {}


def score_channel(q, locus: Locus, cset_key: str, channel: str):
    key = (id(locus), cset_key, quartet_key(q), channel)
    if key in _SCORE_CACHE:
        return _SCORE_CACHE[key]
    geno = locus.geno_c0 if cset_key == "c0" else locus.geno_c1
    stack = locus.stack_c0 if cset_key == "c0" else locus.stack_c1
    if channel == "em":
        v = em_residual(q, locus.fractions, locus.chi_r)
    elif channel == "af":
        v, _ = af_l1(q, geno, locus.obs_af, locus.weight, locus.chi_r)
    elif channel == "disc":
        v, _ = af_disc(q, stack, geno, locus.obs_af, locus.weight, locus.chi_r)
    elif channel == "phase":
        v, _ = phase_score(q, locus.sites, locus.pair_counts, geno, locus.chi_r)
    else:
        v = None
    _SCORE_CACHE[key] = v
    return v


def score_quartet(q, locus: Locus, cset_key: str, channels: Sequence[str]):
    return {ch: score_channel(q, locus, cset_key, ch) for ch in channels}


def search_alts(strat, locus: Locus):
    cands = locus.c0 if strat["cset"] == "c0" else locus.c1
    if strat["search"] == "full":
        qs = enum_full(cands)
    elif strat["search"] == "swap":
        qs = enum_swap(locus.cur)
        # discard any swap that introduces a non-candidate (shouldn't happen)
        qs = [q for q in qs if all(a in cands for a in q)]
    elif strat["search"] == "one":
        qs = enum_one(locus.cur, cands)
    else:
        qs = []
    seen = {quartet_key(locus.cur)}
    out = []
    for q in qs:
        k = quartet_key(q)
        if k in seen:
            continue
        seen.add(k)
        out.append(q)
    return out


def evaluate_strategy(strat, locus: Locus):
    cur_score = score_quartet(locus.cur, locus, strat["cset"], strat["channels"])
    # Reject locus for this strategy if any required current channel is None.
    for ch in strat["channels"]:
        if cur_score.get(ch) is None:
            return None
    best = None
    best_total_gain = 0.0
    for q in search_alts(strat, locus):
        sc = score_quartet(q, locus, strat["cset"], strat["channels"])
        ok = True
        total_gain = 0.0
        for ch in strat["channels"]:
            v = sc.get(ch)
            if v is None:
                ok = False
                break
            d = cur_score[ch] - v
            if d < EPS[ch]:
                ok = False
                break
            total_gain += d
        if ok and total_gain > best_total_gain:
            best_total_gain = total_gain
            best = (q, sc)
    return dict(cur=cur_score, alt=best[1] if best else None,
                alt_q=list(best[0]) if best else None,
                gain=best_total_gain if best else 0.0)


# ---------------- truth scoring ----------------

def truth_score(quartet, truth_p, truth_d, gmap) -> int:
    p = sorted([quartet[0], quartet[1]])
    d = sorted([quartet[2], quartet[3]])
    tp = normalize_for_display(truth_p, "2field", gmap)
    td = normalize_for_display(truth_d, "2field", gmap)
    return overlap(tp, p) + overlap(td, d)


# ---------------- driver ----------------

def discover_samples(fq_root: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for fq1 in sorted(fq_root.glob("*/*_R1_001.fastq.gz")):
        sample = fq1.name.replace("_R1_001.fastq.gz", "")
        parent = fq1.parent.name
        label = parent.replace(" ", "-").replace("_", "-").lower()
        if label.startswith("set-") and len(label) >= 5:
            out[sample] = label
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--out-root", required=True, type=Path)
    ap.add_argument("--fq-root", required=True, type=Path)
    ap.add_argument("--truth-dir", required=True, type=Path)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--imgt", required=True, type=Path)
    ap.add_argument("--report-prefix", required=True, type=Path)
    ap.add_argument("--top-k", type=int, default=10)
    args = ap.parse_args()

    samples = discover_samples(args.fq_root)
    truths = {}
    for label in {v for v in samples.values()}:
        tp = args.truth_dir / f"truth_typing-{label}.tsv"
        if tp.exists():
            truths[label] = load_truth(tp)
    gmap_path = SCRIPT_DIR / "resources" / "spechla" / "db" / "HLA" / "hla_nom_g.txt"
    gmap = load_g_group(gmap_path) if gmap_path.exists() else {}

    rows: List[Dict[str, str]] = []
    cur_total = 0
    for sample in sorted(samples):
        label = samples[sample]
        truth = truths.get(label)
        for gene in GENES:
            _SCORE_CACHE.clear()
            locus = build_locus(sample, gene, label, args.asm_root, args.out_root,
                                args.ref, args.bed, args.imgt, gmap, args.top_k)
            if locus is None:
                continue
            tp = truth["PATIENT"][gene] if truth else []
            td = truth["DONOR"][gene] if truth else []
            cs = truth_score(locus.cur, tp, td, gmap) if truth else 0
            cur_total += cs
            print(f"[{sample}\t{gene}\tcur={','.join(locus.cur)}\tcur_score={cs}/4]",
                  flush=True)
            for strat in STRATEGIES:
                r = evaluate_strategy(strat, locus)
                if r is None:
                    rows.append(dict(strategy=strat["name"], sample=sample,
                                     gene=gene, set=label,
                                     proposed="", verdict="skip",
                                     cur_score=str(cs), alt_score="",
                                     delta_score="0"))
                    continue
                if r["alt_q"] is None:
                    rows.append(dict(strategy=strat["name"], sample=sample,
                                     gene=gene, set=label,
                                     proposed="", verdict="hold",
                                     cur_score=str(cs), alt_score="",
                                     delta_score="0"))
                    continue
                asc = truth_score(r["alt_q"], tp, td, gmap) if truth else 0
                delta = asc - cs
                if delta > 0:
                    verdict = "improve"
                elif delta < 0:
                    verdict = "regress"
                else:
                    verdict = "neutral"
                rows.append(dict(strategy=strat["name"], sample=sample,
                                 gene=gene, set=label,
                                 proposed=",".join(r["alt_q"]),
                                 verdict=verdict,
                                 cur_score=str(cs), alt_score=str(asc),
                                 delta_score=str(delta)))

    args.report_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_tsv = args.report_prefix.with_suffix(".tsv")
    cols = ["strategy", "set", "sample", "gene", "proposed", "verdict",
            "cur_score", "alt_score", "delta_score"]
    with out_tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in cols) + "\n")

    # Per-strategy summary.
    summary_path = args.report_prefix.with_suffix(".summary")
    by = defaultdict(lambda: dict(prop=0, imp=0, reg=0, neu=0, hold=0,
                                  skip=0, net=0,
                                  by_set=defaultdict(lambda: [0, 0, 0, 0]),
                                  by_gene=defaultdict(lambda: [0, 0, 0, 0])))
    for r in rows:
        s = by[r["strategy"]]
        v = r["verdict"]
        if v == "skip":
            s["skip"] += 1
            continue
        if v == "hold":
            s["hold"] += 1
            continue
        s["prop"] += 1
        if v == "improve":
            s["imp"] += 1
        elif v == "regress":
            s["reg"] += 1
        elif v == "neutral":
            s["neu"] += 1
        try:
            s["net"] += int(r["delta_score"])
        except ValueError:
            pass
        sset = s["by_set"][r["set"]]
        sgen = s["by_gene"][r["gene"]]
        for arr in (sset, sgen):
            if v == "improve":
                arr[0] += 1
            elif v == "regress":
                arr[1] += 1
            elif v == "neutral":
                arr[2] += 1
            try:
                arr[3] += int(r["delta_score"])
            except ValueError:
                pass

    with summary_path.open("w") as fh:
        fh.write(f"# baseline 2field score (no rerank) = {cur_total}/360\n")
        fh.write("# strategy\tprop\timp\treg\tneu\thold\tskip\tnet\tafter\n")
        for name in [s["name"] for s in STRATEGIES]:
            s = by[name]
            after = cur_total + s["net"]
            fh.write(f"{name}\t{s['prop']}\t{s['imp']}\t{s['reg']}\t{s['neu']}"
                     f"\t{s['hold']}\t{s['skip']}\t{s['net']}\t{after}/360\n")
        fh.write("\n# per-strategy by_set imp/reg/neu/net\n")
        for name in [s["name"] for s in STRATEGIES]:
            s = by[name]
            line = name
            for k in sorted(s["by_set"]):
                v = s["by_set"][k]
                line += f"\t{k}={v[0]}/{v[1]}/{v[2]}/{v[3]}"
            fh.write(line + "\n")
        fh.write("\n# per-strategy by_gene imp/reg/neu/net\n")
        for name in [s["name"] for s in STRATEGIES]:
            s = by[name]
            line = name
            for k in sorted(s["by_gene"]):
                v = s["by_gene"][k]
                line += f"\t{k}={v[0]}/{v[1]}/{v[2]}/{v[3]}"
            fh.write(line + "\n")

    print()
    print(f"# baseline 2field score (no rerank) = {cur_total}/360")
    print("# strategy\tprop\timp\treg\tneu\thold\tnet\tafter")
    for name in [s["name"] for s in STRATEGIES]:
        s = by[name]
        after = cur_total + s["net"]
        print(f"{name}\t{s['prop']}\t{s['imp']}\t{s['reg']}\t{s['neu']}"
              f"\t{s['hold']}\t{s['net']}\t{after}/360")


if __name__ == "__main__":
    main()
