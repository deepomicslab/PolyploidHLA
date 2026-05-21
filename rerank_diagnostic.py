#!/usr/bin/env python3
"""Truth-free dual-evidence rerank diagnostic.

For every (sample, gene) in the current rescue run this tool:

  1. Builds a *restricted* candidate allele set for the gene from
       final_calls 2-field alleles (R/D)
       baseline calls.tsv 2-field alleles (R/D, if available)
       EM tf_counts top-K alleles
     No new alleles are introduced from outside.

  2. Enumerates all (R1<=R2, D1<=D2) quartets over that small set.

  3. Scores every quartet with two **orthogonal** channels:
       EM channel: L1 residual of expected vs observed 2-field fractions
                   under the per-gene fitted chi_R (same metric used by
                   fit_4hap / em_refine_gate.quartet_residual).
       AF channel: weighted L1 residual of expected vs observed allele
                   frequency on pooled-continuous biallelic sites
                   (caller_free_4hap mechanism, independent of read EM).

  4. Proposes an override only when an alternative quartet improves
     BOTH channels relative to the current final, by at least a fixed
     a-priori margin (no tuning on truth).

  5. Compares the proposal against truth ONLY for post-hoc evaluation.
     Truth is never consulted to decide overrides.

The script writes a TSV report and prints a summary stratified by gene
and by source set (set-a / set-b / set-c) so that any apparent gain or
regression can be inspected in a leave-one-set-out spirit.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import os
import re
import sys
from collections import Counter, defaultdict
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
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
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


def af_residual(quartet: Sequence[str], geno: Dict[str, np.ndarray],
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


# ---------------- candidate set / enumeration ----------------

def build_candidates(final_q: Sequence[str], baseline_q: Optional[Sequence[str]],
                     tf_top: Sequence[Tuple[str, float]], top_k: int
                     ) -> List[str]:
    seen: List[str] = []
    def add(a: str) -> None:
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


def enumerate_quartets(cands: Sequence[str]) -> List[Tuple[str, str, str, str]]:
    out: List[Tuple[str, str, str, str]] = []
    pairs = []
    for i in range(len(cands)):
        for j in range(i, len(cands)):
            pairs.append((cands[i], cands[j]))
    for r1, r2 in pairs:
        for d1, d2 in pairs:
            out.append((r1, r2, d1, d2))
    return out


def quartet_key(q: Sequence[str]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    return (tuple(sorted(q[:2])), tuple(sorted(q[2:4])))


# ---------------- per-(sample, gene) processing ----------------

def process_locus(sample: str, gene: str, asm_root: Path, out_root: Path,
                  ref_path: Path, bed_path: Path, imgt_path: Path,
                  top_k: int, eps_em: float, eps_af: float
                  ) -> Optional[Dict[str, str]]:
    final_path = asm_root / sample / f"{sample}.final_calls.tsv"
    finals = read_final(final_path)
    if gene not in finals:
        return None
    cur_q = finals[gene]

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

    cands = build_candidates(cur_q, baseline_q, tf_rows, top_k)
    if len(cands) < 2:
        return None

    contig = gene.replace("-", "_")
    bed_se = parse_bed(str(bed_path), contig)
    if bed_se[0] is None:
        return None
    s, e = bed_se
    pc_vcf = out_root / sample / f"{sample}.pooled_continuous.vcf.gz"
    if not pc_vcf.exists():
        return None
    obs = collect_obs_af(str(pc_vcf), contig, (s, e))
    sites = sorted(obs.keys())
    if len(sites) < 30:
        af_ok = False
    else:
        af_ok = True
    obs_af = np.asarray([obs[k][0] for k in sites], dtype=float)
    obs_dp = np.asarray([obs[k][1] for k in sites], dtype=float)
    weight = np.minimum(obs_dp, 200.0)

    geno: Dict[str, np.ndarray] = {}
    if af_ok:
        ref_seq = pysam.FastaFile(str(ref_path)).fetch(contig, s, e).upper()
        prefix = gene.replace("HLA-", "") + "*"
        imgt = load_imgt(str(imgt_path), prefix)
        # Restrict the IMGT representative pool to the candidate 2-field
        # families to keep alignment cost bounded and predictable.
        keep = set(cands)
        imgt_sub = {n: s for n, s in imgt.items() if two_field(n) in keep}
        if not imgt_sub:
            af_ok = False
        else:
            geno_2field = imgt_genotypes_at_sites(imgt_sub, ref_seq, s, sites)
            geno = {a: geno_2field[a] for a in cands if a in geno_2field}

    quartets = enumerate_quartets(cands)
    em_scores: Dict[Tuple, float] = {}
    af_scores: Dict[Tuple, Tuple[Optional[float], int]] = {}
    for q in quartets:
        k = quartet_key(q)
        if k not in em_scores:
            em_scores[k] = em_residual(q, fractions, chi_r)
        if af_ok and k not in af_scores:
            af_scores[k] = af_residual(q, geno, obs_af, weight, chi_r)
        elif not af_ok:
            af_scores[k] = (None, 0)

    cur_key = quartet_key(cur_q)
    cur_em = em_scores.get(cur_key, em_residual(cur_q, fractions, chi_r))
    cur_af, cur_af_n = af_scores.get(cur_key, (None, 0))
    if af_ok and cur_key not in af_scores:
        cur_af, cur_af_n = af_residual(cur_q, geno, obs_af, weight, chi_r)

    best_alt = None
    best_gain = 0.0
    if af_ok and cur_af is not None:
        for k, em_s in em_scores.items():
            if k == cur_key:
                continue
            af_pair = af_scores.get(k, (None, 0))
            if af_pair[0] is None:
                continue
            d_em = cur_em - em_s
            d_af = cur_af - af_pair[0]
            if d_em >= eps_em and d_af >= eps_af:
                gain = d_em + d_af  # symmetric joint margin
                if gain > best_gain:
                    best_gain = gain
                    best_alt = (k, em_s, af_pair[0], af_pair[1])

    proposed_q: Optional[List[str]] = None
    if best_alt is not None:
        rp, dp = best_alt[0]
        proposed_q = list(rp) + list(dp)

    return {
        "sample": sample,
        "gene": gene,
        "chi_r": f"{chi_r:.4f}",
        "n_cands": str(len(cands)),
        "n_sites": str(len(sites)),
        "af_ok": "1" if af_ok else "0",
        "current": ",".join(cur_q),
        "cur_em": f"{cur_em:.5f}",
        "cur_af": "" if cur_af is None else f"{cur_af:.5f}",
        "cur_af_n": str(cur_af_n),
        "proposed": ",".join(proposed_q) if proposed_q else "",
        "alt_em": f"{best_alt[1]:.5f}" if best_alt else "",
        "alt_af": f"{best_alt[2]:.5f}" if best_alt else "",
        "alt_af_n": str(best_alt[3]) if best_alt else "",
        "joint_gain": f"{best_gain:.5f}" if best_alt else "",
    }


# ---------------- truth-side evaluation ----------------

def truth_score(quartet: Sequence[str], truth_p: Sequence[str],
                truth_d: Sequence[str], gmap) -> int:
    p = sorted([quartet[0], quartet[1]])
    d = sorted([quartet[2], quartet[3]])
    tp = normalize_for_display(truth_p, "2field", gmap)
    td = normalize_for_display(truth_d, "2field", gmap)
    return overlap(tp, p) + overlap(td, d)


# ---------------- driver ----------------

def discover_samples(fq_root: Path) -> Dict[str, str]:
    """Return {sample_name: set_label} where set_label is 'set-a/b/c'."""
    out: Dict[str, str] = {}
    for fq1 in sorted(fq_root.glob("*/*_R1_001.fastq.gz")):
        sample = fq1.name.replace("_R1_001.fastq.gz", "")
        parent = fq1.parent.name
        label = parent.replace(" ", "-").replace("_", "-").lower()
        if label.startswith("set-") and len(label) >= 5:
            out[sample] = label
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--out-root", required=True, type=Path)
    ap.add_argument("--fq-root", required=True, type=Path)
    ap.add_argument("--truth-dir", required=True, type=Path)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--imgt", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--top-k", type=int, default=10,
                    help="EM top-K alleles added to candidate set")
    ap.add_argument("--eps-em", type=float, default=0.01,
                    help="minimum absolute EM L1 residual reduction "
                         "required to consider an override (generic, "
                         "not tuned on truth)")
    ap.add_argument("--eps-af", type=float, default=0.005,
                    help="minimum absolute AF L1 residual reduction "
                         "required to consider an override (generic, "
                         "not tuned on truth)")
    ap.add_argument("--genes", nargs="+", default=GENES)
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
    for sample in sorted(samples):
        label = samples[sample]
        truth = truths.get(label)
        for gene in args.genes:
            r = process_locus(sample, gene, args.asm_root, args.out_root,
                              args.ref, args.bed, args.imgt,
                              args.top_k, args.eps_em, args.eps_af)
            if r is None:
                continue
            r["set"] = label
            r["truth_p"] = ""
            r["truth_d"] = ""
            r["cur_score"] = ""
            r["alt_score"] = ""
            r["delta_score"] = ""
            r["verdict"] = ""
            if truth is not None:
                tp = truth["PATIENT"][gene]
                td = truth["DONOR"][gene]
                cur_q = r["current"].split(",")
                cs = truth_score(cur_q, tp, td, gmap)
                r["truth_p"] = ",".join(normalize_for_display(tp, "2field", gmap))
                r["truth_d"] = ",".join(normalize_for_display(td, "2field", gmap))
                r["cur_score"] = str(cs)
                if r["proposed"]:
                    pq = r["proposed"].split(",")
                    asc = truth_score(pq, tp, td, gmap)
                    r["alt_score"] = str(asc)
                    delta = asc - cs
                    r["delta_score"] = str(delta)
                    if delta > 0:
                        r["verdict"] = "improve"
                    elif delta < 0:
                        r["verdict"] = "regress"
                    else:
                        r["verdict"] = "neutral"
            rows.append(r)
            print(f"[{sample}\t{gene}\tcur={r['current']}"
                  f"\tcur_em={r['cur_em']}\tcur_af={r['cur_af']}"
                  f"\tprop={r['proposed']}\tverdict={r['verdict']}",
                  flush=True)

    cols = ["set", "sample", "gene", "chi_r", "n_cands", "n_sites", "af_ok",
            "current", "cur_em", "cur_af", "cur_af_n",
            "proposed", "alt_em", "alt_af", "alt_af_n", "joint_gain",
            "truth_p", "truth_d", "cur_score", "alt_score",
            "delta_score", "verdict"]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in cols) + "\n")

    proposed = [r for r in rows if r["proposed"]]
    improve = [r for r in proposed if r["verdict"] == "improve"]
    regress = [r for r in proposed if r["verdict"] == "regress"]
    neutral = [r for r in proposed if r["verdict"] == "neutral"]
    net = sum(int(r["delta_score"]) for r in proposed if r["delta_score"])
    cur_total = sum(int(r["cur_score"]) for r in rows if r["cur_score"])
    print()
    print(f"# samples: {len(samples)}; loci scored: {len(rows)}; "
          f"proposed overrides: {len(proposed)}")
    print(f"# improve={len(improve)} regress={len(regress)} "
          f"neutral={len(neutral)}  net_delta_score={net}")
    print(f"# baseline 2field score (no rerank): {cur_total}/360")
    print(f"# rerank-applied 2field score:     {cur_total + net}/360")
    by_set = defaultdict(lambda: [0, 0, 0, 0])  # imp, reg, neu, net
    by_gene = defaultdict(lambda: [0, 0, 0, 0])
    for r in proposed:
        for d, k in ((by_set, r["set"]), (by_gene, r["gene"])):
            v = d[k]
            if r["verdict"] == "improve":
                v[0] += 1
            elif r["verdict"] == "regress":
                v[1] += 1
            elif r["verdict"] == "neutral":
                v[2] += 1
            try:
                v[3] += int(r["delta_score"]) if r["delta_score"] else 0
            except ValueError:
                pass
    print("# overrides by set (imp/reg/neu/net_delta):")
    for k in sorted(by_set):
        v = by_set[k]
        print(f"#   {k}\t{v[0]}/{v[1]}/{v[2]}\tnet={v[3]}")
    print("# overrides by gene (imp/reg/neu/net_delta):")
    for k in sorted(by_gene):
        v = by_gene[k]
        print(f"#   {k}\t{v[0]}/{v[1]}/{v[2]}\tnet={v[3]}")


if __name__ == "__main__":
    main()
