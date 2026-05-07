#!/usr/bin/env python
"""Estimate per-gene and global chi_R from a freebayes --pooled-continuous VCF.

Pooled-continuous mode reports raw observed AF per site (no GT). For a true
chimeric tetraploid signal, AF clusters at:
  {chi_R/2, chi_D/2, 1-chi_D/2, 1-chi_R/2}                    (dose 1 and 3)
  {chi_R, chi_D, 1-chi_D, 1-chi_R}                            (dose 2 = unique to one side)
  ...all symmetric around 0.5.

Strategy: fold AF -> min(AF, 1-AF) to [0, 0.5]. Then look at the LEFT mode
(<0.25) -> represents the smaller of (chi_R/2, chi_D/2), i.e. chi_R/2 if
chi_R < 0.5. Mode found via simple histogram peak after dropping noise floor.
"""
import sys, gzip, argparse
from collections import defaultdict
import numpy as np


def parse_vcf(path):
    """Return list of (chrom, pos, af, dp) from --pooled-continuous output."""
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, qual, flt, info, fmt = f[:9]
            sample = f[9]
            kv = dict(zip(fmt.split(":"), sample.split(":")))
            try:
                ro = int(kv.get("RO", "0"))
                ao = int(kv["AO"].split(",")[0])
            except (ValueError, KeyError):
                continue
            dp = ro + ao
            if dp < 30: continue
            if "," in alt: continue  # skip multiallelic
            af = ao / dp
            if af < 0.005 or af > 0.995: continue
            out.append((chrom, int(pos), af, dp))
    return out


def estimate_chi_from_af(afs, fold=True):
    """Pick the smallest mode. Use KDE-like histogram on folded AF."""
    if not afs: return None
    arr = np.asarray(afs)
    if fold:
        arr = np.where(arr > 0.5, 1 - arr, arr)
    # We want the small-amplitude mode, representing chi_R/2 (low-dose recipient)
    # or chi_R (dose-2 recipient unique). Bins:
    bins = np.arange(0.0, 0.5 + 1e-9, 0.01)  # 1% bins
    h, edges = np.histogram(arr, bins=bins)
    # smooth (3-bin moving average)
    hs = np.convolve(h, np.ones(3)/3.0, mode='same')
    # find peaks with min height = max/4 and away from 0
    peaks = []
    for i in range(2, len(hs)-1):
        if hs[i] >= hs[i-1] and hs[i] >= hs[i+1] and hs[i] >= max(hs)*0.15:
            center = (edges[i] + edges[i+1]) / 2
            if center >= 0.015:  # ignore <1.5% noise floor
                peaks.append((center, hs[i]))
    if not peaks:
        return None
    # smallest-AF peak corresponds to chi_R/2 (recipient low-dose)
    peaks.sort()
    chi_r_over2 = peaks[0][0]
    chi_r = 2 * chi_r_over2
    return {
        "chi_r": chi_r,
        "chi_r_over2_peak": chi_r_over2,
        "n_af": len(arr),
        "all_peaks": peaks[:6],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("--per-gene", action="store_true")
    args = ap.parse_args()
    rows = parse_vcf(args.vcf)
    print(f"# {len(rows)} biallelic AF rows after filters", file=sys.stderr)
    res = estimate_chi_from_af([r[2] for r in rows])
    if res is None:
        print(f"GLOBAL  chi_R=NA  n={len(rows)}  peaks=[]")
        return
    print(f"GLOBAL  chi_R={res['chi_r']:.4f}  n={res['n_af']}  peaks={res['all_peaks']}")
    if args.per_gene:
        per = defaultdict(list)
        for chrom, pos, af, dp in rows:
            per[chrom].append(af)
        for chrom, afs in sorted(per.items()):
            r = estimate_chi_from_af(afs)
            if r is None:
                print(f"  {chrom:14s} n={len(afs):4d}  (no peak)")
            else:
                print(f"  {chrom:14s} n={len(afs):4d}  chi_R={r['chi_r']:.4f}  "
                      f"peaks={r['all_peaks'][:3]}")


if __name__ == "__main__":
    main()
