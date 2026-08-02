#!/usr/bin/env python
"""Estimate per-gene and global chi_R from a freebayes --pooled-continuous VCF.

Pooled-continuous mode reports raw observed AF per site (no GT). For a true
chimeric tetraploid signal, AF clusters at:
  {chi_R/2, chi_D/2, 1-chi_D/2, 1-chi_R/2}                    (dose 1 and 3)
  {chi_R, chi_D, 1-chi_D, 1-chi_R}                            (dose 2 = unique to one side)
  ...all symmetric around 0.5.

Strategy: fold AF -> min(AF, 1-AF) to [0, 0.5]. Then look at the left mode,
which represents the smaller of (chi_R/2, chi_D/2). The production model is
validated for mixture fractions of at least 10%, so peaks below 4% folded AF
are treated as noise rather than eligible low-component modes.
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


def estimate_chi_from_af(
    afs, fold=True, min_chi=0.10, prior_chi=None, prior_max_af_rows=2000,
):
    """Estimate the low mixture component from the left folded-AF mode."""
    if afs is None or len(afs) == 0:
        return None
    arr = np.asarray(afs)
    if fold:
        arr = np.where(arr > 0.5, 1 - arr, arr)
    # We want the small-amplitude mode, representing chi_R/2 (low-dose recipient)
    # or chi_R (dose-2 recipient unique). Bins:
    bins = np.arange(0.0, 0.5 + 1e-9, 0.01)  # 1% bins
    h, edges = np.histogram(arr, bins=bins)
    # smooth (3-bin moving average)
    hs = np.convolve(h, np.ones(3)/3.0, mode='same')
    minimum_peak = max(0.015, min_chi / 2.0 - 0.01)
    peaks = []
    for i in range(2, len(hs)-1):
        if hs[i] >= hs[i-1] and hs[i] >= hs[i+1] and hs[i] >= max(hs)*0.15:
            center = (edges[i] + edges[i+1]) / 2
            if minimum_peak <= center <= 0.255:
                peaks.append((center, hs[i]))
    if not peaks:
        return None
    peaks.sort()
    chi_r_over2 = peaks[0][0]
    selection = "leftmost"
    if prior_chi is not None and len(arr) < prior_max_af_rows:
        folded_prior = min(float(prior_chi), 1.0 - float(prior_chi))
        leftmost_fraction = 2.0 * chi_r_over2
        chi_r_over2 = min(
            peaks,
            key=lambda peak: (
                abs(2.0 * peak[0] - folded_prior)
                + 0.5 * (2.0 * peak[0] - leftmost_fraction),
                peak[0],
            ),
        )[0]
        selection = "prior_guided"
    chi_r = 2 * chi_r_over2
    return {
        "chi_r": chi_r,
        "chi_r_over2_peak": chi_r_over2,
        "n_af": len(arr),
        "all_peaks": peaks[:6],
        "minimum_peak": minimum_peak,
        "selection": selection,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("--per-gene", action="store_true")
    ap.add_argument(
        "--min-chi", type=float, default=0.10,
        help="minimum supported mixture fraction (default: 0.10)",
    )
    ap.add_argument(
        "--prior-chi", type=float,
        help="optional GT-based chi used to choose among low-evidence AF peaks",
    )
    args = ap.parse_args()
    rows = parse_vcf(args.vcf)
    print(f"# {len(rows)} biallelic AF rows after filters", file=sys.stderr)
    res = estimate_chi_from_af(
        [r[2] for r in rows], min_chi=args.min_chi, prior_chi=args.prior_chi,
    )
    if res is None:
        print(f"GLOBAL  chi_R=NA  n={len(rows)}  peaks=[]")
        return
    print(f"GLOBAL  chi_R={res['chi_r']:.4f}  n={res['n_af']}  peaks={res['all_peaks']}")
    if args.per_gene:
        per = defaultdict(list)
        for chrom, pos, af, dp in rows:
            per[chrom].append(af)
        for chrom, afs in sorted(per.items()):
            r = estimate_chi_from_af(
                afs, min_chi=args.min_chi, prior_chi=args.prior_chi,
            )
            if r is None:
                print(f"  {chrom:14s} n={len(afs):4d}  (no peak)")
            else:
                print(f"  {chrom:14s} n={len(afs):4d}  chi_R={r['chi_r']:.4f}  "
                      f"peaks={r['all_peaks'][:3]}")


if __name__ == "__main__":
    main()
