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


QC_GENE_MIN_AF = 30
# Development-derived gates; freeze before independent validation.
MODEL_MISMATCH_GENE_DELTA = 0.15
MODEL_MISMATCH_PEAK_RATIO = 0.67
MODEL_MISMATCH_RESIDUAL_MARGIN = 0.02
LOW_CONFIDENCE_CI_WIDTH = 0.15
LOW_CONFIDENCE_BOOTSTRAP_FINITE = 0.90


def parse_vcf(path, include_contigs=None):
    """Return list of (chrom, pos, af, dp) from --pooled-continuous output."""
    include_contigs = set(include_contigs or [])
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, qual, flt, info, fmt = f[:9]
            if include_contigs and chrom not in include_contigs:
                continue
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
    dps=None,
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
    selected_peak = min(peaks, key=lambda peak: abs(peak[0] - chi_r_over2))
    alternatives = [peak for peak in peaks if peak != selected_peak]
    best_alternative_height = max((peak[1] for peak in alternatives), default=None)
    peak_ratio = (
        selected_peak[1] / best_alternative_height
        if best_alternative_height not in (None, 0) else None
    )

    weights = np.asarray(dps, dtype=float) if dps is not None else np.ones(len(arr))

    def weighted_residual(candidate):
        grid = np.asarray([
            candidate / 2.0, candidate, (1.0 - candidate) / 2.0, 0.5,
        ])
        distances = np.min(np.abs(arr[:, None] - grid[None, :]), axis=1)
        return float(np.average(distances, weights=weights))

    residual = weighted_residual(chi_r)
    alternative_residuals = [weighted_residual(2.0 * peak[0]) for peak in alternatives]
    residual_margin = (
        residual - min(alternative_residuals) if alternative_residuals else None
    )
    return {
        "chi_r": chi_r,
        "chi_r_over2_peak": chi_r_over2,
        "n_af": len(arr),
        "all_peaks": peaks[:6],
        "minimum_peak": minimum_peak,
        "selection": selection,
        "selected_peak_height": float(selected_peak[1]),
        "peak_ratio": peak_ratio,
        "weighted_residual": residual,
        "residual_margin": residual_margin,
    }


def bootstrap_chi(afs, dps, min_chi, prior_chi, replicates, seed):
    """Return a deterministic row-bootstrap interval and empirical peak odds."""
    if replicates <= 0:
        return None
    afs = np.asarray(afs, dtype=float)
    dps = np.asarray(dps, dtype=float)
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(replicates):
        indices = rng.integers(0, len(afs), size=len(afs))
        result = estimate_chi_from_af(
            afs[indices], min_chi=min_chi, prior_chi=prior_chi,
            dps=dps[indices],
        )
        estimates.append(np.nan if result is None else result["chi_r"])
    values = np.asarray(estimates, dtype=float)
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return {
            "ci_low": None, "ci_high": None, "finite_fraction": 0.0,
            "peak_support": 0.0, "peak_odds": 0.0,
        }
    rounded = np.round(finite / 0.02) * 0.02
    _, counts = np.unique(rounded, return_counts=True)
    counts = np.sort(counts)[::-1]
    best = int(counts[0])
    second = int(counts[1]) if len(counts) > 1 else 0
    return {
        "ci_low": float(np.quantile(finite, 0.025)),
        "ci_high": float(np.quantile(finite, 0.975)),
        "finite_fraction": len(finite) / replicates,
        "peak_support": best / len(finite),
        "peak_odds": (best + 0.5) / (second + 0.5),
    }


def fmt(value, digits=4):
    return "NA" if value is None else f"{value:.{digits}f}"


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
    ap.add_argument(
        "--bootstrap", type=int, default=200,
        help="row-bootstrap replicates for the chi interval (default: 200)",
    )
    ap.add_argument("--bootstrap-seed", type=int, default=20260802)
    ap.add_argument(
        "--include-contigs", nargs="+",
        help="restrict the global and per-gene estimates to these contigs",
    )
    args = ap.parse_args()
    rows = parse_vcf(args.vcf, args.include_contigs)
    print(f"# {len(rows)} biallelic AF rows after filters", file=sys.stderr)
    res = estimate_chi_from_af(
        [r[2] for r in rows], min_chi=args.min_chi, prior_chi=args.prior_chi,
        dps=[r[3] for r in rows],
    )
    if res is None:
        print(f"GLOBAL  chi_R=NA  raw_chi_R=NA  status=LOW_CONFIDENCE  "
              f"reasons=no_eligible_peak  n_af={len(rows)}  peaks=[]")
        return

    per = defaultdict(list)
    for chrom, pos, af, dp in rows:
        per[chrom].append((af, dp))
    gene_results = {}
    for chrom, observations in sorted(per.items()):
        if len(observations) < QC_GENE_MIN_AF:
            continue
        gene_result = estimate_chi_from_af(
            [item[0] for item in observations], min_chi=args.min_chi,
            prior_chi=args.prior_chi, dps=[item[1] for item in observations],
        )
        if gene_result is not None:
            gene_results[chrom] = gene_result

    gene_values = np.asarray([result["chi_r"] for result in gene_results.values()])
    gene_median = float(np.median(gene_values)) if len(gene_values) else None
    gene_mad = (
        float(np.median(np.abs(gene_values - gene_median)))
        if gene_median is not None else None
    )
    gene_delta = (
        abs(res["chi_r"] - gene_median) if gene_median is not None else None
    )
    genes_agreeing = int(np.sum(np.abs(gene_values - res["chi_r"]) <= 0.05))

    bootstrap = bootstrap_chi(
        [r[2] for r in rows], [r[3] for r in rows], args.min_chi,
        args.prior_chi, args.bootstrap, args.bootstrap_seed,
    )
    reasons = []
    mismatch = []
    if gene_delta is not None and gene_delta >= MODEL_MISMATCH_GENE_DELTA:
        mismatch.append("cross_gene_delta")
    if res["peak_ratio"] is not None and res["peak_ratio"] <= MODEL_MISMATCH_PEAK_RATIO:
        mismatch.append("weak_selected_peak")
    if (res["residual_margin"] is not None
            and res["residual_margin"] >= MODEL_MISMATCH_RESIDUAL_MARGIN):
        mismatch.append("better_alternative_fit")
    if mismatch:
        status = "MODEL_MISMATCH"
        reasons.extend(mismatch)
    else:
        low_confidence = []
        if bootstrap is None:
            low_confidence.append("bootstrap_disabled")
        else:
            ci_width = bootstrap["ci_high"] - bootstrap["ci_low"]
            if bootstrap["finite_fraction"] < LOW_CONFIDENCE_BOOTSTRAP_FINITE:
                low_confidence.append("bootstrap_failures")
            if ci_width > LOW_CONFIDENCE_CI_WIDTH:
                low_confidence.append("wide_interval")
        status = "LOW_CONFIDENCE" if low_confidence else "PASS"
        reasons.extend(low_confidence)

    reported_chi = res["chi_r"] if status == "PASS" else None
    dps = np.asarray([r[3] for r in rows], dtype=float)
    ci_low = bootstrap["ci_low"] if bootstrap else None
    ci_high = bootstrap["ci_high"] if bootstrap else None
    print(
        f"GLOBAL  chi_R={fmt(reported_chi)}  raw_chi_R={fmt(res['chi_r'])}  "
        f"status={status}  reasons={','.join(reasons) if reasons else 'none'}  "
        f"ci95={fmt(ci_low)}-{fmt(ci_high)}  "
        f"bootstrap_peak_odds={fmt(bootstrap['peak_odds'] if bootstrap else None, 2)}  "
        f"peak_support={fmt(bootstrap['peak_support'] if bootstrap else None, 3)}  "
        f"peak_ratio={fmt(res['peak_ratio'], 3)}  "
        f"gene_median={fmt(gene_median)}  gene_mad={fmt(gene_mad)}  "
        f"gene_delta={fmt(gene_delta)}  genes_valid={len(gene_values)}  "
        f"genes_agreeing={genes_agreeing}  n_af={res['n_af']}  "
        f"median_dp={fmt(float(np.median(dps)), 1)}  "
        f"total_dp={int(np.sum(dps))}  residual={fmt(res['weighted_residual'], 5)}  "
        f"residual_margin={fmt(res['residual_margin'], 5)}  "
        f"selection={res['selection']}  peaks={res['all_peaks']}"
    )
    if args.per_gene:
        for chrom, observations in sorted(per.items()):
            gene_result = gene_results.get(chrom)
            if gene_result is None:
                print(f"  {chrom:14s} n={len(observations):4d}  (no reliable peak)")
            else:
                print(f"  {chrom:14s} n={len(observations):4d}  "
                      f"chi_R={gene_result['chi_r']:.4f}  "
                      f"residual={gene_result['weighted_residual']:.5f}  "
                      f"peaks={gene_result['all_peaks'][:3]}")


if __name__ == "__main__":
    main()
