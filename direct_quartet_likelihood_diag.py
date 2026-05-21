#!/usr/bin/env python
"""Diagnostic direct read-level quartet likelihood on existing EM outputs.

This script reuses <sample>/em_refine/<gene>.tf_counts.tsv and
<gene>.aug.fa, remaps the gene FASTQs to the augmented reference, and scores
quartets directly from read-level alignment likelihoods.
"""
import argparse
import csv
import os
import re
import subprocess
import tempfile
from collections import Counter

from iterative_remap_em import (
    bwa_mem_all,
    fit_4hap_read_likelihood,
    parse_sam_to_reads,
    safe,
    two_field,
)


GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]


def read_tf_counts(path):
    counts = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            counts[row["allele_2field"]] = float(row["em_weight"])
    return counts


def read_safe2name(fasta_path):
    out = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            name = line[1:].split()[0]
            if "*" not in name and "_" in name:
                parts = name.split("_")
                name = parts[0] + "*" + ":".join(parts[1:])
            out[safe(name)] = name
    return out


def read_calls_quartet(path, gene=None):
    if not path or not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            if gene and row.get("gene") and row.get("gene") != gene:
                continue
            if all(row.get(k) for k in ("R1_2field", "R2_2field", "D1_2field", "D2_2field")):
                return tuple(row[k] for k in ("R1_2field", "R2_2field", "D1_2field", "D2_2field"))
            allele = row.get("allele") or row.get("allele_2field") or row.get("call")
            if not allele:
                continue
            rows.append((row.get("assignment", ""), two_field(allele)))
    if len(rows) < 4:
        return None
    recip = [a for side, a in rows if side == "R"]
    donor = [a for side, a in rows if side == "D"]
    if len(recip) == 2 and len(donor) == 2:
        return tuple(recip + donor)
    return tuple(a for _, a in rows[:4])


def parse_chi(path):
    with open(path) as fh:
        text = fh.read()
    m = re.search(r"chi_R=([0-9.]+)", text)
    if not m:
        raise ValueError(f"could not parse chi_R from {path}")
    return float(m.group(1))


def parse_summary_chi(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return None
    try:
        return float(rows[0].get("chi_r_fit", ""))
    except ValueError:
        return None


def score_pair(pred, truth):
    cp = Counter(pred)
    ct = Counter(truth)
    return sum(min(cp[k], ct[k]) for k in cp)


def quartet_score(pred, truth_r, truth_d):
    return score_pair(pred[:2], truth_r) + score_pair(pred[2:], truth_d)


def split_alleles(text):
    return [x for x in text.split(",") if x]


def load_summary(path, samples, genes, max_rows):
    rows = []
    samples = set(samples or [])
    genes = set(genes or [])
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            if samples and row["sample"] not in samples:
                continue
            if genes and row["gene"] not in genes:
                continue
            rows.append(row)
            if max_rows and len(rows) >= max_rows:
                break
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True,
                    help="quartet_summary TSV with sample,set,gene,truth_R,truth_D")
    ap.add_argument("--spechla-root", required=True,
                    help="root containing per-sample gene FASTQs and em_refine outputs")
    ap.add_argument("--asm-root", required=True,
                    help="root containing per-sample final_calls.tsv and baseline gene calls")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--sample", action="append", default=[])
    ap.add_argument("--gene", action="append", default=[])
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-as-frac", type=float, default=0.95)
    ap.add_argument("--em-T", type=float, default=2.0)
    ap.add_argument("--direct-top-n", type=int, default=10)
    ap.add_argument("--direct-min-frac", type=float, default=0.002)
    ap.add_argument("--direct-log-floor", type=float, default=-25.0)
    ap.add_argument("--direct-per-gene-chi", action="store_true")
    ap.add_argument("--chi-lo", type=float, default=0.005)
    ap.add_argument("--chi-hi", type=float, default=0.5)
    ap.add_argument("--direct-chi-step", type=float, default=0.01)
    ap.add_argument("--direct-chi-prior", type=float, default=0.0)
    ap.add_argument("--direct-dose-prior", type=float, default=0.0)
    ap.add_argument("--direct-family-agg", choices=["max", "logsum"], default="max")
    ap.add_argument("--sam-cache-dir", default=None)
    ap.add_argument("--keep-sam", action="store_true")
    args = ap.parse_args()

    todo = load_summary(args.summary, args.sample, args.gene or GENES, args.max_rows)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    if args.sam_cache_dir:
        os.makedirs(args.sam_cache_dir, exist_ok=True)

    fields = [
        "sample", "set", "gene", "current_score", "direct_score",
        "delta", "direct_score_swap", "n_reads", "nll", "gap", "chi_direct",
        "direct_R", "direct_D", "truth_R", "truth_D", "current_quartet",
        "baseline_quartet",
    ]
    with open(args.out_tsv, "w") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for i, row in enumerate(todo, 1):
            sample = row["sample"]
            gene = row["gene"]
            short = gene.split("-")[1]
            sample_spechla = os.path.join(args.spechla_root, sample)
            sample_asm = os.path.join(args.asm_root, sample)
            em_dir = os.path.join(sample_spechla, "em_refine")
            ref_fa = os.path.join(em_dir, f"{gene}.aug.fa")
            tf_path = os.path.join(em_dir, f"{gene}.tf_counts.tsv")
            summary_path = os.path.join(em_dir, f"{gene}.summary.tsv")
            fq1 = os.path.join(sample_spechla, f"{short}.R1.fq.gz")
            fq2 = os.path.join(sample_spechla, f"{short}.R2.fq.gz")
            chi_path = os.path.join(sample_spechla, f"{sample}.chi_pooled.txt")
            if not all(os.path.exists(p) for p in (ref_fa, tf_path, fq1, fq2, chi_path)):
                print(f"skip missing input: {sample} {gene}", flush=True)
                continue
            counts = read_tf_counts(tf_path)
            safe2name = read_safe2name(ref_fa)
            contig_set = set(safe2name)
            chi_r = parse_summary_chi(summary_path)
            if chi_r is None:
                chi_r = parse_chi(chi_path)
            final_path = os.path.join(sample_asm, f"{sample}.final_calls.tsv")
            baseline_path = os.path.join(sample_asm, gene.lower(), gene, "calls.tsv")
            current = read_calls_quartet(final_path, gene=gene)
            baseline = read_calls_quartet(baseline_path)
            force_names = set(current or []) | set(baseline or [])
            force_quartets = [q for q in (current, baseline) if q]
            if args.sam_cache_dir:
                sam = os.path.join(args.sam_cache_dir, f"{sample}.{gene}.aug.sam")
            else:
                tmp = tempfile.NamedTemporaryFile(prefix=f"{sample}.{gene}.", suffix=".sam", delete=False)
                sam = tmp.name
                tmp.close()
            if not os.path.exists(sam) or os.path.getsize(sam) == 0:
                bwa_mem_all(ref_fa, fq1, fq2, sample, args.threads, sam)
            reads = parse_sam_to_reads(sam, contig_set, args.min_as_frac)
            if not args.keep_sam and not args.sam_cache_dir:
                os.unlink(sam)
            direct, nll, chi_direct, gap = fit_4hap_read_likelihood(
                reads, safe2name, counts, chi_r,
                top_n=args.direct_top_n,
                min_frac=args.direct_min_frac,
                per_gene_chi=args.direct_per_gene_chi,
                chi_lo=args.chi_lo,
                chi_hi=args.chi_hi,
                chi_step=args.direct_chi_step,
                chi_prior_lambda=args.direct_chi_prior,
                dose_prior_lambda=args.direct_dose_prior,
                T=args.em_T,
                log_floor=args.direct_log_floor,
                force_names=force_names,
                force_quartets=force_quartets,
                family_agg=args.direct_family_agg,
            )
            if not direct:
                print(f"direct failed: {sample} {gene}", flush=True)
                continue
            truth_r = split_alleles(row["truth_R"])
            truth_d = split_alleles(row["truth_D"])
            direct_score = quartet_score(direct, truth_r, truth_d)
            direct_score_swap = quartet_score(direct, truth_d, truth_r)
            current_score = int(row["score2"])
            w.writerow({
                "sample": sample,
                "set": row["set"],
                "gene": gene,
                "current_score": current_score,
                "direct_score": direct_score,
                "delta": direct_score - current_score,
                "direct_score_swap": direct_score_swap,
                "n_reads": len(reads),
                "nll": f"{nll:.3f}",
                "gap": f"{gap:.3f}",
                "chi_direct": f"{chi_direct:.4f}",
                "direct_R": ",".join(direct[:2]),
                "direct_D": ",".join(direct[2:]),
                "truth_R": row["truth_R"],
                "truth_D": row["truth_D"],
                "current_quartet": ",".join(current or []),
                "baseline_quartet": ",".join(baseline or []),
            })
            out.flush()
            print(f"[{i}/{len(todo)}] {sample} {gene}: {current_score}->{direct_score} "
                  f"gap={gap:.2f} nreads={len(reads)}", flush=True)


if __name__ == "__main__":
    main()