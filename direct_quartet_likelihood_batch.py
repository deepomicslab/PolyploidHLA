#!/usr/bin/env python3
"""Batch direct read-level quartet likelihood for selected summary rows."""
from __future__ import annotations

import argparse
import csv
import sys
import traceback
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from direct_quartet_likelihood_diag import (  # noqa: E402
    parse_chi,
    parse_summary_chi,
    quartet_score,
    read_calls_quartet,
    read_safe2name,
    read_tf_counts,
    split_alleles,
)
from iterative_remap_em import (  # noqa: E402
    bwa_mem_all,
    fit_4hap_read_likelihood,
    parse_sam_to_reads,
)


def load_rows(path: Path, genes, samples):
    gene_set = set(genes or [])
    sample_set = set(samples or [])
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if gene_set and row["gene"] not in gene_set:
                continue
            if sample_set and row["sample"] not in sample_set:
                continue
            yield row


def make_runtime_rows(genes, samples):
    if not genes:
        raise SystemExit("--gene is required when --summary is not provided")
    if not samples:
        raise SystemExit("--sample is required when --summary is not provided")
    for sample in samples:
        for gene in genes:
            yield {
                "sample": sample,
                "set": "",
                "gene": gene,
                "score2": "",
                "truth_R": "",
                "truth_D": "",
            }


def load_done_keys(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return set()
    with path.open() as handle:
        return {(row["sample"], row["gene"]) for row in csv.DictReader(handle, delimiter="\t")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=None,
                        help="Optional diagnostic summary with truth/score columns. "
                             "If omitted, --sample and --gene define truth-free runtime rows.")
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--sam-cache-dir", required=True, type=Path)
    parser.add_argument("--gene", action="append", default=[])
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-as-frac", type=float, default=0.95)
    parser.add_argument("--em-T", type=float, default=2.0)
    parser.add_argument("--direct-top-n", type=int, default=8)
    parser.add_argument("--direct-min-frac", type=float, default=0.002)
    parser.add_argument("--direct-log-floor", type=float, default=-25.0)
    parser.add_argument("--direct-family-agg", choices=("max", "logsum"), default="max")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fail-tsv", type=Path, default=None)
    args = parser.parse_args()

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.sam_cache_dir.mkdir(parents=True, exist_ok=True)
    if args.summary:
        rows = list(load_rows(args.summary, args.gene, args.sample))
    else:
        rows = list(make_runtime_rows(args.gene, args.sample))
    fields = [
        "sample", "set", "gene", "current_score", "direct_score", "delta",
        "direct_score_swap", "n_reads", "nll", "gap", "chi_direct",
        "direct_R", "direct_D", "truth_R", "truth_D", "current_quartet",
        "baseline_quartet",
    ]
    fail_fields = ["sample", "set", "gene", "error"]
    done = load_done_keys(args.out_tsv) if args.resume else set()
    mode = "a" if args.resume and args.out_tsv.exists() and args.out_tsv.stat().st_size > 0 else "w"
    with args.out_tsv.open(mode) as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        if mode == "w":
            writer.writeheader()
        for index, row in enumerate(rows, 1):
            sample = row["sample"]
            gene = row["gene"]
            if (sample, gene) in done:
                print(f"[{index}/{len(rows)}] skip done {sample} {gene}", flush=True)
                continue
            short = gene.split("-")[1]
            sample_spechla = args.spechla_root / sample
            sample_asm = args.asm_root / sample
            em_dir = sample_spechla / "em_refine"
            ref_fa = em_dir / f"{gene}.aug.fa"
            tf_path = em_dir / f"{gene}.tf_counts.tsv"
            summary_path = em_dir / f"{gene}.summary.tsv"
            fq1 = sample_spechla / f"{short}.R1.fq.gz"
            fq2 = sample_spechla / f"{short}.R2.fq.gz"
            chi_path = sample_spechla / f"{sample}.chi_pooled.txt"
            if not all(input_path.exists() for input_path in (ref_fa, tf_path, fq1, fq2, chi_path)):
                print(f"[{index}/{len(rows)}] skip missing input: {sample} {gene}", flush=True)
                continue
            try:
                counts = read_tf_counts(str(tf_path))
                safe2name = read_safe2name(str(ref_fa))
                contig_set = set(safe2name)
                chi_r = parse_summary_chi(str(summary_path))
                if chi_r is None:
                    chi_r = parse_chi(str(chi_path))
                final_path = sample_asm / f"{sample}.final_calls.tsv"
                baseline_path = sample_asm / gene.lower() / gene / "calls.tsv"
                current = read_calls_quartet(str(final_path), gene=gene)
                baseline = read_calls_quartet(str(baseline_path))
                force_names = set(current or []) | set(baseline or [])
                force_quartets = [quartet for quartet in (current, baseline) if quartet]
                sam = args.sam_cache_dir / f"{sample}.{gene}.aug.sam"
                if not sam.exists() or sam.stat().st_size == 0:
                    bwa_mem_all(str(ref_fa), str(fq1), str(fq2), sample, args.threads, str(sam))
                reads = parse_sam_to_reads(str(sam), contig_set, args.min_as_frac)
                direct, nll, chi_direct, gap = fit_4hap_read_likelihood(
                    reads,
                    safe2name,
                    counts,
                    chi_r,
                    top_n=args.direct_top_n,
                    min_frac=args.direct_min_frac,
                    T=args.em_T,
                    log_floor=args.direct_log_floor,
                    force_names=force_names,
                    force_quartets=force_quartets,
                    family_agg=args.direct_family_agg,
                )
                if not direct:
                    print(f"[{index}/{len(rows)}] direct failed: {sample} {gene}", flush=True)
                    continue
                truth_r = split_alleles(row.get("truth_R", ""))
                truth_d = split_alleles(row.get("truth_D", ""))
                score2_text = row.get("score2", "")
                has_truth = bool(truth_r and truth_d and score2_text and score2_text != "NA")
                if has_truth:
                    direct_score = quartet_score(direct, truth_r, truth_d)
                    direct_score_swap = quartet_score(direct, truth_d, truth_r)
                    current_score = int(row["score2"])
                    delta = direct_score - current_score
                else:
                    direct_score = "NA"
                    direct_score_swap = "NA"
                    current_score = "NA"
                    delta = "NA"
                writer.writerow({
                    "sample": sample,
                    "set": row.get("set", ""),
                    "gene": gene,
                    "current_score": current_score,
                    "direct_score": direct_score,
                    "delta": delta,
                    "direct_score_swap": direct_score_swap,
                    "n_reads": len(reads),
                    "nll": f"{nll:.3f}",
                    "gap": f"{gap:.3f}",
                    "chi_direct": f"{chi_direct:.4f}",
                    "direct_R": ",".join(direct[:2]),
                    "direct_D": ",".join(direct[2:]),
                    "truth_R": row.get("truth_R", ""),
                    "truth_D": row.get("truth_D", ""),
                    "current_quartet": ",".join(current or []),
                    "baseline_quartet": ",".join(baseline or []),
                })
                out.flush()
                print(
                    f"[{index}/{len(rows)}] {sample} {gene}: {current_score}->{direct_score} "
                    f"gap={gap:.2f} nreads={len(reads)}",
                    flush=True,
                )
            except Exception as exc:  # diagnostic batch should keep going across loci
                message = f"{type(exc).__name__}: {exc}"
                print(f"[{index}/{len(rows)}] failed {sample} {gene}: {message}", flush=True)
                if args.fail_tsv:
                    fail_mode = "a" if args.fail_tsv.exists() and args.fail_tsv.stat().st_size > 0 else "w"
                    with args.fail_tsv.open(fail_mode) as fail_handle:
                        fail_writer = csv.DictWriter(fail_handle, delimiter="\t", fieldnames=fail_fields)
                        if fail_mode == "w":
                            fail_writer.writeheader()
                        fail_writer.writerow({
                            "sample": sample,
                            "set": row["set"],
                            "gene": gene,
                            "error": message,
                        })
                traceback.print_exc()
    print(f"wrote {args.out_tsv}")


if __name__ == "__main__":
    main()
