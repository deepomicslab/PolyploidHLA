#!/usr/bin/env python3
"""Check whether truth-allele variants are called in existing pipeline VCFs.

For each sample and truth allele, align a representative IMGT allele sequence to
SpecHLA's HLA_REF gene segment, derive truth SNVs relative to the pipeline
reference, then check whether freebayes/regt/phased VCFs contain those ALT bases.

Truth is used only for diagnostics/scoring, not for production calling.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from hla_polyphase_assemble import _align_allele_to_ref, load_imgt_alleles  # noqa: E402

DEFAULT_HLA_REF = SCRIPT_DIR / "resources" / "spechla" / "db" / "ref" / "hla.ref.extend.fa"
DEFAULT_IMGT = SCRIPT_DIR / "resources" / "spechla" / "db" / "ref" / "hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
DEFAULT_GENE_BED = SCRIPT_DIR / "gene.spechla.bed"
DEFAULT_TRUTH_DIR = SCRIPT_DIR.parent / "truth"
DEFAULT_SPECHLA_ROOT = Path("/data2/wangxuedong/polyploid-hla-realsets/spechla_out_gendx_amp_abc_20260520")

GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]
SIDES = ("PATIENT", "DONOR")


def two_field(allele: str) -> str:
    allele = allele.replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    rest = rest.rstrip("GP")
    fields = rest.split(":")
    return f"{gene}*{':'.join(fields[:2])}" if len(fields) >= 2 else f"{gene}*{fields[0]}"


def strip_suffix(allele: str) -> str:
    allele = allele.replace("HLA-", "").strip()
    if allele and allele[-1].isalpha() and allele[-1] != "G":
        allele = allele[:-1]
    if allele.endswith("G") or allele.endswith("P"):
        allele = allele[:-1]
    return allele


def normalize_truth_allele(gene: str, allele: str) -> str:
    short = gene.replace("HLA-", "")
    allele = strip_suffix(allele)
    if "*" in allele:
        return allele
    return f"{short}*{allele}"


def set_label(sample: str) -> str:
    first = sample[0].upper()
    if first == "A":
        return "set-a"
    if first == "B":
        return "set-b"
    if first == "C":
        return "set-c"
    raise ValueError(f"cannot infer set for {sample}")


def tag_for_gene(gene: str) -> str:
    return gene.lower()


def bam_name_for_gene(gene: str) -> str:
    return gene.replace("HLA-", "") + ".bam"


def load_gene_bed(path: Path) -> dict[str, tuple[str, int, int]]:
    out = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gene, *_ = line.split()
            out[gene] = (chrom, int(start), int(end))
    return out


def load_truth(path: Path) -> dict[str, dict[str, list[str]]]:
    with path.open() as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    header = rows[0][1:]
    truth: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in rows[1:]:
        side = row[0]
        for gene_short, allele in zip(header, row[1:]):
            gene = f"HLA-{gene_short}"
            if gene in GENES:
                truth[side][gene].append(normalize_truth_allele(gene, allele))
    return truth


def choose_imgt_allele(allele: str, imgt: dict[str, str]) -> tuple[str | None, str]:
    allele = strip_suffix(allele)
    candidates = []
    if allele in imgt:
        return allele, "exact"
    for name in imgt:
        if name == allele or name.startswith(allele + ":"):
            candidates.append(name)
    if candidates:
        return sorted(candidates, key=lambda x: (len(x), x))[0], "prefix"
    tf = two_field(allele)
    for name in imgt:
        if two_field(name) == tf:
            candidates.append(name)
    if candidates:
        return sorted(candidates, key=lambda x: (len(x), x))[0], "two_field"
    return None, "missing_imgt"


def iter_truth_snvs(ref_aln: str, allele_aln: str, chrom: str, gene_start: int) -> Iterable[dict[str, object]]:
    ref_index = 0
    insertion_count = 0
    deletion_count = 0
    for ref_base, allele_base in zip(ref_aln, allele_aln):
        if ref_base != "-" and allele_base != "-":
            ref_index += 1
            if ref_base.upper() in "ACGT" and allele_base.upper() in "ACGT" and ref_base.upper() != allele_base.upper():
                yield {
                    "kind": "SNV",
                    "chrom": chrom,
                    "pos": gene_start + ref_index,
                    "ref": ref_base.upper(),
                    "alt": allele_base.upper(),
                }
        elif ref_base != "-" and allele_base == "-":
            ref_index += 1
            deletion_count += 1
        elif ref_base == "-" and allele_base != "-":
            insertion_count += 1
    yield {"kind": "INDEL_SUMMARY", "insertions": insertion_count, "deletions": deletion_count}


def vcf_key_stats(vcf_path: Path) -> dict[tuple[str, int, str, str], dict[str, object]]:
    out = {}
    if not vcf_path.exists():
        return out
    try:
        vcf = pysam.VariantFile(str(vcf_path))
    except Exception as exc:
        print(f"WARN: cannot read VCF {vcf_path}: {exc}", file=sys.stderr)
        return out
    for rec in vcf:
        if len(rec.ref) != 1 or not rec.alts:
            continue
        sample = next(iter(rec.samples.values())) if rec.samples else None
        gt = sample.get("GT") if sample else None
        ad = sample.get("AD") if sample else None
        for idx, alt in enumerate(rec.alts, 1):
            if len(alt) != 1:
                continue
            key = (rec.chrom, rec.pos, rec.ref.upper(), alt.upper())
            alt_count = "NA"
            total_depth = "NA"
            alt_fraction = "NA"
            if ad and len(ad) > idx and ad[idx] is not None:
                alt_count = ad[idx]
                total = sum(x for x in ad if x is not None)
                total_depth = total
                alt_fraction = f"{(ad[idx] / total):.5f}" if total else "NA"
            out[key] = {
                "gt": "/".join("." if x is None else str(x) for x in gt) if gt else "NA",
                "gt_has_alt": bool(gt and idx in gt),
                "alt_count": alt_count,
                "total_depth": total_depth,
                "alt_fraction": alt_fraction,
            }
    return out


def summarize_hits(snvs: list[dict[str, object]], calls: dict[tuple[str, int, str, str], dict[str, object]]) -> tuple[int, int, int, str]:
    expected = len(snvs)
    present = 0
    in_gt = 0
    examples = []
    for snv in snvs:
        key = (str(snv["chrom"]), int(snv["pos"]), str(snv["ref"]), str(snv["alt"]))
        hit = calls.get(key)
        if hit:
            present += 1
            if hit.get("gt_has_alt"):
                in_gt += 1
            if len(examples) < 5:
                examples.append(f"{key[0]}:{key[1]}:{key[2]}>{key[3]}:AF={hit.get('alt_fraction')}:GT={hit.get('gt')}")
    return expected, present, in_gt, ";".join(examples)


def bam_depth_array(bam_path: Path, chrom: str, start: int, end: int) -> tuple[list[int] | None, str]:
    if not bam_path.exists():
        return None, "missing_bam"
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        cov = bam.count_coverage(chrom, start, end, quality_threshold=0)
    return [sum(base_cov[i] for base_cov in cov) for i in range(end - start)], "ok"


def covered_snvs(
    snvs: list[dict[str, object]],
    bam_path: Path,
    min_depth: int,
    gene_start: int | None = None,
    depth_arr: list[int] | None = None,
) -> tuple[list[dict[str, object]], str]:
    if min_depth <= 0:
        return snvs, "NA"
    if depth_arr is None or gene_start is None:
        if not bam_path.exists():
            return [], "missing_bam"
        if not snvs:
            return [], "0"
        chrom = str(snvs[0]["chrom"])
        start0 = min(int(snv["pos"]) - 1 for snv in snvs)
        end0 = max(int(snv["pos"]) for snv in snvs)
        depth_arr, status = bam_depth_array(bam_path, chrom, start0, end0)
        if depth_arr is None:
            return [], status
        gene_start = start0
    if depth_arr is None:
        return [], "missing_bam"
    kept = []
    depths = []
    for snv in snvs:
        idx = int(snv["pos"]) - 1 - gene_start
        depth = depth_arr[idx] if 0 <= idx < len(depth_arr) else 0
        if depth >= min_depth:
            item = dict(snv)
            item["depth"] = depth
            kept.append(item)
            depths.append(depth)
    if not depths:
        return kept, "0"
    return kept, f"min={min(depths)};median={sorted(depths)[len(depths)//2]};max={max(depths)}"


def write_summary(raw_path: Path, summary_path: Path) -> None:
    rows = list(csv.DictReader(raw_path.open(), delimiter="\t"))
    agg = defaultdict(lambda: defaultdict(int))
    for row in rows:
        key = (row["sample"], row["gene"])
        for field in (
            "truth_snvs", "covered_snvs", "freebayes_present", "freebayes_gt",
            "regt_present", "regt_gt", "phased_present", "phased_gt",
            "freebayes_present_covered", "freebayes_gt_covered",
            "regt_gt_covered", "phased_gt_covered",
        ):
            value = row.get(field) or "0"
            if value != "NA":
                agg[key][field] += int(value)
    fields = [
        "sample", "gene", "truth_snvs", "covered_snvs",
        "freebayes_present_rate", "freebayes_gt_rate", "regt_gt_rate", "phased_gt_rate",
        "freebayes_present_covered_rate", "freebayes_gt_covered_rate",
        "regt_gt_covered_rate", "phased_gt_covered_rate",
    ]
    with summary_path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for (sample, gene), vals in sorted(agg.items()):
            truth_snvs = vals["truth_snvs"]
            covered = vals["covered_snvs"]

            def rate(num: int, denom: int) -> str:
                return "NA" if denom == 0 else f"{num}/{denom} ({num / denom:.3f})"

            writer.writerow({
                "sample": sample,
                "gene": gene,
                "truth_snvs": truth_snvs,
                "covered_snvs": covered,
                "freebayes_present_rate": rate(vals["freebayes_present"], truth_snvs),
                "freebayes_gt_rate": rate(vals["freebayes_gt"], truth_snvs),
                "regt_gt_rate": rate(vals["regt_gt"], truth_snvs),
                "phased_gt_rate": rate(vals["phased_gt"], truth_snvs),
                "freebayes_present_covered_rate": rate(vals["freebayes_present_covered"], covered),
                "freebayes_gt_covered_rate": rate(vals["freebayes_gt_covered"], covered),
                "regt_gt_covered_rate": rate(vals["regt_gt_covered"], covered),
                "phased_gt_covered_rate": rate(vals["phased_gt_covered"], covered),
            })


def filter_raw(raw_in: Path, raw_out: Path, exclude_samples: set[str]) -> None:
    with raw_in.open() as in_handle, raw_out.open("w") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row.get("sample") not in exclude_samples:
                writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="+", default=[])
    parser.add_argument("--discover-samples", action="store_true",
                        help="use sample directories under --spechla-root")
    parser.add_argument("--exclude-samples", nargs="+", default=[],
                        help="sample IDs to skip even if listed in --samples")
    parser.add_argument("--spechla-root", type=Path, default=DEFAULT_SPECHLA_ROOT)
    parser.add_argument("--truth-dir", type=Path, default=DEFAULT_TRUTH_DIR)
    parser.add_argument("--hla-ref", type=Path, default=DEFAULT_HLA_REF)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--gene-bed", type=Path, default=DEFAULT_GENE_BED)
    parser.add_argument("--genes", nargs="+", default=GENES)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, default=None)
    parser.add_argument("--raw-in", type=Path, default=None,
                        help="existing raw TSV to filter/summarize instead of recomputing")
    parser.add_argument("--min-depth", type=int, default=10,
                        help="minimum per-position BAM depth for covered-SNV summaries")
    parser.add_argument("--aligner", choices=("parasail", "mappy"), default="parasail")
    args = parser.parse_args()

    exclude_samples = set(args.exclude_samples)
    if args.raw_in:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        filter_raw(args.raw_in, args.out, exclude_samples)
        print(f"wrote {args.out}")
        if args.summary_out:
            write_summary(args.out, args.summary_out)
            print(f"wrote {args.summary_out}")
        return

    samples = list(args.samples)
    if args.discover_samples:
        samples.extend(path.name for path in args.spechla_root.iterdir() if path.is_dir())
    samples = sorted(set(samples))
    if not samples:
        parser.error("provide --samples or --discover-samples")
    samples = [sample for sample in samples if sample not in exclude_samples]

    genes = [gene for gene in args.genes if gene in GENES]
    bed = load_gene_bed(args.gene_bed)
    imgt = load_imgt_alleles(str(args.imgt))
    ref = pysam.FastaFile(str(args.hla_ref))
    truth_cache = {label: load_truth(args.truth_dir / f"truth_typing-{label}.tsv") for label in ("set-a", "set-b", "set-c")}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample", "set", "side", "gene", "truth_allele", "truth_2field", "imgt_allele", "imgt_match",
        "truth_snvs", "covered_snvs", "covered_depth_summary", "freebayes_present", "freebayes_gt",
        "regt_present", "regt_gt", "phased_present", "phased_gt", "freebayes_present_covered",
        "freebayes_gt_covered", "regt_gt_covered", "phased_gt_covered", "freebayes_examples",
    ]
    with args.out.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for sample in samples:
            label = set_label(sample)
            truth = truth_cache[label]
            vcf_cache = {}
            for gene in genes:
                chrom, start, end = bed[gene]
                gene_ref = ref.fetch(chrom, start, end).upper()
                tag = tag_for_gene(gene)
                vcf_cache[(gene, "freebayes")] = vcf_key_stats(args.spechla_root / sample / f"{sample}.freebayes.{tag}.vcf.gz")
                vcf_cache[(gene, "regt")] = vcf_key_stats(args.spechla_root / sample / f"{sample}.freebayes_regt.{tag}.vcf.gz")
                vcf_cache[(gene, "phased")] = vcf_key_stats(args.spechla_root / sample / f"{sample}.phased.{tag}.vcf.gz")
                bam_path = args.spechla_root / sample / bam_name_for_gene(gene)
                depth_arr, depth_status = bam_depth_array(bam_path, chrom, start, end) if args.min_depth > 0 else (None, "NA")
                for side in SIDES:
                    for allele in sorted(set(truth[side][gene])):
                        imgt_allele, match = choose_imgt_allele(allele, imgt)
                        if imgt_allele is None:
                            writer.writerow({
                                "sample": sample, "set": label, "side": side, "gene": gene,
                                "truth_allele": allele, "truth_2field": two_field(allele),
                                "imgt_allele": "NA", "imgt_match": match,
                            })
                            continue
                        amap = _align_allele_to_ref(imgt[imgt_allele], gene_ref, args.aligner)
                        snvs = []
                        for item in iter_truth_snvs(amap[0], amap[1], chrom, start):
                            if item["kind"] == "SNV":
                                snvs.append(item)
                        cov_snvs, depth_summary = covered_snvs(snvs, bam_path, args.min_depth, start, depth_arr)
                        if depth_status != "ok" and args.min_depth > 0:
                            depth_summary = depth_status
                        fb_expected, fb_present, fb_gt, fb_examples = summarize_hits(snvs, vcf_cache[(gene, "freebayes")])
                        _r_expected, regt_present, regt_gt, _regt_examples = summarize_hits(snvs, vcf_cache[(gene, "regt")])
                        _p_expected, phased_present, phased_gt, _phased_examples = summarize_hits(snvs, vcf_cache[(gene, "phased")])
                        _cfb_expected, fb_present_cov, fb_gt_cov, _fb_cov_examples = summarize_hits(cov_snvs, vcf_cache[(gene, "freebayes")])
                        _cr_expected, _regt_present_cov, regt_gt_cov, _regt_cov_examples = summarize_hits(cov_snvs, vcf_cache[(gene, "regt")])
                        _cp_expected, _phased_present_cov, phased_gt_cov, _phased_cov_examples = summarize_hits(cov_snvs, vcf_cache[(gene, "phased")])
                        writer.writerow({
                            "sample": sample,
                            "set": label,
                            "side": side,
                            "gene": gene,
                            "truth_allele": allele,
                            "truth_2field": two_field(allele),
                            "imgt_allele": imgt_allele,
                            "imgt_match": match,
                            "truth_snvs": fb_expected,
                            "covered_snvs": len(cov_snvs),
                            "covered_depth_summary": depth_summary,
                            "freebayes_present": fb_present,
                            "freebayes_gt": fb_gt,
                            "regt_present": regt_present,
                            "regt_gt": regt_gt,
                            "phased_present": phased_present,
                            "phased_gt": phased_gt,
                            "freebayes_present_covered": fb_present_cov,
                            "freebayes_gt_covered": fb_gt_cov,
                            "regt_gt_covered": regt_gt_cov,
                            "phased_gt_covered": phased_gt_cov,
                            "freebayes_examples": fb_examples,
                        })
    print(f"wrote {args.out}")
    if args.summary_out:
        write_summary(args.out, args.summary_out)
        print(f"wrote {args.summary_out}")


if __name__ == "__main__":
    main()
