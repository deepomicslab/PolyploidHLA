#!/usr/bin/env python3
"""Build per-truth-allele VCFs and compare them with pipeline VCFs.

Each truth allele sequence is aligned to the pipeline HLA reference with
minimap2. Variants are called from that alignment with bcftools and then
compared to the existing freebayes/freebayes_regt/phased VCFs.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from diagnose_truth_variants_called import (  # noqa: E402
    DEFAULT_GENE_BED,
    DEFAULT_HLA_REF,
    DEFAULT_IMGT,
    DEFAULT_SPECHLA_ROOT,
    DEFAULT_TRUTH_DIR,
    GENES,
    SIDES,
    bam_depth_array,
    bam_name_for_gene,
    choose_imgt_allele,
    load_gene_bed,
    load_truth,
    set_label,
    tag_for_gene,
    two_field,
)
from hla_polyphase_assemble import load_imgt_alleles  # noqa: E402

STAGES = ("freebayes", "regt", "phased")


def run_command(command: str) -> None:
    subprocess.run(command, shell=True, check=True, executable="/bin/bash")


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"required tool not found in PATH: {name}")
    return path


def safe_name(name: str) -> str:
    return name.replace("*", "_star_").replace(":", "_").replace("/", "_")


def write_fasta(path: Path, name: str, seq: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for offset in range(0, len(seq), 80):
            handle.write(seq[offset:offset + 80] + "\n")


def write_fastq(path: Path, name: str, seq: str) -> None:
    qualities = "I" * len(seq)
    with path.open("w") as handle:
        handle.write(f"@{name}\n{seq}\n+\n{qualities}\n")


def write_alignment_variants_vcf(bam_path: Path, ref_path: Path, out_vcf: Path) -> None:
    ref = pysam.FastaFile(str(ref_path))
    chrom = ref.references[0]
    ref_seq = ref.fetch(chrom).upper()
    variants: set[tuple[str, int, str, str, str]] = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            if record.is_unmapped or record.is_secondary:
                continue
            query_seq = (record.query_sequence or "").upper()
            ref_pos = record.reference_start
            query_pos = 0
            for operation, length in record.cigartuples or []:
                if operation in (0, 7, 8):
                    for offset in range(length):
                        ref_base = ref_seq[ref_pos + offset]
                        query_base = query_seq[query_pos + offset]
                        if ref_base in "ACGT" and query_base in "ACGT" and ref_base != query_base:
                            variants.add((chrom, ref_pos + offset + 1, ref_base, query_base, "SNV"))
                    ref_pos += length
                    query_pos += length
                elif operation == 1:
                    inserted = query_seq[query_pos:query_pos + length]
                    if ref_pos > 0 and inserted:
                        anchor = ref_seq[ref_pos - 1]
                        variants.add((chrom, ref_pos, anchor, anchor + inserted, "INDEL"))
                    query_pos += length
                elif operation == 2:
                    deleted = ref_seq[ref_pos:ref_pos + length]
                    if ref_pos > 0 and deleted:
                        anchor = ref_seq[ref_pos - 1]
                        variants.add((chrom, ref_pos, anchor + deleted, anchor, "INDEL"))
                    ref_pos += length
                elif operation in (4, 5):
                    if operation == 4:
                        query_pos += length
                elif operation == 3:
                    ref_pos += length

    with out_vcf.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID={chrom},length={len(ref_seq)}>\n")
        handle.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Alignment variant type\">\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref_base, alt_base, var_type in sorted(variants, key=lambda item: (item[0], item[1], item[2], item[3])):
            handle.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t60\tPASS\tTYPE={var_type}\n")


def ensure_gene_reference(ref: pysam.FastaFile, chrom: str, cache_dir: Path) -> Path:
    ref_dir = cache_dir / "gene_refs"
    ref_dir.mkdir(parents=True, exist_ok=True)
    path = ref_dir / f"{chrom}.fa"
    if not path.exists():
        write_fasta(path, chrom, ref.fetch(chrom).upper())
    if not Path(str(path) + ".fai").exists():
        run_command(f"samtools faidx {path}")
    return path


def build_truth_vcf(
    gene_ref: Path,
    gene: str,
    imgt_allele: str,
    allele_seq: str,
    cache_dir: Path,
    force: bool,
) -> Path:
    truth_dir = cache_dir / "truth_vcfs" / gene
    truth_dir.mkdir(parents=True, exist_ok=True)
    stem = safe_name(imgt_allele)
    fastq_path = truth_dir / f"{stem}.fq"
    bam_path = truth_dir / f"{stem}.bam"
    raw_vcf = truth_dir / f"{stem}.raw.vcf"
    out_vcf = truth_dir / f"{stem}.truth.vcf.gz"
    if out_vcf.exists() and Path(str(out_vcf) + ".tbi").exists() and not force:
        return out_vcf
    write_fastq(fastq_path, safe_name(imgt_allele), allele_seq.upper())
    run_command(
        f"minimap2 -a -x asm5 --eqx {gene_ref} {fastq_path} "
        f"| samtools sort -o {bam_path} -"
    )
    run_command(f"samtools index {bam_path}")
    write_alignment_variants_vcf(bam_path, gene_ref, raw_vcf)
    run_command(f"bcftools norm -f {gene_ref} -a -m -any -Oz -o {out_vcf} {raw_vcf}")
    run_command(f"bcftools index -t {out_vcf}")
    return out_vcf


def normalize_pipeline_vcf(vcf_path: Path, gene_ref: Path, cache_dir: Path, force: bool) -> Path | None:
    if not vcf_path.exists():
        return None
    out_dir = cache_dir / "pipeline_norm"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_vcf = out_dir / (vcf_path.name.replace(".vcf.gz", ".norm.vcf.gz"))
    if out_vcf.exists() and Path(str(out_vcf) + ".tbi").exists() and not force:
        return out_vcf
    try:
        run_command(f"bcftools norm -f {gene_ref} -a -m -any -Oz -o {out_vcf} {vcf_path}")
        run_command(f"bcftools index -t {out_vcf}")
    except subprocess.CalledProcessError:
        return None
    return out_vcf


def variant_type(ref_base: str, alt_base: str) -> str:
    if len(ref_base) == 1 and len(alt_base) == 1:
        return "SNV"
    return "INDEL"


def read_vcf_keys(path: Path | None) -> tuple[set[tuple[str, int, str, str]], set[tuple[str, int, str, str]], dict[tuple[str, int, str, str], str]]:
    present: set[tuple[str, int, str, str]] = set()
    gt_has_alt: set[tuple[str, int, str, str]] = set()
    types: dict[tuple[str, int, str, str], str] = {}
    if path is None or not path.exists():
        return present, gt_has_alt, types
    try:
        vcf = pysam.VariantFile(str(path))
    except Exception:
        return present, gt_has_alt, types
    for record in vcf:
        sample = next(iter(record.samples.values())) if record.samples else None
        genotype = sample.get("GT") if sample else None
        for alt_index, alt_base in enumerate(record.alts or (), 1):
            key = (record.chrom, record.pos, record.ref.upper(), alt_base.upper())
            present.add(key)
            types[key] = variant_type(record.ref, alt_base)
            if genotype and alt_index in genotype:
                gt_has_alt.add(key)
    return present, gt_has_alt, types


def coverage_filter(
    keys: set[tuple[str, int, str, str]],
    depth_arr: list[int] | None,
    min_depth: int,
) -> set[tuple[str, int, str, str]]:
    if min_depth <= 0 or depth_arr is None:
        return set(keys)
    kept = set()
    for chrom, pos, ref_base, alt_base in keys:
        depth = depth_arr[pos - 1] if 0 <= pos - 1 < len(depth_arr) else 0
        if depth >= min_depth:
            kept.add((chrom, pos, ref_base, alt_base))
    return kept


def count_hits(truth_keys: set[tuple[str, int, str, str]], call_keys: set[tuple[str, int, str, str]]) -> int:
    return len(truth_keys & call_keys)


def write_gene_summary(sample_gene_path: Path, gene_summary_path: Path) -> None:
    aggregate: dict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
    numeric_fields = [
        "truth_variants", "truth_snvs", "truth_indels", "covered_variants",
        "covered_snvs", "covered_indels", "freebayes_present_covered",
        "freebayes_gt_covered", "regt_gt_covered", "phased_gt_covered",
    ]
    with sample_gene_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values = aggregate[row["gene"]]
            for field in numeric_fields:
                values[field] += int(row[field])
    fields = [
        "gene", "truth_variants", "truth_snvs", "truth_indels", "covered_variants",
        "covered_snvs", "covered_indels", "freebayes_present_covered_rate",
        "freebayes_gt_covered_rate", "regt_gt_covered_rate", "phased_gt_covered_rate",
    ]
    with gene_summary_path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for gene, values in sorted(aggregate.items()):
            covered = values["covered_variants"]

            def rate(count: int) -> str:
                return "NA" if covered == 0 else f"{count}/{covered} ({count / covered:.3f})"

            writer.writerow({
                "gene": gene,
                "truth_variants": values["truth_variants"],
                "truth_snvs": values["truth_snvs"],
                "truth_indels": values["truth_indels"],
                "covered_variants": covered,
                "covered_snvs": values["covered_snvs"],
                "covered_indels": values["covered_indels"],
                "freebayes_present_covered_rate": rate(values["freebayes_present_covered"]),
                "freebayes_gt_covered_rate": rate(values["freebayes_gt_covered"]),
                "regt_gt_covered_rate": rate(values["regt_gt_covered"]),
                "phased_gt_covered_rate": rate(values["phased_gt_covered"]),
            })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="+", default=[])
    parser.add_argument("--discover-samples", action="store_true")
    parser.add_argument("--exclude-samples", nargs="+", default=[])
    parser.add_argument("--spechla-root", type=Path, default=DEFAULT_SPECHLA_ROOT)
    parser.add_argument("--truth-dir", type=Path, default=DEFAULT_TRUTH_DIR)
    parser.add_argument("--hla-ref", type=Path, default=DEFAULT_HLA_REF)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--gene-bed", type=Path, default=DEFAULT_GENE_BED)
    parser.add_argument("--genes", nargs="+", default=GENES)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=Path("diagnostics/truth_vcf_cache"))
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    for tool_name in ("minimap2", "samtools", "bcftools"):
        require_tool(tool_name)

    samples = list(args.samples)
    if args.discover_samples:
        samples.extend(path.name for path in args.spechla_root.iterdir() if path.is_dir())
    samples = sorted(set(samples) - set(args.exclude_samples))
    if not samples:
        parser.error("provide --samples or --discover-samples")

    genes = [gene for gene in args.genes if gene in GENES]
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    bed = load_gene_bed(args.gene_bed)
    imgt = load_imgt_alleles(str(args.imgt))
    ref = pysam.FastaFile(str(args.hla_ref))
    truth_cache = {label: load_truth(args.truth_dir / f"truth_typing-{label}.tsv") for label in ("set-a", "set-b", "set-c")}

    allele_path = args.out_prefix.with_suffix(".allele.tsv")
    sample_gene_path = args.out_prefix.with_suffix(".sample_gene_unique.tsv")
    gene_summary_path = args.out_prefix.with_suffix(".gene_unique_summary.tsv")

    allele_fields = [
        "sample", "set", "side", "gene", "truth_allele", "truth_2field", "imgt_allele",
        "imgt_match", "truth_variants", "truth_snvs", "truth_indels", "covered_variants",
        "covered_snvs", "covered_indels", "freebayes_present", "freebayes_gt", "regt_gt",
        "phased_gt", "freebayes_present_covered", "freebayes_gt_covered",
        "regt_gt_covered", "phased_gt_covered",
    ]
    sample_gene_fields = [
        "sample", "set", "gene", "truth_variants", "truth_snvs", "truth_indels",
        "covered_variants", "covered_snvs", "covered_indels", "freebayes_present",
        "freebayes_gt", "regt_gt", "phased_gt", "freebayes_present_covered",
        "freebayes_gt_covered", "regt_gt_covered", "phased_gt_covered",
    ]

    truth_vcf_cache: dict[tuple[str, str], tuple[set[tuple[str, int, str, str]], dict[tuple[str, int, str, str], str]]] = {}

    with allele_path.open("w") as allele_handle, sample_gene_path.open("w") as sample_gene_handle:
        allele_writer = csv.DictWriter(allele_handle, delimiter="\t", fieldnames=allele_fields)
        sample_gene_writer = csv.DictWriter(sample_gene_handle, delimiter="\t", fieldnames=sample_gene_fields)
        allele_writer.writeheader()
        sample_gene_writer.writeheader()

        for sample in samples:
            label = set_label(sample)
            truth = truth_cache[label]
            for gene in genes:
                chrom, _start, _end = bed[gene]
                gene_ref = ensure_gene_reference(ref, chrom, args.cache_dir)
                tag = tag_for_gene(gene)
                pipeline_paths = {
                    "freebayes": args.spechla_root / sample / f"{sample}.freebayes.{tag}.vcf.gz",
                    "regt": args.spechla_root / sample / f"{sample}.freebayes_regt.{tag}.vcf.gz",
                    "phased": args.spechla_root / sample / f"{sample}.phased.{tag}.vcf.gz",
                }
                pipeline_calls = {}
                for stage in STAGES:
                    norm_path = normalize_pipeline_vcf(pipeline_paths[stage], gene_ref, args.cache_dir, args.force)
                    pipeline_calls[stage] = read_vcf_keys(norm_path)
                depth_arr = None
                if args.min_depth > 0:
                    bam_path = args.spechla_root / sample / bam_name_for_gene(gene)
                    depth_arr, _depth_status = bam_depth_array(bam_path, chrom, 0, len(ref.fetch(chrom)))

                sample_gene_truth: set[tuple[str, int, str, str]] = set()
                sample_gene_types: dict[tuple[str, int, str, str], str] = {}
                for side in SIDES:
                    for truth_allele in sorted(set(truth[side][gene])):
                        imgt_allele, match = choose_imgt_allele(truth_allele, imgt)
                        if imgt_allele is None:
                            allele_writer.writerow({
                                "sample": sample,
                                "set": label,
                                "side": side,
                                "gene": gene,
                                "truth_allele": truth_allele,
                                "truth_2field": two_field(truth_allele),
                                "imgt_allele": "NA",
                                "imgt_match": match,
                            })
                            continue
                        cache_key = (gene, imgt_allele)
                        if cache_key not in truth_vcf_cache:
                            truth_vcf = build_truth_vcf(
                                gene_ref, gene, imgt_allele, imgt[imgt_allele], args.cache_dir, args.force
                            )
                            truth_keys, _truth_gt, truth_types = read_vcf_keys(truth_vcf)
                            truth_vcf_cache[cache_key] = (truth_keys, truth_types)
                        truth_keys, truth_types = truth_vcf_cache[cache_key]
                        covered_keys = coverage_filter(truth_keys, depth_arr, args.min_depth)
                        sample_gene_truth.update(truth_keys)
                        sample_gene_types.update(truth_types)

                        row = {
                            "sample": sample,
                            "set": label,
                            "side": side,
                            "gene": gene,
                            "truth_allele": truth_allele,
                            "truth_2field": two_field(truth_allele),
                            "imgt_allele": imgt_allele,
                            "imgt_match": match,
                            "truth_variants": len(truth_keys),
                            "truth_snvs": sum(1 for key in truth_keys if truth_types.get(key) == "SNV"),
                            "truth_indels": sum(1 for key in truth_keys if truth_types.get(key) == "INDEL"),
                            "covered_variants": len(covered_keys),
                            "covered_snvs": sum(1 for key in covered_keys if truth_types.get(key) == "SNV"),
                            "covered_indels": sum(1 for key in covered_keys if truth_types.get(key) == "INDEL"),
                        }
                        for stage in STAGES:
                            present_keys, gt_keys, _types = pipeline_calls[stage]
                            if stage == "freebayes":
                                row["freebayes_present"] = count_hits(truth_keys, present_keys)
                                row["freebayes_present_covered"] = count_hits(covered_keys, present_keys)
                            row[f"{stage}_gt"] = count_hits(truth_keys, gt_keys)
                            row[f"{stage}_gt_covered"] = count_hits(covered_keys, gt_keys)
                        allele_writer.writerow(row)

                covered_sample_gene = coverage_filter(sample_gene_truth, depth_arr, args.min_depth)
                sample_gene_row = {
                    "sample": sample,
                    "set": label,
                    "gene": gene,
                    "truth_variants": len(sample_gene_truth),
                    "truth_snvs": sum(1 for key in sample_gene_truth if sample_gene_types.get(key) == "SNV"),
                    "truth_indels": sum(1 for key in sample_gene_truth if sample_gene_types.get(key) == "INDEL"),
                    "covered_variants": len(covered_sample_gene),
                    "covered_snvs": sum(1 for key in covered_sample_gene if sample_gene_types.get(key) == "SNV"),
                    "covered_indels": sum(1 for key in covered_sample_gene if sample_gene_types.get(key) == "INDEL"),
                }
                for stage in STAGES:
                    present_keys, gt_keys, _types = pipeline_calls[stage]
                    if stage == "freebayes":
                        sample_gene_row["freebayes_present"] = count_hits(sample_gene_truth, present_keys)
                        sample_gene_row["freebayes_present_covered"] = count_hits(covered_sample_gene, present_keys)
                    sample_gene_row[f"{stage}_gt"] = count_hits(sample_gene_truth, gt_keys)
                    sample_gene_row[f"{stage}_gt_covered"] = count_hits(covered_sample_gene, gt_keys)
                sample_gene_writer.writerow(sample_gene_row)

    write_gene_summary(sample_gene_path, gene_summary_path)
    print(f"wrote {allele_path}")
    print(f"wrote {sample_gene_path}")
    print(f"wrote {gene_summary_path}")


if __name__ == "__main__":
    main()