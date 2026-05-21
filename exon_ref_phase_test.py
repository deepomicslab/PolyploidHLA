#!/usr/bin/env python3
"""Experimental exon-reference variant calling and phasing for class-II HLA.

This script is intentionally separate from the main PolyploidHLA workflow. It
keeps class-I calls unchanged by copying the existing per-gene outputs, then
for class-II genes it remaps binned per-gene reads to an exon-DB allele chosen
as the local reference, calls variants, phases the exon contig, builds four
exon haplotype sequences, and types each haplotype against the same exon DB.

The output is an asm-root-like directory that can be evaluated with the usual
aggregate/evaluate scripts:

  <out-root>/<sample>/<gene-lower>/<gene>/calls.tsv
  <out-root>/<sample>/<gene-lower>/<gene>/hap{1..4}.fa
  <out-root>/<sample>/<sample>.final_calls.tsv

Truth is never used by this script.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BUNDLED_SPECHLA = SCRIPT_DIR / "resources" / "spechla"
DEFAULT_LEGACY_SPECHLA = SCRIPT_DIR.parent / "SpecHLA"
DEFAULT_SPECHLA = Path(
    os.environ.get(
        "SPECHLA",
        DEFAULT_BUNDLED_SPECHLA if DEFAULT_BUNDLED_SPECHLA.exists() else DEFAULT_LEGACY_SPECHLA,
    )
)
DEFAULT_EXON_DB = DEFAULT_SPECHLA / "db" / "HLA" / "exon"
DEFAULT_G_GROUP = DEFAULT_SPECHLA / "db" / "HLA" / "hla_nom_g.txt"
DEFAULT_CLASS_I = ["HLA-A", "HLA-B", "HLA-C"]
DEFAULT_CLASS_II = ["HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]


def run(cmd: list[str], log_path: Optional[Path] = None, timeout: Optional[int] = None) -> None:
    if log_path:
        with log_path.open("a") as log:
            log.write("$ " + " ".join(cmd) + "\n")
            try:
                proc = subprocess.run(
                    cmd,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=timeout,
                )
            except subprocess.TimeoutExpired as exc:
                log.write(f"[warn] command timed out after {timeout}s\n")
                raise RuntimeError(f"command timed out after {timeout}s: {' '.join(cmd)}") from exc
    else:
        try:
            proc = subprocess.run(cmd, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(f"command timed out after {timeout}s: {' '.join(cmd)}") from exc
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")


def run_pipe(cmds: list[list[str]], out_path: Path, log_path: Path) -> None:
    with log_path.open("a") as log:
        log.write("$ " + " | ".join(" ".join(cmd) for cmd in cmds) + f" > {out_path}\n")
        procs = []
        prev_stdout = None
        for idx, cmd in enumerate(cmds):
            stdout = subprocess.PIPE if idx < len(cmds) - 1 else out_path.open("wb")
            proc = subprocess.Popen(
                cmd,
                stdin=prev_stdout,
                stdout=stdout,
                stderr=log,
            )
            if prev_stdout is not None:
                prev_stdout.close()
            if idx < len(cmds) - 1:
                prev_stdout = proc.stdout
            procs.append((proc, stdout))
        for proc, stdout in procs:
            rc = proc.wait()
            if stdout is not subprocess.PIPE and hasattr(stdout, "close"):
                stdout.close()
            if rc != 0:
                raise RuntimeError(f"pipeline command failed ({rc}): {' '.join(proc.args)}")


def read_fasta_records(path: Path) -> list[tuple[str, str, str]]:
    records = []
    header = ""
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header and seq_parts:
                    allele = allele_from_header(header)
                    if allele:
                        records.append((header, allele, "".join(seq_parts).upper().replace("-", "")))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header and seq_parts:
        allele = allele_from_header(header)
        if allele:
            records.append((header, allele, "".join(seq_parts).upper().replace("-", "")))
    return records


def allele_from_header(header: str) -> Optional[str]:
    match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", header)
    return match.group(1) if match else None


def strip_expr_suffix(field: str) -> str:
    return field[:-1] if field and field[-1].isalpha() and field[-1] != "G" else field


def clean_allele(allele: str) -> str:
    if not allele or allele == "NA" or "*" not in allele:
        return allele or "NA"
    gene, rest = allele.replace("HLA-", "").split("*", 1)
    fields = rest.split(":")
    fields[-1] = strip_expr_suffix(fields[-1])
    return f"{gene}*{':'.join(fields)}"


def two_field(allele: str) -> str:
    allele = clean_allele(allele)
    if allele == "NA" or "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    fields = rest.replace("G", "").split(":")
    return f"{gene}*{':'.join(fields[:2])}" if len(fields) >= 2 else f"{gene}*{fields[0]}"


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def iter_kmers(seq: str, k: int) -> Iterable[str]:
    seq = seq.upper()
    for idx in range(0, len(seq) - k + 1):
        kmer = seq[idx:idx + k]
        if "N" not in kmer:
            yield kmer


def iter_fastq_sequences(path: Path) -> Iterable[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            if seq:
                yield seq


def build_family_unique_kmers(records: list[tuple[str, str, str]], k: int) -> tuple[dict[str, str], Counter[str]]:
    kmer_families: dict[str, set[str]] = defaultdict(set)
    for _header, allele, seq in records:
        family = two_field(allele)
        seen = set(iter_kmers(seq, k)) | set(iter_kmers(reverse_complement(seq), k))
        for kmer in seen:
            kmer_families[kmer].add(family)
    unique_map: dict[str, str] = {}
    unique_counts: Counter[str] = Counter()
    for kmer, families in kmer_families.items():
        if len(families) == 1:
            family = next(iter(families))
            unique_map[kmer] = family
            unique_counts[family] += 1
    return unique_map, unique_counts


def count_family_read_support(fastqs: list[Path], unique_map: dict[str, str], k: int) -> Counter[str]:
    read_hits: Counter[str] = Counter()
    for fastq in fastqs:
        if not fastq.exists():
            continue
        for seq in iter_fastq_sequences(fastq):
            families = {unique_map[kmer] for kmer in iter_kmers(seq, k) if kmer in unique_map}
            for family in families:
                read_hits[family] += 1
    return read_hits


def choose_exon_reference(records: list[tuple[str, str, str]], fastqs: list[Path], mode: str, k: int):
    if not records:
        raise ValueError("empty exon DB records")
    if mode == "longest":
        header, allele, seq = max(records, key=lambda row: (len(row[2]), row[1]))
        return header, allele, seq, "longest", 0
    unique_map, _unique_counts = build_family_unique_kmers(records, k)
    read_hits = count_family_read_support(fastqs, unique_map, k)
    if not read_hits:
        header, allele, seq = max(records, key=lambda row: (len(row[2]), row[1]))
        return header, allele, seq, "top-read-support-fallback-longest", 0
    top_family, support = max(read_hits.items(), key=lambda row: (row[1], row[0]))
    family_records = [row for row in records if two_field(row[1]) == top_family]
    header, allele, seq = max(family_records, key=lambda row: (len(row[2]), row[1]))
    return header, allele, seq, f"top-read-support:{top_family}", support


def write_reference(path: Path, contig: str, allele: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(f">{contig} allele={allele}\n")
        for idx in range(0, len(seq), 80):
            handle.write(seq[idx:idx + 80] + "\n")


def index_reference(ref_fa: Path, bwa: str, samtools: str, log: Path) -> None:
    if not (ref_fa.with_suffix(ref_fa.suffix + ".bwt")).exists():
        run([bwa, "index", str(ref_fa)], log)
    run([samtools, "faidx", str(ref_fa)], log)


def copy_indexed_vcf(src_vcf: Path, dst_vcf: Path) -> None:
    dst_tbi = Path(str(dst_vcf) + ".tbi")
    for path in (dst_vcf, dst_tbi):
        if path.exists():
            path.unlink()
    shutil.copy2(src_vcf, dst_vcf)
    shutil.copy2(Path(str(src_vcf) + ".tbi"), dst_tbi)


def parse_chi(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    patterns = [r"^GLOBAL\s+.*?chi_R=([0-9.]+)", r"chi_R=([0-9.]+)"]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            try:
                value = float(match.group(1))
            except ValueError:
                continue
            if 0.0 < value < 1.0:
                return value
    return None


def phase_assignment_from_vcf(vcf_path: Path, chi_r: Optional[float]) -> dict[int, str]:
    if chi_r is None or not vcf_path.exists():
        return {1: "R", 2: "R", 3: "D", 4: "D"}
    obs = []
    try:
        vcf = pysam.VariantFile(str(vcf_path))
        for rec in vcf:
            sample = rec.samples[0]
            gt = sample.get("GT")
            if gt is None or len(gt) != 4 or any(value is None for value in gt):
                continue
            if not rec.alts or len(rec.alts) != 1 or not set(gt) <= {0, 1}:
                continue
            ad = sample.get("AD")
            if not ad or len(ad) < 2 or any(value is None for value in ad) or sum(ad) < 10:
                continue
            obs.append((tuple(gt), ad[1] / sum(ad)))
    except (OSError, ValueError):
        obs = []
    if len(obs) < 3:
        return {1: "R", 2: "R", 3: "D", 4: "D"}
    best = None
    for recip_zero_based in itertools.combinations(range(4), 2):
        recip = set(recip_zero_based)
        weights = [chi_r / 2.0 if idx in recip else (1.0 - chi_r) / 2.0 for idx in range(4)]
        residual = 0.0
        for gt, af in obs:
            exp = sum(weights[idx] * gt[idx] for idx in range(4))
            residual += abs(af - exp)
        if best is None or residual < best[0]:
            best = (residual, recip)
    recip = best[1] if best else {0, 1}
    return {idx + 1: ("R" if idx in recip else "D") for idx in range(4)}


def fasta_sequence(path: Path) -> str:
    seq = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def coverage_mask(bam_path: Path, contig: str, length: int, min_depth: int) -> list[bool]:
    if min_depth <= 0:
        return [False] * length
    cov = [0] * length
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for col in bam.pileup(contig, 0, length, truncate=True, stepper="nofilter"):
        depth = 0
        for pileup_read in col.pileups:
            read = pileup_read.alignment
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            depth += 1
        if 0 <= col.reference_pos < length:
            cov[col.reference_pos] = depth
    bam.close()
    return [depth < min_depth for depth in cov]


def consensus_hap(ref_fa: Path, vcf_gz: Path, contig: str, hap_idx: int, sample: str,
                  fallback_seq: str, bcftools: str, samtools: str, log: Path) -> str:
    if not vcf_gz.exists() or not Path(str(vcf_gz) + ".tbi").exists():
        return fallback_seq
    region = f"{contig}:1-{len(fallback_seq)}"
    cmd1 = [samtools, "faidx", str(ref_fa), region]
    cmd2 = [bcftools, "consensus", "-H", str(hap_idx), "-s", sample, str(vcf_gz)]
    with log.open("a") as log_handle:
        log_handle.write("$ " + " ".join(cmd1) + " | " + " ".join(cmd2) + "\n")
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=log_handle)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=log_handle, text=True)
        p1.stdout.close()
        out, _err = p2.communicate()
        rc1 = p1.wait()
        if rc1 != 0 or p2.returncode != 0:
            log_handle.write("[warn] consensus failed; using reference sequence\n")
            return fallback_seq
    seq = "".join(line.strip().upper() for line in out.splitlines() if not line.startswith(">"))
    return seq or fallback_seq


def apply_mask(seq: str, mask: list[bool]) -> str:
    return "".join("N" if idx < len(mask) and mask[idx] else base for idx, base in enumerate(seq))


def build_exon_aligner(records: list[tuple[str, str, str]], prefilter_top: int, backend: str):
    try:
        from hla_polyphase_assemble import BaseLevelAligner
    except ImportError as exc:
        raise RuntimeError("cannot import BaseLevelAligner from hla_polyphase_assemble.py") from exc
    tmpdir = tempfile.TemporaryDirectory(prefix="exon_ref_type_")
    selected_backend = backend
    if selected_backend == "auto":
        try:
            import parasail  # noqa: F401
            selected_backend = "parasail"
        except ImportError:
            selected_backend = "mappy"
    aligner = BaseLevelAligner(backend=selected_backend, prefilter_top=prefilter_top)
    allele_seqs = {allele: record_seq for _header, allele, record_seq in records}
    aligner.index_alleles(allele_seqs, tmpdir.name)
    return aligner, tmpdir


def type_haplotype(seq: str, aligner) -> tuple[str, float, int]:
    query = seq.replace("N", "")
    if not query:
        return "NA", 0.0, 0
    scores = aligner.score_against_all(query, top_k=1, use_prefilter=True)
    if not scores:
        return "NA", 0.0, 0
    allele, score = max(scores.items(), key=lambda row: row[1])
    return allele, float(score), len(query)


def call_and_phase_gene(args, sample: str, gene: str, sample_out: Path) -> None:
    short = gene.split("-", 1)[1]
    contig = gene.replace("-", "_")
    sample_spechla = args.spechla_root / sample
    fastqs = [sample_spechla / f"{short}.R1.fq.gz", sample_spechla / f"{short}.R2.fq.gz"]
    if not all(path.exists() and path.stat().st_size > 0 for path in fastqs):
        raise FileNotFoundError(f"missing binned FASTQs for {sample} {gene}: {fastqs}")
    gene_dir = sample_out / gene.lower() / gene
    work_dir = gene_dir / "exon_ref_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    log = work_dir / "exon_ref_phase.log"

    records = read_fasta_records(args.exon_dir / f"{gene.replace('-', '_')}.fasta")
    _header, ref_allele, ref_seq, ref_reason, ref_support = choose_exon_reference(
        records, fastqs, args.ref_choice, args.ref_k
    )
    ref_fa = work_dir / f"{contig}.exon_ref.fa"
    write_reference(ref_fa, contig, ref_allele, ref_seq)
    index_reference(ref_fa, args.bwa, args.samtools, log)

    bam = work_dir / f"{gene}.exon_ref.bam"
    if args.force or not bam.exists():
        read_group = f"@RG\tID:{sample}\tSM:{sample}"
        cmds = [
            [
                args.bwa, "mem", "-t", str(args.threads), "-U", "10000",
                "-L", "10000,10000", "-R", read_group,
                str(ref_fa), str(fastqs[0]), str(fastqs[1]),
            ],
            [args.samtools, "view", "-@", str(args.samtools_threads), "-bS", "-F", "0x800", "-"],
            [args.samtools, "sort", "-@", str(args.samtools_threads), "-o", str(bam), "-"],
        ]
        run_pipe(cmds, Path(os.devnull), log)
    run([args.samtools, "index", "-@", str(args.samtools_threads), str(bam)], log)

    raw_vcf = work_dir / f"{gene}.freebayes.raw.vcf.gz"
    vcf = work_dir / f"{gene}.freebayes.vcf.gz"
    vcf_tbi = Path(str(vcf) + ".tbi")
    if args.force or not vcf.exists() or not vcf_tbi.exists():
        for partial in (raw_vcf, vcf, vcf_tbi):
            if partial.exists():
                partial.unlink()
        cmds = [
            [
                args.freebayes,
                "-p", str(args.ploidy),
                "--min-alternate-fraction", str(args.min_af),
                "--min-alternate-count", str(args.min_ac),
                "--min-base-quality", str(args.min_bq),
                "--min-mapping-quality", str(args.min_mq),
                "--min-coverage", str(args.min_cov),
                "--haplotype-length", "0",
                "--use-best-n-alleles", "4",
                "-f", str(ref_fa), str(bam),
            ],
            [args.bcftools, "norm", "-f", str(ref_fa), "-a", "-m", "-any", "-Oz", "-o", str(vcf)],
        ]
        run_pipe(cmds, raw_vcf, log)
        if raw_vcf.exists():
            raw_vcf.unlink()
    run([args.tabix, "-f", "-p", "vcf", str(vcf)], log)

    phased_vcf = work_dir / f"{gene}.phased.vcf.gz"
    phased_tbi = Path(str(phased_vcf) + ".tbi")
    if args.skip_whatshap:
        for partial in (phased_vcf, phased_tbi):
            if partial.exists():
                partial.unlink()
        with log.open("a") as handle:
            handle.write("[warn] --skip-whatshap set; using unphased freebayes VCF fallback\n")
        copy_indexed_vcf(vcf, phased_vcf)
    elif args.force or not phased_vcf.exists() or not phased_tbi.exists():
        for partial in (phased_vcf, phased_tbi):
            if partial.exists():
                partial.unlink()
        try:
            run([
                args.whatshap, "polyphase",
                "--ploidy", str(args.ploidy),
                "--reference", str(ref_fa),
                "--threads", str(args.threads),
                "--ignore-read-groups",
                "--output", str(phased_vcf),
                str(vcf), str(bam),
            ], log, timeout=args.whatshap_timeout if args.whatshap_timeout > 0 else None)
            run([args.tabix, "-f", "-p", "vcf", str(phased_vcf)], log)
        except RuntimeError:
            with log.open("a") as handle:
                handle.write("[warn] whatshap failed or timed out; using unphased freebayes VCF fallback\n")
            copy_indexed_vcf(vcf, phased_vcf)

    chi = parse_chi(sample_spechla / f"{sample}.chi_pooled.txt")
    if chi is None:
        chi = parse_chi(sample_spechla / f"{sample}.chimerism.txt")
    assignment = phase_assignment_from_vcf(phased_vcf, chi)
    mask = coverage_mask(bam, contig, len(ref_seq), args.mask_min_depth)
    aligner, aligner_tmpdir = build_exon_aligner(records, args.typing_prefilter_top, args.typing_aligner)

    try:
        with (gene_dir / "calls.tsv").open("w") as calls_out:
            calls_out.write("global_hap\tassignment\tallele\ttotal_assembly_score\n")
            for hap_idx in range(1, 5):
                seq = consensus_hap(ref_fa, phased_vcf, contig, hap_idx, sample, ref_seq,
                                    args.bcftools, args.samtools, log)
                seq = apply_mask(seq, mask)
                allele, score, _aligned_len = type_haplotype(seq, aligner)
                side = assignment[hap_idx]
                with (gene_dir / f"hap{hap_idx}.fa").open("w") as hap_out:
                    hap_out.write(
                        f">{gene}_hap{hap_idx} assignment={side} allele={allele} "
                        f"ref={ref_allele} reason={ref_reason} chi_R={chi if chi is not None else 'NA'}\n"
                    )
                    for offset in range(0, len(seq), 80):
                        hap_out.write(seq[offset:offset + 80] + "\n")
                calls_out.write(f"{hap_idx}\t{side}\t{allele}\t{score:.2f}\n")
    finally:
        aligner_tmpdir.cleanup()

    with (work_dir / "reference_choice.tsv").open("w") as handle:
        handle.write("sample\tgene\tref_choice\tref_allele\tref_len\tref_reason\tref_support\tchi_R\n")
        handle.write(
            f"{sample}\t{gene}\t{args.ref_choice}\t{ref_allele}\t{len(ref_seq)}\t"
            f"{ref_reason}\t{ref_support}\t{chi if chi is not None else 'NA'}\n"
        )


def copy_class_i_gene(baseline_asm_root: Path, sample: str, gene: str, sample_out: Path) -> None:
    src = baseline_asm_root / sample / gene.lower() / gene
    dst = sample_out / gene.lower() / gene
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def aggregate(sample_out_root: Path, sample: str, genes: list[str], args) -> Path:
    out = sample_out_root / sample / f"{sample}.final_calls.tsv"
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "aggregate_calls.py"),
        "--asm-root", str(sample_out_root),
        "--sample", sample,
        "--genes", *genes,
        "--out", str(out),
        "--g-group", str(args.g_group),
    ]
    run(cmd)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--spechla-root", required=True, type=Path,
                        help="root containing <sample>/<gene-short>.R1/R2.fq.gz after read binning")
    parser.add_argument("--baseline-asm-root", required=True, type=Path,
                        help="existing asm root used to keep class-I calls unchanged")
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--exon-dir", type=Path, default=DEFAULT_EXON_DB)
    parser.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP)
    parser.add_argument("--class-i-genes", nargs="+", default=DEFAULT_CLASS_I)
    parser.add_argument("--class-ii-genes", nargs="+", default=DEFAULT_CLASS_II)
    parser.add_argument("--genes", nargs="+", default=None,
                        help="optional subset of class-II genes to run")
    parser.add_argument("--ref-choice", choices=("longest", "top-read-support"), default="longest")
    parser.add_argument("--ref-k", type=int, default=51)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--samtools-threads", type=int, default=0)
    parser.add_argument("--ploidy", type=int, default=4)
    parser.add_argument("--min-af", type=float, default=0.01)
    parser.add_argument("--min-ac", type=int, default=2)
    parser.add_argument("--min-bq", type=int, default=13)
    parser.add_argument("--min-mq", type=int, default=10)
    parser.add_argument("--min-cov", type=int, default=5)
    parser.add_argument("--mask-min-depth", type=int, default=0)
    parser.add_argument("--typing-prefilter-top", type=int, default=200,
                        help="mappy prefilter top-K alleles before parasail exon typing; 0 disables")
    parser.add_argument("--typing-aligner", choices=("auto", "parasail", "mappy"), default="auto")
    parser.add_argument("--whatshap-timeout", type=int, default=300,
                        help="seconds allowed for whatshap polyphase per gene; <=0 disables timeout")
    parser.add_argument("--skip-whatshap", action="store_true",
                        help="skip polyphase and use the indexed freebayes VCF as a fallback input")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--bwa", default=os.environ.get("BWA", "bwa"))
    parser.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    parser.add_argument("--freebayes", default=os.environ.get("FREEBAYES", "freebayes"))
    parser.add_argument("--whatshap", default=os.environ.get("WHATSHAP", "whatshap"))
    parser.add_argument("--tabix", default=os.environ.get("TABIX", "tabix"))
    args = parser.parse_args()
    if args.samtools_threads <= 0:
        args.samtools_threads = args.threads
    return args


def main() -> None:
    args = parse_args()
    for tool in (args.bwa, args.samtools, args.bcftools, args.freebayes, args.whatshap, args.tabix):
        if shutil.which(tool) is None:
            raise SystemExit(f"required tool not found on PATH: {tool}")

    sample_out = args.out_root / args.sample
    sample_out.mkdir(parents=True, exist_ok=True)
    class_ii = args.genes if args.genes else args.class_ii_genes
    all_genes = list(args.class_i_genes) + list(class_ii)

    for gene in args.class_i_genes:
        print(f"[class-I copy] {args.sample} {gene}", flush=True)
        copy_class_i_gene(args.baseline_asm_root, args.sample, gene, sample_out)

    for gene in class_ii:
        print(f"[exon-ref] {args.sample} {gene} ref_choice={args.ref_choice}", flush=True)
        call_and_phase_gene(args, args.sample, gene, sample_out)

    final_calls = aggregate(args.out_root, args.sample, all_genes, args)
    print(f"[done] wrote {final_calls}", flush=True)


if __name__ == "__main__":
    main()
