#!/usr/bin/env python3
"""Generate reproducible mixed-source HLA FASTQs with wgsim."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, NamedTuple, Sequence, Tuple


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BENCH_ROOT = DEFAULT_PROJECT_ROOT.parent / "PolyploidHLA_simulation"
DEFAULT_IMGT_FASTA = Path(
    "/data2/wangxuedong/compare_hla_db/IMGTHLA-3.62.0-alpha/hla_gen.fasta"
)
DEFAULT_CALLER_ALLELES = (
    DEFAULT_PROJECT_ROOT
    / "resources/spechla/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
)
DEFAULT_GENES = ("A", "B", "C", "DRB1", "DQB1", "DPB1")
COPY_LAYOUT = (("recipient", "R1"), ("recipient", "R2"), ("graft", "G1"), ("graft", "G2"))
EXPRESSION_SUFFIX = re.compile(r"[NLSCAQ]$")


class AlleleRecord(NamedTuple):
    sequence_id: str
    gene: str
    allele_full: str
    allele_2field: str
    sequence: str


class CopyTruth(NamedTuple):
    source: str
    haplotype: str
    allele: AlleleRecord


class FastqRecord(NamedTuple):
    name: str
    sequence: str
    plus: str
    quality: str


def stable_seed(master_seed: int, *parts: object) -> int:
    payload = "\x1f".join([str(master_seed)] + [str(part) for part in parts])
    value = int.from_bytes(hashlib.sha256(payload.encode("utf-8")).digest()[:8], "big")
    return value % 2_147_483_646 + 1


def parse_allele_header(header: str) -> Tuple[str, str, str, str]:
    fields = header.split()
    allele_index = next((index for index, field in enumerate(fields) if "*" in field), None)
    if allele_index is None:
        raise ValueError("not an IMGT/HLA genomic allele header")
    allele_token = fields[allele_index]
    gene, allele_fields = allele_token.split("*", 1)
    resolution = allele_fields.split(":")
    if len(resolution) < 2:
        raise ValueError("allele does not have 2-field resolution")
    allele_full = f"{gene}*{allele_fields}"
    allele_2field = f"{gene}*{resolution[0]}:{resolution[1]}"
    sequence_id = fields[0] if fields[0].startswith("HLA:") else allele_full
    return sequence_id, gene, allele_full, allele_2field


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    chunks: List[str] = []
    with path.open("rt", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"sequence before first header in {path}")
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def load_candidates(
    fasta: Path,
    genes: Sequence[str],
    min_length: int,
    max_ambiguous_fraction: float,
    callable_alleles: set,
) -> Dict[str, Dict[str, List[AlleleRecord]]]:
    requested = set(genes)
    candidates: Dict[str, Dict[str, List[AlleleRecord]]] = {
        gene: defaultdict(list) for gene in genes
    }
    for header, sequence in read_fasta(fasta):
        try:
            sequence_id, gene, allele_full, allele_2field = parse_allele_header(header)
        except ValueError:
            continue
        if (
            gene not in requested
            or allele_full not in callable_alleles
            or EXPRESSION_SUFFIX.search(allele_full)
            or len(sequence) < min_length
        ):
            continue
        if set(sequence) - set("ACGTN"):
            continue
        if sequence.count("N") / len(sequence) > max_ambiguous_fraction:
            continue
        candidates[gene][allele_2field].append(
            AlleleRecord(sequence_id, gene, allele_full, allele_2field, sequence)
        )
    missing = [gene for gene in genes if len(candidates[gene]) < 4]
    if missing:
        details = ", ".join(f"{gene}={len(candidates[gene])}" for gene in missing)
        raise ValueError(f"fewer than four eligible 2-field families: {details}")
    return candidates


def load_callable_alleles(path: Path) -> set:
    alleles = set()
    for header, _sequence in read_fasta(path):
        token = next((field for field in header.split() if "*" in field), None)
        if token:
            alleles.add(token)
    if not alleles:
        raise ValueError(f"no callable allele names found in {path}")
    return alleles


def choose_record(records: Sequence[AlleleRecord], rng: random.Random) -> AlleleRecord:
    return records[rng.randrange(len(records))]


def select_gene_truth(
    families: Dict[str, List[AlleleRecord]], scenario: str, rng: random.Random
) -> List[CopyTruth]:
    required = 4 if scenario == "distinct4" else 3
    selected_families = rng.sample(sorted(families), required)
    selected = [choose_record(families[family], rng) for family in selected_families]
    if scenario == "distinct4":
        alleles = selected
    elif scenario == "shared1":
        alleles = [selected[0], selected[1], selected[0], selected[2]]
    elif scenario == "homozygous_graft":
        alleles = [selected[0], selected[1], selected[2], selected[2]]
    else:
        raise ValueError(f"unsupported scenario: {scenario}")
    return [
        CopyTruth(source, haplotype, allele)
        for (source, haplotype), allele in zip(COPY_LAYOUT, alleles)
    ]


def select_genotypes(
    candidates: Dict[str, Dict[str, List[AlleleRecord]]],
    genes: Sequence[str],
    scenario: str,
    individuals: int,
    sample_start: int,
    master_seed: int,
) -> Dict[str, Dict[str, List[CopyTruth]]]:
    genotypes: Dict[str, Dict[str, List[CopyTruth]]] = {}
    for index in range(sample_start, sample_start + individuals):
        sample_id = f"SIM{index:04d}"
        genotypes[sample_id] = {}
        for gene in genes:
            rng = random.Random(stable_seed(master_seed, "genotype", sample_id, gene, scenario))
            genotypes[sample_id][gene] = select_gene_truth(candidates[gene], scenario, rng)
    return genotypes


def expected_fraction(copy_truth: CopyTruth, graft_fraction: float) -> float:
    if copy_truth.source == "graft":
        return graft_fraction / 2.0
    return (1.0 - graft_fraction) / 2.0


def read_pairs_for_copy(
    total_coverage: float, fraction: float, sequence_length: int, read_length: int
) -> int:
    return max(1, int(math.floor(total_coverage * fraction * sequence_length / (2 * read_length) + 0.5)))


def write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("wt", encoding="ascii") as handle:
        handle.write(f">{name}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset : offset + 80] + "\n")


def fastq_records(path: Path) -> Iterator[FastqRecord]:
    with path.open("rt", encoding="ascii") as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ record in {path}")
            yield FastqRecord(name.rstrip(), sequence.rstrip(), plus.rstrip(), quality.rstrip())


def normalized_read_name(name: str) -> str:
    token = name.split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def read_fastq_pairs(r1_path: Path, r2_path: Path) -> Iterator[Tuple[FastqRecord, FastqRecord]]:
    r1_iter = fastq_records(r1_path)
    r2_iter = fastq_records(r2_path)
    while True:
        try:
            r1 = next(r1_iter)
        except StopIteration:
            try:
                next(r2_iter)
            except StopIteration:
                return
            raise ValueError(f"R2 has more records than R1: {r2_path}")
        try:
            r2 = next(r2_iter)
        except StopIteration:
            raise ValueError(f"R1 has more records than R2: {r1_path}")
        if normalized_read_name(r1.name) != normalized_read_name(r2.name):
            raise ValueError(f"mate name mismatch: {r1.name} != {r2.name}")
        yield r1, r2


def run_wgsim(
    executable: str,
    fasta: Path,
    r1_path: Path,
    r2_path: Path,
    read_pairs: int,
    read_length: int,
    insert_mean: int,
    insert_sd: int,
    error_rate: float,
    seed: int,
    log_handle: object,
) -> None:
    command = [
        executable,
        "-N", str(read_pairs),
        "-1", str(read_length),
        "-2", str(read_length),
        "-d", str(insert_mean),
        "-s", str(insert_sd),
        "-e", str(error_rate),
        "-r", "0",
        "-R", "0",
        "-X", "0",
        "-S", str(seed),
        str(fasta),
        str(r1_path),
        str(r2_path),
    ]
    subprocess.run(command, check=True, stdout=log_handle, stderr=log_handle)


def write_merged_fastqs(
    inputs: Sequence[Tuple[Path, Path, str]],
    output_r1: Path,
    output_r2: Path,
    name_map_path: Path,
    sample_id: str,
    merge_seed: int,
) -> int:
    pairs: List[Tuple[FastqRecord, FastqRecord, str]] = []
    for r1_path, r2_path, copy_key in inputs:
        for r1, r2 in read_fastq_pairs(r1_path, r2_path):
            pairs.append((r1, r2, copy_key))
    random.Random(merge_seed).shuffle(pairs)
    output_r1.parent.mkdir(parents=True, exist_ok=True)
    name_map_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_r1, "wt", encoding="ascii") as r1_handle, gzip.open(
        output_r2, "wt", encoding="ascii"
    ) as r2_handle, gzip.open(name_map_path, "wt", encoding="ascii", newline="") as map_handle:
        map_writer = csv.writer(map_handle, delimiter="\t", lineterminator="\n")
        map_writer.writerow(("blinded_read_name", "original_read_name", "copy_key"))
        for index, (r1, r2, copy_key) in enumerate(pairs, start=1):
            blinded = f"{sample_id}_{index:09d}"
            r1_handle.write(f"@{blinded}/1\n{r1.sequence}\n+\n{r1.quality}\n")
            r2_handle.write(f"@{blinded}/2\n{r2.sequence}\n+\n{r2.quality}\n")
            map_writer.writerow((blinded, normalized_read_name(r1.name), copy_key))
    return len(pairs)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def condition_name(graft_fraction: float, coverage: float) -> str:
    graft_percent = int(round(graft_fraction * 100))
    coverage_label = f"{coverage:g}".replace(".", "p")
    return f"graft{graft_percent:02d}_cov{coverage_label}x"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench-root", type=Path, default=DEFAULT_BENCH_ROOT)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT_ROOT)
    parser.add_argument("--imgt-fasta", type=Path, default=DEFAULT_IMGT_FASTA)
    parser.add_argument("--caller-alleles", type=Path, default=DEFAULT_CALLER_ALLELES)
    parser.add_argument("--experiment", default="smoke")
    parser.add_argument("--scenario", choices=("distinct4", "shared1", "homozygous_graft"), default="distinct4")
    parser.add_argument("--genes", nargs="+", default=list(DEFAULT_GENES))
    parser.add_argument("--individuals", type=int, default=1)
    parser.add_argument("--sample-start", type=int, default=1)
    parser.add_argument("--graft-fractions", nargs="+", type=float, default=[0.10])
    parser.add_argument("--coverages", nargs="+", type=float, default=[50.0])
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--insert-mean", type=int, default=350)
    parser.add_argument("--insert-sd", type=int, default=50)
    parser.add_argument("--error-rate", type=float, default=0.001)
    parser.add_argument("--min-sequence-length", type=int, default=2500)
    parser.add_argument("--max-ambiguous-fraction", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=20260728)
    parser.add_argument("--wgsim", default=shutil.which("wgsim") or "wgsim")
    parser.add_argument("--dry-run", action="store_true", help="write design/truth manifests without FASTQs")
    parser.add_argument("--run-caller", action="store_true", help="run polyphase_v2.sh after each FASTQ pair")
    parser.add_argument(
        "--batch-results",
        type=Path,
        help="shared main-pipeline result TSV updated after each completed sample",
    )
    parser.add_argument("--keep-gene-reads", action="store_true")
    parser.add_argument(
        "--cleanup-caller-intermediates",
        action="store_true",
        help="after a successful caller run, remove regenerable FASTQ/BAM and augmented-index files",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.individuals < 1:
        raise ValueError("--individuals must be at least 1")
    if args.sample_start < 1:
        raise ValueError("--sample-start must be at least 1")
    if not args.imgt_fasta.is_file():
        raise FileNotFoundError(f"missing IMGT/HLA genomic FASTA: {args.imgt_fasta}")
    if not args.caller_alleles.is_file():
        raise FileNotFoundError(f"missing production caller allele FASTA: {args.caller_alleles}")
    if args.read_length < 1 or args.insert_mean < 2 * args.read_length:
        raise ValueError("insert mean must be at least twice the read length")
    if args.insert_sd < 0 or not 0 <= args.error_rate < 1:
        raise ValueError("invalid insert SD or error rate")
    if not 0 <= args.max_ambiguous_fraction <= 1:
        raise ValueError("--max-ambiguous-fraction must be in [0, 1]")
    for fraction in args.graft_fractions:
        if not 0.10 <= fraction <= 0.50:
            raise ValueError("graft fractions must be in [0.10, 0.50]")
    if any(coverage <= 0 for coverage in args.coverages):
        raise ValueError("coverages must be positive")
    if args.run_caller and args.dry_run:
        raise ValueError("--run-caller cannot be combined with --dry-run")
    if args.run_caller and not (args.project_root / "polyphase_v2.sh").is_file():
        raise FileNotFoundError(f"missing caller: {args.project_root / 'polyphase_v2.sh'}")
    if args.run_caller:
        required = ("python", "whatshap", "freebayes", "bowtie2", "bwa", "samtools", "bcftools")
        missing = [executable for executable in required if shutil.which(executable) is None]
        if missing:
            raise FileNotFoundError(
                "missing caller executables on PATH: " + ", ".join(missing)
            )
    if not args.dry_run and shutil.which(args.wgsim) is None and not Path(args.wgsim).is_file():
        raise FileNotFoundError(f"wgsim executable not found: {args.wgsim}")


def run_caller(
    args: argparse.Namespace,
    sample_id: str,
    condition: str,
    graft_fraction: float,
    coverage: float,
    r1: Path,
    r2: Path,
) -> int:
    work_dir = args.bench_root / "runs" / args.experiment / condition / sample_id
    log_path = args.bench_root / "logs" / args.experiment / condition / f"{sample_id}.polyphase.log"
    work_dir.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(
        {
            "FQ1": str(r1),
            "FQ2": str(r2),
            "SAMPLE": sample_id,
            "RECIPIENT_MAJOR": "1",
            "FB_MIN_AF": "0.03",
            "GT_DROP_FP_AF": "0.05",
            "CLASS2_DPB1_RARE_COLLAPSE": "0",
            "WORK_DIR": str(work_dir),
            "RESULT_EXPERIMENT": args.experiment,
            "RESULT_CONDITION": condition,
            "RESULT_SCENARIO": args.scenario,
            "RESULT_GRAFT_FRACTION": f"{graft_fraction:.6f}",
            "RESULT_TOTAL_COVERAGE": f"{coverage:g}",
            "RESULT_READ_LENGTH": str(args.read_length),
            "RESULT_INSERT_MEAN": str(args.insert_mean),
            "RESULT_INSERT_SD": str(args.insert_sd),
            "RESULT_ERROR_RATE": str(args.error_rate),
            "RESULT_MASTER_SEED": str(args.seed),
        }
    )
    if args.batch_results is not None:
        environment["BATCH_RESULTS_FILE"] = str(args.batch_results)
    with log_path.open("wt", encoding="utf-8") as log_handle:
        completed = subprocess.run(
            ["bash", str(args.project_root / "polyphase_v2.sh")],
            check=False,
            cwd=str(args.project_root),
            env=environment,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )
    return completed.returncode


def cleanup_caller_intermediates(work_dir: Path, sample_id: str, r1: Path, r2: Path) -> int:
    sample_root = work_dir / "spechla_out" / sample_id
    patterns = (
        "*.fq.gz",
        "*.bam",
        "*.bam.bai",
        "em_refine/*.aug.fa*",
        "drb345/*.aug.fa*",
    )
    removed = 0
    for pattern in patterns:
        for path in sample_root.glob(pattern):
            if path.is_file():
                path.unlink()
                removed += 1
    for path in (r1, r2):
        if path.is_file():
            path.unlink()
            removed += 1
    return removed


def main(argv: Sequence[str] = None) -> int:
    args = build_parser().parse_args(argv)
    validate_args(args)
    genes = tuple(
        dict.fromkeys(gene[4:] if gene.startswith("HLA-") else gene for gene in args.genes)
    )
    callable_alleles = load_callable_alleles(args.caller_alleles)
    candidates = load_candidates(
        args.imgt_fasta,
        genes,
        args.min_sequence_length,
        args.max_ambiguous_fraction,
        callable_alleles,
    )
    genotypes = select_genotypes(
        candidates, genes, args.scenario, args.individuals, args.sample_start, args.seed
    )
    experiment_truth = args.bench_root / "truth" / args.experiment
    config_root = args.bench_root / "config" / args.experiment
    config_root.mkdir(parents=True, exist_ok=True)
    experiment_truth.mkdir(parents=True, exist_ok=True)

    candidate_rows = [
        {
            "gene": gene,
            "eligible_2field_families": len(candidates[gene]),
            "eligible_full_records": sum(len(records) for records in candidates[gene].values()),
        }
        for gene in genes
    ]
    write_tsv(
        config_root / "candidate_counts.tsv",
        ("gene", "eligible_2field_families", "eligible_full_records"),
        candidate_rows,
    )
    design = {
        "experiment": args.experiment,
        "scenario": args.scenario,
        "individuals": args.individuals,
        "sample_start": args.sample_start,
        "genes": genes,
        "graft_fractions": args.graft_fractions,
        "coverages": args.coverages,
        "read_length": args.read_length,
        "insert_mean": args.insert_mean,
        "insert_sd": args.insert_sd,
        "error_rate": args.error_rate,
        "min_sequence_length": args.min_sequence_length,
        "max_ambiguous_fraction": args.max_ambiguous_fraction,
        "master_seed": args.seed,
        "imgt_fasta": str(args.imgt_fasta.resolve()),
        "imgt_fasta_sha256": hashlib.sha256(args.imgt_fasta.read_bytes()).hexdigest(),
        "caller_alleles": str(args.caller_alleles.resolve()),
        "caller_alleles_sha256": hashlib.sha256(args.caller_alleles.read_bytes()).hexdigest(),
        "wgsim": args.wgsim,
        "dry_run": args.dry_run,
    }
    with (config_root / "design.json").open("wt", encoding="utf-8") as handle:
        json.dump(design, handle, indent=2, sort_keys=True)
        handle.write("\n")

    sample_rows: List[Dict[str, object]] = []
    copy_rows: List[Dict[str, object]] = []
    caller_failures = 0
    for graft_fraction in args.graft_fractions:
        for coverage in args.coverages:
            condition = condition_name(graft_fraction, coverage)
            for sample_id, gene_truth in genotypes.items():
                final_r1 = args.bench_root / "reads" / args.experiment / condition / f"{sample_id}.R1.fastq.gz"
                final_r2 = args.bench_root / "reads" / args.experiment / condition / f"{sample_id}.R2.fastq.gz"
                merge_seed = stable_seed(args.seed, "merge", args.experiment, condition, sample_id)
                sample_row = {
                    "sample_id": sample_id,
                    "experiment": args.experiment,
                    "condition": condition,
                    "scenario": args.scenario,
                    "graft_fraction": f"{graft_fraction:.6f}",
                    "total_coverage": f"{coverage:g}",
                    "read_length": args.read_length,
                    "insert_mean": args.insert_mean,
                    "insert_sd": args.insert_sd,
                    "error_rate": args.error_rate,
                    "merge_seed": merge_seed,
                    "recipient_major": 1,
                    "fq1": str(final_r1),
                    "fq2": str(final_r2),
                    "read_pairs": "",
                    "caller_exit_code": "",
                    "caller_intermediates_cleaned": "",
                    "status": "designed" if args.dry_run else "pending",
                }
                inputs: List[Tuple[Path, Path, str]] = []
                for gene in genes:
                    for truth in gene_truth[gene]:
                        fraction = expected_fraction(truth, graft_fraction)
                        pair_count = read_pairs_for_copy(
                            coverage, fraction, len(truth.allele.sequence), args.read_length
                        )
                        copy_seed = stable_seed(
                            args.seed, "wgsim", args.experiment, condition, sample_id, gene, truth.haplotype
                        )
                        copy_key = f"{sample_id}|{gene}|{truth.source}|{truth.haplotype}"
                        copy_rows.append(
                            {
                                "sample_id": sample_id,
                                "experiment": args.experiment,
                                "condition": condition,
                                "gene": f"HLA-{gene}",
                                "source": truth.source,
                                "haplotype": truth.haplotype,
                                "allele_full": truth.allele.allele_full,
                                "allele_2field": truth.allele.allele_2field,
                                "sequence_id": truth.allele.sequence_id,
                                "sequence_length": len(truth.allele.sequence),
                                "expected_fraction": f"{fraction:.6f}",
                                "expected_depth": f"{coverage * fraction:.6f}",
                                "read_pairs": pair_count,
                                "copy_seed": copy_seed,
                            }
                        )
                        if args.dry_run:
                            continue
                        copy_root = args.bench_root / "gene_reads" / args.experiment / condition / sample_id / gene / truth.haplotype
                        fasta_path = copy_root / "truth.fa"
                        r1_path = copy_root / "copy.R1.fastq"
                        r2_path = copy_root / "copy.R2.fastq"
                        if copy_root.exists() and not args.overwrite:
                            raise FileExistsError(f"temporary output already exists: {copy_root}")
                        copy_root.mkdir(parents=True, exist_ok=True)
                        write_fasta(fasta_path, "truth_copy", truth.allele.sequence)
                        log_path = args.bench_root / "logs" / args.experiment / condition / f"{sample_id}.wgsim.log"
                        log_path.parent.mkdir(parents=True, exist_ok=True)
                        with log_path.open("at", encoding="utf-8") as log_handle:
                            run_wgsim(
                                args.wgsim,
                                fasta_path,
                                r1_path,
                                r2_path,
                                pair_count,
                                args.read_length,
                                args.insert_mean,
                                args.insert_sd,
                                args.error_rate,
                                copy_seed,
                                log_handle,
                            )
                        inputs.append((r1_path, r2_path, copy_key))
                if not args.dry_run:
                    if (final_r1.exists() or final_r2.exists()) and not args.overwrite:
                        raise FileExistsError(f"final FASTQ already exists for {sample_id} {condition}")
                    name_map = experiment_truth / "read_name_maps" / condition / f"{sample_id}.tsv.gz"
                    pair_total = write_merged_fastqs(
                        inputs, final_r1, final_r2, name_map, sample_id, merge_seed
                    )
                    sample_row["read_pairs"] = pair_total
                    sample_row["status"] = "fastq_complete"
                    if not args.keep_gene_reads:
                        shutil.rmtree(
                            args.bench_root / "gene_reads" / args.experiment / condition / sample_id
                        )
                    if args.run_caller:
                        exit_code = run_caller(
                            args, sample_id, condition, graft_fraction, coverage, final_r1, final_r2
                        )
                        sample_row["caller_exit_code"] = exit_code
                        if exit_code == 0:
                            sample_row["status"] = "caller_complete"
                            if args.cleanup_caller_intermediates:
                                removed = cleanup_caller_intermediates(
                                    args.bench_root / "runs" / args.experiment / condition / sample_id,
                                    sample_id,
                                    final_r1,
                                    final_r2,
                                )
                                sample_row["caller_intermediates_cleaned"] = removed
                        else:
                            sample_row["status"] = "caller_failed"
                            caller_failures += 1
                sample_rows.append(sample_row)

    sample_fields = (
        "sample_id", "experiment", "condition", "scenario", "graft_fraction",
        "total_coverage", "read_length", "insert_mean", "insert_sd", "error_rate",
        "merge_seed", "recipient_major", "fq1", "fq2", "read_pairs",
        "caller_exit_code", "caller_intermediates_cleaned", "status",
    )
    copy_fields = (
        "sample_id", "experiment", "condition", "gene", "source", "haplotype",
        "allele_full", "allele_2field", "sequence_id", "sequence_length",
        "expected_fraction", "expected_depth", "read_pairs", "copy_seed",
    )
    write_tsv(config_root / "samples.tsv", sample_fields, sample_rows)
    write_tsv(experiment_truth / "copies.tsv", copy_fields, copy_rows)
    print(f"[simulation] experiment={args.experiment} samples={len(sample_rows)} dry_run={args.dry_run}")
    print(f"[simulation] config={config_root}")
    if not args.dry_run:
        print(f"[simulation] reads={args.bench_root / 'reads' / args.experiment}")
    if caller_failures:
        print(f"[simulation] ERROR: {caller_failures} caller run(s) failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (FileNotFoundError, FileExistsError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[simulation] ERROR: {error}", file=sys.stderr)
        sys.exit(2)