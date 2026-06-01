#!/usr/bin/env python3
"""Diagnose raw-read exon support for HLA allele families.

This truth-evaluation helper asks whether alleles missed by the current
assembled-haplotype caller have direct support in per-gene FASTQs. It builds
2-field-family-specific k-mers from the SpecHLA exon FASTA and counts those
k-mers in raw per-gene reads.

Truth is used only to report ranks for known missed alleles; no production call
is made here.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


DEFAULT_GENES = ["HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]


def two_field(allele: str) -> str:
    allele = allele.strip().replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    rest = rest.rstrip("GP")
    fields = rest.split(":")
    return f"{gene}*{fields[0]}:{fields[1]}" if len(fields) >= 2 else f"{gene}*{fields[0]}"


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def iter_kmers(seq: str, k: int) -> Iterable[str]:
    seq = seq.upper()
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" not in kmer:
            yield kmer


def read_exon_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and seq:
                    records.append((name, "".join(seq).upper()))
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name and seq:
        records.append((name, "".join(seq).upper()))
    return records


def allele_from_header(header: str) -> str | None:
    match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", header)
    return match.group(1) if match else None


def build_unique_kmer_map(exon_fasta: Path, k: int) -> tuple[dict[str, str], Counter[str]]:
    kmer_families: dict[str, set[str]] = defaultdict(set)
    for header, seq in read_exon_records(exon_fasta):
        allele = allele_from_header(header)
        if not allele:
            continue
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


def count_fastq_support(paths: list[Path], unique_map: dict[str, str], k: int) -> tuple[Counter[str], Counter[str], int]:
    read_hits: Counter[str] = Counter()
    kmer_hits: Counter[str] = Counter()
    reads = 0
    for path in paths:
        if not path.exists():
            continue
        for seq in iter_fastq_sequences(path):
            reads += 1
            families_in_read: set[str] = set()
            for kmer in iter_kmers(seq, k):
                family = unique_map.get(kmer)
                if family is None:
                    continue
                families_in_read.add(family)
                kmer_hits[family] += 1
            for family in families_in_read:
                read_hits[family] += 1
    return read_hits, kmer_hits, reads


def load_targets(path: Path, genes: set[str]) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("gene") in genes]
    return rows


def rank_family(family: str, read_hits: Counter[str], unique_counts: Counter[str]) -> tuple[int | str, int, int, str]:
    ranked = sorted(
        ((fam, reads, unique_counts.get(fam, 0)) for fam, reads in read_hits.items()),
        key=lambda row: (-row[1], -row[2], row[0]),
    )
    rank = "NA"
    for idx, (fam, _reads, _unique) in enumerate(ranked, 1):
        if fam == family:
            rank = idx
            break
    top5 = ";".join(f"{fam}:{reads}" for fam, reads, _unique in ranked[:5])
    return rank, read_hits.get(family, 0), unique_counts.get(family, 0), top5


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--errors", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--exon-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--k", type=int, default=31)
    args = parser.parse_args()

    genes = set(args.genes)
    targets = load_targets(args.errors, genes)
    by_gene_sample: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in targets:
        by_gene_sample[(row["gene"], row["sample"])].append(row)

    kmer_maps: dict[str, tuple[dict[str, str], Counter[str]]] = {}
    for gene in genes:
        exon_fasta = args.exon_dir / f"{gene.replace('-', '_')}.fasta"
        kmer_maps[gene] = build_unique_kmer_map(exon_fasta, args.k)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "set", "sample", "gene", "side", "pos", "truth", "truth_2field",
        "truth_in_pool_2f", "score2", "em_rank_truth", "fastq_reads",
        "exon_unique_kmers", "truth_support_reads", "truth_exon_rank", "top5_exon_support",
    ]
    with args.out.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for (gene, sample), rows in sorted(by_gene_sample.items()):
            unique_map, unique_counts = kmer_maps[gene]
            short = gene.split("-", 1)[1]
            sample_dir = args.spechla_root / sample
            fastqs = [sample_dir / f"{short}.R1.fq.gz", sample_dir / f"{short}.R2.fq.gz"]
            read_hits, _kmer_hits, fastq_reads = count_fastq_support(fastqs, unique_map, args.k)
            for row in rows:
                truth_2field = two_field(row["truth"])
                rank, support_reads, unique_kmers, top5 = rank_family(truth_2field, read_hits, unique_counts)
                writer.writerow({
                    "set": row.get("set", ""),
                    "sample": sample,
                    "gene": gene,
                    "side": row.get("side", ""),
                    "pos": row.get("pos", ""),
                    "truth": row.get("truth", ""),
                    "truth_2field": truth_2field,
                    "truth_in_pool_2f": row.get("truth_in_pool_2f", ""),
                    "score2": row.get("score2", ""),
                    "em_rank_truth": row.get("em_rank_truth", ""),
                    "fastq_reads": fastq_reads,
                    "exon_unique_kmers": unique_kmers,
                    "truth_support_reads": support_reads,
                    "truth_exon_rank": rank,
                    "top5_exon_support": top5,
                })
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
