#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Iterator, Optional, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_EXON_DIR = SCRIPT_DIR / "resources" / "spechla" / "db" / "HLA" / "exon"
DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]


def normalize_read_name(header: str) -> str:
    name = header.strip().split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def normalize_allele(allele: str) -> str:
    allele = allele.strip().replace("HLA-", "").replace("G", "").rstrip("P")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    parts = rest.split(":")
    if parts and parts[-1].isalpha():
        parts[-1] = parts[-1][:-1]
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def allele_number(allele: str) -> int:
    match = re.search(r"\*(\d+)", allele)
    return int(match.group(1)) if match else 9999


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def iter_kmers(seq: str, k: int) -> Iterator[str]:
    seq = seq.upper()
    for start in range(max(0, len(seq) - k + 1)):
        kmer = seq[start : start + k]
        if "N" not in kmer:
            yield kmer


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name = ""
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and seq:
                    yield name, "".join(seq).upper().replace("-", "")
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name and seq:
        yield name, "".join(seq).upper().replace("-", "")


def allele_from_header(header: str) -> str | None:
    match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", header)
    return normalize_allele(match.group(1)) if match else None


def build_gene_unique_kmers(exon_dir: Path, genes: list[str], k: int) -> dict[str, set[str]]:
    gene_kmers: dict[str, set[str]] = {gene: set() for gene in genes}
    for gene in genes:
        fasta = exon_dir / f"{gene.replace('-', '_')}.fasta"
        for _header, seq in read_fasta(fasta):
            gene_kmers[gene].update(iter_kmers(seq, k))
            gene_kmers[gene].update(iter_kmers(revcomp(seq), k))
    owners: dict[str, list[str]] = defaultdict(list)
    for gene, kmers in gene_kmers.items():
        for kmer in kmers:
            owners[kmer].append(gene)
    return {gene: {kmer for kmer in kmers if owners[kmer] == [gene]} for gene, kmers in gene_kmers.items()}


def build_dpb1_family_unique_kmers(exon_dir: Path, k: int, gene_unique: Optional[set[str]]) -> tuple[dict[str, set[str]], dict[str, str]]:
    fasta = exon_dir / "HLA_DPB1.fasta"
    family_kmers: dict[str, set[str]] = defaultdict(set)
    for header, seq in read_fasta(fasta):
        family = allele_from_header(header)
        if not family:
            continue
        kmers = set(iter_kmers(seq, k)) | set(iter_kmers(revcomp(seq), k))
        family_kmers[family].update(kmers & gene_unique if gene_unique is not None else kmers)
    owners: dict[str, set[str]] = defaultdict(set)
    for family, kmers in family_kmers.items():
        for kmer in kmers:
            owners[kmer].add(family)
    unique_by_family = {
        family: {kmer for kmer in kmers if owners[kmer] == {family}}
        for family, kmers in family_kmers.items()
    }
    kmer_to_family = {
        kmer: family
        for family, kmers in unique_by_family.items()
        for kmer in kmers
    }
    return unique_by_family, kmer_to_family


FastqRecord = Tuple[str, str, str, str]


def open_fastq(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def iter_fastq(path: Path) -> Iterator[tuple[str, FastqRecord]]:
    with open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            qual = handle.readline().rstrip("\n")
            if not qual:
                return
            yield normalize_read_name(header), (header.rstrip("\n"), seq, plus, qual)


def read_names(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return {name for name, _record in iter_fastq(path)}


def record_families(seq1: str, seq2: str, k: int, kmer_to_family: dict[str, str]) -> set[str]:
    families = set()
    for seq in (seq1, seq2):
        for kmer in set(iter_kmers(seq, k)):
            family = kmer_to_family.get(kmer)
            if family:
                families.add(family)
    return families


def backup_path(path: Path, suffix: str) -> Path:
    return path.with_name(path.name.replace(".fq.gz", f".{suffix}.fq.gz"))


def ensure_backup(path: Path, suffix: str) -> Path:
    backup = backup_path(path, suffix)
    if path.exists() and not backup.exists():
        shutil.copy2(path, backup)
    return backup if backup.exists() else path


def append_records(target: Path, source: Path, records: list[FastqRecord]) -> None:
    tmp = target.with_name(f"{target.name}.tmp")
    if source.exists():
        shutil.copy2(source, tmp)
        mode = "at"
    else:
        mode = "wt"
    with gzip.open(tmp, mode) as handle:
        for record in records:
            handle.write("\n".join(record) + "\n")
    tmp.replace(target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq1", required=True, type=Path)
    parser.add_argument("--fq2", required=True, type=Path)
    parser.add_argument("--fq-dir", required=True, type=Path)
    parser.add_argument("--exon-dir", type=Path, default=DEFAULT_EXON_DIR)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--min-full-support", type=int, default=2)
    parser.add_argument("--min-missed-pairs", type=int, default=1)
    parser.add_argument("--max-retained-fraction", type=float, default=0.50)
    parser.add_argument("--min-unique-kmers", type=int, default=100)
    parser.add_argument("--max-allele-number", type=int, default=100)
    parser.add_argument("--max-add-pairs", type=int, default=5000)
    parser.add_argument("--use-gene-unique-filter", action="store_true")
    parser.add_argument("--backup-suffix", default="pre_dpb1_family_rescue")
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.use_gene_unique_filter:
        gene_unique = build_gene_unique_kmers(args.exon_dir, DEFAULT_GENES, args.k)["HLA-DPB1"]
    else:
        gene_unique = None
    unique_by_family, kmer_to_family = build_dpb1_family_unique_kmers(args.exon_dir, args.k, gene_unique)
    dpb1_r1 = args.fq_dir / "DPB1.R1.fq.gz"
    dpb1_r2 = args.fq_dir / "DPB1.R2.fq.gz"
    existing = read_names(dpb1_r1)

    full_support: dict[str, set[str]] = defaultdict(set)
    binned_support: dict[str, set[str]] = defaultdict(set)
    rescued_records: dict[str, tuple[FastqRecord, FastqRecord]] = {}
    read_families: dict[str, set[str]] = {}
    pairs_scanned = 0

    for (name1, rec1), (name2, rec2) in zip(iter_fastq(args.fq1), iter_fastq(args.fq2)):
        pairs_scanned += 1
        name = name1
        families = record_families(rec1[1], rec2[1], args.k, kmer_to_family)
        if not families:
            continue
        read_families[name] = families
        for family in families:
            full_support[family].add(name)
        if name not in existing:
            rescued_records[name] = (rec1, rec2)

    for name, rec in iter_fastq(dpb1_r1):
        families = record_families(rec[1], "", args.k, kmer_to_family)
        for family in families:
            binned_support[family].add(name)
    for name, rec in iter_fastq(dpb1_r2):
        families = record_families(rec[1], "", args.k, kmer_to_family)
        for family in families:
            binned_support[family].add(name)

    candidate_families = set()
    family_rows = []
    for family in sorted(full_support):
        full_count = len(full_support[family])
        binned_count = len(binned_support.get(family, set()))
        missed = len(full_support[family] - binned_support.get(family, set()))
        retained = binned_count / full_count if full_count else 0.0
        unique_count = len(unique_by_family.get(family, set()))
        eligible = (
            full_count >= args.min_full_support
            and missed >= args.min_missed_pairs
            and retained <= args.max_retained_fraction
            and unique_count >= args.min_unique_kmers
            and allele_number(family) <= args.max_allele_number
        )
        if eligible:
            candidate_families.add(family)
        family_rows.append({
            "family": family,
            "unique_kmers": unique_count,
            "full_support_pairs": full_count,
            "binned_support_pairs": binned_count,
            "missed_pairs": missed,
            "retained_fraction": f"{retained:.6f}",
            "candidate": int(eligible),
        })

    add_names = [
        name for name, families in read_families.items()
        if name not in existing and families & candidate_families
    ]
    add_names = add_names[: args.max_add_pairs]
    add_r1 = [rescued_records[name][0] for name in add_names if name in rescued_records]
    add_r2 = [rescued_records[name][1] for name in add_names if name in rescued_records]
    status = "dry_run" if args.dry_run else "written"
    if not add_names:
        status = "no_rescues"
    elif not args.dry_run:
        source_r1 = ensure_backup(dpb1_r1, args.backup_suffix)
        source_r2 = ensure_backup(dpb1_r2, args.backup_suffix)
        append_records(dpb1_r1, source_r1, add_r1)
        append_records(dpb1_r2, source_r2, add_r2)

    manifest = args.manifest or (args.fq_dir / "dpb1_family_rescue_manifest.tsv")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w") as handle:
        handle.write(f"# pairs_scanned\t{pairs_scanned}\n")
        handle.write(f"# existing_dpb1_pairs\t{len(existing)}\n")
        handle.write(f"# candidate_families\t{','.join(sorted(candidate_families))}\n")
        handle.write(f"# add_pairs\t{len(add_names)}\n")
        handle.write(f"# status\t{status}\n")
        fields = ["family", "unique_kmers", "full_support_pairs", "binned_support_pairs", "missed_pairs", "retained_fraction", "candidate"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(family_rows)
    print(f"wrote {manifest}")
    print(f"status={status} add_pairs={len(add_names)} candidate_families={','.join(sorted(candidate_families))}")


if __name__ == "__main__":
    main()