#!/usr/bin/env python3
"""Rescue gene-informative read pairs missed by strict read binning.

This truth-free helper scans the full sample FASTQs for k-mers that are unique
to one HLA gene across the bundled exon references, then appends unbinned read
pairs to that gene's per-gene FASTQs. It is intended to run after
assign_reads_to_genes.py and before per-gene alignment.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_EXON_DIR = SCRIPT_DIR / "resources" / "spechla" / "db" / "HLA" / "exon"
DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]
GENE_TO_SHORT = {
    "HLA-A": "A",
    "HLA-B": "B",
    "HLA-C": "C",
    "HLA-DRB1": "DRB1",
    "HLA-DPB1": "DPB1",
    "HLA-DQB1": "DQB1",
}


def normalize_gene(gene: str) -> str:
    gene = gene.strip()
    if not gene:
        raise ValueError("empty gene")
    if gene.startswith("HLA-"):
        return gene
    return f"HLA-{gene}"


def short_gene(gene: str) -> str:
    gene = normalize_gene(gene)
    return GENE_TO_SHORT.get(gene, gene.replace("HLA-", ""))


def normalize_read_name(header: str) -> str:
    name = header.strip().split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def iter_kmers(seq: str, k: int) -> Iterator[str]:
    seq = seq.upper()
    if len(seq) < k:
        return
    for start in range(len(seq) - k + 1):
        kmer = seq[start : start + k]
        if "N" not in kmer:
            yield kmer


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode)


def read_fasta(path: Path) -> Iterator[str]:
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    yield "".join(seq_parts).upper().replace("-", "")
                seq_parts = []
            else:
                seq_parts.append(line)
    if seq_parts:
        yield "".join(seq_parts).upper().replace("-", "")


def exon_fasta(exon_dir: Path, gene: str) -> Path:
    return exon_dir / f"{normalize_gene(gene).replace('-', '_')}.fasta"


def build_gene_unique_kmers(exon_dir: Path, genes: list[str], k: int) -> tuple[dict[str, set[str]], dict[str, str]]:
    gene_kmers: dict[str, set[str]] = {gene: set() for gene in genes}
    for gene in genes:
        fasta = exon_fasta(exon_dir, gene)
        if not fasta.exists():
            raise FileNotFoundError(fasta)
        for seq in read_fasta(fasta):
            gene_kmers[gene].update(iter_kmers(seq, k))
            gene_kmers[gene].update(iter_kmers(revcomp(seq), k))

    kmer_genes: dict[str, list[str]] = defaultdict(list)
    for gene, kmers in gene_kmers.items():
        for kmer in kmers:
            kmer_genes[kmer].append(gene)

    unique_by_gene: dict[str, set[str]] = {gene: set() for gene in genes}
    kmer_to_gene: dict[str, str] = {}
    for kmer, owners in kmer_genes.items():
        if len(owners) == 1:
            gene = owners[0]
            unique_by_gene[gene].add(kmer)
            kmer_to_gene[kmer] = gene
    return unique_by_gene, kmer_to_gene


FastqRecord = Tuple[str, str, str, str]


def iter_fastq(path: Path) -> Iterator[tuple[str, FastqRecord]]:
    with open_text(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                return
            record = (header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n"))
            yield normalize_read_name(header), record


def write_record(handle, record: FastqRecord) -> None:
    handle.write("\n".join(record) + "\n")


def read_fastq_names(path: Path) -> set[str]:
    names: set[str] = set()
    if not path.exists():
        return names
    for name, _record in iter_fastq(path):
        names.add(name)
    return names


def backup_path(path: Path, suffix: str) -> Path:
    name = path.name
    match = re.match(r"(.+)\.R([12])\.fq\.gz$", name)
    if match:
        return path.with_name(f"{match.group(1)}.R{match.group(2)}.{suffix}.fq.gz")
    return path.with_name(f"{name}.{suffix}")


def source_path(path: Path, suffix: str) -> Path:
    backup = backup_path(path, suffix)
    return backup if backup.exists() else path


def ensure_backup(path: Path, suffix: str) -> Path:
    backup = backup_path(path, suffix)
    if path.exists() and not backup.exists():
        shutil.copy2(path, backup)
    return backup if backup.exists() else path


def count_gene_hits_for_seq(seq: str, k: int, kmer_to_gene: dict[str, str]) -> Counter[str]:
    hits: Counter[str] = Counter()
    for kmer in set(iter_kmers(seq, k)):
        gene = kmer_to_gene.get(kmer)
        if gene:
            hits[gene] += 1
    return hits


def count_gene_hits(seq1: str, seq2: str, k: int, kmer_to_gene: dict[str, str]) -> tuple[Counter[str], Counter[str], Counter[str]]:
    hits1 = count_gene_hits_for_seq(seq1, k, kmer_to_gene)
    hits2 = count_gene_hits_for_seq(seq2, k, kmer_to_gene)
    seen_kmers = set(iter_kmers(seq1, k))
    seen_kmers.update(iter_kmers(seq2, k))
    hits: Counter[str] = Counter()
    for kmer in seen_kmers:
        gene = kmer_to_gene.get(kmer)
        if gene:
            hits[gene] += 1
    return hits, hits1, hits2


def choose_gene(
    hits: Counter[str],
    min_hits: int,
    min_margin: int,
    hits1: Counter[str] | None = None,
    hits2: Counter[str] | None = None,
    require_both_mates: bool = False,
    min_mate_hits: int = 1,
) -> tuple[str | None, str, int, int]:
    if not hits:
        return None, "no_gene_unique_kmer", 0, 0
    ranked = hits.most_common()
    top_gene, top_hits = ranked[0]
    second_hits = ranked[1][1] if len(ranked) > 1 else 0
    if top_hits < min_hits:
        return None, "below_min_hits", top_hits, second_hits
    if top_hits - second_hits < min_margin:
        return None, "ambiguous_gene_hits", top_hits, second_hits
    if require_both_mates:
        mate1_hits = hits1.get(top_gene, 0) if hits1 else 0
        mate2_hits = hits2.get(top_gene, 0) if hits2 else 0
        if mate1_hits < min_mate_hits or mate2_hits < min_mate_hits:
            return None, "missing_mate_support", top_hits, second_hits
    return top_gene, "rescued", top_hits, second_hits


def append_records(target: Path, source: Path, records: list[FastqRecord]) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(f"{target.name}.tmp")
    if source.exists():
        shutil.copy2(source, tmp)
        mode = "at"
    else:
        mode = "wt"
    with gzip.open(tmp, mode) as handle:
        for record in records:
            write_record(handle, record)
    tmp.replace(target)


def main() -> None:
    parser = argparse.ArgumentParser(description="Rescue missed HLA gene-binned read pairs from full FASTQs")
    parser.add_argument("--fq1", required=True, type=Path)
    parser.add_argument("--fq2", required=True, type=Path)
    parser.add_argument("--fq-dir", required=True, type=Path, help="Directory containing <gene>.R1/R2.fq.gz")
    parser.add_argument("--exon-dir", type=Path, default=DEFAULT_EXON_DIR)
    parser.add_argument("--gene", action="append", default=[],
                        help="Target gene to rescue. May be supplied more than once. Defaults to all supported genes.")
    parser.add_argument("--background-gene", action="append", default=[],
                        help="Gene set used to decide gene-unique k-mers. Defaults to all supported genes.")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--min-hits", type=int, default=1)
    parser.add_argument("--min-margin", type=int, default=1)
    parser.add_argument("--require-both-mates", action="store_true",
                        help="Require both R1 and R2 to contain target gene-unique k-mers")
    parser.add_argument("--min-mate-hits", type=int, default=1,
                        help="Minimum target gene-unique k-mers required on each mate when --require-both-mates is set")
    parser.add_argument("--max-rescue-fraction", type=float, default=0.25)
    parser.add_argument("--max-rescue-pairs", type=int, default=100000)
    parser.add_argument("--retention-gate", action="store_true",
                        help="Only rescue genes with abnormal loss of full-FASTQ gene-informative read pairs")
    parser.add_argument("--retention-min-full-pairs", type=int, default=50)
    parser.add_argument("--retention-max-retained-fraction", type=float, default=0.10)
    parser.add_argument("--retention-min-missing-fraction", type=float, default=0.30)
    parser.add_argument("--retention-min-rescue-pairs", type=int, default=50)
    parser.add_argument("--backup-suffix", default="pre_read_rescue")
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    background_genes = [normalize_gene(gene) for gene in (args.background_gene or DEFAULT_GENES)]
    target_genes = [normalize_gene(gene) for gene in (args.gene or background_genes)]
    unknown = [gene for gene in background_genes + target_genes if gene not in GENE_TO_SHORT]
    if unknown:
        raise SystemExit(f"unsupported genes: {','.join(unknown)}")
    missing_background = [gene for gene in target_genes if gene not in background_genes]
    if missing_background:
        raise SystemExit(
            "target genes must be included in the background gene set: " + ",".join(missing_background)
        )

    unique_by_gene, kmer_to_gene = build_gene_unique_kmers(args.exon_dir, background_genes, args.k)
    print("gene_unique_kmers", {gene: len(unique_by_gene[gene]) for gene in target_genes}, flush=True)

    existing_names: dict[str, set[str]] = {}
    original_existing_names: dict[str, set[str]] = {}
    r1_paths: dict[str, Path] = {}
    r2_paths: dict[str, Path] = {}
    source_r1_paths: dict[str, Path] = {}
    source_r2_paths: dict[str, Path] = {}
    original_pair_counts: dict[str, int] = {}
    for gene in target_genes:
        short = short_gene(gene)
        r1_path = args.fq_dir / f"{short}.R1.fq.gz"
        r2_path = args.fq_dir / f"{short}.R2.fq.gz"
        r1_paths[gene] = r1_path
        r2_paths[gene] = r2_path
        source_r1 = source_path(r1_path, args.backup_suffix)
        source_r2 = source_path(r2_path, args.backup_suffix)
        source_r1_paths[gene] = source_r1
        source_r2_paths[gene] = source_r2
        names = read_fastq_names(source_r1)
        original_existing_names[gene] = set(names)
        existing_names[gene] = set(names)
        original_pair_counts[gene] = len(names)

    rescued: dict[str, list[tuple[str, FastqRecord, FastqRecord, int, int]]] = {gene: [] for gene in target_genes}
    full_gene_pairs: Counter[str] = Counter()
    retained_gene_pairs: Counter[str] = Counter()
    counts: Counter[str] = Counter()

    for (name1, rec1), (name2, rec2) in zip(iter_fastq(args.fq1), iter_fastq(args.fq2)):
        counts["pairs_scanned"] += 1
        if name1 != name2:
            counts["name_mismatch_pairs"] += 1
        read_name = name1
        hits, hits1, hits2 = count_gene_hits(rec1[1], rec2[1], args.k, kmer_to_gene)
        gene, reason, top_hits, second_hits = choose_gene(
            hits,
            args.min_hits,
            args.min_margin,
            hits1=hits1,
            hits2=hits2,
            require_both_mates=args.require_both_mates,
            min_mate_hits=args.min_mate_hits,
        )
        counts[reason] += 1
        if gene is None:
            continue
        if gene not in rescued:
            counts["assigned_to_non_target_gene"] += 1
            continue
        full_gene_pairs[gene] += 1
        if read_name in original_existing_names[gene]:
            retained_gene_pairs[gene] += 1
            counts["already_in_target_gene"] += 1
            continue
        if read_name in existing_names[gene]:
            counts["duplicate_rescue_candidate"] += 1
            continue
        existing_names[gene].add(read_name)
        rescued[gene].append((read_name, rec1, rec2, top_hits, second_hits))

    rows = []
    for gene in target_genes:
        original_pairs = original_pair_counts[gene]
        rescue_pairs = len(rescued[gene])
        eligible_pairs = full_gene_pairs[gene]
        retained_pairs = retained_gene_pairs[gene]
        retained_fraction = retained_pairs / eligible_pairs if eligible_pairs else 0.0
        missing_fraction = rescue_pairs / eligible_pairs if eligible_pairs else 0.0
        max_by_fraction = int(original_pairs * args.max_rescue_fraction) if original_pairs else args.max_rescue_pairs
        allowed_pairs = min(args.max_rescue_pairs, max_by_fraction if max_by_fraction > 0 else args.max_rescue_pairs)
        status = "dry_run" if args.dry_run else "written"
        retention_reason = "not_checked"
        if rescue_pairs > allowed_pairs:
            status = "skipped_too_many_rescues"
        if args.retention_gate and rescue_pairs > 0 and status not in {"skipped_too_many_rescues"}:
            loss_anomaly = (
                retained_fraction <= args.retention_max_retained_fraction
                or missing_fraction >= args.retention_min_missing_fraction
            )
            enough_support = eligible_pairs >= args.retention_min_full_pairs
            enough_rescues = rescue_pairs >= args.retention_min_rescue_pairs
            if not enough_support:
                status = "skipped_low_full_gene_support"
                retention_reason = f"eligible_pairs<{args.retention_min_full_pairs}"
            elif not enough_rescues:
                status = "skipped_low_rescue_pairs"
                retention_reason = f"rescued_pairs<{args.retention_min_rescue_pairs}"
            elif not loss_anomaly:
                status = "skipped_retention_ok"
                retention_reason = (
                    f"retained_fraction>{args.retention_max_retained_fraction}"
                    f";missing_fraction<{args.retention_min_missing_fraction}"
                )
            else:
                retention_reason = "passed"
        if rescue_pairs == 0:
            status = "no_rescues" if not args.dry_run else "dry_run"
            retention_reason = "no_rescues"
        rows.append({
            "gene": gene,
            "short_gene": short_gene(gene),
            "k": args.k,
            "unique_kmers": len(unique_by_gene[gene]),
            "original_pairs": original_pairs,
            "full_gene_pairs": eligible_pairs,
            "retained_gene_pairs": retained_pairs,
            "retained_fraction": f"{retained_fraction:.6f}",
            "missing_fraction": f"{missing_fraction:.6f}",
            "rescued_pairs": rescue_pairs,
            "allowed_pairs": allowed_pairs,
            "retention_gate": int(args.retention_gate),
            "retention_reason": retention_reason,
            "status": status,
        })
        if args.dry_run or status != "written":
            continue
        r1_path = r1_paths[gene]
        r2_path = r2_paths[gene]
        source_r1 = ensure_backup(r1_path, args.backup_suffix)
        source_r2 = ensure_backup(r2_path, args.backup_suffix)
        append_records(r1_path, source_r1, [rec1 for _name, rec1, _rec2, _top, _second in rescued[gene]])
        append_records(r2_path, source_r2, [rec2 for _name, _rec1, rec2, _top, _second in rescued[gene]])

    manifest = args.manifest or (args.fq_dir / "read_bin_rescue_manifest.tsv")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w") as handle:
        fieldnames = [
            "gene", "short_gene", "k", "unique_kmers", "original_pairs",
            "full_gene_pairs", "retained_gene_pairs", "retained_fraction", "missing_fraction",
            "rescued_pairs", "allowed_pairs", "retention_gate", "retention_reason", "status",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        handle.write("# counters\n")
        for key, value in counts.most_common():
            handle.write(f"# {key}\t{value}\n")

    print(f"wrote {manifest}")
    for row in rows:
        print(
            row["gene"],
            "original=", row["original_pairs"],
            "rescued=", row["rescued_pairs"],
            "status=", row["status"],
            flush=True,
        )


if __name__ == "__main__":
    main()