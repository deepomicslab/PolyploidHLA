#!/usr/bin/env python3
"""Diagnose raw FASTQ support for remaining GenDx missing truth alleles.

This is a diagnostic only. It can optionally apply a proposal manifest in memory
to ask what errors remain after a guarded rescue layer, but it never rewrites
pipeline outputs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_EXON_DIR = SCRIPT_DIR / "resources" / "spechla" / "db" / "HLA" / "exon"
GENE_SHORT = {
    "HLA-A": "A",
    "HLA-B": "B",
    "HLA-C": "C",
    "HLA-DRB1": "DRB1",
    "HLA-DQB1": "DQB1",
    "HLA-DPB1": "DPB1",
}


def read_tsv(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_alleles(text: str) -> list[str]:
    return [item for item in (text or "").split(",") if item]


def normalize_allele(allele: str) -> str:
    allele = (allele or "").strip().replace("HLA-", "").replace("G", "").rstrip("P")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    parts = rest.split(":")
    if parts and parts[-1].isalpha():
        parts[-1] = parts[-1][:-1]
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def normalize_read_name(header: str) -> str:
    name = header.strip().split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def iter_kmers(seq: str, k: int) -> Iterator[str]:
    seq = seq.upper()
    if len(seq) < k:
        return
    for start in range(len(seq) - k + 1):
        kmer = seq[start:start + k]
        if "N" not in kmer:
            yield kmer


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name = ""
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and seq_parts:
                    yield name, "".join(seq_parts).upper().replace("-", "")
                name = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name and seq_parts:
        yield name, "".join(seq_parts).upper().replace("-", "")


def allele_from_header(header: str) -> str | None:
    match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", header)
    return normalize_allele(match.group(1)) if match else None


def build_unique_target_kmers(exon_dir: Path, gene: str, alleles: set[str], k: int):
    fasta = exon_dir / f"{gene.replace('-', '_')}.fasta"
    target_kmers = {allele: set() for allele in alleles}
    if not fasta.exists():
        return target_kmers
    for header, seq in read_fasta(fasta):
        allele = allele_from_header(header)
        if allele not in target_kmers:
            continue
        target_kmers[allele].update(iter_kmers(seq, k))
        target_kmers[allele].update(iter_kmers(revcomp(seq), k))
    target_union = set().union(*target_kmers.values()) if target_kmers else set()
    owners = defaultdict(set)
    if target_union:
        for header, seq in read_fasta(fasta):
            allele = allele_from_header(header)
            if not allele:
                continue
            seen = (set(iter_kmers(seq, k)) | set(iter_kmers(revcomp(seq), k))) & target_union
            for kmer in seen:
                owners[kmer].add(allele)
    return {
        allele: {kmer for kmer in kmers if owners.get(kmer) == {allele}}
        for allele, kmers in target_kmers.items()
    }


def iter_fastq_pairs(fq1: Path, fq2: Path):
    with gzip.open(fq1, "rt") as h1, gzip.open(fq2, "rt") as h2:
        while True:
            r1_header = h1.readline()
            r2_header = h2.readline()
            if not r1_header or not r2_header:
                return
            r1_seq = h1.readline().strip().upper()
            r2_seq = h2.readline().strip().upper()
            h1.readline(); h1.readline()
            h2.readline(); h2.readline()
            yield normalize_read_name(r1_header), r1_seq, r2_seq


def scan_pairs(fq1: Path, fq2: Path, allele_kmers: dict[str, set[str]], k: int):
    support = {allele: set() for allele in allele_kmers}
    kmer_to_alleles = defaultdict(list)
    for allele, kmers in allele_kmers.items():
        for kmer in kmers:
            kmer_to_alleles[kmer].append(allele)
    if not fq1.exists() or not fq2.exists():
        return 0, support, [str(fq1), str(fq2)]
    total = 0
    for read_name, seq1, seq2 in iter_fastq_pairs(fq1, fq2):
        total += 1
        hits = set()
        for seq in (seq1, seq2):
            for kmer in iter_kmers(seq, k):
                hits.update(kmer_to_alleles.get(kmer, ()))
        for allele in hits:
            support[allele].add(read_name)
    return total, support, []


def tf_counts_for(spechla_root: Path, sample: str, gene: str):
    path = spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
    rows = read_tsv(path)
    out = {}
    for index, row in enumerate(rows, 1):
        allele = normalize_allele(row.get("allele_2field", ""))
        if not allele:
            continue
        try:
            frac = float(row.get("fraction", "0") or 0)
        except ValueError:
            frac = 0.0
        out[allele] = (index, frac)
    return out


def consume_missing(truth_vals: list[str], pred_vals: list[str]):
    remaining = Counter(pred_vals)
    missing = []
    for allele in truth_vals:
        if remaining[allele] > 0:
            remaining[allele] -= 1
        else:
            missing.append(allele)
    return missing


def classify(full_count: int, binned_count: int, unique_count: int, full_min: int, retained: str) -> str:
    if unique_count == 0:
        return "no_unique_kmers"
    if full_count == 0:
        return "no_raw_support_by_unique_kmers"
    if full_count < full_min:
        return "low_raw_support"
    retained_value = 0.0 if retained == "NA" else float(retained)
    if binned_count == 0 or retained_value < 0.25:
        return "binning_loss"
    return "in_bin_not_selected"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quartet-summary", required=True, type=Path)
    parser.add_argument("--combined-manifest", type=Path, default=None)
    parser.add_argument("--fq-root", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--exon-dir", type=Path, default=DEFAULT_EXON_DIR)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--k", action="append", type=int, default=[])
    parser.add_argument("--min-full-support", type=int, default=5)
    args = parser.parse_args()

    ks = args.k or [31, 51]
    proposals = {}
    if args.combined_manifest:
        for row in read_tsv(args.combined_manifest):
            proposals[(row["sample"], row["gene"])] = split_alleles(row["proposed_quartet"])

    missing_rows = []
    targets_by_gene = defaultdict(set)
    for row in read_tsv(args.quartet_summary):
        sample, gene = row["sample"], row["gene"]
        pred = proposals.get((sample, gene), split_alleles(row["pred_R"]) + split_alleles(row["pred_D"]))
        truth_r = split_alleles(row["truth_R"])
        truth_d = split_alleles(row["truth_D"])
        for side, allele in [("R", a) for a in consume_missing(truth_r, pred[:2])] + [("D", a) for a in consume_missing(truth_d, pred[2:4])]:
            norm = normalize_allele(allele)
            missing_rows.append({
                "set": row["set"],
                "sample": sample,
                "gene": gene,
                "side": side,
                "missing_allele": norm,
                "post_pred_R": ",".join(pred[:2]),
                "post_pred_D": ",".join(pred[2:4]),
                "original_score": row["score2"],
            })
            targets_by_gene[gene].add(norm)

    kmer_cache = {}
    for gene, alleles in targets_by_gene.items():
        for k in ks:
            kmer_cache[(gene, k)] = build_unique_target_kmers(args.exon_dir, gene, alleles, k)

    scan_cache = {}
    out_rows = []
    for row in missing_rows:
        sample, gene, allele = row["sample"], row["gene"], row["missing_allele"]
        fq1 = args.fq_root / f"{sample}_R1_001.fastq.gz"
        fq2 = args.fq_root / f"{sample}_R2_001.fastq.gz"
        short = GENE_SHORT[gene]
        b1 = args.spechla_root / sample / f"{short}.R1.fq.gz"
        b2 = args.spechla_root / sample / f"{short}.R2.fq.gz"
        tf = tf_counts_for(args.spechla_root, sample, gene)
        out = dict(row)
        for k in ks:
            allele_kmers = kmer_cache[(gene, k)]
            full_key = (sample, gene, k, "full")
            bin_key = (sample, gene, k, "binned")
            if full_key not in scan_cache:
                scan_cache[full_key] = scan_pairs(fq1, fq2, allele_kmers, k)
            if bin_key not in scan_cache:
                scan_cache[bin_key] = scan_pairs(b1, b2, allele_kmers, k)
            full_total, full_support, full_missing_paths = scan_cache[full_key]
            binned_total, binned_support, binned_missing_paths = scan_cache[bin_key]
            full_count = len(full_support.get(allele, set()))
            binned_count = len(binned_support.get(allele, set()))
            retained = "NA" if full_count == 0 else f"{binned_count / full_count:.6f}"
            out[f"unique_kmers_k{k}"] = len(allele_kmers.get(allele, set()))
            out[f"full_total_pairs_k{k}"] = full_total
            out[f"binned_total_pairs_k{k}"] = binned_total
            out[f"full_support_pairs_k{k}"] = full_count
            out[f"binned_support_pairs_k{k}"] = binned_count
            out[f"missed_by_binning_pairs_k{k}"] = max(0, full_count - binned_count)
            out[f"retained_fraction_k{k}"] = retained
            out[f"missing_paths_k{k}"] = ";".join(full_missing_paths + binned_missing_paths)
        rank, frac = tf.get(allele, ("NA", 0.0))
        out["em_rank"] = rank
        out["em_frac"] = f"{frac:.8f}" if isinstance(frac, float) else frac
        primary_k = ks[0]
        out["raw_support_class"] = classify(
            int(out[f"full_support_pairs_k{primary_k}"]),
            int(out[f"binned_support_pairs_k{primary_k}"]),
            int(out[f"unique_kmers_k{primary_k}"]),
            args.min_full_support,
            out[f"retained_fraction_k{primary_k}"],
        )
        out_rows.append(out)

    fields = ["set", "sample", "gene", "side", "missing_allele", "original_score", "post_pred_R", "post_pred_D"]
    for k in ks:
        fields.extend([
            f"unique_kmers_k{k}", f"full_total_pairs_k{k}", f"binned_total_pairs_k{k}",
            f"full_support_pairs_k{k}", f"binned_support_pairs_k{k}",
            f"missed_by_binning_pairs_k{k}", f"retained_fraction_k{k}", f"missing_paths_k{k}",
        ])
    fields.extend(["em_rank", "em_frac", "raw_support_class"])
    write_tsv(args.out, fields, out_rows)

    summary = defaultdict(lambda: {"missing": 0, "full_support": 0, "binned_support": 0})
    primary_k = ks[0]
    for row in out_rows:
        key = (row["gene"], row["raw_support_class"])
        summary[key]["missing"] += 1
        summary[key]["full_support"] += int(row[f"full_support_pairs_k{primary_k}"])
        summary[key]["binned_support"] += int(row[f"binned_support_pairs_k{primary_k}"])
    summary_rows = []
    for (gene, klass), stats in sorted(summary.items()):
        summary_rows.append({
            "gene": gene,
            "raw_support_class": klass,
            "missing_allele_copies": stats["missing"],
            "full_support_pairs_k" + str(primary_k): stats["full_support"],
            "binned_support_pairs_k" + str(primary_k): stats["binned_support"],
        })
    write_tsv(args.summary, ["gene", "raw_support_class", "missing_allele_copies", "full_support_pairs_k" + str(primary_k), "binned_support_pairs_k" + str(primary_k)], summary_rows)
    print(f"missing_copies\t{len(out_rows)}")
    print(f"wrote\t{args.out}")
    print(f"wrote\t{args.summary}")


if __name__ == "__main__":
    main()