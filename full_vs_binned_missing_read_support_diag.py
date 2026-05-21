#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


FULL_ROOT = Path("fqs")
SPECHLA_ROOT = Path("/data2/wangxuedong/polyploid-hla-realsets/spechla_out_abc_realsets_rescue_20260512")
EXON_DIR = Path("scripts/resources/spechla/db/HLA/exon")
INPUT_TSV = Path("diagnostics/remaining_missing_raw_read_support_common_minor_20260513.tsv")
DEFAULT_OUTPUT_TSV = Path("diagnostics/full_vs_binned_missing_read_support_20260513.tsv")
DEFAULT_SUMMARY = Path("diagnostics/full_vs_binned_missing_read_support_20260513.summary")

SET_DIR = {"set-a": "set_a", "set-b": "set B", "set-c": "set C"}


def normalize_read_name(header):
    name = header.strip().split()[0]
    if name.startswith("@"):
        name = name[1:]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def normalize_allele(allele):
    allele = (allele or "").strip().replace("HLA-", "").replace("G", "").rstrip("P")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    parts = rest.split(":")
    if parts and parts[-1].isalpha():
        parts[-1] = parts[-1][:-1]
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def iter_kmers(seq, k):
    seq = seq.upper()
    if len(seq) < k:
        return
    for start in range(len(seq) - k + 1):
        kmer = seq[start : start + k]
        if "N" not in kmer:
            yield kmer


def read_fasta(path):
    name = None
    seq = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and seq:
                    yield name, "".join(seq).upper()
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name and seq:
        yield name, "".join(seq).upper()


def allele_from_header(header):
    match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", header)
    return normalize_allele(match.group(1)) if match else None


def build_unique_target_kmers(exon_dir, gene, k, target_alleles):
    fasta = exon_dir / f"{gene.replace('-', '_')}.fasta"
    if not fasta.exists():
        raise FileNotFoundError(fasta)

    target = set(target_alleles)
    target_kmers = {allele: set() for allele in target}
    for header, seq in read_fasta(fasta):
        family = allele_from_header(header)
        if not family:
            continue
        if family not in target:
            continue
        target_kmers[family].update(iter_kmers(seq, k))
        target_kmers[family].update(iter_kmers(revcomp(seq), k))

    target_union = set()
    for kmers in target_kmers.values():
        target_union.update(kmers)
    if not target_union:
        return target_kmers

    owners = defaultdict(set)
    for header, seq in read_fasta(fasta):
        family = allele_from_header(header)
        if not family:
            continue
        seen = (set(iter_kmers(seq, k)) | set(iter_kmers(revcomp(seq), k))) & target_union
        for kmer in seen:
            owners[kmer].add(family)

    unique = {allele: set() for allele in target}
    for allele, kmers in target_kmers.items():
        unique[allele] = {kmer for kmer in kmers if owners.get(kmer) == {allele}}
    return unique


def iter_fastq_records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            if seq:
                yield normalize_read_name(header), seq


def full_fastqs(full_root, set_id, sample):
    folder = full_root / SET_DIR[set_id]
    return [folder / f"{sample}_R1_001.fastq.gz", folder / f"{sample}_R2_001.fastq.gz"]


def gene_fastqs(spechla_root, sample, gene):
    short = gene.split("-", 1)[1]
    return [spechla_root / sample / f"{short}.R1.fq.gz", spechla_root / sample / f"{short}.R2.fq.gz"]


def scan(paths, allele_kmers, k):
    support = {allele: set() for allele in allele_kmers}
    kmer_to_alleles = defaultdict(list)
    for allele, kmers in allele_kmers.items():
        for kmer in kmers:
            kmer_to_alleles[kmer].append(allele)

    total = 0
    missing = []
    for path in paths:
        if not path.exists():
            missing.append(str(path))
            continue
        for read_name, seq in iter_fastq_records(path):
            total += 1
            hit_alleles = set()
            for kmer in iter_kmers(seq, k):
                for allele in kmer_to_alleles.get(kmer, ()):
                    hit_alleles.add(allele)
            for allele in hit_alleles:
                support[allele].add(read_name)
    return total, support, missing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", action="append", default=[])
    parser.add_argument("--side", action="append", default=[])
    parser.add_argument("--k", action="append", type=int, default=[])
    parser.add_argument("--input", type=Path, default=INPUT_TSV)
    parser.add_argument("--full-root", type=Path, default=FULL_ROOT)
    parser.add_argument("--spechla-root", type=Path, default=SPECHLA_ROOT)
    parser.add_argument("--exon-dir", type=Path, default=EXON_DIR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT_TSV)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    args = parser.parse_args()

    genes = set(args.gene)
    sides = set(args.side)
    ks = tuple(args.k or [31, 51])
    rows = []
    with args.input.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if genes and row["gene"] not in genes:
                continue
            if sides and row["side"] not in sides:
                continue
            rows.append(row)
    if not rows:
        raise SystemExit("no rows selected")

    targets_by_gene = defaultdict(set)
    for row in rows:
        targets_by_gene[row["gene"]].add(row["missing_allele"])

    kmers = {}
    for gene, alleles in sorted(targets_by_gene.items()):
        for k in ks:
            kmers[(gene, k)] = build_unique_target_kmers(args.exon_dir, gene, k, alleles)
            print("built", gene, k, {allele: len(vals) for allele, vals in kmers[(gene, k)].items()})

    full_cache = {}
    gene_cache = {}
    all_missing_paths = []
    for row in rows:
        sample, set_id, gene = row["sample"], row["set"], row["gene"]
        for k in ks:
            key = (sample, gene, k)
            if key not in full_cache:
                full_cache[key] = scan(full_fastqs(args.full_root, set_id, sample), kmers[(gene, k)], k)
                all_missing_paths.extend(full_cache[key][2])
                print("scanned full", key, "reads", full_cache[key][0])
            row_spechla_root = Path(row.get("spechla_root") or args.spechla_root)
            gene_key = (str(row_spechla_root), sample, gene, k)
            if gene_key not in gene_cache:
                gene_cache[gene_key] = scan(gene_fastqs(row_spechla_root, sample, gene), kmers[(gene, k)], k)
                all_missing_paths.extend(gene_cache[gene_key][2])
                print("scanned gene", gene_key, "reads", gene_cache[gene_key][0])

    fields = [
        "set",
        "sample",
        "gene",
        "side",
        "missing_allele",
        "mean_mask",
        "chi_r_fit",
        "expected_side_single_fraction",
        "em_rank",
        "em_frac",
        "full_total_reads_k31",
        "binned_total_reads_k31",
        "unique_kmers_k31",
        "full_support_reads_k31",
        "binned_support_reads_k31",
        "missed_by_binning_reads_k31",
        "retained_fraction_k31",
        "full_total_reads_k51",
        "binned_total_reads_k51",
        "unique_kmers_k51",
        "full_support_reads_k51",
        "binned_support_reads_k51",
        "missed_by_binning_reads_k51",
        "retained_fraction_k51",
        "binning_diagnosis",
    ]
    output_rows = []
    for row in rows:
        sample, gene, allele = row["sample"], row["gene"], row["missing_allele"]
        out = {field: row.get(field, "") for field in fields[:10]}
        for k in (31, 51):
            out[f"full_total_reads_k{k}"] = 0
            out[f"binned_total_reads_k{k}"] = 0
            out[f"unique_kmers_k{k}"] = 0
            out[f"full_support_reads_k{k}"] = 0
            out[f"binned_support_reads_k{k}"] = 0
            out[f"missed_by_binning_reads_k{k}"] = 0
            out[f"retained_fraction_k{k}"] = "NA"
        for k in ks:
            full_total, full_support, _ = full_cache[(sample, gene, k)]
            row_spechla_root = Path(row.get("spechla_root") or args.spechla_root)
            binned_total, binned_support, _ = gene_cache[(str(row_spechla_root), sample, gene, k)]
            full_names = full_support.get(allele, set())
            binned_names = binned_support.get(allele, set())
            missed = full_names - binned_names
            retained = (len(full_names & binned_names) / len(full_names)) if full_names else None
            out[f"full_total_reads_k{k}"] = full_total
            out[f"binned_total_reads_k{k}"] = binned_total
            out[f"unique_kmers_k{k}"] = len(kmers[(gene, k)].get(allele, set()))
            out[f"full_support_reads_k{k}"] = len(full_names)
            out[f"binned_support_reads_k{k}"] = len(binned_names)
            out[f"missed_by_binning_reads_k{k}"] = len(missed)
            out[f"retained_fraction_k{k}"] = "NA" if retained is None else f"{retained:.4f}"

        full31 = int(out["full_support_reads_k31"])
        full51 = int(out["full_support_reads_k51"])
        binned31 = int(out["binned_support_reads_k31"])
        binned51 = int(out["binned_support_reads_k51"])
        missed31 = int(out["missed_by_binning_reads_k31"])
        missed51 = int(out["missed_by_binning_reads_k51"])
        if full31 == 0 and full51 == 0:
            diagnosis = "not_in_full_fastq_by_unique_kmers"
        elif (missed31 > 0 or missed51 > 0) and (binned31 == 0 and binned51 == 0):
            diagnosis = "present_in_full_but_lost_by_binning"
        elif missed31 > 0 or missed51 > 0:
            diagnosis = "partly_lost_by_binning"
        else:
            diagnosis = "retained_by_binning_if_detectable"
        out["binning_diagnosis"] = diagnosis
        output_rows.append(out)

    args.out.parent.mkdir(exist_ok=True)
    with args.out.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)

    by_diag = Counter(row["binning_diagnosis"] for row in output_rows)
    by_side_diag = defaultdict(Counter)
    by_gene = Counter()
    by_gene_diag = defaultdict(Counter)
    for row in output_rows:
        by_side_diag[row["side"]][row["binning_diagnosis"]] += 1
        by_gene[row["gene"]] += 1
        by_gene_diag[row["gene"]][row["binning_diagnosis"]] += 1

    with args.summary.open("w") as handle:
        handle.write(f"records\t{len(output_rows)}\n")
        handle.write(f"missing_paths\t{len(set(all_missing_paths))}\n")
        for path in sorted(set(all_missing_paths))[:20]:
            handle.write(f"missing_path\t{path}\n")
        handle.write("by_binning_diagnosis\n")
        for diagnosis, count in by_diag.most_common():
            handle.write(f"{diagnosis}\t{count}\n")
        handle.write("by_side_binning_diagnosis\n")
        for side, counter in sorted(by_side_diag.items()):
            handle.write(f"{side}\t{dict(counter)}\n")
        handle.write("by_gene\n")
        for gene, count in by_gene.most_common():
            handle.write(f"{gene}\t{count}\t{dict(by_gene_diag[gene])}\n")
        handle.write("rows_with_full_support_lost_or_partial\n")
        for row in output_rows:
            if row["binning_diagnosis"] in ("present_in_full_but_lost_by_binning", "partly_lost_by_binning"):
                handle.write(
                    "\t".join(
                        str(row[field])
                        for field in [
                            "set",
                            "sample",
                            "gene",
                            "side",
                            "missing_allele",
                            "full_support_reads_k31",
                            "binned_support_reads_k31",
                            "missed_by_binning_reads_k31",
                            "full_support_reads_k51",
                            "binned_support_reads_k51",
                            "missed_by_binning_reads_k51",
                            "binning_diagnosis",
                        ]
                    )
                    + "\n"
                )

    print("wrote", args.out)
    print("wrote", args.summary)
    print(args.summary.read_text())


if __name__ == "__main__":
    main()