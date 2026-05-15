#!/usr/bin/env python3
"""Truth-free DRB3/4/5 typing linked to the current DRB1 haplotypes.

DRB345 is a reporting convenience, not a single locus. This script keeps it as
an independent diagnostic/add-on path: extract read pairs with competitive DB
alignments to DRB3/DRB4/DRB5, EM-remap them to a combined DRB3/4/5 allele
reference, then choose one DRB345 allele per R1/R2/D1/D2 haplotype constrained
by the current DRB1-linked locus expectation.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Optional

from aggregate_calls import allele_2field, allele_g_group, clean_allele, load_g_group
from iterative_remap_em import (
    build_aug_ref,
    bwa_mem_all,
    load_imgt,
    parse_sam_to_reads,
    run_em,
    safe,
    two_field,
)


TARGET_PREFIXES = ("DRB3*", "DRB4*", "DRB5*")
DISABLE_MIN_AS = -100000000
DRB345_LOCI = ("DRB3", "DRB4", "DRB5")
HAP_SLOTS = (
    ("R1", "R"),
    ("R2", "R"),
    ("D1", "D"),
    ("D2", "D"),
)
DEFAULT_GENES = ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5")
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SPECHLA = Path(os.environ.get("SPECHLA", SCRIPT_DIR / "resources" / "spechla"))
DEFAULT_IMGT = DEFAULT_SPECHLA / "db" / "ref" / "hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
DEFAULT_G_GROUP = DEFAULT_SPECHLA / "db" / "HLA" / "hla_nom_g.txt"


def alignment_score(aln) -> int:
    try:
        return int(aln.get_tag("AS"))
    except KeyError:
        return int(aln.query_alignment_length or 0)


def passes_relative_score_gate(target_score: int, best_score: int, min_as_frac: float) -> bool:
    if min_as_frac <= 0:
        return True
    if best_score > 0:
        return target_score >= best_score * min_as_frac
    if best_score < 0:
        return target_score >= best_score / min_as_frac
    return target_score >= 0


def collect_drb345_read_names(db_bam: Path, min_as_frac: float, min_as: int) -> set[str]:
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit("type_drb345.py requires pysam to read the DB BAM") from exc

    best_by_read: dict[str, int] = defaultdict(lambda: -10**9)
    target_best_by_read: dict[str, int] = defaultdict(lambda: -10**9)
    with pysam.AlignmentFile(str(db_bam), "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_id < 0:
                continue
            qname = aln.query_name
            score = alignment_score(aln)
            if score > best_by_read[qname]:
                best_by_read[qname] = score
            ref = bam.get_reference_name(aln.reference_id) or ""
            if ref.startswith(TARGET_PREFIXES) and score > target_best_by_read[qname]:
                target_best_by_read[qname] = score

    keep = set()
    for qname, target_score in target_best_by_read.items():
        best = best_by_read.get(qname, target_score)
        if target_score >= min_as and passes_relative_score_gate(target_score, best, min_as_frac):
            keep.add(qname)
    return keep


def fastq_name(header: str) -> str:
    name = header[1:].split()[0]
    return name[:-2] if name.endswith(("/1", "/2")) else name


def write_selected_fastq(src: Path, dst: Path, wanted: set[str]) -> int:
    opener = gzip.open if src.suffix == ".gz" else open
    n = 0
    with opener(src, "rt") as inp, gzip.open(dst, "wt") as out:
        while True:
            h = inp.readline()
            if not h:
                break
            s = inp.readline()
            p = inp.readline()
            q = inp.readline()
            if not q:
                break
            if fastq_name(h.rstrip("\n")) in wanted:
                out.write(h)
                out.write(s)
                out.write(p)
                out.write(q)
                n += 1
    return n


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def iter_kmers(seq: str, k: int):
    seq = seq.upper()
    for idx in range(0, len(seq) - k + 1):
        kmer = seq[idx:idx + k]
        if "N" not in kmer:
            yield kmer


def locus_unique_kmers(imgt: Path, k: int) -> dict[str, set[str]]:
    db = load_imgt(str(imgt))
    owners: dict[str, Optional[str]] = {}
    for name, seq in db.items():
        locus = locus_of_allele(name)
        if locus not in DRB345_LOCI:
            continue
        for kmer in iter_kmers(seq, k):
            prev = owners.get(kmer)
            owners[kmer] = locus if prev is None else (prev if prev == locus else "multi")
    unique = {locus: set() for locus in DRB345_LOCI}
    for kmer, locus in owners.items():
        if locus in unique:
            unique[locus].add(kmer)
            unique[locus].add(revcomp(kmer))
    return unique


def locus_unique_pair_support(imgt: Path, fq1: Path, fq2: Path, k: int) -> dict[str, object]:
    unique = locus_unique_kmers(imgt, k)
    kmer_to_locus = {}
    for locus, kmers in unique.items():
        for kmer in kmers:
            kmer_to_locus[kmer] = locus
    counts = {locus: 0 for locus in DRB345_LOCI}
    total = 0
    multi = 0
    opener1 = gzip.open if fq1.suffix == ".gz" else open
    opener2 = gzip.open if fq2.suffix == ".gz" else open
    with opener1(fq1, "rt") as inp1, opener2(fq2, "rt") as inp2:
        while True:
            h1 = inp1.readline()
            h2 = inp2.readline()
            if not h1 or not h2:
                break
            s1 = inp1.readline().strip()
            s2 = inp2.readline().strip()
            inp1.readline(); inp1.readline()
            inp2.readline(); inp2.readline()
            total += 1
            hit_loci = set()
            for seq in (s1, s2):
                for kmer in iter_kmers(seq, k):
                    locus = kmer_to_locus.get(kmer)
                    if locus:
                        hit_loci.add(locus)
            for locus in hit_loci:
                counts[locus] += 1
            if len(hit_loci) > 1:
                multi += 1
    return {
        "total_pairs": total,
        "counts": counts,
        "fractions": {locus: (counts[locus] / total if total else 0.0) for locus in DRB345_LOCI},
        "unique_kmers": {locus: len(unique[locus]) for locus in DRB345_LOCI},
        "multi_pairs": multi,
    }


def drb1_linked_locus(allele: str) -> str:
    allele = clean_allele(allele)
    if "*" not in allele:
        return "absent"
    gene, rest = allele.split("*", 1)
    if gene != "DRB1":
        return "absent"
    family = rest.split(":", 1)[0]
    if family in {"03", "11", "12", "13", "14"}:
        return "DRB3"
    if family in {"04", "07", "09"}:
        return "DRB4"
    if family in {"15", "16"}:
        return "DRB5"
    return "absent"


def read_drb1_haplotypes(final_calls: Path) -> dict[str, str]:
    rows = read_tsv(final_calls)
    for row in rows:
        if row.get("gene") == "HLA-DRB1":
            return {slot: row.get(f"{slot}_full", "NA") for slot, _side in HAP_SLOTS}
    return {slot: "NA" for slot, _side in HAP_SLOTS}


def read_drb1_row(final_calls: Path) -> dict[str, str]:
    for row in read_tsv(final_calls):
        if row.get("gene") == "HLA-DRB1":
            return row
    return {}


def drb1_is_untrusted(row: dict[str, str], mask_threshold: float) -> bool:
    warning = row.get("warning", "")
    if "high_mask" in warning:
        return True
    try:
        if float(row.get("mean_mask_fraction", "0") or 0) >= mask_threshold:
            return True
    except ValueError:
        pass
    return False


def read_chi_r(fq_dir: Path, sample: str) -> Optional[float]:
    pooled = fq_dir / f"{sample}.chi_pooled.txt"
    if pooled.exists():
        for line in pooled.read_text().splitlines():
            if line.startswith("GLOBAL") and "chi_R=" in line:
                for item in line.split():
                    if item.startswith("chi_R="):
                        try:
                            return float(item.split("=", 1)[1])
                        except ValueError:
                            pass
    chimerism = fq_dir / f"{sample}.chimerism.txt"
    if chimerism.exists():
        for line in chimerism.read_text().splitlines():
            if "chi_R=" not in line:
                continue
            for item in line.split():
                if item.startswith("chi_R="):
                    try:
                        return float(item.split("=", 1)[1])
                    except ValueError:
                        pass
    return None


def build_drb345_candidates(imgt: Path, subs_per_2field: int) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    db = load_imgt(str(imgt))
    by_tf: dict[str, list[str]] = defaultdict(list)
    for name in db:
        if name.startswith(TARGET_PREFIXES):
            by_tf[two_field(name)].append(name)
    chosen = []
    for names in by_tf.values():
        names.sort(key=lambda n: -len(db[n]))
        chosen.extend(names[:subs_per_2field])
    contigs = {safe(name): db[name] for name in chosen}
    safe2name = {safe(name): name for name in chosen}
    tf_to_safe: dict[str, str] = {}
    for contig, name in safe2name.items():
        tf = two_field(name)
        if tf not in tf_to_safe or len(contigs[contig]) > len(contigs[tf_to_safe[tf]]):
            tf_to_safe[tf] = contig
    return contigs, safe2name, tf_to_safe


def locus_of_allele(allele: str) -> str:
    return allele.split("*", 1)[0] if "*" in allele else "absent"


def candidates_by_locus(tf_counts: dict[str, float], per_locus: int, allowed_loci: Optional[set[str]] = None) -> dict[str, list[str]]:
    grouped: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for allele, count in tf_counts.items():
        locus = locus_of_allele(allele)
        if allowed_loci is not None and locus not in allowed_loci:
            continue
        grouped[locus].append((allele, count))
    out = {}
    for locus in DRB345_LOCI:
        values = sorted(grouped.get(locus, []), key=lambda item: (-item[1], item[0]))[:per_locus]
        out[locus] = [allele for allele, _count in values]
    return out


def score_quartet(quartet, slots, observed: dict[str, float]) -> float:
    expected = defaultdict(float)
    for allele, (_slot, _side, _locus, weight, _candidates) in zip(quartet, slots):
        if allele != "NA":
            expected[allele] += weight
    scale = sum(expected.values())
    if scale > 0:
        expected = defaultdict(float, {allele: value / scale for allele, value in expected.items()})
    names = set(observed) | set(expected)
    return sum(abs(observed.get(name, 0.0) - expected.get(name, 0.0)) for name in names)


def fit_linked_calls(tf_counts: dict[str, float], drb1: dict[str, str], chi_r: float, per_locus: int):
    total = sum(tf_counts.values()) or 1.0
    observed = {allele: count / total for allele, count in tf_counts.items()}
    locus_candidates = candidates_by_locus(tf_counts, per_locus)
    slots = []
    for slot, side in HAP_SLOTS:
        locus = drb1_linked_locus(drb1.get(slot, "NA"))
        if locus == "absent":
            candidates = ["NA"]
        else:
            candidates = locus_candidates.get(locus) or ["NA"]
        weight = chi_r / 2.0 if side == "R" else (1.0 - chi_r) / 2.0
        slots.append((slot, side, locus, weight, candidates))

    best = None
    for quartet in product(*(slot[-1] for slot in slots)):
        score = score_quartet(quartet, slots, observed)
        if best is None or score < best[0]:
            best = (score, quartet)
    score, quartet = best if best is not None else (math.inf, ("NA", "NA", "NA", "NA"))
    return slots, quartet, score


def fit_evidence_calls(tf_counts: dict[str, float], chi_r: float, per_locus: int, allowed_loci: set[str]):
    filtered = {allele: count for allele, count in tf_counts.items() if locus_of_allele(allele) in allowed_loci}
    total = sum(filtered.values()) or 1.0
    observed = {allele: count / total for allele, count in filtered.items()}
    locus_candidates = candidates_by_locus(filtered, per_locus, allowed_loci)
    candidate_pool = sorted({allele for values in locus_candidates.values() for allele in values}, key=lambda a: (-filtered.get(a, 0.0), a))
    if not candidate_pool:
        candidate_pool = ["NA"]
    slots = []
    for slot, side in HAP_SLOTS:
        weight = chi_r / 2.0 if side == "R" else (1.0 - chi_r) / 2.0
        slots.append((slot, side, "evidence", weight, candidate_pool))

    best = None
    for quartet in product(*(slot[-1] for slot in slots)):
        score = score_quartet(quartet, slots, observed)
        if best is None or score < best[0]:
            best = (score, quartet)
    score, quartet = best if best is not None else (math.inf, ("NA", "NA", "NA", "NA"))
    return slots, quartet, score, filtered


def confident_loci(unique_support: dict[str, object], chi_r: float, min_fraction: float) -> tuple[set[str], float]:
    threshold = min_fraction
    if threshold < 0:
        threshold = max(0.02, min(0.05, (chi_r / 2.0) * 0.5))
    fractions = unique_support.get("fractions", {})
    loci = {locus for locus in DRB345_LOCI if fractions.get(locus, 0.0) >= threshold}
    return loci, threshold


def format_fraction(value: float) -> str:
    return f"{value:.6f}"


def format_count(value: float) -> str:
    return f"{value:.2f}"


def write_calls(call_dir: Path, slots, quartet, tf_counts, total: float, tf_to_safe, safe2name) -> list[dict[str, str]]:
    call_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for idx, (allele_2f, (slot, side, locus, hap_fraction, _candidates)) in enumerate(zip(quartet, slots), 1):
        if allele_2f == "NA":
            full = "NA"
            read_count = 0.0
        else:
            full = safe2name.get(tf_to_safe.get(allele_2f, ""), allele_2f)
            read_count = tf_counts.get(allele_2f, 0.0)
        read_fraction = read_count / total if total else 0.0
        rows.append({
            "global_hap": str(idx),
            "assignment": side,
            "allele": full,
            "drb1_linked_locus": locus_of_allele(allele_2f) if locus == "evidence" else locus,
            "hap_fraction": format_fraction(hap_fraction),
            "allele_read_fraction": format_fraction(read_fraction),
            "allele_read_count": format_count(read_count),
            "em_weight": format_count(read_count),
        })
    fields = ["global_hap", "assignment", "allele", "drb1_linked_locus", "hap_fraction", "allele_read_fraction", "allele_read_count", "em_weight"]
    write_tsv(call_dir / "calls.tsv", fields, rows)
    return rows


def update_final_calls(final_calls: Path, rows: list[dict[str, str]], gmap, source: str, warning: str) -> None:
    existing = read_tsv(final_calls)
    fields = list(existing[0].keys()) if existing else [
        "sample", "gene", "R1_full", "R2_full", "D1_full", "D2_full",
        "R1_2field", "R2_2field", "D1_2field", "D2_2field",
        "R1_g_group", "R2_g_group", "D1_g_group", "D2_g_group",
        "R1_report", "R2_report", "D1_report", "D2_report",
        "R1_fraction", "R2_fraction", "D1_fraction", "D2_fraction",
        "R1_read_fraction", "R2_read_fraction", "D1_read_fraction", "D2_read_fraction",
        "R1_read_count", "R2_read_count", "D1_read_count", "D2_read_count",
        "source", "mean_mask_fraction", "report_level", "warning",
    ]
    by_hap = {row["global_hap"]: row for row in rows}
    alleles = [by_hap.get(str(i), {}).get("allele", "NA") for i in range(1, 5)]
    out = {
        "sample": final_calls.stem.replace(".final_calls", ""),
        "gene": "HLA-DRB345",
        "source": source,
        "mean_mask_fraction": "NA",
        "report_level": "linked-em",
        "warning": warning,
    }
    for slot, allele in zip(("R1", "R2", "D1", "D2"), alleles):
        out[f"{slot}_full"] = allele
        out[f"{slot}_2field"] = allele_2field(allele)
        out[f"{slot}_g_group"] = allele_g_group(allele, gmap)
        out[f"{slot}_report"] = allele_2field(allele)
    for idx, slot in enumerate(("R1", "R2", "D1", "D2"), 1):
        row = by_hap.get(str(idx), {})
        out[f"{slot}_fraction"] = row.get("hap_fraction", "NA")
        out[f"{slot}_read_fraction"] = row.get("allele_read_fraction", "NA")
        out[f"{slot}_read_count"] = row.get("allele_read_count", "NA")
    kept = [row for row in existing if row.get("gene") != "HLA-DRB345"]
    kept.append(out)
    write_tsv(final_calls, fields, kept)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq-dir", required=True, type=Path)
    parser.add_argument("--db-bam", required=True, type=Path)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--final-calls", required=True, type=Path)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--chi-r", type=float, default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--subs-per-2field", type=int, default=5)
    parser.add_argument("--top-per-locus", type=int, default=8)
    parser.add_argument("--db-min-as-frac", type=float, default=0.90)
    parser.add_argument("--db-min-as", type=int, default=DISABLE_MIN_AS)
    parser.add_argument("--remap-min-as-frac", type=float, default=0.95)
    parser.add_argument("--evidence-k", type=int, default=71)
    parser.add_argument("--min-locus-unique-frac", type=float, default=-1.0)
    parser.add_argument("--drb1-untrusted-mask", type=float, default=0.50)
    parser.add_argument("--em-iter", type=int, default=300)
    parser.add_argument("--em-T", type=float, default=2.0)
    args = parser.parse_args()

    out_dir = args.out_dir or (args.fq_dir / "drb345")
    out_dir.mkdir(parents=True, exist_ok=True)
    chi_r = args.chi_r if args.chi_r is not None else read_chi_r(args.fq_dir, args.sample)
    if chi_r is None:
        raise SystemExit("could not infer chi_R; pass --chi-r or keep chimerism logs in fq-dir")

    drb1_row = read_drb1_row(args.final_calls)
    drb1 = {slot: drb1_row.get(f"{slot}_full", "NA") for slot, _side in HAP_SLOTS} if drb1_row else {slot: "NA" for slot, _side in HAP_SLOTS}
    read_names = collect_drb345_read_names(args.db_bam, args.db_min_as_frac, args.db_min_as)
    fq1 = args.fq_dir / f"{args.sample}.uniq.R1.fq.gz"
    fq2 = args.fq_dir / f"{args.sample}.uniq.R2.fq.gz"
    drb_fq1 = out_dir / "DRB345.R1.fq.gz"
    drb_fq2 = out_dir / "DRB345.R2.fq.gz"
    n1 = write_selected_fastq(fq1, drb_fq1, read_names)
    n2 = write_selected_fastq(fq2, drb_fq2, read_names)

    contigs, safe2name, tf_to_safe = build_drb345_candidates(args.imgt, args.subs_per_2field)
    ref_fa = out_dir / "HLA-DRB345.aug.fa"
    build_aug_ref(str(ref_fa), contigs)
    sam = out_dir / "HLA-DRB345.aug.sam"
    if n1 and n2:
        bwa_mem_all(str(ref_fa), str(drb_fq1), str(drb_fq2), args.sample, args.threads, str(sam))
        reads = parse_sam_to_reads(str(sam), set(contigs), args.remap_min_as_frac)
        _theta, counts, iters = run_em(reads, list(contigs), n_iter=args.em_iter, T=args.em_T)
    else:
        reads = {}
        counts = {}
        iters = 0
    tf_counts: dict[str, float] = defaultdict(float)
    for contig, count in counts.items():
        tf_counts[two_field(safe2name[contig])] += count
    total = sum(tf_counts.values())

    tf_path = out_dir / "HLA-DRB345.tf_counts.tsv"
    with tf_path.open("w") as handle:
        handle.write("allele_2field\tlocus\tem_weight\tfraction\n")
        for allele, count in sorted(tf_counts.items(), key=lambda item: (-item[1], item[0])):
            handle.write(f"{allele}\t{locus_of_allele(allele)}\t{count:.4f}\t{(count / total if total else 0.0):.8f}\n")

    unique_support = locus_unique_pair_support(args.imgt, drb_fq1, drb_fq2, args.evidence_k) if n1 and n2 else {
        "total_pairs": 0,
        "counts": {locus: 0 for locus in DRB345_LOCI},
        "fractions": {locus: 0.0 for locus in DRB345_LOCI},
        "unique_kmers": {locus: 0 for locus in DRB345_LOCI},
        "multi_pairs": 0,
    }
    allowed_loci, locus_threshold = confident_loci(unique_support, chi_r, args.min_locus_unique_frac)
    mode = "linked"
    score_counts = dict(tf_counts)
    if drb1_is_untrusted(drb1_row, args.drb1_untrusted_mask) and allowed_loci:
        slots, quartet, score, score_counts = fit_evidence_calls(dict(tf_counts), chi_r, args.top_per_locus, allowed_loci)
        mode = "evidence"
    else:
        slots, quartet, score = fit_linked_calls(dict(tf_counts), drb1, chi_r, args.top_per_locus)
    call_dir = args.asm_root / args.sample / "hla-drb345" / "HLA-DRB345"
    call_rows = write_calls(call_dir, slots, quartet, score_counts, sum(score_counts.values()), tf_to_safe, safe2name)
    gmap = load_g_group(args.g_group)
    source = "drb345-evidence-em" if mode == "evidence" else "drb345-linked-em"
    warning = "drb1_untrusted_drb345_evidence" if mode == "evidence" else "drb1_linked_drb345"
    update_final_calls(args.final_calls, call_rows, gmap, source, warning)

    manifest = out_dir / "HLA-DRB345.manifest.tsv"
    manifest_fields = ["sample", "chi_r", "mode", "db_target_read_names", "fq1_pairs", "fq2_pairs", "em_reads", "em_iters", "sum_abs_diff", "drb1_untrusted", "locus_unique_k", "locus_unique_threshold", "allowed_loci", "unique_pairs_DRB3", "unique_pairs_DRB4", "unique_pairs_DRB5", "unique_frac_DRB3", "unique_frac_DRB4", "unique_frac_DRB5", "drb1_R1", "drb1_R2", "drb1_D1", "drb1_D2", "quartet"]
    unique_counts = unique_support.get("counts", {})
    unique_fracs = unique_support.get("fractions", {})
    manifest_row = {
        "sample": args.sample,
        "chi_r": f"{chi_r:.6f}",
        "mode": mode,
        "db_target_read_names": str(len(read_names)),
        "fq1_pairs": str(n1),
        "fq2_pairs": str(n2),
        "em_reads": str(len(reads)),
        "em_iters": str(iters),
        "sum_abs_diff": f"{score:.6f}" if math.isfinite(score) else "NA",
        "drb1_untrusted": "1" if drb1_is_untrusted(drb1_row, args.drb1_untrusted_mask) else "0",
        "locus_unique_k": str(args.evidence_k),
        "locus_unique_threshold": f"{locus_threshold:.6f}",
        "allowed_loci": ",".join(sorted(allowed_loci)) if allowed_loci else "NA",
        "unique_pairs_DRB3": str(unique_counts.get("DRB3", 0)),
        "unique_pairs_DRB4": str(unique_counts.get("DRB4", 0)),
        "unique_pairs_DRB5": str(unique_counts.get("DRB5", 0)),
        "unique_frac_DRB3": f"{unique_fracs.get('DRB3', 0.0):.6f}",
        "unique_frac_DRB4": f"{unique_fracs.get('DRB4', 0.0):.6f}",
        "unique_frac_DRB5": f"{unique_fracs.get('DRB5', 0.0):.6f}",
        "drb1_R1": drb1.get("R1", "NA"),
        "drb1_R2": drb1.get("R2", "NA"),
        "drb1_D1": drb1.get("D1", "NA"),
        "drb1_D2": drb1.get("D2", "NA"),
        "quartet": ",".join(quartet),
    }
    write_tsv(manifest, manifest_fields, [manifest_row])
    if sam.exists():
        sam.unlink()
    print(f"[drb345] wrote {call_dir / 'calls.tsv'} and updated {args.final_calls}")


if __name__ == "__main__":
    main()