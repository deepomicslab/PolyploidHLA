#!/usr/bin/env python3
"""Direct read evidence gate for competing HLA quartets.

This offline validator compares a baseline quartet against a validation/rescue
quartet using per-gene read pairs and allele-sequence k-mer evidence. It does
not call variants, phase haplotypes, or use truth for the decision. Truth-derived
status columns from a score TSV are used only to summarize whether the likelihood
direction would have accepted improvements and rejected regressions.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from hla_polyphase_assemble import load_imgt_alleles  # noqa: E402


DEFAULT_IMGT = SCRIPT_ROOT / "resources" / "spechla" / "db" / "ref" / "hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
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
    return gene if gene.startswith("HLA-") else f"HLA-{gene}"


def short_gene(gene: str) -> str:
    return GENE_TO_SHORT.get(normalize_gene(gene), normalize_gene(gene).replace("HLA-", ""))


def strip_expr_suffix(value: str) -> str:
    return value[:-1] if value and value[-1].isalpha() and value[-1] != "G" else value


def clean_allele(allele: str) -> str:
    allele = (allele or "").strip().replace("HLA-", "")
    if not allele or allele == "NA" or "*" not in allele:
        return ""
    gene, fields = allele.split("*", 1)
    parts = fields.split(":")
    parts[-1] = strip_expr_suffix(parts[-1])
    return f"{gene}*{':'.join(parts)}"


def allele_2field(allele: str) -> str:
    allele = clean_allele(allele)
    if not allele:
        return ""
    gene, fields = allele.split("*", 1)
    parts = fields.replace("G", "").split(":")
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def split_alleles(text: str) -> list[str]:
    return [allele_2field(item) for item in (text or "").split(",") if allele_2field(item)]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def iter_fastq_pairs(fq1: Path, fq2: Path):
    with open_text(fq1) as handle1, open_text(fq2) as handle2:
        while True:
            h1 = handle1.readline()
            h2 = handle2.readline()
            if not h1 or not h2:
                return
            s1 = handle1.readline().strip().upper()
            s2 = handle2.readline().strip().upper()
            handle1.readline()
            handle2.readline()
            q1 = handle1.readline()
            q2 = handle2.readline()
            if not q1 or not q2:
                return
            yield s1, s2


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def iter_kmers(seq: str, k: int) -> Iterable[str]:
    if len(seq) < k:
        return
    for start in range(len(seq) - k + 1):
        kmer = seq[start:start + k]
        if "N" not in kmer:
            yield kmer


def allele_prefix_regex(gene: str, allele: str) -> re.Pattern[str]:
    gene_prefix = normalize_gene(gene).replace("HLA-", "")
    _allele_gene, fields = allele.split("*", 1)
    return re.compile(rf"^{re.escape(gene_prefix)}\*{re.escape(fields)}(?::|$)")


@lru_cache(maxsize=4)
def load_imgt_cached(imgt: str) -> dict[str, str]:
    return load_imgt_alleles(imgt)


def load_candidate_sequences(imgt: Path, gene: str, candidates: list[str]) -> dict[str, list[str]]:
    alleles = load_imgt_cached(str(imgt))
    by_candidate: dict[str, list[str]] = {candidate: [] for candidate in candidates}
    patterns = {candidate: allele_prefix_regex(gene, candidate) for candidate in candidates}
    for full_name, seq in alleles.items():
        clean_name = clean_allele(full_name)
        if not clean_name:
            continue
        for candidate, pattern in patterns.items():
            if pattern.match(clean_name):
                by_candidate[candidate].append(seq.upper())
    return by_candidate


def build_informative_kmers(
    seqs_by_allele: dict[str, list[str]],
    k: int,
    max_full_alleles_per_2field: int,
    max_owner_fraction: float,
) -> dict[str, tuple[str, ...]]:
    kmer_owners: dict[str, set[str]] = defaultdict(set)
    candidates = sorted(seqs_by_allele)
    for allele, seqs in seqs_by_allele.items():
        for seq in seqs[:max_full_alleles_per_2field]:
            for source in (seq, revcomp(seq)):
                for kmer in iter_kmers(source, k):
                    kmer_owners[kmer].add(allele)
    max_owners = max(1, math.floor(len(candidates) * max_owner_fraction))
    out = {}
    for kmer, owners in kmer_owners.items():
        if len(owners) == 0 or len(owners) == len(candidates):
            continue
        if len(owners) > max_owners:
            continue
        out[kmer] = tuple(sorted(owners))
    return out


def read_chi_r(spechla_root: Path, sample: str) -> float:
    sample_dir = spechla_root / sample
    pooled = sample_dir / f"{sample}.chi_pooled.txt"
    if pooled.exists():
        for line in pooled.read_text().splitlines():
            if not line.startswith("GLOBAL") or "chi_R=" not in line:
                continue
            for item in line.split():
                if item.startswith("chi_R="):
                    try:
                        value = float(item.split("=", 1)[1])
                    except ValueError:
                        continue
                    if 0.02 <= value <= 0.98:
                        return value
    chimerism = sample_dir / f"{sample}.chimerism.txt"
    if chimerism.exists():
        for line in chimerism.read_text().splitlines():
            if "chi_R=" not in line:
                continue
            for item in line.split():
                if item.startswith("chi_R="):
                    try:
                        value = float(item.split("=", 1)[1])
                    except ValueError:
                        continue
                    if 0.02 <= value <= 0.98:
                        return value
    return 0.5


def gene_fastqs(spechla_root: Path, sample: str, gene: str) -> tuple[Path, Path]:
    base = spechla_root / sample / short_gene(gene)
    return Path(f"{base}.R1.fq.gz"), Path(f"{base}.R2.fq.gz")


def effective_weights(quartet: list[str], chi_r: float) -> dict[str, float]:
    slot_weights = [chi_r / 2.0, chi_r / 2.0, (1.0 - chi_r) / 2.0, (1.0 - chi_r) / 2.0]
    weights: dict[str, float] = defaultdict(float)
    for allele, weight in zip(quartet, slot_weights):
        weights[allele] += weight
    return dict(weights)


def logsumexp(values: list[float]) -> float:
    max_value = max(values)
    return max_value + math.log(sum(math.exp(value - max_value) for value in values))


def quartet_pair_loglik(allele_counts: Counter[str], quartet: list[str], chi_r: float, score_scale: float) -> float:
    weights = effective_weights(quartet, chi_r)
    values = []
    for allele, weight in weights.items():
        values.append(math.log(max(weight, 1e-12)) + score_scale * allele_counts.get(allele, 0.0))
    return logsumexp(values)


def score_row(row: dict[str, str], args: argparse.Namespace) -> dict[str, object]:
    sample = row["sample"]
    gene = normalize_gene(row["gene"])
    baseline = split_alleles(row["baseline_pred_R"]) + split_alleles(row["baseline_pred_D"])
    validation = split_alleles(row["validation_pred_R"]) + split_alleles(row["validation_pred_D"])
    candidates = sorted(set(baseline) | set(validation))
    result = {
        "sample": sample,
        "gene": gene,
        "status": row.get("status", ""),
        "score_delta": row.get("delta", ""),
        "baseline_quartet": ",".join(baseline),
        "validation_quartet": ",".join(validation),
        "candidate_alleles": ",".join(candidates),
    }
    if len(baseline) != 4 or len(validation) != 4 or not candidates:
        result.update({"decision": "missing_quartet", "error": "missing baseline/validation quartet"})
        return result

    expected = "validation" if row.get("status") == "improved" else "baseline" if row.get("status") == "regressed" else "same"
    if baseline == validation:
        result.update({
            "decision": "baseline",
            "expected_for_changed_row": expected,
            "direction_correct": 1,
            "chi_r": "",
            "total_pairs": 0,
            "informative_pairs": 0,
            "discriminating_pairs": 0,
            "validation_pair_wins": 0,
            "baseline_pair_wins": 0,
            "informative_kmers": 0,
            "mean_pair_informative_kmers": "0.000",
            "baseline_loglik": "0.0000",
            "validation_loglik": "0.0000",
            "validation_minus_baseline": "0.0000",
            "error": "",
        })
        return result

    seqs_by_allele = load_candidate_sequences(args.imgt, gene, candidates)
    missing = [allele for allele, seqs in seqs_by_allele.items() if not seqs]
    if missing:
        result.update({"decision": "missing_allele_sequence", "error": ",".join(missing)})
        return result

    kmer_owners = build_informative_kmers(
        seqs_by_allele,
        args.k,
        args.max_full_alleles_per_2field,
        args.max_owner_fraction,
    )
    fq1, fq2 = gene_fastqs(args.spechla_root, sample, gene)
    if not fq1.exists() or not fq2.exists():
        result.update({"decision": "missing_fastq", "error": f"{fq1},{fq2}"})
        return result

    chi_r = read_chi_r(args.spechla_root, sample)
    baseline_ll = 0.0
    validation_ll = 0.0
    informative_pairs = 0
    discriminating_pairs = 0
    validation_pair_wins = 0
    baseline_pair_wins = 0
    total_pairs = 0
    total_informative_kmers = 0.0
    for seq1, seq2 in iter_fastq_pairs(fq1, fq2):
        total_pairs += 1
        if args.max_pairs > 0 and total_pairs > args.max_pairs:
            break
        allele_counts: Counter[str] = Counter()
        for kmer in set(iter_kmers(seq1, args.k)) | set(iter_kmers(seq2, args.k)):
            owners = kmer_owners.get(kmer)
            if not owners:
                continue
            inc = 1.0 / len(owners)
            for allele in owners:
                allele_counts[allele] += inc
        informative = sum(allele_counts.values())
        if informative < args.min_pair_informative_kmers:
            continue
        informative_pairs += 1
        total_informative_kmers += informative
        pair_baseline = quartet_pair_loglik(allele_counts, baseline, chi_r, args.score_scale)
        pair_validation = quartet_pair_loglik(allele_counts, validation, chi_r, args.score_scale)
        baseline_ll += pair_baseline
        validation_ll += pair_validation
        if abs(pair_validation - pair_baseline) >= args.pair_margin:
            discriminating_pairs += 1
            if pair_validation > pair_baseline:
                validation_pair_wins += 1
            else:
                baseline_pair_wins += 1

    margin = validation_ll - baseline_ll
    decision = "validation" if margin >= args.margin and discriminating_pairs >= args.min_discriminating_pairs else "baseline"
    result.update({
        "decision": decision,
        "expected_for_changed_row": expected,
        "direction_correct": int((expected == "validation" and decision == "validation") or (expected == "baseline" and decision == "baseline") or expected == "same"),
        "chi_r": f"{chi_r:.4f}",
        "total_pairs": total_pairs,
        "informative_pairs": informative_pairs,
        "discriminating_pairs": discriminating_pairs,
        "validation_pair_wins": validation_pair_wins,
        "baseline_pair_wins": baseline_pair_wins,
        "informative_kmers": len(kmer_owners),
        "mean_pair_informative_kmers": f"{(total_informative_kmers / informative_pairs) if informative_pairs else 0.0:.3f}",
        "baseline_loglik": f"{baseline_ll:.4f}",
        "validation_loglik": f"{validation_ll:.4f}",
        "validation_minus_baseline": f"{margin:.4f}",
        "error": "",
    })
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare baseline vs validation HLA quartets from direct read evidence")
    parser.add_argument("--score-tsv", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path,
                        help="root containing <sample>/<gene-short>.R1/R2.fq.gz")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--statuses", default="improved,regressed",
                        help="comma-separated score row statuses to evaluate; use ALL for every row")
    parser.add_argument("--genes", default="",
                        help="optional comma-separated gene filter")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-pairs", type=int, default=0)
    parser.add_argument("--max-full-alleles-per-2field", type=int, default=25)
    parser.add_argument("--max-owner-fraction", type=float, default=0.75)
    parser.add_argument("--score-scale", type=float, default=0.35)
    parser.add_argument("--margin", type=float, default=5.0)
    parser.add_argument("--pair-margin", type=float, default=0.05)
    parser.add_argument("--min-discriminating-pairs", type=int, default=3)
    parser.add_argument("--min-pair-informative-kmers", type=float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    status_filter = None if args.statuses.upper() == "ALL" else {item.strip() for item in args.statuses.split(",") if item.strip()}
    gene_filter = {normalize_gene(item.strip()) for item in args.genes.split(",") if item.strip()}
    rows = []
    for row in read_tsv(args.score_tsv):
        if status_filter is not None and row.get("status") not in status_filter:
            continue
        if gene_filter and normalize_gene(row.get("gene", "")) not in gene_filter:
            continue
        rows.append(score_row(row, args))

    fields = [
        "sample", "gene", "status", "score_delta", "decision", "expected_for_changed_row",
        "direction_correct", "chi_r", "total_pairs", "informative_pairs",
        "discriminating_pairs", "validation_pair_wins", "baseline_pair_wins",
        "informative_kmers", "mean_pair_informative_kmers", "baseline_loglik",
        "validation_loglik", "validation_minus_baseline", "baseline_quartet",
        "validation_quartet", "candidate_alleles", "error",
    ]
    write_tsv(args.out, fields, rows)

    summary_path = args.summary or args.out.with_suffix(".summary.tsv")
    total = len(rows)
    decided = [row for row in rows if row.get("decision") in {"baseline", "validation"}]
    correct = sum(int(row.get("direction_correct", 0)) for row in decided)
    improved = [row for row in decided if row.get("status") == "improved"]
    regressed = [row for row in decided if row.get("status") == "regressed"]
    summary_rows = [{
        "rows": total,
        "decided_rows": len(decided),
        "direction_correct": correct,
        "direction_accuracy": f"{(correct / len(decided)) if decided else 0.0:.4f}",
        "improved_accept_validation": sum(1 for row in improved if row.get("decision") == "validation"),
        "improved_total": len(improved),
        "regressed_keep_baseline": sum(1 for row in regressed if row.get("decision") == "baseline"),
        "regressed_total": len(regressed),
        "errors": sum(1 for row in rows if row.get("error")),
    }]
    write_tsv(summary_path, list(summary_rows[0]), summary_rows)
    print(f"wrote {args.out}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()