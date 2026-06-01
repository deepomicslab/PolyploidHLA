#!/usr/bin/env python3
"""Prototype raw-read to allele-database quartet caller.

This script is an offline experiment. It keeps chimerism from the existing
pipeline outputs, but derives allele candidates and quartets directly from read
pairs and the IMGT allele database. It does not use baseline predictions for
candidate generation or scoring, and it does not call variants or phase VCFs.
Truth columns in the input score TSV are used only for validation metrics.
"""
from __future__ import annotations

import argparse
import math
import sys
import subprocess
import tempfile
from collections import Counter, defaultdict
from itertools import combinations_with_replacement, permutations
from pathlib import Path
from typing import Callable, Mapping

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from direct_read_quartet_likelihood import (  # noqa: E402
    DEFAULT_IMGT,
    allele_2field,
    clean_allele,
    effective_weights,
    gene_fastqs,
    iter_fastq_pairs,
    iter_kmers,
    load_imgt_cached,
    logsumexp,
    normalize_gene,
    open_text,
    quartet_pair_loglik,
    read_chi_r,
    read_tsv,
    revcomp,
    split_alleles,
    write_tsv,
)
from score_gendx_validation_root import score_quartet  # noqa: E402


BACKGROUND_GENES = {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"}

BOOL_TRUE = {"1", "true", "yes", "y", "on"}
BOOL_FALSE = {"0", "false", "no", "n", "off"}


def parse_bool_config(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in BOOL_TRUE:
        return True
    if normalized in BOOL_FALSE:
        return False
    raise ValueError(f"expected boolean value, got {value!r}")


def parse_optional_float_config(value: str) -> float | None:
    normalized = value.strip().lower()
    if normalized in {"", "none", "null", "na", "nan"}:
        return None
    return float(value)


GENE_CONFIG_FIELDS: dict[str, tuple[str, Callable[[str], object]]] = {
    "gene_action": ("gene_action", str),
    "selector_iterations": ("selector_iterations", int),
    "evidence_backend": ("evidence_backend", str),
    "max_candidates": ("max_candidates", int),
    "min_candidate_support": ("min_candidate_support", float),
    "max_candidate_pairs": ("max_candidate_pairs", int),
    "max_evidence_pairs": ("max_evidence_pairs", int),
    "alignment_candidate_ranker": ("alignment_candidate_ranker", str),
    "alignment_count_prior_weight": ("alignment_count_prior_weight", float),
    "alignment_top_per_pair": ("alignment_top_per_pair", int),
    "alignment_top_delta": ("alignment_top_delta", float),
    "alignment_full_allele_collapse": ("alignment_full_allele_collapse", parse_bool_config),
    "alignment_missing_penalty": ("alignment_missing_penalty", float),
    "quartet_selector": ("quartet_selector", str),
    "ratio_count_weight": ("ratio_count_weight", float),
    "ratio_count_epsilon": ("ratio_count_epsilon", float),
    "ratio_count_use_ambiguity_groups": ("ratio_count_use_ambiguity_groups", parse_bool_config),
    "constrained_ratio_min_dropped_group_observed": ("constrained_ratio_min_dropped_group_observed", float),
    "constrained_ratio_max_dropped_group_copies": ("constrained_ratio_max_dropped_group_copies", int),
    "baseline_gate_margin": ("baseline_gate_margin", parse_optional_float_config),
    "baseline_confidence_max_loglik_per_pair": ("baseline_confidence_max_loglik_per_pair", parse_optional_float_config),
    "baseline_multiset_gate": ("baseline_multiset_gate", parse_bool_config),
    "side_copy_gate": ("side_copy_gate", str),
    "reject_balanced_het_to_hom_side": ("reject_balanced_het_to_hom_side", parse_bool_config),
    "ambiguity_collapse_gate": ("ambiguity_collapse_gate", str),
    "ambiguity_collapse_min_pairs": ("ambiguity_collapse_min_pairs", int),
    "ambiguity_collapse_min_shared_fraction": ("ambiguity_collapse_min_shared_fraction", float),
    "ambiguity_collapse_max_mean_delta": ("ambiguity_collapse_max_mean_delta", float),
    "ambiguity_collapse_max_pair_delta": ("ambiguity_collapse_max_pair_delta", float),
    "ambiguity_collapse_min_within_pair_delta_fraction": ("ambiguity_collapse_min_within_pair_delta_fraction", float),
    "ambiguity_collapse_adaptive_max_mean_fraction": ("ambiguity_collapse_adaptive_max_mean_fraction", float),
    "ambiguity_collapse_adaptive_max_pair_fraction": ("ambiguity_collapse_adaptive_max_pair_fraction", float),
    "ambiguity_collapse_adaptive_min_mean_gap": ("ambiguity_collapse_adaptive_min_mean_gap", float),
    "replacement_gate": ("replacement_gate", str),
    "replacement_min_loglik_per_pair": ("replacement_min_loglik_per_pair", float),
    "replacement_min_support": ("replacement_min_support", float),
    "replacement_min_em_weight": ("replacement_min_em_weight", float),
    "replacement_max_new_copies": ("replacement_max_new_copies", int),
    "replacement_allow_same_ambiguity_multiset": ("replacement_allow_same_ambiguity_multiset", parse_bool_config),
    "replacement_allow_balanced_group_dosage_change": ("replacement_allow_balanced_group_dosage_change", parse_bool_config),
    "residual_gate": ("residual_gate", str),
    "residual_min_improvement_per_pair": ("residual_min_improvement_per_pair", float),
    "residual_missing_penalty": ("residual_missing_penalty", float),
    "residual_indel_weight": ("residual_indel_weight", float),
    "residual_model": ("residual_model", str),
    "residual_scale": ("residual_scale", float),
    "raw_gene_enrich": ("raw_gene_enrich", parse_bool_config),
    "raw_gene_enrich_gene_unique": ("raw_gene_enrich_gene_unique", parse_bool_config),
    "raw_gene_enrich_min_kmers": ("raw_gene_enrich_min_kmers", float),
    "raw_gene_enrich_scan_pairs": ("raw_gene_enrich_scan_pairs", int),
    "em_iterations": ("em_iterations", int),
    "em_min_weight": ("em_min_weight", float),
    "max_quartets": ("max_quartets", int),
}


GENE_CONFIG_CHOICES: dict[str, set[str]] = {
    "gene_action": {"call", "baseline"},
    "evidence_backend": {"kmer", "alignment"},
    "alignment_candidate_ranker": {"em", "support", "support-count-prior"},
    "quartet_selector": {"read-likelihood", "ratio-count", "constrained-ratio", "support-count", "multiset-likelihood"},
    "side_copy_gate": {"off", "multiset", "no-new-2field"},
    "ambiguity_collapse_gate": {"off", "group-side", "group-multiset", "adaptive-side", "adaptive-multiset"},
    "replacement_gate": {"off", "likelihood"},
    "residual_gate": {"off", "not-worse", "improvement"},
    "residual_model": {"min", "weighted"},
}


def gene_prefix(gene: str) -> str:
    return normalize_gene(gene).replace("HLA-", "") + "*"


def load_gene_representatives(imgt: Path, gene: str) -> dict[str, str]:
    prefix = gene_prefix(gene)
    representatives: dict[str, str] = {}
    for full_name, seq in load_imgt_cached(str(imgt)).items():
        clean_name = clean_allele(full_name)
        if not clean_name.startswith(prefix):
            continue
        allele = allele_2field(clean_name)
        if not allele:
            continue
        seq = seq.upper().replace("-", "")
        if not seq:
            continue
        if allele not in representatives or len(seq) > len(representatives[allele]):
            representatives[allele] = seq
    return representatives


def load_gene_full_allele_representatives(imgt: Path, gene: str) -> dict[str, str]:
    prefix = gene_prefix(gene)
    representatives: dict[str, str] = {}
    for full_name, seq in load_imgt_cached(str(imgt)).items():
        clean_name = clean_allele(full_name)
        if not clean_name.startswith(prefix):
            continue
        seq = seq.upper().replace("-", "")
        if seq:
            representatives[clean_name] = seq
    return representatives


def full_allele_record_counts(full_representatives: Mapping[str, str]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for full_allele in full_representatives:
        allele = allele_2field(full_allele)
        if allele:
            counts[allele] += 1
    return dict(counts)


def load_all_gene_representatives(imgt: Path, genes: set[str]) -> dict[str, dict[str, str]]:
    prefixes = {normalize_gene(gene).replace("HLA-", ""): normalize_gene(gene) for gene in genes}
    representatives: dict[str, dict[str, str]] = {normalize_gene(gene): {} for gene in genes}
    for full_name, seq in load_imgt_cached(str(imgt)).items():
        clean_name = clean_allele(full_name)
        if "*" not in clean_name:
            continue
        name_gene = clean_name.split("*", 1)[0]
        gene = prefixes.get(name_gene)
        if gene is None:
            continue
        allele = allele_2field(clean_name)
        seq = seq.upper().replace("-", "")
        if not allele or not seq:
            continue
        gene_reps = representatives[gene]
        if allele not in gene_reps or len(seq) > len(gene_reps[allele]):
            gene_reps[allele] = seq
    return representatives


def canonical_config_key(column: str) -> str:
    return column.strip().lstrip("-").replace("-", "_")


def serialize_gene_config(overrides: Mapping[str, object]) -> str:
    return ";".join(f"{key}={value}" for key, value in sorted(overrides.items()))


def load_gene_config_profiles(path: Path | None) -> dict[str, dict[str, object]]:
    if path is None:
        return {}
    profiles: dict[str, dict[str, object]] = {}
    for row in read_tsv(path):
        gene = normalize_gene(row.get("gene", ""))
        if not gene:
            raise ValueError(f"gene profile row in {path} is missing a gene")
        overrides: dict[str, object] = {}
        for raw_column, raw_value in row.items():
            if raw_column is None:
                continue
            column = canonical_config_key(raw_column)
            if column in {"gene", "profile", "note", "notes", "comment"}:
                continue
            value = (raw_value or "").strip()
            if value == "":
                continue
            if column not in GENE_CONFIG_FIELDS:
                allowed = ", ".join(sorted(GENE_CONFIG_FIELDS))
                raise ValueError(f"unknown gene config column {raw_column!r} in {path}; allowed columns: {allowed}")
            attr, converter = GENE_CONFIG_FIELDS[column]
            converted = converter(value)
            allowed_choices = GENE_CONFIG_CHOICES.get(attr)
            if allowed_choices is not None and str(converted) not in allowed_choices:
                choices = ", ".join(sorted(allowed_choices))
                raise ValueError(f"invalid {raw_column}={value!r} for {gene}; allowed: {choices}")
            overrides[attr] = converted
        if gene in profiles:
            raise ValueError(f"duplicate gene profile for {gene} in {path}")
        profiles[gene] = overrides
    return profiles


def args_for_gene(args: argparse.Namespace, gene: str) -> argparse.Namespace:
    overrides = args.gene_config_profiles.get(normalize_gene(gene), {})
    if not overrides:
        return args
    gene_args = argparse.Namespace(**vars(args))
    for attr, value in overrides.items():
        setattr(gene_args, attr, value)
    gene_args.gene_config_overrides = serialize_gene_config(overrides)
    return gene_args


def build_gene_kmer_owners(
    representatives: dict[str, str],
    other_gene_representatives: list[dict[str, str]],
    k: int,
    max_owner_alleles: int,
    gene_unique: bool,
) -> dict[str, tuple[str, ...]]:
    owner_sets: dict[str, set[str]] = defaultdict(set)
    for allele, seq in representatives.items():
        for source in (seq, revcomp(seq)):
            for kmer in iter_kmers(source, k):
                owners = owner_sets[kmer]
                if len(owners) <= max_owner_alleles:
                    owners.add(allele)
    blocked_kmers: set[str] = set()
    if gene_unique:
        target_kmers = set(owner_sets)
        for other_reps in other_gene_representatives:
            for seq in other_reps.values():
                for source in (seq, revcomp(seq)):
                    for kmer in iter_kmers(source, k):
                        if kmer in target_kmers:
                            blocked_kmers.add(kmer)
    return {
        kmer: tuple(sorted(owners))
        for kmer, owners in owner_sets.items()
        if kmer not in blocked_kmers and 0 < len(owners) <= max_owner_alleles
    }


def build_candidate_kmer_owners(
    representatives: dict[str, str],
    candidates: list[str],
    k: int,
    max_owner_fraction: float,
) -> dict[str, tuple[str, ...]]:
    owner_sets: dict[str, set[str]] = defaultdict(set)
    candidate_set = set(candidates)
    for allele in candidates:
        seq = representatives.get(allele)
        if not seq:
            continue
        for source in (seq, revcomp(seq)):
            for kmer in iter_kmers(source, k):
                owner_sets[kmer].add(allele)
    max_owners = max(1, math.floor(len(candidate_set) * max_owner_fraction))
    return {
        kmer: tuple(sorted(owners))
        for kmer, owners in owner_sets.items()
        if 0 < len(owners) < len(candidate_set) and len(owners) <= max_owners
    }


def build_seed_kmer_owners(
    representatives: dict[str, str],
    k: int,
    max_owner_alleles: int,
) -> dict[str, tuple[str, ...]]:
    owner_sets: dict[str, set[str]] = defaultdict(set)
    for allele, seq in representatives.items():
        for source in (seq, revcomp(seq)):
            for kmer in iter_kmers(source, k):
                owners = owner_sets[kmer]
                if len(owners) <= max_owner_alleles:
                    owners.add(allele)
    return {
        kmer: tuple(sorted(owners))
        for kmer, owners in owner_sets.items()
        if 0 < len(owners) <= max_owner_alleles
    }


def raw_fastqs(fq_root: Path, sample: str) -> tuple[Path, Path]:
    return fq_root / f"{sample}_R1_001.fastq.gz", fq_root / f"{sample}_R2_001.fastq.gz"


def select_fastqs(args: argparse.Namespace, sample: str, gene: str) -> tuple[Path, Path, str]:
    if args.fq_root is not None:
        fq1, fq2 = raw_fastqs(args.fq_root, sample)
        return fq1, fq2, "raw"
    fq1, fq2 = gene_fastqs(args.spechla_root, sample, gene)
    return fq1, fq2, "gene_bin"


def scan_candidate_support(
    fq1: Path,
    fq2: Path,
    kmer_owners: dict[str, tuple[str, ...]],
    k: int,
    max_pairs: int,
) -> tuple[int, Counter[str]]:
    support: Counter[str] = Counter()
    total_pairs = 0
    for seq1, seq2 in iter_fastq_pairs(fq1, fq2):
        total_pairs += 1
        if max_pairs > 0 and total_pairs > max_pairs:
            break
        pair_alleles = Counter()
        for kmer in set(iter_kmers(seq1, k)) | set(iter_kmers(seq2, k)):
            owners = kmer_owners.get(kmer)
            if not owners:
                continue
            increment = 1.0 / len(owners)
            for allele in owners:
                pair_alleles[allele] += increment
        if not pair_alleles:
            continue
        for allele, value in pair_alleles.items():
            support[allele] += value
    return total_pairs, support


def collect_evidence_pairs(
    fq1: Path,
    fq2: Path,
    kmer_owners: dict[str, tuple[str, ...]],
    k: int,
    max_pairs: int,
    min_pair_informative_kmers: float,
) -> tuple[int, list[Counter[str]]]:
    total_pairs = 0
    evidence: list[Counter[str]] = []
    for seq1, seq2 in iter_fastq_pairs(fq1, fq2):
        total_pairs += 1
        if max_pairs > 0 and total_pairs > max_pairs:
            break
        allele_counts: Counter[str] = Counter()
        for kmer in set(iter_kmers(seq1, k)) | set(iter_kmers(seq2, k)):
            owners = kmer_owners.get(kmer)
            if not owners:
                continue
            increment = 1.0 / len(owners)
            for allele in owners:
                allele_counts[allele] += increment
        if sum(allele_counts.values()) >= min_pair_informative_kmers:
            evidence.append(allele_counts)
    return total_pairs, evidence


class ReadPairAligner:
    def __init__(self) -> None:
        import parasail

        self.parasail = parasail
        self.matrix = parasail.matrix_create("ACGT", 2, -3)
        self.gap_open = 5
        self.gap_extend = 2
        self._score_cache: dict[tuple[str, str], int] = {}

    def score_read(self, read: str, allele: str, ref: str) -> int:
        cache_key = (read, allele)
        cached = self._score_cache.get(cache_key)
        if cached is not None:
            return cached
        score = max(
            self.parasail.sg_dx_trace_striped_16(read, ref, self.gap_open, self.gap_extend, self.matrix).score,
            self.parasail.sg_dx_trace_striped_16(revcomp(read), ref, self.gap_open, self.gap_extend, self.matrix).score,
        )
        self._score_cache[cache_key] = score
        return score

    def score_pair(self, seq1: str, seq2: str, allele: str, ref: str) -> int:
        return self.score_read(seq1, allele, ref) + self.score_read(seq2, allele, ref)


class MappyReadPairAligner:
    def __init__(self, representatives: Mapping[str, str], best_n: int) -> None:
        import mappy

        self.mappy = mappy
        self._tmpdir = tempfile.TemporaryDirectory(prefix="raw_read_allele_mappy_")
        fasta = Path(self._tmpdir.name) / "alleles.fa"
        with fasta.open("w") as handle:
            for allele, seq in representatives.items():
                handle.write(f">{allele}\n{seq}\n")
        self.aligner = mappy.Aligner(str(fasta), preset="sr", best_n=max(1, best_n))
        if not self.aligner:
            raise RuntimeError(f"failed to build mappy index for {fasta}")

    def close(self) -> None:
        self._tmpdir.cleanup()

    def score_read(self, read: str) -> dict[str, int]:
        scores: dict[str, int] = {}
        for hit in self.aligner.map(read):
            score = hit.mlen - 2 * hit.NM
            if hit.ctg not in scores or score > scores[hit.ctg]:
                scores[hit.ctg] = score
        return scores

    def score_pair(self, seq1: str, seq2: str) -> dict[str, int]:
        scores = self.score_read(seq1)
        for allele, score in self.score_read(seq2).items():
            scores[allele] = scores.get(allele, 0) + score
        return scores


def top_alignment_loglikes(
    seq1: str,
    seq2: str,
    representatives: Mapping[str, str],
    seed_kmer_owners: Mapping[str, tuple[str, ...]],
    aligner: ReadPairAligner,
    seed_k: int,
    prefilter_alleles: int,
    top_per_pair: int,
    top_delta: float,
    score_scale: float,
) -> tuple[dict[str, float], int]:
    seed_support: Counter[str] = Counter()
    for kmer in set(iter_kmers(seq1, seed_k)) | set(iter_kmers(seq2, seed_k)):
        owners = seed_kmer_owners.get(kmer)
        if not owners:
            continue
        increment = 1.0 / len(owners)
        for allele in owners:
            seed_support[allele] += increment
    candidate_alleles = [allele for allele, _ in seed_support.most_common(prefilter_alleles)]
    if not candidate_alleles:
        return {}, 0
    scores = []
    for allele in candidate_alleles:
        ref = representatives.get(allele)
        if ref is None:
            continue
        scores.append((allele, aligner.score_pair(seq1, seq2, allele, ref)))
    if not scores:
        return {}, 0
    scores.sort(key=lambda item: -item[1])
    best = scores[0][1]
    kept: dict[str, float] = {}
    for allele, score in scores[:top_per_pair]:
        if best - score > top_delta:
            continue
        kept[allele] = (score - best) / max(score_scale, 1e-6)
    return kept, best


def collect_alignment_evidence(
    fq1: Path,
    fq2: Path,
    representatives: Mapping[str, str],
    seed_kmer_owners: Mapping[str, tuple[str, ...]],
    max_pairs: int,
    max_evidence: int,
    seed_k: int,
    prefilter_alleles: int,
    top_per_pair: int,
    top_delta: float,
    score_scale: float,
    min_score_per_base: float,
) -> tuple[int, list[dict[str, float]], Counter[str]]:
    aligner = ReadPairAligner()
    total_pairs = 0
    evidence: list[dict[str, float]] = []
    support: Counter[str] = Counter()
    for seq1, seq2 in iter_fastq_pairs(fq1, fq2):
        total_pairs += 1
        if max_pairs > 0 and total_pairs > max_pairs:
            break
        loglikes, best_score = top_alignment_loglikes(
            seq1,
            seq2,
            representatives,
            seed_kmer_owners,
            aligner,
            seed_k,
            prefilter_alleles,
            top_per_pair,
            top_delta,
            score_scale,
        )
        min_score = min_score_per_base * max(1, len(seq1) + len(seq2))
        if not loglikes or best_score < min_score:
            continue
        best_loglike = max(loglikes.values())
        weights = {allele: math.exp(value - best_loglike) for allele, value in loglikes.items()}
        total_weight = sum(weights.values()) or 1.0
        for allele, weight in weights.items():
            support[allele] += weight / total_weight
        evidence.append(loglikes)
        if max_evidence > 0 and len(evidence) >= max_evidence:
            break
    return total_pairs, evidence, support


def write_allele_fasta(representatives: Mapping[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for allele, seq in representatives.items():
            handle.write(f">{allele}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def fastq_qname(header: str) -> str:
    return header.strip()[1:].split()[0]


def pair_kmer_support(seq1: str, seq2: str, kmer_owners: Mapping[str, tuple[str, ...]], k: int) -> float:
    support = 0.0
    for kmer in set(iter_kmers(seq1, k)) | set(iter_kmers(seq2, k)):
        owners = kmer_owners.get(kmer)
        if owners:
            support += 1.0 / len(owners)
    return support


def write_fastq_subset(
    fq1: Path,
    fq2: Path,
    out1: Path,
    out2: Path,
    max_pairs: int,
    enrich_kmer_owners: Mapping[str, tuple[str, ...]] | None = None,
    enrich_k: int = 31,
    min_enrich_kmers: float = 1.0,
    scan_pairs: int = 0,
) -> tuple[int, int, dict[str, int]]:
    selected_pairs = 0
    scanned_pairs = 0
    read_lengths: dict[str, int] = {}
    with open_text(fq1) as handle1, open_text(fq2) as handle2, out1.open("w") as out_handle1, out2.open("w") as out_handle2:
        while True:
            record1 = [handle1.readline() for _ in range(4)]
            record2 = [handle2.readline() for _ in range(4)]
            if not record1[0] or not record2[0]:
                break
            if any(not item for item in record1 + record2):
                break
            scanned_pairs += 1
            if scan_pairs > 0 and scanned_pairs > scan_pairs:
                scanned_pairs -= 1
                break
            if enrich_kmer_owners is not None:
                seq1 = record1[1].strip().upper()
                seq2 = record2[1].strip().upper()
                if pair_kmer_support(seq1, seq2, enrich_kmer_owners, enrich_k) < min_enrich_kmers:
                    continue
            selected_pairs += 1
            if max_pairs > 0 and selected_pairs > max_pairs:
                selected_pairs -= 1
                break
            out_handle1.writelines(record1)
            out_handle2.writelines(record2)
            qname = fastq_qname(record1[0])
            read_lengths[qname] = len(record1[1].strip()) + len(record2[1].strip())
    return selected_pairs, scanned_pairs, read_lengths


def sam_tag(fields: list[str], tag: str) -> int | None:
    prefix = f"{tag}:i:"
    for field in fields[11:]:
        if field.startswith(prefix):
            try:
                return int(field[len(prefix):])
            except ValueError:
                return None
    return None


def cigar_indel_bases(cigar: str) -> int:
    total = 0
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        length = int(number) if number else 0
        if char in {"I", "D"}:
            total += length
        number = ""
    return total


def alignment_residual(nm: int, cigar: str, indel_weight: float) -> float:
    indel_bases = cigar_indel_bases(cigar)
    mismatch_bases = max(0, nm - indel_bases)
    return float(mismatch_bases) + max(0.0, indel_weight) * float(indel_bases)


def collect_bwa_alignment_evidence(
    fq1: Path,
    fq2: Path,
    representatives: Mapping[str, str],
    max_pairs: int,
    max_evidence: int,
    top_per_pair: int,
    top_delta: float,
    score_scale: float,
    min_score_per_base: float,
    threads: int,
    bwa_path: str,
    all_alignments: bool,
    enrich_kmer_owners: Mapping[str, tuple[str, ...]] | None = None,
    enrich_k: int = 31,
    min_enrich_kmers: float = 1.0,
    enrich_scan_pairs: int = 0,
    residual_indel_weight: float = 2.0,
    collapse_2field: bool = False,
) -> tuple[int, int, list[dict[str, float]], list[dict[str, float]], Counter[str]]:
    pair_limit = max_pairs if max_pairs > 0 else 2000
    with tempfile.TemporaryDirectory(prefix="raw_read_bwa_") as tmp:
        tmpdir = Path(tmp)
        allele_fa = tmpdir / "alleles.fa"
        sub1 = tmpdir / "subset.R1.fq"
        sub2 = tmpdir / "subset.R2.fq"
        write_allele_fasta(representatives, allele_fa)
        total_pairs, scanned_pairs, read_lengths = write_fastq_subset(
            fq1,
            fq2,
            sub1,
            sub2,
            pair_limit,
            enrich_kmer_owners,
            enrich_k,
            min_enrich_kmers,
            enrich_scan_pairs,
        )
        subprocess.run([bwa_path, "index", str(allele_fa)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        command = [bwa_path, "mem", "-T", "0", "-t", str(max(1, threads))]
        if all_alignments:
            command.append("-a")
        command.extend([str(allele_fa), str(sub1), str(sub2)])
        proc = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    read_metrics: dict[tuple[str, int], dict[str, tuple[int, float]]] = defaultdict(dict)
    for line in proc.stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & 4:
            continue
        allele = fields[2]
        if allele == "*" or allele not in representatives:
            continue
        mate = 1 if flag & 64 else 2 if flag & 128 else 0
        score = sam_tag(fields, "AS")
        nm = sam_tag(fields, "NM")
        if score is None or nm is None:
            continue
        residual = alignment_residual(nm, fields[5], residual_indel_weight)
        key = (fields[0], mate)
        current = read_metrics[key].get(allele)
        if current is None or score > current[0] or (score == current[0] and residual < current[1]):
            read_metrics[key][allele] = (score, residual)

    pair_scores: dict[str, Counter[str]] = defaultdict(Counter)
    pair_residuals: dict[str, Counter[str]] = defaultdict(Counter)
    for (qname, _), allele_scores in read_metrics.items():
        for allele, (score, residual) in allele_scores.items():
            pair_scores[qname][allele] += score
            pair_residuals[qname][allele] += residual

    evidence: list[dict[str, float]] = []
    residual_evidence: list[dict[str, float]] = []
    support: Counter[str] = Counter()
    for qname, scores in pair_scores.items():
        if not scores:
            continue
        best = max(scores.values())
        min_score = min_score_per_base * max(1, read_lengths.get(qname, 0))
        if best < min_score:
            continue
        kept: dict[str, float] = {}
        kept_residuals: dict[str, float] = {}
        for allele, score in scores.most_common(top_per_pair):
            if best - score > top_delta:
                continue
            evidence_allele = allele_2field(allele) if collapse_2field else allele
            if not evidence_allele:
                continue
            loglike = (score - best) / max(score_scale, 1e-6)
            if evidence_allele in kept:
                kept[evidence_allele] = logsumexp([kept[evidence_allele], loglike])
                kept_residuals[evidence_allele] = min(kept_residuals[evidence_allele], float(pair_residuals[qname].get(allele, 0.0)))
            else:
                kept[evidence_allele] = loglike
                kept_residuals[evidence_allele] = float(pair_residuals[qname].get(allele, 0.0))
        if not kept:
            continue
        best_loglike = max(kept.values())
        weights = {allele: math.exp(value - best_loglike) for allele, value in kept.items()}
        total_weight = sum(weights.values()) or 1.0
        for allele, weight in weights.items():
            support[allele] += weight / total_weight
        evidence.append(kept)
        residual_evidence.append(kept_residuals)
        if max_evidence > 0 and len(evidence) >= max_evidence:
            break
    return total_pairs, scanned_pairs, evidence, residual_evidence, support


def collect_mappy_alignment_evidence(
    fq1: Path,
    fq2: Path,
    representatives: Mapping[str, str],
    max_pairs: int,
    max_evidence: int,
    top_per_pair: int,
    top_delta: float,
    score_scale: float,
    min_score_per_base: float,
) -> tuple[int, list[dict[str, float]], Counter[str]]:
    aligner = MappyReadPairAligner(representatives, top_per_pair)
    total_pairs = 0
    evidence: list[dict[str, float]] = []
    support: Counter[str] = Counter()
    try:
        for seq1, seq2 in iter_fastq_pairs(fq1, fq2):
            total_pairs += 1
            if max_pairs > 0 and total_pairs > max_pairs:
                break
            scores = aligner.score_pair(seq1, seq2)
            if not scores:
                continue
            ranked = sorted(scores.items(), key=lambda item: -item[1])
            best_score = ranked[0][1]
            min_score = min_score_per_base * max(1, len(seq1) + len(seq2))
            if best_score < min_score:
                continue
            loglikes: dict[str, float] = {}
            for allele, score in ranked[:top_per_pair]:
                if best_score - score > top_delta:
                    continue
                loglikes[allele] = (score - best_score) / max(score_scale, 1e-6)
            if not loglikes:
                continue
            best_loglike = max(loglikes.values())
            weights = {allele: math.exp(value - best_loglike) for allele, value in loglikes.items()}
            total_weight = sum(weights.values()) or 1.0
            for allele, weight in weights.items():
                support[allele] += weight / total_weight
            evidence.append(loglikes)
            if max_evidence > 0 and len(evidence) >= max_evidence:
                break
    finally:
        aligner.close()
    return total_pairs, evidence, support


def em_fit_weights(
    evidence: list[dict[str, float]],
    candidates: list[str],
    iterations: int,
    min_weight: float,
) -> tuple[dict[str, float], float]:
    if not evidence or not candidates:
        return {}, -math.inf
    weights = {allele: 1.0 / len(candidates) for allele in candidates}
    candidate_set = set(candidates)
    total_loglik = -math.inf
    for _ in range(max(1, iterations)):
        counts = {allele: min_weight for allele in candidates}
        total_loglik = 0.0
        for loglikes in evidence:
            values = []
            alleles = []
            for allele in candidate_set & loglikes.keys():
                values.append(math.log(max(weights[allele], 1e-12)) + loglikes[allele])
                alleles.append(allele)
            if not values:
                continue
            denom = logsumexp(values)
            total_loglik += denom
            for allele, value in zip(alleles, values):
                counts[allele] += math.exp(value - denom)
        total = sum(counts.values()) or 1.0
        weights = {allele: value / total for allele, value in counts.items()}
    return weights, total_loglik


def adjusted_candidate_support(
    support: Mapping[str, float],
    record_counts: Mapping[str, int] | None,
    count_prior_weight: float,
) -> dict[str, float]:
    if not record_counts or count_prior_weight <= 0:
        return {allele: float(value) for allele, value in support.items()}
    adjusted: dict[str, float] = {}
    for allele, value in support.items():
        prior = 1.0 + count_prior_weight * math.log1p(max(0, int(record_counts.get(allele, 1))))
        adjusted[allele] = float(value) * prior
    return adjusted


def alignment_candidates(
    evidence: list[dict[str, float]],
    support: Counter[str],
    max_candidates: int,
    em_iterations: int,
    em_min_weight: float,
    ranker: str = "em",
    record_counts: Mapping[str, int] | None = None,
    count_prior_weight: float = 0.0,
) -> tuple[list[str], dict[str, float], float]:
    adjusted_support = adjusted_candidate_support(support, record_counts, count_prior_weight)
    ranked_by_support = sorted(support, key=lambda allele: (-adjusted_support.get(allele, 0.0), -support.get(allele, 0.0), allele))
    initial = ranked_by_support[:max(max_candidates * 4, max_candidates)]
    if not initial:
        return [], {}, -math.inf
    if ranker in {"support", "support-count-prior"}:
        candidates = initial[:max_candidates]
        final_weights, final_loglik = em_fit_weights(evidence, candidates, em_iterations, em_min_weight)
        ranked = sorted(candidates, key=lambda allele: (-adjusted_support.get(allele, 0.0), -support.get(allele, 0.0), -final_weights.get(allele, 0.0), allele))
        return ranked, final_weights, final_loglik
    em_weights, _ = em_fit_weights(evidence, initial, em_iterations, em_min_weight)
    ranked = sorted(initial, key=lambda allele: (-em_weights.get(allele, 0.0), -support.get(allele, 0.0), allele))
    candidates = ranked[:max_candidates]
    final_weights, final_loglik = em_fit_weights(evidence, candidates, em_iterations, em_min_weight)
    ranked = sorted(candidates, key=lambda allele: (-final_weights.get(allele, 0.0), -support.get(allele, 0.0), allele))
    return ranked, final_weights, final_loglik


def quartet_pair_alignment_loglik(loglikes: dict[str, float], quartet: list[str], chi_r: float, missing_penalty: float) -> float:
    weights = effective_weights(quartet, chi_r)
    if loglikes:
        floor = max(loglikes.values()) - missing_penalty
    else:
        floor = -missing_penalty
    values = []
    for allele, weight in weights.items():
        values.append(math.log(max(weight, 1e-12)) + loglikes.get(allele, floor))
    return logsumexp(values)


def multiset_pair_alignment_loglik(loglikes: dict[str, float], multiset: list[str], missing_penalty: float) -> float:
    if loglikes:
        floor = max(loglikes.values()) - missing_penalty
    else:
        floor = -missing_penalty
    values = []
    for allele, copies in Counter(multiset).items():
        weight = copies / max(1, len(multiset))
        values.append(math.log(max(weight, 1e-12)) + loglikes.get(allele, floor))
    return logsumexp(values)


def score_multiset_alignment_evidence(evidence: list[dict[str, float]], multiset: list[str], missing_penalty: float) -> float:
    if len(multiset) != 4:
        return -math.inf
    return sum(multiset_pair_alignment_loglik(loglikes, multiset, missing_penalty) for loglikes in evidence)


def multiset_search_alignment(
    candidates: list[str],
    evidence: list[dict[str, float]],
    missing_penalty: float,
    max_multisets: int,
) -> tuple[list[str], float, int]:
    best_multiset: list[str] = []
    best_score = -math.inf
    evaluated = 0
    for multiset_tuple in combinations_with_replacement(candidates, 4):
        evaluated += 1
        if max_multisets > 0 and evaluated > max_multisets:
            return best_multiset, best_score, evaluated - 1
        multiset = list(multiset_tuple)
        score = score_multiset_alignment_evidence(evidence, multiset, missing_penalty)
        if score > best_score:
            best_score = score
            best_multiset = multiset
    return best_multiset, best_score, evaluated


def quartet_search_alignment(
    candidates: list[str],
    evidence: list[dict[str, float]],
    chi_r: float,
    missing_penalty: float,
    max_quartets: int,
) -> tuple[list[str], float, int]:
    diploids = list(combinations_with_replacement(candidates, 2))
    best_quartet: list[str] = []
    best_score = -math.inf
    evaluated = 0
    for r_diploid in diploids:
        for d_diploid in diploids:
            evaluated += 1
            if max_quartets > 0 and evaluated > max_quartets:
                return best_quartet, best_score, evaluated - 1
            quartet = [r_diploid[0], r_diploid[1], d_diploid[0], d_diploid[1]]
            score = 0.0
            for loglikes in evidence:
                score += quartet_pair_alignment_loglik(loglikes, quartet, chi_r, missing_penalty)
            if score > best_score:
                best_score = score
                best_quartet = quartet
    return best_quartet, best_score, evaluated


def quartet_pair_residual(
    residuals: Mapping[str, float],
    quartet: list[str],
    chi_r: float,
    missing_penalty: float,
    model: str,
    scale: float,
) -> float:
    if not residuals:
        return missing_penalty
    floor = max(float(value) for value in residuals.values()) + missing_penalty
    if model == "weighted":
        values = []
        for allele, weight in effective_weights(quartet, chi_r).items():
            residual = float(residuals.get(allele, floor))
            values.append(math.log(max(weight, 1e-12)) - residual / max(scale, 1e-6))
        return -max(scale, 1e-6) * logsumexp(values)
    present = [float(residuals[allele]) for allele in set(quartet) if allele in residuals]
    return min(present) if present else floor


def score_quartet_residual_evidence(
    residual_evidence: list[dict[str, float]],
    quartet: list[str],
    chi_r: float,
    missing_penalty: float,
    model: str,
    scale: float,
) -> float:
    if len(quartet) != 4 or not residual_evidence:
        return math.inf
    return sum(quartet_pair_residual(residuals, quartet, chi_r, missing_penalty, model, scale) for residuals in residual_evidence)


def posterior_abundance(
    evidence: list[dict[str, float]],
    alleles: set[str],
    allele_to_group: Mapping[str, str] | None = None,
) -> dict[str, float]:
    counts: Counter[str] = Counter()
    total = 0.0
    for loglikes in evidence:
        kept = {allele: value for allele, value in loglikes.items() if allele in alleles}
        if not kept:
            continue
        best = max(kept.values())
        weights = {allele: math.exp(value - best) for allele, value in kept.items()}
        denom = sum(weights.values()) or 1.0
        for allele, weight in weights.items():
            key = allele_to_group.get(allele, allele) if allele_to_group else allele
            counts[key] += weight / denom
            total += weight / denom
    if total <= 0:
        return {}
    return {allele: value / total for allele, value in counts.items()}


def collapsed_expected_weights(quartet: list[str], chi_r: float, allele_to_group: Mapping[str, str] | None = None) -> dict[str, float]:
    expected: Counter[str] = Counter()
    for allele, weight in effective_weights(quartet, chi_r).items():
        key = allele_to_group.get(allele, allele) if allele_to_group else allele
        expected[key] += weight
    return dict(expected)


def abundance_fit_loglik(observed: Mapping[str, float], expected: Mapping[str, float], evidence_pairs: int, epsilon: float) -> float:
    if not observed or evidence_pairs <= 0:
        return 0.0
    keys = set(observed) | set(expected)
    smoothed_total = 1.0 + epsilon * len(keys)
    score = 0.0
    for key in keys:
        observed_count = float(observed.get(key, 0.0)) * evidence_pairs
        expected_fraction = (float(expected.get(key, 0.0)) + epsilon) / smoothed_total
        score += observed_count * math.log(max(expected_fraction, 1e-12))
    return score


def support_count_abundance(
    support: Mapping[str, float],
    candidates: list[str],
    record_counts: Mapping[str, int] | None,
    count_prior_weight: float,
) -> dict[str, float]:
    adjusted = adjusted_candidate_support(support, record_counts, count_prior_weight)
    total = sum(max(0.0, adjusted.get(allele, 0.0)) for allele in candidates)
    if total <= 0:
        return {}
    return {allele: max(0.0, adjusted.get(allele, 0.0)) / total for allele in candidates}


def score_quartet_support_count(
    quartet: list[str],
    chi_r: float,
    observed_abundance: Mapping[str, float],
    evidence_pairs: int,
    abundance_epsilon: float,
) -> float:
    expected = collapsed_expected_weights(quartet, chi_r)
    return abundance_fit_loglik(observed_abundance, expected, evidence_pairs, abundance_epsilon)


def quartet_search_support_count(
    candidates: list[str],
    support: Mapping[str, float],
    chi_r: float,
    evidence_pairs: int,
    max_quartets: int,
    abundance_epsilon: float,
    record_counts: Mapping[str, int] | None,
    count_prior_weight: float,
    top_quartets: int = 1,
) -> tuple[list[tuple[list[str], float]], int, dict[str, float]]:
    observed = support_count_abundance(support, candidates, record_counts, count_prior_weight)
    diploids = list(combinations_with_replacement(candidates, 2))
    ranked: list[tuple[list[str], float]] = []
    evaluated = 0
    for r_diploid in diploids:
        for d_diploid in diploids:
            evaluated += 1
            if max_quartets > 0 and evaluated > max_quartets:
                ranked.sort(key=lambda item: item[1], reverse=True)
                return ranked[:max(1, top_quartets)], evaluated - 1, observed
            quartet = [r_diploid[0], r_diploid[1], d_diploid[0], d_diploid[1]]
            score = score_quartet_support_count(quartet, chi_r, observed, evidence_pairs, abundance_epsilon)
            ranked.append((quartet, score))
    ranked.sort(key=lambda item: item[1], reverse=True)
    return ranked[:max(1, top_quartets)], evaluated, observed


def score_quartet_alignment_ratio_count(
    evidence: list[dict[str, float]],
    quartet: list[str],
    chi_r: float,
    missing_penalty: float,
    observed_abundance: Mapping[str, float],
    allele_to_group: Mapping[str, str] | None,
    count_weight: float,
    abundance_epsilon: float,
) -> float:
    read_score = score_quartet_alignment_evidence(evidence, quartet, chi_r, missing_penalty)
    expected = collapsed_expected_weights(quartet, chi_r, allele_to_group)
    count_score = abundance_fit_loglik(observed_abundance, expected, len(evidence), abundance_epsilon)
    return read_score + count_weight * count_score


def quartet_search_alignment_ratio_count(
    candidates: list[str],
    evidence: list[dict[str, float]],
    chi_r: float,
    missing_penalty: float,
    max_quartets: int,
    allele_to_group: Mapping[str, str] | None,
    count_weight: float,
    abundance_epsilon: float,
) -> tuple[list[str], float, int, dict[str, float]]:
    observed = posterior_abundance(evidence, set(candidates), allele_to_group)
    diploids = list(combinations_with_replacement(candidates, 2))
    best_quartet: list[str] = []
    best_score = -math.inf
    evaluated = 0
    for r_diploid in diploids:
        for d_diploid in diploids:
            evaluated += 1
            if max_quartets > 0 and evaluated > max_quartets:
                return best_quartet, best_score, evaluated - 1, observed
            quartet = [r_diploid[0], r_diploid[1], d_diploid[0], d_diploid[1]]
            score = score_quartet_alignment_ratio_count(
                evidence,
                quartet,
                chi_r,
                missing_penalty,
                observed,
                allele_to_group,
                count_weight,
                abundance_epsilon,
            )
            if score > best_score:
                best_score = score
                best_quartet = quartet
    return best_quartet, best_score, evaluated, observed


def ranked_quartets_alignment_ratio_count(
    candidates: list[str],
    evidence: list[dict[str, float]],
    chi_r: float,
    missing_penalty: float,
    max_quartets: int,
    allele_to_group: Mapping[str, str] | None,
    count_weight: float,
    abundance_epsilon: float,
    top_quartets: int,
) -> tuple[list[tuple[list[str], float]], int, dict[str, float]]:
    observed = posterior_abundance(evidence, set(candidates), allele_to_group)
    diploids = list(combinations_with_replacement(candidates, 2))
    ranked: list[tuple[list[str], float]] = []
    evaluated = 0
    for r_diploid in diploids:
        for d_diploid in diploids:
            evaluated += 1
            if max_quartets > 0 and evaluated > max_quartets:
                ranked.sort(key=lambda item: item[1], reverse=True)
                return ranked[:max(1, top_quartets)], evaluated - 1, observed
            quartet = [r_diploid[0], r_diploid[1], d_diploid[0], d_diploid[1]]
            score = score_quartet_alignment_ratio_count(
                evidence,
                quartet,
                chi_r,
                missing_penalty,
                observed,
                allele_to_group,
                count_weight,
                abundance_epsilon,
            )
            ranked.append((quartet, score))
    ranked.sort(key=lambda item: item[1], reverse=True)
    return ranked[:max(1, top_quartets)], evaluated, observed


def baseline_preferred_group_relabel(
    quartet: list[str],
    baseline: list[str],
    allele_to_group: Mapping[str, str],
) -> tuple[list[str], str]:
    if len(quartet) != 4 or len(baseline) != 4 or not allele_to_group:
        return quartet, ""
    relabeled = list(quartet)
    notes = []
    groups = sorted(set(collapse_quartet_groups(quartet, allele_to_group)) | set(collapse_quartet_groups(baseline, allele_to_group)))
    for group in groups:
        baseline_members = [allele for allele in baseline if allele_to_group.get(allele, allele) == group]
        quartet_members = [allele for allele in relabeled if allele_to_group.get(allele, allele) == group]
        missing = Counter(baseline_members)
        missing.subtract(Counter(quartet_members))
        replacement_slots = [
            index for index, allele in enumerate(relabeled)
            if allele_to_group.get(allele, allele) == group and Counter(baseline_members)[allele] < Counter(quartet_members)[allele]
        ]
        for allele, count in sorted((item, value) for item, value in missing.items() if value > 0):
            for _ in range(count):
                if not replacement_slots:
                    break
                index = replacement_slots.pop(0)
                old = relabeled[index]
                relabeled[index] = allele
                notes.append(f"{old}->{allele}")
    return relabeled, ";".join(notes)


def low_observed_dropped_group_debug(
    baseline: list[str],
    quartet: list[str],
    allele_to_group: Mapping[str, str],
    observed_abundance: Mapping[str, float],
    min_observed: float,
    max_dropped_copies: int,
) -> str:
    if len(baseline) != 4 or len(quartet) != 4 or not allele_to_group:
        return ""
    baseline_counts = Counter(collapse_quartet_groups(baseline, allele_to_group))
    quartet_counts = Counter(collapse_quartet_groups(quartet, allele_to_group))
    weak = []
    for group, baseline_count in sorted(baseline_counts.items()):
        dropped = baseline_count - quartet_counts.get(group, 0)
        if dropped <= 0:
            continue
        if max_dropped_copies >= 0 and dropped > max_dropped_copies:
            weak.append(f"{group}:dropped={dropped},max={max_dropped_copies}")
            continue
        observed = float(observed_abundance.get(group, 0.0))
        if observed < min_observed:
            weak.append(f"{group}:dropped={dropped},observed={observed:.4f}")
    return "low_observed_dropped_group=" + "|".join(weak) if weak else ""


def unique_quartet_permutations(quartet: list[str]) -> list[list[str]]:
    if len(quartet) != 4:
        return []
    return [list(items) for items in sorted(set(permutations(quartet, 4)))]


def score_fixed_quartets_alignment(
    quartets: list[list[str]],
    evidence: list[dict[str, float]],
    chi_r: float,
    missing_penalty: float,
) -> tuple[list[str], float, int]:
    best_quartet: list[str] = []
    best_score = -math.inf
    evaluated = 0
    for quartet in quartets:
        evaluated += 1
        score = score_quartet_alignment_evidence(evidence, quartet, chi_r, missing_penalty)
        if score > best_score:
            best_score = score
            best_quartet = quartet
    return best_quartet, best_score, evaluated


def quartet_search_multiset_likelihood(
    candidates: list[str],
    evidence: list[dict[str, float]],
    chi_r: float,
    missing_penalty: float,
    max_multisets: int,
) -> tuple[list[str], float, int, float]:
    multiset, multiset_score, evaluated = multiset_search_alignment(candidates, evidence, missing_penalty, max_multisets)
    quartet, quartet_score, _ = score_fixed_quartets_alignment(unique_quartet_permutations(multiset), evidence, chi_r, missing_penalty)
    return quartet, quartet_score, evaluated, multiset_score


def allele_pair_likelihood_delta(
    evidence: list[dict[str, float]],
    allele_a: str,
    allele_b: str,
    missing_penalty: float,
) -> tuple[int, int, float, list[float]]:
    observed = 0
    shared = 0
    total_delta = 0.0
    deltas: list[float] = []
    for loglikes in evidence:
        has_a = allele_a in loglikes
        has_b = allele_b in loglikes
        if not has_a and not has_b:
            continue
        observed += 1
        if has_a and has_b:
            shared += 1
        floor = (max(loglikes.values()) if loglikes else 0.0) - missing_penalty
        delta = abs(loglikes.get(allele_a, floor) - loglikes.get(allele_b, floor))
        total_delta += delta
        deltas.append(delta)
    mean_delta = total_delta / observed if observed else math.inf
    return observed, shared, mean_delta, deltas


def build_ambiguity_groups(
    alleles: list[str],
    evidence: list[dict[str, float]],
    missing_penalty: float,
    min_observed_pairs: int,
    min_shared_fraction: float,
    max_mean_delta: float,
    max_pair_delta: float,
    min_within_pair_delta_fraction: float,
) -> dict[str, str]:
    unique_alleles = sorted(set(alleles))
    parent = {allele: allele for allele in unique_alleles}

    def find(allele: str) -> str:
        while parent[allele] != allele:
            parent[allele] = parent[parent[allele]]
            allele = parent[allele]
        return allele

    def union(left: str, right: str) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[max(root_left, root_right)] = min(root_left, root_right)

    for index, allele_a in enumerate(unique_alleles):
        for allele_b in unique_alleles[index + 1:]:
            observed, shared, mean_delta, deltas = allele_pair_likelihood_delta(
                evidence,
                allele_a,
                allele_b,
                missing_penalty,
            )
            if observed < min_observed_pairs:
                continue
            shared_fraction = shared / observed if observed else 0.0
            if shared_fraction < min_shared_fraction:
                continue
            within_fraction = sum(1 for delta in deltas if delta <= max_pair_delta) / len(deltas) if deltas else 0.0
            if mean_delta <= max_mean_delta and within_fraction >= min_within_pair_delta_fraction:
                union(allele_a, allele_b)

    groups: dict[str, list[str]] = defaultdict(list)
    for allele in unique_alleles:
        groups[find(allele)].append(allele)
    allele_to_group: dict[str, str] = {}
    for members in groups.values():
        group_name = "/".join(sorted(members))
        for allele in members:
            allele_to_group[allele] = group_name
    return allele_to_group


def build_adaptive_ambiguity_groups(
    alleles: list[str],
    evidence: list[dict[str, float]],
    missing_penalty: float,
    min_observed_pairs: int,
    min_shared_fraction: float,
    max_mean_delta_fraction: float,
    max_pair_delta_fraction: float,
    min_within_pair_delta_fraction: float,
    min_mean_gap: float,
) -> tuple[dict[str, str], str]:
    unique_alleles = sorted(set(alleles))
    parent = {allele: allele for allele in unique_alleles}

    def find(allele: str) -> str:
        while parent[allele] != allele:
            parent[allele] = parent[parent[allele]]
            allele = parent[allele]
        return allele

    def union(left: str, right: str) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[max(root_left, root_right)] = min(root_left, root_right)

    max_mean_delta = max(0.0, missing_penalty * max_mean_delta_fraction)
    max_pair_delta = max(0.0, missing_penalty * max_pair_delta_fraction)
    pair_rows: list[dict[str, object]] = []
    for index, allele_a in enumerate(unique_alleles):
        for allele_b in unique_alleles[index + 1:]:
            observed, shared, mean_delta, deltas = allele_pair_likelihood_delta(
                evidence,
                allele_a,
                allele_b,
                missing_penalty,
            )
            if observed < min_observed_pairs or not math.isfinite(mean_delta):
                continue
            shared_fraction = shared / observed if observed else 0.0
            if shared_fraction < min_shared_fraction:
                continue
            shared_deltas = []
            for loglikes in evidence:
                if allele_a in loglikes and allele_b in loglikes:
                    shared_deltas.append(abs(loglikes[allele_a] - loglikes[allele_b]))
            within_fraction = sum(1 for delta in shared_deltas if delta <= max_pair_delta) / len(shared_deltas) if shared_deltas else 0.0
            if within_fraction < min_within_pair_delta_fraction:
                continue
            pair_rows.append({
                "allele_a": allele_a,
                "allele_b": allele_b,
                "mean_delta": mean_delta,
                "shared_fraction": shared_fraction,
                "within_fraction": within_fraction,
            })

    sorted_means = sorted(float(row["mean_delta"]) for row in pair_rows if float(row["mean_delta"]) <= max_mean_delta)
    adaptive_mean_delta = max_mean_delta if sorted_means else -math.inf
    if sorted_means and min_mean_gap > 0:
        adaptive_mean_delta = sorted_means[0]
        for left, right in zip(sorted_means, sorted_means[1:]):
            if right - left >= min_mean_gap:
                adaptive_mean_delta = left
                break
            adaptive_mean_delta = right
        adaptive_mean_delta = min(adaptive_mean_delta, max_mean_delta)

    collapsed_pairs = 0
    for row in pair_rows:
        if float(row["mean_delta"]) <= adaptive_mean_delta:
            union(str(row["allele_a"]), str(row["allele_b"]))
            collapsed_pairs += 1

    groups: dict[str, list[str]] = defaultdict(list)
    for allele in unique_alleles:
        groups[find(allele)].append(allele)
    allele_to_group: dict[str, str] = {}
    for members in groups.values():
        group_name = "/".join(sorted(members))
        for allele in members:
            allele_to_group[allele] = group_name
    debug = (
        f"adaptive_pairs={len(pair_rows)},collapsed_pairs={collapsed_pairs},"
        f"mean_delta={adaptive_mean_delta:.4f},max_mean_delta={max_mean_delta:.4f},"
        f"max_pair_delta={max_pair_delta:.4f}"
    )
    return allele_to_group, debug


def collapse_quartet_groups(quartet: list[str], allele_to_group: Mapping[str, str]) -> list[str]:
    return [allele_to_group.get(allele, allele) for allele in quartet]


def same_collapsed_side_counts(baseline: list[str], quartet: list[str], allele_to_group: Mapping[str, str]) -> bool:
    baseline_groups = collapse_quartet_groups(baseline, allele_to_group)
    quartet_groups = collapse_quartet_groups(quartet, allele_to_group)
    return Counter(baseline_groups[:2]) == Counter(quartet_groups[:2]) and Counter(baseline_groups[2:]) == Counter(quartet_groups[2:])


def same_collapsed_multiset(baseline: list[str], quartet: list[str], allele_to_group: Mapping[str, str]) -> bool:
    return Counter(collapse_quartet_groups(baseline, allele_to_group)) == Counter(collapse_quartet_groups(quartet, allele_to_group))


def serialize_ambiguity_groups(allele_to_group: Mapping[str, str], alleles: list[str]) -> str:
    groups: dict[str, list[str]] = defaultdict(list)
    for allele in sorted(set(alleles)):
        groups[allele_to_group.get(allele, allele)].append(allele)
    return ";".join(
        "/".join(sorted(members))
        for members in sorted(groups.values(), key=lambda members: (members[0], len(members)))
    )


def is_balanced_het_to_hom_side_flip(baseline: list[str], quartet: list[str]) -> bool:
    if len(baseline) != 4 or len(quartet) != 4:
        return False
    baseline_r = Counter(baseline[:2])
    baseline_d = Counter(baseline[2:])
    if baseline_r != baseline_d or len(baseline_r) != 2:
        return False
    direct_r = Counter(quartet[:2])
    direct_d = Counter(quartet[2:])
    return len(direct_r) == 1 and len(direct_d) == 1 and direct_r != direct_d


def replacement_gate_debug(
    baseline: list[str],
    quartet: list[str],
    baseline_loglik: float,
    quartet_loglik: float,
    evidence_pairs: int,
    support: Mapping[str, float],
    em_weights: Mapping[str, float],
    allele_to_group: Mapping[str, str],
    min_loglik_per_pair: float,
    min_support: float,
    min_em_weight: float,
    max_new_copies: int,
    reject_same_ambiguity_multiset: bool,
    reject_balanced_group_dosage_change: bool,
) -> tuple[bool, str]:
    if len(baseline) != 4 or len(quartet) != 4:
        return False, "invalid_quartet"
    if Counter(quartet) == Counter(baseline):
        return False, "same_multiset"
    if not math.isfinite(baseline_loglik) or not math.isfinite(quartet_loglik):
        return False, "nonfinite_loglik"
    if reject_same_ambiguity_multiset and allele_to_group and same_collapsed_multiset(baseline, quartet, allele_to_group):
        return False, "same_ambiguity_multiset"
    if reject_balanced_group_dosage_change and allele_to_group:
        baseline_groups = collapse_quartet_groups(baseline, allele_to_group)
        quartet_groups = collapse_quartet_groups(quartet, allele_to_group)
        baseline_r = Counter(baseline_groups[:2])
        baseline_d = Counter(baseline_groups[2:])
        if baseline_r == baseline_d and len(baseline_r) >= 2 and Counter(baseline_groups) != Counter(quartet_groups):
            return False, "balanced_group_dosage_change"

    new_copies = Counter(quartet)
    new_copies.subtract(Counter(baseline))
    new_copies = Counter({allele: count for allele, count in new_copies.items() if count > 0})
    total_new_copies = sum(new_copies.values())
    if max_new_copies >= 0 and total_new_copies > max_new_copies:
        return False, f"too_many_new_copies={total_new_copies}"

    min_total_delta = min_loglik_per_pair * max(1, evidence_pairs)
    loglik_delta = quartet_loglik - baseline_loglik
    if loglik_delta < min_total_delta:
        return False, f"low_delta={loglik_delta:.4f}<min={min_total_delta:.4f}"

    weak = []
    for allele in sorted(new_copies):
        allele_support = float(support.get(allele, 0.0))
        allele_em = float(em_weights.get(allele, 0.0))
        if allele_support < min_support or allele_em < min_em_weight:
            weak.append(f"{allele}:support={allele_support:.2f},em={allele_em:.4f}")
    if weak:
        return False, "weak_new_allele=" + "|".join(weak)

    new_desc = ",".join(f"{allele}x{count}" for allele, count in sorted(new_copies.items()))
    return True, f"pass_delta={loglik_delta:.4f},min={min_total_delta:.4f},new={new_desc}"


def score_quartet_alignment_evidence(
    evidence: list[dict[str, float]],
    quartet: list[str],
    chi_r: float,
    missing_penalty: float,
) -> float:
    if len(quartet) != 4:
        return -math.inf
    return sum(quartet_pair_alignment_loglik(loglikes, quartet, chi_r, missing_penalty) for loglikes in evidence)


def top_candidates(support: Counter[str], max_candidates: int, min_support: float) -> list[str]:
    ranked = [allele for allele, value in support.most_common() if value >= min_support]
    return ranked[:max_candidates]


def quartet_search(
    candidates: list[str],
    evidence: list[Counter[str]],
    chi_r: float,
    score_scale: float,
    max_quartets: int,
) -> tuple[list[str], float, int]:
    diploids = list(combinations_with_replacement(candidates, 2))
    best_quartet: list[str] = []
    best_score = -math.inf
    evaluated = 0
    for r_diploid in diploids:
        for d_diploid in diploids:
            evaluated += 1
            if max_quartets > 0 and evaluated > max_quartets:
                return best_quartet, best_score, evaluated - 1
            quartet = [r_diploid[0], r_diploid[1], d_diploid[0], d_diploid[1]]
            score = 0.0
            for allele_counts in evidence:
                score += quartet_pair_loglik(allele_counts, quartet, chi_r, score_scale)
            if score > best_score:
                best_score = score
                best_quartet = quartet
    return best_quartet, best_score, evaluated


def serialize_support(support: Counter[str], candidates: list[str]) -> str:
    return ",".join(f"{allele}:{support[allele]:.2f}" for allele in candidates)


def serialize_weights(weights: Mapping[str, float], candidates: list[str]) -> str:
    return ",".join(f"{allele}:{weights.get(allele, 0.0):.4f}" for allele in candidates)

def serialize_float_mapping(values: Mapping[str, float]) -> str:
    return ",".join(f"{key}:{values[key]:.4f}" for key in sorted(values))


def load_independent_sidecopy_allow(path: Path | None) -> set[tuple[str, str]]:
    if path is None:
        return set()
    allowed: set[tuple[str, str]] = set()
    for row in read_tsv(path):
        sample = row.get("sample", "")
        gene = normalize_gene(row.get("gene", ""))
        if not sample or not gene:
            continue
        allow_value = row.get("allow", row.get("accept", "1")).strip().lower()
        if allow_value in {"1", "true", "yes", "y", "accept", "accepted", "allow", "allowed"}:
            allowed.add((sample, gene))
    return allowed


def quartet_multiset_score(quartet: list[str], truth_r: list[str], truth_d: list[str]) -> int:
    truth_counts = Counter(truth_r + truth_d)
    pred_counts = Counter(quartet)
    return sum(min(truth_counts[allele], pred_counts[allele]) for allele in truth_counts)


def score_row(
    row: dict[str, str],
    args: argparse.Namespace,
    representatives: dict[str, str],
    gene_kmer_owners: dict[str, tuple[str, ...]],
    alignment_seed_owners: dict[str, tuple[str, ...]],
    alignment_representatives: dict[str, str] | None = None,
    allele_record_counts: Mapping[str, int] | None = None,
) -> dict[str, object]:
    sample = row["sample"]
    gene = normalize_gene(row["gene"])
    baseline = split_alleles(row["baseline_pred_R"]) + split_alleles(row["baseline_pred_D"])
    truth_r = split_alleles(row["truth_R"])
    truth_d = split_alleles(row["truth_D"])
    baseline_multiset_score = quartet_multiset_score(baseline, truth_r, truth_d) if len(baseline) == 4 else 0
    result: dict[str, object] = {
        "sample": sample,
        "set": row.get("set", ""),
        "gene": gene,
        "gene_config": getattr(args, "gene_config_overrides", ""),
        "baseline_score": row.get("baseline_score", ""),
        "baseline_multiset_score": baseline_multiset_score,
        "truth_R": row.get("truth_R", ""),
        "truth_D": row.get("truth_D", ""),
        "baseline_pred_R": row.get("baseline_pred_R", ""),
        "baseline_pred_D": row.get("baseline_pred_D", ""),
    }
    if args.gene_action == "baseline":
        baseline_score = int(row.get("baseline_score", 0))
        result.update({
            "direct_score": baseline_score,
            "delta_vs_baseline": 0,
            "direct_multiset_score": baseline_multiset_score,
            "delta_multiset_vs_baseline": 0,
            "status": "same",
            "direct_pred_R": row.get("baseline_pred_R", ""),
            "direct_pred_D": row.get("baseline_pred_D", ""),
            "evidence_backend": args.evidence_backend,
            "alignment_engine": args.alignment_engine if args.evidence_backend == "alignment" else "",
            "gate_decision": "baseline_gene_config",
            "baseline_quartet": ",".join(baseline),
            "error": "",
        })
        return result
    fq1, fq2, read_source = select_fastqs(args, sample, gene)
    result["read_source"] = read_source
    if not fq1.exists() or not fq2.exists():
        result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "missing_fastq", "error": f"{fq1},{fq2}"})
        return result

    chi_r = read_chi_r(args.spechla_root, sample)
    em_weights: dict[str, float] = {}
    em_loglik = -math.inf
    candidate_kmer_owners: dict[str, tuple[str, ...]] = {}
    raw_scan_pairs: int | str = ""
    gate_decision = "direct"
    baseline_loglik = -math.inf
    ambiguity_groups = ""
    baseline_group_quartet = ""
    direct_group_quartet = ""
    ambiguity_collapse_debug = ""
    replacement_allele_to_group: dict[str, str] = {}
    unrestricted_quartet: list[str] = []
    unrestricted_loglik = -math.inf
    unrestricted_residual = math.inf
    ratio_count_observed = ""
    ratio_count_groups = ""
    ratio_count_observed_values: dict[str, float] = {}
    replacement_debug = ""
    selector_attempts_debug = ""
    baseline_residual = math.inf
    best_residual = math.inf
    residual_debug = ""
    if args.evidence_backend == "kmer":
        total_pairs, support = scan_candidate_support(fq1, fq2, gene_kmer_owners, args.k, args.max_candidate_pairs)
        candidates = top_candidates(support, args.max_candidates, args.min_candidate_support)
        if len(candidates) < 1:
            result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "no_candidates", "error": "no supported allele candidates"})
            return result

        candidate_kmer_owners = build_candidate_kmer_owners(representatives, candidates, args.k, args.max_owner_fraction)
        _, evidence = collect_evidence_pairs(
            fq1,
            fq2,
            candidate_kmer_owners,
            args.k,
            args.max_evidence_pairs,
            args.min_pair_informative_kmers,
        )
        if not evidence:
            result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "no_evidence", "error": "no informative read pairs"})
            return result
        quartet, loglik, evaluated = quartet_search(candidates, evidence, chi_r, args.score_scale, args.max_quartets)
        evidence_pairs = len(evidence)
    else:
        try:
            if args.alignment_engine == "bwa":
                enrich_owners = gene_kmer_owners if read_source == "raw" and args.raw_gene_enrich else None
                total_pairs, raw_scan_pairs, alignment_evidence, alignment_residual_evidence, support = collect_bwa_alignment_evidence(
                    fq1,
                    fq2,
                    alignment_representatives or representatives,
                    args.max_candidate_pairs,
                    args.max_evidence_pairs,
                    args.alignment_top_per_pair,
                    args.alignment_top_delta,
                    args.alignment_score_scale,
                    args.alignment_min_score_per_base,
                    args.bwa_threads,
                    args.bwa_path,
                    not args.bwa_primary_only,
                    enrich_owners,
                    args.k,
                    args.raw_gene_enrich_min_kmers,
                    args.raw_gene_enrich_scan_pairs,
                    args.residual_indel_weight,
                    bool(args.alignment_full_allele_collapse),
                )
            else:
                total_pairs, alignment_evidence, support = collect_alignment_evidence(
                    fq1,
                    fq2,
                    representatives,
                    alignment_seed_owners,
                    args.max_candidate_pairs,
                    args.max_evidence_pairs,
                    args.alignment_seed_k,
                    args.alignment_prefilter_alleles,
                    args.alignment_top_per_pair,
                    args.alignment_top_delta,
                    args.alignment_score_scale,
                    args.alignment_min_score_per_base,
                )
                alignment_residual_evidence = []
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "alignment_error", "error": str(exc)})
            return result
        if not alignment_evidence:
            result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "no_evidence", "error": "no passing alignment read pairs"})
            return result
        candidates, em_weights, em_loglik = alignment_candidates(
            alignment_evidence,
            support,
            args.max_candidates,
            args.em_iterations,
            args.em_min_weight,
            args.alignment_candidate_ranker,
            allele_record_counts,
            args.alignment_count_prior_weight,
        )
        if args.side_copy_gate != "off" or args.baseline_multiset_gate:
            for allele in baseline:
                if allele in representatives and allele not in candidates:
                    candidates.append(allele)
        if len(candidates) < 1:
            result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "no_candidates", "error": "no alignment-supported allele candidates"})
            return result
        candidate_set = set(candidates)
        evidence = []
        residual_evidence: list[dict[str, float]] = []
        for index, loglikes in enumerate(alignment_evidence):
            if not (candidate_set & loglikes.keys()):
                continue
            evidence.append({allele: value for allele, value in loglikes.items() if allele in candidate_set})
            if index < len(alignment_residual_evidence):
                residual_evidence.append({allele: value for allele, value in alignment_residual_evidence[index].items() if allele in candidate_set})
            else:
                residual_evidence.append({})
        if not evidence:
            result.update({"direct_score": 0, "delta_vs_baseline": -int(row.get("baseline_score", 0)), "direct_multiset_score": 0, "delta_multiset_vs_baseline": -baseline_multiset_score, "status": "no_evidence", "error": "no candidate-overlapping alignment read pairs"})
            return result
        selector_allele_to_group: dict[str, str] = {}
        if args.quartet_selector in {"ratio-count", "constrained-ratio"} and args.ratio_count_use_ambiguity_groups:
            selector_alleles = sorted(set(candidates) | set(baseline))
            selector_allele_to_group, selector_debug = build_adaptive_ambiguity_groups(
                selector_alleles,
                evidence,
                args.alignment_missing_penalty,
                args.ambiguity_collapse_min_pairs,
                args.ambiguity_collapse_min_shared_fraction,
                args.ambiguity_collapse_adaptive_max_mean_fraction,
                args.ambiguity_collapse_adaptive_max_pair_fraction,
                args.ambiguity_collapse_min_within_pair_delta_fraction,
                args.ambiguity_collapse_adaptive_min_mean_gap,
            )
            ratio_count_groups = serialize_ambiguity_groups(selector_allele_to_group, selector_alleles)
            ambiguity_collapse_debug = selector_debug
        selector_ranked_quartets: list[tuple[list[str], float]] = []
        if args.quartet_selector in {"ratio-count", "constrained-ratio"}:
            if args.selector_iterations > 1:
                selector_ranked_quartets, evaluated, observed_abundance = ranked_quartets_alignment_ratio_count(
                    candidates,
                    evidence,
                    chi_r,
                    args.alignment_missing_penalty,
                    args.max_quartets,
                    selector_allele_to_group if selector_allele_to_group else None,
                    args.ratio_count_weight,
                    args.ratio_count_epsilon,
                    args.selector_iterations,
                )
                quartet, loglik = selector_ranked_quartets[0] if selector_ranked_quartets else ([], -math.inf)
            else:
                quartet, loglik, evaluated, observed_abundance = quartet_search_alignment_ratio_count(
                    candidates,
                    evidence,
                    chi_r,
                    args.alignment_missing_penalty,
                    args.max_quartets,
                    selector_allele_to_group if selector_allele_to_group else None,
                    args.ratio_count_weight,
                    args.ratio_count_epsilon,
                )
                selector_ranked_quartets = [(list(quartet), loglik)] if quartet else []
            ratio_count_observed_values = dict(observed_abundance)
            ratio_count_observed = serialize_float_mapping(observed_abundance)
        elif args.quartet_selector == "support-count":
            selector_ranked_quartets, evaluated, observed_abundance = quartet_search_support_count(
                candidates,
                support,
                chi_r,
                len(evidence),
                args.max_quartets,
                args.ratio_count_epsilon,
                allele_record_counts,
                args.alignment_count_prior_weight,
                args.selector_iterations,
            )
            quartet, loglik = selector_ranked_quartets[0] if selector_ranked_quartets else ([], -math.inf)
            ratio_count_observed_values = dict(observed_abundance)
            ratio_count_observed = serialize_float_mapping(observed_abundance)
        elif args.quartet_selector == "multiset-likelihood":
            quartet, loglik, evaluated, multiset_loglik = quartet_search_multiset_likelihood(
                candidates,
                evidence,
                chi_r,
                args.alignment_missing_penalty,
                args.max_quartets,
            )
            selector_ranked_quartets = [(list(quartet), loglik)] if quartet else []
            ratio_count_observed = f"multiset_loglik:{multiset_loglik:.4f}"
        else:
            quartet, loglik, evaluated = quartet_search_alignment(candidates, evidence, chi_r, args.alignment_missing_penalty, args.max_quartets)
        baseline_loglik = score_quartet_alignment_evidence(evidence, baseline, chi_r, args.alignment_missing_penalty)
        baseline_residual = score_quartet_residual_evidence(
            residual_evidence, baseline, chi_r, args.residual_missing_penalty, args.residual_model, args.residual_scale
        )
        unrestricted_quartet = list(quartet)
        unrestricted_loglik = loglik
        unrestricted_residual = score_quartet_residual_evidence(
            residual_evidence, unrestricted_quartet, chi_r, args.residual_missing_penalty, args.residual_model, args.residual_scale
        )
        constrained_ratio_debug = ""
        if args.quartet_selector == "constrained-ratio" and selector_allele_to_group:
            relabeled_quartet, relabel_debug = baseline_preferred_group_relabel(unrestricted_quartet, baseline, selector_allele_to_group)
            if relabeled_quartet != unrestricted_quartet:
                unrestricted_quartet = relabeled_quartet
                unrestricted_loglik = score_quartet_alignment_ratio_count(
                    evidence,
                    unrestricted_quartet,
                    chi_r,
                    args.alignment_missing_penalty,
                    ratio_count_observed_values,
                    selector_allele_to_group,
                    args.ratio_count_weight,
                    args.ratio_count_epsilon,
                )
                unrestricted_residual = score_quartet_residual_evidence(
                    residual_evidence, unrestricted_quartet, chi_r, args.residual_missing_penalty, args.residual_model, args.residual_scale
                )
                constrained_ratio_debug = f"baseline_group_relabel={relabel_debug}"
        if args.side_copy_gate == "multiset":
            quartet, loglik, evaluated = score_fixed_quartets_alignment(
                unique_quartet_permutations(baseline),
                evidence,
                chi_r,
                args.alignment_missing_penalty,
            )
            gate_decision = "sidecopy_multiset"
        elif args.side_copy_gate == "no-new-2field":
            allowed = sorted(set(allele for allele in baseline if allele in representatives))
            quartet, loglik, evaluated = quartet_search_alignment(
                allowed,
                evidence,
                chi_r,
                args.alignment_missing_penalty,
                args.max_quartets,
            )
            gate_decision = "sidecopy_no_new_2field"
        if args.baseline_gate_margin is not None and math.isfinite(baseline_loglik):
            if loglik - baseline_loglik < args.baseline_gate_margin:
                quartet = baseline
                loglik = baseline_loglik
                gate_decision = "baseline"
        if args.baseline_confidence_max_loglik_per_pair is not None and math.isfinite(baseline_loglik):
            baseline_loglik_per_pair = baseline_loglik / max(1, len(evidence))
            if baseline_loglik_per_pair >= args.baseline_confidence_max_loglik_per_pair:
                quartet = baseline
                loglik = baseline_loglik
                gate_decision = f"baseline_confidence={baseline_loglik_per_pair:.4f}"
        if args.ambiguity_collapse_gate != "off":
            collapse_alleles = sorted(set(candidates) | set(baseline) | set(quartet))
            if args.ambiguity_collapse_gate.startswith("adaptive-"):
                allele_to_group, ambiguity_collapse_debug = build_adaptive_ambiguity_groups(
                    collapse_alleles,
                    evidence,
                    args.alignment_missing_penalty,
                    args.ambiguity_collapse_min_pairs,
                    args.ambiguity_collapse_min_shared_fraction,
                    args.ambiguity_collapse_adaptive_max_mean_fraction,
                    args.ambiguity_collapse_adaptive_max_pair_fraction,
                    args.ambiguity_collapse_min_within_pair_delta_fraction,
                    args.ambiguity_collapse_adaptive_min_mean_gap,
                )
            else:
                allele_to_group = build_ambiguity_groups(
                    collapse_alleles,
                    evidence,
                    args.alignment_missing_penalty,
                    args.ambiguity_collapse_min_pairs,
                    args.ambiguity_collapse_min_shared_fraction,
                    args.ambiguity_collapse_max_mean_delta,
                    args.ambiguity_collapse_max_pair_delta,
                    args.ambiguity_collapse_min_within_pair_delta_fraction,
                )
            replacement_allele_to_group = dict(allele_to_group)
            ambiguity_groups = serialize_ambiguity_groups(allele_to_group, collapse_alleles)
            baseline_group_quartet = ",".join(collapse_quartet_groups(baseline, allele_to_group))
            direct_group_quartet = ",".join(collapse_quartet_groups(quartet, allele_to_group))
            if args.ambiguity_collapse_gate in {"group-side", "adaptive-side"} and not same_collapsed_side_counts(baseline, quartet, allele_to_group):
                independent_allowed = (sample, gene) in args.independent_sidecopy_allow
                if independent_allowed:
                    gate_decision = f"{gate_decision}_independent_sidecopy"
                else:
                    quartet = baseline
                    loglik = baseline_loglik if math.isfinite(baseline_loglik) else loglik
                    gate_decision = "baseline_ambiguity_group_side"
                    direct_group_quartet = baseline_group_quartet
            elif args.ambiguity_collapse_gate in {"group-multiset", "adaptive-multiset"} and not same_collapsed_multiset(baseline, quartet, allele_to_group):
                quartet = baseline
                loglik = baseline_loglik if math.isfinite(baseline_loglik) else loglik
                gate_decision = "baseline_ambiguity_group_multiset"
                direct_group_quartet = baseline_group_quartet
        if args.replacement_gate == "likelihood":
            replacement_ok = False
            attempt_notes = []
            proposal_source = selector_ranked_quartets if selector_ranked_quartets else [(unrestricted_quartet, unrestricted_loglik)]
            for rank, (candidate_quartet, candidate_loglik) in enumerate(proposal_source[:max(1, args.selector_iterations)], 1):
                proposal_quartet = list(candidate_quartet)
                proposal_loglik = candidate_loglik
                proposal_prefix = ""
                if args.quartet_selector == "constrained-ratio" and replacement_allele_to_group:
                    relabeled_quartet, relabel_debug = baseline_preferred_group_relabel(proposal_quartet, baseline, replacement_allele_to_group)
                    if relabeled_quartet != proposal_quartet:
                        proposal_quartet = relabeled_quartet
                        proposal_loglik = score_quartet_alignment_ratio_count(
                            evidence,
                            proposal_quartet,
                            chi_r,
                            args.alignment_missing_penalty,
                            ratio_count_observed_values,
                            replacement_allele_to_group,
                            args.ratio_count_weight,
                            args.ratio_count_epsilon,
                        )
                        proposal_prefix = f"baseline_group_relabel={relabel_debug}"
                    proposal_residual = score_quartet_residual_evidence(
                        residual_evidence, proposal_quartet, chi_r, args.residual_missing_penalty, args.residual_model, args.residual_scale
                    )
                dropped_group_debug = ""
                if args.quartet_selector == "constrained-ratio" and replacement_allele_to_group:
                    dropped_group_debug = low_observed_dropped_group_debug(
                        baseline,
                        proposal_quartet,
                        replacement_allele_to_group,
                        ratio_count_observed_values,
                        args.constrained_ratio_min_dropped_group_observed,
                        args.constrained_ratio_max_dropped_group_copies,
                    )
                if dropped_group_debug:
                    candidate_ok = False
                    candidate_debug = dropped_group_debug
                else:
                    candidate_ok, candidate_debug = replacement_gate_debug(
                        baseline,
                        proposal_quartet,
                        baseline_loglik,
                        proposal_loglik,
                        len(evidence),
                        support,
                        em_weights,
                        replacement_allele_to_group,
                        args.replacement_min_loglik_per_pair,
                        args.replacement_min_support,
                        args.replacement_min_em_weight,
                        args.replacement_max_new_copies,
                        not args.replacement_allow_same_ambiguity_multiset,
                        not args.replacement_allow_balanced_group_dosage_change,
                    )
                if candidate_ok and args.residual_gate != "off" and math.isfinite(baseline_residual) and math.isfinite(proposal_residual):
                    residual_improvement = baseline_residual - proposal_residual
                    min_residual_improvement = args.residual_min_improvement_per_pair * max(1, len(residual_evidence))
                    if args.residual_gate == "not-worse" and residual_improvement < -min_residual_improvement:
                        candidate_ok = False
                        candidate_debug = f"residual_worse={residual_improvement:.4f}<min={-min_residual_improvement:.4f}"
                    elif args.residual_gate == "improvement" and residual_improvement < min_residual_improvement:
                        candidate_ok = False
                        candidate_debug = f"residual_low_gain={residual_improvement:.4f}<min={min_residual_improvement:.4f}"
                    else:
                        candidate_debug = f"{candidate_debug};residual_gain={residual_improvement:.4f},min={min_residual_improvement:.4f}"
                if proposal_prefix:
                    candidate_debug = f"{proposal_prefix};{candidate_debug}" if candidate_debug else proposal_prefix
                attempt_notes.append(f"rank{rank}:{candidate_debug}")
                if candidate_ok:
                    replacement_ok = True
                    replacement_debug = candidate_debug
                    unrestricted_quartet = proposal_quartet
                    unrestricted_loglik = proposal_loglik
                    unrestricted_residual = proposal_residual
                    break
            selector_attempts_debug = "|".join(attempt_notes[:max(1, args.selector_iterations)])
            if not replacement_debug and attempt_notes:
                replacement_debug = attempt_notes[-1].split(":", 1)[1]
            if replacement_ok:
                quartet = unrestricted_quartet
                loglik = unrestricted_loglik
                gate_decision = "replacement_likelihood"
                if baseline_group_quartet:
                    direct_group_quartet = ""
        if args.reject_balanced_het_to_hom_side and is_balanced_het_to_hom_side_flip(baseline, quartet):
            quartet = baseline
            loglik = baseline_loglik if math.isfinite(baseline_loglik) else loglik
            gate_decision = "baseline_balanced_het_guard"
            if baseline_group_quartet:
                direct_group_quartet = baseline_group_quartet
        if args.baseline_multiset_gate and Counter(quartet) != Counter(baseline):
            quartet = baseline
            loglik = baseline_loglik if math.isfinite(baseline_loglik) else loglik
            gate_decision = "baseline_multiset"
        evidence_pairs = len(evidence)
        best_residual = score_quartet_residual_evidence(
            residual_evidence, quartet, chi_r, args.residual_missing_penalty, args.residual_model, args.residual_scale
        )
        if math.isfinite(baseline_residual) and math.isfinite(best_residual):
            residual_debug = f"gain={baseline_residual - best_residual:.4f}"
    direct_score = score_quartet(quartet, truth_r, truth_d) if len(quartet) == 4 else 0
    direct_multiset_score = quartet_multiset_score(quartet, truth_r, truth_d) if len(quartet) == 4 else 0
    baseline_score = int(row.get("baseline_score", 0))
    delta = direct_score - baseline_score
    delta_multiset = direct_multiset_score - baseline_multiset_score
    status = "improved" if delta > 0 else "regressed" if delta < 0 else "same"
    result.update({
        "direct_score": direct_score,
        "delta_vs_baseline": delta,
        "direct_multiset_score": direct_multiset_score,
        "delta_multiset_vs_baseline": delta_multiset,
        "status": status,
        "direct_pred_R": ",".join(quartet[:2]),
        "direct_pred_D": ",".join(quartet[2:4]),
        "chi_r": f"{chi_r:.4f}",
        "evidence_backend": args.evidence_backend,
        "alignment_engine": args.alignment_engine if args.evidence_backend == "alignment" else "",
        "raw_or_gene_total_pairs": total_pairs,
        "raw_scan_pairs": raw_scan_pairs,
        "evidence_pairs": evidence_pairs,
        "gene_db_alleles": len(representatives),
        "gene_informative_kmers": len(gene_kmer_owners),
        "candidate_informative_kmers": len(candidate_kmer_owners),
        "candidates": ",".join(candidates),
        "candidate_support": serialize_support(support, candidates),
        "em_weights": serialize_weights(em_weights, candidates),
        "em_loglik": f"{em_loglik:.4f}" if math.isfinite(em_loglik) else "",
        "ambiguity_groups": ambiguity_groups,
        "baseline_group_quartet": baseline_group_quartet,
        "direct_group_quartet": direct_group_quartet,
        "ambiguity_collapse_debug": ambiguity_collapse_debug,
        "unrestricted_quartet": ",".join(unrestricted_quartet),
        "unrestricted_loglik": f"{unrestricted_loglik:.4f}" if math.isfinite(unrestricted_loglik) else "",
        "unrestricted_residual": f"{unrestricted_residual:.4f}" if math.isfinite(unrestricted_residual) else "",
        "replacement_debug": replacement_debug,
        "selector_attempts": selector_attempts_debug,
        "baseline_residual": f"{baseline_residual:.4f}" if math.isfinite(baseline_residual) else "",
        "best_residual": f"{best_residual:.4f}" if math.isfinite(best_residual) else "",
        "residual_debug": residual_debug,
        "ratio_count_observed": ratio_count_observed,
        "ratio_count_groups": ratio_count_groups,
        "quartets_evaluated": evaluated,
        "baseline_loglik": f"{baseline_loglik:.4f}" if math.isfinite(baseline_loglik) else "",
        "gate_decision": gate_decision,
        "best_loglik": f"{loglik:.4f}",
        "baseline_quartet": ",".join(baseline),
        "error": "",
    })
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Direct raw-read allele DB quartet caller prototype")
    parser.add_argument("--score-tsv", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path, help="existing output root used only for chi_R and gene-bin fallback")
    parser.add_argument("--fq-root", type=Path, default=None, help="raw FASTQ root with <sample>_R1/R2_001.fastq.gz")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--genes", default="")
    parser.add_argument("--samples", default="", help="optional comma-separated sample filter")
    parser.add_argument("--statuses", default="ALL", help="comma-separated input row statuses, or ALL")
    parser.add_argument("--gene-config-tsv", type=Path, default=None,
                        help="optional TSV with one row per gene and algorithm/threshold columns overriding global defaults")
    parser.add_argument("--gene-action", choices=("call", "baseline"), default="call",
                        help="call the gene or keep the input baseline quartet; intended mainly for per-gene config profiles")
    parser.add_argument("--selector-iterations", type=int, default=1,
                        help="number of ranked unrestricted quartet proposals to try through replacement gates")
    parser.add_argument("--evidence-backend", choices=("kmer", "alignment"), default="alignment")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-owner-alleles", type=int, default=8)
    parser.add_argument("--no-global-gene-unique", action="store_true", help="do not remove target-gene k-mers also present in other HLA genes")
    parser.add_argument("--max-owner-fraction", type=float, default=0.75)
    parser.add_argument("--max-candidates", type=int, default=8)
    parser.add_argument("--min-candidate-support", type=float, default=5.0)
    parser.add_argument("--max-candidate-pairs", type=int, default=0)
    parser.add_argument("--max-evidence-pairs", type=int, default=2500)
    parser.add_argument("--alignment-candidate-ranker", choices=("em", "support", "support-count-prior"), default="em")
    parser.add_argument("--alignment-count-prior-weight", type=float, default=0.0,
                        help="default-off weight for boosting 2-field alleles represented by many full IMGT records")
    parser.add_argument("--min-pair-informative-kmers", type=float, default=1.0)
    parser.add_argument("--score-scale", type=float, default=0.35)
    parser.add_argument("--alignment-score-scale", type=float, default=12.0)
    parser.add_argument("--alignment-engine", choices=("bwa", "parasail"), default="bwa")
    parser.add_argument("--alignment-seed-k", type=int, default=23)
    parser.add_argument("--alignment-seed-max-owner-alleles", type=int, default=32)
    parser.add_argument("--alignment-prefilter-alleles", type=int, default=48)
    parser.add_argument("--alignment-top-per-pair", type=int, default=16)
    parser.add_argument("--alignment-top-delta", type=float, default=36.0)
    parser.add_argument("--alignment-full-allele-collapse", action="store_true",
                        help="BWA-align against every full IMGT allele record, then collapse per-read evidence to 2-field alleles")
    parser.add_argument("--alignment-min-score-per-base", type=float, default=0.55)
    parser.add_argument("--alignment-missing-penalty", type=float, default=8.0)
    parser.add_argument("--quartet-selector", choices=("read-likelihood", "ratio-count", "constrained-ratio", "support-count", "multiset-likelihood"), default="read-likelihood",
                        help="quartet search objective for unrestricted candidate selection")
    parser.add_argument("--ratio-count-weight", type=float, default=0.20,
                        help="weight of chi_R-aware global abundance fit in --quartet-selector ratio-count")
    parser.add_argument("--ratio-count-epsilon", type=float, default=0.01,
                        help="smoothing added to expected allele/group abundance fractions")
    parser.add_argument("--no-ratio-count-ambiguity-groups", dest="ratio_count_use_ambiguity_groups", action="store_false",
                        help="score ratio-count abundance at exact allele level instead of adaptive ambiguity-group level")
    parser.set_defaults(ratio_count_use_ambiguity_groups=True)
    parser.add_argument("--constrained-ratio-min-dropped-group-observed", type=float, default=0.05,
                        help="with --quartet-selector constrained-ratio, reject replacement that drops a baseline ambiguity group below this observed posterior abundance")
    parser.add_argument("--constrained-ratio-max-dropped-group-copies", type=int, default=1,
                        help="with --quartet-selector constrained-ratio, reject replacement that removes more than this many copies from any baseline ambiguity group; -1 disables")
    parser.add_argument("--baseline-gate-margin", type=float, default=None, help="accept alignment proposal only if proposal log-likelihood beats baseline by this margin")
    parser.add_argument("--baseline-confidence-max-loglik-per-pair", type=float, default=None,
                        help="default off; reject replacements when baseline alignment log-likelihood per evidence pair is at or above this value")
    parser.add_argument("--baseline-multiset-gate", action="store_true", help="accept only R/D reassignments that preserve the baseline 4-copy allele multiset")
    parser.add_argument("--side-copy-gate", choices=("off", "multiset", "no-new-2field"), default="off",
                        help="constrain BWA/EM output to side/copy correction: 'multiset' searches only permutations of the baseline 4-copy multiset; 'no-new-2field' searches quartets using only baseline 2-field alleles")
    parser.add_argument("--reject-balanced-het-to-hom-side", action="store_true",
                        help="reject moves from balanced R/D heterozygotes (ab/ab) to homozygous side splits (aa/bb), a short-read ambiguity pattern")
    parser.add_argument("--ambiguity-collapse-gate", choices=("off", "group-side", "group-multiset", "adaptive-side", "adaptive-multiset"), default="off",
                        help="collapse read-likelihood-indistinguishable alleles and reject proposals that move outside collapsed group constraints")
    parser.add_argument("--ambiguity-collapse-min-pairs", type=int, default=20,
                        help="minimum read pairs where either allele appears before two alleles can be collapsed")
    parser.add_argument("--ambiguity-collapse-min-shared-fraction", type=float, default=0.50,
                        help="minimum fraction of compared read pairs where both alleles appear before collapsing")
    parser.add_argument("--ambiguity-collapse-max-mean-delta", type=float, default=0.35,
                        help="maximum mean per-read log-likelihood delta for collapsed alleles")
    parser.add_argument("--ambiguity-collapse-max-pair-delta", type=float, default=1.25,
                        help="per-read-pair log-likelihood delta threshold used by the robust within-fraction collapse test")
    parser.add_argument("--ambiguity-collapse-min-within-pair-delta-fraction", type=float, default=0.90,
                        help="minimum fraction of compared read pairs with delta <= --ambiguity-collapse-max-pair-delta")
    parser.add_argument("--ambiguity-collapse-adaptive-max-mean-fraction", type=float, default=0.35,
                        help="adaptive collapse cap for mean pair delta, as a fraction of --alignment-missing-penalty")
    parser.add_argument("--ambiguity-collapse-adaptive-max-pair-fraction", type=float, default=0.40,
                        help="adaptive collapse per-read-pair delta threshold, as a fraction of --alignment-missing-penalty")
    parser.add_argument("--ambiguity-collapse-adaptive-min-mean-gap", type=float, default=0.0,
                        help="optional gap in sorted candidate-pair mean deltas used to stop the adaptive low-delta cluster; 0 disables gap splitting")
    parser.add_argument("--independent-sidecopy-allow-tsv", type=Path, default=None,
                        help="optional truth-free TSV with sample,gene[,allow] rows that approve side/copy moves outside collapsed ambiguity groups")
    parser.add_argument("--replacement-gate", choices=("off", "likelihood"), default="off",
                        help="default-off diagnostic allele replacement gate using unrestricted candidate-quartet likelihood and support evidence")
    parser.add_argument("--replacement-min-loglik-per-pair", type=float, default=0.05,
                        help="minimum unrestricted-vs-baseline log-likelihood improvement per evidence pair before allele replacement is accepted")
    parser.add_argument("--replacement-min-support", type=float, default=1.0,
                        help="minimum BWA read-pair support for every new allele introduced by --replacement-gate likelihood")
    parser.add_argument("--replacement-min-em-weight", type=float, default=0.02,
                        help="minimum EM weight for every new allele introduced by --replacement-gate likelihood")
    parser.add_argument("--replacement-max-new-copies", type=int, default=2,
                        help="maximum number of allele copies not present in the baseline quartet; -1 disables this cap")
    parser.add_argument("--replacement-allow-same-ambiguity-multiset", action="store_true",
                        help="allow full allele replacement proposals that preserve the collapsed ambiguity-group multiset")
    parser.add_argument("--replacement-allow-balanced-group-dosage-change", action="store_true",
                        help="allow full allele replacement to alter collapsed group dosage when baseline R/D have balanced group counts")
    parser.add_argument("--residual-gate", choices=("off", "not-worse", "improvement"), default="off",
                        help="gate replacement proposals by mismatch/indel residual against the baseline quartet")
    parser.add_argument("--residual-min-improvement-per-pair", type=float, default=0.0,
                        help="minimum per-read-pair residual improvement required by --residual-gate improvement")
    parser.add_argument("--residual-missing-penalty", type=float, default=20.0,
                        help="residual penalty when a read pair has no alignment to any quartet allele")
    parser.add_argument("--residual-indel-weight", type=float, default=2.0,
                        help="weight for CIGAR insertion/deletion bases in mismatch/indel residual")
    parser.add_argument("--residual-model", choices=("min", "weighted"), default="weighted",
                        help="residual quartet model: min uses the best allele per read pair; weighted uses chi_R/copy mixture weights")
    parser.add_argument("--residual-scale", type=float, default=4.0,
                        help="softmin scale for --residual-model weighted")
    parser.add_argument("--bwa-path", default="bwa")
    parser.add_argument("--bwa-threads", type=int, default=2)
    parser.add_argument("--bwa-primary-only", action="store_true", help="do not ask bwa mem for all secondary alignments")
    parser.add_argument("--raw-gene-enrich", action="store_true", help="for raw FASTQs, select BWA input read pairs by target-gene informative k-mers")
    parser.add_argument("--raw-gene-enrich-gene-unique", dest="raw_gene_enrich_gene_unique", action="store_true",
                        help="when enriching raw FASTQs, require target-gene k-mers to be absent from other configured HLA genes")
    parser.add_argument("--no-raw-gene-enrich-gene-unique", dest="raw_gene_enrich_gene_unique", action="store_false",
                        help="when enriching raw FASTQs, keep target-gene k-mers even if they appear in other configured HLA genes")
    parser.set_defaults(raw_gene_enrich_gene_unique=True)
    parser.add_argument("--raw-gene-enrich-min-kmers", type=float, default=1.0)
    parser.add_argument("--raw-gene-enrich-scan-pairs", type=int, default=0, help="maximum raw read pairs to scan while selecting BWA input; 0 scans until enough selected or EOF")
    parser.add_argument("--em-iterations", type=int, default=40)
    parser.add_argument("--em-min-weight", type=float, default=1e-6)
    parser.add_argument("--max-quartets", type=int, default=0)
    parser.add_argument("--quiet", action="store_true", help="suppress per-row progress messages")
    return parser.parse_args()


def summarize(rows: list[dict[str, object]], summary_path: Path) -> None:
    summary: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        gene = str(row["gene"])
        direct_score = int(row.get("direct_score", 0) or 0)
        baseline_score = int(row.get("baseline_score", 0) or 0)
        direct_multiset_score = int(row.get("direct_multiset_score", 0) or 0)
        baseline_multiset_score = int(row.get("baseline_multiset_score", 0) or 0)
        for label in ("ALL", gene):
            summary[label]["rows"] += 1
            summary[label]["baseline_score"] += baseline_score
            summary[label]["direct_score"] += direct_score
            summary[label]["baseline_multiset_score"] += baseline_multiset_score
            summary[label]["direct_multiset_score"] += direct_multiset_score
            summary[label][str(row.get("status", "unknown"))] += 1
            if row.get("error"):
                summary[label]["errors"] += 1
    out_rows = []
    for label in ["ALL"] + sorted(gene for gene in summary if gene != "ALL"):
        stats = summary[label]
        out_rows.append({
            "gene": label,
            "rows": stats["rows"],
            "baseline_score": stats["baseline_score"],
            "direct_score": stats["direct_score"],
            "delta": stats["direct_score"] - stats["baseline_score"],
            "accuracy": f"{(stats['direct_score'] / (stats['rows'] * 4)) if stats['rows'] else 0.0:.4f}",
            "baseline_multiset_score": stats["baseline_multiset_score"],
            "direct_multiset_score": stats["direct_multiset_score"],
            "delta_multiset": stats["direct_multiset_score"] - stats["baseline_multiset_score"],
            "multiset_accuracy": f"{(stats['direct_multiset_score'] / (stats['rows'] * 4)) if stats['rows'] else 0.0:.4f}",
            "improved_rows": stats["improved"],
            "regressed_rows": stats["regressed"],
            "same_rows": stats["same"],
            "error_rows": stats["errors"],
        })
    write_tsv(summary_path, [
        "gene", "rows", "baseline_score", "direct_score", "delta", "accuracy",
        "baseline_multiset_score", "direct_multiset_score", "delta_multiset", "multiset_accuracy",
        "improved_rows", "regressed_rows", "same_rows", "error_rows",
    ], out_rows)


def main() -> None:
    args = parse_args()
    args.independent_sidecopy_allow = load_independent_sidecopy_allow(args.independent_sidecopy_allow_tsv)
    args.gene_config_profiles = load_gene_config_profiles(args.gene_config_tsv)
    args.gene_config_overrides = ""
    status_filter = None if args.statuses.upper() == "ALL" else {item.strip() for item in args.statuses.split(",") if item.strip()}
    gene_filter = {normalize_gene(item.strip()) for item in args.genes.split(",") if item.strip()}
    sample_filter = {item.strip() for item in args.samples.split(",") if item.strip()}
    input_rows = []
    for row in read_tsv(args.score_tsv):
        gene = normalize_gene(row.get("gene", ""))
        if sample_filter and row.get("sample") not in sample_filter:
            continue
        if gene_filter and gene not in gene_filter:
            continue
        if status_filter is not None and row.get("status") not in status_filter:
            continue
        input_rows.append(row)

    rows: list[dict[str, object]] = []
    genes = {normalize_gene(row["gene"]) for row in input_rows}
    all_representatives = load_all_gene_representatives(args.imgt, genes | BACKGROUND_GENES)
    for gene in sorted(genes):
        gene_args = args_for_gene(args, gene)
        representatives = all_representatives[gene]
        alignment_representatives = None
        allele_record_counts: dict[str, int] | None = None
        if gene_args.evidence_backend == "alignment" and gene_args.alignment_engine == "bwa" and gene_args.alignment_full_allele_collapse:
            alignment_representatives = load_gene_full_allele_representatives(args.imgt, gene)
            allele_record_counts = full_allele_record_counts(alignment_representatives)
        other_gene_representatives = [reps for other_gene, reps in all_representatives.items() if other_gene != gene]
        gene_kmer_owners = build_gene_kmer_owners(
            representatives,
            other_gene_representatives,
            gene_args.k,
            gene_args.max_owner_alleles,
            (not gene_args.no_global_gene_unique) and gene_args.raw_gene_enrich_gene_unique,
        )
        alignment_seed_owners = build_seed_kmer_owners(
            representatives,
            gene_args.alignment_seed_k,
            gene_args.alignment_seed_max_owner_alleles,
        )
        gene_rows = [item for item in input_rows if normalize_gene(item["gene"]) == gene]
        for index, row in enumerate(gene_rows, 1):
            if not args.quiet:
                print(f"[row] {gene} {index}/{len(gene_rows)} {row.get('sample', '')}", file=sys.stderr, flush=True)
            result = score_row(row, gene_args, representatives, gene_kmer_owners, alignment_seed_owners, alignment_representatives, allele_record_counts)
            rows.append(result)
            if not args.quiet:
                print(
                    f"[done] {gene} {index}/{len(gene_rows)} {row.get('sample', '')} "
                    f"baseline={result.get('baseline_score', '')} direct={result.get('direct_score', '')} "
                    f"status={result.get('status', '')} gate={result.get('gate_decision', '')} error={result.get('error', '')}",
                    file=sys.stderr,
                    flush=True,
                )

    fields = [
        "sample", "set", "gene", "gene_config", "baseline_score", "direct_score", "delta_vs_baseline",
        "baseline_multiset_score", "direct_multiset_score", "delta_multiset_vs_baseline", "status",
        "truth_R", "truth_D", "baseline_pred_R", "baseline_pred_D", "direct_pred_R", "direct_pred_D",
        "chi_r", "evidence_backend", "alignment_engine", "read_source", "raw_or_gene_total_pairs", "raw_scan_pairs", "evidence_pairs", "gene_db_alleles",
        "gene_informative_kmers", "candidate_informative_kmers", "candidates", "candidate_support", "em_weights", "em_loglik",
        "ambiguity_groups", "baseline_group_quartet", "direct_group_quartet", "ambiguity_collapse_debug",
        "unrestricted_quartet", "unrestricted_loglik", "unrestricted_residual", "replacement_debug", "selector_attempts",
        "baseline_residual", "best_residual", "residual_debug", "ratio_count_observed", "ratio_count_groups",
        "quartets_evaluated", "baseline_loglik", "gate_decision", "best_loglik", "baseline_quartet", "error",
    ]
    write_tsv(args.out, fields, rows)
    summary_path = args.summary or args.out.with_suffix(".summary.tsv")
    summarize(rows, summary_path)
    print(f"wrote {args.out}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()