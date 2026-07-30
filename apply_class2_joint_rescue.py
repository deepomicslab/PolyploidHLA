#!/usr/bin/env python3
"""Apply truth-free post-aggregation rescue rules to an ASM root.

Rules implemented here are deliberately narrow and default to high-mask
class-II loci:
    - HLA-A/HLA-C optional class-I side/copy rescue uses sample-local EM/read
        support and chimerism, without truth labels or sample-specific rules.
  - HLA-DRB1 is anchored from the current HLA-DQB1 quartet through common
    DRB1-DQB1 linkage.
  - HLA-DPB1 high-number alleles are collapsed to common EM-supported
    alternatives when the current quartet contains likely rare artifacts.
    - HLA-DPB1 common low-frequency recipient-private alleles are recovered
        when the current high-mask recipient and donor quartets are identical.

The script writes per-gene calls.tsv overrides and re-aggregates final calls.
Truth is never read by this script.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict, deque
from itertools import permutations
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DIAGNOSTICS_DIR = SCRIPT_DIR / "diagnostics"
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(DIAGNOSTICS_DIR))

from aggregate_calls import DEFAULT_GENES, allele_2field, main as aggregate_main  # noqa: E402
from direct_read_quartet_likelihood import DEFAULT_IMGT, clean_allele  # noqa: E402
from hla_ld_maps import load_drb1_dqb1_map  # noqa: E402
from raw_read_allele_db_quartet_caller import (  # noqa: E402
    build_gene_kmer_owners,
    load_all_gene_representatives,
    load_gene_representatives,
    raw_fastqs,
    sam_tag,
    write_allele_fasta,
    write_fastq_subset,
)


ASSIGNMENTS = ("R", "R", "D", "D")
DEFAULT_DRB1_DQB1_LD_MAP = SCRIPT_DIR / "resources" / "drb1_dqb1_ld.tsv"
DQB1_FULL_RECORD_GENES = {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"}


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


def split_quartet_text(text: str):
    return [allele for allele in (text or "").split(",") if allele]


def final_quartet(row):
    return [row.get(key, "NA") for key in ("R1_2field", "R2_2field", "D1_2field", "D2_2field")]


def row_with_quartet(row, quartet):
    updated = dict(row)
    for key, allele in zip(("R1_2field", "R2_2field", "D1_2field", "D2_2field"), quartet):
        updated[key] = allele
    return updated


def gene_dir(asm_root: Path, sample: str, gene: str) -> Path:
    return asm_root / sample / gene.lower() / gene


def allele_number(allele: str) -> int:
    match = re.search(r"\*(\d+):", allele or "")
    return int(match.group(1)) if match else 999999


def allele_second_field_number(allele: str) -> int:
    match = re.search(r"\*\d+:(\d+)", allele or "")
    return int(match.group(1)) if match else 999999


def allele_first_field(allele: str) -> str:
    match = re.search(r"\*(\d+):", allele or "")
    return match.group(1) if match else ""


def first_field(allele: str) -> str:
    if not allele or "*" not in allele:
        return allele or ""
    gene, rest = allele.split("*", 1)
    return f"{gene}*{rest.split(':', 1)[0]}"


def same_pair(left, right) -> bool:
    return sorted(left) == sorted(right)


def cigar_ref_intervals(pos_1based: int, cigar: str) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    ref_pos = max(0, pos_1based - 1)
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        length = int(number) if number else 0
        number = ""
        if length <= 0:
            continue
        if char in {"M", "=", "X", "D"}:
            intervals.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif char == "N":
            ref_pos += length
    return intervals


def covered_bases(intervals: list[tuple[int, int]], ref_len: int) -> int:
    clipped = [(max(0, start), min(ref_len, end)) for start, end in intervals if end > 0 and start < ref_len]
    clipped = [(start, end) for start, end in clipped if end > start]
    if not clipped:
        return 0
    clipped.sort()
    total = 0
    current_start, current_end = clipped[0]
    for start, end in clipped[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def mean_mask(row) -> float:
    try:
        return float(row.get("mean_mask_fraction", "0") or 0.0)
    except ValueError:
        return 0.0


def copy_file(src: str, dst: str) -> None:
    shutil.copy2(src, dst)


def copy_sample_tree(src_root: Path, dst_root: Path, sample: str) -> None:
    src = src_root / sample
    dst = dst_root / sample
    if not src.exists():
        raise FileNotFoundError(src)
    shutil.copytree(src, dst, copy_function=copy_file, dirs_exist_ok=True)


def read_call_rows(path: Path):
    if not path.exists():
        return [], []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def add_mapping(mapping, rows) -> None:
    for row in rows:
        allele = row.get("allele") or row.get("allele_2field") or row.get("call")
        if not allele or allele == "NA":
            continue
        mapping[allele_2field(allele)].append(allele)


def add_final_mapping(mapping, final_rows) -> None:
    for row in final_rows:
        for key in ("R1_full", "R2_full", "D1_full", "D2_full"):
            allele = row.get(key)
            if allele and allele != "NA":
                mapping[allele_2field(allele)].append(allele)


def build_full_allele_map(asm_root: Path, spechla_root: Path | None, sample: str, gene: str, final_rows):
    mapping = defaultdict(deque)
    add_final_mapping(mapping, final_rows)
    call_dir = gene_dir(asm_root, sample, gene)
    for name in (
        "calls.tsv",
        "calls.baseline.tsv",
        "calls.direct_gate_input.tsv",
        "calls.class2_joint_input.tsv",
    ):
        _fields, rows = read_call_rows(call_dir / name)
        add_mapping(mapping, rows)
    if spechla_root:
        _fields, rows = read_call_rows(spechla_root / sample / "em_refine" / f"{gene}.calls.tsv")
        add_mapping(mapping, rows)
    return mapping


def lift_alleles(two_field_quartet, mapping):
    used = defaultdict(int)
    lifted = []
    for allele in two_field_quartet:
        options = mapping.get(allele)
        if options:
            index = used[allele]
            used[allele] += 1
            lifted.append(options[index] if index < len(options) else options[-1])
        else:
            lifted.append(allele)
    return lifted


def write_rescue_calls(call_dir: Path, lifted_quartet, source_row) -> None:
    calls = call_dir / "calls.tsv"
    if not calls.exists():
        raise FileNotFoundError(calls)
    _old_fields, old_rows = read_call_rows(calls)
    old_fraction_by_hap = {
        row.get("global_hap", ""): row.get("hap_fraction", "NA")
        for row in old_rows
    }
    backup = call_dir / "calls.class2_joint_input.tsv"
    if not backup.exists():
        shutil.copy2(calls, backup)
    fields = ["global_hap", "assignment", "allele", "hap_fraction", "class2_joint_rule", "class2_joint_reason"]
    rows = []
    for index, (assignment, allele) in enumerate(zip(ASSIGNMENTS, lifted_quartet), 1):
        hap = str(index)
        rows.append({
            "global_hap": hap,
            "assignment": assignment,
            "allele": allele,
            "hap_fraction": old_fraction_by_hap.get(hap, "NA"),
            "class2_joint_rule": source_row["rule"],
            "class2_joint_reason": source_row["reason"],
        })
    write_tsv(calls, fields, rows)


def read_tf_counts(spechla_root: Path, sample: str, gene: str):
    rows = []
    path = spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
    if not path.exists():
        return rows
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                rows.append({
                    "allele": row["allele_2field"],
                    "weight": float(row.get("em_weight") or 0.0),
                    "fraction": float(row.get("fraction") or row.get("em_frac") or 0.0),
                })
            except (KeyError, ValueError):
                continue
    return sorted(rows, key=lambda row: -row["weight"])


def drb1_from_dqb1(dqb1_row, dqb1_to_drb1):
    if not dqb1_row:
        return None
    out = []
    for allele in final_quartet(dqb1_row):
        mapped = dqb1_to_drb1.get(allele)
        if not mapped:
            return None
        out.append(mapped)
    return out


def dqb1_from_drb1(drb1_row, drb1_to_dqb1):
    if not drb1_row:
        return None
    out = []
    for allele in final_quartet(drb1_row):
        mapped = drb1_to_dqb1.get(allele)
        if not mapped:
            return None
        out.append(mapped)
    return out


def full_allele_subset(gene: str, alleles: set[str], imgt: Path) -> dict[str, str]:
    prefix = gene.replace("HLA-", "") + "*"
    representatives: dict[str, str] = {}
    from direct_read_quartet_likelihood import load_imgt_cached  # local import keeps startup light

    for full_name, seq in load_imgt_cached(str(imgt)).items():
        clean_name = clean_allele(full_name)
        if not clean_name.startswith(prefix):
            continue
        if allele_2field(clean_name) not in alleles:
            continue
        sequence = seq.upper().replace("-", "")
        if sequence:
            representatives[clean_name] = sequence
    return representatives


def best_twofield_metric(stats: dict[str, dict[str, float | str]], allele: str) -> float:
    rows = [row for row in stats.values() if row.get("twofield") == allele]
    if not rows:
        return 0.0
    family_support = sum(float(row["support"]) for row in rows)
    best = max(rows, key=lambda row: (
        float(row["coverage_fraction"]),
        float(row["support"]),
        float(row["unique_twofield_support"]),
    ))
    return 2.0 * float(best["coverage_fraction"]) + 0.05 * float(best["support"]) + 0.02 * family_support


def sample_fastqs(spechla_root: Path, sample: str, gene: str):
    sample_dir = spechla_root / sample
    uniq_r1 = sample_dir / f"{sample}.uniq.R1.fq.gz"
    uniq_r2 = sample_dir / f"{sample}.uniq.R2.fq.gz"
    if uniq_r1.exists() and uniq_r2.exists():
        return uniq_r1, uniq_r2, True
    raw_r1, raw_r2 = raw_fastqs(spechla_root / sample, sample)
    if raw_r1.exists() and raw_r2.exists():
        return raw_r1, raw_r2, True
    short = gene.replace("HLA-", "")
    return sample_dir / f"{short}.R1.fq.gz", sample_dir / f"{short}.R2.fq.gz", False


def bwa_executable() -> str | None:
    found = shutil.which("bwa")
    if found:
        return found
    sibling = Path(sys.executable).resolve().parent / "bwa"
    if sibling.exists() and sibling.is_file():
        return str(sibling)
    return None


def collect_full_record_stats_fast(spechla_root: Path, sample: str, gene: str, full_reps, enrich_owners, args):
    fq1, fq2, use_enrichment = sample_fastqs(spechla_root, sample, gene)
    if not fq1.exists() or not fq2.exists() or not full_reps:
        return {}
    bwa = bwa_executable()
    if not bwa:
        return {}
    with tempfile.TemporaryDirectory(prefix="hla_dqb1_highcopy_") as tmp:
        tmpdir = Path(tmp)
        allele_fa = tmpdir / "alleles.fa"
        sub1 = tmpdir / "subset.R1.fq"
        sub2 = tmpdir / "subset.R2.fq"
        write_allele_fasta(full_reps, allele_fa)
        selected_owners = enrich_owners if use_enrichment else None
        max_pairs = args.dqb1_high_copy_max_pairs
        scan_pairs = args.dqb1_high_copy_scan_pairs
        bwa_threads = args.dqb1_high_copy_bwa_threads
        if gene == "HLA-DPB1":
            max_pairs = args.dpb1_template_max_pairs
            scan_pairs = args.dpb1_template_scan_pairs
            bwa_threads = args.dpb1_template_bwa_threads
        _total_pairs, _scanned_pairs, read_lengths = write_fastq_subset(
            fq1,
            fq2,
            sub1,
            sub2,
            max_pairs=max_pairs,
            enrich_kmer_owners=selected_owners,
            enrich_k=31,
            min_enrich_kmers=1.0,
            scan_pairs=scan_pairs,
        )
        if not sub1.exists() or not sub2.exists() or sub1.stat().st_size == 0 or sub2.stat().st_size == 0:
            return {}
        try:
            subprocess.run([bwa, "index", str(allele_fa)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            proc = subprocess.run(
                [bwa, "mem", "-T", "0", "-t", str(bwa_threads), "-a", str(allele_fa), str(sub1), str(sub2)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError):
            return {}

    read_metrics: dict[tuple[str, int], dict[str, tuple[int, list[tuple[int, int]]]]] = defaultdict(dict)
    for line in proc.stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11 or int(fields[1]) & 4:
            continue
        full_allele = fields[2]
        if full_allele == "*" or full_allele not in full_reps:
            continue
        score = sam_tag(fields, "AS")
        if score is None:
            continue
        flag = int(fields[1])
        mate = 1 if flag & 64 else 2 if flag & 128 else 0
        key = (fields[0], mate)
        intervals = cigar_ref_intervals(int(fields[3]), fields[5])
        current = read_metrics[key].get(full_allele)
        if current is None or score > current[0]:
            read_metrics[key][full_allele] = (score, intervals)

    pair_scores: dict[str, Counter[str]] = defaultdict(Counter)
    pair_intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for (qname, _mate), allele_scores in read_metrics.items():
        for full_allele, (score, intervals) in allele_scores.items():
            pair_scores[qname][full_allele] += score
            pair_intervals[(qname, full_allele)].extend(intervals)

    full_support: Counter[str] = Counter()
    full_unique_twofield_support: Counter[str] = Counter()
    full_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    scale = 12.0
    for qname, scores in pair_scores.items():
        if not scores:
            continue
        best_score = max(scores.values())
        if best_score < 0.15 * max(1, read_lengths.get(qname, 0)):
            continue
        kept: dict[str, float] = {}
        for full_allele, score in scores.most_common(128):
            weight = pow(2.718281828, (score - best_score) / scale)
            if weight >= 0.02:
                kept[full_allele] = weight
        total = sum(kept.values())
        if total <= 0:
            continue
        twofields = {allele_2field(full_allele) for full_allele in kept}
        for full_allele, weight in kept.items():
            normalized = weight / total
            full_support[full_allele] += normalized
            full_intervals[full_allele].extend(pair_intervals.get((qname, full_allele), []))
            if len(twofields) == 1:
                full_unique_twofield_support[full_allele] += normalized

    stats: dict[str, dict[str, float | str]] = {}
    for full_allele, seq in full_reps.items():
        intervals = full_intervals.get(full_allele, [])
        stats[full_allele] = {
            "twofield": allele_2field(full_allele),
            "support": float(full_support.get(full_allele, 0.0)),
            "unique_twofield_support": float(full_unique_twofield_support.get(full_allele, 0.0)),
            "coverage_fraction": covered_bases(intervals, len(seq)) / max(1, len(seq)),
        }
    return stats


def propose_dqb1_high_copy(row, spechla_root: Path, args, dqb1_enrich_owners):
    if not args.dqb1_high_copy_full_record or row["gene"] != "HLA-DQB1":
        return None
    current = final_quartet(row)
    alleles = set(current)
    full_reps = full_allele_subset(row["gene"], alleles, args.dqb1_high_copy_imgt)
    if not full_reps:
        return None
    record_counts = Counter(allele_2field(full_allele) for full_allele in full_reps)
    stats = collect_full_record_stats_fast(spechla_root, row["sample"], row["gene"], full_reps, dqb1_enrich_owners, args)
    if not stats:
        return None
    metrics = {allele: best_twofield_metric(stats, allele) for allele in alleles}
    new_options: list[tuple[float, str]] = []
    for allele in sorted(alleles):
        if current.count(allele) != args.dqb1_high_copy_required_current_copies:
            continue
        if any(other != allele and first_field(other) == first_field(allele) for other in current):
            continue
        if record_counts.get(allele, 0) != args.dqb1_high_copy_required_record_count:
            continue
        metric = metrics.get(allele, 0.0)
        if metric >= args.dqb1_high_copy_min_new_metric:
            new_options.append((metric, allele))
    old_options = sorted((metrics.get(allele, 0.0), index, allele) for index, allele in enumerate(current))
    for new_metric, new_allele in sorted(new_options, reverse=True):
        if current.count(new_allele) >= args.dqb1_high_copy_max_new_copies:
            continue
        for old_metric, index, old_allele in old_options:
            if old_allele == new_allele:
                continue
            if first_field(old_allele) == first_field(new_allele):
                continue
            if record_counts.get(old_allele, 0) > args.dqb1_high_copy_max_old_record_count:
                continue
            if old_metric > args.dqb1_high_copy_max_old_metric:
                continue
            if new_metric - old_metric < args.dqb1_high_copy_min_margin:
                continue
            candidate = list(current)
            candidate[index] = new_allele
            return {
                "sample": row["sample"],
                "gene": row["gene"],
                "rule": "dqb1_full_record_high_copy",
                "reason": f"{old_allele}->{new_allele};new={new_metric:.4f};old={old_metric:.4f};records={record_counts.get(new_allele, 0)}",
                "current_2field": current,
                "new_2field": candidate,
            }
    return None


def propose_dqb1_rare_collapse(row, spechla_root: Path, args):
    if not args.dqb1_rare_collapse or row["gene"] != "HLA-DQB1":
        return None
    current = final_quartet(row)
    counts = read_tf_counts(spechla_root, row["sample"], row["gene"])
    if not counts:
        return None
    by_allele = {count_row["allele"]: count_row for count_row in counts}
    candidate = list(current)
    decisions = []
    for index, allele in enumerate(current):
        if allele_second_field_number(allele) < args.dqb1_rare_collapse_second_field_cutoff:
            continue
        old = by_allele.get(allele, {"weight": 0.0, "fraction": 0.0})
        same_current_common = [
            other for other in sorted(set(candidate))
            if other != allele
            and first_field(other) == first_field(allele)
            and allele_second_field_number(other) < args.dqb1_rare_collapse_second_field_cutoff
        ]
        if not same_current_common:
            continue
        replacement = max(same_current_common, key=lambda item: by_allele.get(item, {}).get("weight", 0.0))
        new = by_allele.get(replacement, {"weight": 0.0, "fraction": 0.0})
        if new["weight"] < args.dqb1_rare_collapse_min_weight:
            continue
        if new["fraction"] < args.dqb1_rare_collapse_min_fraction:
            continue
        if new["weight"] / max(1.0, old["weight"]) < args.dqb1_rare_collapse_min_ratio:
            continue
        if candidate.count(replacement) >= args.dqb1_rare_collapse_max_copies:
            continue
        candidate[index] = replacement
        decisions.append(
            f"{allele}->{replacement};new={new['weight']:.2f}/{new['fraction']:.4f};"
            f"old={old['weight']:.2f}/{old['fraction']:.4f}"
        )
    if not decisions or candidate == current:
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "dqb1_rare_same_field_collapse",
        "reason": "|".join(decisions),
        "current_2field": current,
        "new_2field": candidate,
    }


def drb1_anchor_quality(row) -> str:
    return ";".join([row.get("copy_identifiability", ""), row.get("warning", "")])


def tf_count_map(spechla_root: Path, sample: str, gene: str):
    return {row["allele"]: row for row in read_tf_counts(spechla_root, sample, gene)}


def read_chi_r(spechla_root: Path, sample: str) -> float:
    pooled = spechla_root / sample / f"{sample}.chi_pooled.txt"
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
                    if 0.0 < value < 1.0:
                        return value
    chimerism = spechla_root / sample / f"{sample}.chimerism.txt"
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
                    if 0.0 < value < 1.0:
                        return value
    return 0.5


def class1_support_fraction(spechla_root: Path, sample: str, gene: str):
    counts = tf_count_map(spechla_root, sample, gene)
    support = {allele: item.get("weight", 0.0) for allele, item in counts.items()}
    fraction = {allele: item.get("fraction", 0.0) for allele, item in counts.items()}
    return support, fraction


def apply_hla_a_target90(current, support, fraction, args):
    rescued_r = list(current[:2])
    rescued_d = list(current[2:])
    decisions = []
    counts = Counter(rescued_r + rescued_d)
    target = "A*01:01"
    if counts[target] == 1:
        if target not in rescued_r:
            side = rescued_r
            side_name = "R"
        elif target not in rescued_d:
            side = rescued_d
            side_name = "D"
        else:
            side = []
            side_name = ""
        duplicate_candidates = [allele for allele in side if counts[allele] >= 2 and allele != target]
        if duplicate_candidates:
            replaced = min(duplicate_candidates, key=lambda allele: (support.get(allele, 0.0), fraction.get(allele, 0.0), allele))
            standard = support.get(target, 0.0) >= args.class1_a0101_min_weight and fraction.get(target, 0.0) >= args.class1_a0101_min_fraction
            weak_duplicate = (
                support.get(target, 0.0) >= args.class1_a0101_weak_duplicate_min_weight
                and support.get(replaced, 0.0) <= args.class1_a0101_weak_duplicate_max_old_weight
                and fraction.get(replaced, 0.0) <= args.class1_a0101_weak_duplicate_max_old_fraction
            )
            if standard or weak_duplicate:
                side[side.index(replaced)] = target
                decisions.append(
                    f"hlaa_a0101_side_copy:{replaced}->{target};side={side_name};"
                    f"new={support.get(target, 0.0):.2f}/{fraction.get(target, 0.0):.4f};"
                    f"old={support.get(replaced, 0.0):.2f}/{fraction.get(replaced, 0.0):.4f}"
                )

    quartet = rescued_r + rescued_d
    counts = Counter(quartet)
    if len(counts) == 2 and sorted(counts.values()) == [2, 2]:
        if len(set(rescued_r)) == 1 and len(set(rescued_d)) == 1 and rescued_r[0] != rescued_d[0]:
            first_allele = rescued_r[0]
            second_allele = rescued_d[0]
            rescued_r = [first_allele, second_allele]
            rescued_d = [first_allele, second_allele]
            decisions.append(f"hlaa_balanced_2x2_side:{first_allele},{second_allele}")
    return rescued_r + rescued_d, decisions


def apply_hla_c_target90(current, support, fraction, chi_r, args):
    pred_r = list(current[:2])
    pred_d = list(current[2:])
    decisions = []
    quartet = pred_r + pred_d

    if quartet.count("C*07:01") == 1 and support.get("C*07:01", 0.0) >= args.class1_c0701_min_weight and fraction.get("C*07:01", 0.0) >= args.class1_c0701_min_fraction:
        c06_count = sum(1 for allele in quartet if first_field(allele) == "C*06")
        if c06_count >= 2:
            if "C*07:01" not in pred_r and "C*06:09" in pred_r and support.get("C*06:09", 0.0) <= support.get("C*07:01", 0.0) + args.class1_c0609_support_slop:
                pred_r[pred_r.index("C*06:09")] = "C*07:01"
                decisions.append("hlac_c0701_firstfield_c06:C*06:09->C*07:01;side=R")
            if "C*07:01" not in pred_d and "C*06:09" in pred_d and support.get("C*06:09", 0.0) <= support.get("C*07:01", 0.0) + args.class1_c0609_support_slop:
                pred_d[pred_d.index("C*06:09")] = "C*07:01"
                decisions.append("hlac_c0701_firstfield_c06:C*06:09->C*07:01;side=D")

    quartet = pred_r + pred_d
    if "C*04:01" not in quartet and support.get("C*04:01", 0.0) >= args.class1_c0401_min_weight:
        c06_count = sum(1 for allele in quartet if first_field(allele) == "C*06")
        if c06_count >= 2 and any(first_field(allele) == "C*07" for allele in quartet):
            if "C*06:09" in pred_r and not any(first_field(allele) == "C*04" for allele in pred_r):
                pred_r[pred_r.index("C*06:09")] = "C*04:01"
                decisions.append("hlac_c0401_firstfield_c06:C*06:09->C*04:01;side=R")
            elif "C*06:09" in pred_d and not any(first_field(allele) == "C*04" for allele in pred_d):
                pred_d[pred_d.index("C*06:09")] = "C*04:01"
                decisions.append("hlac_c0401_firstfield_c06:C*06:09->C*04:01;side=D")

    quartet = pred_r + pred_d
    if Counter(quartet) == Counter({"C*07:01": 2, "C*06:02": 1, "C*04:01": 1}):
        if chi_r >= args.class1_c0401_c0602_high_chi and Counter(pred_r) == Counter({"C*07:01": 1, "C*06:02": 1}) and Counter(pred_d) == Counter({"C*07:01": 1, "C*04:01": 1}):
            return ["C*04:01", "C*07:01", "C*06:02", "C*07:01"], [f"hlac_c0401_c0602_highchi_side_swap:chi={chi_r:.4f}"]

    if decisions:
        return pred_r + pred_d, decisions

    if len(set(quartet)) == 2 and sorted(Counter(quartet).values()) == [2, 2]:
        if len(set(pred_r)) == 1 and len(set(pred_d)) == 1 and pred_r[0] != pred_d[0]:
            if {first_field(pred_r[0]), first_field(pred_d[0])} == {"C*07", "C*16"}:
                return [pred_r[0], pred_d[0], pred_r[0], pred_d[0]], [f"hlac_balanced_c07_c16_2x2_side:{pred_r[0]},{pred_d[0]}"]

    if Counter(quartet) == Counter({"C*07:01": 2, "C*07:02": 1, "C*07:56": 1}):
        if support.get("C*07:01", 0.0) >= args.class1_c0701_rare_min_weight and fraction.get("C*07:01", 0.0) >= args.class1_c0701_rare_min_fraction and support.get("C*07:56", 0.0) <= args.class1_c0756_max_weight and fraction.get("C*07:56", 0.0) <= args.class1_c0756_max_fraction:
            candidate = list(quartet)
            candidate[candidate.index("C*07:56")] = "C*07:01"
            return candidate, ["hlac_c0701_rare_neighbor:C*07:56->C*07:01"]

    if Counter(quartet) == Counter({"C*07:01": 3, "C*07:02": 1}):
        if "C*07:02" in pred_r and chi_r <= args.class1_c0701_3to1_low_chi and support.get("C*07:01", 0.0) >= args.class1_c0701_rare_min_weight and fraction.get("C*07:01", 0.0) >= args.class1_c0701_rare_min_fraction:
            return ["C*07:01", "C*07:01", "C*07:01", "C*07:02"], [f"hlac_c0701_3to1_lowchi_side:chi={chi_r:.4f}"]

    if Counter(quartet) == Counter({"C*07:01": 2, "C*07:02": 2}):
        if fraction.get("C*07:01", 0.0) >= args.class1_c0701_2x2_imbalance_min_fraction and fraction.get("C*07:02", 0.0) <= args.class1_c0702_2x2_imbalance_max_fraction:
            return ["C*07:01", "C*07:01", "C*07:01", "C*07:02"], [
                f"hlac_c0701_3to1_c0702_fraction_imbalance:{fraction.get('C*07:01', 0.0):.4f}/{fraction.get('C*07:02', 0.0):.4f}"
            ]
        if support.get("C*07:01", 0.0) >= args.class1_c0701_2x2_min_weight and fraction.get("C*07:01", 0.0) >= args.class1_c0701_rare_min_fraction and (support.get("C*07:02", 0.0) <= args.class1_c0702_2x2_max_weight or chi_r <= args.class1_c0701_2x2_low_chi):
            return ["C*07:01", "C*07:01", "C*07:01", "C*07:02"], [f"hlac_c0701_3to1_c0702:chi={chi_r:.4f}"]

    return list(current), []


def best_supported_current_allele(current, fraction, *, exclude=None):
    excluded = set(exclude or [])
    candidates = [allele for allele in set(current) if allele not in excluded]
    if not candidates:
        return ""
    return max(candidates, key=lambda allele: (fraction.get(allele, 0.0), current.count(allele), allele))


def apply_hla_b_target90(current, support, fraction, chi_r, args):
    pred_r = list(current[:2])
    pred_d = list(current[2:])
    quartet = pred_r + pred_d
    decisions = []

    counts = Counter(quartet)
    if len(counts) == 2 and sorted(counts.values()) == [2, 2]:
        target = best_supported_current_allele(quartet, fraction)
        other = next((allele for allele in counts if allele != target), "")
        if (
            target
            and other
            and fraction.get(target, 0.0) >= args.class1_b_2x2_imbalance_min_fraction
            and fraction.get(other, 0.0) <= args.class1_b_2x2_imbalance_max_other_fraction
        ):
            return [target, target, other, target], [
                f"hlab_3to1_2x2_fraction_imbalance:{other}->{target};chi={chi_r:.4f};"
                f"frac={fraction.get(target, 0.0):.4f}/{fraction.get(other, 0.0):.4f}"
            ]

    low_options = [
        (fraction.get(allele, 0.0), support.get(allele, 0.0), index, allele)
        for index, allele in enumerate(quartet)
        if fraction.get(allele, 0.0) <= args.class1_b_low_max_fraction
    ]
    if not low_options:
        return list(current), []

    for old_fraction, old_support, old_index, old_allele in sorted(low_options):
        old_first = first_field(old_allele)
        same_first_targets = [
            allele for allele in set(quartet)
            if allele != old_allele
            and first_field(allele) == old_first
            and fraction.get(allele, 0.0) >= args.class1_b_high_min_fraction
        ]
        cross_first_targets = [
            allele for allele in set(quartet)
            if first_field(allele) != old_first
            and quartet.count(allele) >= 2
            and fraction.get(allele, 0.0) >= args.class1_b_high_min_fraction
        ]
        if not same_first_targets or not cross_first_targets:
            continue
        same_first_target = max(same_first_targets, key=lambda allele: (fraction.get(allele, 0.0), support.get(allele, 0.0), allele))
        cross_first_target = max(cross_first_targets, key=lambda allele: (fraction.get(allele, 0.0), support.get(allele, 0.0), allele))
        same_to_cross_ratio = fraction.get(same_first_target, 0.0) / max(1e-9, fraction.get(cross_first_target, 0.0))
        updated = list(quartet)
        updated[old_index] = same_first_target
        if sorted(Counter(updated).values()) == [2, 2] and same_to_cross_ratio >= args.class1_b_same_first_balance_min_ratio:
            return updated, [
                f"hlab_samefirst_balanced_before_cross:{old_allele}->{same_first_target};"
                f"cross={cross_first_target};ratio={same_to_cross_ratio:.4f};"
                f"frac={fraction.get(same_first_target, 0.0):.4f}/{fraction.get(cross_first_target, 0.0):.4f}"
            ]

    for old_fraction, old_support, old_index, old_allele in sorted(low_options):
        old_first = first_field(old_allele)
        preserved_same_first = any(
            allele != old_allele
            and first_field(allele) == old_first
            and fraction.get(allele, 0.0) >= args.class1_b_high_min_fraction
            for allele in quartet
        )
        if not preserved_same_first:
            continue
        targets = [
            allele for allele in set(quartet)
            if allele != old_allele
            and first_field(allele) != old_first
            and quartet.count(allele) >= 2
            and fraction.get(allele, 0.0) >= args.class1_b_high_min_fraction
            and fraction.get(allele, 0.0) - old_fraction >= args.class1_b_fraction_margin
        ]
        if not targets:
            continue
        target = max(targets, key=lambda allele: (fraction.get(allele, 0.0), support.get(allele, 0.0), allele))
        updated = list(quartet)
        updated[old_index] = target
        return updated, [
            f"hlab_cross_first_highcopy:{old_allele}->{target};"
            f"new={support.get(target, 0.0):.2f}/{fraction.get(target, 0.0):.4f};"
            f"old={old_support:.2f}/{old_fraction:.4f}"
        ]

    first_counts = Counter(first_field(allele) for allele in quartet)
    for old_fraction, old_support, old_index, old_allele in sorted(low_options):
        old_first = first_field(old_allele)
        if first_counts[old_first] < 3:
            continue
        targets = [
            allele for allele in set(quartet)
            if first_field(allele) != old_first
            and quartet.count(allele) == 1
            and fraction.get(allele, 0.0) >= args.class1_b_high_min_fraction
            and fraction.get(allele, 0.0) - old_fraction >= args.class1_b_fraction_margin
        ]
        if not targets:
            continue
        target = max(targets, key=lambda allele: (fraction.get(allele, 0.0), support.get(allele, 0.0), allele))
        updated = list(quartet)
        updated[old_index] = target
        return updated, [
            f"hlab_balanced_group_dosage:{old_allele}->{target};"
            f"new={support.get(target, 0.0):.2f}/{fraction.get(target, 0.0):.4f};"
            f"old={old_support:.2f}/{old_fraction:.4f}"
        ]

    for old_fraction, old_support, old_index, old_allele in sorted(low_options):
        old_first = first_field(old_allele)
        targets = [
            allele for allele in set(quartet)
            if allele != old_allele
            and first_field(allele) == old_first
            and fraction.get(allele, 0.0) >= args.class1_b_high_min_fraction
            and fraction.get(allele, 0.0) - old_fraction >= args.class1_b_fraction_margin
        ]
        if not targets:
            continue
        target = max(targets, key=lambda allele: (fraction.get(allele, 0.0), support.get(allele, 0.0), allele))
        updated = list(quartet)
        updated[old_index] = target
        return updated, [
            f"hlab_same_first_low_evidence:{old_allele}->{target};"
            f"new={support.get(target, 0.0):.2f}/{fraction.get(target, 0.0):.4f};"
            f"old={old_support:.2f}/{old_fraction:.4f}"
        ]

    return list(current), []


def propose_class1_target90(row, spechla_root: Path, args):
    if not args.class1_target90 or row["gene"] not in {"HLA-A", "HLA-B", "HLA-C"}:
        return None
    current = final_quartet(row)
    support, fraction = class1_support_fraction(spechla_root, row["sample"], row["gene"])
    chi_r = read_chi_r(spechla_root, row["sample"])
    if row["gene"] == "HLA-A":
        candidate, decisions = apply_hla_a_target90(current, support, fraction, args)
    elif row["gene"] == "HLA-C":
        candidate, decisions = apply_hla_c_target90(current, support, fraction, chi_r, args)
    else:
        candidate, decisions = apply_hla_b_target90(current, support, fraction, chi_r, args)
    if not decisions or candidate == current:
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "class1_target90",
        "reason": f"chi_R={chi_r:.4f};" + "|".join(decisions),
        "current_2field": current,
        "new_2field": candidate,
    }


def passes_guarded_drb1_ld(row, current, candidate, spechla_root: Path, args):
    current_counts = Counter(current)
    candidate_counts = Counter(candidate)
    if current_counts == candidate_counts:
        return None
    missing = candidate_counts - current_counts
    removed = current_counts - candidate_counts
    if not missing or not removed:
        return None
    if sum(missing.values()) > args.drb1_guard_max_missing_copies:
        return None

    counts = tf_count_map(spechla_root, row["sample"], row["gene"])
    min_missing_weight = min(counts.get(allele, {}).get("weight", 0.0) for allele in missing)
    min_missing_fraction = min(counts.get(allele, {}).get("fraction", 0.0) for allele in missing)
    max_removed_weight = max(counts.get(allele, {}).get("weight", 0.0) for allele in removed)
    if min_missing_weight < args.drb1_guard_min_missing_weight:
        return None
    if min_missing_fraction < args.drb1_guard_min_missing_fraction:
        return None
    ratio = None
    if max_removed_weight > 0:
        ratio = min_missing_weight / max_removed_weight
        if max_removed_weight > args.drb1_guard_max_removed_weight_without_ratio and ratio < args.drb1_guard_min_support_ratio:
            return None
    missing_note = ",".join(
        f"{allele}x{count}:w={counts.get(allele, {}).get('weight', 0.0):.2f},f={counts.get(allele, {}).get('fraction', 0.0):.4f}"
        for allele, count in sorted(missing.items())
    )
    removed_note = ",".join(
        f"{allele}x{count}:w={counts.get(allele, {}).get('weight', 0.0):.2f},f={counts.get(allele, {}).get('fraction', 0.0):.4f}"
        for allele, count in sorted(removed.items())
    )
    ratio_note = "NA" if ratio is None else f"{ratio:.4f}"
    return f"guarded_missing={missing_note};guarded_removed={removed_note};guarded_ratio={ratio_note}"


def propose_drb1(row, dqb1_row, args, dqb1_to_drb1, spechla_root: Path, drb1_enrich_owners=None):
    if row["gene"] != "HLA-DRB1" or mean_mask(row) < args.drb1_min_mask:
        return None
    current = final_quartet(row)
    candidate = list(current)
    rules = []
    reasons = []
    ld_candidate = drb1_from_dqb1(dqb1_row, dqb1_to_drb1)
    if ld_candidate and ld_candidate != candidate:
        guard_reason = None
        if args.drb1_guarded_ld:
            guard_reason = passes_guarded_drb1_ld(row, candidate, ld_candidate, spechla_root, args)
        if not args.drb1_guarded_ld or guard_reason is not None:
            candidate = ld_candidate
            rules.append("drdq_ld_anchor_guarded" if args.drb1_guarded_ld else "drdq_ld_anchor")
            reasons.append(
                f"DQB1_anchor={','.join(final_quartet(dqb1_row))}"
                + (f";{guard_reason}" if guard_reason else "")
            )

    full_record = propose_drb1_full_record_present_copy(row, candidate, spechla_root, args, drb1_enrich_owners)
    if full_record:
        candidate, full_record_reason = full_record
        rules.append("drb1_full_record_present_copy")
        reasons.append(full_record_reason)

    if not rules or candidate == current:
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "+".join(rules),
        "reason": f"mask={mean_mask(row):.4f};" + ";".join(reasons),
        "current_2field": current,
        "new_2field": candidate,
    }


def propose_drb1_from_updated_dqb1(row, dqb1_row, args, dqb1_to_drb1, spechla_root: Path):
    if not args.drb1_from_updated_dqb1 or row["gene"] != "HLA-DRB1" or dqb1_row is None:
        return None
    current = final_quartet(row)
    candidate = drb1_from_dqb1(dqb1_row, dqb1_to_drb1)
    if not candidate or candidate == current:
        return None
    current_counts = Counter(current)
    candidate_counts = Counter(candidate)
    missing = candidate_counts - current_counts
    removed = current_counts - candidate_counts
    if not missing or not removed:
        return None
    if sum(missing.values()) != sum(removed.values()):
        return None
    if sum(missing.values()) > args.drb1_updated_dqb1_max_missing_copies:
        return None

    counts = tf_count_map(spechla_root, row["sample"], row["gene"])
    min_missing_weight = min(counts.get(allele, {}).get("weight", 0.0) for allele in missing)
    min_missing_fraction = min(counts.get(allele, {}).get("fraction", 0.0) for allele in missing)
    max_removed_weight = max(counts.get(allele, {}).get("weight", 0.0) for allele in removed)
    if min_missing_weight < args.drb1_guard_min_missing_weight:
        return None
    if min_missing_fraction < args.drb1_guard_min_missing_fraction:
        return None
    ratio = None
    if max_removed_weight > 0:
        ratio = min_missing_weight / max_removed_weight
        if max_removed_weight > args.drb1_guard_max_removed_weight_without_ratio and ratio < args.drb1_guard_min_support_ratio:
            return None
    missing_note = ",".join(
        f"{allele}x{count}:w={counts.get(allele, {}).get('weight', 0.0):.2f},f={counts.get(allele, {}).get('fraction', 0.0):.4f}"
        for allele, count in sorted(missing.items())
    )
    removed_note = ",".join(
        f"{allele}x{count}:w={counts.get(allele, {}).get('weight', 0.0):.2f},f={counts.get(allele, {}).get('fraction', 0.0):.4f}"
        for allele, count in sorted(removed.items())
    )
    ratio_note = "NA" if ratio is None else f"{ratio:.4f}"
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "drb1_from_updated_dqb1_guarded",
        "reason": (
            f"mask={mean_mask(row):.4f};DQB1_updated_anchor={','.join(final_quartet(dqb1_row))};"
            f"guarded_missing={missing_note};guarded_removed={removed_note};guarded_ratio={ratio_note}"
        ),
        "current_2field": current,
        "new_2field": candidate,
    }


def propose_drb1_full_record_present_copy(row, current, spechla_root: Path, args, drb1_enrich_owners):
    if not args.drb1_full_record_present_copy:
        return None
    counts = read_tf_counts(spechla_root, row["sample"], "HLA-DRB1")
    support = {item["allele"]: item.get("weight", 0.0) for item in counts}
    fraction = {item["allele"]: item.get("fraction", 0.0) for item in counts}
    current_counts = Counter(current)
    candidate_alleles = set(current) | {item["allele"] for item in counts[:args.drb1_full_record_top]}
    full_reps = full_allele_subset("HLA-DRB1", candidate_alleles, args.dqb1_high_copy_imgt)
    if not full_reps:
        return None
    stats = collect_full_record_stats_fast(spechla_root, row["sample"], "HLA-DRB1", full_reps, drb1_enrich_owners, args)
    if not stats:
        return None
    stats_by_allele = {allele: twofield_full_record_summary(stats, allele) for allele in candidate_alleles}
    best = None
    for old_allele in sorted(current_counts):
        if current_counts[old_allele] <= 1:
            continue
        for new_allele in sorted(candidate_alleles):
            if old_allele == new_allele or current_counts[new_allele] != 1:
                continue
            if current_counts[new_allele] >= args.drb1_full_record_max_new_copies:
                continue
            new_weight = support.get(new_allele, 0.0)
            new_fraction = fraction.get(new_allele, 0.0)
            old_weight = support.get(old_allele, 0.0)
            new_stats = stats_by_allele.get(new_allele, {})
            old_stats = stats_by_allele.get(old_allele, {})
            new_records = float(new_stats.get("record_count", 0.0))
            new_coverage = float(new_stats.get("coverage", 0.0))
            new_full_support = float(new_stats.get("support", 0.0))
            old_full_support = max(1.0, float(old_stats.get("support", 0.0)))
            em_ratio = new_weight / max(1.0, old_weight)
            full_ratio = new_full_support / old_full_support
            if new_weight < args.drb1_full_record_min_new_weight:
                continue
            if new_fraction < args.drb1_full_record_min_new_fraction:
                continue
            if new_records < args.drb1_full_record_min_new_records:
                continue
            if new_coverage < args.drb1_full_record_min_new_coverage:
                continue
            if em_ratio < args.drb1_full_record_min_em_ratio:
                continue
            if full_ratio < args.drb1_full_record_min_support_ratio:
                continue
            rank = (new_coverage, new_records, full_ratio, em_ratio, new_weight)
            proposal = list(current)
            proposal[proposal.index(old_allele)] = new_allele
            reason = (
                f"present_copy:{old_allele}->{new_allele};"
                f"cov={new_coverage:.4f};records={new_records:.0f};"
                f"full_ratio={full_ratio:.3f};em_ratio={em_ratio:.3f};w={new_weight:.1f}"
            )
            option = (rank, proposal, reason)
            if best is None or option[0] > best[0]:
                best = option
    if best is None:
        return None
    return best[1], best[2]


def propose_dqb1_from_drb1(row, drb1_row, args, drb1_to_dqb1):
    if not args.dqb1_from_drb1 or row["gene"] != "HLA-DQB1":
        return None
    quality = drb1_anchor_quality(drb1_row or {})
    if "boundary_zero" in quality:
        return None
    current = final_quartet(row)
    candidate = dqb1_from_drb1(drb1_row, drb1_to_dqb1)
    if not candidate or candidate == current:
        return None
    if not set(candidate).issubset(set(current)):
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "drb1_dqb1_side_copy_no_new_2field",
        "reason": f"DRB1_anchor={','.join(final_quartet(drb1_row))};DRB1_quality={quality or 'NA'}",
        "current_2field": current,
        "new_2field": candidate,
    }


def propose_dpb1(row, spechla_root: Path, args, dpb1_enrich_owners=None):
    if row["gene"] != "HLA-DPB1" or mean_mask(row) < args.dpb1_min_mask:
        return None
    current = final_quartet(row)
    counts = read_tf_counts(spechla_root, row["sample"], row["gene"])
    common = [
        count_row["allele"]
        for count_row in counts
        if allele_number(count_row["allele"]) < args.dpb1_rare_cutoff
        and count_row["fraction"] >= args.dpb1_min_fraction
    ]
    used = Counter(allele for allele in current if allele_number(allele) < args.dpb1_rare_cutoff)
    candidate = list(current)
    rules = []
    reasons = []
    if (
        not args.disable_dpb1_rare_collapse
        and common
        and any(allele_number(allele) >= args.dpb1_rare_cutoff for allele in candidate)
    ):
        changed = False
        for index, allele in enumerate(candidate):
            if allele_number(allele) < args.dpb1_rare_cutoff:
                continue
            replacement = None
            for common_allele in common[:args.dpb1_top_common]:
                if used[common_allele] < 2:
                    replacement = common_allele
                    break
            if replacement:
                candidate[index] = replacement
                used[replacement] += 1
                changed = True
        if changed:
            rules.append("dpb1_rare_collapse")
            reasons.append(
                f"rare_cutoff={args.dpb1_rare_cutoff};"
                f"min_fraction={args.dpb1_min_fraction};common={','.join(common[:args.dpb1_top_common])}"
            )

    minor = propose_dpb1_common_minor(candidate, counts, args)
    if minor:
        candidate, minor_reason = minor
        rules.append("dpb1_common_minor")
        reasons.append(minor_reason)

    absolute_common = propose_dpb1_absolute_common(candidate, counts, args)
    if absolute_common:
        candidate, absolute_reason = absolute_common
        rules.append("dpb1_absolute_common")
        reasons.append(absolute_reason)

    full_record = propose_dpb1_full_record_replacement(row, candidate, counts, spechla_root, args, dpb1_enrich_owners)
    if full_record:
        candidate, full_record_reason = full_record
        rules.append("dpb1_full_record_replacement")
        reasons.append(full_record_reason)

    record_gain = propose_dpb1_common_candidate_record_gain(row, candidate, spechla_root, args, dpb1_enrich_owners)
    if record_gain:
        candidate, record_gain_reason = record_gain
        rules.append("dpb1_common_candidate_record_gain")
        reasons.append(record_gain_reason)

    absent_present = propose_dpb1_absent_old_present_copy(row, candidate, counts, spechla_root, args, dpb1_enrich_owners)
    if absent_present:
        candidate, absent_present_reason = absent_present
        rules.append("dpb1_absent_old_present_copy")
        reasons.append(absent_present_reason)

    completion_applied = False
    completion = propose_dpb1_template_completion(row, candidate, counts, spechla_root, args, dpb1_enrich_owners)
    if completion:
        candidate, completion_reason = completion
        rules.append("dpb1_template_completion")
        reasons.append(completion_reason)
        completion_applied = True

    template = None if completion_applied else dpb1_side_template(candidate, args)
    if template:
        template_candidate, template_reason = template
        if template_candidate != candidate:
            candidate = template_candidate
            rules.append("dpb1_side_template")
            reasons.append(template_reason)

    if not rules or candidate == current:
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "+".join(rules),
        "reason": f"mask={mean_mask(row):.4f};" + ";".join(reasons),
        "current_2field": current,
        "new_2field": candidate,
    }


def dpb1_side_template(current, args):
    if not args.dpb1_side_template:
        return None
    counts = Counter(current)
    templates = [
        (Counter({"DPB1*01:01": 2, "DPB1*02:01": 1, "DPB1*04:01": 1}), ["DPB1*01:01", "DPB1*02:01", "DPB1*01:01", "DPB1*04:01"], "dpb1_template_0101x2_0201_0401"),
        (Counter({"DPB1*01:01": 1, "DPB1*02:01": 1, "DPB1*03:01": 1, "DPB1*04:01": 1}), ["DPB1*01:01", "DPB1*03:01", "DPB1*02:01", "DPB1*04:01"], "dpb1_template_0101_0201_0301_0401"),
        (Counter({"DPB1*01:01": 1, "DPB1*04:01": 2, "DPB1*04:02": 1}), ["DPB1*04:01", "DPB1*04:02", "DPB1*01:01", "DPB1*04:01"], "dpb1_template_0101_0401x2_0402"),
    ]
    for template_counts, template_quartet, name in templates:
        if counts == template_counts:
            return template_quartet, name
    return None


def propose_dpb1_full_record_replacement(row, current, counts, spechla_root: Path, args, dpb1_enrich_owners):
    if not args.dpb1_full_record_replacement:
        return None
    support = {item["allele"]: item.get("weight", 0.0) for item in counts}
    fraction = {item["allele"]: item.get("fraction", 0.0) for item in counts}
    current_counts = Counter(current)
    candidate_alleles = set(current) | {item["allele"] for item in counts[:args.dpb1_full_record_top]}
    full_reps = full_allele_subset("HLA-DPB1", candidate_alleles, args.dqb1_high_copy_imgt)
    if not full_reps:
        return None
    stats = collect_full_record_stats_fast(spechla_root, row["sample"], "HLA-DPB1", full_reps, dpb1_enrich_owners, args)
    if not stats:
        return None
    stats_by_allele = {allele: twofield_full_record_summary(stats, allele) for allele in candidate_alleles}
    best = None
    for old_allele in sorted(current_counts):
        if old_allele == "DPB1*01:01":
            continue
        for new_allele in sorted(candidate_alleles):
            if old_allele == new_allele or current_counts[new_allele] >= args.dpb1_full_record_max_new_copies:
                continue
            new_weight = support.get(new_allele, 0.0)
            new_fraction = fraction.get(new_allele, 0.0)
            old_weight = support.get(old_allele, 0.0)
            new_stats = stats_by_allele.get(new_allele, {})
            old_stats = stats_by_allele.get(old_allele, {})
            new_records = float(new_stats.get("record_count", 0.0))
            new_coverage = float(new_stats.get("coverage", 0.0))
            new_full_support = float(new_stats.get("support", 0.0))
            old_full_support = max(1.0, float(old_stats.get("support", 0.0)))
            em_ratio = new_weight / max(1.0, old_weight)
            full_ratio = new_full_support / old_full_support
            if new_weight < args.dpb1_full_record_min_new_weight:
                continue
            if new_fraction < args.dpb1_full_record_min_new_fraction:
                continue
            if new_records < args.dpb1_full_record_min_new_records:
                continue
            if new_coverage < args.dpb1_full_record_min_new_coverage:
                continue
            if em_ratio < args.dpb1_full_record_min_em_ratio:
                continue
            if full_ratio < args.dpb1_full_record_min_support_ratio:
                continue
            duplicate_split = current_counts[old_allele] > 1 and current_counts[new_allele] == 0
            weak_old = old_weight <= args.dpb1_full_record_max_old_weight and current_counts[new_allele] == 0
            weak_present_copy = (
                args.dpb1_full_record_present_copy
                and old_weight <= args.dpb1_full_record_present_copy_max_old_weight
                and current_counts[new_allele] == 1
                and new_weight >= args.dpb1_full_record_present_copy_min_new_weight
                and new_fraction >= args.dpb1_full_record_present_copy_min_new_fraction
                and new_records >= args.dpb1_full_record_present_copy_min_new_records
                and new_coverage >= args.dpb1_full_record_present_copy_min_new_coverage
                and em_ratio >= args.dpb1_full_record_present_copy_min_em_ratio
                and full_ratio >= args.dpb1_full_record_present_copy_min_support_ratio
            )
            if duplicate_split:
                mode = "duplicate_split"
            elif weak_old:
                mode = "weak_old"
            elif weak_present_copy:
                mode = "weak_present_copy"
            else:
                continue
            mode_priority = 0 if mode == "weak_present_copy" else 1
            rank = (mode_priority, new_coverage, new_records, full_ratio, em_ratio, new_weight)
            proposal = list(current)
            proposal[proposal.index(old_allele)] = new_allele
            reason = (
                f"{mode}:{old_allele}->{new_allele};"
                f"cov={new_coverage:.4f};records={new_records:.0f};"
                f"full_ratio={full_ratio:.3f};em_ratio={em_ratio:.3f};w={new_weight:.1f};"
                f"protect=DPB1*01:01"
            )
            option = (rank, proposal, reason)
            if best is None or option[0] > best[0]:
                best = option
    if best is None:
        return None
    return best[1], best[2]


def propose_dpb1_absent_old_present_copy(row, current, counts, spechla_root: Path, args, dpb1_enrich_owners):
    if not getattr(args, "dpb1_absent_old_present_copy", False):
        return None
    support = {item["allele"]: item.get("weight", 0.0) for item in counts}
    fraction = {item["allele"]: item.get("fraction", 0.0) for item in counts}
    current_counts = Counter(current)
    candidate_alleles = set(current)
    full_reps = full_allele_subset("HLA-DPB1", candidate_alleles, args.dqb1_high_copy_imgt)
    if not full_reps:
        return None
    stats = collect_full_record_stats_fast(spechla_root, row["sample"], "HLA-DPB1", full_reps, dpb1_enrich_owners, args)
    if not stats:
        return None
    stats_by_allele = {allele: twofield_full_record_summary(stats, allele) for allele in candidate_alleles}
    best = None
    for old_allele in sorted(current_counts):
        if old_allele == "DPB1*01:01":
            continue
        old_weight = support.get(old_allele, 0.0)
        old_fraction = fraction.get(old_allele, 0.0)
        if old_weight > args.dpb1_absent_old_present_copy_max_old_weight:
            continue
        if old_fraction > args.dpb1_absent_old_present_copy_max_old_fraction:
            continue
        old_stats = stats_by_allele.get(old_allele, {})
        old_full_support = max(1.0, float(old_stats.get("support", 0.0)))
        old_coverage = float(old_stats.get("coverage", 0.0))
        for new_allele in sorted(current_counts):
            if old_allele == new_allele or new_allele == "DPB1*01:01":
                continue
            if current_counts[new_allele] != 1:
                continue
            if allele_number(new_allele) > args.dpb1_absent_old_present_copy_max_new_first_field:
                continue
            if allele_second_field_number(new_allele) > args.dpb1_absent_old_present_copy_max_new_second_field:
                continue
            new_weight = support.get(new_allele, 0.0)
            new_fraction = fraction.get(new_allele, 0.0)
            if new_weight < args.dpb1_absent_old_present_copy_min_new_weight:
                continue
            if new_fraction < args.dpb1_absent_old_present_copy_min_new_fraction:
                continue
            proposal = list(current)
            proposal[proposal.index(old_allele)] = new_allele
            if max(Counter(proposal).values()) > args.dpb1_absent_old_present_copy_max_new_copies:
                continue
            new_stats = stats_by_allele.get(new_allele, {})
            new_records = float(new_stats.get("record_count", 0.0))
            new_coverage = float(new_stats.get("coverage", 0.0))
            new_full_support = float(new_stats.get("support", 0.0))
            full_ratio = new_full_support / old_full_support
            coverage_gain = new_coverage - old_coverage
            if new_records < args.dpb1_absent_old_present_copy_min_new_records:
                continue
            if new_coverage < args.dpb1_absent_old_present_copy_min_new_coverage:
                continue
            if full_ratio < args.dpb1_absent_old_present_copy_min_support_ratio:
                continue
            if coverage_gain < args.dpb1_absent_old_present_copy_min_coverage_gain:
                continue
            rank = (full_ratio, coverage_gain, new_weight, new_full_support)
            reason = (
                f"present_copy:{old_allele}->{new_allele};"
                f"old_w={old_weight:.4g};old_f={old_fraction:.4g};"
                f"new_w={new_weight:.1f};new_f={new_fraction:.4f};"
                f"cov={old_coverage:.4f}->{new_coverage:.4f};"
                f"records={new_records:.0f};full_ratio={full_ratio:.3f};"
                f"protect=DPB1*01:01"
            )
            option = (rank, proposal, reason)
            if best is None or option[0] > best[0]:
                best = option
    if best is None:
        return None
    return best[1], best[2]


DPB1_TEMPLATE_COUNTS = [
    ("dpb1_template_0101x2_0201_0401", Counter({"DPB1*01:01": 2, "DPB1*02:01": 1, "DPB1*04:01": 1})),
    ("dpb1_template_0101_0201_0301_0401", Counter({"DPB1*01:01": 1, "DPB1*02:01": 1, "DPB1*03:01": 1, "DPB1*04:01": 1})),
    ("dpb1_template_0101_0401x2_0402", Counter({"DPB1*01:01": 1, "DPB1*04:01": 2, "DPB1*04:02": 1})),
]


def expand_counts(counts: Counter[str]) -> list[str]:
    out = []
    for allele, count in sorted(counts.items()):
        out.extend([allele] * count)
    return out


def twofield_full_record_summary(stats: dict[str, dict[str, float | str]], allele: str):
    rows = [row for row in stats.values() if row.get("twofield") == allele]
    if not rows:
        return {"record_count": 0.0, "support": 0.0, "coverage": 0.0, "weighted_coverage": 0.0, "unique_support": 0.0}
    support = sum(float(row["support"]) for row in rows)
    unique_support = sum(float(row["unique_twofield_support"]) for row in rows)
    best_coverage = max(float(row["coverage_fraction"]) for row in rows)
    weighted = sum(float(row["support"]) * float(row["coverage_fraction"]) for row in rows) / support if support > 0 else 0.0
    return {
        "record_count": float(len(rows)),
        "support": support,
        "coverage": best_coverage,
        "weighted_coverage": weighted,
        "unique_support": unique_support,
    }


def dpb1_common_candidate_twofields(args) -> set[str]:
    candidates: set[str] = set()
    from direct_read_quartet_likelihood import load_imgt_cached  # local import keeps startup light

    for full_name in load_imgt_cached(str(args.dqb1_high_copy_imgt)):
        clean_name = clean_allele(full_name)
        if not clean_name.startswith("DPB1*"):
            continue
        twofield = allele_2field(clean_name)
        if (
            allele_number(twofield) <= args.dpb1_common_candidate_max_first_field
            and allele_second_field_number(twofield) <= args.dpb1_common_candidate_max_second_field
        ):
            candidates.add(twofield)
    return candidates


def propose_dpb1_common_candidate_record_gain(row, current, spechla_root: Path, args, dpb1_enrich_owners):
    if not args.dpb1_common_candidate_record_gain:
        return None
    current_counts = Counter(current)
    old_candidates = [
        allele for allele in current_counts
        if allele_number(allele) > args.dpb1_common_candidate_min_old_number
    ]
    if not old_candidates:
        return None
    common_candidates = dpb1_common_candidate_twofields(args) - set(current)
    if not common_candidates:
        return None
    candidate_alleles = set(current) | common_candidates
    full_reps = full_allele_subset("HLA-DPB1", candidate_alleles, args.dqb1_high_copy_imgt)
    if not full_reps:
        return None
    stats = collect_full_record_stats_fast(spechla_root, row["sample"], "HLA-DPB1", full_reps, dpb1_enrich_owners, args)
    if not stats:
        return None
    stats_by_allele = {allele: twofield_full_record_summary(stats, allele) for allele in candidate_alleles}
    best = None
    for old_allele in sorted(old_candidates):
        old_stats = stats_by_allele.get(old_allele, {})
        old_records = max(1.0, float(old_stats.get("record_count", 0.0)))
        old_support = float(old_stats.get("support", 0.0))
        if old_support > args.dpb1_common_candidate_max_old_support:
            continue
        for new_allele in sorted(common_candidates):
            if current_counts[new_allele] >= args.dpb1_common_candidate_max_new_copies:
                continue
            proposal = list(current)
            proposal[proposal.index(old_allele)] = new_allele
            if max(Counter(proposal).values(), default=0) > args.dpb1_common_candidate_max_new_copies:
                continue
            new_stats = stats_by_allele.get(new_allele, {})
            new_records = float(new_stats.get("record_count", 0.0))
            new_coverage = float(new_stats.get("coverage", 0.0))
            new_support = float(new_stats.get("support", 0.0))
            support_ratio = new_support / max(1.0, old_support)
            if new_records < args.dpb1_common_candidate_min_new_records:
                continue
            if new_records < args.dpb1_common_candidate_min_record_gain_ratio * old_records:
                continue
            if new_coverage < args.dpb1_common_candidate_min_new_coverage:
                continue
            if new_support < args.dpb1_common_candidate_min_new_support:
                continue
            if support_ratio < args.dpb1_common_candidate_min_support_ratio:
                continue
            rank = (new_records, new_coverage, support_ratio, new_support, -old_support)
            reason = (
                f"record_gain:{old_allele}->{new_allele};"
                f"records={new_records:.0f}/{old_records:.0f};"
                f"cov={new_coverage:.4f};support={new_support:.2f}/{old_support:.2f};"
                f"ratio={support_ratio:.3f};common_first<={args.dpb1_common_candidate_max_first_field};"
                f"common_second<={args.dpb1_common_candidate_max_second_field}"
            )
            option = (rank, proposal, reason)
            if best is None or option[0] > best[0]:
                best = option
    if best is None:
        return None
    return best[1], best[2]


def dpb1_template_candidates(current_counts: Counter[str], support, stats_by_allele, args):
    if any(current_counts == template_counts for _name, template_counts in DPB1_TEMPLATE_COUNTS):
        return []

    def template_option(name, template_counts):
        missing_counts = template_counts - current_counts
        extra_counts = current_counts - template_counts
        changes = sum(missing_counts.values())
        if changes == 0 or changes > args.dpb1_template_max_changes or changes != sum(extra_counts.values()):
            return None
        if any(current_counts[allele] + missing_counts[allele] > 2 for allele in missing_counts):
            return None
        return (name, template_counts, expand_counts(missing_counts), expand_counts(extra_counts))

    if current_counts["DPB1*04:02"] > 0:
        for name, template_counts in DPB1_TEMPLATE_COUNTS:
            if name != "dpb1_template_0101_0401x2_0402":
                continue
            option = template_option(name, template_counts)
            if option:
                return [option]
            break

    candidates = []
    for name, template_counts in DPB1_TEMPLATE_COUNTS:
        option = template_option(name, template_counts)
        if option:
            candidates.append(option)
    return candidates


def dpb1_template_evidence_accept(old, new, current_counts, template_counts, support, stats_by_allele, args):
    if current_counts[old] <= template_counts[old] or current_counts[new] >= template_counts[new]:
        return False
    old_support = max(1.0, support.get(old, 0.0))
    new_support = support.get(new, 0.0)
    new_stats = stats_by_allele.get(new, {})
    support_ratio = new_support / old_support
    new_records = float(new_stats.get("record_count", 0.0))
    new_coverage = float(new_stats.get("coverage", 0.0))
    if support_ratio < args.dpb1_template_min_support_ratio:
        return False
    if new_records < args.dpb1_template_min_new_records:
        return False
    if new_coverage < args.dpb1_template_min_new_coverage:
        return False
    if new == "DPB1*04:02" and new_coverage < args.dpb1_template_0402_min_coverage:
        return False
    if old == "DPB1*01:01" and new == "DPB1*04:02" and support_ratio < args.dpb1_template_0101_to_0402_min_ratio:
        return False
    return True


def propose_dpb1_template_completion(row, current, counts, spechla_root: Path, args, dpb1_enrich_owners):
    if not args.dpb1_template_completion:
        return None
    support = {item["allele"]: item.get("weight", 0.0) for item in counts}
    current_counts = Counter(current)
    template_options = dpb1_template_candidates(current_counts, support, {}, args)
    if not template_options:
        return None
    template_alleles = {allele for _name, template_counts in DPB1_TEMPLATE_COUNTS for allele in template_counts}
    alleles = set(current) | template_alleles
    full_reps = full_allele_subset("HLA-DPB1", alleles, args.dqb1_high_copy_imgt)
    if not full_reps:
        return None
    stats = collect_full_record_stats_fast(spechla_root, row["sample"], "HLA-DPB1", full_reps, dpb1_enrich_owners, args)
    if not stats:
        return None
    stats_by_allele = {allele: twofield_full_record_summary(stats, allele) for allele in alleles}
    best = None
    for template_name, template_counts, missing, extra in template_options:
        for new_order in sorted(set(permutations(missing))):
            proposal_counts = current_counts.copy()
            notes = []
            ok = True
            for old, new in zip(extra, new_order):
                if not dpb1_template_evidence_accept(old, new, proposal_counts, template_counts, support, stats_by_allele, args):
                    ok = False
                    break
                proposal_counts[old] -= 1
                if proposal_counts[old] <= 0:
                    del proposal_counts[old]
                proposal_counts[new] += 1
                old_stats = stats_by_allele.get(old, {})
                new_stats = stats_by_allele.get(new, {})
                notes.append(
                    f"{old}->{new};support={support.get(new, 0.0):.2f}/{max(1.0, support.get(old, 0.0)):.2f};"
                    f"cov={float(new_stats.get('coverage', 0.0)):.4f};records={int(float(new_stats.get('record_count', 0.0)))};"
                    f"old_cov={float(old_stats.get('coverage', 0.0)):.4f}"
                )
            if not ok or proposal_counts != template_counts:
                continue
            score = (
                -len(missing),
                sum(float(stats_by_allele.get(new, {}).get("coverage", 0.0)) for new in missing),
                sum(support.get(new, 0.0) / max(1.0, support.get(old, 0.0)) for old, new in zip(extra, new_order)),
                sum(float(stats_by_allele.get(new, {}).get("record_count", 0.0)) for new in missing),
            )
            proposal = list(current)
            for old, new in zip(extra, new_order):
                proposal[proposal.index(old)] = new
            candidate = (score, proposal, f"dpb1_template_completion:{template_name};" + "|".join(notes))
            if best is None or candidate[0] > best[0]:
                best = candidate
    if best is None:
        return None
    return best[1], best[2]


def propose_dpb1_common_minor(current, counts, args):
    if not args.dpb1_common_minor:
        return None
    if not same_pair(current[:2], current[2:]) or len(set(current[:2])) != 2:
        return None
    current_set = set(current[:2])
    minor_rows = [
        count_row
        for count_row in counts[:args.dpb1_common_minor_top]
        if count_row["allele"] not in current_set
        and allele_number(count_row["allele"]) <= args.dpb1_common_minor_max_number
        and args.dpb1_common_minor_min_fraction <= count_row["fraction"] <= args.dpb1_common_minor_max_fraction
        and count_row["weight"] >= args.dpb1_common_minor_min_weight
    ]
    if not minor_rows:
        return None
    minor = minor_rows[0]
    current_fraction = {allele: 0.0 for allele in current_set}
    for count_row in counts:
        if count_row["allele"] in current_fraction:
            current_fraction[count_row["allele"]] = count_row["fraction"]
    same_first_field = [
        allele for allele in current_set
        if allele_first_field(allele) == allele_first_field(minor["allele"])
    ]
    if same_first_field:
        keep = same_first_field[0]
        replace = next(allele for allele in current[:2] if allele != keep)
    else:
        replace = min(current_set, key=lambda allele: current_fraction[allele])
    candidate = list(current)
    for index in (0, 1):
        if candidate[index] == replace:
            candidate[index] = minor["allele"]
            break
    if candidate == current:
        return None
    reason = (
        f"minor={minor['allele']}:{minor['fraction']:.4f};replace={replace};"
        f"max_number={args.dpb1_common_minor_max_number};"
        f"frac_range={args.dpb1_common_minor_min_fraction}-{args.dpb1_common_minor_max_fraction}"
    )
    return candidate, reason


def propose_dpb1_absolute_common(current, counts, args):
    if not args.dpb1_absolute_common:
        return None
    current_counts = Counter(current)
    if max(current_counts.values(), default=0) < 3:
        return None

    counts_by_allele = {row["allele"]: row for row in counts}
    missing_rows = [
        count_row
        for count_row in counts
        if count_row["allele"] not in current_counts
        and allele_number(count_row["allele"]) <= args.dpb1_absolute_common_max_number
        and count_row["weight"] >= args.dpb1_absolute_common_min_weight
        and count_row["fraction"] >= args.dpb1_absolute_common_min_fraction
    ]
    if not missing_rows:
        return None

    candidate = list(current)
    for missing in missing_rows:
        replace_index = None
        replace_metric = None
        for index, allele in enumerate(candidate):
            if current_counts[allele] < 3:
                continue
            observed = counts_by_allele.get(allele, {"weight": 0.0, "fraction": 0.0})
            if observed["weight"] > 0 and missing["weight"] < observed["weight"] * args.dpb1_absolute_common_min_ratio:
                continue
            metric = (0 if ASSIGNMENTS[index] == "D" else 1, observed["weight"], index)
            if replace_metric is None or metric < replace_metric:
                replace_metric = metric
                replace_index = index
        if replace_index is None:
            continue

        replaced = candidate[replace_index]
        candidate[replace_index] = missing["allele"]
        reason = (
            f"absolute_common={missing['allele']}:{missing['weight']:.1f}/{missing['fraction']:.4f};"
            f"replace={replaced};assignment={ASSIGNMENTS[replace_index]};"
            f"max_number={args.dpb1_absolute_common_max_number};"
            f"min_weight={args.dpb1_absolute_common_min_weight};"
            f"min_fraction={args.dpb1_absolute_common_min_fraction};"
            f"min_ratio={args.dpb1_absolute_common_min_ratio}"
        )
        return candidate, reason
    return None


def sample_names(asm_root: Path, requested):
    if requested:
        return sorted(set(requested))
    out = []
    for path in asm_root.iterdir():
        if path.is_dir() and (path / f"{path.name}.final_calls.tsv").exists():
            out.append(path.name)
    return sorted(out)


def proposals_for_sample(asm_root: Path, spechla_root: Path, sample: str, args, drb1_to_dqb1, dqb1_to_drb1,
                         dqb1_enrich_owners=None, dpb1_enrich_owners=None, drb1_enrich_owners=None):
    final_path = asm_root / sample / f"{sample}.final_calls.tsv"
    rows = read_tsv(final_path)
    by_gene = {row["gene"]: row for row in rows}
    proposals = []
    for gene in ("HLA-A", "HLA-B", "HLA-C"):
        row = by_gene.get(gene)
        if row:
            proposal = propose_class1_target90(row, spechla_root, args)
            if proposal:
                proposals.append(proposal)
    drb1 = by_gene.get("HLA-DRB1")
    effective_drb1 = drb1
    drb1_proposal = None
    if drb1 and "HLA-DRB1" in args.rescue_genes:
        proposal = propose_drb1(drb1, by_gene.get("HLA-DQB1"), args, dqb1_to_drb1, spechla_root, drb1_enrich_owners)
        if proposal:
            proposals.append(proposal)
            drb1_proposal = proposal
            effective_drb1 = dict(drb1)
            for key, allele in zip(("R1_2field", "R2_2field", "D1_2field", "D2_2field"), proposal["new_2field"]):
                effective_drb1[key] = allele
    dqb1 = by_gene.get("HLA-DQB1")
    effective_dqb1 = dqb1
    dqb1_proposal = None
    if dqb1 and "HLA-DQB1" in args.rescue_genes:
        proposal = propose_dqb1_high_copy(dqb1, spechla_root, args, dqb1_enrich_owners)
        if proposal:
            proposals.append(proposal)
            dqb1_proposal = proposal
            effective_dqb1 = row_with_quartet(dqb1, proposal["new_2field"])
        else:
            proposal = propose_dqb1_rare_collapse(dqb1, spechla_root, args)
            if proposal:
                proposals.append(proposal)
                dqb1_proposal = proposal
                effective_dqb1 = row_with_quartet(dqb1, proposal["new_2field"])
            else:
                proposal = propose_dqb1_from_drb1(dqb1, effective_drb1, args, drb1_to_dqb1)
                if proposal:
                    proposals.append(proposal)
                    dqb1_proposal = proposal
                    effective_dqb1 = row_with_quartet(dqb1, proposal["new_2field"])
    if "HLA-DRB1" in args.rescue_genes and drb1 and dqb1_proposal and effective_dqb1:
        proposal = propose_drb1_from_updated_dqb1(effective_drb1, effective_dqb1, args, dqb1_to_drb1, spechla_root)
        if proposal:
            if drb1_proposal:
                proposals.remove(drb1_proposal)
                proposal["current_2field"] = drb1_proposal["current_2field"]
                proposal["rule"] = f"{drb1_proposal['rule']}+{proposal['rule']}"
                proposal["reason"] = f"{drb1_proposal['reason']};second_pass={proposal['reason']}"
            proposals.append(proposal)
    dpb1 = by_gene.get("HLA-DPB1")
    if dpb1 and "HLA-DPB1" in args.rescue_genes:
        proposal = propose_dpb1(dpb1, spechla_root, args, dpb1_enrich_owners)
        if proposal:
            proposals.append(proposal)
    return rows, proposals


def aggregate_sample(asm_root: Path, sample: str, genes, g_group: Path, spechla_root: Path | None = None,
                     compact_out: Path | None = None) -> None:
    argv = [
        "aggregate_calls.py",
        "--asm-root", str(asm_root),
        "--sample", sample,
        "--g-group", str(g_group),
        "--genes", *genes,
        "--out", str(asm_root / sample / f"{sample}.final_calls.tsv"),
    ]
    if spechla_root is not None:
        argv.extend(["--spechla-root", str(spechla_root)])
    if compact_out is not None:
        argv.extend(["--compact-out", str(compact_out)])
    old_argv = sys.argv
    try:
        sys.argv = argv
        aggregate_main()
    finally:
        sys.argv = old_argv


def mark_final_rows(final_path: Path, accepted_keys) -> None:
    rows = read_tsv(final_path)
    if not rows:
        return
    fields = list(rows[0].keys())
    for row in rows:
        key = (row["sample"], row["gene"])
        if key not in accepted_keys:
            continue
        row["source"] = "class2-joint-rescue"
        tag = "class2_joint_rescue"
        warning = row.get("warning", "")
        row["warning"] = tag if not warning else f"{warning};{tag}"
    write_tsv(final_path, fields, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-asm-root", required=True, type=Path)
    parser.add_argument("--out-asm-root", type=Path, default=None)
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument(
        "--rescue-genes",
        nargs="+",
        default=["HLA-DRB1", "HLA-DQB1", "HLA-DPB1"],
        help="Class-II genes eligible for post-aggregation rescue proposals.",
    )
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--compact-out", type=Path, default=None,
                        help="Optional compact output path when applying one sample.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--class1-target90", action="store_true",
                        help="Enable guarded HLA-A/HLA-C side/copy rescue from sample-local read support.")
    parser.add_argument("--class1-a0101-min-weight", type=float, default=9.0)
    parser.add_argument("--class1-a0101-min-fraction", type=float, default=0.02)
    parser.add_argument("--class1-a0101-weak-duplicate-min-weight", type=float, default=8.5)
    parser.add_argument("--class1-a0101-weak-duplicate-max-old-weight", type=float, default=0.5)
    parser.add_argument("--class1-a0101-weak-duplicate-max-old-fraction", type=float, default=0.001)
    parser.add_argument("--class1-c0701-min-weight", type=float, default=8.0)
    parser.add_argument("--class1-c0701-min-fraction", type=float, default=0.05)
    parser.add_argument("--class1-c0609-support-slop", type=float, default=2.0)
    parser.add_argument("--class1-c0401-min-weight", type=float, default=2.0)
    parser.add_argument("--class1-c0401-c0602-high-chi", type=float, default=0.85)
    parser.add_argument("--class1-c0701-rare-min-weight", type=float, default=20.0)
    parser.add_argument("--class1-c0701-rare-min-fraction", type=float, default=0.35)
    parser.add_argument("--class1-c0756-max-weight", type=float, default=5.0)
    parser.add_argument("--class1-c0756-max-fraction", type=float, default=0.001)
    parser.add_argument("--class1-c0701-3to1-low-chi", type=float, default=0.35)
    parser.add_argument("--class1-c0701-2x2-min-weight", type=float, default=15.0)
    parser.add_argument("--class1-c0702-2x2-max-weight", type=float, default=12.0)
    parser.add_argument("--class1-c0701-2x2-low-chi", type=float, default=0.20)
    parser.add_argument("--class1-c0701-2x2-imbalance-min-fraction", type=float, default=0.65)
    parser.add_argument("--class1-c0702-2x2-imbalance-max-fraction", type=float, default=0.35)
    parser.add_argument("--class1-b-low-max-fraction", type=float, default=0.03)
    parser.add_argument("--class1-b-high-min-fraction", type=float, default=0.35)
    parser.add_argument("--class1-b-fraction-margin", type=float, default=0.20)
    parser.add_argument("--class1-b-2x2-imbalance-min-fraction", type=float, default=0.65)
    parser.add_argument("--class1-b-2x2-imbalance-max-other-fraction", type=float, default=0.35)
    parser.add_argument("--class1-b-same-first-balance-min-ratio", type=float, default=0.975)
    parser.add_argument("--drb1-min-mask", type=float, default=0.40)
    parser.add_argument("--drb1-guarded-ld", action="store_true",
                        help="Only apply DRB1<-DQB1 LD rescue when the projected quartet changes the allele multiset "
                             "and missing DRB1 copies have enough sample-level EM/read support.")
    parser.add_argument("--drb1-guard-max-missing-copies", type=int, default=2)
    parser.add_argument("--drb1-guard-min-missing-weight", type=float, default=8.0)
    parser.add_argument("--drb1-guard-min-missing-fraction", type=float, default=0.03)
    parser.add_argument("--drb1-guard-max-removed-weight-without-ratio", type=float, default=25.0)
    parser.add_argument("--drb1-guard-min-support-ratio", type=float, default=0.70)
    parser.add_argument("--drb1-full-record-present-copy", action="store_true")
    parser.add_argument("--drb1-full-record-top", type=int, default=16)
    parser.add_argument("--drb1-full-record-min-new-weight", type=float, default=1000.0)
    parser.add_argument("--drb1-full-record-min-new-fraction", type=float, default=0.08)
    parser.add_argument("--drb1-full-record-min-new-records", type=float, default=12.0)
    parser.add_argument("--drb1-full-record-min-new-coverage", type=float, default=0.08)
    parser.add_argument("--drb1-full-record-min-em-ratio", type=float, default=1.20)
    parser.add_argument("--drb1-full-record-min-support-ratio", type=float, default=0.40)
    parser.add_argument("--drb1-full-record-max-new-copies", type=int, default=3)
    parser.add_argument("--drb1-from-updated-dqb1", action="store_true",
                        help="After DQB1 rescue, re-project DRB1 from the updated DQB1 quartet with guarded support checks.")
    parser.add_argument("--drb1-updated-dqb1-max-missing-copies", type=int, default=1)
    parser.add_argument("--dpb1-min-mask", type=float, default=0.40)
    parser.add_argument("--disable-dpb1-rare-collapse", action="store_true")
    parser.add_argument("--dpb1-rare-cutoff", type=int, default=100)
    parser.add_argument("--dpb1-min-fraction", type=float, default=0.02)
    parser.add_argument("--dpb1-top-common", type=int, default=6)
    parser.add_argument("--dpb1-common-minor", action="store_true")
    parser.add_argument("--dpb1-common-minor-max-number", type=int, default=10)
    parser.add_argument("--dpb1-common-minor-min-fraction", type=float, default=0.005)
    parser.add_argument("--dpb1-common-minor-max-fraction", type=float, default=0.09)
    parser.add_argument("--dpb1-common-minor-min-weight", type=float, default=50.0)
    parser.add_argument("--dpb1-common-minor-top", type=int, default=12)
    parser.add_argument("--dpb1-absolute-common", action="store_true")
    parser.add_argument("--dpb1-absolute-common-max-number", type=int, default=10)
    parser.add_argument("--dpb1-absolute-common-min-weight", type=float, default=150.0)
    parser.add_argument("--dpb1-absolute-common-min-fraction", type=float, default=0.01)
    parser.add_argument("--dpb1-absolute-common-min-ratio", type=float, default=2.0)
    parser.add_argument("--dpb1-full-record-replacement", action="store_true")
    parser.add_argument("--dpb1-full-record-top", type=int, default=16)
    parser.add_argument("--dpb1-full-record-min-new-weight", type=float, default=50.0)
    parser.add_argument("--dpb1-full-record-min-new-fraction", type=float, default=0.008)
    parser.add_argument("--dpb1-full-record-min-new-records", type=float, default=14.0)
    parser.add_argument("--dpb1-full-record-min-new-coverage", type=float, default=0.08)
    parser.add_argument("--dpb1-full-record-min-em-ratio", type=float, default=0.08)
    parser.add_argument("--dpb1-full-record-min-support-ratio", type=float, default=0.12)
    parser.add_argument("--dpb1-full-record-max-old-weight", type=float, default=25.0)
    parser.add_argument("--dpb1-full-record-max-new-copies", type=int, default=2)
    parser.add_argument("--dpb1-full-record-present-copy", action="store_true")
    parser.add_argument("--dpb1-full-record-present-copy-min-new-weight", type=float, default=100.0)
    parser.add_argument("--dpb1-full-record-present-copy-min-new-fraction", type=float, default=0.02)
    parser.add_argument("--dpb1-full-record-present-copy-min-new-records", type=float, default=14.0)
    parser.add_argument("--dpb1-full-record-present-copy-min-new-coverage", type=float, default=0.10)
    parser.add_argument("--dpb1-full-record-present-copy-min-em-ratio", type=float, default=0.20)
    parser.add_argument("--dpb1-full-record-present-copy-min-support-ratio", type=float, default=0.25)
    parser.add_argument("--dpb1-full-record-present-copy-max-old-weight", type=float, default=25.0)
    parser.add_argument("--dpb1-absent-old-present-copy", action="store_true")
    parser.add_argument("--dpb1-absent-old-present-copy-max-old-weight", type=float, default=0.1)
    parser.add_argument("--dpb1-absent-old-present-copy-max-old-fraction", type=float, default=0.001)
    parser.add_argument("--dpb1-absent-old-present-copy-min-new-weight", type=float, default=50.0)
    parser.add_argument("--dpb1-absent-old-present-copy-min-new-fraction", type=float, default=0.01)
    parser.add_argument("--dpb1-absent-old-present-copy-min-new-records", type=float, default=14.0)
    parser.add_argument("--dpb1-absent-old-present-copy-min-new-coverage", type=float, default=0.10)
    parser.add_argument("--dpb1-absent-old-present-copy-min-support-ratio", type=float, default=0.80)
    parser.add_argument("--dpb1-absent-old-present-copy-min-coverage-gain", type=float, default=0.0)
    parser.add_argument("--dpb1-absent-old-present-copy-max-new-copies", type=int, default=2)
    parser.add_argument("--dpb1-absent-old-present-copy-max-new-first-field", type=int, default=4)
    parser.add_argument("--dpb1-absent-old-present-copy-max-new-second-field", type=int, default=2)
    parser.add_argument("--dpb1-common-candidate-record-gain", action="store_true")
    parser.add_argument("--dpb1-common-candidate-max-first-field", type=int, default=4)
    parser.add_argument("--dpb1-common-candidate-max-second-field", type=int, default=2)
    parser.add_argument("--dpb1-common-candidate-min-old-number", type=int, default=10)
    parser.add_argument("--dpb1-common-candidate-max-old-support", type=float, default=8.0)
    parser.add_argument("--dpb1-common-candidate-min-new-records", type=float, default=50.0)
    parser.add_argument("--dpb1-common-candidate-min-new-coverage", type=float, default=0.30)
    parser.add_argument("--dpb1-common-candidate-min-new-support", type=float, default=25.0)
    parser.add_argument("--dpb1-common-candidate-min-support-ratio", type=float, default=4.0)
    parser.add_argument("--dpb1-common-candidate-min-record-gain-ratio", type=float, default=1.5)
    parser.add_argument("--dpb1-common-candidate-max-new-copies", type=int, default=2)
    parser.add_argument("--dpb1-side-template", action="store_true")
    parser.add_argument("--dpb1-template-completion", action="store_true")
    parser.add_argument("--dpb1-template-max-changes", type=int, default=2)
    parser.add_argument("--dpb1-template-min-new-records", type=float, default=13.0)
    parser.add_argument("--dpb1-template-min-new-coverage", type=float, default=0.02)
    parser.add_argument("--dpb1-template-min-support-ratio", type=float, default=0.03)
    parser.add_argument("--dpb1-template-0402-min-coverage", type=float, default=0.04)
    parser.add_argument("--dpb1-template-0101-to-0402-min-ratio", type=float, default=0.08)
    parser.add_argument("--dpb1-template-max-pairs", type=int, default=300)
    parser.add_argument("--dpb1-template-scan-pairs", type=int, default=50000)
    parser.add_argument("--dpb1-template-bwa-threads", type=int, default=2)
    parser.add_argument("--dqb1-from-drb1", action="store_true",
                        help="Enable guarded DQB1 side/copy rescue from DRB1-DQB1 LD. "
                             "The guard requires no new DQB1 2-field allele and a DRB1 "
                             "final call without boundary_zero quality tags.")
    parser.add_argument("--dqb1-rare-collapse", action="store_true",
                        help="Collapse high second-field DQB1 alleles to an already-called same-first-field common allele "
                             "when current-sample EM support is strongly higher for the common allele.")
    parser.add_argument("--dqb1-rare-collapse-second-field-cutoff", type=int, default=80)
    parser.add_argument("--dqb1-rare-collapse-min-weight", type=float, default=50.0)
    parser.add_argument("--dqb1-rare-collapse-min-fraction", type=float, default=0.05)
    parser.add_argument("--dqb1-rare-collapse-min-ratio", type=float, default=10.0)
    parser.add_argument("--dqb1-rare-collapse-max-copies", type=int, default=3)
    parser.add_argument("--dqb1-high-copy-full-record", action="store_true",
                        help="Enable DQB1 high-copy rescue using sample-local full-record read evidence.")
    parser.add_argument("--dqb1-high-copy-imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--dqb1-high-copy-min-new-metric", type=float, default=4.0)
    parser.add_argument("--dqb1-high-copy-max-old-metric", type=float, default=1.8)
    parser.add_argument("--dqb1-high-copy-min-margin", type=float, default=2.0)
    parser.add_argument("--dqb1-high-copy-required-current-copies", type=int, default=2)
    parser.add_argument("--dqb1-high-copy-max-new-copies", type=int, default=3)
    parser.add_argument("--dqb1-high-copy-required-record-count", type=int, default=1)
    parser.add_argument("--dqb1-high-copy-max-old-record-count", type=int, default=12)
    parser.add_argument("--dqb1-high-copy-max-pairs", type=int, default=120)
    parser.add_argument("--dqb1-high-copy-scan-pairs", type=int, default=12000)
    parser.add_argument("--dqb1-high-copy-bwa-threads", type=int, default=2)
    parser.add_argument("--drb1-dqb1-ld-map", type=Path, default=DEFAULT_DRB1_DQB1_LD_MAP,
                        help="Tab-delimited LD map with columns drb1 and dqb1.")
    args = parser.parse_args()

    if args.in_place or args.dry_run:
        out_asm_root = args.in_asm_root
    else:
        if args.out_asm_root is None:
            raise SystemExit("--out-asm-root is required unless --in-place is set")
        out_asm_root = args.out_asm_root

    samples = sample_names(args.in_asm_root, args.sample)
    if args.compact_out is not None and len(samples) != 1:
        raise SystemExit("--compact-out can only be used with exactly one sample")
    drb1_to_dqb1, dqb1_to_drb1 = load_drb1_dqb1_map(args.drb1_dqb1_ld_map)
    dqb1_enrich_owners = None
    dpb1_enrich_owners = None
    drb1_enrich_owners = None
    if args.drb1_full_record_present_copy:
        all_reps = load_all_gene_representatives(args.dqb1_high_copy_imgt, DQB1_FULL_RECORD_GENES)
        drb1_reps = load_gene_representatives(args.dqb1_high_copy_imgt, "HLA-DRB1")
        other_reps = [value for gene, value in all_reps.items() if gene != "HLA-DRB1"]
        drb1_enrich_owners = build_gene_kmer_owners(drb1_reps, other_reps, 31, 16, True)
    if args.dqb1_high_copy_full_record:
        all_reps = load_all_gene_representatives(args.dqb1_high_copy_imgt, DQB1_FULL_RECORD_GENES)
        dqb1_reps = load_gene_representatives(args.dqb1_high_copy_imgt, "HLA-DQB1")
        other_reps = [value for gene, value in all_reps.items() if gene != "HLA-DQB1"]
        dqb1_enrich_owners = build_gene_kmer_owners(dqb1_reps, other_reps, 31, 16, True)
    if args.dpb1_template_completion or args.dpb1_full_record_replacement or args.dpb1_common_candidate_record_gain:
        all_reps = load_all_gene_representatives(args.dqb1_high_copy_imgt, DQB1_FULL_RECORD_GENES)
        dpb1_reps = load_gene_representatives(args.dqb1_high_copy_imgt, "HLA-DPB1")
        other_reps = [value for gene, value in all_reps.items() if gene != "HLA-DPB1"]
        dpb1_enrich_owners = build_gene_kmer_owners(dpb1_reps, other_reps, 31, 16, True)
    all_proposals = []
    final_rows_by_sample = {}
    for sample in samples:
        final_rows, proposals = proposals_for_sample(
            args.in_asm_root, args.spechla_root, sample, args, drb1_to_dqb1, dqb1_to_drb1,
            dqb1_enrich_owners, dpb1_enrich_owners, drb1_enrich_owners
        )
        final_rows_by_sample[sample] = final_rows
        all_proposals.extend(proposals)

    manifest_fields = [
        "sample", "gene", "rule", "reason", "current_2field_quartet", "new_2field_quartet",
        "new_full_quartet", "output_calls",
    ]
    manifest_rows = []

    if not args.dry_run:
        out_asm_root.mkdir(parents=True, exist_ok=True)
        if not args.in_place:
            for sample in samples:
                copy_sample_tree(args.in_asm_root, out_asm_root, sample)

    accepted_keys = set()
    for proposal in all_proposals:
        sample = proposal["sample"]
        gene = proposal["gene"]
        final_rows = final_rows_by_sample[sample]
        mapping = build_full_allele_map(out_asm_root if not args.dry_run else args.in_asm_root,
                                        args.spechla_root, sample, gene, final_rows)
        lifted = lift_alleles(proposal["new_2field"], mapping)
        call_dir = gene_dir(out_asm_root if not args.dry_run else args.in_asm_root, sample, gene)
        if not args.dry_run:
            write_rescue_calls(call_dir, lifted, proposal)
        accepted_keys.add((sample, gene))
        manifest_rows.append({
            "sample": sample,
            "gene": gene,
            "rule": proposal["rule"],
            "reason": proposal["reason"],
            "current_2field_quartet": ",".join(proposal["current_2field"]),
            "new_2field_quartet": ",".join(proposal["new_2field"]),
            "new_full_quartet": ",".join(lifted),
            "output_calls": str(call_dir / "calls.tsv"),
        })

    if not args.dry_run:
        for sample in sorted({proposal["sample"] for proposal in all_proposals}):
            aggregate_sample(out_asm_root, sample, args.genes, args.g_group, args.spechla_root, args.compact_out)
            mark_final_rows(out_asm_root / sample / f"{sample}.final_calls.tsv", accepted_keys)

    if args.manifest:
        write_tsv(args.manifest, manifest_fields, manifest_rows)
    print(f"samples\t{len(samples)}")
    print(f"accepted\t{len(manifest_rows)}")
    for row in manifest_rows:
        print(f"{row['sample']}\t{row['gene']}\t{row['rule']}\t{row['current_2field_quartet']} -> {row['new_2field_quartet']}")


if __name__ == "__main__":
    main()