#!/usr/bin/env python3
"""Apply the frozen source-agnostic quartet optimizers to one sample."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
import tempfile
from collections import Counter, defaultdict, deque
from itertools import permutations
from pathlib import Path
from types import SimpleNamespace

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from aggregate_calls import DEFAULT_GENES, allele_2field, main as aggregate_main  # noqa: E402
from diagnostics.direct_read_quartet_likelihood import (  # noqa: E402
    DEFAULT_IMGT,
    build_informative_kmers,
    load_candidate_sequences,
)
from diagnostics.offline_class_i_hybrid_quartet import (  # noqa: E402
    compare_quartets,
    iter_fastq_pairs,
    pair_evidence,
    proposal_is_supported,
)
from diagnostics.offline_joint_quartet_posterior import (  # noqa: E402
    call_quartet,
    read_counts,
    read_major_fraction_prior,
)

CLASS_I_GENES = ("HLA-A", "HLA-B", "HLA-C")
CLASS_II_GENES = ("HLA-DRB1", "HLA-DPB1", "HLA-DQB1")
ASSIGNMENTS = ("R", "R", "D", "D")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return [], []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def gene_dir(asm_root: Path, sample: str, gene: str) -> Path:
    return asm_root / sample / gene.lower() / gene


def baseline_rows(call_dir: Path) -> tuple[list[str], list[dict[str, str]]]:
    backup = call_dir / "calls.quartet_optimization_input.tsv"
    return read_tsv(backup if backup.exists() else call_dir / "calls.tsv")


def baseline_quartet(rows: list[dict[str, str]]) -> tuple[str, ...]:
    ordered = sorted(
        rows,
        key=lambda row: int(row.get("global_hap", "0")) if row.get("global_hap", "").isdigit() else 0,
    )
    return tuple(allele_2field(row.get("allele", "")) for row in ordered)


def read_chi_path(spechla_root: Path, sample: str) -> Path | None:
    sample_dir = spechla_root / sample
    for path in (sample_dir / f"{sample}.chi_pooled.txt", sample_dir / f"{sample}.chimerism.txt"):
        if path.exists():
            return path
    return None


def full_allele_map(
    asm_root: Path,
    spechla_root: Path,
    sample: str,
    gene: str,
) -> dict[str, deque[str]]:
    mapping: dict[str, deque[str]] = defaultdict(deque)
    call_dir = gene_dir(asm_root, sample, gene)
    paths = [
        call_dir / "calls.tsv",
        call_dir / "calls.baseline.tsv",
        call_dir / "calls.class2_joint_input.tsv",
        call_dir / "calls.pre_private_rescue.tsv",
        call_dir / "calls.quartet_optimization_input.tsv",
        spechla_root / sample / "em_refine" / f"{gene}.calls.tsv",
    ]
    for path in paths:
        _fields, rows = read_tsv(path)
        for row in rows:
            allele = row.get("allele") or row.get("allele_2field") or row.get("call")
            if allele and allele != "NA":
                mapping[allele_2field(allele)].append(allele)
    return mapping


def lift_quartet(quartet: tuple[str, ...], mapping: dict[str, deque[str]]) -> list[str] | None:
    used: Counter[str] = Counter()
    lifted = []
    for allele in quartet:
        options = mapping.get(allele)
        if not options:
            return None
        index = used[allele]
        used[allele] += 1
        lifted.append(options[index] if index < len(options) else options[-1])
    return lifted


def order_pair(pair: tuple[str, str], baseline: tuple[str, str]) -> tuple[str, str]:
    return min(set(permutations(pair)), key=lambda item: (sum(a != b for a, b in zip(item, baseline)), item))


def slot_quartet(
    major_group: tuple[str, str],
    minor_group: tuple[str, str],
    chi_r: float,
    baseline: tuple[str, ...],
) -> tuple[str, ...]:
    recipient, donor = (major_group, minor_group) if chi_r >= 0.5 else (minor_group, major_group)
    return order_pair(recipient, baseline[:2]) + order_pair(donor, baseline[2:])


def class_i_read_gate(
    args: argparse.Namespace,
    gene: str,
    baseline: tuple[str, ...],
    proposal: tuple[str, ...],
) -> tuple[bool, dict[str, object]]:
    if tuple(sorted(baseline)) == tuple(sorted(proposal)):
        return False, {"decision": "same", "reason": "proposal_equals_baseline"}
    short = gene.removeprefix("HLA-")
    fq1 = args.spechla_root / args.sample / f"{short}.R1.fq.gz"
    fq2 = args.spechla_root / args.sample / f"{short}.R2.fq.gz"
    if not fq1.exists() or not fq2.exists():
        return False, {"decision": "fallback", "reason": "missing_gene_fastq"}
    candidates = sorted(set(baseline) | set(proposal))
    sequences = load_candidate_sequences(args.imgt, gene, candidates)
    missing = sorted(allele for allele, records in sequences.items() if not records)
    if missing:
        return False, {"decision": "fallback", "reason": f"missing_sequences={','.join(missing)}"}
    owners = build_informative_kmers(sequences, args.k, args.max_full_alleles, args.max_owner_fraction)
    evidence_rows = []
    private_counts: Counter[str] = Counter()
    total_pairs = 0
    for sequence1, quality1, sequence2, quality2 in iter_fastq_pairs(fq1, fq2):
        total_pairs += 1
        evidence, private = pair_evidence(
            sequence1, quality1, sequence2, quality2, owners, args.k, args.concordance_bonus
        )
        if sum(evidence.values()) < args.min_pair_evidence:
            continue
        evidence_rows.append(evidence)
        private_counts.update(private)
    if not evidence_rows:
        return False, {"decision": "fallback", "reason": "no_informative_pairs", "total_pairs": total_pairs}
    comparison = compare_quartets(evidence_rows, tuple(sorted(baseline)), tuple(sorted(proposal)), args.score_scale)
    comparison["informative_pairs"] = len(evidence_rows)
    baseline_private = sum(private_counts[allele] for allele in set(baseline) - set(proposal))
    proposal_private = sum(private_counts[allele] for allele in set(proposal) - set(baseline))
    gate_args = SimpleNamespace(
        min_log_bayes_factor=args.min_log_bayes_factor,
        min_log_bayes_factor_per_informative_pair=args.min_normalized_bf,
        min_discriminating_pairs=args.min_discriminating_pairs,
        min_proposal_private_pairs=args.min_proposal_private_pairs,
        min_private_pair_ratio=args.min_private_pair_ratio,
        private_pair_slack=0,
    )
    accepted = proposal_is_supported(comparison, baseline_private, proposal_private, gate_args)
    log_bf = float(comparison["log_bayes_factor"])
    return accepted, {
        "decision": "proposal" if accepted else "baseline",
        "reason": "frozen_normalized_v1_gate",
        "total_pairs": total_pairs,
        "informative_pairs": len(evidence_rows),
        "log_bayes_factor": f"{log_bf:.6f}",
        "normalized_bf": f"{log_bf / len(evidence_rows):.6f}",
        "discriminating_pairs": comparison["discriminating_pairs"],
        "baseline_private_pairs": baseline_private,
        "proposal_private_pairs": proposal_private,
        "private_pair_ratio": f"{(proposal_private + 1) / (baseline_private + 1):.6f}",
    }


def write_calls(
    call_dir: Path,
    fields: list[str],
    source_rows: list[dict[str, str]],
    lifted: list[str],
    rule: str,
) -> None:
    calls = call_dir / "calls.tsv"
    backup = call_dir / "calls.quartet_optimization_input.tsv"
    if not backup.exists():
        shutil.copy2(calls, backup)
    output_fields = list(fields)
    for field in ("quartet_optimization_rule", "quartet_optimization_profile"):
        if field not in output_fields:
            output_fields.append(field)
    ordered = sorted(
        source_rows,
        key=lambda row: int(row.get("global_hap", "0")) if row.get("global_hap", "").isdigit() else 0,
    )
    output_rows = []
    for index, (assignment, allele) in enumerate(zip(ASSIGNMENTS, lifted), 1):
        row = dict(ordered[index - 1])
        row["global_hap"] = str(index)
        row["assignment"] = assignment
        row["allele"] = allele
        for field in ("allele_read_count", "allele_read_fraction", "read_count", "read_fraction"):
            if field in row:
                row[field] = ""
        row["quartet_optimization_rule"] = rule
        row["quartet_optimization_profile"] = "normalized_joint_v1"
        output_rows.append(row)
    write_tsv(calls, output_fields, output_rows)


def aggregate(args: argparse.Namespace) -> None:
    argv = [
        "aggregate_calls.py", "--asm-root", str(args.asm_root), "--sample", args.sample,
        "--genes", *args.genes, "--spechla-root", str(args.spechla_root),
        "--g-group", str(args.g_group), "--out", str(args.asm_root / args.sample / f"{args.sample}.final_calls.tsv"),
    ]
    if args.compact_out:
        argv.extend(["--compact-out", str(args.compact_out)])
    old_argv = sys.argv
    try:
        sys.argv = argv
        aggregate_main()
    finally:
        sys.argv = old_argv


def optimize_gene(args: argparse.Namespace, gene: str, chi_r: float) -> dict[str, object]:
    call_dir = gene_dir(args.asm_root, args.sample, gene)
    fields, rows = baseline_rows(call_dir)
    baseline = baseline_quartet(rows)
    audit: dict[str, object] = {
        "sample": args.sample, "gene": gene, "profile": args.profile,
        "baseline_2field": ",".join(baseline), "proposal_2field": "", "selected_2field": ",".join(baseline),
        "decision": "fallback", "reason": "", "applied": "0",
    }
    if len(rows) != 4 or len(baseline) != 4 or any(not allele or allele == "NA" for allele in baseline):
        audit["reason"] = "invalid_baseline"
        return audit
    counts_path = args.spechla_root / args.sample / "em_refine" / f"{gene}.tf_counts.tsv"
    if not counts_path.exists():
        audit["reason"] = "missing_tf_counts"
        return audit
    try:
        result = call_quartet(read_counts(counts_path), tuple(sorted(baseline)), max(chi_r, 1.0 - chi_r))
    except (OSError, ValueError, ArithmeticError) as error:
        audit["reason"] = f"joint_error:{type(error).__name__}"
        return audit
    proposal = tuple(result["quartet"])
    ordered_proposal = slot_quartet(result["major_group"], result["minor_group"], chi_r, baseline)
    audit.update({
        "proposal_2field": ",".join(proposal),
        "posterior_gap": f"{float(result['posterior_gap']):.6f}",
        "fitted_major_fraction": f"{float(result['major_fraction']):.4f}",
    })
    if gene in CLASS_I_GENES:
        accepted, gate = class_i_read_gate(args, gene, baseline, proposal)
        audit.update(gate)
        for key in (
            "total_pairs", "informative_pairs", "log_bayes_factor", "normalized_bf",
            "discriminating_pairs", "baseline_private_pairs", "proposal_private_pairs", "private_pair_ratio",
        ):
            audit[key] = gate.get(key, "")
        if not accepted:
            return audit
        rule = "class_i_normalized_v1"
    else:
        if tuple(sorted(baseline)) == proposal:
            audit.update({"decision": "same", "reason": "proposal_equals_baseline"})
            return audit
        audit.update({"decision": "proposal", "reason": "class_ii_joint_v2"})
        rule = "class_ii_joint_v2"
    mapping = full_allele_map(args.asm_root, args.spechla_root, args.sample, gene)
    lifted = lift_quartet(ordered_proposal, mapping)
    if lifted is None:
        audit.update({"decision": "fallback", "reason": "missing_full_allele_mapping"})
        return audit
    audit["selected_2field"] = ",".join(ordered_proposal)
    audit["selected_full"] = ",".join(lifted)
    if args.profile == "normalized_joint_v1":
        write_calls(call_dir, fields, rows, lifted, rule)
        audit["applied"] = "1"
    return audit


def safely_optimize_gene(args: argparse.Namespace, gene: str, chi_r: float) -> dict[str, object]:
    try:
        return optimize_gene(args, gene, chi_r)
    except Exception as error:
        return {
            "sample": args.sample,
            "gene": gene,
            "profile": args.profile,
            "baseline_2field": "",
            "proposal_2field": "",
            "selected_2field": "",
            "decision": "fallback",
            "reason": f"unexpected_error:{type(error).__name__}",
            "applied": "0",
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--profile", choices=("shadow", "normalized_joint_v1"), required=True)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--compact-out", type=Path)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-full-alleles", type=int, default=25)
    parser.add_argument("--max-owner-fraction", type=float, default=0.75)
    parser.add_argument("--score-scale", type=float, default=0.35)
    parser.add_argument("--concordance-bonus", type=float, default=0.5)
    parser.add_argument("--min-pair-evidence", type=float, default=1.0)
    parser.add_argument("--min-log-bayes-factor", type=float, default=5.0)
    parser.add_argument("--min-normalized-bf", type=float, default=0.10)
    parser.add_argument("--min-discriminating-pairs", type=int, default=3)
    parser.add_argument("--min-proposal-private-pairs", type=int, default=10)
    parser.add_argument("--min-private-pair-ratio", type=float, default=5.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    chi_path = read_chi_path(args.spechla_root, args.sample)
    if chi_path is None:
        raise SystemExit(f"missing chimerism evidence for {args.sample}")
    major_fraction = read_major_fraction_prior(chi_path)
    text = chi_path.read_text()
    chi_r = major_fraction
    import re
    match = re.search(r"chi_R=([0-9.]+)", text)
    if match:
        chi_r = float(match.group(1))
    audits = [safely_optimize_gene(args, gene, chi_r) for gene in args.genes if gene in CLASS_I_GENES + CLASS_II_GENES]
    fields = [
        "sample", "gene", "profile", "baseline_2field", "proposal_2field", "selected_2field", "selected_full",
        "decision", "reason", "applied", "posterior_gap", "fitted_major_fraction", "total_pairs",
        "informative_pairs", "log_bayes_factor", "normalized_bf", "discriminating_pairs",
        "baseline_private_pairs", "proposal_private_pairs", "private_pair_ratio",
    ]
    write_tsv(args.manifest, fields, audits)
    if args.profile == "normalized_joint_v1" and any(row["applied"] == "1" for row in audits):
        aggregate(args)
    selected = sum(row["decision"] == "proposal" for row in audits)
    applied = sum(row["applied"] == "1" for row in audits)
    print(f"[quartet-optimization:{args.profile}] sample={args.sample} selected={selected} applied={applied} manifest={args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())