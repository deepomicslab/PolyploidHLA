#!/usr/bin/env python3
"""Apply truth-free private-read allele rescue to one PolyploidHLA sample."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path

from diagnostics.direct_read_quartet_likelihood import DEFAULT_IMGT
from offline_class_i_private_rescue import (
    cached_private_pair_support,
    normalize_2field,
    propose_rescue,
    read_em_counts,
)


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_rows_atomic(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def apply_gene(args: argparse.Namespace, gene: str) -> dict[str, str]:
    call_dir = args.asm_root / args.sample / gene.lower() / gene
    calls_path = call_dir / "calls.tsv"
    backup_path = call_dir / "calls.pre_private_rescue.tsv"
    if not calls_path.exists():
        return {"sample": args.sample, "gene": gene, "decision": "missing_calls"}

    source_path = backup_path if backup_path.exists() else calls_path
    fields, rows = read_rows(source_path)
    baseline = [normalize_2field(row.get("allele", "")) for row in rows]
    sample_root = args.spechla_root / args.sample
    counts = read_em_counts(sample_root / "em_refine" / f"{gene}.tf_counts.tsv")
    em_gap_override = gene in args.em_gap_genes
    preliminary = propose_rescue(
        baseline, counts, {}, True, args.top_n, em_gap_override, args.top_n,
        args.min_em_fraction, args.min_em_count, args.min_em_gap, 0,
        args.weak_singleton_max, args.candidate_weak_ratio,
    )
    evidence_candidates = preliminary[-1]
    comparison_support: dict[str, Counter[str]] = {}
    if evidence_candidates:
        short = gene.removeprefix("HLA-")
        comparison_support = cached_private_pair_support(
            args.support_cache,
            sample_root / f"{short}.R1.fq.gz",
            sample_root / f"{short}.R2.fq.gz",
            gene,
            baseline,
            evidence_candidates,
            args.imgt,
            args.k,
            args.max_full_alleles,
        )
    rescued, decision, candidate, replaced, candidate_fraction, _ = propose_rescue(
        baseline, counts, comparison_support, True, args.top_n,
        em_gap_override, args.top_n, args.min_em_fraction, args.min_em_count,
        args.min_em_gap, args.min_candidate_private, args.weak_singleton_max,
        args.candidate_weak_ratio,
    )
    support = comparison_support.get(candidate, Counter())
    second_decision = "not_run"
    second_candidate = ""
    second_replaced = ""
    second_candidate_fraction = 0.0
    second_support: Counter[str] = Counter()
    if gene in args.second_pass_genes and decision == "rescue":
        second_preliminary = propose_rescue(
            rescued, counts, {}, True, args.top_n, em_gap_override, args.top_n,
            args.min_em_fraction, args.min_em_count, args.min_em_gap, 0,
            args.second_weak_singleton_max, args.second_candidate_weak_ratio,
            blocked_candidates={replaced}, require_weak_singleton=True,
        )
        second_evidence_candidates = second_preliminary[-1]
        second_comparison_support: dict[str, Counter[str]] = {}
        if second_evidence_candidates:
            short = gene.removeprefix("HLA-")
            second_comparison_support = cached_private_pair_support(
                args.support_cache,
                sample_root / f"{short}.R1.fq.gz",
                sample_root / f"{short}.R2.fq.gz",
                gene,
                rescued,
                second_evidence_candidates,
                args.imgt,
                args.k,
                args.max_full_alleles,
            )
        (
            second_rescued,
            second_decision,
            second_candidate,
            second_replaced,
            second_candidate_fraction,
            _,
        ) = propose_rescue(
            rescued, counts, second_comparison_support, True, args.top_n,
            em_gap_override, args.top_n, args.min_em_fraction, args.min_em_count,
            args.min_em_gap, args.second_min_candidate_private,
            args.second_weak_singleton_max, args.second_candidate_weak_ratio,
            blocked_candidates={replaced}, require_weak_singleton=True,
        )
        second_support = second_comparison_support.get(second_candidate, Counter())
        if second_decision == "rescue":
            rescued = second_rescued
            decision = "rescue_twice"
    audit = {
        "sample": args.sample,
        "gene": gene,
        "decision": decision,
        "candidate": candidate or ".",
        "replaced": replaced or ".",
        "candidate_em_fraction": f"{candidate_fraction:.6f}",
        "candidate_private_pairs": str(support[candidate] if candidate else 0),
        "replaced_private_pairs": str(support[replaced] if replaced else 0),
        "second_decision": second_decision,
        "second_candidate": second_candidate or ".",
        "second_replaced": second_replaced or ".",
        "second_candidate_em_fraction": f"{second_candidate_fraction:.6f}",
        "second_candidate_private_pairs": str(second_support[second_candidate] if second_candidate else 0),
        "second_replaced_private_pairs": str(second_support[second_replaced] if second_replaced else 0),
        "baseline_2field": ",".join(baseline),
        "rescued_2field": ",".join(rescued),
    }
    if decision not in {"rescue", "rescue_twice"}:
        return audit

    if not backup_path.exists():
        shutil.copy2(calls_path, backup_path)
    working = list(baseline)
    replacements = [(candidate, replaced, candidate_fraction)]
    if decision == "rescue_twice":
        replacements.append((second_candidate, second_replaced, second_candidate_fraction))
    for incoming, outgoing, fraction in replacements:
        replaced_index = working.index(outgoing)
        working[replaced_index] = incoming
        rows[replaced_index]["allele"] = incoming
        for field in ("allele_read_count", "em_weight"):
            if field in fields:
                rows[replaced_index][field] = f"{counts[incoming]:.2f}"
        if "allele_read_fraction" in fields:
            rows[replaced_index]["allele_read_fraction"] = f"{fraction:.6f}"
    write_rows_atomic(calls_path, fields, rows)
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genes", nargs="+", default=["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"])
    parser.add_argument("--em-gap-genes", nargs="+", default=["HLA-A", "HLA-B", "HLA-C"])
    parser.add_argument("--second-pass-genes", nargs="*", default=["HLA-A", "HLA-B", "HLA-C"])
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--support-cache", type=Path)
    parser.add_argument("--imgt", type=Path, default=DEFAULT_IMGT)
    parser.add_argument("--top-n", type=int, default=8)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--max-full-alleles", type=int, default=25)
    parser.add_argument("--min-em-fraction", type=float, default=0.005)
    parser.add_argument("--min-em-count", type=float, default=20.0)
    parser.add_argument("--min-em-gap", type=float, default=1.5)
    parser.add_argument("--min-candidate-private", type=int, default=30)
    parser.add_argument("--weak-singleton-max", type=int, default=10)
    parser.add_argument("--candidate-weak-ratio", type=float, default=3.0)
    parser.add_argument("--second-min-candidate-private", type=int, default=50)
    parser.add_argument("--second-weak-singleton-max", type=int, default=5)
    parser.add_argument("--second-candidate-weak-ratio", type=float, default=5.0)
    args = parser.parse_args()

    audits = [apply_gene(args, gene) for gene in args.genes]
    fields = [
        "sample", "gene", "decision", "candidate", "replaced",
        "candidate_em_fraction", "candidate_private_pairs", "replaced_private_pairs",
        "second_decision", "second_candidate", "second_replaced",
        "second_candidate_em_fraction", "second_candidate_private_pairs", "second_replaced_private_pairs",
        "baseline_2field", "rescued_2field",
    ]
    write_rows_atomic(args.manifest, fields, audits)
    rescued_count = sum(row["decision"] in {"rescue", "rescue_twice"} for row in audits)
    print(f"[private-read-rescue] sample={args.sample} rescued={rescued_count} manifest={args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())