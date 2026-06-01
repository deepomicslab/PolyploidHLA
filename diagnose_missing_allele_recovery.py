#!/usr/bin/env python3
"""Classify where missed truth alleles remain recoverable in existing outputs."""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_calls import load_truth, norm_allele  # noqa: E402
from rerank_multi_strategy import read_calls_quartet, read_final, read_tf_counts, two_field  # noqa: E402

GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"]
SIDES = ("PATIENT", "DONOR")
SIDE_SLOTS = {"PATIENT": (0, 1), "DONOR": (2, 3)}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sample_set(sample: str) -> str:
    prefix = sample[:1].lower()
    if prefix in {"a", "b", "c", "d"}:
        return f"set-{prefix}"
    raise ValueError(f"cannot infer set from sample name: {sample}")


def discover_samples(asm_root: Path, exclude: set[str]) -> list[str]:
    samples = []
    for path in sorted(asm_root.iterdir()):
        if path.name in exclude or not path.is_dir():
            continue
        final = path / f"{path.name}.final_calls.tsv"
        if final.exists():
            samples.append(path.name)
    return samples


def normalize_values(values: list[str]) -> list[str]:
    return [norm_allele(value, "2field") for value in values]


def missing_with_extras(truth_values: list[str], pred_values: list[str]) -> tuple[list[str], list[str]]:
    truth_counter = Counter(truth_values)
    pred_counter = Counter(pred_values)
    matched = Counter()
    for allele in list(truth_counter):
        matched[allele] = min(truth_counter[allele], pred_counter.get(allele, 0))
    missing = []
    extra = []
    for allele, count in truth_counter.items():
        missing.extend([allele] * (count - matched[allele]))
    for allele, count in pred_counter.items():
        extra.extend([allele] * (count - matched[allele]))
    return missing, extra


def read_match_score_candidates(path: Path, top_n: int) -> set[str]:
    rows = read_tsv(path)
    out: set[str] = set()
    by_block_hap: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_block_hap[(row.get("block", ""), row.get("local_hap", ""))].append(row)
    for group_rows in by_block_hap.values():
        group_rows.sort(key=lambda row: float(row.get("score", "0") or 0), reverse=True)
        for row in group_rows[:top_n]:
            out.add(two_field(row.get("allele", "")))
    return out


def counts_text(counter: Counter[str]) -> str:
    return ";".join(f"{key}:{counter[key]}" for key in sorted(counter))


def classify_rows(args: argparse.Namespace) -> list[dict[str, object]]:
    truths = {
        label: load_truth(args.truth_dir / f"truth_typing-{label}.tsv")
        for label in ("set-a", "set-b", "set-c")
    }
    rows = []
    for sample in discover_samples(args.asm_root, set(args.exclude_samples)):
        label = sample_set(sample)
        truth = truths[label]
        final_path = args.asm_root / sample / f"{sample}.final_calls.tsv"
        finals = read_final(final_path)
        for gene in GENES:
            if gene not in finals:
                continue
            current = finals[gene]
            call_dir = args.asm_root / sample / gene.lower() / gene
            baseline = read_calls_quartet(call_dir / "calls.baseline.tsv") or []
            calls = read_calls_quartet(call_dir / "calls.tsv") or []
            em_top = [allele for allele, _frac in read_tf_counts(args.spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv")[:args.top_n]]
            match_top = read_match_score_candidates(call_dir / "match_scores.tsv", args.top_n)
            for side in SIDES:
                slots = SIDE_SLOTS[side]
                truth_side = normalize_values(truth[side][gene])
                pred_side = [current[index] for index in slots]
                missing, extra = missing_with_extras(truth_side, pred_side)
                for allele in missing:
                    rows.append({
                        "sample": sample,
                        "set": label,
                        "gene": gene,
                        "side": side,
                        "missing_allele": allele,
                        "extra_alleles": ",".join(extra),
                        "truth_side": ",".join(truth_side),
                        "pred_side": ",".join(pred_side),
                        "current_quartet": ",".join(current),
                        "baseline_quartet": ",".join(baseline),
                        "calls_quartet": ",".join(calls),
                        "in_current_quartet": allele in current,
                        "in_opposite_side": allele in [current[index] for index in ({0, 1, 2, 3} - set(slots))],
                        "in_baseline_quartet": allele in baseline,
                        "in_calls_quartet": allele in calls,
                        "in_em_top": allele in em_top,
                        "in_match_scores_top": allele in match_top,
                        "recoverability": recoverability(allele, current, slots, baseline, calls, em_top, match_top),
                    })
    return rows


def recoverability(allele: str, current: list[str], slots: tuple[int, int], baseline: list[str], calls: list[str],
                   em_top: list[str], match_top: set[str]) -> str:
    same_side = [current[index] for index in slots]
    other_side = [current[index] for index in ({0, 1, 2, 3} - set(slots))]
    if allele in same_side:
        return "copy_underreported_same_side"
    if allele in other_side:
        return "role_swap_or_cross_side"
    if allele in current:
        return "current_quartet_other"
    if allele in baseline:
        return "baseline_candidate"
    if allele in calls:
        return "calls_candidate"
    if allele in em_top and allele in match_top:
        return "em_and_match_candidate"
    if allele in em_top:
        return "em_candidate"
    if allele in match_top:
        return "match_candidate"
    return "not_in_top_candidates"


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for row in rows:
        grouped[(str(row["gene"]), str(row["recoverability"]))][str(row["missing_allele"])] += 1
    out = []
    for (gene, recoverability), counter in sorted(grouped.items()):
        out.append({
            "gene": gene,
            "recoverability": recoverability,
            "missing_count": sum(counter.values()),
            "allele_counts": counts_text(counter),
        })
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--truth-dir", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--exclude-samples", nargs="*", default=[])
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    rows = classify_rows(args)
    fields = [
        "sample", "set", "gene", "side", "missing_allele", "extra_alleles",
        "truth_side", "pred_side", "current_quartet", "baseline_quartet", "calls_quartet",
        "in_current_quartet", "in_opposite_side", "in_baseline_quartet", "in_calls_quartet",
        "in_em_top", "in_match_scores_top", "recoverability",
    ]
    write_tsv(args.out_prefix.with_suffix(".tsv"), fields, rows)
    write_tsv(args.out_prefix.with_suffix(".summary.tsv"), ["gene", "recoverability", "missing_count", "allele_counts"], summarize(rows))
    print(f"wrote {args.out_prefix.with_suffix('.tsv')}")
    print(f"wrote {args.out_prefix.with_suffix('.summary.tsv')}")


if __name__ == "__main__":
    main()
