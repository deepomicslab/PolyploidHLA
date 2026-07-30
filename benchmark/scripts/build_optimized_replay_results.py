#!/usr/bin/env python3
"""Build a non-destructive optimized replay table from completed v4 outputs."""

from __future__ import annotations

import argparse
import csv
import glob
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

MAIN_GENES = {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"}
PRIVATE_GENES = {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1"}


def normalize_2field(allele: str) -> str:
    allele = (allele or "NA").replace("HLA-", "").replace("G", "")
    if allele == "NA" or "*" not in allele:
        return allele
    gene, fields = allele.split("*", 1)
    return f"{gene}*{':'.join(fields.split(':')[:2])}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_private_replay(path: Path) -> dict[tuple[str, str, str, str], dict[str, str]]:
    return {
        (row["experiment"], row["condition"], row["sample_id"], row["gene"]): row
        for row in read_tsv(path)
    }


def load_truth(bench_root: Path, experiment_glob: str) -> dict[tuple[str, str, str, str], list[str]]:
    truth: dict[tuple[str, str, str, str], list[str]] = defaultdict(list)
    for path in sorted((bench_root / "truth").glob(f"{experiment_glob}/copies.tsv")):
        experiment = path.parent.name
        for row in read_tsv(path):
            if row["gene"] in MAIN_GENES:
                truth[(experiment, row["condition"], row["sample_id"], row["gene"])].append(
                    normalize_2field(row["allele_2field"])
                )
    return truth


def load_safe_class2(bench_root: Path, experiment_glob: str) -> dict[tuple[str, str, str, str], list[str]]:
    calls = {}
    pattern = bench_root / "runs" / experiment_glob / "*" / "SIM*" / "asm_v2" / "SIM*" / "hla-*" / "HLA-*" / "calls.tsv"
    for path_text in glob.glob(str(pattern)):
        path = Path(path_text)
        gene = path.parent.name
        if gene not in {"HLA-DQB1", "HLA-DPB1"}:
            continue
        experiment = path.parents[6].name
        condition = path.parents[5].name
        sample = path.parents[4].name
        backup = path.with_name("calls.class2_joint_input.tsv")
        source = backup if backup.exists() else path
        calls[(experiment, condition, sample, gene)] = [
            normalize_2field(row["allele"]) for row in read_tsv(source)
        ]
    return calls


def multiset_delta(old: list[str], new: list[str]) -> tuple[str, str]:
    removed = list((Counter(old) - Counter(new)).elements())
    added = list((Counter(new) - Counter(old)).elements())
    return ",".join(sorted(added)) or ".", ",".join(sorted(removed)) or "."


def write_atomic(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
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


def summarize(rows: list[dict[str, object]], profile: str) -> list[dict[str, object]]:
    scored = [row for row in rows if row["gene"] in MAIN_GENES]
    groups = [("overall", "all", scored)]
    groups.extend(("gene", gene, [row for row in scored if row["gene"] == gene]) for gene in sorted(MAIN_GENES))
    groups.extend(
        ("condition", condition, [row for row in scored if row["condition"] == condition])
        for condition in sorted({str(row["condition"]) for row in scored})
    )
    summary = []
    for scope, name, selected in groups:
        baseline_correct = 0
        optimized_correct = 0
        baseline_exact = 0
        optimized_exact = 0
        changed = 0
        for row in selected:
            expected = row["truth_allele_2field_multiset"].split(",")
            baseline = row["baseline_allele_2field_multiset"].split(",")
            optimized = row["optimized_allele_2field_multiset"].split(",")
            old_hits = sum((Counter(expected) & Counter(baseline)).values())
            new_hits = sum((Counter(expected) & Counter(optimized)).values())
            baseline_correct += old_hits
            optimized_correct += new_hits
            baseline_exact += old_hits == 4
            optimized_exact += new_hits == 4
            changed += Counter(baseline) != Counter(optimized)
        copies = 4 * len(selected)
        summary.append({
            "replay_profile": profile,
            "scope": scope,
            "name": name,
            "sample_loci": len(selected),
            "baseline_correct_copies": baseline_correct,
            "optimized_correct_copies": optimized_correct,
            "delta_correct_copies": optimized_correct - baseline_correct,
            "truth_copies": copies,
            "baseline_recall": f"{baseline_correct / copies:.6f}",
            "optimized_recall": f"{optimized_correct / copies:.6f}",
            "baseline_exact_quartets": baseline_exact,
            "optimized_exact_quartets": optimized_exact,
            "optimized_exact_quartet_accuracy": f"{optimized_exact / len(selected):.6f}",
            "changed_loci": changed,
        })
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench-root", required=True, type=Path)
    parser.add_argument("--baseline-results", required=True, type=Path)
    parser.add_argument("--class-i-replay", required=True, type=Path)
    parser.add_argument("--drb1-replay", required=True, type=Path)
    parser.add_argument("--experiment-glob", default="accuracy_main_v4_shard*")
    parser.add_argument("--profile", default="class_i_v8b+drb1_v1+class2_safe")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"refusing to overwrite existing replay: {args.output}")

    private = load_private_replay(args.class_i_replay)
    private.update(load_private_replay(args.drb1_replay))
    safe_class2 = load_safe_class2(args.bench_root, args.experiment_glob)
    truth = load_truth(args.bench_root, args.experiment_glob)
    output_rows: list[dict[str, object]] = []
    for baseline_row in read_tsv(args.baseline_results):
        key = (
            baseline_row["experiment"], baseline_row["condition"],
            baseline_row["sample"], baseline_row["gene"],
        )
        baseline = [normalize_2field(value) for value in baseline_row["allele_2field_multiset"].split(",")]
        optimized = list(baseline)
        decision = "unchanged"
        candidate = "."
        replaced = "."
        source = "baseline"
        if key in private:
            replay = private[key]
            optimized = [normalize_2field(value) for value in replay["private_v2"].split(",")]
            decision = replay["decision"]
            candidate = replay.get("candidate", ".")
            replaced = replay.get("replaced", ".")
            source = "class_i_v8b" if key[3] != "HLA-DRB1" else "drb1_private_v1"
        elif key in safe_class2:
            optimized = safe_class2[key]
            candidate, replaced = multiset_delta(baseline, optimized)
            decision = "restore_pre_class2" if Counter(optimized) != Counter(baseline) else "unchanged"
            source = "pre_class2_safe"

        expected = truth.get(key, [])
        correct = sum((Counter(expected) & Counter(optimized)).values()) if expected else ""
        output_rows.append({
            "replay_profile": args.profile,
            "per_locus_decision_uses_truth": 0,
            "profile_selected_on_development_truth": 1,
            **{field: baseline_row[field] for field in (
                "experiment", "condition", "scenario", "graft_fraction", "total_coverage",
                "read_length", "insert_mean", "insert_sd", "error_rate", "master_seed",
                "sample", "gene",
            )},
            "replay_source": source,
            "decision": decision,
            "candidate_2field": candidate,
            "replaced_2field": replaced,
            "baseline_allele_2field_multiset": ",".join(baseline),
            "optimized_allele_2field_multiset": ",".join(optimized),
            "truth_allele_2field_multiset": ",".join(sorted(expected)),
            "correct_copies": correct,
            "exact_quartet": int(correct == 4) if expected else "",
        })

    expected_private = 1600
    expected_class2 = 800
    if len(private) != expected_private:
        raise ValueError(f"expected {expected_private} private replay loci, found {len(private)}")
    if len(safe_class2) != expected_class2:
        raise ValueError(f"expected {expected_class2} safe class-II loci, found {len(safe_class2)}")
    if len(output_rows) != 2800:
        raise ValueError(f"expected 2800 output rows, found {len(output_rows)}")
    fields = list(output_rows[0])
    write_atomic(args.output, fields, output_rows)
    if args.summary_output is not None:
        summary = summarize(output_rows, args.profile)
        write_atomic(args.summary_output, list(summary[0]), summary)
    scored = [row for row in output_rows if row["gene"] in MAIN_GENES]
    correct = sum(int(row["correct_copies"]) for row in scored)
    exact = sum(int(row["exact_quartet"]) for row in scored)
    print(f"[optimized-replay] rows={len(output_rows)} correct={correct}/9600 exact={exact}/2400")
    print(f"[optimized-replay] output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())