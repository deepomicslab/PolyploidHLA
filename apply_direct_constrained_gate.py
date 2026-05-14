#!/usr/bin/env python3
"""Apply constrained direct-likelihood quartet overrides to an ASM root.

The gate is truth-free: accept a direct quartet only when it differs from the
current quartet and its likelihood gap is at least the configured threshold.
Accepted 2-field alleles are lifted back to full allele strings by reusing full
alleles already present in current/baseline/EM calls for the same 2-field name.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from collections import defaultdict, deque
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from aggregate_calls import DEFAULT_GENES, allele_2field, main as aggregate_main  # noqa: E402


ASSIGNMENTS = ("R", "R", "D", "D")


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_quartet(row):
    return [x for x in row["direct_R"].split(",") if x] + [x for x in row["direct_D"].split(",") if x]


def current_quartet(row):
    return [x for x in row["current_quartet"].split(",") if x]


def gene_dir(asm_root: Path, sample: str, gene: str) -> Path:
    return asm_root / sample / gene.lower() / gene


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
        fields = reader.fieldnames or []
        return fields, list(reader)


def add_mapping(mapping, rows) -> None:
    for row in rows:
        allele = row.get("allele") or row.get("allele_2field") or row.get("call")
        if not allele or allele == "NA":
            continue
        mapping[allele_2field(allele)].append(allele)


def build_full_allele_map(call_dir: Path, spechla_root: Path | None, sample: str, gene: str):
    mapping = defaultdict(deque)
    for name in ("calls.tsv", "calls.baseline.tsv", "calls.direct_gate_input.tsv"):
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


def write_direct_calls(call_dir: Path, lifted_quartet, source_row) -> None:
    calls = call_dir / "calls.tsv"
    if not calls.exists():
        raise FileNotFoundError(calls)
    backup = call_dir / "calls.direct_gate_input.tsv"
    if not backup.exists():
        shutil.copy2(calls, backup)
    fields = ["global_hap", "assignment", "allele", "hap_fraction", "direct_likelihood_gap", "direct_chi", "direct_n_reads"]
    rows = []
    chi = float(source_row["chi_direct"])
    for index, (assignment, allele) in enumerate(zip(ASSIGNMENTS, lifted_quartet), 1):
        hap_fraction = chi / 2.0 if assignment == "R" else (1.0 - chi) / 2.0
        rows.append({
            "global_hap": str(index),
            "assignment": assignment,
            "allele": allele,
            "hap_fraction": f"{hap_fraction:.6f}",
            "direct_likelihood_gap": source_row["gap"],
            "direct_chi": source_row["chi_direct"],
            "direct_n_reads": source_row["n_reads"],
        })
    write_tsv(calls, fields, rows)


def accepted_direct_rows(path: Path, gap_threshold: float):
    accepted = []
    for row in read_tsv(path):
        if float(row["gap"]) < gap_threshold:
            continue
        if split_quartet(row) == current_quartet(row):
            continue
        accepted.append(row)
    return accepted


def aggregate_sample(asm_root: Path, sample: str, genes, g_group: Path) -> None:
    argv = [
        "aggregate_calls.py",
        "--asm-root", str(asm_root),
        "--sample", sample,
        "--g-group", str(g_group),
        "--genes", *genes,
        "--out", str(asm_root / sample / f"{sample}.final_calls.tsv"),
    ]
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
        row["source"] = "direct-constrained"
        warning = row.get("warning", "")
        tag = "direct_constrained_gap_gate"
        row["warning"] = tag if not warning else f"{warning};{tag}"
    write_tsv(final_path, fields, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct-tsv", required=True, type=Path)
    parser.add_argument("--in-asm-root", required=True, type=Path)
    parser.add_argument("--out-asm-root", type=Path, default=None)
    parser.add_argument("--spechla-root", type=Path, default=None)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--gap-threshold", type=float, default=150.0)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--in-place", action="store_true",
                        help="Apply accepted overrides directly under --in-asm-root. "
                             "Each edited calls.tsv is backed up as calls.direct_gate_input.tsv.")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.in_place:
        out_asm_root = args.in_asm_root
    else:
        if args.out_asm_root is None:
            raise SystemExit("--out-asm-root is required unless --in-place is set")
        out_asm_root = args.out_asm_root

    accepted = accepted_direct_rows(args.direct_tsv, args.gap_threshold)
    samples = sorted({row["sample"] for row in read_tsv(args.direct_tsv)})
    accepted_keys = {(row["sample"], row["gene"]) for row in accepted}
    manifest_fields = [
        "sample", "gene", "gap", "n_reads", "chi_direct", "current_quartet",
        "direct_2field_quartet", "direct_full_quartet", "output_calls",
    ]
    manifest_rows = []

    if not args.dry_run:
        out_asm_root.mkdir(parents=True, exist_ok=True)
        if not args.in_place:
            for sample in samples:
                copy_sample_tree(args.in_asm_root, out_asm_root, sample)

    for row in accepted:
        sample = row["sample"]
        gene = row["gene"]
        call_dir = gene_dir(out_asm_root if not args.dry_run else args.in_asm_root, sample, gene)
        direct_2field = split_quartet(row)
        mapping = build_full_allele_map(call_dir, args.spechla_root, sample, gene)
        direct_full = lift_alleles(direct_2field, mapping)
        if not args.dry_run:
            write_direct_calls(call_dir, direct_full, row)
        manifest_rows.append({
            "sample": sample,
            "gene": gene,
            "gap": row["gap"],
            "n_reads": row["n_reads"],
            "chi_direct": row["chi_direct"],
            "current_quartet": row["current_quartet"],
            "direct_2field_quartet": ",".join(direct_2field),
            "direct_full_quartet": ",".join(direct_full),
            "output_calls": str(call_dir / "calls.tsv"),
        })
        print(f"[direct-gate] {sample} {gene}: gap={row['gap']} -> override", flush=True)

    if not args.dry_run:
        for sample in samples:
            aggregate_sample(out_asm_root, sample, args.genes, args.g_group)
            mark_final_rows(out_asm_root / sample / f"{sample}.final_calls.tsv", accepted_keys)

    manifest = args.manifest or (out_asm_root / "direct_constrained_gate_manifest.tsv")
    if not args.dry_run:
        write_tsv(manifest, manifest_fields, manifest_rows)
    print(
        f"accepted={len(accepted)} samples={len(samples)} gap_threshold={args.gap_threshold} "
        f"out={out_asm_root}",
        flush=True,
    )


if __name__ == "__main__":
    main()
