#!/usr/bin/env python3
"""Create R/D-free four-copy HLA output from aggregate_calls compact TSV.

The main aggregate table keeps legacy R1/R2/D1/D2 slot names. For mixed donor /
recipient samples those side labels can be biologically ambiguous, so this script
emits copy1-copy4 records sorted by estimated copy proportion without implying
which two copies belong to recipient or donor.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

SLOTS = ("R1", "R2", "D1", "D2")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_fraction(value: str) -> float | None:
    if value in {"", "NA", None}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def format_fraction(value: float) -> str:
    if abs(value) < 1e-4:
        return f"{value:.3e}"
    return f"{value:.6f}"


def allele_2field(allele: str) -> str:
    if not allele or allele == "NA" or "*" not in allele:
        return allele or "NA"
    gene, fields = allele.replace("HLA-", "").split("*", 1)
    parts = fields.replace("G", "").split(":")
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def normalize_copy_fractions(copies: list[dict[str, str]]) -> str:
    parsed = [parse_fraction(copy["raw_copy_fraction"]) for copy in copies]
    if all(value is not None for value in parsed):
        total = sum(value for value in parsed if value is not None)
        if total > 0:
            for copy, value in zip(copies, parsed):
                copy["copy_fraction"] = format_fraction((value or 0.0) / total)
            return "copy_fraction_fit"
    fallback = 1.0 / len(copies) if copies else 0.0
    for copy in copies:
        copy["copy_fraction"] = format_fraction(fallback)
    return "equal_fraction_fallback"


def copy_rows_from_gene(row: dict[str, str]) -> tuple[list[dict[str, str]], dict[str, str]]:
    copies: list[dict[str, str]] = []
    for slot_index, slot in enumerate(SLOTS, 1):
        allele = row.get(f"{slot}_allele", "NA") or "NA"
        copies.append({
            "sample": row.get("sample", ""),
            "gene": row.get("gene", ""),
            "legacy_slot": slot,
            "legacy_slot_index": str(slot_index),
            "allele": allele,
            "allele_2field": allele_2field(allele),
            "raw_copy_fraction": row.get(f"{slot}_copy_fraction", "NA") or "NA",
            "copy_fraction": "NA",
            "copy_identifiability": row.get("copy_identifiability", ""),
            "copy_fit_error": row.get("copy_fit_error", ""),
        })
    source = normalize_copy_fractions(copies)
    for copy in copies:
        copy["proportion_source"] = source
    copies.sort(key=lambda item: (-parse_fraction(item["copy_fraction"]), item["allele"], item["legacy_slot_index"]))
    for rank, copy in enumerate(copies, 1):
        copy["copy_id"] = f"copy{rank}"
        copy["copy_rank"] = str(rank)
    compact = {
        "sample": row.get("sample", ""),
        "gene": row.get("gene", ""),
        "allele_multiset": ",".join(copy["allele"] for copy in copies),
        "allele_2field_multiset": ",".join(copy["allele_2field"] for copy in copies),
        "copy_fractions": ",".join(copy["copy_fraction"] for copy in copies),
        "proportion_source": source,
        "copy_identifiability": row.get("copy_identifiability", ""),
        "copy_fit_error": row.get("copy_fit_error", ""),
    }
    return copies, compact


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create R/D-free HLA copy call outputs from final_calls.compact.tsv")
    parser.add_argument("--compact", required=True, type=Path, help="input <sample>.final_calls.compact.tsv")
    parser.add_argument("--out", required=True, type=Path, help="long-format copy output TSV")
    parser.add_argument("--compact-out", type=Path, default=None, help="one-row-per-gene copy output TSV")
    args = parser.parse_args()

    long_rows: list[dict[str, str]] = []
    compact_rows: list[dict[str, str]] = []
    for row in read_tsv(args.compact):
        copies, compact = copy_rows_from_gene(row)
        long_rows.extend(copies)
        compact_rows.append(compact)

    long_fields = [
        "sample", "gene", "copy_id", "copy_rank", "allele", "allele_2field", "copy_fraction",
        "proportion_source", "legacy_slot", "raw_copy_fraction", "copy_identifiability", "copy_fit_error",
    ]
    compact_fields = [
        "sample", "gene", "allele_multiset", "allele_2field_multiset", "copy_fractions",
        "proportion_source", "copy_identifiability", "copy_fit_error",
    ]
    write_tsv(args.out, long_fields, long_rows)
    if args.compact_out is not None:
        write_tsv(args.compact_out, compact_fields, compact_rows)


if __name__ == "__main__":
    main()
