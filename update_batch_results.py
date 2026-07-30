#!/usr/bin/env python3
"""Merge one sample's compact copy calls into a shared batch result TSV."""

from __future__ import annotations

import argparse
import csv
import fcntl
import os
import tempfile
from pathlib import Path

METADATA_FIELDS = [
    "experiment", "condition", "scenario", "graft_fraction", "total_coverage",
    "read_length", "insert_mean", "insert_sd", "error_rate", "master_seed",
]
RESULT_FIELDS = [
    "sample", "gene", "allele_multiset", "allele_2field_multiset", "copy_fractions",
    "allele_read_counts", "copy_read_counts", "proportion_source",
    "copy_identifiability", "copy_fit_error",
]
FIELDS = METADATA_FIELDS + RESULT_FIELDS


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if list(reader.fieldnames or []) != FIELDS:
            raise ValueError(f"unexpected batch result header in {path}")
        return list(reader)


def parse_metadata(values: list[str]) -> dict[str, str]:
    metadata = {field: "" for field in METADATA_FIELDS}
    for value in values:
        if "=" not in value:
            raise ValueError(f"metadata must use KEY=VALUE syntax: {value!r}")
        key, item = value.split("=", 1)
        if key not in metadata:
            raise ValueError(f"unsupported metadata key: {key}")
        metadata[key] = item
    return metadata


def write_atomic(path: Path, rows: list[dict[str, str]]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-result", required=True, type=Path)
    parser.add_argument("--batch-result", required=True, type=Path)
    parser.add_argument("--metadata", action="append", default=[], metavar="KEY=VALUE")
    args = parser.parse_args()

    metadata = parse_metadata(args.metadata)
    with args.sample_result.open() as handle:
        sample_rows = list(csv.DictReader(handle, delimiter="\t"))
    incoming = [{**metadata, **{field: row.get(field, "") for field in RESULT_FIELDS}} for row in sample_rows]
    if not incoming:
        raise ValueError(f"no result rows found in {args.sample_result}")

    args.batch_result.parent.mkdir(parents=True, exist_ok=True)
    lock_path = args.batch_result.with_suffix(args.batch_result.suffix + ".lock")
    with lock_path.open("a+") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        existing = read_tsv(args.batch_result)
        incoming_keys = {
            (row["experiment"], row["condition"], row["sample"], row["gene"])
            for row in incoming
        }
        combined = [
            row for row in existing
            if (row["experiment"], row["condition"], row["sample"], row["gene"]) not in incoming_keys
        ]
        combined.extend(incoming)
        combined.sort(key=lambda row: (row["experiment"], row["condition"], row["sample"], row["gene"]))
        write_atomic(args.batch_result, combined)
        fcntl.flock(lock_handle, fcntl.LOCK_UN)
    print(f"[batch-results] updated {args.batch_result} rows={len(combined)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())