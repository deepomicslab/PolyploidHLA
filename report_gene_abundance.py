#!/usr/bin/env python3
"""Write per-gene mixture and conditional copy-abundance diagnostics."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from estimate_chi_pooled import QC_GENE_MIN_AF, estimate_chi_from_af, parse_vcf


FIXED_DIPLOID_GENES = (
    "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1",
    "HLA-E", "HLA-F", "HLA-G", "HLA-H", "MICA", "MICB",
)
DRB_GENES = ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5")
FIELDS = (
    "sample", "gene", "model", "global_chi", "global_chi_source",
    "global_qc_status", "global_ci95", "local_chi", "n_af", "median_dp",
    "residual", "delta_from_global", "low_source_copies",
    "high_source_copies", "low_source_fraction", "high_source_fraction",
    "expected_gene_abundance", "observed_low_fraction",
    "observed_read_fraction", "called_copies", "status", "reasons",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_key_values(line: str) -> dict[str, str]:
    values = {}
    for token in line.split():
        if "=" in token:
            key, value = token.split("=", 1)
            values[key] = value
    return values


def finite_float(value: str | None) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def read_global_chi(pooled_log: Path, chimerism_log: Path) -> dict[str, object]:
    pooled = {}
    if pooled_log.exists():
        for line in pooled_log.read_text().splitlines():
            if line.startswith("GLOBAL"):
                pooled = parse_key_values(line)
                break
    pooled_chi = finite_float(pooled.get("chi_R"))
    raw_chi = finite_float(pooled.get("raw_chi_R"))
    if pooled.get("status") == "PASS" and pooled_chi is not None:
        return {
            "chi": pooled_chi, "source": "pooled_pass", "status": "PASS",
            "ci95": pooled.get("ci95", "NA"),
        }

    fallback = None
    if chimerism_log.exists():
        for line in chimerism_log.read_text().splitlines():
            values = parse_key_values(line)
            candidate = finite_float(values.get("chi_R"))
            if candidate is not None:
                fallback = candidate
    if fallback is not None:
        fallback = min(fallback, 1.0 - fallback)
        return {
            "chi": fallback, "source": "gt_fallback",
            "status": pooled.get("status", "NOT_AVAILABLE"),
            "ci95": pooled.get("ci95", "NA"),
        }
    return {
        "chi": raw_chi, "source": "unaccepted_raw" if raw_chi is not None else "none",
        "status": pooled.get("status", "NOT_AVAILABLE"),
        "ci95": pooled.get("ci95", "NA"),
    }


def read_gene_contigs(gene_bed: Path) -> dict[str, str]:
    mapping = {}
    with gene_bed.open() as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 4 and not line.startswith("#"):
                mapping[parts[3]] = parts[0]
    return mapping


def fmt(value: float | None, digits: int = 4) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def observed_minor_fraction(row: dict[str, str] | None) -> float | None:
    if row is None:
        return None
    values = sorted(
        value for value in (
            finite_float(row.get(f"{slot}_copy_fraction"))
            for slot in ("R1", "R2", "D1", "D2")
        )
        if value is not None
    )
    return sum(values[:2]) if len(values) == 4 else None


def fixed_gene_rows(
    sample: str,
    vcf: Path,
    gene_bed: Path,
    compact: Path,
    global_info: dict[str, object],
) -> list[dict[str, str]]:
    contigs = read_gene_contigs(gene_bed)
    observations: dict[str, list[tuple[float, int]]] = {}
    for chrom, _pos, af, dp in parse_vcf(str(vcf)):
        observations.setdefault(chrom, []).append((af, dp))
    calls = {row.get("gene", ""): row for row in read_tsv(compact)}
    global_chi = global_info["chi"]
    rows = []
    for gene in FIXED_DIPLOID_GENES:
        call = calls.get(gene)
        contig = contigs.get(gene)
        values = observations.get(contig, []) if contig else []
        local = None
        if len(values) >= QC_GENE_MIN_AF:
            local = estimate_chi_from_af(
                [item[0] for item in values], dps=[item[1] for item in values],
                prior_chi=global_chi,
            )
        local_chi = local["chi_r"] if local else None
        delta = (
            abs(local_chi - global_chi)
            if local_chi is not None and global_chi is not None else None
        )
        reasons = []
        if call is None and contig is None:
            status = "NOT_ENABLED"
            reasons.append("gene_not_enabled")
        elif call is None:
            status = "LOW_CONFIDENCE"
            reasons.append("missing_final_call")
        elif global_info["source"] != "pooled_pass":
            status = "LOW_CONFIDENCE"
            reasons.append(
                "global_gt_fallback"
                if global_info["source"] == "gt_fallback"
                else "missing_accepted_global_chi"
            )
        elif local is None:
            status = "LOW_CONFIDENCE"
            reasons.append("no_reliable_local_peak")
        elif delta is not None and delta >= 0.15:
            status = "MODEL_MISMATCH"
            reasons.append("local_global_delta")
        else:
            status = "PASS"
        observed_low = observed_minor_fraction(call)
        dps = [item[1] for item in values]
        rows.append({
            "sample": sample,
            "gene": gene,
            "model": "fixed_diploid",
            "global_chi": fmt(global_chi),
            "global_chi_source": str(global_info["source"]),
            "global_qc_status": str(global_info["status"]),
            "global_ci95": str(global_info["ci95"]),
            "local_chi": fmt(local_chi),
            "n_af": str(len(values)),
            "median_dp": fmt(float(np.median(dps)) if dps else None, 1),
            "residual": fmt(local["weighted_residual"] if local else None, 5),
            "delta_from_global": fmt(delta),
            "low_source_copies": "2",
            "high_source_copies": "2",
            "low_source_fraction": fmt(global_chi),
            "high_source_fraction": fmt(1.0 - global_chi if global_chi is not None else None),
            "expected_gene_abundance": "1.0000" if global_chi is not None else "NA",
            "observed_low_fraction": fmt(observed_low),
            "observed_read_fraction": "NA",
            "called_copies": "4" if call else "0",
            "status": status,
            "reasons": ",".join(reasons) if reasons else "none",
        })
    return rows


def drb_gene_rows(
    sample: str,
    calls_path: Path,
    tf_counts_path: Path,
    global_info: dict[str, object],
) -> list[dict[str, str]]:
    calls = read_tsv(calls_path)
    tf_counts = read_tsv(tf_counts_path)
    global_chi = global_info["chi"]
    total_read_count = sum(finite_float(row.get("em_weight")) or 0.0 for row in tf_counts)
    rows = []
    for gene in DRB_GENES:
        locus = gene.replace("HLA-", "")
        linked = [row for row in calls if row.get("drb1_linked_locus") == locus]
        low = sum(
            1 for row in linked
            if (finite_float(row.get("hap_fraction")) or 0.0) <= 0.25
        )
        high = len(linked) - low
        called = [
            row for row in linked
            if row.get("allele", "NA") not in {"", "NA"}
        ]
        read_count = sum(
            finite_float(row.get("em_weight")) or 0.0
            for row in tf_counts if row.get("locus") == locus
        )
        expected = (
            (low * global_chi + high * (1.0 - global_chi)) / 2.0
            if global_chi is not None else None
        )
        reasons = []
        if not calls:
            status = "NOT_RUN"
            reasons.append("missing_drb345_calls")
        elif not linked:
            status = "PASS"
            reasons.append("structural_absence")
        elif len(called) < len(linked):
            status = "LOW_CONFIDENCE"
            reasons.append("linked_copy_without_allele")
        elif read_count <= 0:
            status = "LOW_CONFIDENCE"
            reasons.append("no_read_support")
        elif global_info["source"] != "pooled_pass":
            status = "LOW_CONFIDENCE"
            reasons.append(
                "global_gt_fallback"
                if global_info["source"] == "gt_fallback"
                else "missing_accepted_global_chi"
            )
        else:
            status = "PASS"
        rows.append({
            "sample": sample,
            "gene": gene,
            "model": "drb1_linked_conditional_copy",
            "global_chi": fmt(global_chi),
            "global_chi_source": str(global_info["source"]),
            "global_qc_status": str(global_info["status"]),
            "global_ci95": str(global_info["ci95"]),
            "local_chi": "NA",
            "n_af": "NA",
            "median_dp": "NA",
            "residual": "NA",
            "delta_from_global": "NA",
            "low_source_copies": str(low),
            "high_source_copies": str(high),
            "low_source_fraction": fmt(global_chi),
            "high_source_fraction": fmt(1.0 - global_chi if global_chi is not None else None),
            "expected_gene_abundance": fmt(expected),
            "observed_low_fraction": "NA",
            "observed_read_fraction": fmt(read_count / total_read_count if total_read_count else None),
            "called_copies": str(len(called)),
            "status": status,
            "reasons": ",".join(reasons) if reasons else "none",
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pooled-vcf", required=True, type=Path)
    parser.add_argument("--pooled-log", required=True, type=Path)
    parser.add_argument("--chimerism-log", required=True, type=Path)
    parser.add_argument("--gene-bed", required=True, type=Path)
    parser.add_argument("--compact-calls", required=True, type=Path)
    parser.add_argument("--drb345-calls", required=True, type=Path)
    parser.add_argument("--drb345-tf-counts", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    global_info = read_global_chi(args.pooled_log, args.chimerism_log)
    rows = fixed_gene_rows(
        args.sample, args.pooled_vcf, args.gene_bed, args.compact_calls,
        global_info,
    )
    rows.extend(drb_gene_rows(
        args.sample, args.drb345_calls, args.drb345_tf_counts, global_info,
    ))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[gene-abundance] wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()