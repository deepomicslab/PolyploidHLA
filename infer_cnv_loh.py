#!/usr/bin/env python3
"""Infer source-resolved HLA CNV/LOH states with a joint dosage MILP."""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix


SLOTS = ("R1", "R2", "D1", "D2")
CORE_GENES = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")
FIELDS = (
    "sample", "gene", "chi", "raw_median_depth", "normalized_depth", "breadth",
    "R1_dosage", "R2_dosage", "D1_dosage", "D2_dosage", "total_copies",
    "event", "confidence", "best_objective", "second_objective", "objective_gap",
    "normal_objective", "event_support", "allele_groups", "reasons",
)


@dataclass(frozen=True)
class DepthEvidence:
    median: float
    mean: float
    breadth: float


@dataclass(frozen=True)
class MilpResult:
    state: tuple[int, int, int, int]
    objective: float


def finite_float(value: str | None) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_intervals(path: Path) -> dict[str, tuple[str, int, int]]:
    intervals = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, gene, *_ = line.split()
            intervals[gene] = (chrom, int(start), int(end))
    return intervals


def bam_depth(bam: pysam.AlignmentFile, chrom: str, start: int, end: int,
              breadth_min_depth: int) -> DepthEvidence | None:
    if chrom not in bam.references or end <= start:
        return None
    coverage = np.zeros(end - start, dtype=np.int32)
    for column in bam.pileup(
        chrom, start, end, truncate=True, stepper="nofilter",
        min_base_quality=0, min_mapping_quality=0,
    ):
        depth = 0
        for pileup_read in column.pileups:
            alignment = pileup_read.alignment
            if (alignment.is_unmapped or alignment.is_secondary or
                    alignment.is_supplementary or alignment.is_duplicate or
                    pileup_read.is_del or pileup_read.is_refskip):
                continue
            depth += 1
        coverage[column.reference_pos - start] = depth
    return DepthEvidence(
        median=float(np.median(coverage)),
        mean=float(np.mean(coverage)),
        breadth=float(np.mean(coverage >= breadth_min_depth)),
    )


def allele_groups(row: dict[str, str]) -> list[tuple[str, tuple[int, ...], float]]:
    grouped: dict[str, dict[str, object]] = {}
    for index, slot in enumerate(SLOTS):
        allele = row.get(f"{slot}_allele", "NA")
        fraction = finite_float(row.get(f"{slot}_copy_fraction"))
        if not allele or allele == "NA" or fraction is None:
            continue
        entry = grouped.setdefault(allele, {"indices": [], "fraction": 0.0})
        entry["indices"].append(index)
        entry["fraction"] = float(entry["fraction"]) + max(0.0, fraction)
    total = sum(float(entry["fraction"]) for entry in grouped.values())
    if total <= 0:
        return []
    return [
        (allele, tuple(entry["indices"]), float(entry["fraction"]) / total)
        for allele, entry in grouped.items()
    ]


def _state_from_bits(values: np.ndarray) -> tuple[int, int, int, int]:
    return tuple(int(round(values[2 * index])) + 2 * int(round(values[2 * index + 1]))
                 for index in range(4))


def solve_dosage_milp(
    normalized_depth: float,
    chi: float,
    groups: list[tuple[str, tuple[int, ...], float]],
    prior_weight: float = 0.06,
    r_is_low: bool = True,
    exclude_bits: tuple[int, ...] | None = None,
    fixed_state: tuple[int, int, int, int] | None = None,
) -> MilpResult:
    r_weight, d_weight = (chi, 1.0 - chi) if r_is_low else (1.0 - chi, chi)
    weights = (r_weight, r_weight, d_weight, d_weight)
    observations = [(indices, normalized_depth * fraction) for _, indices, fraction in groups]
    observations.append((tuple(range(4)), normalized_depth))
    observation_weights = [1.0] * len(groups) + [1.5]

    bit_count = 8
    residual_offset = bit_count
    prior_offset = residual_offset + 2 * len(observations)
    variable_count = prior_offset + 8
    objective = np.zeros(variable_count)
    lower = np.zeros(variable_count)
    upper = np.full(variable_count, np.inf)
    upper[:bit_count] = 1.0
    integrality = np.zeros(variable_count, dtype=int)
    integrality[:bit_count] = 1

    rows = len(observations) + 4 + (1 if exclude_bits is not None else 0)
    matrix = lil_matrix((rows, variable_count), dtype=float)
    constraint_lower = np.zeros(rows)
    constraint_upper = np.zeros(rows)

    for obs_index, ((indices, observed), obs_weight) in enumerate(
        zip(observations, observation_weights)
    ):
        for slot_index in indices:
            matrix[obs_index, 2 * slot_index] = weights[slot_index] / 2.0
            matrix[obs_index, 2 * slot_index + 1] = weights[slot_index]
        positive = residual_offset + 2 * obs_index
        negative = positive + 1
        matrix[obs_index, positive] = 1.0
        matrix[obs_index, negative] = -1.0
        objective[positive] = obs_weight
        objective[negative] = obs_weight
        constraint_lower[obs_index] = observed
        constraint_upper[obs_index] = observed

    for slot_index in range(4):
        row_index = len(observations) + slot_index
        matrix[row_index, 2 * slot_index] = 1.0
        matrix[row_index, 2 * slot_index + 1] = 2.0
        positive = prior_offset + 2 * slot_index
        negative = positive + 1
        matrix[row_index, positive] = 1.0
        matrix[row_index, negative] = -1.0
        objective[positive] = prior_weight
        objective[negative] = prior_weight
        constraint_lower[row_index] = 1.0
        constraint_upper[row_index] = 1.0

    if exclude_bits is not None:
        row_index = rows - 1
        one_count = 0
        for bit_index, bit in enumerate(exclude_bits):
            if bit:
                matrix[row_index, bit_index] = -1.0
                one_count += 1
            else:
                matrix[row_index, bit_index] = 1.0
        constraint_lower[row_index] = 1.0 - one_count
        constraint_upper[row_index] = np.inf

    if fixed_state is not None:
        for slot_index, dosage in enumerate(fixed_state):
            lower[2 * slot_index] = upper[2 * slot_index] = dosage & 1
            lower[2 * slot_index + 1] = upper[2 * slot_index + 1] = dosage >> 1

    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(lower, upper),
        constraints=LinearConstraint(matrix.tocsr(), constraint_lower, constraint_upper),
        options={"presolve": True},
    )
    if not result.success or result.x is None:
        raise RuntimeError(f"MILP failed: {result.message}")
    return MilpResult(_state_from_bits(result.x), float(result.fun))


def classify_state(state: tuple[int, int, int, int]) -> str:
    if state == (1, 1, 1, 1):
        return "NORMAL"
    if sum(state) != 4:
        return "CNV"
    if any(sorted(state[offset:offset + 2]) == [0, 2] for offset in (0, 2)):
        return "COPY_NEUTRAL_LOH"
    return "COPY_REMODELING"


def infer_gene(row: dict[str, str], depth: DepthEvidence, depth_scale: float,
               chi: float, min_gap: float) -> dict[str, str]:
    normalized_depth = depth.median / depth_scale
    groups = allele_groups(row)
    if not groups:
        raise ValueError("no usable allele fractions")
    slot_fractions = [finite_float(row.get(f"{slot}_copy_fraction")) for slot in SLOTS]
    fraction_total = sum(value for value in slot_fractions if value is not None)
    r_total = sum(value or 0.0 for value in slot_fractions[:2]) / fraction_total
    r_is_low = abs(r_total - chi) <= abs(r_total - (1.0 - chi))
    best = solve_dosage_milp(normalized_depth, chi, groups, r_is_low=r_is_low)
    best_bits = tuple(bit for dosage in best.state for bit in (dosage & 1, dosage >> 1))
    second = solve_dosage_milp(
        normalized_depth, chi, groups, r_is_low=r_is_low, exclude_bits=best_bits,
    )
    normal = solve_dosage_milp(
        normalized_depth, chi, groups, r_is_low=r_is_low, fixed_state=(1, 1, 1, 1),
    )
    gap = second.objective - best.objective
    event_support = normal.objective - best.objective
    event = classify_state(best.state)
    confidence = "HIGH" if gap >= min_gap else "AMBIGUOUS"
    reasons = []
    if depth.breadth < 0.8:
        confidence = "LOW_EVIDENCE"
        reasons.append("low_breadth")
    cross_source_shared = any(
        any(index < 2 for index in indices) and any(index >= 2 for index in indices)
        for _, indices, _ in groups
    )
    if event == "COPY_NEUTRAL_LOH" and cross_source_shared:
        confidence = "ASSIGNMENT_AMBIGUOUS"
        reasons.append("loh_allele_shared_across_sources")
    if len(groups) == 1 and event in {"COPY_NEUTRAL_LOH", "COPY_REMODELING"}:
        confidence = "UNIDENTIFIABLE"
        reasons.append("all_slots_same_allele")
    if event != "NORMAL" and event_support <= 0:
        confidence = "AMBIGUOUS"
        reasons.append("no_improvement_over_normal")
    return {
        "sample": row.get("sample", ""), "gene": row.get("gene", ""),
        "chi": f"{chi:.6f}", "raw_median_depth": f"{depth.median:.3f}",
        "normalized_depth": f"{normalized_depth:.5f}", "breadth": f"{depth.breadth:.5f}",
        **{f"{slot}_dosage": str(best.state[index]) for index, slot in enumerate(SLOTS)},
        "total_copies": str(sum(best.state)), "event": event, "confidence": confidence,
        "best_objective": f"{best.objective:.6f}",
        "second_objective": f"{second.objective:.6f}", "objective_gap": f"{gap:.6f}",
        "normal_objective": f"{normal.objective:.6f}", "event_support": f"{event_support:.6f}",
        "allele_groups": str(len(groups)), "reasons": ";".join(reasons),
    }


def global_chi(abundance_rows: list[dict[str, str]]) -> float:
    for row in abundance_rows:
        value = finite_float(row.get("global_chi"))
        if value is not None:
            return value
    raise ValueError("gene abundance report has no usable global_chi")


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--gene-abundance", required=True, type=Path)
    parser.add_argument("--compact-calls", required=True, type=Path)
    parser.add_argument("--merged-bam", required=True, type=Path)
    parser.add_argument("--gene-bed", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--breadth-min-depth", type=int, default=10)
    parser.add_argument("--min-gap", type=float, default=0.08)
    args = parser.parse_args()

    abundance_rows = read_tsv(args.gene_abundance)
    chi = global_chi(abundance_rows)
    calls = {row["gene"]: row for row in read_tsv(args.compact_calls)}
    intervals = read_intervals(args.gene_bed)
    depths = {}
    with pysam.AlignmentFile(args.merged_bam, "rb") as bam:
        for gene in calls.keys() & intervals.keys():
            evidence = bam_depth(bam, *intervals[gene], args.breadth_min_depth)
            if evidence is not None:
                depths[gene] = evidence
    reference_depths = [depths[gene].median for gene in CORE_GENES
                        if gene in depths and depths[gene].breadth >= 0.8 and depths[gene].median > 0]
    if len(reference_depths) < 3:
        raise SystemExit("fewer than three covered core genes; cannot normalize HLA depth")
    depth_scale = float(np.median(reference_depths))
    rows = []
    for gene, row in calls.items():
        if gene not in depths or gene not in intervals:
            continue
        try:
            rows.append(infer_gene(row, depths[gene], depth_scale, chi, args.min_gap))
        except ValueError:
            continue
    write_rows(args.out, rows)
    print(f"[CNV/LOH] wrote {args.out} ({len(rows)} genes; depth_scale={depth_scale:.3f})")


if __name__ == "__main__":
    main()