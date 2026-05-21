#!/usr/bin/env python3
"""Truth-free offline reranking experiments over existing PolyploidHLA outputs.

This script deliberately reuses existing per-gene calls/support files and writes
only lightweight trial final_calls.tsv files. It does not rerun read binning,
deduplication, database mapping, assembly, or EM remapping.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from aggregate_calls import DEFAULT_G_GROUP, SLOTS, allele_2field, allele_g_group, load_g_group  # noqa: E402
from evaluate_calls import load_truth, normalize_for_display, norm_allele, overlap, overlap_g_group_truth_resolution  # noqa: E402
from generate_summaries import read_chi_r as summary_read_chi_r  # noqa: E402


OFFICIAL = Path("/data2/wangxuedong/polyploid-hla-realsets/asm_v2_abc_realsets_class2_joint_common_minor_safe_20260513")
HYBRID = Path("/data2/wangxuedong/polyploid-hla-realsets/asm_v2_abc_realsets_dpb1_readbin_rescue_reuse_smoke_abscommon_20260514")
SETD = Path("/data2/wangxuedong/polyploid-hla-realsets/asm_v2_set_d_current_20260514")
REALSETS = Path("/data2/wangxuedong/polyploid-hla-realsets")

SUPPORT_ROOTS = [
    REALSETS / "spechla_out_abc_realsets_dpb1_readbin_rescue_reuse_smoke_20260514",
    REALSETS / "spechla_out_abc_realsets_class2_readbin_rescue_reuse_20260513",
    REALSETS / "spechla_out_abc_realsets_rescue_20260512",
    REALSETS / "spechla_out_abc_realsets_local_20260511",
    REALSETS / "spechla_out_set_d_current_20260514",
    REALSETS / "proportion_support_rerun_20260518",
]

TRUTH = {
    "set-a": ROOT / "truth" / "truth_typing-set-a.tsv",
    "set-b": ROOT / "truth" / "truth_typing-set-b.tsv",
    "set-c": ROOT / "truth" / "truth_typing-set-c.tsv",
    "set-d": ROOT / "truth" / "truth_typing-set-d.tsv",
}

CLASS2 = {"HLA-DRB1", "HLA-DQB1", "HLA-DPB1"}
STRATEGIES = (
    "current",
    "class2_input",
    "class2_baseline",
    "dpb1_input",
    "dpb1_baseline",
    "dose_fit",
    "dose_fit_guarded",
)


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
    prefix = sample.split("-", 1)[0]
    if prefix == "248444":
        return "set-d"
    value = int(prefix)
    if 267015 <= value <= 267019:
        return "set-a"
    if 267020 <= value <= 267024:
        return "set-b"
    if 267025 <= value <= 267029:
        return "set-c"
    return "unknown"


def final_path(root: Path, sample: str) -> Path:
    return root / sample / f"{sample}.final_calls.tsv"


def gene_dir(root: Path, sample: str, gene: str) -> Path:
    return root / sample / gene.lower() / gene


def input_samples() -> list[tuple[str, Path, str]]:
    hybrid_samples = {
        path.name for path in HYBRID.iterdir()
        if path.is_dir() and final_path(HYBRID, path.name).exists()
    }
    out: list[tuple[str, Path, str]] = []
    for path in sorted(OFFICIAL.iterdir()):
        if path.is_dir() and final_path(OFFICIAL, path.name).exists():
            root = HYBRID if path.name in hybrid_samples else OFFICIAL
            label = "abc_hybrid" if root == HYBRID else "abc_official"
            out.append((path.name, root, label))
    for path in sorted(SETD.iterdir()):
        if path.is_dir() and final_path(SETD, path.name).exists():
            out.append((path.name, SETD, "setd"))
    return out


def quartet_from_final(row: dict[str, str]) -> list[str]:
    return [row.get(f"{slot}_full") or row.get(f"{slot}_2field") or "NA" for slot in SLOTS]


def quartet_2field(quartet: Iterable[str]) -> list[str]:
    return [allele_2field(allele) for allele in quartet]


def quartet_from_call_file(path: Path) -> Optional[list[str]]:
    rows = read_tsv(path)
    if not rows:
        return None
    rows.sort(key=lambda row: int(row.get("global_hap", "0") or 0))
    recipient = [row.get("allele") or row.get("allele_2field") or "NA" for row in rows if row.get("assignment") == "R"]
    donor = [row.get("allele") or row.get("allele_2field") or "NA" for row in rows if row.get("assignment") == "D"]
    if len(recipient) < 2 or len(donor) < 2:
        return None
    return recipient[:2] + donor[:2]


def float_or_none(value: object) -> Optional[float]:
    try:
        if value in (None, "", "NA"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def mean_mask(row: dict[str, str]) -> float:
    return float_or_none(row.get("mean_mask_fraction")) or 0.0


def support_path(sample: str, gene: str) -> Optional[Path]:
    for root in SUPPORT_ROOTS:
        for path in (
            root / sample / "em_refine" / f"{gene}.tf_counts.tsv",
            root / sample / f"{gene}.tf_counts.tsv",
            root / sample / gene.lower().replace("hla-", "") / f"{gene}.tf_counts.tsv",
        ):
            if path.exists():
                return path
    return None


def read_support(sample: str, gene: str) -> dict[str, float]:
    path = support_path(sample, gene)
    counts: dict[str, float] = {}
    if path is None:
        return counts
    for row in read_tsv(path):
        allele = row.get("allele_2field") or allele_2field(row.get("allele", ""))
        value = float_or_none(row.get("em_weight"))
        if value is None:
            value = float_or_none(row.get("fraction"))
        if allele and allele != "NA" and value is not None:
            counts[allele] = max(counts.get(allele, 0.0), value)
    return counts


def add_callfile_candidates(candidates: set[str], root: Path, sample: str, gene: str) -> None:
    call_dir = gene_dir(root, sample, gene)
    for name in ("calls.tsv", "calls.class2_joint_input.tsv", "calls.baseline.tsv", "calls.direct_gate_input.tsv"):
        quartet = quartet_from_call_file(call_dir / name)
        if quartet:
            candidates.update(quartet_2field(quartet))


def callfile_quartets(root: Path, sample: str, gene: str) -> list[tuple[str, str, str, str]]:
    out = []
    call_dir = gene_dir(root, sample, gene)
    for name in ("calls.tsv", "calls.class2_joint_input.tsv", "calls.baseline.tsv", "calls.direct_gate_input.tsv"):
        quartet = quartet_from_call_file(call_dir / name)
        if quartet:
            out.append(tuple(quartet_2field(quartet)))
    return out


def expected_fraction(quartet: tuple[str, str, str, str], allele: str, chi_r: float) -> float:
    return (quartet[:2].count(allele) * chi_r / 2.0) + (quartet[2:].count(allele) * (1.0 - chi_r) / 2.0)


def residual(quartet: tuple[str, str, str, str], fractions: dict[str, float], chi_r: float) -> float:
    keys = set(fractions) | set(quartet)
    return sum(abs(fractions.get(allele, 0.0) - expected_fraction(quartet, allele, chi_r)) for allele in keys)


def best_dose_quartet(row: dict[str, str], root: Path, sample: str, gene: str) -> tuple[Optional[list[str]], str]:
    chi_r = summary_read_chi_r(sample)
    if chi_r is None or not 0.0 < chi_r < 1.0:
        return None, "no_chi_r"
    if gene not in CLASS2 and chi_r > 0.20:
        return None, f"skip_nonclass2_chi={chi_r:.4f}"
    if gene in CLASS2 and mean_mask(row) < 0.15 and chi_r > 0.20:
        return None, f"skip_lowmask_chi={chi_r:.4f}"
    support = read_support(sample, gene)
    if not support:
        return None, "no_support"
    total = sum(support.values())
    if total <= 0:
        return None, "zero_support"
    fractions = {allele: value / total for allele, value in support.items()}
    current = tuple(quartet_2field(quartet_from_final(row)))
    candidates = {allele for allele, frac in sorted(fractions.items(), key=lambda item: -item[1])[:8] if frac >= 0.003}
    candidates.update(current)
    add_callfile_candidates(candidates, root, sample, gene)
    candidates = {allele for allele in candidates if allele and allele != "NA"}
    if len(candidates) < 2 or len(candidates) > 10:
        return None, f"candidate_count={len(candidates)}"
    names = sorted(candidates)
    candidate_quartets: set[tuple[str, str, str, str]] = {current}
    candidate_quartets.update(callfile_quartets(root, sample, gene))
    current_list = list(current)
    for left, right in ((0, 2), (0, 3), (1, 2), (1, 3), (0, 1), (2, 3)):
        trial = list(current_list)
        trial[left], trial[right] = trial[right], trial[left]
        candidate_quartets.add(tuple(trial))
    for index in range(4):
        for allele in names:
            trial = list(current_list)
            trial[index] = allele
            candidate_quartets.add(tuple(trial))

    best: Optional[tuple[tuple[str, str, str, str], float]] = None
    for quartet in candidate_quartets:
        score = residual(quartet, fractions, chi_r)
        # Avoid fully erasing the minor source unless the current call already does.
        if len(set(quartet[:2])) == 1 and quartet[0] in quartet[2:]:
            score += 0.02
        if best is None or score < best[1]:
            best = (quartet, score)
    if best is None:
        return None, "no_best"
    current_score = residual(current, fractions, chi_r)
    return list(best[0]), f"dose_fit;current={current_score:.4f};best={best[1]:.4f};n={len(candidates)};q={len(candidate_quartets)}"


def set_quartet(row: dict[str, str], quartet: list[str], gmap, source_tag: str, reason: str) -> dict[str, str]:
    out = dict(row)
    high_mask = mean_mask(out) >= 0.15 or out.get("report_level") == "2-field"
    for slot, allele in zip(SLOTS, quartet):
        out[f"{slot}_full"] = allele
        out[f"{slot}_2field"] = allele_2field(allele)
        out[f"{slot}_g_group"] = allele_g_group(allele, gmap)
        out[f"{slot}_report"] = allele_2field(allele) if high_mask else allele
    out["source"] = source_tag
    warning = out.get("warning", "")
    extra = f"offline_rerank:{reason}"
    out["warning"] = extra if not warning else f"{warning};{extra}"
    return out


def apply_strategy(row: dict[str, str], root: Path, sample: str, strategy: str, gmap) -> tuple[dict[str, str], str]:
    gene = row["gene"]
    if strategy == "current" or gene == "HLA-DRB345":
        return dict(row), "current"
    call_dir = gene_dir(root, sample, gene)
    if strategy == "class2_input" and gene in CLASS2 and mean_mask(row) >= 0.15:
        quartet = quartet_from_call_file(call_dir / "calls.class2_joint_input.tsv")
        if quartet:
            return set_quartet(row, quartet, gmap, "offline-class2-input", "class2_input"), "class2_input"
    if strategy == "class2_baseline" and gene in CLASS2 and mean_mask(row) >= 0.15:
        quartet = quartet_from_call_file(call_dir / "calls.baseline.tsv")
        if quartet:
            return set_quartet(row, quartet, gmap, "offline-class2-baseline", "class2_baseline"), "class2_baseline"
    if strategy == "dpb1_input" and gene == "HLA-DPB1" and mean_mask(row) >= 0.15:
        quartet = quartet_from_call_file(call_dir / "calls.class2_joint_input.tsv")
        if quartet:
            return set_quartet(row, quartet, gmap, "offline-dpb1-input", "dpb1_input"), "dpb1_input"
    if strategy == "dpb1_baseline" and gene == "HLA-DPB1" and mean_mask(row) >= 0.15:
        quartet = quartet_from_call_file(call_dir / "calls.baseline.tsv")
        if quartet:
            return set_quartet(row, quartet, gmap, "offline-dpb1-baseline", "dpb1_baseline"), "dpb1_baseline"
    if strategy in {"dose_fit", "dose_fit_guarded"}:
        quartet, reason = best_dose_quartet(row, root, sample, gene)
        if quartet:
            current = quartet_2field(quartet_from_final(row))
            changed = quartet != current
            if strategy == "dose_fit_guarded":
                # Conservative gate: only use clear improvements and avoid touching high-confidence full calls.
                parts = dict(item.split("=", 1) for item in reason.split(";") if "=" in item)
                current_score = float(parts.get("current", "0"))
                best_score = float(parts.get("best", "0"))
                if not changed or current_score - best_score < 0.08 or (mean_mask(row) < 0.15 and gene not in CLASS2):
                    return dict(row), "dose_fit_guarded_reject:" + reason
            return set_quartet(row, quartet, gmap, "offline-dose-fit", reason), reason
        return dict(row), reason
    return dict(row), "no_change"


def pred_values(row: dict[str, str], side: str) -> list[str]:
    prefix = ("R1", "R2") if side == "PATIENT" else ("D1", "D2")
    return [row.get(f"{slot}_report") or row.get(f"{slot}_full") or "NA" for slot in prefix]


def score_rows(rows_by_gene: dict[str, dict[str, str]], truth, gmap) -> dict[str, tuple[int, int]]:
    scores = {"2field": [0, 0], "g_group": [0, 0]}
    for side in ("PATIENT", "DONOR"):
        for gene, truth_vals in truth[side].items():
            row = rows_by_gene.get(gene)
            if row is None:
                continue
            preds = pred_values(row, side)
            truth_2 = normalize_for_display(truth_vals, "2field", gmap)
            pred_2 = sorted(norm_allele(value, "2field", gmap) for value in preds)
            scores["2field"][0] += overlap(truth_2, pred_2)
            scores["2field"][1] += len(truth_2)
            scores["g_group"][0] += overlap_g_group_truth_resolution(truth_vals, preds, gmap)
            scores["g_group"][1] += len(truth_vals)
    return {key: (value[0], value[1]) for key, value in scores.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=ROOT / "diagnostics" / "offline_rerank_existing_20260518")
    parser.add_argument("--strategy", action="append", choices=STRATEGIES, default=[])
    parser.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP)
    args = parser.parse_args()

    strategies = args.strategy or list(STRATEGIES)
    gmap = load_g_group(args.g_group)
    truth = {name: load_truth(path) for name, path in TRUTH.items()}
    fields: Optional[list[str]] = None
    manifest_rows = []
    summary_rows = []
    per_strategy_totals = {strategy: {"2field": [0, 0], "g_group": [0, 0], "changed": 0} for strategy in strategies}

    for strategy in strategies:
        for sample, root, root_label in input_samples():
            final_rows = read_tsv(final_path(root, sample))
            if not final_rows:
                continue
            fields = fields or list(final_rows[0].keys())
            new_rows = []
            changed = 0
            for row in final_rows:
                new_row, reason = apply_strategy(row, root, sample, strategy, gmap)
                if quartet_2field(quartet_from_final(new_row)) != quartet_2field(quartet_from_final(row)):
                    changed += 1
                    manifest_rows.append({
                        "strategy": strategy,
                        "sample": sample,
                        "sample_set": sample_set(sample),
                        "gene": row["gene"],
                        "root": root_label,
                        "old_quartet": ",".join(quartet_2field(quartet_from_final(row))),
                        "new_quartet": ",".join(quartet_2field(quartet_from_final(new_row))),
                        "reason": reason,
                    })
                new_rows.append(new_row)
            out_path = args.out_dir / strategy / sample / f"{sample}.final_calls.tsv"
            write_tsv(out_path, fields, new_rows)
            rows_by_gene = {row["gene"]: row for row in new_rows}
            scores = score_rows(rows_by_gene, truth[sample_set(sample)], gmap)
            for level in ("2field", "g_group"):
                per_strategy_totals[strategy][level][0] += scores[level][0]
                per_strategy_totals[strategy][level][1] += scores[level][1]
            per_strategy_totals[strategy]["changed"] += changed

    current_ok = None
    for strategy in strategies:
        total = per_strategy_totals[strategy]
        if strategy == "current":
            current_ok = total["2field"][0]
        delta = "NA" if current_ok is None else str(total["2field"][0] - current_ok)
        summary_rows.append({
            "strategy": strategy,
            "changed_gene_rows": total["changed"],
            "two_field": f"{total['2field'][0]}/{total['2field'][1]}",
            "two_field_delta_vs_current": delta,
            "g_group": f"{total['g_group'][0]}/{total['g_group'][1]}",
        })

    write_tsv(args.out_dir / "strategy_summary.tsv", ["strategy", "changed_gene_rows", "two_field", "two_field_delta_vs_current", "g_group"], summary_rows)
    write_tsv(args.out_dir / "manifest.tsv", ["strategy", "sample", "sample_set", "gene", "root", "old_quartet", "new_quartet", "reason"], manifest_rows)
    for row in summary_rows:
        print("\t".join(str(row[key]) for key in ("strategy", "changed_gene_rows", "two_field", "two_field_delta_vs_current", "g_group")))


if __name__ == "__main__":
    main()