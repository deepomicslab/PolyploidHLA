#!/usr/bin/env python3
"""Apply truth-free class-II joint rescue rules to an ASM root.

Rules implemented here are deliberately narrow and default to high-mask
class-II loci:
  - HLA-DRB1 is anchored from the current HLA-DQB1 quartet through common
    DRB1-DQB1 linkage.
  - HLA-DPB1 high-number alleles are collapsed to common EM-supported
    alternatives when the current quartet contains likely rare artifacts.
    - HLA-DPB1 common low-frequency recipient-private alleles are recovered
        when the current high-mask recipient and donor quartets are identical.

The script writes per-gene calls.tsv overrides and re-aggregates final calls.
Truth is never read by this script.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from aggregate_calls import DEFAULT_GENES, allele_2field, main as aggregate_main  # noqa: E402


ASSIGNMENTS = ("R", "R", "D", "D")
DRB1_DQB1_LD = {
    "DQB1*02:01": "DRB1*03:01",
    "DQB1*02:02": "DRB1*07:01",
    "DQB1*02:82": "DRB1*07:01",
    "DQB1*02:109": "DRB1*03:01",
    "DQB1*03:01": "DRB1*04:01",
    "DQB1*06:02": "DRB1*15:01",
}


def read_tsv(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_quartet_text(text: str):
    return [allele for allele in (text or "").split(",") if allele]


def final_quartet(row):
    return [row.get(key, "NA") for key in ("R1_2field", "R2_2field", "D1_2field", "D2_2field")]


def gene_dir(asm_root: Path, sample: str, gene: str) -> Path:
    return asm_root / sample / gene.lower() / gene


def allele_number(allele: str) -> int:
    match = re.search(r"\*(\d+):", allele or "")
    return int(match.group(1)) if match else 999999


def allele_first_field(allele: str) -> str:
    match = re.search(r"\*(\d+):", allele or "")
    return match.group(1) if match else ""


def same_pair(left, right) -> bool:
    return sorted(left) == sorted(right)


def mean_mask(row) -> float:
    try:
        return float(row.get("mean_mask_fraction", "0") or 0.0)
    except ValueError:
        return 0.0


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
        return reader.fieldnames or [], list(reader)


def add_mapping(mapping, rows) -> None:
    for row in rows:
        allele = row.get("allele") or row.get("allele_2field") or row.get("call")
        if not allele or allele == "NA":
            continue
        mapping[allele_2field(allele)].append(allele)


def add_final_mapping(mapping, final_rows) -> None:
    for row in final_rows:
        for key in ("R1_full", "R2_full", "D1_full", "D2_full"):
            allele = row.get(key)
            if allele and allele != "NA":
                mapping[allele_2field(allele)].append(allele)


def build_full_allele_map(asm_root: Path, spechla_root: Path | None, sample: str, gene: str, final_rows):
    mapping = defaultdict(deque)
    add_final_mapping(mapping, final_rows)
    call_dir = gene_dir(asm_root, sample, gene)
    for name in (
        "calls.tsv",
        "calls.baseline.tsv",
        "calls.direct_gate_input.tsv",
        "calls.class2_joint_input.tsv",
    ):
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


def write_rescue_calls(call_dir: Path, lifted_quartet, source_row) -> None:
    calls = call_dir / "calls.tsv"
    if not calls.exists():
        raise FileNotFoundError(calls)
    _old_fields, old_rows = read_call_rows(calls)
    old_fraction_by_hap = {
        row.get("global_hap", ""): row.get("hap_fraction", "NA")
        for row in old_rows
    }
    backup = call_dir / "calls.class2_joint_input.tsv"
    if not backup.exists():
        shutil.copy2(calls, backup)
    fields = ["global_hap", "assignment", "allele", "hap_fraction", "class2_joint_rule", "class2_joint_reason"]
    rows = []
    for index, (assignment, allele) in enumerate(zip(ASSIGNMENTS, lifted_quartet), 1):
        hap = str(index)
        rows.append({
            "global_hap": hap,
            "assignment": assignment,
            "allele": allele,
            "hap_fraction": old_fraction_by_hap.get(hap, "NA"),
            "class2_joint_rule": source_row["rule"],
            "class2_joint_reason": source_row["reason"],
        })
    write_tsv(calls, fields, rows)


def read_tf_counts(spechla_root: Path, sample: str, gene: str):
    rows = []
    path = spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
    if not path.exists():
        return rows
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                rows.append({
                    "allele": row["allele_2field"],
                    "weight": float(row.get("em_weight") or 0.0),
                    "fraction": float(row.get("fraction") or row.get("em_frac") or 0.0),
                })
            except (KeyError, ValueError):
                continue
    return sorted(rows, key=lambda row: -row["weight"])


def drb1_from_dqb1(dqb1_row):
    if not dqb1_row:
        return None
    out = []
    for allele in final_quartet(dqb1_row):
        mapped = DRB1_DQB1_LD.get(allele)
        if not mapped:
            return None
        out.append(mapped)
    return out


def propose_drb1(row, dqb1_row, args):
    if row["gene"] != "HLA-DRB1" or mean_mask(row) < args.drb1_min_mask:
        return None
    current = final_quartet(row)
    candidate = drb1_from_dqb1(dqb1_row)
    if not candidate or candidate == current:
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "drdq_ld_anchor",
        "reason": f"mask={mean_mask(row):.4f};DQB1_anchor={','.join(final_quartet(dqb1_row))}",
        "current_2field": current,
        "new_2field": candidate,
    }


def propose_dpb1(row, spechla_root: Path, args):
    if row["gene"] != "HLA-DPB1" or mean_mask(row) < args.dpb1_min_mask:
        return None
    current = final_quartet(row)
    counts = read_tf_counts(spechla_root, row["sample"], row["gene"])
    common = [
        count_row["allele"]
        for count_row in counts
        if allele_number(count_row["allele"]) < args.dpb1_rare_cutoff
        and count_row["fraction"] >= args.dpb1_min_fraction
    ]
    used = Counter(allele for allele in current if allele_number(allele) < args.dpb1_rare_cutoff)
    candidate = list(current)
    rules = []
    reasons = []
    if common and any(allele_number(allele) >= args.dpb1_rare_cutoff for allele in candidate):
        changed = False
        for index, allele in enumerate(candidate):
            if allele_number(allele) < args.dpb1_rare_cutoff:
                continue
            replacement = None
            for common_allele in common[:args.dpb1_top_common]:
                if used[common_allele] < 2:
                    replacement = common_allele
                    break
            if replacement:
                candidate[index] = replacement
                used[replacement] += 1
                changed = True
        if changed:
            rules.append("dpb1_rare_collapse")
            reasons.append(
                f"rare_cutoff={args.dpb1_rare_cutoff};"
                f"min_fraction={args.dpb1_min_fraction};common={','.join(common[:args.dpb1_top_common])}"
            )

    minor = propose_dpb1_common_minor(candidate, counts, args)
    if minor:
        candidate, minor_reason = minor
        rules.append("dpb1_common_minor")
        reasons.append(minor_reason)

    absolute_common = propose_dpb1_absolute_common(candidate, counts, args)
    if absolute_common:
        candidate, absolute_reason = absolute_common
        rules.append("dpb1_absolute_common")
        reasons.append(absolute_reason)

    if not rules or candidate == current:
        return None
    return {
        "sample": row["sample"],
        "gene": row["gene"],
        "rule": "+".join(rules),
        "reason": f"mask={mean_mask(row):.4f};" + ";".join(reasons),
        "current_2field": current,
        "new_2field": candidate,
    }


def propose_dpb1_common_minor(current, counts, args):
    if not args.dpb1_common_minor:
        return None
    if not same_pair(current[:2], current[2:]) or len(set(current[:2])) != 2:
        return None
    current_set = set(current[:2])
    minor_rows = [
        count_row
        for count_row in counts[:args.dpb1_common_minor_top]
        if count_row["allele"] not in current_set
        and allele_number(count_row["allele"]) <= args.dpb1_common_minor_max_number
        and args.dpb1_common_minor_min_fraction <= count_row["fraction"] <= args.dpb1_common_minor_max_fraction
        and count_row["weight"] >= args.dpb1_common_minor_min_weight
    ]
    if not minor_rows:
        return None
    minor = minor_rows[0]
    current_fraction = {allele: 0.0 for allele in current_set}
    for count_row in counts:
        if count_row["allele"] in current_fraction:
            current_fraction[count_row["allele"]] = count_row["fraction"]
    same_first_field = [
        allele for allele in current_set
        if allele_first_field(allele) == allele_first_field(minor["allele"])
    ]
    if same_first_field:
        keep = same_first_field[0]
        replace = next(allele for allele in current[:2] if allele != keep)
    else:
        replace = min(current_set, key=lambda allele: current_fraction[allele])
    candidate = list(current)
    for index in (0, 1):
        if candidate[index] == replace:
            candidate[index] = minor["allele"]
            break
    if candidate == current:
        return None
    reason = (
        f"minor={minor['allele']}:{minor['fraction']:.4f};replace={replace};"
        f"max_number={args.dpb1_common_minor_max_number};"
        f"frac_range={args.dpb1_common_minor_min_fraction}-{args.dpb1_common_minor_max_fraction}"
    )
    return candidate, reason


def propose_dpb1_absolute_common(current, counts, args):
    if not args.dpb1_absolute_common:
        return None
    current_counts = Counter(current)
    if max(current_counts.values(), default=0) < 3:
        return None

    counts_by_allele = {row["allele"]: row for row in counts}
    missing_rows = [
        count_row
        for count_row in counts
        if count_row["allele"] not in current_counts
        and allele_number(count_row["allele"]) <= args.dpb1_absolute_common_max_number
        and count_row["weight"] >= args.dpb1_absolute_common_min_weight
        and count_row["fraction"] >= args.dpb1_absolute_common_min_fraction
    ]
    if not missing_rows:
        return None

    candidate = list(current)
    for missing in missing_rows:
        replace_index = None
        replace_metric = None
        for index, allele in enumerate(candidate):
            if current_counts[allele] < 3:
                continue
            observed = counts_by_allele.get(allele, {"weight": 0.0, "fraction": 0.0})
            if observed["weight"] > 0 and missing["weight"] < observed["weight"] * args.dpb1_absolute_common_min_ratio:
                continue
            metric = (0 if ASSIGNMENTS[index] == "D" else 1, observed["weight"], index)
            if replace_metric is None or metric < replace_metric:
                replace_metric = metric
                replace_index = index
        if replace_index is None:
            continue

        replaced = candidate[replace_index]
        candidate[replace_index] = missing["allele"]
        reason = (
            f"absolute_common={missing['allele']}:{missing['weight']:.1f}/{missing['fraction']:.4f};"
            f"replace={replaced};assignment={ASSIGNMENTS[replace_index]};"
            f"max_number={args.dpb1_absolute_common_max_number};"
            f"min_weight={args.dpb1_absolute_common_min_weight};"
            f"min_fraction={args.dpb1_absolute_common_min_fraction};"
            f"min_ratio={args.dpb1_absolute_common_min_ratio}"
        )
        return candidate, reason
    return None


def sample_names(asm_root: Path, requested):
    if requested:
        return sorted(set(requested))
    out = []
    for path in asm_root.iterdir():
        if path.is_dir() and (path / f"{path.name}.final_calls.tsv").exists():
            out.append(path.name)
    return sorted(out)


def proposals_for_sample(asm_root: Path, spechla_root: Path, sample: str, args):
    final_path = asm_root / sample / f"{sample}.final_calls.tsv"
    rows = read_tsv(final_path)
    by_gene = {row["gene"]: row for row in rows}
    proposals = []
    drb1 = by_gene.get("HLA-DRB1")
    if drb1:
        proposal = propose_drb1(drb1, by_gene.get("HLA-DQB1"), args)
        if proposal:
            proposals.append(proposal)
    dpb1 = by_gene.get("HLA-DPB1")
    if dpb1:
        proposal = propose_dpb1(dpb1, spechla_root, args)
        if proposal:
            proposals.append(proposal)
    return rows, proposals


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
        row["source"] = "class2-joint-rescue"
        tag = "class2_joint_rescue"
        warning = row.get("warning", "")
        row["warning"] = tag if not warning else f"{warning};{tag}"
    write_tsv(final_path, fields, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-asm-root", required=True, type=Path)
    parser.add_argument("--out-asm-root", type=Path, default=None)
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--g-group", required=True, type=Path)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--drb1-min-mask", type=float, default=0.40)
    parser.add_argument("--dpb1-min-mask", type=float, default=0.40)
    parser.add_argument("--dpb1-rare-cutoff", type=int, default=100)
    parser.add_argument("--dpb1-min-fraction", type=float, default=0.02)
    parser.add_argument("--dpb1-top-common", type=int, default=6)
    parser.add_argument("--dpb1-common-minor", action="store_true")
    parser.add_argument("--dpb1-common-minor-max-number", type=int, default=10)
    parser.add_argument("--dpb1-common-minor-min-fraction", type=float, default=0.005)
    parser.add_argument("--dpb1-common-minor-max-fraction", type=float, default=0.09)
    parser.add_argument("--dpb1-common-minor-min-weight", type=float, default=50.0)
    parser.add_argument("--dpb1-common-minor-top", type=int, default=12)
    parser.add_argument("--dpb1-absolute-common", action="store_true")
    parser.add_argument("--dpb1-absolute-common-max-number", type=int, default=10)
    parser.add_argument("--dpb1-absolute-common-min-weight", type=float, default=150.0)
    parser.add_argument("--dpb1-absolute-common-min-fraction", type=float, default=0.01)
    parser.add_argument("--dpb1-absolute-common-min-ratio", type=float, default=2.0)
    args = parser.parse_args()

    if args.in_place or args.dry_run:
        out_asm_root = args.in_asm_root
    else:
        if args.out_asm_root is None:
            raise SystemExit("--out-asm-root is required unless --in-place is set")
        out_asm_root = args.out_asm_root

    samples = sample_names(args.in_asm_root, args.sample)
    all_proposals = []
    final_rows_by_sample = {}
    for sample in samples:
        final_rows, proposals = proposals_for_sample(args.in_asm_root, args.spechla_root, sample, args)
        final_rows_by_sample[sample] = final_rows
        all_proposals.extend(proposals)

    manifest_fields = [
        "sample", "gene", "rule", "reason", "current_2field_quartet", "new_2field_quartet",
        "new_full_quartet", "output_calls",
    ]
    manifest_rows = []

    if not args.dry_run:
        out_asm_root.mkdir(parents=True, exist_ok=True)
        if not args.in_place:
            for sample in samples:
                copy_sample_tree(args.in_asm_root, out_asm_root, sample)

    accepted_keys = set()
    for proposal in all_proposals:
        sample = proposal["sample"]
        gene = proposal["gene"]
        final_rows = final_rows_by_sample[sample]
        mapping = build_full_allele_map(out_asm_root if not args.dry_run else args.in_asm_root,
                                        args.spechla_root, sample, gene, final_rows)
        lifted = lift_alleles(proposal["new_2field"], mapping)
        call_dir = gene_dir(out_asm_root if not args.dry_run else args.in_asm_root, sample, gene)
        if not args.dry_run:
            write_rescue_calls(call_dir, lifted, proposal)
        accepted_keys.add((sample, gene))
        manifest_rows.append({
            "sample": sample,
            "gene": gene,
            "rule": proposal["rule"],
            "reason": proposal["reason"],
            "current_2field_quartet": ",".join(proposal["current_2field"]),
            "new_2field_quartet": ",".join(proposal["new_2field"]),
            "new_full_quartet": ",".join(lifted),
            "output_calls": str(call_dir / "calls.tsv"),
        })

    if not args.dry_run:
        for sample in sorted({proposal["sample"] for proposal in all_proposals}):
            aggregate_sample(out_asm_root, sample, args.genes, args.g_group)
            mark_final_rows(out_asm_root / sample / f"{sample}.final_calls.tsv", accepted_keys)

    if args.manifest:
        write_tsv(args.manifest, manifest_fields, manifest_rows)
    print(f"samples\t{len(samples)}")
    print(f"accepted\t{len(manifest_rows)}")
    for row in manifest_rows:
        print(f"{row['sample']}\t{row['gene']}\t{row['rule']}\t{row['current_2field_quartet']} -> {row['new_2field_quartet']}")


if __name__ == "__main__":
    main()