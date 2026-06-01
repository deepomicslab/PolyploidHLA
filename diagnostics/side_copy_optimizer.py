#!/usr/bin/env python3
"""Truth-free conservative side/copy optimizer prototype.

This is a standalone proposal producer. It does not require truth and does not
modify pipeline outputs. When --quartet-summary is provided, truth columns are
used only for post-hoc validation of the emitted proposals.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_calls import load_g_group, normalize_for_display, overlap  # noqa: E402
from caller_free_4hap import parse_bed  # noqa: E402
from hla_polyphase_assemble import (  # noqa: E402
    PLOIDY,
    _score_combo,
    collect_blocks,
    compute_chim_penalty_vaf,
)
from rerank_multi_strategy import (  # noqa: E402
    EPS,
    STRATEGIES,
    build_locus,
    enum_one,
    enum_permute,
    evaluate_strategy,
    infer_set_label,
    quartet_key,
    score_quartet,
    two_field,
)

DEFAULT_PROFILE = (
    "C0_one_EM_AF_PHASE_REBAL",
)
GENES = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1")
CHIM_STRATEGIES = {"CHIM_perm", "CHIM_rebal"}
STRUCTURAL_STRATEGIES = {"DUP_interleave"}


def read_tsv(path: Path):
    if not path or not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_alleles(text: str):
    return [item for item in (text or "").split(",") if item]


def discover_samples_from_asm(asm_root: Path):
    out = {}
    for sample_dir in sorted(asm_root.iterdir()):
        if sample_dir.is_dir() and (sample_dir / f"{sample_dir.name}.final_calls.tsv").exists():
            out[sample_dir.name] = infer_set_label(sample_dir.name, sample_dir.parent.name)
    return out


def score_quartet_side_aware(quartet, truth_r, truth_d, gmap) -> int:
    pred_r = sorted([quartet[0], quartet[1]])
    pred_d = sorted([quartet[2], quartet[3]])
    truth_r_norm = normalize_for_display(truth_r, "2field", gmap)
    truth_d_norm = normalize_for_display(truth_d, "2field", gmap)
    return overlap(truth_r_norm, pred_r) + overlap(truth_d_norm, pred_d)


def current_quartet(row):
    return split_alleles(row["pred_R"]) + split_alleles(row["pred_D"])


def read_match_blocks(path: Path):
    if not path.exists():
        return []
    blocks = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                block_index = int(row["block"])
                local_hap = int(row["local_hap"]) - 1
                score = float(row["score"])
            except (KeyError, ValueError):
                continue
            if local_hap < 0 or local_hap >= PLOIDY:
                continue
            allele = two_field(row.get("allele", ""))
            if not allele:
                continue
            per_hap = blocks.setdefault(block_index, [dict() for _ in range(PLOIDY)])
            old = per_hap[local_hap].get(allele)
            if old is None or score > old:
                per_hap[local_hap][allele] = score
    return [blocks[index] for index in sorted(blocks)]


def assembly_chim_evidence(sample: str, gene: str, asm_root: Path, spechla_root: Path,
                           bed_path: Path, chim_weight_vaf: float):
    gene_lc = gene.lower()
    match_blocks = read_match_blocks(asm_root / sample / gene_lc / gene / "match_scores.tsv")
    if not match_blocks:
        return None
    contig = gene.replace("-", "_")
    start, end = parse_bed(str(bed_path), contig)
    if start is None:
        return None
    vcf_path = spechla_root / sample / f"{sample}.phased.{gene_lc}.vcf.gz"
    if not vcf_path.exists():
        return None
    try:
        blocks_info = collect_blocks(str(vcf_path), contig, start, end)
    except Exception:
        return None
    if not blocks_info:
        return None
    n_blocks = min(len(match_blocks), len(blocks_info))
    if n_blocks == 0:
        return None
    blocks_obs = [block[3] for block in blocks_info[:n_blocks]]
    perms24 = list(itertools.permutations(range(PLOIDY)))
    return {
        "match_blocks": match_blocks[:n_blocks],
        "blocks_obs": blocks_obs,
        "perms24": perms24,
        "chim_weight_vaf": chim_weight_vaf,
    }


def assembly_chim_score(quartet, evidence, chi_r: float):
    penalty_vaf = compute_chim_penalty_vaf(evidence["blocks_obs"], evidence["perms24"], chi_r)
    penalty = []
    for row in penalty_vaf:
        penalty.append([evidence["chim_weight_vaf"] * value for value in row])
    score, _perms = _score_combo(evidence["match_blocks"], tuple(quartet), evidence["perms24"], penalty)
    return score


def channel_guard(current_scores, alt_scores, channels, require_any_gain: bool):
    any_gain = False
    for channel in channels:
        cur_value = current_scores.get(channel)
        alt_value = alt_scores.get(channel)
        if cur_value is None or alt_value is None:
            continue
        if alt_value > cur_value + EPS[channel]:
            return False, any_gain
        if cur_value - alt_value >= EPS[channel]:
            any_gain = True
    if require_any_gain and not any_gain:
        return False, any_gain
    return True, any_gain


def evaluate_structural_strategy(name, locus):
    if name != "DUP_interleave":
        return None
    r1, r2, d1, d2 = locus.cur
    if r1 != r2 or d1 != d2 or r1 == d1:
        return None
    return {
        "alt_q": [r1, d1, r2, d2],
        "gain": 0.0,
        "any_channel_gain": "shape",
    }


def evaluate_chim_strategy(name, locus, args):
    evidence = assembly_chim_evidence(locus.sample, locus.gene, args.asm_root, args.spechla_root,
                                      args.bed, args.chim_weight_vaf)
    if evidence is None:
        return None
    if name == "CHIM_perm":
        candidates = enum_permute(locus.cur)
        min_gain = args.chim_perm_min_gain
        require_any_gain = False
    elif name == "CHIM_rebal":
        if locus.gene == "HLA-DPB1":
            return None
        candidates = [q for q in enum_one(locus.cur, locus.cur) if set(q).issubset(set(locus.cur))]
        min_gain = args.chim_rebal_min_gain
        require_any_gain = True
    else:
        return None
    current_assembly = assembly_chim_score(locus.cur, evidence, locus.chi_r)
    channels = ["em", "af", "phase"]
    current_scores = score_quartet(locus.cur, locus, "c0", channels)
    best = None
    seen = {quartet_key(locus.cur)}
    for quartet in candidates:
        key = quartet_key(quartet)
        if key in seen:
            continue
        seen.add(key)
        alt_assembly = assembly_chim_score(quartet, evidence, locus.chi_r)
        assembly_gain = alt_assembly - current_assembly
        if assembly_gain < min_gain:
            continue
        alt_scores = score_quartet(quartet, locus, "c0", channels)
        ok, any_gain = channel_guard(current_scores, alt_scores, channels, require_any_gain)
        if not ok:
            continue
        rank_gain = assembly_gain + sum(
            max(0.0, (current_scores.get(channel) or 0.0) - (alt_scores.get(channel) or 0.0))
            for channel in channels
        )
        if best is None or rank_gain > best["rank_gain"]:
            best = {
                "alt_q": list(quartet),
                "gain": rank_gain,
                "rank_gain": rank_gain,
                "assembly_current": current_assembly,
                "assembly_alt": alt_assembly,
                "assembly_gain": assembly_gain,
                "any_channel_gain": any_gain,
            }
    return best


def validate_manifest(summary_rows, proposals, gmap):
    by_key = {(row["sample"], row["gene"]): row for row in summary_rows}
    proposal_by_key = {(row["sample"], row["gene"]): split_alleles(row["proposed_quartet"]) for row in proposals}
    rows = []
    stats = defaultdict(Counter)
    base_total = 0
    new_total = 0
    for key, summary in sorted(by_key.items()):
        truth_r = split_alleles(summary["truth_R"])
        truth_d = split_alleles(summary["truth_D"])
        current = current_quartet(summary)
        proposed = proposal_by_key.get(key, current)
        base = score_quartet_side_aware(current, truth_r, truth_d, gmap)
        new = score_quartet_side_aware(proposed, truth_r, truth_d, gmap)
        base_total += base
        new_total += new
        if key not in proposal_by_key:
            continue
        delta = new - base
        verdict = "improve" if delta > 0 else "regress" if delta < 0 else "neutral"
        stats[("ALL", "ALL")][verdict] += 1
        stats[("ALL", "ALL")]["rows"] += 1
        stats[("ALL", "ALL")]["delta"] += delta
        stats[("gene", key[1])][verdict] += 1
        stats[("gene", key[1])]["rows"] += 1
        stats[("gene", key[1])]["delta"] += delta
        rows.append({
            "sample": key[0],
            "gene": key[1],
            "current_score": base,
            "proposed_score": new,
            "delta": delta,
            "verdict": verdict,
        })
    summary_out = [{
        "scope": "ALL",
        "name": "side_copy_optimizer",
        "rows": stats[("ALL", "ALL")]["rows"],
        "improve": stats[("ALL", "ALL")]["improve"],
        "regress": stats[("ALL", "ALL")]["regress"],
        "neutral": stats[("ALL", "ALL")]["neutral"],
        "delta": stats[("ALL", "ALL")]["delta"],
        "baseline": f"{base_total}/{len(summary_rows) * 4}",
        "optimized": f"{new_total}/{len(summary_rows) * 4}",
    }]
    for (scope, name), counter in sorted(stats.items()):
        if scope == "ALL":
            continue
        summary_out.append({
            "scope": scope,
            "name": name,
            "rows": counter["rows"],
            "improve": counter["improve"],
            "regress": counter["regress"],
            "neutral": counter["neutral"],
            "delta": counter["delta"],
            "baseline": "",
            "optimized": "",
        })
    return rows, summary_out


def choose_proposals_for_locus(locus, strategies_by_name, strategy_names, args):
    accepted = []
    for name in strategy_names:
        if name in CHIM_STRATEGIES:
            strategy = {"name": name, "channels": ["assembly_chim", "em", "af", "phase"], "search": name.replace("CHIM_", "")}
            result = evaluate_chim_strategy(name, locus, args)
        elif name in STRUCTURAL_STRATEGIES:
            strategy = {"name": name, "channels": ["shape"], "search": "duplicate_side_interleave"}
            result = evaluate_structural_strategy(name, locus)
        else:
            strategy = strategies_by_name[name]
            result = evaluate_strategy(strategy, locus)
        if not result or not result.get("alt_q"):
            continue
        proposed = result["alt_q"]
        accepted.append({
            "sample": locus.sample,
            "set": locus.set_label,
            "gene": locus.gene,
            "strategy": name,
            "current_quartet": ",".join(locus.cur),
            "proposed_quartet": ",".join(proposed),
            "channels": ",".join(strategy.get("channels", [])),
            "search": strategy.get("search", ""),
            "gain": f"{result.get('gain', 0.0):.6f}",
            "assembly_current": f"{result.get('assembly_current', '')}",
            "assembly_alt": f"{result.get('assembly_alt', '')}",
            "assembly_gain": f"{result.get('assembly_gain', '')}",
            "any_channel_gain": f"{result.get('any_channel_gain', '')}",
        })
        break
    return accepted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm-root", required=True, type=Path)
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--ref", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--imgt", required=True, type=Path)
    parser.add_argument("--out-manifest", required=True, type=Path)
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--gene", action="append", default=[])
    parser.add_argument("--strategy", action="append", default=[])
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--chim-weight-vaf", type=float, default=200.0)
    parser.add_argument("--chim-perm-min-gain", type=float, default=20.0)
    parser.add_argument("--chim-rebal-min-gain", type=float, default=0.0)
    parser.add_argument("--g-group", type=Path, default=SCRIPT_ROOT / "resources" / "spechla" / "db" / "HLA" / "hla_nom_g.txt")
    parser.add_argument("--quartet-summary", type=Path, default=None,
                        help="Optional validation summary; truth used only for post-hoc scoring.")
    parser.add_argument("--out-validation", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    args = parser.parse_args()

    strategies_by_name = {strategy["name"]: strategy for strategy in STRATEGIES}
    strategy_names = args.strategy or list(DEFAULT_PROFILE)
    unknown = [
        name for name in strategy_names
        if name not in strategies_by_name and name not in CHIM_STRATEGIES and name not in STRUCTURAL_STRATEGIES
    ]
    if unknown:
        raise SystemExit(f"unknown strategies: {','.join(unknown)}")

    samples = discover_samples_from_asm(args.asm_root)
    if args.sample:
        sample_set = set(args.sample)
        samples = {sample: label for sample, label in samples.items() if sample in sample_set}
    genes = args.gene or list(GENES)
    gmap = load_g_group(args.g_group) if args.g_group.exists() else {}

    proposal_rows = []
    skipped = Counter()
    for sample, set_label in sorted(samples.items()):
        for gene in genes:
            locus = build_locus(sample, gene, set_label, args.asm_root, args.spechla_root,
                                args.ref, args.bed, args.imgt, gmap, args.top_k)
            if locus is None:
                skipped[gene] += 1
                continue
            proposal_rows.extend(choose_proposals_for_locus(locus, strategies_by_name, strategy_names, args))

    fields = [
        "sample", "set", "gene", "strategy", "current_quartet", "proposed_quartet",
        "channels", "search", "gain", "assembly_current", "assembly_alt", "assembly_gain",
        "any_channel_gain",
    ]
    write_tsv(args.out_manifest, fields, proposal_rows)

    if args.quartet_summary:
        summary_rows = read_tsv(args.quartet_summary)
        if args.sample:
            sample_set = set(args.sample)
            summary_rows = [row for row in summary_rows if row["sample"] in sample_set]
        if args.gene:
            gene_set = set(args.gene)
            summary_rows = [row for row in summary_rows if row["gene"] in gene_set]
        validation_rows, summary_out = validate_manifest(summary_rows, proposal_rows, gmap)
        if args.out_validation:
            write_tsv(args.out_validation, ["sample", "gene", "current_score", "proposed_score", "delta", "verdict"], validation_rows)
        if args.out_summary:
            write_tsv(args.out_summary, ["scope", "name", "rows", "improve", "regress", "neutral", "delta", "baseline", "optimized"], summary_out)

    print(f"samples\t{len(samples)}")
    print(f"proposals\t{len(proposal_rows)}")
    for gene in sorted(skipped):
        print(f"skipped_locus\t{gene}\t{skipped[gene]}")
    if args.out_summary and args.out_summary.exists():
        for row in read_tsv(args.out_summary):
            if row.get("scope") == "ALL":
                print("summary\t" + "\t".join(row.get(col, "") for col in ["rows", "improve", "regress", "neutral", "delta", "baseline", "optimized"]))


if __name__ == "__main__":
    main()
