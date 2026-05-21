#!/usr/bin/env python3
"""Run DRB345 add-on typing for ABC real sample sets.

This is intentionally a thin driver around type_drb345.py. It reuses the
existing polyphase outputs for each sample, especially <sample>.map_database.bam,
so the competitive DB alignment step is not repeated.
"""
from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REALSETS = Path("/data2/wangxuedong/polyploid-hla-realsets")
DEFAULT_SUMMARY = Path(
    "/data6/wangxuedong/polyploid_hla/real_sample_summaries/real_sample_simple_gene_copy_fraction_fit_20260518.tsv"
)
DEFAULT_SPECHLA_ROOT = REALSETS / "spechla_out_abc_realsets_rescue_20260512"
DEFAULT_ASM_ROOT = REALSETS / "asm_v2_abc_realsets_class2_joint_common_minor_safe_20260513"
DEFAULT_SETS = ("set-a", "set-b", "set-c")


def read_samples(summary: Path, wanted_sets: set[str], wanted_samples: set[str] | None) -> list[tuple[str, str]]:
    seen: set[tuple[str, str]] = set()
    out: list[tuple[str, str]] = []
    with summary.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_set = row.get("sample_set", "")
            sample = row.get("sample", "")
            if sample_set not in wanted_sets:
                continue
            if wanted_samples is not None and sample not in wanted_samples:
                continue
            key = (sample_set, sample)
            if key not in seen:
                seen.add(key)
                out.append(key)
    return out


def final_has_drb345(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return any(row.get("gene") == "HLA-DRB345" for row in reader)


def quote_cmd(cmd: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def build_command(args: argparse.Namespace, sample: str) -> list[str]:
    fq_dir = args.spechla_root / sample
    return [
        args.python,
        str(args.type_drb345),
        "--sample",
        sample,
        "--fq-dir",
        str(fq_dir),
        "--db-bam",
        str(fq_dir / f"{sample}.map_database.bam"),
        "--asm-root",
        str(args.asm_root),
        "--final-calls",
        str(args.asm_root / sample / f"{sample}.final_calls.tsv"),
        "--threads",
        str(args.threads),
        "--subs-per-2field",
        str(args.subs_per_2field),
        "--top-per-locus",
        str(args.top_per_locus),
        "--db-min-as-frac",
        str(args.db_min_as_frac),
        "--remap-min-as-frac",
        str(args.remap_min_as_frac),
        "--evidence-k",
        str(args.evidence_k),
        "--min-locus-unique-frac",
        str(args.min_locus_unique_frac),
        "--drb1-untrusted-mask",
        str(args.drb1_untrusted_mask),
    ]


def missing_inputs(args: argparse.Namespace, sample: str) -> list[Path]:
    fq_dir = args.spechla_root / sample
    return [
        path
        for path in (
            fq_dir / f"{sample}.uniq.R1.fq.gz",
            fq_dir / f"{sample}.uniq.R2.fq.gz",
            fq_dir / f"{sample}.map_database.bam",
            args.asm_root / sample / f"{sample}.final_calls.tsv",
        )
        if not path.exists()
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--spechla-root", type=Path, default=DEFAULT_SPECHLA_ROOT)
    parser.add_argument("--asm-root", type=Path, default=DEFAULT_ASM_ROOT)
    parser.add_argument("--sets", nargs="+", default=list(DEFAULT_SETS), help="sample sets to run, default: set-a set-b set-c")
    parser.add_argument("--samples", nargs="+", default=None, help="optional explicit sample names")
    parser.add_argument("--type-drb345", type=Path, default=SCRIPT_DIR / "type_drb345.py")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--subs-per-2field", type=int, default=5)
    parser.add_argument("--top-per-locus", type=int, default=8)
    parser.add_argument("--db-min-as-frac", type=float, default=0.90)
    parser.add_argument("--remap-min-as-frac", type=float, default=0.95)
    parser.add_argument("--evidence-k", type=int, default=71)
    parser.add_argument("--min-locus-unique-frac", type=float, default=-1.0)
    parser.add_argument("--drb1-untrusted-mask", type=float, default=0.50)
    parser.add_argument("--force", action="store_true", help="rerun even when calls.tsv and final HLA-DRB345 row already exist")
    parser.add_argument("--dry-run", action="store_true", help="print commands without running")
    parser.add_argument("--keep-going", action="store_true", help="continue after a sample fails")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    wanted_samples = set(args.samples) if args.samples else None
    samples = read_samples(args.summary, set(args.sets), wanted_samples)
    if not samples:
        print("[drb345-batch] no samples selected", file=sys.stderr)
        return 2

    print(f"[drb345-batch] selected {len(samples)} samples from {','.join(args.sets)}")
    print(f"[drb345-batch] spechla_root={args.spechla_root}")
    print(f"[drb345-batch] asm_root={args.asm_root}")

    skipped = 0
    failed = 0
    for sample_set, sample in samples:
        call_path = args.asm_root / sample / "hla-drb345" / "HLA-DRB345" / "calls.tsv"
        final_path = args.asm_root / sample / f"{sample}.final_calls.tsv"
        if not args.force and call_path.exists() and final_has_drb345(final_path):
            skipped += 1
            print(f"[skip] {sample_set} {sample}: DRB345 already present")
            continue

        missing = missing_inputs(args, sample)
        if missing:
            failed += 1
            print(f"[missing] {sample_set} {sample}: " + ", ".join(str(path) for path in missing), file=sys.stderr)
            if not args.keep_going:
                return 1
            continue

        cmd = build_command(args, sample)
        if args.dry_run:
            print(quote_cmd(cmd))
            continue

        print(f"[run] {sample_set} {sample}")
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            failed += 1
            print(f"[failed] {sample}: exit={result.returncode}", file=sys.stderr)
            if not args.keep_going:
                return result.returncode

    print(f"[drb345-batch] done: selected={len(samples)} skipped={skipped} failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())