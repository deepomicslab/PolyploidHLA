#!/usr/bin/env python3
"""Offline diagnostics for three rescue ideas on current real-sample calls.

The script reuses existing final calls, per-gene FASTQs, tf_counts, and DRB345
outputs. It does not rerun binning, mapping, assembly, or EM.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(ROOT))

import evaluate_calls as ev  # noqa: E402
from aggregate_calls import DEFAULT_G_GROUP, allele_2field  # noqa: E402
from generate_summaries import read_chi_r  # noqa: E402


REALSETS = Path("/data2/wangxuedong/polyploid-hla-realsets")
OFFICIAL = REALSETS / "asm_v2_abc_realsets_class2_joint_common_minor_safe_20260513"
HYBRID = REALSETS / "asm_v2_abc_realsets_dpb1_readbin_rescue_reuse_smoke_abscommon_20260514"
SETD = REALSETS / "asm_v2_set_d_current_20260514"

ASM_ROOTS = [
    REALSETS / "asm_v2_abc_realsets_local_20260511",
    REALSETS / "asm_v2_abc_realsets_rescue_20260512",
    REALSETS / "asm_v2_abc_realsets_direct_gate_20260513",
    REALSETS / "asm_v2_abc_realsets_direct_gate_safe_20260513",
    REALSETS / "asm_v2_abc_realsets_class2_joint_safe_20260513",
    REALSETS / "asm_v2_abc_realsets_class2_joint_common_minor_safe_20260513",
    REALSETS / "asm_v2_abc_realsets_dpb1_readbin_rescue_reuse_smoke_abscommon_20260514",
]

SPECHLA_ROOTS = [
    REALSETS / "spechla_out_abc_realsets_dpb1_readbin_rescue_reuse_smoke_20260514",
    REALSETS / "spechla_out_abc_realsets_class2_readbin_rescue_reuse_20260513",
    REALSETS / "spechla_out_abc_realsets_rescue_20260512",
    REALSETS / "spechla_out_abc_realsets_local_20260511",
    REALSETS / "spechla_out_set_d_current_20260514",
]

TRUTH = {
    "set-a": ROOT / "truth" / "truth_typing-set-a.tsv",
    "set-b": ROOT / "truth" / "truth_typing-set-b.tsv",
    "set-c": ROOT / "truth" / "truth_typing-set-c.tsv",
    "set-d": ROOT / "truth" / "truth_typing-set-d.tsv",
}

SLOTS = ("R1", "R2", "D1", "D2")


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


def current_roots() -> list[tuple[str, Path, str]]:
    hybrid_samples = {path.name for path in HYBRID.iterdir() if path.is_dir() and final_path(HYBRID, path.name).exists()}
    out = []
    for path in sorted(OFFICIAL.iterdir()):
        if path.is_dir() and final_path(OFFICIAL, path.name).exists():
            root = HYBRID if path.name in hybrid_samples else OFFICIAL
            out.append((path.name, root, "abc_hybrid" if root == HYBRID else "abc_official"))
    for path in sorted(SETD.iterdir()):
        if path.is_dir() and final_path(SETD, path.name).exists():
            out.append((path.name, SETD, "setd"))
    return out


def load_final(path: Path) -> dict[str, dict[str, str]]:
    return {row["gene"]: row for row in read_tsv(path)}


def quartet(row: dict[str, str]) -> list[str]:
    return [row.get(f"{slot}_report") or row.get(f"{slot}_2field") or row.get(f"{slot}_full") or "NA" for slot in SLOTS]


def side_pred(row: dict[str, str], side: str) -> list[str]:
    slots = ("R1", "R2") if side == "PATIENT" else ("D1", "D2")
    return [row.get(f"{slot}_report") or row.get(f"{slot}_2field") or row.get(f"{slot}_full") or "NA" for slot in slots]


def score_gene(row: dict[str, str], truth_vals: list[str], side: str, gmap) -> int:
    truth_2 = ev.normalize_for_display(truth_vals, "2field", gmap)
    pred_2 = sorted(ev.norm_allele(value, "2field", gmap) for value in side_pred(row, side))
    return ev.overlap(truth_2, pred_2)


def score_sample(rows: dict[str, dict[str, str]], truth, gmap) -> tuple[int, int, int, int]:
    ok2 = total2 = okg = totalg = 0
    for side in ("PATIENT", "DONOR"):
        for gene, truth_vals in truth[side].items():
            row = rows.get(gene)
            if row is None:
                continue
            ok2 += score_gene(row, truth_vals, side, gmap)
            total2 += len(truth_vals)
            okg += ev.overlap_g_group_truth_resolution(truth_vals, side_pred(row, side), gmap)
            totalg += len(truth_vals)
    return ok2, total2, okg, totalg


def score_all(sample_rows: dict[str, dict[str, dict[str, str]]], truth_by_set, gmap) -> tuple[int, int, int, int]:
    total = [0, 0, 0, 0]
    for sample, rows in sample_rows.items():
        scores = score_sample(rows, truth_by_set[sample_set(sample)], gmap)
        for idx, value in enumerate(scores):
            total[idx] += value
    return tuple(total)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def iter_fastq(path: Path) -> Iterable[str]:
    if not path.exists():
        return
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip().upper()
            handle.readline(); handle.readline()
            if seq:
                yield seq


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def iter_kmers(seq: str, k: int):
    seq = seq.upper()
    for start in range(0, len(seq) - k + 1):
        kmer = seq[start:start + k]
        if "N" not in kmer:
            yield kmer


def read_fasta(path: Path):
    name = ""
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name and seq:
                    yield name, "".join(seq).upper().replace("-", "")
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name and seq:
        yield name, "".join(seq).upper().replace("-", "")


def allele_from_header(header: str) -> Optional[str]:
    match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", header)
    return allele_2field(match.group(1)) if match else None


def build_family_unique_kmers(exon_fasta: Path, k: int) -> tuple[dict[str, str], Counter[str]]:
    owners: dict[str, set[str]] = defaultdict(set)
    for header, seq in read_fasta(exon_fasta):
        family = allele_from_header(header)
        if not family:
            continue
        for kmer in set(iter_kmers(seq, k)) | set(iter_kmers(revcomp(seq), k)):
            owners[kmer].add(family)
    unique = {}
    counts = Counter()
    for kmer, families in owners.items():
        if len(families) == 1:
            family = next(iter(families))
            unique[kmer] = family
            counts[family] += 1
    return unique, counts


def spechla_dir(sample: str) -> Optional[Path]:
    for root in SPECHLA_ROOTS:
        path = root / sample
        if path.exists():
            return path
    return None


def fastq_support(sample: str, gene: str, unique_map: dict[str, str], k: int) -> Counter[str]:
    sample_dir = spechla_dir(sample)
    if sample_dir is None:
        return Counter()
    short = gene.split("-", 1)[1]
    hits = Counter()
    for fq in (sample_dir / f"{short}.R1.fq.gz", sample_dir / f"{short}.R2.fq.gz"):
        for seq in iter_fastq(fq):
            families = {unique_map[kmer] for kmer in iter_kmers(seq, k) if kmer in unique_map}
            for family in families:
                hits[family] += 1
    return hits


def expected(quartet_2f: list[str], allele: str, chi: float) -> float:
    return quartet_2f[:2].count(allele) * chi / 2.0 + quartet_2f[2:].count(allele) * (1.0 - chi) / 2.0


def support_residual(quartet_2f: list[str], fractions: dict[str, float], chi: float) -> float:
    keys = set(quartet_2f) | set(fractions)
    return sum(abs(fractions.get(allele, 0.0) - expected(quartet_2f, allele, chi)) for allele in keys)


def with_quartet(row: dict[str, str], new_quartet: list[str], tag: str) -> dict[str, str]:
    out = dict(row)
    for slot, allele in zip(SLOTS, new_quartet):
        out[f"{slot}_full"] = allele
        out[f"{slot}_2field"] = allele_2field(allele)
        out[f"{slot}_report"] = allele_2field(allele)
    out["source"] = tag
    return out


def run_existing_root_test(current_rows, truth_by_set, gmap, out_dir: Path) -> None:
    rows = []
    current_score = score_all(current_rows, truth_by_set, gmap)
    for root in ASM_ROOTS:
        if not root.exists():
            continue
        trial = {sample: dict(rows_by_gene) for sample, rows_by_gene in current_rows.items()}
        changed = 0
        for sample in list(trial):
            alt = load_final(final_path(root, sample))
            if not alt:
                continue
            for gene, row in alt.items():
                if gene in trial[sample]:
                    if quartet(row) != quartet(trial[sample][gene]):
                        changed += 1
                    trial[sample][gene] = row
        score = score_all(trial, truth_by_set, gmap)
        rows.append({
            "root": root.name,
            "changed_gene_rows": changed,
            "two_field": f"{score[0]}/{score[1]}",
            "two_field_delta": score[0] - current_score[0],
            "g_group": f"{score[2]}/{score[3]}",
        })
    write_tsv(out_dir / "existing_read_rescue_roots.tsv", ["root", "changed_gene_rows", "two_field", "two_field_delta", "g_group"], rows)


def run_dpb1_kmer_test(current_rows, truth_by_set, gmap, out_dir: Path) -> None:
    unique_map, unique_counts = build_family_unique_kmers(SCRIPT_DIR / "resources" / "spechla" / "db" / "HLA" / "exon" / "HLA_DPB1.fasta", 31)
    trial = {sample: dict(rows_by_gene) for sample, rows_by_gene in current_rows.items()}
    manifest = []
    for sample, rows_by_gene in current_rows.items():
        row = rows_by_gene.get("HLA-DPB1")
        if row is None:
            continue
        hits = fastq_support(sample, "HLA-DPB1", unique_map, 31)
        filtered = {allele: count for allele, count in hits.items() if count >= 3 and unique_counts[allele] >= 100 and allele_number(allele) <= 100}
        if not filtered:
            continue
        total = sum(filtered.values())
        fractions = {allele: count / total for allele, count in filtered.items()}
        chi = read_chi_r(sample)
        if chi is None:
            continue
        current_q = [allele_2field(value) for value in quartet(row)]
        candidates = sorted(set(current_q) | set(sorted(filtered, key=lambda a: -filtered[a])[:8]))
        neighbor_quartets = {tuple(current_q)}
        for idx in range(4):
            for allele in candidates:
                trial_q = list(current_q)
                trial_q[idx] = allele
                neighbor_quartets.add(tuple(trial_q))
        best = min(neighbor_quartets, key=lambda q: support_residual(list(q), fractions, chi))
        current_res = support_residual(current_q, fractions, chi)
        best_res = support_residual(list(best), fractions, chi)
        if list(best) != current_q and current_res - best_res >= 0.05:
            trial[sample]["HLA-DPB1"] = with_quartet(row, list(best), "offline-dpb1-kmer")
            manifest.append({
                "sample": sample,
                "old": ",".join(current_q),
                "new": ",".join(best),
                "current_residual": f"{current_res:.6f}",
                "best_residual": f"{best_res:.6f}",
                "top_support": ";".join(f"{a}:{filtered[a]}" for a in sorted(filtered, key=lambda a: -filtered[a])[:8]),
            })
    score = score_all(trial, truth_by_set, gmap)
    current_score = score_all(current_rows, truth_by_set, gmap)
    write_tsv(out_dir / "dpb1_unique_kmer_manifest.tsv", ["sample", "old", "new", "current_residual", "best_residual", "top_support"], manifest)
    write_tsv(out_dir / "dpb1_unique_kmer_summary.tsv", ["strategy", "changed_gene_rows", "two_field", "two_field_delta", "g_group"], [{
        "strategy": "dpb1_unique_kmer_neighbor",
        "changed_gene_rows": len(manifest),
        "two_field": f"{score[0]}/{score[1]}",
        "two_field_delta": score[0] - current_score[0],
        "g_group": f"{score[2]}/{score[3]}",
    }])


def allele_number(allele: str) -> int:
    match = re.search(r"\*(\d+):", allele or "")
    return int(match.group(1)) if match else 999999


def drb1_locus(allele: str) -> str:
    family = allele_2field(allele).split("*", 1)[1].split(":", 1)[0] if "*" in allele_2field(allele) else ""
    if family in {"03", "11", "12", "13", "14"}:
        return "DRB3"
    if family in {"04", "07", "09"}:
        return "DRB4"
    if family in {"15", "16"}:
        return "DRB5"
    return "absent"


def allele_locus(allele: str) -> str:
    clean = allele.replace("HLA-", "")
    return clean.split("*", 1)[0] if "*" in clean else "absent"


def read_tf_counts(sample: str, gene: str) -> Counter[str]:
    counts = Counter()
    for root in SPECHLA_ROOTS:
        path = root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
        if not path.exists():
            continue
        for row in read_tsv(path):
            allele = row.get("allele_2field") or allele_2field(row.get("allele", ""))
            try:
                counts[allele] = max(counts[allele], float(row.get("em_weight") or row.get("fraction") or 0))
            except ValueError:
                pass
        if counts:
            break
    return counts


def run_drb345_consistency_test(current_rows, truth_by_set, gmap, out_dir: Path) -> None:
    trial = {sample: dict(rows_by_gene) for sample, rows_by_gene in current_rows.items()}
    manifest = []
    for sample, rows_by_gene in current_rows.items():
        drb1 = rows_by_gene.get("HLA-DRB1")
        drb345 = rows_by_gene.get("HLA-DRB345")
        if drb1 is None or drb345 is None:
            continue
        counts = read_tf_counts(sample, "HLA-DRB1")
        if not counts:
            continue
        by_locus = defaultdict(list)
        for allele, count in counts.items():
            by_locus[drb1_locus(allele)].append((allele, count))
        for values in by_locus.values():
            values.sort(key=lambda item: (-item[1], item[0]))
        old = [allele_2field(value) for value in quartet(drb1)]
        desired = [allele_locus(value) for value in quartet(drb345)]
        new = list(old)
        changed = False
        for idx, locus in enumerate(desired):
            if locus not in {"DRB3", "DRB4", "DRB5"}:
                continue
            if drb1_locus(new[idx]) == locus:
                continue
            replacement = next((allele for allele, _count in by_locus.get(locus, []) if allele not in new), None)
            if replacement:
                new[idx] = replacement
                changed = True
        if changed:
            trial[sample]["HLA-DRB1"] = with_quartet(drb1, new, "offline-drb345-consistency")
            manifest.append({"sample": sample, "old": ",".join(old), "drb345_loci": ",".join(desired), "new": ",".join(new)})
    score = score_all(trial, truth_by_set, gmap)
    current_score = score_all(current_rows, truth_by_set, gmap)
    write_tsv(out_dir / "drb345_consistency_manifest.tsv", ["sample", "old", "drb345_loci", "new"], manifest)
    write_tsv(out_dir / "drb345_consistency_summary.tsv", ["strategy", "changed_gene_rows", "two_field", "two_field_delta", "g_group"], [{
        "strategy": "drb345_locus_consistency",
        "changed_gene_rows": len(manifest),
        "two_field": f"{score[0]}/{score[1]}",
        "two_field_delta": score[0] - current_score[0],
        "g_group": f"{score[2]}/{score[3]}",
    }])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=ROOT / "diagnostics" / "three_rescue_ideas_20260518")
    parser.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP)
    args = parser.parse_args()

    gmap = ev.load_g_group(args.g_group)
    truth_by_set = {name: ev.load_truth(path) for name, path in TRUTH.items()}
    current = {sample: load_final(final_path(root, sample)) for sample, root, _label in current_roots()}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    current_score = score_all(current, truth_by_set, gmap)
    write_tsv(args.out_dir / "current_summary.tsv", ["strategy", "two_field", "g_group"], [{
        "strategy": "current",
        "two_field": f"{current_score[0]}/{current_score[1]}",
        "g_group": f"{current_score[2]}/{current_score[3]}",
    }])
    run_existing_root_test(current, truth_by_set, gmap, args.out_dir)
    run_dpb1_kmer_test(current, truth_by_set, gmap, args.out_dir)
    run_drb345_consistency_test(current, truth_by_set, gmap, args.out_dir)
    for path in sorted(args.out_dir.glob("*summary.tsv")):
        print(path)
        print(path.read_text().strip())


if __name__ == "__main__":
    main()