#!/usr/bin/env python3
"""Aggregate per-gene calls.tsv files into a single summary table.

Reads <ASM_ROOT>/<SAMPLE>/<gene_lc>/<TAG>/calls.tsv for each gene and writes
a single tab-separated file. The full allele calls are preserved, and extra
2-field, G group, and report columns are emitted for low-resolution / G group
truth comparison and low-coverage genes.

    sample  gene  R1_full  R2_full  D1_full  D2_full  R1_fraction  R1_read_fraction ...

`R1_fraction`/`R2_fraction`/`D1_fraction`/`D2_fraction` are the modelled
haplotype proportions for each reported call. `R1_read_fraction` etc. are the
allele-family read support fractions from EM `tf_counts.tsv` when available.
`R1_copy_fraction_fit` etc. are fitted copy-level fractions constrained to sum
to 1 when read support is available. A concise companion file keeps only the
reported alleles and copy fractions.
`source` is `em-refined` if a `calls.baseline.tsv` sibling exists (meaning the
EM stage overrode the baseline), otherwise `baseline`.

Usage:
    aggregate_calls.py --asm-root asm_v2 --sample mySample \\
        [--genes HLA-A HLA-B ...] [--out final_calls.tsv]

Defaults to the 6 typed genes and writes
<asm-root>/<sample>/<sample>.final_calls.tsv plus
<asm-root>/<sample>/<sample>.final_calls.compact.tsv.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np

DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]
SLOTS = ("R1", "R2", "D1", "D2")
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BUNDLED_SPECHLA = SCRIPT_DIR / "resources" / "spechla"
DEFAULT_LEGACY_SPECHLA = SCRIPT_DIR.parent / "SpecHLA"
DEFAULT_SPECHLA = Path(os.environ.get("SPECHLA", DEFAULT_BUNDLED_SPECHLA if DEFAULT_BUNDLED_SPECHLA.exists() else DEFAULT_LEGACY_SPECHLA))
DEFAULT_G_GROUP = DEFAULT_SPECHLA / "db" / "HLA" / "hla_nom_g.txt"


def load_g_group(path: Path):
    gmap = {}
    if not path.exists():
        return gmap
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(";")
            if len(parts) < 3:
                continue
            gene = parts[0].rstrip("*")
            members = parts[-2].split("/") if parts[-2] else []
            group = parts[-1] or parts[-2]
            if not group:
                continue
            group_name = f"{gene}*{group}"
            for member in members:
                gmap[f"{gene}*{member}"] = group_name
            gmap[group_name] = group_name
    return gmap


def strip_expr_suffix(field: str) -> str:
    return field[:-1] if field and field[-1].isalpha() and field[-1] != "G" else field


def clean_allele(allele: str) -> str:
    if not allele or allele == "NA" or "*" not in allele:
        return allele or "NA"
    gene, fields = allele.replace("HLA-", "").split("*", 1)
    parts = fields.split(":")
    parts[-1] = strip_expr_suffix(parts[-1])
    return f"{gene}*{':'.join(parts)}"


def allele_2field(allele: str) -> str:
    allele = clean_allele(allele)
    if allele == "NA" or "*" not in allele:
        return allele
    gene, fields = allele.split("*", 1)
    parts = fields.replace("G", "").split(":")
    return f"{gene}*{':'.join(parts[:2])}" if len(parts) >= 2 else f"{gene}*{parts[0]}"


def allele_g_group(allele: str, gmap) -> str:
    allele = clean_allele(allele)
    if allele == "NA" or "*" not in allele:
        return allele
    if allele.endswith("G"):
        return allele
    gene, rest = allele.split("*", 1)
    candidates = [allele]
    fields = rest.split(":")
    if len(fields) == 2:
        candidates.extend([f"{gene}*{rest}:01", f"{gene}*{rest}:01:01"])
    elif len(fields) == 3:
        candidates.append(f"{gene}*{rest}:01")
    for candidate in candidates:
        if candidate in gmap:
            return gmap[candidate]
    return allele


def fasta_n_fraction(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    seq = []
    with path.open() as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    s = "".join(seq)
    return None if not s else s.count("N") / len(s)


def read_calls(path: Path):
    """Return per-haplotype call rows sorted by global_hap."""
    rows = []
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_h = header.index("global_hap")
            i_a = header.index("assignment")
            i_l = header.index("allele")
        except ValueError:
            sys.stderr.write(f"[warn] {path}: unexpected header {header}\n")
            return []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_h, i_a, i_l):
                continue
            row = {name: parts[idx] if idx < len(parts) else "" for idx, name in enumerate(header)}
            rows.append(row)
    rows.sort(key=lambda r: int(r.get("global_hap", "")) if r.get("global_hap", "").isdigit() else r.get("global_hap", ""))
    return rows


def read_tsv_dicts(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header == [""]:
            return []
        rows = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rows.append({name: parts[idx] if idx < len(parts) else "" for idx, name in enumerate(header)})
        return rows


def load_read_support(spechla_root: Optional[Path], sample: str, gene: str) -> tuple[dict[str, dict[str, str]], bool]:
    if spechla_root is None:
        return {}, False
    path = spechla_root / sample / "em_refine" / f"{gene}.tf_counts.tsv"
    support = {}
    for row in read_tsv_dicts(path):
        allele = row.get("allele_2field", "")
        if not allele:
            continue
        support[allele] = {
            "read_count": row.get("em_weight", "NA"),
            "read_fraction": row.get("fraction", "NA"),
        }
    return support, path.exists() and bool(support)


def read_chi_r(spechla_root: Optional[Path], sample: str) -> Optional[float]:
    if spechla_root is None:
        return None
    sample_dir = spechla_root / sample
    pooled = sample_dir / f"{sample}.chi_pooled.txt"
    if pooled.exists():
        for line in pooled.read_text().splitlines():
            if not line.startswith("GLOBAL") or "chi_R=" not in line:
                continue
            for item in line.split():
                if item.startswith("chi_R="):
                    try:
                        value = float(item.split("=", 1)[1])
                    except ValueError:
                        continue
                    if 0.02 <= value <= 0.49:
                        return value
    chimerism = sample_dir / f"{sample}.chimerism.txt"
    fallback = None
    if chimerism.exists():
        for line in chimerism.read_text().splitlines():
            if "chi_R=" not in line:
                continue
            for item in line.split():
                if item.startswith("chi_R="):
                    try:
                        fallback = float(item.split("=", 1)[1])
                    except ValueError:
                        pass
    return fallback if fallback is not None and 0.0 < fallback < 1.0 else None


def call_value(row: Optional[dict], key: str, default: str = "NA") -> str:
    if not row:
        return default
    value = row.get(key, default)
    return value if value not in {"", None} else default


def call_fraction(row: Optional[dict]) -> str:
    if not row:
        return "NA"
    for key in ("hap_fraction", "haplotype_fraction", "fraction", "ratio"):
        value = row.get(key)
        if value not in {None, "", "NA"}:
            try:
                return f"{float(value):.6f}"
            except ValueError:
                return str(value)
    return "NA"


def format_float(value: str, digits: int) -> str:
    if value in {"", "NA", None}:
        return "NA"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def call_read_support(row: Optional[dict], allele: str, support: dict[str, dict[str, str]], complete_support: bool) -> tuple[str, str]:
    if row:
        row_fraction = row.get("allele_read_fraction") or row.get("read_fraction")
        row_count = row.get("allele_read_count") or row.get("read_count")
        if row_fraction or row_count:
            return format_float(row_fraction or "NA", 6), format_float(row_count or "NA", 2)
    allele_support = support.get(allele_2field(allele), {})
    if not allele_support and complete_support and allele_2field(allele) != "NA":
        return "0.000000", "0.00"
    return (
        format_float(allele_support.get("read_fraction", "NA"), 6),
        format_float(allele_support.get("read_count", "NA"), 2),
    )


def float_or_none(value) -> Optional[float]:
    try:
        if value in (None, "", "NA"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def chi_from_slot_fractions(slot_fractions: list[str], fallback: Optional[float]) -> Optional[float]:
    values = [float_or_none(value) for value in slot_fractions]
    if values[0] is not None and values[1] is not None:
        chi_r = values[0] + values[1]
        if 0.0 < chi_r < 1.0:
            return chi_r
    return fallback


def fit_copy_fractions(slot_alleles: list[str], slot_read_fractions: list[str], chi_r: Optional[float]) -> dict[str, str]:
    alleles_by_slot = [allele_2field(allele) for allele in slot_alleles]
    support_by_allele: dict[str, float] = {}
    for allele, fraction in zip(alleles_by_slot, slot_read_fractions):
        value = float_or_none(fraction)
        if allele == "NA" or value is None:
            continue
        support_by_allele[allele] = max(support_by_allele.get(allele, 0.0), value)
    support_sum = sum(support_by_allele.values())
    if support_sum <= 0:
        return {
            "R1_copy_fraction_fit": "NA", "R2_copy_fraction_fit": "NA",
            "D1_copy_fraction_fit": "NA", "D2_copy_fraction_fit": "NA",
            "copy_fit_error": "NA", "copy_identifiability": "no_read_support",
            "copy_chi_r": "NA", "allele_support_fraction_sum": "0.000000",
        }

    alleles = sorted(support_by_allele)
    y = np.array([support_by_allele[allele] / support_sum for allele in alleles], dtype=float)
    matrix = np.zeros((len(alleles), len(SLOTS)), dtype=float)
    for row_idx, allele in enumerate(alleles):
        for col_idx, slot_allele in enumerate(alleles_by_slot):
            if slot_allele == allele:
                matrix[row_idx, col_idx] = 1.0

    design = matrix.copy()
    target = y.copy()
    chi_used = chi_r if chi_r is not None and 0.0 < chi_r < 1.0 else None
    if chi_used is not None:
        weight = 0.25
        design = np.vstack([
            design,
            np.sqrt(weight) * np.array([1.0, 1.0, 0.0, 0.0]),
            np.sqrt(weight) * np.array([0.0, 0.0, 1.0, 1.0]),
        ])
        target = np.concatenate([target, np.sqrt(weight) * np.array([chi_used, 1.0 - chi_used])])

    best_x = None
    best_score = None
    for mask in range(1, 1 << len(SLOTS)):
        active = [idx for idx in range(len(SLOTS)) if mask & (1 << idx)]
        sub_design = design[:, active]
        gram = 2.0 * sub_design.T @ sub_design + 1e-9 * np.eye(len(active))
        rhs = 2.0 * sub_design.T @ target
        equality = np.ones((1, len(active)))
        kkt = np.block([[gram, equality.T], [equality, np.zeros((1, 1))]])
        kkt_rhs = np.concatenate([rhs, np.array([1.0])])
        try:
            solution = np.linalg.lstsq(kkt, kkt_rhs, rcond=None)[0][:len(active)]
        except np.linalg.LinAlgError:
            continue
        if np.any(solution < -1e-7):
            continue
        x = np.zeros(len(SLOTS), dtype=float)
        for idx, value in zip(active, solution):
            x[idx] = max(0.0, value)
        total = x.sum()
        if total <= 0:
            continue
        x = x / total
        score = float(np.sum((matrix @ x - y) ** 2))
        if best_score is None or score < best_score:
            best_score = score
            best_x = x

    if best_x is None:
        best_x = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)
        best_score = float(np.sum((matrix @ best_x - y) ** 2))

    identifiable = "identifiable" if np.linalg.matrix_rank(matrix) == len(SLOTS) else "underdetermined"
    if identifiable == "underdetermined" and chi_used is not None:
        identifiable = "underdetermined_chi_regularized"
    elif identifiable == "underdetermined":
        identifiable = "underdetermined_min_norm"
    if best_score is not None and best_score > 0.0025:
        identifiable += ";high_fit_error"
    if np.any(best_x <= 1e-6):
        identifiable += ";boundary_zero"
    return {
        "R1_copy_fraction_fit": f"{best_x[0]:.6f}",
        "R2_copy_fraction_fit": f"{best_x[1]:.6f}",
        "D1_copy_fraction_fit": f"{best_x[2]:.6f}",
        "D2_copy_fraction_fit": f"{best_x[3]:.6f}",
        "copy_fit_error": f"{best_score:.8f}",
        "copy_identifiability": identifiable,
        "copy_chi_r": f"{chi_used:.6f}" if chi_used is not None else "NA",
        "allele_support_fraction_sum": f"{support_sum:.6f}",
    }


def collect(asm_root: Path, sample: str, genes, mask_warn: float, gmap, spechla_root: Optional[Path] = None):
    out_rows = []
    sample_chi_r = read_chi_r(spechla_root, sample)
    for gene in genes:
        gene_lc = gene.lower()
        d = asm_root / sample / gene_lc / gene
        calls = d / "calls.tsv"
        mask_values = [fasta_n_fraction(d / f"hap{i}.fa") for i in range(1, 5)]
        mask_values = [x for x in mask_values if x is not None]
        mean_mask = sum(mask_values) / len(mask_values) if mask_values else None
        if not calls.exists():
            out_rows.append({
                "sample": sample, "gene": gene,
                "R1_full": "NA", "R2_full": "NA", "D1_full": "NA", "D2_full": "NA",
                "R1_fraction": "NA", "R2_fraction": "NA", "D1_fraction": "NA", "D2_fraction": "NA",
                "R1_read_fraction": "NA", "R2_read_fraction": "NA",
                "D1_read_fraction": "NA", "D2_read_fraction": "NA",
                "R1_read_count": "NA", "R2_read_count": "NA",
                "D1_read_count": "NA", "D2_read_count": "NA",
                "R1_copy_fraction_fit": "NA", "R2_copy_fraction_fit": "NA",
                "D1_copy_fraction_fit": "NA", "D2_copy_fraction_fit": "NA",
                "copy_fit_error": "NA", "copy_identifiability": "missing_calls_tsv",
                "copy_chi_r": "NA", "allele_support_fraction_sum": "0.000000",
                "source": "missing", "mean_mask_fraction": "NA",
                "report_level": "missing", "warning": "missing_calls_tsv",
            })
            continue
        read_support, complete_support = load_read_support(spechla_root, sample, gene)
        rows = read_calls(calls)
        r_rows = [row for row in rows if row.get("assignment") == "R"]
        d_rows = [row for row in rows if row.get("assignment") == "D"]
        r_rows = (r_rows + [None, None])[:2]
        d_rows = (d_rows + [None, None])[:2]
        rs = [call_value(row, "allele") for row in r_rows]
        ds = [call_value(row, "allele") for row in d_rows]
        rf = [call_fraction(row) for row in r_rows]
        df = [call_fraction(row) for row in d_rows]
        r_support = [call_read_support(row, allele, read_support, complete_support) for row, allele in zip(r_rows, rs)]
        d_support = [call_read_support(row, allele, read_support, complete_support) for row, allele in zip(d_rows, ds)]
        slot_fractions = [rf[0], rf[1], df[0], df[1]]
        slot_read_fractions = [r_support[0][0], r_support[1][0], d_support[0][0], d_support[1][0]]
        copy_fit = fit_copy_fractions(rs + ds, slot_read_fractions, chi_from_slot_fractions(slot_fractions, sample_chi_r))
        source = "em-refined" if (d / "calls.baseline.tsv").exists() else "baseline"
        high_mask = mean_mask is not None and mean_mask >= mask_warn
        report_level = "2-field" if high_mask else "full"
        warning = "high_mask_report_2field" if high_mask else ""
        out_rows.append({
            "sample": sample, "gene": gene,
            "R1_full": rs[0], "R2_full": rs[1], "D1_full": ds[0], "D2_full": ds[1],
            "R1_2field": allele_2field(rs[0]), "R2_2field": allele_2field(rs[1]),
            "D1_2field": allele_2field(ds[0]), "D2_2field": allele_2field(ds[1]),
            "R1_g_group": allele_g_group(rs[0], gmap), "R2_g_group": allele_g_group(rs[1], gmap),
            "D1_g_group": allele_g_group(ds[0], gmap), "D2_g_group": allele_g_group(ds[1], gmap),
            "R1_report": allele_2field(rs[0]) if high_mask else rs[0],
            "R2_report": allele_2field(rs[1]) if high_mask else rs[1],
            "D1_report": allele_2field(ds[0]) if high_mask else ds[0],
            "D2_report": allele_2field(ds[1]) if high_mask else ds[1],
            "R1_fraction": rf[0], "R2_fraction": rf[1],
            "D1_fraction": df[0], "D2_fraction": df[1],
            "R1_read_fraction": r_support[0][0], "R2_read_fraction": r_support[1][0],
            "D1_read_fraction": d_support[0][0], "D2_read_fraction": d_support[1][0],
            "R1_read_count": r_support[0][1], "R2_read_count": r_support[1][1],
            "D1_read_count": d_support[0][1], "D2_read_count": d_support[1][1],
            "R1_copy_fraction_fit": copy_fit["R1_copy_fraction_fit"],
            "R2_copy_fraction_fit": copy_fit["R2_copy_fraction_fit"],
            "D1_copy_fraction_fit": copy_fit["D1_copy_fraction_fit"],
            "D2_copy_fraction_fit": copy_fit["D2_copy_fraction_fit"],
            "copy_fit_error": copy_fit["copy_fit_error"],
            "copy_identifiability": copy_fit["copy_identifiability"],
            "copy_chi_r": copy_fit["copy_chi_r"],
            "allele_support_fraction_sum": copy_fit["allele_support_fraction_sum"],
            "source": source,
            "mean_mask_fraction": "NA" if mean_mask is None else f"{mean_mask:.4f}",
            "report_level": report_level,
            "warning": warning,
        })
    return out_rows


def default_compact_path(out_path: Path) -> Path:
    if out_path.name.endswith(".final_calls.tsv"):
        return out_path.with_name(out_path.name.replace(".final_calls.tsv", ".final_calls.compact.tsv"))
    return out_path.with_name(f"{out_path.stem}.compact.tsv")


def write_rows(path: Path, cols: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(cols) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in cols) + "\n")


def compact_row(row: dict[str, str]) -> dict[str, str]:
    out = {"sample": row.get("sample", ""), "gene": row.get("gene", "")}
    for slot in SLOTS:
        out[f"{slot}_allele"] = row.get(f"{slot}_report") or row.get(f"{slot}_2field") or row.get(f"{slot}_full", "NA")
        out[f"{slot}_copy_fraction"] = row.get(f"{slot}_copy_fraction_fit", "NA")
    out["copy_identifiability"] = row.get("copy_identifiability", "")
    out["copy_fit_error"] = row.get("copy_fit_error", "")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--mask-warn", type=float, default=0.15,
                    help="mean hap FASTA N fraction above which report columns "
                    "are downgraded to 2-field")
    ap.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP,
                    help="WMDA hla_nom_g.txt used for G group conversion")
    ap.add_argument("--spechla-root", type=Path, default=None,
                    help="per-sample SpecHLA output root used to read EM tf_counts.tsv for allele read support")
    ap.add_argument("--compact-out", type=Path, default=None,
                    help="write concise allele/copy-fraction result table (default: <sample>.final_calls.compact.tsv)")
    ap.add_argument("--no-compact", action="store_true",
                    help="do not write the concise compact result table")
    args = ap.parse_args()

    gmap = load_g_group(args.g_group)
    rows = collect(args.asm_root, args.sample, args.genes, args.mask_warn, gmap, args.spechla_root)
    out_path = args.out or (args.asm_root / args.sample / f"{args.sample}.final_calls.tsv")
    cols = [
        "sample", "gene",
        "R1_full", "R2_full", "D1_full", "D2_full",
        "R1_2field", "R2_2field", "D1_2field", "D2_2field",
        "R1_g_group", "R2_g_group", "D1_g_group", "D2_g_group",
        "R1_report", "R2_report", "D1_report", "D2_report",
        "R1_fraction", "R2_fraction", "D1_fraction", "D2_fraction",
        "R1_read_fraction", "R2_read_fraction", "D1_read_fraction", "D2_read_fraction",
        "R1_read_count", "R2_read_count", "D1_read_count", "D2_read_count",
        "R1_copy_fraction_fit", "R2_copy_fraction_fit", "D1_copy_fraction_fit", "D2_copy_fraction_fit",
        "copy_fit_error", "copy_identifiability", "copy_chi_r", "allele_support_fraction_sum",
        "source", "mean_mask_fraction", "report_level", "warning",
    ]
    write_rows(out_path, cols, rows)
    sys.stderr.write(f"[aggregate] wrote {out_path} ({len(rows)} genes)\n")
    if not args.no_compact:
        compact_path = args.compact_out or default_compact_path(out_path)
        compact_cols = [
            "sample", "gene",
            "R1_allele", "R1_copy_fraction", "R2_allele", "R2_copy_fraction",
            "D1_allele", "D1_copy_fraction", "D2_allele", "D2_copy_fraction",
            "copy_identifiability", "copy_fit_error",
        ]
        write_rows(compact_path, compact_cols, [compact_row(row) for row in rows])
        sys.stderr.write(f"[aggregate] wrote {compact_path} ({len(rows)} genes, compact)\n")


if __name__ == "__main__":
    main()
