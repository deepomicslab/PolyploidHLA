#!/usr/bin/env python3
"""Aggregate per-gene calls.tsv files into a single summary table.

Reads <ASM_ROOT>/<SAMPLE>/<gene_lc>/<TAG>/calls.tsv for each gene and writes
a single tab-separated file. The full allele calls are preserved, and extra
2-field, G group, and report columns are emitted for low-resolution / G group
truth comparison and low-coverage genes.

    sample  gene  R1_full  R2_full  D1_full  D2_full  R1_fraction ...

`R1_fraction`/`R2_fraction`/`D1_fraction`/`D2_fraction` are the modelled
haplotype proportions for each reported call. `source` is `em-refined` if a
`calls.baseline.tsv` sibling exists (meaning the EM stage overrode the
baseline), otherwise `baseline`.

Usage:
    aggregate_calls.py --asm-root asm_v2 --sample mySample \\
        [--genes HLA-A HLA-B ...] [--out final_calls.tsv]

Defaults to the 6 typed genes and writes
<asm-root>/<sample>/<sample>.final_calls.tsv.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]
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


def collect(asm_root: Path, sample: str, genes, mask_warn: float, gmap):
    out_rows = []
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
                "source": "missing", "mean_mask_fraction": "NA",
                "report_level": "missing", "warning": "missing_calls_tsv",
            })
            continue
        rows = read_calls(calls)
        r_rows = [row for row in rows if row.get("assignment") == "R"]
        d_rows = [row for row in rows if row.get("assignment") == "D"]
        r_rows = (r_rows + [None, None])[:2]
        d_rows = (d_rows + [None, None])[:2]
        rs = [call_value(row, "allele") for row in r_rows]
        ds = [call_value(row, "allele") for row in d_rows]
        rf = [call_fraction(row) for row in r_rows]
        df = [call_fraction(row) for row in d_rows]
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
            "source": source,
            "mean_mask_fraction": "NA" if mean_mask is None else f"{mean_mask:.4f}",
            "report_level": report_level,
            "warning": warning,
        })
    return out_rows


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
    args = ap.parse_args()

    gmap = load_g_group(args.g_group)
    rows = collect(args.asm_root, args.sample, args.genes, args.mask_warn, gmap)
    out_path = args.out or (args.asm_root / args.sample / f"{args.sample}.final_calls.tsv")
    cols = [
        "sample", "gene",
        "R1_full", "R2_full", "D1_full", "D2_full",
        "R1_2field", "R2_2field", "D1_2field", "D2_2field",
        "R1_g_group", "R2_g_group", "D1_g_group", "D2_g_group",
        "R1_report", "R2_report", "D1_report", "D2_report",
        "R1_fraction", "R2_fraction", "D1_fraction", "D2_fraction",
        "source", "mean_mask_fraction", "report_level", "warning",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    sys.stderr.write(f"[aggregate] wrote {out_path} ({len(rows)} genes)\n")


if __name__ == "__main__":
    main()
