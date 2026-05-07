#!/usr/bin/env python3
"""Exon-level G group typing from assembled hap FASTAs.

This is a fallback/diagnostic path for genes with high masked full-length
sequence. It follows the SpecHLA idea: compare inferred hap sequences against
per-gene exon databases, then collapse the best exon hit to G group using
hla_nom_g.txt.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SPECHLA = Path(os.environ.get("SPECHLA", SCRIPT_DIR.parent / "SpecHLA"))
DEFAULT_EXON_DB = DEFAULT_SPECHLA / "db" / "HLA" / "exon"
DEFAULT_G_GROUP = DEFAULT_SPECHLA / "db" / "HLA" / "hla_nom_g.txt"
DEFAULT_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1"]


def load_g_group(path: Path):
    gmap = {}
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


def clean_allele(allele: str) -> str:
    allele = allele.replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    fields = rest.split(":")
    if fields[-1] and fields[-1][-1].isalpha() and fields[-1][-1] != "G":
        fields[-1] = fields[-1][:-1]
    return f"{gene}*{':'.join(fields)}"


def to_2field(allele: str) -> str:
    allele = clean_allele(allele)
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    fields = rest.replace("G", "").split(":")
    return f"{gene}*{':'.join(fields[:2])}"


def to_g_group(allele: str, gmap) -> str:
    allele = clean_allele(allele)
    if allele.endswith("G"):
        return allele
    if "*" not in allele:
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


def parse_exon_db_headers(path: Path):
    id_to_allele = {}
    with path.open() as f:
        for line in f:
            if not line.startswith(">"):
                continue
            text = line[1:].strip()
            seq_id = text.split()[0]
            match = re.search(r"([A-Z0-9]+\*[0-9:]+[A-Z]?)", text)
            if match:
                id_to_allele[seq_id] = match.group(1)
    return id_to_allele


def n_fraction(path: Path) -> float:
    seq = []
    with path.open() as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    s = "".join(seq)
    return 1.0 if not s else s.count("N") / len(s)


def best_exon_hit(hap_fa: Path, exon_db_prefix: Path, id_to_allele):
    outfmt = "6 qseqid sseqid pident length bitscore"
    with tempfile.NamedTemporaryFile("w+t", suffix=".blast") as tmp:
        cmd = [
            "blastn", "-query", str(hap_fa), "-db", str(exon_db_prefix),
            "-outfmt", outfmt, "-max_target_seqs", "20000", "-dust", "no",
            "-out", tmp.name,
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        scores = {}
        with open(tmp.name) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                _, sid, pident, length, bitscore = parts[:5]
                allele = id_to_allele.get(sid)
                if not allele:
                    continue
                score = float(bitscore)
                ident_len = float(pident) * int(length)
                prev = scores.get(allele)
                if prev is None:
                    scores[allele] = [score, ident_len, int(length)]
                else:
                    prev[0] += score
                    prev[1] += ident_len
                    prev[2] += int(length)
        if not scores:
            return "no_match", 0.0, 0
        allele, vals = max(scores.items(), key=lambda kv: (kv[1][0], kv[1][1], kv[1][2]))
        return allele, vals[0], vals[2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-root", required=True, type=Path)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--exon-db", type=Path, default=DEFAULT_EXON_DB)
    ap.add_argument("--g-group", type=Path, default=DEFAULT_G_GROUP)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    gmap = load_g_group(args.g_group)
    out_path = args.out or (args.asm_root / args.sample / f"{args.sample}.exon_calls.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        cols = ["sample", "gene", "global_hap", "assignment", "exon_best", "exon_g_group", "exon_2field", "n_fraction", "blast_score", "aligned_len"]
        out.write("\t".join(cols) + "\n")
        for gene in args.genes:
            d = args.asm_root / args.sample / gene.lower() / gene
            exon_prefix = args.exon_db / f"{gene.replace('-', '_')}.fasta"
            if not exon_prefix.exists():
                continue
            id_to_allele = parse_exon_db_headers(exon_prefix)
            for hap in range(1, 5):
                hap_fa = d / f"hap{hap}.fa"
                if not hap_fa.exists():
                    continue
                allele, score, aln_len = best_exon_hit(hap_fa, exon_prefix, id_to_allele)
                assignment = "R" if hap <= 2 else "D"
                row = {
                    "sample": args.sample,
                    "gene": gene,
                    "global_hap": str(hap),
                    "assignment": assignment,
                    "exon_best": allele,
                    "exon_g_group": to_g_group(allele, gmap) if allele != "no_match" else "no_match",
                    "exon_2field": to_2field(allele) if allele != "no_match" else "no_match",
                    "n_fraction": f"{n_fraction(hap_fa):.4f}",
                    "blast_score": f"{score:.2f}",
                    "aligned_len": str(aln_len),
                }
                out.write("\t".join(row[c] for c in cols) + "\n")
    print(f"[exon-typing] wrote {out_path}")


if __name__ == "__main__":
    main()
