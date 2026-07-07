#!/usr/bin/env python3
"""Prepare optional non-classical HLA/MIC typing resources.

The bundled SpecHLA reference contains the validated classical HLA contigs, while
the IMGT-style allele FASTA also contains additional loci such as HLA-E/F/G/H and
MICA/MICB. This helper builds a small augmented reference package for those loci
without modifying the bundled resource directory.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def gene_prefix(gene: str) -> str:
    return gene.replace("HLA-", "") + "*"


def ref_tag(gene: str) -> str:
    return gene.replace("-", "_")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper().replace("-", "")
                tokens = line[1:].split()
                name = next((token for token in tokens if "*" in token), tokens[0])
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper().replace("-", "")
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for offset in range(0, len(seq), 80):
                handle.write(seq[offset:offset + 80] + "\n")


def run_if_missing(outputs: list[Path], command: list[str], force: bool) -> None:
    if not force and all(path.exists() and path.stat().st_size > 0 for path in outputs):
        return
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare augmented resources for optional extra HLA/MIC typing genes")
    parser.add_argument("--spechla-root", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--base-gene-bed", required=True, type=Path)
    parser.add_argument("--genes", nargs="+", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    spechla_db = args.spechla_root / "db"
    imgt = spechla_db / "ref" / "hla_gen.format.filter.extend.DRB.no26789.v2.fasta"
    base_ref = spechla_db / "ref" / "hla.ref.extend.fa"
    if not imgt.exists():
        raise SystemExit(f"missing IMGT FASTA: {imgt}")
    if not base_ref.exists():
        raise SystemExit(f"missing base reference: {base_ref}")

    db = read_fasta(imgt)
    selected: dict[str, tuple[str, str]] = {}
    for gene in args.genes:
        matches = {name: seq for name, seq in db.items() if name.startswith(gene_prefix(gene))}
        if not matches:
            raise SystemExit(f"no IMGT alleles found for {gene} with prefix {gene_prefix(gene)}")
        representative_name, representative_seq = max(matches.items(), key=lambda item: (len(item[1]), item[0]))
        selected[gene] = (representative_name, representative_seq)

    ref_dir = args.work_dir / "db" / "ref"
    hla_dir = args.work_dir / "db" / "HLA"
    augmented_ref = ref_dir / "hla.ref.extend.extra.fa"
    gene_bed = args.work_dir / "gene.extra.spechla.bed"

    if args.force or not augmented_ref.exists():
        ref_dir.mkdir(parents=True, exist_ok=True)
        with augmented_ref.open("w") as out_handle:
            with base_ref.open() as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
            for gene, (_allele, seq) in selected.items():
                out_handle.write(f">{ref_tag(gene)}\n")
                for offset in range(0, len(seq), 80):
                    out_handle.write(seq[offset:offset + 80] + "\n")

    if args.force or not gene_bed.exists():
        args.work_dir.mkdir(parents=True, exist_ok=True)
        with gene_bed.open("w") as out_handle:
            with args.base_gene_bed.open() as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
            for gene, (_allele, seq) in selected.items():
                out_handle.write(f"{ref_tag(gene)}\t0\t{len(seq)}\t{gene}\n")

    for gene, (_allele, seq) in selected.items():
        tag = ref_tag(gene)
        gene_fasta = hla_dir / tag / f"{tag}.fa"
        if args.force or not gene_fasta.exists():
            write_fasta(gene_fasta, {tag: seq})
        run_if_missing([gene_fasta.with_suffix(gene_fasta.suffix + ".fai")], ["samtools", "faidx", str(gene_fasta)], args.force)
        run_if_missing(
            [Path(str(gene_fasta) + suffix) for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]],
            ["bwa", "index", str(gene_fasta)],
            args.force,
        )

    run_if_missing([Path(str(augmented_ref) + ".fai")], ["samtools", "faidx", str(augmented_ref)], args.force)
    run_if_missing(
        [Path(str(augmented_ref) + suffix) for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]],
        ["bwa", "index", str(augmented_ref)],
        args.force,
    )
    dict_out = augmented_ref.with_suffix(".dict")
    if shutil.which("samtools"):
        run_if_missing([dict_out], ["samtools", "dict", "-o", str(dict_out), str(augmented_ref)], args.force)

    print(f"HLA_REF={augmented_ref}")
    print(f"GENE_BED={gene_bed}")
    print(f"EXTRA_HLA_DIR={hla_dir}")


if __name__ == "__main__":
    main()