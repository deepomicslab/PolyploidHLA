#!/usr/bin/env python3
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pysam


SAMPLE = "267018-HLA-20260415EM_S7_L001"
BAM = Path("/data2/wangxuedong/polyploid-hla-realsets/spechla_out_abc_realsets_rescue_20260512") / SAMPLE / f"{SAMPLE}.map_database.bam"
OUT = Path("/data6/wangxuedong/polyploid_hla/diagnostics/267018_raw_db_allele_support.tsv")
TOKENS = [
    "C*04:01",
    "C*06:02",
    "C*06:09",
    "C*07:01",
    "DQB1*02:01",
    "DQB1*02:02",
    "DQB1*02:82",
    "DQB1*03:01",
]


def normalize_ref(ref: str) -> str:
    return ref.replace("HLA-", "").replace("_", "*", 1).replace("_", ":")


def family(ref: str) -> str | None:
    normalized = normalize_ref(ref)
    for token in TOKENS:
        if normalized.startswith(token):
            return token
    return None


def alignment_score(aln) -> int:
    try:
        return int(aln.get_tag("AS"))
    except KeyError:
        return int(aln.query_alignment_length or 0)


def main() -> None:
    summary = {
        token: {
            "alignments": 0,
            "qnames": set(),
            "top_qnames": set(),
            "near_top_qnames": set(),
            "refs": defaultdict(int),
        }
        for token in TOKENS
    }
    per_read = defaultdict(list)
    with pysam.AlignmentFile(str(BAM), "rb") as bam:
        refs = bam.references
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_id < 0:
                continue
            ref = refs[aln.reference_id]
            fam = family(ref)
            if fam is None:
                continue
            score = alignment_score(aln)
            qname = aln.query_name
            summary[fam]["alignments"] += 1
            summary[fam]["qnames"].add(qname)
            summary[fam]["refs"][normalize_ref(ref)] += 1
            per_read[qname].append((fam, score, normalize_ref(ref)))

    for qname, rows in per_read.items():
        best = max(score for _fam, score, _ref in rows)
        best_families = {fam for fam, score, _ref in rows if score == best}
        near_families = {fam for fam, score, _ref in rows if score >= best - 5}
        for fam in best_families:
            summary[fam]["top_qnames"].add(qname)
        for fam in near_families:
            summary[fam]["near_top_qnames"].add(qname)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as handle:
        handle.write("family\talignments\tdistinct_reads\ttop_tie_reads\tnear_top_minus5_reads\ttop_refs\n")
        for fam in TOKENS:
            refs = summary[fam]["refs"]
            top_refs = ";".join(
                f"{ref}:{count}" for ref, count in sorted(refs.items(), key=lambda item: (-item[1], item[0]))[:8]
            )
            handle.write(
                f"{fam}\t{summary[fam]['alignments']}\t{len(summary[fam]['qnames'])}\t"
                f"{len(summary[fam]['top_qnames'])}\t{len(summary[fam]['near_top_qnames'])}\t{top_refs}\n"
            )

        handle.write("\ncomparison\treads_A\treads_B\toverlap\toverlap_frac_of_B\ttop_overlap\tnear_top_overlap\n")
        pairs = [
            ("C*06:02", "C*06:09"),
            ("C*04:01", "C*06:09"),
            ("C*07:01", "C*06:09"),
            ("DQB1*02:01", "DQB1*02:82"),
            ("DQB1*02:02", "DQB1*02:82"),
            ("DQB1*03:01", "DQB1*02:82"),
        ]
        for left, right in pairs:
            left_reads = summary[left]["qnames"]
            right_reads = summary[right]["qnames"]
            overlap = left_reads & right_reads
            top_overlap = summary[left]["top_qnames"] & summary[right]["top_qnames"]
            near_overlap = summary[left]["near_top_qnames"] & summary[right]["near_top_qnames"]
            frac = len(overlap) / len(right_reads) if right_reads else 0.0
            handle.write(
                f"{left} vs {right}\t{len(left_reads)}\t{len(right_reads)}\t{len(overlap)}\t"
                f"{frac:.6f}\t{len(top_overlap)}\t{len(near_overlap)}\n"
            )


if __name__ == "__main__":
    main()