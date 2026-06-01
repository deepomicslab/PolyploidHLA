#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SUMMARY = ROOT / 'diagnostics/gendx_quartet_summary_noA4.tsv'
DEFAULT_GENES = ['HLA-A', 'HLA-C', 'HLA-DQB1']


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def split_alleles(text):
    return [x for x in text.split(',') if x]


def score(q, truth_r, truth_d):
    def overlap(truth, pred):
        counts = Counter(truth)
        hits = 0
        for allele in pred:
            if counts[allele] > 0:
                counts[allele] -= 1
                hits += 1
        return hits
    return overlap(truth_r, q[:2]) + overlap(truth_d, q[2:])


def tf_counts(spechla_root: Path, sample: str, gene: str):
    path = spechla_root / sample / 'em_refine' / f'{gene}.tf_counts.tsv'
    rows = read_rows(path)
    out = {}
    order = []
    for row in rows:
        allele = row['allele_2field']
        frac = float(row['fraction'])
        out[allele] = frac
        order.append(allele)
    return out, order


def expected(q, recipient_weight: float, donor_weight: float):
    values = Counter()
    for allele in q[:2]:
        values[allele] += recipient_weight
    for allele in q[2:]:
        values[allele] += donor_weight
    return values


def fit_error(q, obs, recipient_weight: float, donor_weight: float):
    exp = expected(q, recipient_weight, donor_weight)
    alleles = set(obs) | set(exp)
    return sum(abs(obs.get(a, 0.0) - exp.get(a, 0.0)) for a in alleles)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--spechla-root', required=True, type=Path,
                        help='Root containing <sample>/em_refine/*.tf_counts.tsv')
    parser.add_argument('--quartet-summary', default=DEFAULT_SUMMARY, type=Path)
    parser.add_argument('--out', type=Path, default=None)
    parser.add_argument('--genes', nargs='+', default=DEFAULT_GENES)
    parser.add_argument('--chi-r', type=float, default=0.21)
    args = parser.parse_args()

    recipient_weight = args.chi_r / 2
    donor_weight = (1 - args.chi_r) / 2
    out_path = args.out or (ROOT / 'diagnostics' / f'{args.sample}.dose_rescue_candidates.tsv')

    summary = {(row['sample'], row['gene']): row for row in read_rows(args.quartet_summary)}
    rows = []
    for gene in args.genes:
        row = summary[(args.sample, gene)]
        obs, order = tf_counts(args.spechla_root, args.sample, gene)
        current = split_alleles(row['pred_R']) + split_alleles(row['pred_D'])
        truth_r = split_alleles(row['truth_R'])
        truth_d = split_alleles(row['truth_D'])
        candidate_pool = order[:5]
        # Keep search small and truth-free: top EM alleles only, ordered R/D pairs.
        seen = set()
        ranked = []
        for r_pair in itertools.combinations_with_replacement(candidate_pool, 2):
            for d_pair in itertools.combinations_with_replacement(candidate_pool, 2):
                quartet = tuple(r_pair + d_pair)
                if quartet in seen:
                    continue
                seen.add(quartet)
                ranked.append((fit_error(quartet, obs, recipient_weight, donor_weight), quartet))
        ranked.sort(key=lambda item: item[0])
        current_err = fit_error(tuple(current), obs, recipient_weight, donor_weight)
        truth_q = tuple(truth_r + truth_d)
        truth_err = fit_error(truth_q, obs, recipient_weight, donor_weight)
        rows.append({
            'gene': gene,
            'rank': 'current',
            'fit_error': f'{current_err:.6f}',
            'posthoc_score': score(current, truth_r, truth_d),
            'quartet': ','.join(current),
            'truth_quartet': ','.join(truth_q),
        })
        rows.append({
            'gene': gene,
            'rank': 'truth',
            'fit_error': f'{truth_err:.6f}',
            'posthoc_score': score(list(truth_q), truth_r, truth_d),
            'quartet': ','.join(truth_q),
            'truth_quartet': ','.join(truth_q),
        })
        for idx, (err, quartet) in enumerate(ranked[:20], 1):
            rows.append({
                'gene': gene,
                'rank': idx,
                'fit_error': f'{err:.6f}',
                'posthoc_score': score(list(quartet), truth_r, truth_d),
                'quartet': ','.join(quartet),
                'truth_quartet': ','.join(truth_q),
            })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w') as handle:
        fields = ['gene', 'rank', 'fit_error', 'posthoc_score', 'quartet', 'truth_quartet']
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out_path)


if __name__ == '__main__':
    main()