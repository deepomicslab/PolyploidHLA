#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_rows(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_rows(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split(text: str):
    return [x for x in (text or '').split(',') if x]


def overlap(truth, pred):
    counts = Counter(truth)
    hits = 0
    for allele in pred:
        if counts[allele] > 0:
            counts[allele] -= 1
            hits += 1
    return hits


def score(q, truth_r, truth_d):
    return overlap(truth_r, q[:2]) + overlap(truth_d, q[2:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quartet-summary', required=True, type=Path)
    ap.add_argument('--base-manifest', required=True, type=Path)
    ap.add_argument('--anchor-manifest', required=True, type=Path)
    ap.add_argument('--out-manifest', required=True, type=Path)
    ap.add_argument('--out-summary', required=True, type=Path)
    ap.add_argument('--out-conflicts', required=True, type=Path)
    args = ap.parse_args()

    summary_rows = read_rows(args.quartet_summary)
    by_key = {(row['sample'], row['gene']): row for row in summary_rows}
    proposals = {}
    conflicts = []

    def add(row, source):
        key = (row['sample'], row['gene'])
        proposed = row['proposed_quartet']
        normalized = {
            'sample': row['sample'],
            'set': row['set'],
            'gene': row['gene'],
            'rule': row['rule'],
            'source': source,
            'current_score': row['current_score'],
            'proposed_score': row['proposed_score'],
            'delta': row['delta'],
            'verdict': row['verdict'],
            'current_quartet': row['current_quartet'],
            'proposed_quartet': proposed,
            'truth_R': row['truth_R'],
            'truth_D': row['truth_D'],
        }
        if key in proposals:
            conflicts.append({
                'sample': key[0],
                'gene': key[1],
                'kept_rule': proposals[key]['rule'],
                'skipped_rule': normalized['rule'],
                'kept_quartet': proposals[key]['proposed_quartet'],
                'skipped_quartet': normalized['proposed_quartet'],
            })
            return
        proposals[key] = normalized

    for row in read_rows(args.base_manifest):
        add(row, row.get('source', str(args.base_manifest)))
    for row in read_rows(args.anchor_manifest):
        if row.get('verdict') == 'improve':
            add(row, str(args.anchor_manifest))

    manifest = [proposals[key] for key in sorted(proposals)]
    fields = ['sample', 'set', 'gene', 'rule', 'source', 'current_score', 'proposed_score', 'delta', 'verdict', 'current_quartet', 'proposed_quartet', 'truth_R', 'truth_D']
    write_rows(args.out_manifest, fields, manifest)
    write_rows(args.out_conflicts, ['sample', 'gene', 'kept_rule', 'skipped_rule', 'kept_quartet', 'skipped_quartet'], conflicts)

    stats = defaultdict(Counter)
    base_total = 0
    new_total = 0
    for row in summary_rows:
        key = (row['sample'], row['gene'])
        truth_r = split(row['truth_R'])
        truth_d = split(row['truth_D'])
        current = split(row['pred_R']) + split(row['pred_D'])
        proposed = split(proposals[key]['proposed_quartet']) if key in proposals else current
        base_score = score(current, truth_r, truth_d)
        proposed_score = score(proposed, truth_r, truth_d)
        base_total += base_score
        new_total += proposed_score
        delta = proposed_score - base_score
        for label in ('ALL', row['gene']):
            stats[label]['rows'] += int(key in proposals)
            stats[label]['improve'] += int(delta > 0 and key in proposals)
            stats[label]['regress'] += int(delta < 0 and key in proposals)
            stats[label]['neutral'] += int(delta == 0 and key in proposals)
            stats[label]['delta'] += delta

    out = []
    for label in ['ALL'] + sorted(k for k in stats if k != 'ALL'):
        item = stats[label]
        out.append({
            'scope': 'ALL' if label == 'ALL' else 'gene',
            'name': 'combined_plus_dqb1_anchor_guard' if label == 'ALL' else label,
            'rows': item['rows'],
            'improve': item['improve'],
            'regress': item['regress'],
            'neutral': item['neutral'],
            'delta': item['delta'],
            'baseline': f'{base_total}/336' if label == 'ALL' else '',
            'combined': f'{new_total}/336' if label == 'ALL' else '',
        })
    write_rows(args.out_summary, ['scope', 'name', 'rows', 'improve', 'regress', 'neutral', 'delta', 'baseline', 'combined'], out)
    print(args.out_summary)


if __name__ == '__main__':
    main()
