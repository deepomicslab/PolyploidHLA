#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from hla_ld_maps import load_drb1_dqb1_map  # noqa: E402

DEFAULT_DRB1_DQB1_LD_MAP = SCRIPT_ROOT / 'resources' / 'drb1_dqb1_ld.tsv'


def rows(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path: Path, fieldnames, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


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


def load_final_row(asm_root: Path | None, sample: str, gene: str):
    if asm_root is None:
        return {}
    path = asm_root / sample / f'{sample}.final_calls.tsv'
    if not path.exists():
        return {}
    for row in rows(path):
        if row.get('gene') == gene:
            return row
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quartet-summary', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--summary', required=True, type=Path)
    ap.add_argument('--asm-root', type=Path, default=None)
    ap.add_argument('--require-drb1-no-boundary-zero', action='store_true')
    ap.add_argument('--drb1-dqb1-ld-map', type=Path, default=DEFAULT_DRB1_DQB1_LD_MAP)
    args = ap.parse_args()

    drb1_to_dqb1, _dqb1_to_drb1 = load_drb1_dqb1_map(args.drb1_dqb1_ld_map)
    by_key = {(row['sample'], row['gene']): row for row in rows(args.quartet_summary)}
    out_rows = []
    stats = Counter()
    for (sample, gene), dqb1 in sorted(by_key.items()):
        if gene != 'HLA-DQB1':
            continue
        drb1 = by_key.get((sample, 'HLA-DRB1'))
        if not drb1:
            continue
        drb1_final = load_final_row(args.asm_root, sample, 'HLA-DRB1')
        drb1_quality = ';'.join([
            drb1_final.get('copy_identifiability', ''),
            drb1_final.get('warning', ''),
        ])
        if args.require_drb1_no_boundary_zero and 'boundary_zero' in drb1_quality:
            continue
        drb1_q = split(drb1['pred_R']) + split(drb1['pred_D'])
        if any(allele not in drb1_to_dqb1 for allele in drb1_q):
            continue
        proposed = [drb1_to_dqb1[allele] for allele in drb1_q]
        current = split(dqb1['pred_R']) + split(dqb1['pred_D'])
        if proposed == current:
            continue
        # Guard: do not introduce a DQB1 2-field absent from current calls.
        if not set(proposed).issubset(set(current)):
            continue
        truth_r = split(dqb1['truth_R'])
        truth_d = split(dqb1['truth_D'])
        current_score = score(current, truth_r, truth_d)
        proposed_score = score(proposed, truth_r, truth_d)
        delta = proposed_score - current_score
        verdict = 'improve' if delta > 0 else 'regress' if delta < 0 else 'neutral'
        stats['rows'] += 1
        stats[verdict] += 1
        stats['delta'] += delta
        stats['current_score'] += current_score
        stats['proposed_score'] += proposed_score
        out_rows.append({
            'sample': sample,
            'set': dqb1['set'],
            'gene': 'HLA-DQB1',
            'rule': 'DRB1_DQB1_side_copy_no_new_2field',
            'current_score': current_score,
            'proposed_score': proposed_score,
            'delta': delta,
            'verdict': verdict,
            'drb1_quartet': ','.join(drb1_q),
            'drb1_quality': drb1_quality,
            'current_quartet': ','.join(current),
            'proposed_quartet': ','.join(proposed),
            'truth_R': dqb1['truth_R'],
            'truth_D': dqb1['truth_D'],
        })

    fields = ['sample', 'set', 'gene', 'rule', 'current_score', 'proposed_score', 'delta', 'verdict', 'drb1_quartet', 'drb1_quality', 'current_quartet', 'proposed_quartet', 'truth_R', 'truth_D']
    write(args.out, fields, out_rows)
    write(args.summary, ['rule', 'rows', 'improve', 'regress', 'neutral', 'delta', 'current_score', 'proposed_score'], [{
        'rule': 'DRB1_DQB1_side_copy_no_new_2field',
        'rows': stats['rows'],
        'improve': stats['improve'],
        'regress': stats['regress'],
        'neutral': stats['neutral'],
        'delta': stats['delta'],
        'current_score': stats['current_score'],
        'proposed_score': stats['proposed_score'],
    }])
    print(args.summary)


if __name__ == '__main__':
    main()
