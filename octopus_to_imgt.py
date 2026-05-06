#!/usr/bin/env python3
"""Method 2: Octopus polyclone phased VCF -> 4 reconstructed haplotype
sequences -> nearest IMGT allele per haplotype.

Input:
  --vcf      octopus polyclone VCF (one contig per file expected)
  --ref      hla.ref.extend.fa (combined per-contig ref)
  --gene     HLA-A / HLA-B / ... (used to filter IMGT db)
  --contig   contig name in ref (e.g. HLA_B). If omitted, use first contig in vcf.
  --bed      gene.spechla.bed (used to take typing region)
  --imgt     IMGT fasta path
Output (stdout):
  hap_idx   freq   ref_id_pct   top_allele1   id1   top_allele2   id2

Strategy:
  - Pick the largest phaseset PS in --gene region.
  - For each phased site within PS, take K-clone GT and ALT/REF.
  - Determine K = max ploidy across phased sites in PS.
  - For sites whose called GT-len < K, broadcast (treat missing as ref).
  - Build K haplotype sequences by editing the ref slice [bed_start, bed_end).
  - Map each via mappy preset 'sr' (or 'asm5') against IMGT subset for the gene.
  - Per-clone frequency = mean of MAP_HF[clone_idx] across phased PASS sites
    that have MAP_HF (fall back to AC-based 1/K for sites without MAP_HF).
"""
import argparse
import sys
import pysam
import mappy as mp


def load_imgt(path, prefix):
    out, n, p = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if n is not None and n.startswith(prefix):
                out[n] = "".join(p).upper().replace("-", "")
            tok = line[1:].split()
            n = next((t for t in tok if "*" in t), tok[0])
            p = []
        else:
            p.append(line)
    if n is not None and n.startswith(prefix):
        out[n] = "".join(p).upper().replace("-", "")
    return out


def parse_bed(path, contig):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == contig:
            return int(f[1]), int(f[2])
    return None, None


def get_ref_seq(ref_path, contig, start, end):
    fa = pysam.FastaFile(ref_path)
    return fa.fetch(contig, start, end).upper()


def collect_phased_sites(vcf_path, contig, bed_start, bed_end):
    """Return list of (pos0, ref, alts, gt_tuple, hf_tuple_or_None, ps).
    Only phased sites with all-defined GT in the gene region."""
    sites = []
    v = pysam.VariantFile(vcf_path)
    for rec in v:
        if contig is not None and rec.chrom != contig:
            continue
        if not (bed_start <= rec.pos - 1 < bed_end):
            continue
        if not rec.alts:
            continue
        # FILTER: take PASS plus harmless tags; reject q10 (low quality)
        flt = list(rec.filter.keys())
        if 'q10' in flt:
            continue
        s = rec.samples[0]
        gt = s.get('GT')
        if gt is None or any(a is None for a in gt):
            continue
        # require phased
        phased = s.phased if hasattr(s, 'phased') else False
        if not phased and len(set(gt)) > 1:
            continue
        ps = s.get('PS')
        hf = s.get('MAP_HF')
        sites.append((rec.pos - 1, rec.ref.upper(),
                      [a.upper() for a in rec.alts], tuple(gt), hf, ps))
    return sites


def pick_best_phaseset(sites):
    # Group by PS; pick the PS whose sites span the largest range AND have
    # ploidy K>=2. Sites with PS=None are kept and broadcast.
    from collections import defaultdict
    by_ps = defaultdict(list)
    for s in sites:
        by_ps[s[5]].append(s)
    if not by_ps:
        return None, []
    # choose PS with most multi-allele (K>=2) phased sites
    def score(ps_sites):
        return sum(1 for x in ps_sites
                   if len(x[3]) >= 2 and len(set(x[3])) > 1)
    cand = [(score(v), k, v) for k, v in by_ps.items() if k is not None]
    if not cand:
        return None, sites
    cand.sort(reverse=True)
    best_ps = cand[0][1]
    return best_ps, by_ps[best_ps]


def build_haplotypes(ref_slice, slice_start, sites, K):
    """Apply edits to ref_slice for each of K haplotypes.
    For unphased AC=K (homozygous alt) sites, apply to all.
    For phased sites with len(gt)<K, broadcast: assume remaining clones are ref.
    """
    # collect per-clone edit list (sorted by pos desc to preserve offsets)
    edits = [list() for _ in range(K)]
    for pos, ref, alts, gt, hf, ps in sites:
        rel = pos - slice_start
        if rel < 0 or rel + len(ref) > len(ref_slice):
            continue
        L = len(gt)
        for c in range(K):
            allele = gt[c] if c < L else 0
            if allele == 0:
                continue
            if allele - 1 >= len(alts):
                continue
            alt = alts[allele - 1]
            edits[c].append((rel, len(ref), alt))
    haps = []
    for c in range(K):
        es = sorted(edits[c], reverse=True)
        s = list(ref_slice)
        for rel, ref_len, alt in es:
            s[rel:rel + ref_len] = list(alt)
        haps.append("".join(s))
    return haps


def freq_per_clone(sites, K):
    """Per-clone average frequency from MAP_HF tag (when available)."""
    import numpy as np
    accum = np.zeros(K, dtype=float)
    n = 0
    for pos, ref, alts, gt, hf, ps in sites:
        if hf is None or len(hf) != K:
            continue
        try:
            v = [float(x) for x in hf]
        except (TypeError, ValueError):
            continue
        accum += np.array(v)
        n += 1
    if n == 0:
        return [1.0 / K] * K
    return (accum / n).tolist()


def map_haps_to_imgt(haps, imgt_db):
    """Build a tiny mappy index over IMGT subset, map each hap, return top-1
    per hap as (allele_name, identity, alen)."""
    # Write a temp fasta
    import tempfile, os
    tf = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False)
    for n, s in imgt_db.items():
        tf.write(f">{n}\n{s}\n")
    tf.close()
    aligner = mp.Aligner(tf.name, preset='asm10', best_n=5)
    if not aligner:
        os.unlink(tf.name)
        raise RuntimeError("mappy aligner failed")
    results = []
    for h in haps:
        best = None  # (identity, name, alen, mlen)
        for hit in aligner.map(h):
            if hit.is_primary:
                ident = hit.mlen / max(1, hit.blen)
                if best is None or ident > best[0] or (
                        ident == best[0] and hit.mlen > best[3]):
                    best = (ident, hit.ctg, hit.blen, hit.mlen)
        if best is None:
            # fall back to non-primary
            for hit in aligner.map(h):
                ident = hit.mlen / max(1, hit.blen)
                if best is None or ident > best[0]:
                    best = (ident, hit.ctg, hit.blen, hit.mlen)
        results.append(best)
    os.unlink(tf.name)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--gene', required=True)
    ap.add_argument('--contig', default=None)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--imgt', default='/data6/wangxuedong/polyploid_hla/SpecHLA/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta')
    args = ap.parse_args()

    contig = args.contig or args.gene.replace('-', '_')
    bed_s, bed_e = parse_bed(args.bed, contig)
    if bed_s is None:
        sys.exit(f"contig {contig} not in bed")
    print(f"# contig={contig} region={bed_s}-{bed_e}", flush=True)

    sites = collect_phased_sites(args.vcf, contig, bed_s, bed_e)
    print(f"# {len(sites)} phased sites in region", flush=True)
    if not sites:
        sys.exit("no phased sites")

    ps, ps_sites = pick_best_phaseset(sites)
    print(f"# best PS={ps} n_sites={len(ps_sites)}", flush=True)

    K = max(len(x[3]) for x in ps_sites)
    if K < 2:
        K = 2
    print(f"# K={K}", flush=True)

    ref_slice = get_ref_seq(args.ref, contig, bed_s, bed_e)
    haps = build_haplotypes(ref_slice, bed_s, ps_sites, K)
    freqs = freq_per_clone(ps_sites, K)

    pfx = args.gene.replace('HLA-', '') + '*'
    imgt_db = load_imgt(args.imgt, pfx)
    print(f"# IMGT alleles for {args.gene}: {len(imgt_db)}", flush=True)

    hits = map_haps_to_imgt(haps, imgt_db)
    print("hap\tfreq\tlen\ttop_allele\tidentity\tblen\tmlen")
    order = sorted(range(K), key=lambda i: -freqs[i])
    for rank, c in enumerate(order):
        h = haps[c]
        if hits[c] is None:
            print(f"{c}\t{freqs[c]:.4f}\t{len(h)}\tNA\tNA\tNA\tNA")
        else:
            ident, name, blen, mlen = hits[c]
            print(f"{c}\t{freqs[c]:.4f}\t{len(h)}\t{name}\t{ident:.4f}\t{blen}\t{mlen}")


if __name__ == '__main__':
    main()
