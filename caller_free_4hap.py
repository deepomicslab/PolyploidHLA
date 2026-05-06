#!/usr/bin/env python3
"""Method 3 (caller-free per-site AF -> IMGT 4-hap fit).

Inputs:
  --pc-vcf      pooled-continuous freebayes VCF (e.g. spechla_out/<S>/<S>.pooled_continuous.vcf.gz)
  --gene        e.g. HLA-B
  --contig      contig in pc-vcf and ref (default: HLA_<gene>)
  --ref         hla.ref.extend.fa
  --imgt        IMGT db fasta
  --bed         gene.spechla.bed (optional, for region restriction)
  --chi-r       recipient chimerism (default: read from --pc-log)
  --pc-log      estimate_chi_pooled.py output (to read chi_R automatically)
  --top-n       restrict to top-N IMGT 2-field alleles by mappable-site coverage
  --candidate-allowlist  comma-separated list, restrict R/D candidates (debug)
Output (stdout):
  rank  R1  R2  D1  D2  score  n_inform_sites  R-frac  D-frac
"""
import argparse, sys, os, tempfile
from collections import defaultdict
import numpy as np
import pysam
import mappy as mp


def two_field(name):
    a, b = name.split("*")
    f = b.split(":")
    return f"{a}*{f[0]}:{f[1]}" if len(f) >= 2 else f"{a}*{f[0]}"


def load_imgt(path, gene_prefix):
    out, n, p = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if n is not None and n.startswith(gene_prefix):
                out[n] = "".join(p).upper().replace("-", "")
            tok = line[1:].split()
            n = next((t for t in tok if "*" in t), tok[0])
            p = []
        else:
            p.append(line)
    if n is not None and n.startswith(gene_prefix):
        out[n] = "".join(p).upper().replace("-", "")
    return out


def parse_bed(path, contig):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == contig:
            return int(f[1]), int(f[2])
    return None, None


def collect_obs_af(vcf_path, contig, region):
    """Return dict (pos0, ref, alt) -> AF (AO/(RO+AO)) at biallelic sites."""
    s, e = region
    obs = {}
    v = pysam.VariantFile(vcf_path)
    for rec in v:
        if rec.chrom != contig: continue
        if not (s <= rec.pos - 1 < e): continue
        if not rec.alts or len(rec.alts) != 1: continue
        smp = rec.samples[0]
        ro = smp.get('RO')
        ao = smp.get('AO')
        if ro is None or ao is None: continue
        ao = ao[0] if hasattr(ao, '__len__') else ao
        if ao is None: continue
        ro = 0 if ro is None else ro
        dp = ro + ao
        if dp < 30: continue
        af = ao / dp
        if af <= 0.005 or af >= 0.995: continue
        obs[(rec.pos - 1, rec.ref.upper(), rec.alts[0].upper())] = (af, dp)
    return obs


def imgt_genotypes_at_sites(imgt_seqs, ref_seq, ref_offset, sites):
    """For each IMGT 2-field representative (pick longest seq per 2-field),
    align to ref_seq with mappy, then for each site (pos0_genome, ref, alt)
    determine if the allele has the alt (1) or ref (0) or unknown (-1).

    Returns: {two_field_name: np.array shape (n_sites,) int8 in {0,1,-1}}
    """
    # Pick longest sub-allele per 2-field as representative
    by2f = defaultdict(list)
    for n, s in imgt_seqs.items():
        by2f[two_field(n)].append((len(s), n, s))
    reps = {}
    for tf, lst in by2f.items():
        lst.sort(reverse=True)
        reps[tf] = lst[0][2]

    # Build mappy index over ref_seq
    tf_ref = tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False)
    tf_ref.write(f">REF\n{ref_seq}\n")
    tf_ref.close()
    aligner = mp.Aligner(tf_ref.name, preset='asm10')
    if not aligner:
        os.unlink(tf_ref.name)
        raise RuntimeError("mappy index failed")

    n_sites = len(sites)
    pos_list = [p for (p, r, a) in sites]
    ref_list = [r for (p, r, a) in sites]
    alt_list = [a for (p, r, a) in sites]

    out = {}
    for tf, seq in reps.items():
        gt = np.full(n_sites, -1, dtype=np.int8)
        best_hit = None
        for hit in aligner.map(seq):
            if best_hit is None or hit.mlen > best_hit.mlen:
                best_hit = hit
        if best_hit is None:
            out[tf] = gt
            continue
        # Walk CIGAR to map ref pos -> query pos
        # hit.r_st is 0-based on ref_seq; we need genome pos = ref_offset + r_pos
        r_pos = best_hit.r_st
        q_pos = best_hit.q_st if best_hit.strand > 0 else (len(seq) - best_hit.q_en)
        # Use simpler approach: reconstruct alignment via cs string when available
        # Instead, build pos map by walking cigar
        cigar = best_hit.cigar  # list of (length, op) where op 0=M,1=I,2=D,4=S
        # Build genome_pos -> query allele dict for this allele
        rp = best_hit.r_st
        qp = best_hit.q_st
        ref_to_q = {}
        for length, op in cigar:
            if op == 0 or op == 7 or op == 8:  # M/=/X
                for i in range(length):
                    ref_to_q[rp + i] = qp + i
                rp += length; qp += length
            elif op == 1:  # I (in query relative to ref)
                qp += length
            elif op == 2 or op == 3:  # D / N
                rp += length
            elif op == 4 or op == 5:  # S / H
                pass

        for i, (pos, r_allele, a_allele) in enumerate(sites):
            rp_local = pos - ref_offset
            qpi = ref_to_q.get(rp_local)
            if qpi is None or qpi >= len(seq): continue
            # SNV only
            if len(r_allele) == 1 and len(a_allele) == 1:
                base = seq[qpi]
                if best_hit.strand < 0:
                    # rare: map flipped; bail
                    continue
                if base == r_allele:
                    gt[i] = 0
                elif base == a_allele:
                    gt[i] = 1
                # else: another base, leave -1
            else:
                # indel: skip for now (caller-free indel matching is fragile)
                continue
        out[tf] = gt
    os.unlink(tf_ref.name)
    return out


def fit_4hap(obs_af, obs_dp, geno_mat, allele_names, chi_r,
             top_n=None, return_topk=10):
    """Brute-force enumerate 4-tuples (R1,R2,D1,D2) (ordered as set: R-pair
    + D-pair, both unordered within pair).

    Pred AF at site i: 0.5*(geno[R1,i]+geno[R2,i])*chi_r + 0.5*(geno[D1,i]+geno[D2,i])*chi_d
    Score: weighted sum |obs - pred| over informative sites (geno != -1 for ALL 4 candidates).

    To bound enum: pick top_n alleles by depth-weighted "AF-explained" score.
    """
    chi_d = 1.0 - chi_r
    n_alleles = len(allele_names)
    n_sites = len(obs_af)
    obs_af = np.asarray(obs_af, dtype=float)
    obs_dp = np.asarray(obs_dp, dtype=float)
    G = np.stack([geno_mat[a] for a in allele_names])  # (A, S)
    # Per-allele score: how well the allele's genotype matches observed AF.
    # match = sum_i w_i * (g==1 ? obs_af[i] : g==0 ? 1-obs_af[i] : 0)
    obs_af_arr = np.asarray(obs_af, dtype=float)
    w_pre = np.minimum(np.asarray(obs_dp, dtype=float), 200.0)
    match_score = np.zeros(n_alleles)
    for i in range(n_alleles):
        g = G[i]
        s = np.where(g == 1, obs_af_arr, np.where(g == 0, 1.0 - obs_af_arr, 0.0))
        match_score[i] = (w_pre * s).sum()
    if top_n and top_n < n_alleles:
        top_idx = np.argsort(-match_score)[:top_n]
    else:
        top_idx = np.arange(n_alleles)
    G_top = G[top_idx]
    names_top = [allele_names[i] for i in top_idx]
    A = len(top_idx)
    print(f"# enumerating {A} alleles -> ~{A*A*A*A//4} 4-tuples (with sym reduction)", flush=True)

    # weights: clipped depth (cap at 200 to avoid one site dominating)
    w = np.minimum(obs_dp, 200.0)
    w_sum = w.sum() + 1e-9

    best = []  # list of (score, R1, R2, D1, D2)
    # Symmetry: R1<=R2, D1<=D2; (R-pair vs D-pair) — try both orientations and keep one
    for r1 in range(A):
        gr1 = G_top[r1]
        for r2 in range(r1, A):
            gr2 = G_top[r2]
            r_sum = (gr1 + gr2).astype(float)
            r_valid = (gr1 != -1) & (gr2 != -1)
            for d1 in range(A):
                gd1 = G_top[d1]
                for d2 in range(d1, A):
                    gd2 = G_top[d2]
                    d_sum = (gd1 + gd2).astype(float)
                    valid = r_valid & (gd1 != -1) & (gd2 != -1)
                    if valid.sum() < 5:
                        continue
                    pred = 0.5 * r_sum * chi_r + 0.5 * d_sum * chi_d
                    diff = np.abs(obs_af - pred)
                    score = (w * diff * valid).sum() / (w * valid).sum()
                    if len(best) < return_topk:
                        best.append((score, r1, r2, d1, d2, int(valid.sum())))
                        best.sort()
                    elif score < best[-1][0]:
                        best[-1] = (score, r1, r2, d1, d2, int(valid.sum()))
                        best.sort()
    out = []
    for sc, r1, r2, d1, d2, ns in best:
        out.append((sc, names_top[r1], names_top[r2], names_top[d1],
                    names_top[d2], ns))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pc-vcf', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--gene', required=True)
    ap.add_argument('--contig', default=None)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--imgt', default='/data6/wangxuedong/polyploid_hla/SpecHLA/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta')
    ap.add_argument('--chi-r', type=float, required=True)
    ap.add_argument('--top-n', type=int, default=40)
    ap.add_argument('--allowlist', default=None)
    args = ap.parse_args()

    contig = args.contig or args.gene.replace('-', '_')
    bed_s, bed_e = parse_bed(args.bed, contig)
    if bed_s is None:
        sys.exit(f"contig {contig} not in bed")
    print(f"# contig={contig} region={bed_s}-{bed_e} chi_R={args.chi_r}", flush=True)

    obs = collect_obs_af(args.pc_vcf, contig, (bed_s, bed_e))
    sites = sorted(obs.keys())
    obs_af = [obs[k][0] for k in sites]
    obs_dp = [obs[k][1] for k in sites]
    print(f"# {len(sites)} biallelic AF sites in region", flush=True)
    if len(sites) < 5:
        sys.exit("too few sites")

    fa = pysam.FastaFile(args.ref)
    ref_seq = fa.fetch(contig, bed_s, bed_e).upper()

    pfx = args.gene.replace('HLA-', '') + '*'
    imgt_seqs = load_imgt(args.imgt, pfx)
    print(f"# IMGT alleles: {len(imgt_seqs)}", flush=True)

    geno = imgt_genotypes_at_sites(imgt_seqs, ref_seq, bed_s, sites)
    if args.allowlist:
        allow = set(args.allowlist.split(','))
        geno = {k: v for k, v in geno.items() if k in allow}
    print(f"# 2-field alleles after dedup: {len(geno)}", flush=True)
    names = sorted(geno.keys())

    res = fit_4hap(obs_af, obs_dp, geno, names, args.chi_r,
                   top_n=args.top_n, return_topk=15)
    print("rank\tscore\tR1\tR2\tD1\tD2\tn_sites")
    for i, (sc, r1, r2, d1, d2, ns) in enumerate(res):
        print(f"{i+1}\t{sc:.5f}\t{r1}\t{r2}\t{d1}\t{d2}\t{ns}")


if __name__ == '__main__':
    main()
