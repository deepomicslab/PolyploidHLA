#!/usr/bin/env python3
"""
HLA polyploid haplotype assembly (k=4) from a phased VCF + IMGT allele DB.

Pipeline (per HLA gene):
  1. Slice phased VCF to the gene region.
  2. Group variants by PS (phase set); for each block, build 4 local
     haplotype sequences via `bcftools consensus -H <h>` (correctly
     handles SNPs + InDels + multi-allelics).
  3. Base-level pairwise alignment (parasail Smith-Waterman / semi-global)
     of every block-haplotype to all IMGT alleles of that gene ->
     per-(block, local_hap, allele) match score matrix.
  4. Global assembly with optional 2+2 paired-diploid constraint:
     choose (a_R1, a_R2) for recipient and (a_D1, a_D2) for donor and a
     permutation pi_b for every block to maximize total score.
  5. Emit 4 full gene-level haplotype FASTAs by stitching block sequences
     according to pi_b; gaps between blocks are filled from the chosen
     IMGT allele (base-level aligned to gene reference).

Requires: pysam, bcftools, samtools (in PATH); parasail (pip) or mappy.
"""

import argparse
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict

import pysam

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BUNDLED_SPECHLA = os.path.join(SCRIPT_DIR, "resources", "spechla")
LEGACY_SPECHLA = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "SpecHLA"))
DEFAULT_SPECHLA = os.environ.get(
    "SPECHLA",
    BUNDLED_SPECHLA if os.path.isdir(BUNDLED_SPECHLA) else LEGACY_SPECHLA,
)
DEFAULT_IMGT = os.environ.get(
    "IMGT_HLA_FASTA",
    os.path.join(DEFAULT_SPECHLA, "db", "ref", "hla_gen.format.filter.extend.DRB.no26789.v2.fasta"),
)

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it, **kw):
        return it


PLOIDY = 4


# ---------- I/O helpers ----------

def read_gene_bed(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split()
            out.append((f[0], int(f[1]), int(f[2]), f[3]))
    return out


def fetch_ref(ref, chrom, start, end):
    return ref.fetch(chrom, start, end).upper()


def load_imgt_alleles(fasta_path):
    """Return dict allele_name -> sequence (uppercase, no gaps)."""
    alleles = {}
    name, parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    alleles[name] = "".join(parts).upper().replace("-", "")
                tokens = line[1:].split()
                star = next((t for t in tokens if "*" in t), tokens[0])
                name, parts = star, []
            else:
                parts.append(line)
        if name is not None:
            alleles[name] = "".join(parts).upper().replace("-", "")
    return alleles


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


# ---------- Step 1+2: per-block VCF slice + bcftools consensus ----------

def collect_blocks(vcf_path, chrom, start, end):
    """Return list of (ps, lo, hi, obs) for each PS block in [start,end].

    obs = list of (pos0, gt_tuple, af_obs) for biallelic phased variants where
    FORMAT/AD is present (used by chimerism-aware penalty).
    """
    vcf = pysam.VariantFile(vcf_path)
    blocks = defaultdict(lambda: {"items": [], "obs": []})
    for rec in vcf.fetch(chrom, start, end):
        sample = rec.samples[0]
        gt = sample.get("GT")
        if gt is None or len(gt) != PLOIDY or any(a is None for a in gt):
            continue
        if not sample.phased:
            continue
        ps = sample.get("PS")
        if ps is None:
            ps = rec.pos
        blocks[ps]["items"].append((rec.pos - 1, len(rec.ref)))
        if rec.alts and len(rec.alts) == 1 and set(gt) <= {0, 1}:
            ad = sample.get("AD")
            if ad and len(ad) >= 2 and not any(x is None for x in ad) and sum(ad) > 0:
                af = ad[1] / sum(ad)
                blocks[ps]["obs"].append((rec.pos - 1, tuple(gt), af))
    out = []
    for ps, d in blocks.items():
        items = sorted(d["items"])
        lo = max(start, items[0][0])
        hi = min(end, items[-1][0] + items[-1][1])
        out.append((ps, lo, hi, d["obs"]))
    out.sort(key=lambda x: x[1])
    return out


def collect_haplotag_coverage(bam_path, chrom, blocks_info):
    """For each PS block return [c1,c2,c3,c4] read counts (HP tag 1..4).

    BAM must be produced by `whatshap haplotag --ploidy 4` (sets HP=1..4 and PS).
    Reads outside the block (different PS) are ignored.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = []
    for ps, lo, hi, _obs in blocks_info:
        cov = [0, 0, 0, 0]
        for read in bam.fetch(chrom, lo, hi):
            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            try:
                hp = read.get_tag("HP")
                rps = read.get_tag("PS")
            except KeyError:
                continue
            if rps != ps:
                continue
            if 1 <= hp <= PLOIDY:
                cov[hp - 1] += 1
        out.append(cov)
    return out


def compute_chim_penalty_vaf(blocks_obs, perms24, chi_r):
    """penalty[bi][pi_idx] = sum_v |AF_obs - E[AF | pi, chi_r]| over biallelic vars."""
    chi_d = 1.0 - chi_r
    w = (chi_r / 2.0, chi_r / 2.0, chi_d / 2.0, chi_d / 2.0)
    out = []
    for obs in blocks_obs:
        if not obs:
            out.append([0.0] * len(perms24))
            continue
        pen = []
        for pi in perms24:
            wl = [w[pi[h]] for h in range(PLOIDY)]
            tot = 0.0
            for _pos, gt, af in obs:
                exp = wl[0]*gt[0] + wl[1]*gt[1] + wl[2]*gt[2] + wl[3]*gt[3]
                tot += abs(af - exp)
            pen.append(tot)
        out.append(pen)
    return out


def compute_chim_penalty_cov(blocks_cov, perms24, chi_r):
    """penalty[bi][pi_idx] = sum_g |obs_frac_g - exp_frac_g| over global haps."""
    chi_d = 1.0 - chi_r
    exp = (chi_r / 2.0, chi_r / 2.0, chi_d / 2.0, chi_d / 2.0)
    out = []
    for cov in blocks_cov:
        tot = sum(cov)
        if tot == 0:
            out.append([0.0] * len(perms24))
            continue
        frac = [c / tot for c in cov]
        pen = []
        for pi in perms24:
            # global g gets local h where pi[h]=g
            inv = [0, 0, 0, 0]
            for h, g in enumerate(pi):
                inv[g] = h
            tot_p = abs(frac[inv[0]] - exp[0]) + abs(frac[inv[1]] - exp[1]) \
                  + abs(frac[inv[2]] - exp[2]) + abs(frac[inv[3]] - exp[3])
            pen.append(tot_p)
        out.append(pen)
    return out


def estimate_chimerism_from_vcf(blocks_obs, min_points=8):
    """Estimate the *major* fraction (>=0.5) from biallelic variant AFs.

    Caller decides which population (R or D) corresponds to the major
    fraction (typical priors: allo-HSCT recipient blood -> donor major;
    solid-organ transplant recipient blood -> recipient major).

    For a tetraploid variant with ALT-dosage d in (0,4), the observed AF
    is one of (chi_R/2, chi_D/2) when d=1, and (1-chi_R/2, 1-chi_D/2)
    when d=3. Map both to [0,0.5], 1D k=2 cluster, return major = 2*high.
    """
    pts = []
    for obs in blocks_obs:
        for _pos, gt, af in obs:
            d = sum(gt)
            if d == 1:
                pts.append(af)
            elif d == 3:
                pts.append(1.0 - af)
    pts = [x for x in pts if 0.02 < x < 0.48]
    if len(pts) < min_points:
        return None, len(pts), None
    pts.sort()
    best = (float("inf"), None, None, None)
    for i in range(1, len(pts)):
        a, b = pts[:i], pts[i:]
        ma, mb = sum(a)/len(a), sum(b)/len(b)
        sse = sum((x-ma)**2 for x in a) + sum((x-mb)**2 for x in b)
        if sse < best[0]:
            best = (sse, ma, mb, i)
    _sse, m_low, m_high, split = best
    # Two estimates of the major fraction; trust the one from the high cluster
    # (low cluster sees more reference-bias / dropout noise).
    chi_major = max(2 * m_high, 1 - 2 * m_low)
    chi_major = min(0.99, max(0.51, chi_major))
    return chi_major, len(pts), (m_low, m_high, split)


def slice_vcf_by_ps(vcf_path, chrom, ps, lo, hi, out_vcf):
    """Subset VCF to records in [lo,hi] AND with FORMAT/PS == ps."""
    region = f"{chrom}:{lo+1}-{hi}"
    expr = f"FMT/PS={ps}"
    run(["bcftools", "view", "-r", region, "-i", expr,
         "-Oz", "-o", out_vcf, vcf_path])
    run(["bcftools", "index", "-t", "-f", out_vcf])


def bcftools_consensus_hap(ref_fa, vcf_gz, chrom, lo, hi, hap_idx, sample,
                           mask_bed=None):
    """samtools faidx | bcftools consensus -H <h> -s <sample> [--mask BED].

    SpecHLA-style: when ``mask_bed`` is given, every position listed in the
    BED becomes 'N' in the consensus output. This handles intra-block
    low-coverage holes (e.g. exon-only sequencing where introns are blank).
    """
    region = f"{chrom}:{lo+1}-{hi}"
    p1 = subprocess.Popen(["samtools", "faidx", ref_fa, region],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cons_cmd = ["bcftools", "consensus", "-H", str(hap_idx),
                "-s", sample]
    if mask_bed:
        cons_cmd += ["--mask", mask_bed]
    cons_cmd += [vcf_gz]
    p2 = subprocess.Popen(cons_cmd, stdin=p1.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    out, err = p2.communicate()
    _, err1 = p1.communicate()
    if p2.returncode != 0:
        raise RuntimeError(
            f"bcftools consensus failed (hap {hap_idx}): "
            f"{err.decode(errors='replace')}\n{err1.decode(errors='replace')}")
    seq = []
    for line in out.decode().splitlines():
        if line.startswith(">"):
            continue
        seq.append(line.strip())
    return "".join(seq).upper()


def build_block_haplotypes_via_bcftools(vcf_path, ref_fa, chrom, ps, lo, hi,
                                        sample, tmpdir, mask_bed=None):
    sub_vcf = os.path.join(tmpdir, f"block_{chrom}_{ps}.vcf.gz")
    slice_vcf_by_ps(vcf_path, chrom, ps, lo, hi, sub_vcf)
    seqs = []
    for h in range(1, PLOIDY + 1):
        seqs.append(
            bcftools_consensus_hap(ref_fa, sub_vcf, chrom, lo, hi, h, sample,
                                   mask_bed=mask_bed))
    return seqs


# ---------- Step 3: base-level pairwise alignment ----------

class BaseLevelAligner:
    def __init__(self, backend="parasail", prefilter_top=200):
        self.backend = backend
        self.prefilter_top = prefilter_top
        self._mappy_prefilter = None
        if backend == "parasail":
            import parasail
            self.parasail = parasail
            self.matrix = parasail.matrix_create("ACGT", 2, -3)
            self.gap_open = 5
            self.gap_extend = 2
        elif backend == "mappy":
            import mappy
            self.mappy = mappy
        else:
            raise ValueError(backend)

    def index_alleles(self, allele_dict, tmpdir):
        self._alleles = allele_dict
        if self.backend == "parasail" and self.prefilter_top \
                and len(allele_dict) > self.prefilter_top:
            try:
                import mappy
                fa = os.path.join(tmpdir, "alleles_prefilter.fa")
                with open(fa, "w") as fh:
                    for n, s in allele_dict.items():
                        fh.write(f">{n}\n{s}\n")
                self._mappy_prefilter = mappy.Aligner(
                    fa, preset="asm5", best_n=self.prefilter_top)
                print(f"  [aligner] mappy prefilter built ({len(allele_dict)} alleles, "
                      f"top-{self.prefilter_top} per query)", flush=True)
            except Exception as e:
                print(f"  [aligner] mappy prefilter unavailable ({e}); "
                      f"will run parasail on ALL {len(allele_dict)} alleles",
                      flush=True)
                self._mappy_prefilter = None
        if self.backend == "mappy":
            fa = os.path.join(tmpdir, "alleles.fa")
            with open(fa, "w") as fh:
                for n, s in allele_dict.items():
                    fh.write(f">{n}\n{s}\n")
            self._aln = self.mappy.Aligner(fa, preset="asm5", best_n=50)

    def _candidate_names(self, query):
        if self._mappy_prefilter is None:
            return list(self._alleles.keys())
        seen = {}
        for hit in self._mappy_prefilter.map(query):
            s = hit.mlen - 2 * hit.NM
            if hit.ctg not in seen or s > seen[hit.ctg]:
                seen[hit.ctg] = s
        if not seen:
            return list(self._alleles.keys())
        return [n for n, _ in sorted(seen.items(), key=lambda kv: -kv[1])]

    def score_against_all(self, query, top_k=20, use_prefilter=True):
        if not query:
            return {}
        if self.backend == "parasail":
            if use_prefilter:
                cand = self._candidate_names(query)
            else:
                cand = list(self._alleles.keys())
            scores = {}
            for name in cand:
                ref = self._alleles[name]
                res = self.parasail.sg_dx_trace_striped_16(
                    query, ref, self.gap_open, self.gap_extend, self.matrix)
                scores[name] = res.score
            if top_k and len(scores) > top_k:
                scores = dict(sorted(scores.items(), key=lambda kv: -kv[1])[:top_k])
            return scores
        scores = {}
        for hit in self._aln.map(query):
            s = hit.mlen - 2 * hit.NM
            if hit.ctg not in scores or s > scores[hit.ctg]:
                scores[hit.ctg] = s
        if top_k and len(scores) > top_k:
            scores = dict(sorted(scores.items(), key=lambda kv: -kv[1])[:top_k])
        return scores


# ---------- Step 4: global assembly with 2+2 paired-diploid constraint ----------

def _score_combo(block_scores, combo, perms24, penalty_per_block=None):
    total = 0.0
    chosen_perms = []
    for bi, blk in enumerate(block_scores):
        best_b = float("-inf")
        best_pi = None
        for pi_idx, pi in enumerate(perms24):
            s = sum(blk[h_local].get(combo[pi[h_local]], 0.0)
                    for h_local in range(PLOIDY))
            if penalty_per_block is not None:
                s -= penalty_per_block[bi][pi_idx]
            if s > best_b:
                best_b, best_pi = s, pi
        total += best_b
        chosen_perms.append(best_pi)
    return total, chosen_perms


def assemble_paired_diploid(block_scores, top_n_per_block=10, paired=True,
                            global_pool_cap=30, penalty_per_block=None):
    """
    If paired=True: hap{1,2}=recipient, hap{3,4}=donor; recipient and
        donor each can be homozygous (same allele x2).
    Returns (chosen_alleles, perms_per_block, total_score, assignment).
    """
    # Per-allele global score = sum over blocks of max-over-haps(score)
    # Use this to cap the search pool to the top `global_pool_cap` alleles.
    allele_global = defaultdict(float)
    for blk in block_scores:
        per_allele_max = defaultdict(lambda: float("-inf"))
        for hap_scores in blk:
            for name, val in hap_scores.items():
                if val > per_allele_max[name]:
                    per_allele_max[name] = val
        for name, val in per_allele_max.items():
            allele_global[name] += val

    pool = set()
    for blk in block_scores:
        for hap_scores in blk:
            top = sorted(hap_scores.items(), key=lambda kv: -kv[1])[:top_n_per_block]
            pool.update(name for name, _ in top)
    pool = sorted(pool)
    if len(pool) > global_pool_cap:
        pool = sorted(pool, key=lambda n: -allele_global.get(n, 0.0))[:global_pool_cap]
        pool.sort()
    if len(pool) < 2:
        return None, None, float("-inf"), None

    perms24 = list(itertools.permutations(range(PLOIDY)))
    best = (float("-inf"), None, None)

    if paired:
        diploid_combos = [(a, b) for i, a in enumerate(pool) for b in pool[i:]]
        total_combos = len(diploid_combos) ** 2
        print(f"  [assemble] pool={len(pool)} alleles, "
              f"diploid_pairs={len(diploid_combos)}, "
              f"combos to evaluate={total_combos}", flush=True)
        pbar = tqdm(total=total_combos, desc="  paired-assembly",
                    mininterval=1.0, unit="combo")
        for r in diploid_combos:
            for d in diploid_combos:
                combo = (r[0], r[1], d[0], d[1])
                total, perms = _score_combo(block_scores, combo, perms24,
                                            penalty_per_block)
                if total > best[0]:
                    best = (total, combo, perms)
                pbar.update(1)
        pbar.close()
        assignment = {1: "R", 2: "R", 3: "D", 4: "D"}
    else:
        combos = list(itertools.combinations_with_replacement(pool, PLOIDY))
        print(f"  [assemble] pool={len(pool)} alleles, "
              f"unconstrained combos={len(combos)}", flush=True)
        for combo in tqdm(combos, desc="  assembly",
                          mininterval=1.0, unit="combo"):
            total, perms = _score_combo(block_scores, combo, perms24,
                                        penalty_per_block)
            if total > best[0]:
                best = (total, combo, perms)
        assignment = {h: "?" for h in range(1, PLOIDY + 1)}

    total, combo, perms = best
    return combo, perms, total, assignment


# ---------- Step 5: stitch full haplotypes ----------

def stitch_haplotypes(blocks_meta, block_seqs, perms, alleles_chosen,
                      allele_seqs, gene_start, gene_end, ref, chrom,
                      aligner_backend):
    full = ["" for _ in range(PLOIDY)]
    gene_ref_seq = fetch_ref(ref, chrom, gene_start, gene_end)

    allele_maps = {}
    for name in set(alleles_chosen):
        seq = allele_seqs.get(name, "")
        if not seq:
            allele_maps[name] = None
            continue
        allele_maps[name] = _align_allele_to_ref(seq, gene_ref_seq, aligner_backend)

    def fill(global_h, ref_lo, ref_hi):
        name = alleles_chosen[global_h]
        amap = allele_maps.get(name)
        if amap is None:
            return gene_ref_seq[ref_lo - gene_start: ref_hi - gene_start]
        return _extract_aligned_region(amap, ref_lo - gene_start,
                                       ref_hi - gene_start, gene_ref_seq)

    cursor = gene_start
    for (b_lo, b_hi), seqs, pi in zip(blocks_meta, block_seqs, perms):
        if b_lo > cursor:
            for h in range(PLOIDY):
                full[h] += fill(h, cursor, b_lo)
        inv = [0] * PLOIDY
        for h_local, h_global in enumerate(pi):
            inv[h_global] = h_local
        for h in range(PLOIDY):
            full[h] += seqs[inv[h]]
        cursor = b_hi
    if cursor < gene_end:
        for h in range(PLOIDY):
            full[h] += fill(h, cursor, gene_end)
    return full


def _align_allele_to_ref(allele_seq, ref_seq, backend):
    """Return (ref_aligned, allele_aligned) strings (with '-' gaps)."""
    if backend == "parasail":
        import parasail
        matrix = parasail.matrix_create("ACGT", 2, -3)
        # semi-global: end-gap free on REFERENCE so the allele can be a
        # subregion of the gene span
        res = parasail.sg_dx_trace_striped_16(
            allele_seq, ref_seq, 5, 2, matrix)
        tb = res.traceback
        return (tb.ref, tb.query)
    import mappy
    aln = mappy.Aligner(seq=ref_seq, preset="asm5")
    hit = next(iter(aln.map(allele_seq)), None)
    if hit is None:
        return None
    q = allele_seq if hit.strand == 1 else mappy.revcomp(allele_seq)
    ref_a, all_a = [], []
    rpos, qpos = hit.r_st, hit.q_st
    for length, op in hit.cigar:
        if op in (0, 7, 8):
            ref_a.append(ref_seq[rpos:rpos+length])
            all_a.append(q[qpos:qpos+length])
            rpos += length; qpos += length
        elif op == 1:
            ref_a.append("-" * length); all_a.append(q[qpos:qpos+length])
            qpos += length
        elif op == 2:
            ref_a.append(ref_seq[rpos:rpos+length]); all_a.append("-" * length)
            rpos += length
        elif op in (4, 5):
            qpos += length
    return ("".join(ref_a), "".join(all_a))


def _extract_aligned_region(amap, ref_lo, ref_hi, gene_ref_seq):
    ref_aln, all_aln = amap
    if ref_aln is None:
        return gene_ref_seq[ref_lo:ref_hi]
    out = []
    rcur = 0
    in_region = False
    for rb, ab in zip(ref_aln, all_aln):
        if rb != "-":
            if ref_lo <= rcur < ref_hi:
                in_region = True
                if ab != "-":
                    out.append(ab)
            rcur += 1
        else:
            if in_region and ab != "-":
                out.append(ab)
        if rcur >= ref_hi:
            break
    s = "".join(out)
    return s if s else gene_ref_seq[ref_lo:ref_hi]


def _bam_coverage(bam_path, chrom, start, end):
    """Return numpy-like list of length (end-start) with per-base depth (primary,
    non-dup, non-secondary reads)."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    cov = [0] * (end - start)
    for col in bam.pileup(chrom, start, end, truncate=True,
                          min_base_quality=0, min_mapping_quality=0,
                          stepper="nofilter"):
        # apply our own read filtering
        n = 0
        for r in col.pileups:
            a = r.alignment
            if a.is_unmapped or a.is_secondary or a.is_supplementary or a.is_duplicate:
                continue
            if r.is_del or r.is_refskip:
                continue
            n += 1
        i = col.reference_pos - start
        if 0 <= i < len(cov):
            cov[i] = n
    return cov


def _lowcov_intervals_from_cov(cov, window, min_depth, offset=0):
    """SpecHLA-style sliding-window mask. Returns list of (start, end) tuples
    in absolute coordinates (offset added). Mean depth over `window` bases
    < `min_depth` -> mask.
    """
    n = len(cov)
    intervals = []
    in_low = False
    s = e = 0
    csum = [0] * (n + 1)
    for i in range(n):
        csum[i + 1] = csum[i] + cov[i]
    for i in range(window, n):
        m = (csum[i] - csum[i - window]) / window
        if m < min_depth:
            if not in_low:
                s = i - window
            e = i
            in_low = True
        else:
            if in_low:
                if intervals and s < intervals[-1][1]:
                    intervals[-1] = (intervals[-1][0], e)
                else:
                    intervals.append((s, e))
            in_low = False
    if in_low:
        if intervals and s < intervals[-1][1]:
            intervals[-1] = (intervals[-1][0], e)
        else:
            intervals.append((s, e))
    return [(s + offset, e + offset) for (s, e) in intervals]


def _write_mask_bed(intervals, chrom, path):
    """Write 0-based half-open BED for `bcftools consensus --mask`."""
    with open(path, "w") as fh:
        for s, e in intervals:
            if e > s:
                fh.write(f"{chrom}\t{s}\t{e}\n")


def mask_hap_by_coverage(hap_seq, gene_ref_seq, lowcov, aligner_backend):
    """Replace hap bases whose ref position has cov < K with 'N'.

    lowcov: bool array (len = len(gene_ref_seq)); True means mask.
    Insertion bases (hap base at ref-gap) inherit the mask of the preceding
    ref position.
    """
    if not hap_seq:
        return hap_seq
    amap = _align_allele_to_ref(hap_seq, gene_ref_seq, aligner_backend)
    if amap is None:
        return hap_seq
    ref_aln, hap_aln = amap
    out = []
    rcur = 0
    last_mask = False
    for rb, hb in zip(ref_aln, hap_aln):
        if rb != "-":
            mask = (0 <= rcur < len(lowcov)) and lowcov[rcur]
            last_mask = mask
            rcur += 1
            if hb != "-":
                out.append("N" if mask else hb)
        else:
            if hb != "-":
                out.append("N" if last_mask else hb)
    # if hap had unaligned tail (clipped), keep as-is
    return "".join(out) if out else hap_seq


# ---------- Driver ----------

def gene_to_imgt(gene):
    return gene[4:] if gene.startswith("HLA-") else gene


def process_gene(chrom, gstart, gend, gene, args, ref, sample_name):
    out_gene_dir = os.path.join(args.out, gene)
    os.makedirs(out_gene_dir, exist_ok=True)
    print(f"[{gene}] region {chrom}:{gstart}-{gend}", flush=True)

    blocks_info = collect_blocks(args.vcf, chrom, gstart, gend)
    if not blocks_info:
        print(f"[{gene}] no phased blocks; skip", flush=True)
        return
    print(f"[{gene}] {len(blocks_info)} block(s)", flush=True)

    tmpdir = tempfile.mkdtemp(prefix=f"{gene}_", dir=args.out)
    try:
        # ---- precompute SpecHLA-style low-depth mask BED ----
        mask_bed = None
        lowcov_arr = None
        if args.bam and args.mask_min_depth > 0:
            cov_full = _bam_coverage(args.bam, chrom, gstart, gend)
            intervals = _lowcov_intervals_from_cov(
                cov_full, args.mask_window, args.mask_min_depth, offset=gstart)
            mask_bed = os.path.join(tmpdir, "low_depth.bed")
            _write_mask_bed(intervals, chrom, mask_bed)
            n_low = sum(e - s for s, e in intervals)
            print(f"  [mask] window={args.mask_window} min_depth={args.mask_min_depth} "
                  f"-> {len(intervals)} intervals, {n_low}/{len(cov_full)} ref bp masked",
                  flush=True)
            lowcov_arr = [False] * (gend - gstart)
            for s, e in intervals:
                for i in range(max(0, s - gstart), min(gend - gstart, e - gstart)):
                    lowcov_arr[i] = True

        block_seqs, blocks_meta, blocks_obs = [], [], []
        for bi, (ps, lo, hi, obs) in enumerate(blocks_info):
            try:
                seqs = build_block_haplotypes_via_bcftools(
                    args.vcf, args.ref, chrom, ps, lo, hi, sample_name, tmpdir,
                    mask_bed=mask_bed)
            except Exception as e:
                print(f"[{gene}] block {bi} (PS={ps}) consensus failed: {e}",
                      file=sys.stderr, flush=True)
                continue
            block_seqs.append(seqs)
            blocks_meta.append((lo, hi))
            blocks_obs.append(obs)
            if args.dump_block_fa:
                for h, s in enumerate(seqs):
                    with open(os.path.join(out_gene_dir,
                                           f"block_{bi}_hap{h+1}.fa"), "w") as fh:
                        fh.write(f">block{bi}_hap{h+1} PS={ps} "
                                 f"{chrom}:{lo}-{hi}\n{s}\n")

        if not block_seqs:
            print(f"[{gene}] no usable blocks; skip", flush=True)
            return

        imgt_path = args.imgt
        if not os.path.exists(imgt_path):
            print(f"[{gene}] IMGT db missing: {imgt_path}; skip", flush=True)
            return
        alleles = load_imgt_alleles(imgt_path)
        # Filter to only this gene's alleles. IMGT names look like "A*01:01:01:01",
        # "DRB1*15:01:01:01" etc. gene is e.g. "HLA-A" -> prefix "A*".
        gene_prefix = gene.replace("HLA-", "") + "*"
        alleles = {n: s for n, s in alleles.items() if n.startswith(gene_prefix)}
        print(f"[{gene}] {len(alleles)} IMGT alleles (filtered by prefix '{gene_prefix}')",
              flush=True)
        if not alleles:
            print(f"[{gene}] no alleles match prefix; skip", flush=True)
            return

        aligner = BaseLevelAligner(backend=args.aligner,
                                   prefilter_top=args.prefilter_top)
        aligner.index_alleles(alleles, tmpdir)

        # ---- optional block filtering by length / variant count ----
        keep_idx = list(range(len(block_seqs)))
        if args.min_block_bp > 0 or args.min_block_variants > 0:
            keep_idx = []
            for bi in range(len(block_seqs)):
                lo, hi = blocks_meta[bi]
                if (hi - lo) < args.min_block_bp:
                    continue
                if len(blocks_obs[bi]) < args.min_block_variants:
                    continue
                keep_idx.append(bi)
            dropped = len(block_seqs) - len(keep_idx)
            print(f"  [filter] kept {len(keep_idx)}/{len(block_seqs)} blocks "
                  f"(min_bp={args.min_block_bp}, min_var={args.min_block_variants}); "
                  f"dropped {dropped}", flush=True)
            if not keep_idx:
                print(f"[{gene}] all blocks filtered out; skip", flush=True)
                return
            block_seqs = [block_seqs[i] for i in keep_idx]
            blocks_meta = [blocks_meta[i] for i in keep_idx]
            blocks_obs = [blocks_obs[i] for i in keep_idx]

        block_scores, score_rows = [], []
        n_calls = sum(len(haps) for haps in block_seqs)
        # ref baseline: score(ref_slice, allele) per block; subtracted to
        # remove single-allele HLA_X ref bias.
        ref_baselines = [None] * len(block_seqs)
        if args.ref_baseline != "off":
            n_calls += len(block_seqs)
        t0 = time.time()
        with tqdm(total=n_calls, desc=f"  [{gene}] scoring blocks x haps",
                  unit="hap", mininterval=0.5) as pbar:
            for bi, haps in enumerate(block_seqs):
                if args.ref_baseline != "off":
                    lo, hi = blocks_meta[bi]
                    ref_seq = fetch_ref(ref, chrom, lo, hi)
                    base = aligner.score_against_all(
                        ref_seq, top_k=0, use_prefilter=False)
                    ref_baselines[bi] = base
                    pbar.update(1)
                per_hap = []
                for h, s in enumerate(haps):
                    sc = aligner.score_against_all(s, top_k=args.top_n_per_block * 5)
                    if ref_baselines[bi] is not None:
                        base = ref_baselines[bi]
                        if args.ref_baseline == "sub":
                            sc = {n: v - base.get(n, 0) for n, v in sc.items()}
                        else:  # bonus
                            lam = args.ref_baseline_weight
                            sc = {n: v + lam * max(0.0, v - base.get(n, 0))
                                  for n, v in sc.items()}
                    per_hap.append(sc)
                    for name, val in sc.items():
                        score_rows.append((bi, h + 1, name, val))
                    pbar.update(1)
                block_scores.append(per_hap)
        print(f"  [{gene}] scoring done in {time.time()-t0:.1f}s"
              + (f" (ref-baseline={args.ref_baseline}"
                 + (f", lam={args.ref_baseline_weight}" if args.ref_baseline == 'bonus' else "")
                 + ")" if args.ref_baseline != "off" else ""),
              flush=True)

        with open(os.path.join(out_gene_dir, "match_scores.tsv"), "w") as fh:
            fh.write("block\tlocal_hap\tallele\tscore\n")
            for row in sorted(score_rows, key=lambda r: (r[0], r[1], -r[3])):
                fh.write("\t".join(str(x) for x in row) + "\n")

        # ---- chimerism-aware R/D penalty (paired mode only) ----
        penalty_per_block = None
        chi_r = None
        if args.paired_diploids:
            if args.chimerism is not None:
                try:
                    chi_r = float(args.chimerism)
                except ValueError:
                    chi_r = None
            if chi_r is None or args.chimerism == "auto":
                est, n_pts, info = estimate_chimerism_from_vcf(blocks_obs)
                if est is None:
                    print(f"  [chim] auto-estimate FAILED ({n_pts} usable points; "
                          f"need >=8); supply --chimerism manually", flush=True)
                else:
                    m_low, m_high, split = info
                    # est = major fraction (>=0.5). Map to recipient by prior:
                    if args.recipient_major:
                        chi_r_auto = est
                        prior_msg = "recipient-major (e.g. solid-organ tx)"
                    else:
                        chi_r_auto = 1.0 - est
                        prior_msg = "donor-major (e.g. allo-HSCT recipient blood)"
                    print(f"  [chim] auto-estimate major={est:.3f} "
                          f"(n={n_pts}, cluster means {m_low:.3f}/{m_high:.3f}, "
                          f"split={split}); prior={prior_msg} -> "
                          f"chi_R={chi_r_auto:.3f}", flush=True)
                    if chi_r is None:
                        chi_r = chi_r_auto
        if args.paired_diploids and chi_r is not None:
            if not (0.0 < chi_r < 1.0):
                print(f"[{gene}] WARN: chi_R={chi_r} not in (0,1); ignored",
                      flush=True)
                chi_r = None
        if chi_r is not None:
            perms24 = list(itertools.permutations(range(PLOIDY)))
            pen_v = compute_chim_penalty_vaf(blocks_obs, perms24, chi_r) \
                if args.chim_weight_vaf > 0 else None
            pen_c = None
            if args.haplotag_bam and args.chim_weight_cov > 0:
                blocks_info_4 = [(blocks_info[bi][0], blocks_info[bi][1],
                                  blocks_info[bi][2], blocks_obs[bi])
                                 for bi in range(len(blocks_obs))]
                cov = collect_haplotag_coverage(args.haplotag_bam, chrom,
                                                blocks_info_4)
                pen_c = compute_chim_penalty_cov(cov, perms24, chi_r)
                cov_summary = ", ".join(
                    f"b{bi}:{tuple(c)}" for bi, c in enumerate(cov[:3]))
                print(f"  [chim] haplotag cov (first 3 blocks): {cov_summary}",
                      flush=True)
            if pen_v is not None or pen_c is not None:
                penalty_per_block = []
                for bi in range(len(blocks_obs)):
                    row = [0.0] * len(perms24)
                    if pen_v is not None:
                        for k in range(len(perms24)):
                            row[k] += args.chim_weight_vaf * pen_v[bi][k]
                    if pen_c is not None:
                        for k in range(len(perms24)):
                            row[k] += args.chim_weight_cov * pen_c[bi][k]
                    penalty_per_block.append(row)
                print(f"  [chim] using chi_R={chi_r:.3f} chi_D={1-chi_r:.3f} "
                      f"lam_vaf={args.chim_weight_vaf} "
                      f"lam_cov={args.chim_weight_cov}", flush=True)
                if abs(chi_r - 0.5) < 0.02:
                    print("  [chim] WARN: chi_R~0.5, R/D label still ambiguous",
                          flush=True)

        chosen, perms, total, assignment = assemble_paired_diploid(
            block_scores, top_n_per_block=args.top_n_per_block,
            paired=args.paired_diploids,
            global_pool_cap=args.global_pool_cap,
            penalty_per_block=penalty_per_block)
        if chosen is None:
            print(f"[{gene}] assembly failed", flush=True)
            return

        with open(os.path.join(out_gene_dir, "calls.tsv"), "w") as fh:
            fh.write("global_hap\tassignment\tallele\thap_fraction\ttotal_assembly_score\n")
            for h, a in enumerate(chosen):
                side = assignment[h + 1]
                if chi_r is None or side not in {"R", "D"}:
                    hap_fraction = "NA"
                else:
                    frac = chi_r / 2.0 if side == "R" else (1.0 - chi_r) / 2.0
                    hap_fraction = f"{frac:.6f}"
                fh.write(f"{h+1}\t{side}\t{a}\t{hap_fraction}\t{total:.2f}\n")
        tag = "paired(R/D)" if args.paired_diploids else "unconstrained"
        print(f"[{gene}] {tag}: " +
              ", ".join(f"{assignment[i+1]}:{a}" for i, a in enumerate(chosen)) +
              f" (score={total:.1f})", flush=True)

        full = stitch_haplotypes(blocks_meta, block_seqs, perms, chosen,
                                 alleles, gstart, gend, ref, chrom,
                                 aligner_backend=args.aligner)

        # ---- N-mask between-block fills (block-internal already masked
        #      via `bcftools consensus --mask`).
        if lowcov_arr is not None:
            gene_ref_seq = fetch_ref(ref, chrom, gstart, gend)
            full = [mask_hap_by_coverage(s, gene_ref_seq, lowcov_arr, args.aligner)
                    for s in full]

        for h, s in enumerate(full):
            with open(os.path.join(out_gene_dir, f"hap{h+1}.fa"), "w") as fh:
                fh.write(f">{gene}_hap{h+1} assignment={assignment[h+1]} "
                         f"allele={chosen[h]} {chrom}:{gstart}-{gend}\n{s}\n")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--imgt", default=DEFAULT_IMGT,
                    help="IMGT/HLA allele FASTA used for allele scoring")
    ap.add_argument("--imgt-dir", required=False,
                    help="deprecated; use --imgt")
    ap.add_argument("--gene-bed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default=None,
                    help="sample name in VCF (default: first sample)")
    ap.add_argument("--top-n-per-block", type=int, default=10)
    ap.add_argument("--global-pool-cap", type=int, default=30,
                    help="max number of unique alleles in assembly search pool")
    ap.add_argument("--prefilter-top", type=int, default=200,
                    help="mappy prefilter top-K alleles per query before parasail "
                         "(0 disables)")
    ap.add_argument("--paired-diploids", action="store_true",
                    help="enforce 2+2 recipient/donor split")
    ap.add_argument("--chimerism", default=None,
                    help="recipient fraction in (0,1); or 'auto' to estimate "
                         "from VCF AD. Enables R/D-disambiguating penalties.")
    ap.add_argument("--recipient-major", action="store_true",
                    help="prior for auto-estimate: recipient is the MAJOR "
                         "population (e.g. solid-organ transplant recipient "
                         "blood, donor cfDNA <50%). Default assumes "
                         "recipient is the MINOR population (e.g. allo-HSCT "
                         "recipient blood, donor chimerism >50%).")
    ap.add_argument("--chim-weight-vaf", type=float, default=200.0,
                    help="lambda for VAF-based chimerism penalty (0 disables)")
    ap.add_argument("--chim-weight-cov", type=float, default=200.0,
                    help="lambda for haplotag-coverage chimerism penalty "
                         "(0 disables; requires --haplotag-bam)")
    ap.add_argument("--haplotag-bam", default=None,
                    help="BAM produced by `whatshap haplotag --ploidy 4` with "
                         "HP=1..4 and PS tags; enables coverage-based R/D penalty.")
    ap.add_argument("--bam", default=None,
                    help="BAM aligned to --ref. If set with --mask-min-depth>0, "
                         "haplotype FASTA bases at low-coverage ref positions "
                         "are written as 'N'.")
    ap.add_argument("--mask-min-depth", type=int, default=0,
                    help="min BAM depth (sliding-window mean); ref positions below "
                         "this are masked to N via `bcftools consensus --mask` and "
                         "in between-block fills (0 disables).")
    ap.add_argument("--mask-window", type=int, default=20,
                    help="sliding-window size (bp) for low-depth mask, "
                         "SpecHLA default 20.")
    ap.add_argument("--aligner", choices=["parasail", "mappy"],
                    default="parasail",
                    help="base-level (parasail) or fast (mappy)")
    ap.add_argument("--dump-block-fa", action="store_true")
    ap.add_argument("--genes", default="")
    ap.add_argument("--ref-baseline", choices=["off", "bonus", "sub"],
                    default="bonus",
                    help="Per-block ref-baseline scoring. 'bonus' (default): "
                         "score' = raw + lambda * max(0, raw - baseline) -- "
                         "reward alleles closer to hap than to gene reference, "
                         "never punish. 'sub': raw - baseline (aggressive, may "
                         "hurt non-ref-biased genes). 'off': original raw "
                         "score (subject to single-allele HLA_X ref bias).")
    ap.add_argument("--ref-baseline-weight", type=float, default=1.0,
                    help="lambda weight for the bonus term (only with mode=bonus)")
    ap.add_argument("--min-block-bp", type=int, default=0,
                    help="Drop blocks shorter than this many ref bp from scoring.")
    ap.add_argument("--min-block-variants", type=int, default=0,
                    help="Drop blocks with fewer phased het variants from scoring.")
    args = ap.parse_args()

    for tool in ("bcftools", "samtools"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: '{tool}' not found in PATH")

    os.makedirs(args.out, exist_ok=True)
    ref = pysam.FastaFile(args.ref)

    sample_name = args.sample
    if sample_name is None:
        with pysam.VariantFile(args.vcf) as vf:
            sample_name = list(vf.header.samples)[0]
    print(f"[INFO] sample: {sample_name}", flush=True)

    bed = read_gene_bed(args.gene_bed)
    if args.genes:
        wanted = set(args.genes.split(","))
        bed = [b for b in bed if b[3] in wanted]

    for chrom, gstart, gend, gene in bed:
        try:
            process_gene(chrom, gstart, gend, gene, args, ref, sample_name)
        except Exception as e:
            import traceback
            print(f"[{gene}] ERROR: {e}", file=sys.stderr, flush=True)
            traceback.print_exc()


if __name__ == "__main__":
    main()
