#!/usr/bin/env python
"""Iterative remap with EM-based read re-weighting (Salmon-like).

Pipeline per gene:
 1. Build aug ref: ONE longest sub-allele per 2-field IMGT name in this gene.
 2. bwa mem -a (report ALL alignments per read; secondary kept).
 3. Parse SAM: collect, for each read, all alignments and their alignment
    score (AS tag). Drop alignments far below the read's best AS.
 4. EM: each read's mass is fractionally distributed over its candidate
    contigs proportional to (abundance[c] * exp((AS[c]-AS_max)/T)). Iterate.
 5. From EM-weighted per-contig mass, fit 4-hap multiset under chimerism
    dose model: AF(nR,nD)=(nR/2)*chi_R+(nD/2)*chi_D. Search constrained
    to top-N contigs by mass.

Usage:
  iterative_remap_em.py --asm-dir D --sample S --fq-dir Q --chi-r 0.27 \
      --gene HLA-A [...] --out-dir O
"""
import argparse, os, sys, subprocess, tempfile, time, math, itertools
from collections import defaultdict
import numpy as np

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


def load_imgt(path=DEFAULT_IMGT):
    out, n, p = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if n is not None: out[n] = "".join(p).upper().replace("-","")
            tok = line[1:].split()
            n = next((t for t in tok if "*" in t), tok[0]); p = []
        else: p.append(line)
    if n: out[n] = "".join(p).upper().replace("-","")
    return out


def filter_gene(db, gene):
    pfx = gene.replace("HLA-", "") + "*"
    return {n: s for n, s in db.items() if n.startswith(pfx)}


def two_field(name):
    a, b = name.split("*")
    fields = b.split(":")
    return f"{a}*{fields[0]}:{fields[1]}" if len(fields) >= 2 else f"{a}*{fields[0]}"


def safe(name):
    return name.replace("*", "_").replace(":", "_")


def build_aug_ref(out_path, name2seq):
    with open(out_path, "w") as fh:
        for n, s in name2seq.items():
            fh.write(f">{n}\n{s}\n")
    subprocess.run(["samtools", "faidx", out_path], check=True)
    subprocess.run(["bwa", "index", out_path], check=True, capture_output=True)


def bwa_mem_all(ref_fa, fq1, fq2, sample, threads, out_sam, max_alt=200):
    """bwa mem -a + -h max_alt : report up to max_alt alignments per read
    in the XA tag-style; with -a everything goes as separate records.
    We instead use -h to bound multi-map explosion via XA tag and parse XA.
    Simpler: just keep bwa mem -a but limit -h."""
    rg = f"@RG\\tID:{sample}\\tSM:{sample}"
    cmd = (f"bwa mem -a -h {max_alt} -t {threads} -U 10000 -L 10000,10000 -R '{rg}' "
           f"{ref_fa} {fq1} {fq2} 2>/dev/null > {out_sam}")
    subprocess.run(["bash", "-c", cmd], check=True)


def parse_sam_to_reads(sam_path, contig_set, min_as_frac=0.95):
    """Return dict: read_id -> list[(contig, AS)] keeping only alignments
    whose AS is >= min_as_frac * best_AS for that read."""
    reads = defaultdict(list)
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11: continue
            flag = int(f[1])
            if flag & 0x4: continue          # unmapped
            qname = f[0]
            # mate distinction: include r1/r2 marker
            mate = "1" if (flag & 0x40) else ("2" if (flag & 0x80) else "0")
            rid = f"{qname}/{mate}"
            ctg = f[2]
            if ctg not in contig_set: continue
            AS = 0
            for tag in f[11:]:
                if tag.startswith("AS:i:"):
                    AS = int(tag[5:]); break
            reads[rid].append((ctg, AS))
    # filter by min_as_frac per read
    out = {}
    for rid, lst in reads.items():
        if not lst: continue
        best = max(a for _, a in lst)
        if best <= 0: continue
        keep = [(c, a) for c, a in lst if a >= best * min_as_frac]
        # collapse duplicates per contig (keep best AS)
        bestc = {}
        for c, a in keep:
            if c not in bestc or a > bestc[c]: bestc[c] = a
        out[rid] = list(bestc.items())
    return out


def run_em(reads, contigs, n_iter=300, tol=1e-6, T=2.0):
    """Vectorized EM via numpy.
    reads: dict rid -> list[(contig, AS)]; we convert to flat arrays."""
    cidx = {c: i for i, c in enumerate(contigs)}
    nC = len(contigs)
    # flatten: row_starts[r], row_ends[r] index into ctgs / aswt
    starts = []
    rows_ctg = []
    rows_w = []  # exp((AS-AS_max)/T) per alignment
    for rid, lst in reads.items():
        if not lst: continue
        starts.append(len(rows_ctg))
        mx = max(a for _, a in lst)
        for c, a in lst:
            rows_ctg.append(cidx[c])
            rows_w.append(math.exp((a - mx) / T))
    starts.append(len(rows_ctg))
    rows_ctg = np.asarray(rows_ctg, dtype=np.int32)
    rows_w   = np.asarray(rows_w,   dtype=np.float64)
    starts   = np.asarray(starts,   dtype=np.int64)
    nR = len(starts) - 1
    # row id per alignment
    row_id = np.repeat(np.arange(nR, dtype=np.int64),
                       np.diff(starts).astype(np.int64))

    theta = np.full(nC, 1.0 / nC)
    for it in range(n_iter):
        u = theta[rows_ctg] * rows_w
        # per-read normalization
        denom = np.zeros(nR)
        np.add.at(denom, row_id, u)
        denom[denom == 0] = 1e-30
        u = u / denom[row_id]
        new = np.zeros(nC) + 1e-12
        np.add.at(new, rows_ctg, u)
        new = new / new.sum()
        delta = float(np.abs(new - theta).sum())
        theta = new
        if delta < tol:
            break
    # final expected counts
    u = theta[rows_ctg] * rows_w
    denom = np.zeros(nR)
    np.add.at(denom, row_id, u)
    denom[denom == 0] = 1e-30
    u = u / denom[row_id]
    counts = np.zeros(nC)
    np.add.at(counts, rows_ctg, u)
    return ({c: float(theta[i]) for c, i in cidx.items()},
            {c: float(counts[i]) for c, i in cidx.items()},
            it + 1)


def fit_4hap(counts, chi_r, top_n=14, min_frac=0.005,
             per_gene_chi=False, chi_lo=0.005, chi_hi=0.5,
             chi_step=0.005, chi_prior_lambda=0.0):
    """Search the best (R1,R2,D1,D2) 2-field quartet under a chimerism dose
    model. If per_gene_chi=True, also search chi_r on a grid for each quartet
    and add a soft L1 prior |chi - chi_global|*lambda so the per-gene fit
    cannot drift arbitrarily far from the genome-wide estimate.

    Returns ((R1,R2,D1,D2), score, fitted_chi_r).
    Score is sum |obs - exp| (+ prior penalty when per_gene_chi)."""
    total = sum(counts.values())
    if total == 0: return None, float("inf"), chi_r
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]
    items = [(c, n) for c, n in items if n / total > min_frac]
    if len(items) < 2: return None, float("inf"), chi_r
    obs_frac = {c: n / total for c, n in items}
    names = [c for c, _ in items]
    if per_gene_chi:
        chi_grid = np.arange(chi_lo, chi_hi + 1e-9, chi_step)
    else:
        chi_grid = np.array([chi_r])
    best = None
    for R1, R2, D1, D2 in itertools.product(names, repeat=4):
        if (R1, R2) > (R2, R1): continue
        if (D1, D2) > (D2, D1): continue
        # Build per-candidate (a, b) so that exp[c] = a[c] + b[c]*chi
        nR = defaultdict(int); nD = defaultdict(int)
        for hap, side in ((R1,"R"),(R2,"R"),(D1,"D"),(D2,"D")):
            (nR if side == "R" else nD)[hap] += 1
        a_arr = np.empty(len(names)); b_arr = np.empty(len(names)); o_arr = np.empty(len(names))
        for i, c in enumerate(names):
            o_arr[i] = obs_frac[c]
            a_arr[i] = nD[c] / 2.0
            b_arr[i] = (nR[c] - nD[c]) / 2.0
        # vectorize over chi_grid: diff_matrix[k] = sum_c |o - a - b*chi_k|
        # shape (len(chi),len(names))
        exp_mat = a_arr[None, :] + b_arr[None, :] * chi_grid[:, None]
        diff_vec = np.abs(o_arr[None, :] - exp_mat).sum(axis=1)
        if chi_prior_lambda > 0:
            diff_vec = diff_vec + chi_prior_lambda * np.abs(chi_grid - chi_r)
        k = int(np.argmin(diff_vec))
        score = float(diff_vec[k]); chi_used = float(chi_grid[k])
        if best is None or score < best[1]:
            best = ((R1, R2, D1, D2), score, chi_used)
    return best


def has_expression_suffix(two_field_name):
    return bool(two_field_name and two_field_name[-1].isalpha() and two_field_name[-1] != "G")


def rescue_recipient_minor(counts, winners, chi_r, min_frac=0.001,
                           min_count=20.0, max_frac=0.08):
    """Recover recipient-minor alleles suppressed by donor-major fitting.

    In donor-major mixtures, the L1 dose fit can prefer a symmetric donor-like
    quartet because the recipient-only allele is underweighted by multi-mapping
    EM. If the fitted recipient pair contains no allele absent from the donor
    pair, but a credible low-frequency non-donor allele exists, report it as the
    second recipient haplotype while keeping the donor major pair unchanged.
    """
    if not winners or chi_r <= 0 or chi_r >= 0.45:
        return winners, None
    total = sum(counts.values()) or 1.0
    r_pair = list(winners[:2])
    d_pair = list(winners[2:])
    donor_set = set(d_pair)
    if set(r_pair) - donor_set:
        return winners, None
    candidates = []
    for name, count in counts.items():
        frac = count / total
        if name in donor_set:
            continue
        if has_expression_suffix(name):
            continue
        if count < min_count or frac < min_frac or frac > max_frac:
            continue
        candidates.append((name, count, frac))
    if not candidates:
        return winners, None
    candidates.sort(key=lambda x: (-x[1], x[0]))
    minor, count, frac = candidates[0]
    shared = max(d_pair, key=lambda n: counts.get(n, 0.0))
    rescued = (shared, minor, d_pair[0], d_pair[1])
    detail = f"R2={minor} weight={count:.1f} ({frac * 100:.2f}%)"
    return rescued, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fq-dir", required=True)
    ap.add_argument("--chi-r", type=float, required=True)
    ap.add_argument("--gene", action="append", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--imgt", default=DEFAULT_IMGT,
                    help="IMGT/HLA FASTA used by SpecHLA")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-as-frac", type=float, default=0.95,
                    help="keep alignments with AS >= frac*best_AS")
    ap.add_argument("--em-T", type=float, default=2.0,
                    help="softmax temperature on AS difference")
    ap.add_argument("--em-iter", type=int, default=300)
    ap.add_argument("--subs-per-2field", type=int, default=5,
                    help="number of longest sub-alleles per 2-field to include "
                    "in the augmented reference (>=1)")
    ap.add_argument("--top-n", type=int, default=25,
                    help="top-N 2-fields by EM mass to enter the 4-hap search")
    ap.add_argument("--min-frac", type=float, default=0.002,
                    help="minimum 2-field fraction to enter the 4-hap search")
    ap.add_argument("--per-gene-chi", action="store_true",
                    help="re-fit chi_r per gene/quartet on a grid (recommended "
                    "when global chi_R is small or per-locus dropout is uneven)")
    ap.add_argument("--chi-lo", type=float, default=0.005)
    ap.add_argument("--chi-hi", type=float, default=0.5)
    ap.add_argument("--chi-step", type=float, default=0.005)
    ap.add_argument("--chi-prior", type=float, default=0.5,
                    help="L1 prior penalty weight on |chi_local - chi_global| "
                    "(only used with --per-gene-chi)")
    ap.add_argument("--recipient-minor-rescue", action="store_true",
                    help="recover a low-frequency recipient-only allele when "
                    "donor-major L1 fitting chooses a symmetric donor-like quartet")
    ap.add_argument("--rescue-min-frac", type=float, default=0.001)
    ap.add_argument("--rescue-min-count", type=float, default=20.0)
    ap.add_argument("--rescue-max-frac", type=float, default=0.08)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"loading IMGT db: {args.imgt}", flush=True)
    db = load_imgt(args.imgt)
    print(f"  {len(db)} alleles", flush=True)
    summary = []
    for g in args.gene:
        short = g.split("-")[1]
        sub_db = filter_gene(db, g)
        # Keep up to K longest sub-alleles per 2-field. Single-rep biases EM
        # toward 2-fields with few sub-alleles (their unique sequence draws all
        # matching reads), starving common 2-fields whose ONE chosen rep may
        # not match this sample's sub-allele exactly. Multiple reps per 2-field
        # let EM pick best within each family; tf_counts aggregation recovers
        # full 2-field mass.
        by_2f = defaultdict(list)
        for nm in sub_db:
            by_2f[two_field(nm)].append(nm)
        cands = []
        for tf, names in by_2f.items():
            names.sort(key=lambda n: -len(sub_db[n]))
            cands.extend(names[:args.subs_per_2field])
        print(f"\n=== {g} ===  candidates={len(cands)} "
              f"({len(by_2f)} 2-fields, <={args.subs_per_2field} reps each)",
              flush=True)
        contigs = {safe(n): sub_db[n] for n in cands}
        safe2name = {safe(n): n for n in cands}
        ref_fa = os.path.join(args.out_dir, f"{g}.aug.fa")
        build_aug_ref(ref_fa, contigs)

        fq1 = os.path.join(args.fq_dir, f"{short}.R1.fq.gz")
        fq2 = os.path.join(args.fq_dir, f"{short}.R2.fq.gz")
        if not os.path.exists(fq1):
            print(f"  fq missing: {fq1}", flush=True); continue
        sam = os.path.join(args.out_dir, f"{g}.aug.sam")
        t0 = time.time()
        bwa_mem_all(ref_fa, fq1, fq2, args.sample, args.threads, sam)
        print(f"  bwa mem -a {time.time()-t0:.1f}s", flush=True)
        t0 = time.time()
        reads = parse_sam_to_reads(sam, set(contigs), args.min_as_frac)
        print(f"  parsed {len(reads)} reads, "
              f"avg multi-map={sum(len(v) for v in reads.values())/max(1,len(reads)):.1f}  "
              f"({time.time()-t0:.1f}s)", flush=True)
        t0 = time.time()
        theta, counts, iters = run_em(reads, contigs, n_iter=args.em_iter,
                                      T=args.em_T)
        print(f"  EM converged in {iters} iters ({time.time()-t0:.1f}s)",
              flush=True)
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:12]
        total = sum(counts.values()) or 1
        for c, n in top:
            print(f"    {safe2name[c]}: weight={n:.1f} ({n/total*100:.2f}%)")
        # roll up to 2-field — sibling alleles steal mass from each other
        # because their sequences are nearly identical. The chimerism dose model
        # is naturally a 2-field statement (truth resolves at 2-field anyway).
        tf_counts = defaultdict(float)
        tf_to_safe = {}  # representative safe contig per 2-field
        for c, n in counts.items():
            tf = two_field(safe2name[c])
            tf_counts[tf] += n
            if tf not in tf_to_safe or n > counts.get(tf_to_safe[tf], 0):
                tf_to_safe[tf] = c
        print(f"  --- rolled to 2-field (top 10) ---", flush=True)
        for tf, n in sorted(tf_counts.items(), key=lambda kv: -kv[1])[:10]:
            print(f"    {tf}: weight={n:.1f} ({n/total*100:.2f}%)")
        best, diff, fit_chi = fit_4hap(
            dict(tf_counts), args.chi_r,
            top_n=args.top_n, min_frac=args.min_frac,
            per_gene_chi=args.per_gene_chi,
            chi_lo=args.chi_lo, chi_hi=args.chi_hi, chi_step=args.chi_step,
            chi_prior_lambda=args.chi_prior,
        )
        if best is None:
            print("  fit failed", flush=True); continue
        winners = list(best)  # already 2-field strings
        rescue_detail = None
        if args.recipient_minor_rescue:
            rescued, rescue_detail = rescue_recipient_minor(
                dict(tf_counts), tuple(winners), args.chi_r,
                min_frac=args.rescue_min_frac,
                min_count=args.rescue_min_count,
                max_frac=args.rescue_max_frac,
            )
            winners = list(rescued)
        print(f"  best 4-hap (sumAbsDiff={diff:.3f}, chi_R_fit={fit_chi:.3f}):")
        print(f"    R1={winners[0]}  R2={winners[1]}  "
              f"D1={winners[2]}  D2={winners[3]}", flush=True)
        if rescue_detail:
            print(f"  recipient-minor rescue: {rescue_detail}", flush=True)
        with open(os.path.join(args.out_dir, f"{g}.iterative.tsv"), "w") as fh:
            fh.write("global_hap\tassignment\tallele_2field\tem_weight\n")
            for i, (nm, side) in enumerate(zip(winners, ["R","R","D","D"]), 1):
                fh.write(f"{i}\t{side}\t{nm}\t{tf_counts.get(nm, 0):.2f}\n")
        # also emit a calls.tsv-shaped file matching hla_polyphase_assemble.py
        # output, so polyphase_v2.sh can drop it in as an override.
        # Use the longest-sub-allele representative for each 2-field winner.
        with open(os.path.join(args.out_dir, f"{g}.calls.tsv"), "w") as fh:
            fh.write("global_hap\tassignment\tallele\tem_weight\n")
            for i, (nm, side) in enumerate(zip(winners, ["R","R","D","D"]), 1):
                rep_safe = tf_to_safe.get(nm, safe(nm))
                rep = safe2name.get(rep_safe, nm)
                fh.write(f"{i}\t{side}\t{rep}\t{tf_counts.get(nm, 0):.2f}\n")
        # per-gene summary line for downstream gating
        with open(os.path.join(args.out_dir, f"{g}.summary.tsv"), "w") as fh:
            fh.write("gene\tsum_abs_diff\tn_reads\ttop_frac\tchi_r_fit\n")
            top_frac = max(tf_counts.values()) / total if total else 0.0
            fh.write(f"{g}\t{diff:.4f}\t{len(reads)}\t{top_frac:.4f}\t{fit_chi:.4f}\n")
        # delete sam to save space
        os.unlink(sam)
        summary.append((g, winners, diff))
    print("\n=== summary ===", flush=True)
    for g, w, d in summary:
        print(f"  {g} (sumAbsDiff={d:.3f}): " + ", ".join(
            [f"{['R','R','D','D'][i]}:{w[i]}" for i in range(4)]))


if __name__ == "__main__":
    main()
