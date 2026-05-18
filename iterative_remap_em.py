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
             chi_step=0.005, chi_prior_lambda=0.0,
             force_names=None, force_quartets=None, restrict_names=None):
    """Search the best (R1,R2,D1,D2) 2-field quartet under a chimerism dose
    model. If per_gene_chi=True, also search chi_r on a grid for each quartet
    and add a soft L1 prior |chi - chi_global|*lambda so the per-gene fit
    cannot drift arbitrarily far from the genome-wide estimate.

    Returns ((R1,R2,D1,D2), score, fitted_chi_r).
    Score is sum |obs - exp| (+ prior penalty when per_gene_chi)."""
    total = sum(counts.values())
    if total == 0: return None, float("inf"), chi_r
    force_names = set(force_names or [])
    force_quartets = list(force_quartets or [])
    restrict_names = set(restrict_names or [])
    if restrict_names:
        items_by_name = {name: counts.get(name, 0.0) for name in restrict_names}
    else:
        items_by_name = dict(sorted(counts.items(), key=lambda kv: -kv[1])[:top_n])
    for name in force_names:
        if name in counts:
            items_by_name[name] = counts[name]
        else:
            items_by_name[name] = 0.0
    items = [(c, n) for c, n in items_by_name.items()
             if n / total > min_frac or c in force_names or c in restrict_names]
    if len(items) < 2: return None, float("inf"), chi_r
    obs_frac = {c: n / total for c, n in items}
    names = [c for c, _ in items]
    if per_gene_chi:
        chi_grid = np.arange(chi_lo, chi_hi + 1e-9, chi_step)
    else:
        chi_grid = np.array([chi_r])
    best = None
    quartet_iter = itertools.product(names, repeat=4)
    seen_quartets = set()
    for R1, R2, D1, D2 in itertools.chain(quartet_iter, force_quartets):
        if (R1, R2) > (R2, R1): continue
        if (D1, D2) > (D2, D1): continue
        if any(name not in names for name in (R1, R2, D1, D2)):
            continue
        key = (R1, R2, D1, D2)
        if key in seen_quartets:
            continue
        seen_quartets.add(key)
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


def logsumexp(vals):
    if not vals:
        return -float("inf")
    m = max(vals)
    if not math.isfinite(m):
        return m
    return m + math.log(sum(math.exp(v - m) for v in vals))


def quartet_l1_score(counts, quartet, chi_r):
    total = sum(counts.values())
    if total == 0 or not quartet:
        return float("inf")
    names = set(counts) | set(quartet)
    obs = {c: counts.get(c, 0.0) / total for c in names}
    nR = defaultdict(int); nD = defaultdict(int)
    for hap, side in zip(quartet, ["R", "R", "D", "D"]):
        (nR if side == "R" else nD)[hap] += 1
    diff = 0.0
    for c in names:
        exp = nR[c] * chi_r / 2.0 + nD[c] * (1.0 - chi_r) / 2.0
        diff += abs(obs.get(c, 0.0) - exp)
    return diff


def reads_to_tf_logweights(reads, safe2name, T=2.0, family_agg="max"):
    rows = []
    for lst in reads.values():
        if not lst:
            continue
        mx = max(a for _, a in lst)
        per_tf = defaultdict(list)
        for c, a in lst:
            if c not in safe2name:
                continue
            per_tf[two_field(safe2name[c])].append((a - mx) / T)
        if not per_tf:
            continue
        if family_agg == "logsum":
            rows.append({tf: logsumexp(vs) for tf, vs in per_tf.items()})
        else:
            rows.append({tf: max(vs) for tf, vs in per_tf.items()})
    return rows


def fit_4hap_read_likelihood(reads, safe2name, counts, chi_r, top_n=12,
                             min_frac=0.002, per_gene_chi=False,
                             chi_lo=0.005, chi_hi=0.5, chi_step=0.01,
                             chi_prior_lambda=0.0, dose_prior_lambda=0.0, T=2.0,
                             log_floor=-25.0, force_names=None,
                             force_quartets=None, family_agg="max",
                             restrict_names=None):
    """Search a quartet by direct read likelihood.

    Each read contributes log(sum_h dose_h * P(read|allele_h)), using BWA AS
    differences as relative log-likelihoods. This keeps multi-mapping
    uncertainty at read level instead of first collapsing it into EM fractions.
    Returns ((R1,R2,D1,D2), nll, chi, gap_to_second).
    """
    total = sum(counts.values())
    if total == 0:
        return None, float("inf"), chi_r, 0.0
    force_names = set(force_names or [])
    force_quartets = list(force_quartets or [])
    restrict_names = set(restrict_names or [])
    if restrict_names:
        items_by_name = {name: counts.get(name, 0.0) for name in restrict_names}
    else:
        items_by_name = dict(sorted(counts.items(), key=lambda kv: -kv[1])[:top_n])
    for name in force_names:
        items_by_name.setdefault(name, counts.get(name, 0.0))
    names = [c for c, n in items_by_name.items()
             if n / total > min_frac or c in force_names or c in restrict_names]
    if len(names) < 2:
        return None, float("inf"), chi_r, 0.0
    rows = reads_to_tf_logweights(reads, safe2name, T=T, family_agg=family_agg)
    if not rows:
        return None, float("inf"), chi_r, 0.0
    name_to_idx = {name: i for i, name in enumerate(names)}
    row_mat = np.full((len(rows), len(names)), log_floor, dtype=np.float64)
    for r, row in enumerate(rows):
        for tf, logw in row.items():
            idx = name_to_idx.get(tf)
            if idx is not None:
                row_mat[r, idx] = logw
    chi_grid = (np.arange(chi_lo, chi_hi + 1e-9, chi_step)
                if per_gene_chi else np.array([chi_r]))
    pairs = list(itertools.combinations_with_replacement(names, 2))
    quartet_iter = itertools.chain(
        ((r1, r2, d1, d2) for r1, r2 in pairs for d1, d2 in pairs),
        (q for q in force_quartets if q),
    )
    seen = set()
    best = None
    second_nll = float("inf")
    n_reads = len(rows)
    for quartet in quartet_iter:
        if any(name not in names for name in quartet):
            continue
        if quartet in seen:
            continue
        seen.add(quartet)
        for chi in chi_grid:
            weights = defaultdict(float)
            weights[quartet[0]] += chi / 2.0
            weights[quartet[1]] += chi / 2.0
            weights[quartet[2]] += (1.0 - chi) / 2.0
            weights[quartet[3]] += (1.0 - chi) / 2.0
            idxs = []
            log_weights = []
            for tf, w in weights.items():
                if w > 0:
                    idxs.append(name_to_idx[tf])
                    log_weights.append(math.log(w))
            sub = row_mat[:, idxs] + np.asarray(log_weights, dtype=np.float64)[None, :]
            mx = sub.max(axis=1)
            ll = float(np.sum(mx + np.log(np.exp(sub - mx[:, None]).sum(axis=1))))
            nll = -ll
            if chi_prior_lambda > 0:
                nll += chi_prior_lambda * n_reads * abs(float(chi) - chi_r)
            if dose_prior_lambda > 0:
                nll += dose_prior_lambda * n_reads * quartet_l1_score(
                    counts, quartet, float(chi)
                )
            if best is None or nll < best[1]:
                if best is not None:
                    second_nll = best[1]
                best = (quartet, float(nll), float(chi))
            elif nll < second_nll:
                second_nll = float(nll)
    if best is None:
        return None, float("inf"), chi_r, 0.0
    gap = second_nll - best[1] if math.isfinite(second_nll) else float("inf")
    return best[0], best[1], best[2], gap


def has_expression_suffix(two_field_name):
    return bool(two_field_name and two_field_name[-1].isalpha() and two_field_name[-1] != "G")


def read_baseline_quartet(path):
    """Read a calls.tsv-shaped baseline file and return 2-field R1,R2,D1,D2."""
    if not path or not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_h = header.index("global_hap")
            i_a = header.index("assignment")
            i_l = header.index("allele")
        except ValueError:
            return None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_h, i_a, i_l):
                continue
            rows.append((f[i_h], f[i_a], two_field(f[i_l])))
    rows.sort(key=lambda r: int(r[0]) if str(r[0]).isdigit() else r[0])
    rs = [a for _h, side, a in rows if side == "R"]
    ds = [a for _h, side, a in rows if side == "D"]
    if len(rs) < 2 or len(ds) < 2:
        return None
    return tuple(rs[:2] + ds[:2])

def read_dpb1_rescue_candidates(path):
    """Read DPB1 allele-family rescue candidates from a manifest file."""
    if not path or not os.path.exists(path):
        return set()
    candidates = set()
    rows = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].strip().split("\t", 1)
                if len(parts) == 2 and parts[0] == "candidate_families":
                    for item in parts[1].split(","):
                        item = item.strip()
                        if item:
                            candidates.add(item)
                continue
            rows.append(line)
    if rows:
        header = rows[0].split("\t")
        for line in rows[1:]:
            values = line.split("\t")
            row = dict(zip(header, values))
            if row.get("candidate") == "1" and row.get("family"):
                candidates.add(row["family"])
    return candidates

def dpb1_rescue_manifest_path(args):
    if args.dpb1_rescue_candidate_manifest:
        return args.dpb1_rescue_candidate_manifest
    return os.path.join(args.fq_dir, "dpb1_family_rescue_manifest.tsv")


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


def quartet_residual_from_counts(counts, quartet, chi_r):
    total = sum(counts.values()) or 1.0
    exp = defaultdict(float)
    for allele in quartet[:2]:
        exp[allele] += chi_r / 2.0
    for allele in quartet[2:4]:
        exp[allele] += (1.0 - chi_r) / 2.0
    names = set(counts) | set(exp)
    return sum(abs(counts.get(name, 0.0) / total - exp.get(name, 0.0)) for name in names)


def collapse_low_recipient_private(counts, winners, chi_r, max_frac=0.02,
                                   dose_ratio=0.20, top_n=4):
    """Replace unsupported recipient-private alleles with supported shared alleles.

    Some genes, especially HLA-C in these mixtures, have high-confidence top
    alleles but the dosage fit can spend a recipient slot on a tiny private
    allele to shave residual. A recipient-private allele with far less support
    than a single recipient haplotype dose is unlikely to be real; try replacing
    it with a donor/top allele only when that improves the same dosage residual.
    """
    if not winners or chi_r <= 0 or chi_r >= 0.5:
        return winners, None, None
    total = sum(counts.values()) or 1.0
    recipient_single_dose = chi_r / 2.0
    q = list(winners)
    donor_set = set(q[2:4])
    candidates = []
    for allele in list(q[2:4]) + [name for name, _count in sorted(counts.items(), key=lambda kv: -kv[1])[:top_n]]:
        if allele not in candidates:
            candidates.append(allele)
    changed = []
    current_residual = quartet_residual_from_counts(counts, q, chi_r)
    for idx in (0, 1):
        allele = q[idx]
        if allele in donor_set:
            continue
        frac = counts.get(allele, 0.0) / total
        if frac >= max_frac or frac >= dose_ratio * recipient_single_dose:
            continue
        best_q = None
        best_replacement = None
        best_residual = current_residual
        for replacement in candidates:
            if replacement == allele:
                continue
            trial = list(q)
            trial[idx] = replacement
            residual = quartet_residual_from_counts(counts, trial, chi_r)
            if residual < best_residual - 1e-9:
                best_q = trial
                best_replacement = replacement
                best_residual = residual
        if best_q is None:
            continue
        changed.append((allele, best_replacement, frac, recipient_single_dose))
        q = best_q
        donor_set = set(q[2:4])
        current_residual = best_residual
    if not changed:
        return winners, None, None
    detail = ", ".join(
        f"{old}->{new} frac={frac:.4f} recipient_single={dose:.4f}"
        for old, new, frac, dose in changed
    )
    return tuple(q), detail, current_residual


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
    ap.add_argument("--baseline-root", default=None,
                    help="sample assembly root containing <gene-lower>/<gene>/calls.tsv; "
                    "baseline alleles are forced into the EM quartet candidate pool")
    ap.add_argument("--dpb1-rescue-candidate-mode", choices=["auto", "off"], default="auto",
                    help="for HLA-DPB1, constrain quartet search to baseline alleles plus "
                    "candidate families from dpb1_family_rescue_manifest.tsv when present")
    ap.add_argument("--dpb1-rescue-candidate-manifest", default=None,
                    help="optional DPB1 family rescue manifest path; defaults to "
                    "<fq-dir>/dpb1_family_rescue_manifest.tsv")
    ap.add_argument("--low-recipient-private-rescue", action="store_true",
                    help="collapse very low-support recipient-private alleles to supported shared alleles")
    ap.add_argument("--low-recipient-private-genes", default="HLA-C",
                    help="comma-separated genes where low recipient-private rescue is enabled")
    ap.add_argument("--low-recipient-private-max-frac", type=float, default=0.02)
    ap.add_argument("--low-recipient-private-dose-ratio", type=float, default=0.20)
    ap.add_argument("--direct-quartet-likelihood", action="store_true",
                    help="choose the 4-hap quartet by direct read-level likelihood instead of EM-fraction L1")
    ap.add_argument("--direct-top-n", type=int, default=12,
                    help="top-N 2-field alleles considered by direct read-level likelihood")
    ap.add_argument("--direct-min-frac", type=float, default=-1.0,
                    help="minimum EM fraction for direct likelihood candidates; <0 reuses --min-frac")
    ap.add_argument("--direct-log-floor", type=float, default=-25.0,
                    help="relative log-likelihood assigned when a read lacks an alignment to a quartet allele")
    ap.add_argument("--direct-per-gene-chi", action="store_true",
                    help="grid-search chi_r for direct read-level likelihood")
    ap.add_argument("--direct-chi-step", type=float, default=0.01)
    ap.add_argument("--direct-chi-prior", type=float, default=0.0,
                    help="NLL prior weight per read on |chi_direct - chi_global|")
    ap.add_argument("--direct-dose-prior", type=float, default=0.0,
                    help="NLL prior weight per read on EM dose L1 for direct likelihood quartets")
    ap.add_argument("--direct-family-agg", choices=["max", "logsum"], default="max",
                    help="how to aggregate sub-allele alignments into a 2-field read likelihood")
    args = ap.parse_args()
    low_private_genes = {g.strip() for g in args.low_recipient_private_genes.split(",") if g.strip()}
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
        with open(os.path.join(args.out_dir, f"{g}.tf_counts.tsv"), "w") as fh:
            fh.write("allele_2field\tem_weight\tfraction\n")
            for tf, n in sorted(tf_counts.items(), key=lambda kv: -kv[1]):
                fh.write(f"{tf}\t{n:.4f}\t{n/total:.8f}\n")
        baseline_quartet = None
        if args.baseline_root:
            baseline_path = os.path.join(args.baseline_root, g.lower(), g, "calls.tsv")
            baseline_quartet = read_baseline_quartet(baseline_path)
            if baseline_quartet:
                print("  baseline quartet candidate: " + ", ".join(baseline_quartet), flush=True)
        fit_force_names = set(baseline_quartet or [])
        fit_restrict_names = set()
        if g == "HLA-DPB1" and args.dpb1_rescue_candidate_mode != "off":
            manifest_path = dpb1_rescue_manifest_path(args)
            dpb1_rescue_candidates = read_dpb1_rescue_candidates(manifest_path)
            if dpb1_rescue_candidates:
                fit_restrict_names = set(fit_force_names) | dpb1_rescue_candidates
                fit_force_names.update(dpb1_rescue_candidates)
                print("  DPB1 rescue-constrained quartet candidates: "
                      + ", ".join(sorted(fit_restrict_names))
                      + f"  (manifest={manifest_path})", flush=True)
        best, diff, fit_chi = fit_4hap(
            dict(tf_counts), args.chi_r,
            top_n=args.top_n, min_frac=args.min_frac,
            per_gene_chi=args.per_gene_chi,
            chi_lo=args.chi_lo, chi_hi=args.chi_hi, chi_step=args.chi_step,
            chi_prior_lambda=args.chi_prior,
            force_names=fit_force_names,
            force_quartets=[baseline_quartet] if baseline_quartet else None,
            restrict_names=fit_restrict_names,
        )
        if best is None:
            print("  fit failed", flush=True); continue
        if args.direct_quartet_likelihood:
            direct_min_frac = args.min_frac if args.direct_min_frac < 0 else args.direct_min_frac
            direct_best, direct_nll, direct_chi, direct_gap = fit_4hap_read_likelihood(
                reads, safe2name, dict(tf_counts), args.chi_r,
                top_n=args.direct_top_n, min_frac=direct_min_frac,
                per_gene_chi=args.direct_per_gene_chi,
                chi_lo=args.chi_lo, chi_hi=args.chi_hi,
                chi_step=args.direct_chi_step,
                chi_prior_lambda=args.direct_chi_prior,
                dose_prior_lambda=args.direct_dose_prior,
                T=args.em_T,
                log_floor=args.direct_log_floor,
                force_names=fit_force_names,
                force_quartets=[baseline_quartet] if baseline_quartet else None,
                family_agg=args.direct_family_agg,
                restrict_names=fit_restrict_names,
            )
            if direct_best is not None:
                best = direct_best
                fit_chi = direct_chi
                diff = quartet_l1_score(dict(tf_counts), best, fit_chi)
                print(f"  direct read-likelihood: nll={direct_nll:.2f} "
                      f"gap={direct_gap:.2f} chi_R_fit={fit_chi:.3f} "
                      f"l1={diff:.3f}", flush=True)
            else:
                print("  direct read-likelihood failed; keeping EM-fraction L1 fit", flush=True)
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
        low_private_detail = None
        if args.low_recipient_private_rescue and g in low_private_genes:
            collapsed, low_private_detail, collapsed_diff = collapse_low_recipient_private(
                dict(tf_counts), tuple(winners), fit_chi,
                max_frac=args.low_recipient_private_max_frac,
                dose_ratio=args.low_recipient_private_dose_ratio,
            )
            winners = list(collapsed)
            if collapsed_diff is not None:
                diff = collapsed_diff
        print(f"  best 4-hap (sumAbsDiff={diff:.3f}, chi_R_fit={fit_chi:.3f}):")
        print(f"    R1={winners[0]}  R2={winners[1]}  "
              f"D1={winners[2]}  D2={winners[3]}", flush=True)
        if rescue_detail:
            print(f"  recipient-minor rescue: {rescue_detail}", flush=True)
        if low_private_detail:
            print(f"  low-recipient-private rescue: {low_private_detail}", flush=True)
        with open(os.path.join(args.out_dir, f"{g}.iterative.tsv"), "w") as fh:
            fh.write("global_hap\tassignment\tallele_2field\thap_fraction\tallele_read_fraction\tem_weight\n")
            for i, (nm, side) in enumerate(zip(winners, ["R","R","D","D"]), 1):
                hap_fraction = fit_chi / 2.0 if side == "R" else (1.0 - fit_chi) / 2.0
                allele_read_fraction = tf_counts.get(nm, 0.0) / total if total else 0.0
                fh.write(f"{i}\t{side}\t{nm}\t{hap_fraction:.6f}\t{allele_read_fraction:.6f}\t{tf_counts.get(nm, 0):.2f}\n")
        # also emit a calls.tsv-shaped file matching hla_polyphase_assemble.py
        # output, so polyphase_v2.sh can drop it in as an override.
        # Use the longest-sub-allele representative for each 2-field winner.
        with open(os.path.join(args.out_dir, f"{g}.calls.tsv"), "w") as fh:
            fh.write("global_hap\tassignment\tallele\thap_fraction\tallele_read_fraction\tallele_read_count\tem_weight\n")
            for i, (nm, side) in enumerate(zip(winners, ["R","R","D","D"]), 1):
                rep_safe = tf_to_safe.get(nm, safe(nm))
                rep = safe2name.get(rep_safe, nm)
                hap_fraction = fit_chi / 2.0 if side == "R" else (1.0 - fit_chi) / 2.0
                allele_read_count = tf_counts.get(nm, 0.0)
                allele_read_fraction = allele_read_count / total if total else 0.0
                fh.write(f"{i}\t{side}\t{rep}\t{hap_fraction:.6f}\t{allele_read_fraction:.6f}\t{allele_read_count:.2f}\t{allele_read_count:.2f}\n")
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
