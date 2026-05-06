#!/usr/bin/env python3
"""Reassign 4-ploid GT in a freebayes VCF using observed AD and chimerism prior.

Model: 4 haplotypes split into 2 recipient (R) + 2 donor (D).
With chimerism fraction chi_R (recipient read fraction) and chi_D = 1 - chi_R,
the expected ALT allele frequency for a site with nR alt-carrying R-haps and
nD alt-carrying D-haps is:
    AF(nR, nD) = (nR/2)*chi_R + (nD/2)*chi_D    , nR in {0,1,2}, nD in {0,1,2}

We pick (nR, nD) minimizing |obs_AF - AF(nR, nD)|, breaking ties by:
  (1) prefer larger total alt count  (avoid 0/0/0/0 spurious calls)
  (2) prefer balanced split

We then write GT = unphased "a/a/a/a" with 'nR' alts in R-slots and 'nD' alts
in D-slots; whatshap polyphase will phase the haplotype assignment.

Usage:
    reassign_gt_chimeric.py --vcf in.vcf.gz --chi-r 0.27 --out out.vcf.gz \
        [--min-depth 10] [--min-obs-af 0.05] [--drop-fp-af 0.05]
"""
import argparse, gzip, os, sys, subprocess, tempfile

def best_nrnd(obs_af, chi_r):
    chi_d = 1.0 - chi_r
    best = None  # (delta, -tot_alt, nR, nD)
    for nR in (0, 1, 2):
        for nD in (0, 1, 2):
            exp = 0.5 * nR * chi_r + 0.5 * nD * chi_d
            d = abs(obs_af - exp)
            tot = nR + nD
            key = (d, -tot)
            if best is None or key < best[0]:
                best = (key, nR, nD)
    _, nR, nD = best
    return nR, nD

def gt_string(nR, nD):
    # 4 unphased slots: R1 R2 D1 D2
    slots = ["1"] * nR + ["0"] * (2 - nR) + ["1"] * nD + ["0"] * (2 - nD)
    return "/".join(slots)

def open_in(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def write_vcf(in_path, out_path, chi_r, min_depth, min_obs_af, drop_fp_af):
    chi_d = 1.0 - chi_r
    n_in = n_out = n_drop_lowdepth = n_drop_fp = n_changed = 0
    out_text = []
    with open_in(in_path) as fh:
        for line in fh:
            if line.startswith("##"):
                out_text.append(line)
                continue
            if line.startswith("#CHROM"):
                # add INFO line for our annotation
                out_text.append('##INFO=<ID=CHI_AF,Number=1,Type=Float,Description="Best-fit ALT AF under chimerism model">\n')
                out_text.append('##INFO=<ID=CHI_NR,Number=1,Type=Integer,Description="Recipient ALT-hap count under chimerism model">\n')
                out_text.append('##INFO=<ID=CHI_ND,Number=1,Type=Integer,Description="Donor ALT-hap count under chimerism model">\n')
                out_text.append('##FORMAT=<ID=ORIG_GT,Number=1,Type=String,Description="Original freebayes GT before chimeric reassignment">\n')
                out_text.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                out_text.append(line)
                continue
            n_in += 1
            chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample = f[:10]
            # only act on bi-allelic SNV/indel
            if "," in alt:
                out_text.append(line)
                n_out += 1
                continue
            fmt_keys = fmt.split(":")
            sv = sample.split(":")
            d = dict(zip(fmt_keys, sv))
            ad = d.get("AD", "")
            if not ad or "," not in ad:
                out_text.append(line)
                n_out += 1
                continue
            try:
                ad_list = [int(x) for x in ad.split(",")]
            except ValueError:
                out_text.append(line); n_out += 1; continue
            tot = sum(ad_list)
            if tot < min_depth:
                n_drop_lowdepth += 1
                continue
            obs_af = ad_list[1] / tot if tot else 0.0
            # drop obvious FP (very low AF)
            if obs_af < drop_fp_af:
                n_drop_fp += 1
                continue
            nR, nD = best_nrnd(obs_af, chi_r)
            if nR == 0 and nD == 0:
                # would produce 0/0/0/0; skip as no-call
                n_drop_fp += 1
                continue
            new_gt = gt_string(nR, nD)
            orig_gt = d.get("GT", "./././.")
            # update GT, append ORIG_GT
            d["GT"] = new_gt
            d["ORIG_GT"] = orig_gt
            new_keys = list(fmt_keys)
            if "ORIG_GT" not in new_keys:
                new_keys.append("ORIG_GT")
            new_fmt = ":".join(new_keys)
            new_sv = ":".join(d.get(k, ".") for k in new_keys)
            # annotate INFO
            exp_af = 0.5 * nR * chi_r + 0.5 * nD * chi_d
            extra = f"CHI_AF={exp_af:.3f};CHI_NR={nR};CHI_ND={nD}"
            new_info = info + ";" + extra if info != "." else extra
            new_line = "\t".join([chrom, pos, vid, ref, alt, qual, flt, new_info, new_fmt, new_sv]) + "\n"
            out_text.append(new_line)
            if new_gt != orig_gt:
                n_changed += 1
            n_out += 1

    if out_path.endswith(".gz"):
        # write to a temp uncompressed file in the same dir as output, then bgzip
        out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
        os.makedirs(out_dir, exist_ok=True)
        tmp = tempfile.NamedTemporaryFile("wt", delete=False, suffix=".vcf", dir=out_dir)
        tmp.writelines(out_text); tmp.close()
        subprocess.check_call(["bgzip", "-f", tmp.name])
        os.replace(tmp.name + ".gz", out_path)
        subprocess.check_call(["tabix", "-f", "-p", "vcf", out_path])
    else:
        with open(out_path, "wt") as fo:
            fo.writelines(out_text)
    sys.stderr.write(
        f"[reassign_gt] in={n_in} kept={n_out} changed_GT={n_changed} "
        f"dropped_lowdepth(<{min_depth})={n_drop_lowdepth} dropped_FP(AF<{drop_fp_af})={n_drop_fp}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chi-r", type=float, required=True, help="recipient read fraction (0..1)")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-obs-af", type=float, default=0.05, help="(unused, for symmetry)")
    ap.add_argument("--drop-fp-af", type=float, default=0.05, help="drop sites with obs_AF below this (FP control)")
    args = ap.parse_args()
    write_vcf(args.vcf, args.out, args.chi_r, args.min_depth, args.min_obs_af, args.drop_fp_af)

if __name__ == "__main__":
    main()
