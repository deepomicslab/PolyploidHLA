# Polyploid HLA Typing — `scripts/`

End-to-end pipeline for **chimeric (k=4) HLA typing** of allo-HSCT and
solid-organ transplant samples from short-read FASTQs. Outputs 4 haplotype
sequences per gene tagged `R`(ecipient) / `D`(onor).

- Algorithm reference: [../PIPELINE.md](../PIPELINE.md)
- Installation & environment: [INSTALL.md](INSTALL.md)

---

## 1. Files

| File | Role |
| ---- | ---- |
| `polyphase_v2.sh`           | **driver — run this** |
| `hla_polyphase_assemble.py` | baseline typing engine (whatshap polyphase + IMGT scoring) |
| `reassign_gt_chimeric.py`   | χ-aware GT correction before phasing |
| `estimate_chi_pooled.py`    | pooled-continuous χ_R estimator |
| `iterative_remap_em.py`     | EM refinement (Salmon-style read remap) |
| `aggregate_calls.py`        | merges per-gene `calls.tsv` into one summary table |
| `gene.spechla.bed`          | per-gene typing region on `hla.ref.extend.fa` |
| `environment.yml`           | conda environment spec |
| `octopus_to_imgt.py`, `caller_free_4hap.py` | rejected alternatives, kept for reference |

---

## 2. Install

See [INSTALL.md](INSTALL.md). Short version:

```bash
conda env create -f scripts/environment.yml
conda activate polyploid-hla
# Then place / symlink the SpecHLA database at <repo>/SpecHLA, or
# export SPECHLA=/path/to/SpecHLA before running.
```

---

## 3. Quick start

From the repository root (the directory containing `scripts/`):

```bash
FQ1=/path/sample.R1.fq.gz \
FQ2=/path/sample.R2.fq.gz \
SAMPLE=mySample \
RECIPIENT_MAJOR=0 \
bash scripts/polyphase_v2.sh
```

**Final result — one file, all 6 genes:**

```bash
column -t asm_v2/mySample/mySample.final_calls.tsv
```

Example:

```
sample      gene      R1            R2            D1            D2            source
mySample    HLA-A     A*01:01:01:01 A*23:01:01:01 A*01:01:01:01 A*29:02:01:02 em-refined
mySample    HLA-B     B*08:01:01:01 B*44:03:01:01 B*08:01:01:01 B*45:01:01:01 em-refined
...
```

The per-gene FASTAs (`hap{1..4}.fa`) and raw `calls.tsv` are still kept under
`asm_v2/<SAMPLE>/<gene_lc>/<HLA-X>/` for inspection.

---

## 4. Required inputs

| Var | Meaning |
| --- | ------- |
| `FQ1`, `FQ2`        | paired short-read FASTQs (gz ok) |
| `SAMPLE`            | sample id (used for output dirs) |
| `RECIPIENT_MAJOR`   | `0` = donor major (post-HSCT blood, default); `1` = recipient major (solid-organ / pre-transplant) |

Optional environment / database overrides:

| Var | Default | Meaning |
| --- | ------- | ------- |
| `SPECHLA`  | `<repo>/SpecHLA`           | SpecHLA install root |
| `PYBIN`    | first `python` on PATH     | python binary |
| `WHATSHAP` | first `whatshap` on PATH   | whatshap binary |
| `WORK_DIR` | `<repo>` (parent of scripts/) | base for output dirs |
| `OUT_ROOT` | `${WORK_DIR}/spechla_out`  | per-sample alignments + VCFs |
| `ASM_ROOT` | `${WORK_DIR}/asm_v2`       | typing outputs |

Full env-var list is in [../PIPELINE.md](../PIPELINE.md) §6.

---

## 5. Tuning per sample

The pipeline auto-estimates χ from the data. Defaults work for χ_R in
`[0.05, 0.50]`. For boundary cases:

| Situation | Override |
| --------- | -------- |
| χ_R < 0.03 (very deep chimerism) | `FB_MIN_AF=0.005 FB_MIN_AC=2` |
| χ_R > 0.50 (recipient majority)  | `RECIPIENT_MAJOR=1` |
| Coverage < 50× | `MASK_MIN_DEPTH=3` |
| Coverage > 500× | `MASK_MIN_DEPTH=10` |

See [../PIPELINE.md](../PIPELINE.md) §6.1 for the parameter selection guide.

---

## 6. Outputs

```
asm_v2/<SAMPLE>/
    <SAMPLE>.final_calls.tsv          ★ FINAL aggregated result (one row per gene)
    <gene_lc>/<HLA-X>/
        calls.tsv                     per-gene final 4-hap call (R/D-tagged)
        calls.baseline.tsv            baseline before EM refinement (if overridden)
        hap{1..4}.fa                  per-haplotype masked FASTA

spechla_out/<SAMPLE>/                 intermediate alignments + variants
    <SAMPLE>.merge.bam, .freebayes.vcf.gz, .pooled_continuous.vcf.gz, ...
    <SAMPLE>.chimerism.txt            χ from AD-cluster estimator
    <SAMPLE>.chi_pooled.txt           χ from pooled-continuous (per gene)
    em_refine/<gene>.{calls,summary,iterative}.tsv
```

* `<SAMPLE>.final_calls.tsv` columns:
  `sample | gene | R1 | R2 | D1 | D2 | source` (`source` ∈ {`em-refined`,
  `baseline`, `missing`}).
* Per-gene `calls.tsv` columns:
  `global_hap | assignment(R/D) | allele | em_weight`.

---

## 7. Validation

Two real samples (same patient, two time points), defaults only:

| Sample | χ_R | Field-level accuracy (24 across 6 genes) |
| ------ | --- | ---- |
| 267015 (donor major, χ_R≈0.27) | 0.27 | **24/24** |
| 267016 (deep chimerism, χ_R≈0.05) | 0.05 | **21/24** (3 R2 errors are SNR-limited) |

See [../PIPELINE.md](../PIPELINE.md) §7 for per-gene breakdown.

---

## 8. Re-running idempotently

`SKIP_DONE=1` (default) skips steps whose output already exists. To force a
re-run from a specific step, delete its output and re-invoke the driver:

```bash
rm spechla_out/mySample/mySample.freebayes.vcf.gz   # re-do variant call
rm -r asm_v2/mySample                               # re-do typing
bash scripts/polyphase_v2.sh
```
