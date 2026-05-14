# Polyploid HLA Typing

End-to-end pipeline for **chimeric (k=4) HLA typing** of allo-HSCT and
solid-organ transplant samples from short-read FASTQs. Outputs 4 haplotype
sequences per gene tagged `R`(ecipient) / `D`(onor).

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
| `rescue_gene_binned_reads.py` | conservative read-bin rescue v2 before per-gene alignment |
| `apply_class2_joint_rescue.py` | guarded class-II post-aggregate rescue |
| `em_refine_gate.py`         | per-gene EM override gate logic |
| `aggregate_calls.py`        | merges per-gene `calls.tsv` into one summary table |
| `evaluate_calls.py`         | compares `<SAMPLE>.final_calls.tsv` with `truth_typing.tsv` at 2-field and G group resolution |
| `exon_typing_from_haps.py`  | exon-level G group fallback/diagnostic for high-mask genes |
| `build_resource_indexes.sh` | rebuilds HLA resource indexes when files are missing or a custom resource set is used |
| `gene.spechla.bed`          | per-gene typing region on bundled `hla.ref.extend.fa` |
| `resources/spechla/`        | bundled SpecHLA-derived helper scripts and HLA reference files |
| `environment.yml`           | conda environment spec |
| `octopus_to_imgt.py`, `caller_free_4hap.py` | rejected alternatives, kept for reference |

---

## 2. Install

See [INSTALL.md](INSTALL.md). Short version:

```bash
conda env create -f environment.yml
conda activate polyploid-hla
# HLA reference resources needed by the pipeline are bundled under
# resources/spechla. Set SPECHLA=/path/to/custom/resources only
# when intentionally overriding them.

# Optional: repair/rebuild bundled resource indexes.
bash build_resource_indexes.sh
```

---

## 3. Quick start

From the repository root:

```bash
FQ1=/path/sample.R1.fq.gz \
FQ2=/path/sample.R2.fq.gz \
SAMPLE=mySample \
RECIPIENT_MAJOR=0 \
bash polyphase_v2.sh
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

The file keeps both high-resolution calls and conservative report calls:

* `*_full`: full allele chosen by the pipeline.
* `*_2field`: allele collapsed to 2-field, useful when truth is low resolution.
* `*_g_group`: allele converted through SpecHLA `hla_nom_g.txt`.
* `*_report`: equals `*_full` by default; automatically downgraded to 2-field
  when a gene has high masked sequence fraction.
* `*_fraction`: modelled haplotype proportion for the reported R1/R2/D1/D2
  call. In the standard 2+2 chimerism model, R haplotypes are `chi_R / 2` each
  and D haplotypes are `(1 - chi_R) / 2` each; EM/direct modes use the fitted
  gene-specific chi when available.
* `mean_mask_fraction`, `report_level`, `warning`: explain why a gene was
  reported at full vs. 2-field resolution.

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
| `SPECHLA`  | `resources/spechla` | bundled HLA resource root; override only for a custom database |
| `PYBIN`    | first `python` on PATH     | python binary |
| `WHATSHAP` | first `whatshap` on PATH   | whatshap binary |
| `FREEBAYES` | first `freebayes` on PATH | freebayes binary; use 1.3.6 if newer builds abort |
| `THREADS`  | `8` | threads for bowtie2, BWA, whatshap, EM remap, and samtools helper steps |
| `SAMTOOLS_THREADS` | `$THREADS` | threads for samtools view/sort/index/merge |
| `WORK_DIR` | parent of this repository | base for output dirs |
| `OUT_ROOT` | `${WORK_DIR}/spechla_out`  | per-sample alignments + VCFs |
| `ASM_ROOT` | `${WORK_DIR}/asm_v2`       | typing outputs |
| `EXON_TYPING` | `1` | also write exon-level fallback diagnostics (`<SAMPLE>.exon_calls.tsv`) |
| `BOWTIE2_MODE` | `very-sensitive` | bowtie2 preset for IMGT competitive mapping; use `sensitive` for faster exploratory runs |
| `BOWTIE2_K` | `30` | max alignments reported per read pair during IMGT competitive mapping |
| `ASSEMBLE_ALIGNER` | `parasail` | base-level scorer; `mappy` is faster but less exact |
| `ASSEMBLE_PREFILTER_TOP` | `200` | mappy prefilter size before parasail scoring; smaller is faster |
| `EM_REFINE_PER_GENE_CHI` | `0` | experimental; fixed pooled/global χ is the recommended default |
| `EM_REFINE_RECIPIENT_MINOR_RESCUE` | `1` | recover low-frequency recipient-only alleles when donor-major EM fitting collapses R/D to the donor-like pair |
| `READ_BIN_RESCUE` | `1` | run truth-free read-bin rescue before per-gene alignment |
| `READ_BIN_RESCUE_GENES` | `HLA-DPB1` | default target for conservative rescue; broaden only for validation |
| `READ_BIN_RESCUE_RETENTION_GATE` | `1` | require abnormal full-vs-retained read support loss before rewriting a bin |
| `REUSE_BINNING_ROOT` | empty | seed deduped FASTQs, DB BAM, per-gene FASTQs, and `header.sam` from a prior run |
| `REUSE_BINNING_CLEAN_DOWNSTREAM` | `0` | remove downstream outputs after seeding cache when intentionally recomputing calls |

The options above cover the recommended user-facing settings.

The default read-bin rescue is deliberately conservative: DPB1-only,
paired-mate evidence, and retention-gated. Broad all-gene rescue addresses a
real binning-loss failure mode but can increase phasing cost and should be run
only as an explicit validation mode.

If indexes are missing after copying or replacing the resource directory, run:

```bash
bash build_resource_indexes.sh --resources "${SPECHLA:-resources/spechla}"
```

For exploratory reruns where speed matters more than final reporting, a useful
starting point is:

```bash
THREADS=16 BOWTIE2_MODE=sensitive BOWTIE2_K=15 ASSEMBLE_PREFILTER_TOP=100 \
bash polyphase_v2.sh
```

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

For most samples, keep the defaults and only adjust the overrides above when
the sample falls into one of those boundary cases.

---

## 6. Outputs

```
asm_v2/<SAMPLE>/
    <SAMPLE>.final_calls.tsv          ★ FINAL aggregated result (one row per gene)
    <SAMPLE>.exon_calls.tsv           exon-level G group diagnostic for high-mask genes
    <gene_lc>/<HLA-X>/
        calls.tsv                     per-gene final 4-hap call (R/D-tagged)
        calls.baseline.tsv            baseline before EM refinement (if overridden)
        hap{1..4}.fa                  per-haplotype masked FASTA

spechla_out/<SAMPLE>/                 intermediate alignments + variants
    <SAMPLE>.merge.bam, .freebayes.vcf.gz, .pooled_continuous.vcf.gz, ...
    <SAMPLE>.chimerism.txt            χ from AD-cluster estimator
    <SAMPLE>.chi_pooled.txt           χ from pooled-continuous (per gene)
    read_bin_rescue_manifest.tsv      rescue counts, retention gate, final status
    class2_joint_rescue_manifest.tsv  guarded class-II rescue audit trail
    em_refine/<gene>.{calls,summary,iterative}.tsv
```

* `<SAMPLE>.final_calls.tsv` columns:
  `sample | gene | R1_full | R2_full | D1_full | D2_full | R1_2field | ... |
  R1_g_group | ... | R1_report | ... | R1_fraction | R2_fraction |
  D1_fraction | D2_fraction | source | mean_mask_fraction | report_level |
  warning`.
* Per-gene `calls.tsv` columns:
  `global_hap | assignment(R/D) | allele | hap_fraction | em_weight` for
  EM-refined calls, or `... | hap_fraction | total_assembly_score` for baseline
  assembly calls.

If truth is available, evaluate with:

```bash
python evaluate_calls.py \
  --truth truth/truth_typing.tsv \
  --calls asm_v2/mySample/mySample.final_calls.tsv
```

Evaluation reports only `2field` and `g_group` accuracy. It intentionally does
not score 3-field because many truth entries are 2-field or G group resolution.
For G group scoring, truth alleles that cannot be uniquely mapped through
`hla_nom_g.txt` remain at 2-field resolution instead of being treated as
false mismatches.

---

## 7. Re-running idempotently

`SKIP_DONE=1` (default) skips steps whose output already exists. To force a
re-run from a specific step, delete its output and re-invoke the driver:

```bash
rm spechla_out/mySample/mySample.freebayes.vcf.gz   # re-do variant call
rm -r asm_v2/mySample                               # re-do typing
bash polyphase_v2.sh
```

For expensive real-data replays, reuse the dedup / competitive DB map /
per-gene binning outputs and recompute only downstream steps:

```bash
REUSE_BINNING_ROOT=/path/to/previous/spechla_out \
REUSE_BINNING_CLEAN_DOWNSTREAM=1 \
SKIP_DONE=1 \
bash polyphase_v2.sh
```

When reporting accuracy, regenerate evaluation from the current
`final_calls.tsv`; stale `*.eval.txt` files can describe an older call set.
