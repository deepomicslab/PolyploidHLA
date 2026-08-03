# Polyploid HLA Typing

End-to-end pipeline for **chimeric (k=4) HLA typing** of allo-HSCT and
solid-organ transplant samples from short-read FASTQs. The detailed aggregate
keeps the legacy R1/R2/D1/D2 assembly slots, while the primary mixed-sample
report lists four allele copies per gene sorted by estimated copy proportion
without making recipient/donor assignment part of the final call.

- Installation & environment: [INSTALL.md](INSTALL.md)

---

## 1. Files

| File | Role |
| ---- | ---- |
| `polyphase_v2.sh`           | **driver — run this** |
| `hla_polyphase_assemble.py` | baseline typing engine (whatshap polyphase + IMGT scoring) |
| `reassign_gt_chimeric.py`   | χ-aware GT correction before phasing |
| `estimate_chi_pooled.py`    | pooled-continuous χ_R estimator |
| `report_gene_abundance.py`  | writes the 15-gene mixture/copy-abundance and QC report |
| `iterative_remap_em.py`     | EM refinement (Salmon-style read remap) |
| `apply_quartet_optimization.py` | frozen class-I read-gated and class-II joint quartet optimization |
| `diagnostics/rescue_gene_binned_reads.py` | validation-only read-bin rescue diagnostic/prototype |
| `apply_class2_joint_rescue.py` | guarded class-II post-aggregate rescue |
| `em_refine_gate.py`         | per-gene EM override gate logic |
| `aggregate_calls.py`        | merges per-gene `calls.tsv` into one summary table |
| `evaluate_calls.py`         | compares `<SAMPLE>.final_calls.tsv` with `truth_typing.tsv` at 2-field and G group resolution |
| `exon_typing_from_haps.py`  | exon-level G group fallback/diagnostic for high-mask genes |
| `prepare_extra_typing_resources.py` | builds augmented references for HLA-E/F/G/H and MICA/MICB typing |
| `build_resource_indexes.sh` | rebuilds HLA resource indexes when files are missing or a custom resource set is used |
| `gene.spechla.bed`          | per-gene typing region on bundled `hla.ref.extend.fa` |
| `resources/spechla/`        | bundled SpecHLA-derived helper scripts and HLA reference files |
| `environment.yml`           | conda environment spec |
| `benchmark/` | simulation, scoring, and Slurm benchmark workflows; see [benchmark/README.md](benchmark/README.md) |
| `diagnostics/` | offline diagnostics, validation runners, and rejected alternatives kept for reference |
| `docs/` | experiment protocols and extended project documentation |
| `tests/` | focused unit and integration tests |

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
WORK_DIR=$PWD/run_outputs \
bash polyphase_v2.sh
```

`WORK_DIR` controls where outputs are written. If it is not set, the driver
uses the parent directory of `polyphase_v2.sh` and writes `${WORK_DIR}/asm_v2`
and `${WORK_DIR}/spechla_out`.

**Final result files:**

```bash
column -t run_outputs/asm_v2/mySample/mySample.copy_calls.tsv
column -t run_outputs/asm_v2/mySample/mySample.final_calls.tsv
column -t run_outputs/asm_v2/mySample/mySample.gene_abundance.tsv
```

By default, allele-copy outputs contain the six validated classical genes,
the six extended genes HLA-E/F/G/H, MICA, and MICB, plus the DRB345 linked
add-on row. The abundance report separates DRB345 into HLA-DRB3, HLA-DRB4,
and HLA-DRB5, producing 15 gene rows in total. The extended-gene calls and
abundance rows remain exploratory until independently validated.

Primary side-agnostic copy report:

```
sample    gene   copy_id copy_rank allele          allele_2field copy_fraction allele_read_count copy_read_count proportion_source legacy_slot raw_copy_fraction copy_identifiability        copy_fit_error
mySample  HLA-A  copy1   1         A*29:02:01:02   A*29:02       0.445599      1425.14           1425.14         copy_fraction_fit D2          0.445599          underdetermined_chi_regularized 0.00000000
mySample  HLA-A  copy2   2         A*01:01:01:01   A*01:01       0.344401      1149.35           1100.29         copy_fraction_fit R1          0.344401          underdetermined_chi_regularized 0.00000000
...
```

Detailed R/D-slot aggregate:

```
sample    gene   R1_full        R2_full        D1_full        D2_full        R1_report R2_report D1_report D2_report R1_read_count R2_read_count D1_read_count D2_read_count source
mySample  HLA-A  A*01:01:01:01 A*23:01:01:01 A*01:01:01:01 A*29:02:01:02 A*01:01   A*23:01   A*01:01   A*29:02   1149.35       623.77        1149.35       1425.14       em_refine
...
```

The detailed file keeps both high-resolution calls and conservative report calls:

* `*_full`: full allele chosen by the pipeline.
* `*_2field`: allele collapsed to 2-field, useful when truth is low resolution.
* `*_g_group`: allele converted through SpecHLA `hla_nom_g.txt`.
* `*_report`: equals `*_full` by default; automatically downgraded to 2-field
  when a gene has high masked sequence fraction.
* `*_fraction`: modelled haplotype proportion for the reported R1/R2/D1/D2
  call. In the standard 2+2 chimerism model, R haplotypes are `chi_R / 2` each
  and D haplotypes are `(1 - chi_R) / 2` each; EM/direct modes use the fitted
  gene-specific chi when available.
* `*_read_fraction`: allele-family read support fraction from EM read
  assignment (`tf_counts.tsv`). This is not forced to R1/R2 or D1/D2 being
  1:1 and can differ between two alleles from the same person.
* `*_read_count`: EM-assigned effective read weight for the same 2-field allele
  family.
* `*_copy_fraction_fit`: fitted copy-level fraction for each R1/R2/D1/D2 slot.
  These four values sum to 1 when read support is available. The fit minimizes
  allele-family read-support error under `x >= 0` and `sum(x)=1`, using the
  sample/gene chi as a weak regularizer when duplicated alleles make the copy
  split underdetermined.
  Values smaller than `1e-4` are written in scientific notation, so tiny values
  and exact boundary zeros are not displayed as fixed-width `0.000000`.
* `copy_fit_error`, `copy_identifiability`, `copy_chi_r`, and
  `allele_support_fraction_sum`: diagnostics for the copy-fraction fit.
  `boundary_zero` means the constrained fit placed one slot at zero because it
  lacked independent support; it is not a filled missing value.
* `mean_mask_fraction`, `report_level`, `warning`: explain why a gene was
  reported at full vs. 2-field resolution.

The pipeline also writes a concise companion file,
`<SAMPLE>.final_calls.compact.tsv`, with only the sample, gene, four reported
alleles, four fitted copy fractions, per-allele read counts, and fit diagnostics.
The original `<SAMPLE>.final_calls.tsv` remains the detailed result file.

For mixed donor/recipient samples, the primary side-agnostic report is now
written as `<SAMPLE>.copy_calls.tsv` and `<SAMPLE>.copy_calls.compact.tsv`.
These files list the four allele copies and their estimated proportions without
making R/D assignment part of the main result. The long file carries both
`allele_read_count`, the EM-assigned effective read count supporting each
reported allele family, and `copy_read_count`, the same support apportioned
across repeated copies of that allele according to the fitted copy fractions.
The compact file carries these count vectors in the same order as
`copy_fractions`. `allele_read_count` can be nonzero when a fitted duplicate
copy fraction is zero, but `copy_read_count` will be zero for that row. The
legacy R1/R2/D1/D2 slot is kept only as an annotation in the long file so old
debugging workflows still have a bridge back to the assembly slots.

When `DRB345_TYPING=1` (default), the pipeline also appends an `HLA-DRB345`
row. This is not a seventh ordinary locus: it is a DRB1-linked add-on for the
DRB3/DRB4/DRB5 genes. The add-on extracts read pairs with competitive DB
support for DRB3/4/5, EM-remaps them to a combined DRB345 allele set, and uses
the final DRB1 haplotypes to decide whether each R1/R2/D1/D2 haplotype should
carry DRB3, DRB4, DRB5, or no DRB345 gene. It does not change the classical or
extended fixed-diploid gene calls. DRB345 DB-read extraction accepts bowtie2-style alignment scores,
where the best score is often 0 and imperfect hits are negative. If the DRB1
row is high-mask / low-confidence, the add-on switches to an evidence-first
mode: long locus-unique k-mers decide which DRB3/4/5 loci have credible support,
then EM allele fractions are fit to the R/R/D/D dose model without hard DRB1
linkage. This lets DRB345 remain callable when DRB1 itself is unreliable.

The pipeline also writes `<SAMPLE>.gene_abundance.tsv`. Its fixed-diploid rows
cover the six classical genes and the six default extended genes (HLA-E/F/G/H,
MICA, and MICB), with global mixture, gene-local AF estimate, depth, residual,
and local/global agreement fields. The global mixture remains restricted to the
validated six classical genes. HLA-DRB3/4/5 are reported as three conditional
copy-abundance rows: their low/high source copy counts come from the DRB1-linked
four-haplotype result, and normal structural absence is reported explicitly.
Non-PASS pooled mixture estimates propagate as low-confidence abundance rows
rather than being presented as validated gene-level estimates.

The pipeline then writes `<SAMPLE>.cnv_loh.tsv`. A joint four-haplotype MILP
fits integer R1/R2/D1/D2 dosages in `0..3` against normalized gene depth and
allele-group absolute dosage. It reports the exact best state, a no-good-cut
second optimum, the normal-state counterfactual, and their objective gaps.
Total-copy changes are labeled `CNV`; a source-local `0/2` state at total copy
four is labeled `COPY_NEUTRAL_LOH`. The state remains visible when an allele is
shared across sources, but its confidence is `ASSIGNMENT_AMBIGUOUS` because the
source carrying that dosage cannot be identified from the mixed sample alone.
This diagnostic does not rewrite HLA allele calls.

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
| `EXTRA_TYPING_GENES` | `HLA-E HLA-F HLA-G HLA-H MICA MICB` | extended fixed-diploid loci; set to an explicit empty string to disable |
| `EXTRA_TYPING_RESOURCE_ROOT` | `${WORK_DIR}/extra_typing_resources` | temporary augmented references generated when `EXTRA_TYPING_GENES` is set |
| `GENE_ABUNDANCE_OUTPUT` | `1` | write `<SAMPLE>.gene_abundance.tsv`; set to `0` to disable |
| `CNV_LOH_OUTPUT` | `1` | write the joint integer-dosage `<SAMPLE>.cnv_loh.tsv` diagnostic; set to `0` to disable |
| `POOLED_CHI_CONTIGS` | core six reference contigs | contigs allowed to contribute to the global pooled mixture estimate |
| `EXON_TYPING` | `1` | also write exon-level fallback diagnostics (`<SAMPLE>.exon_calls.tsv`) |
| `BOWTIE2_MODE` | `very-sensitive` | bowtie2 preset for IMGT competitive mapping; use `sensitive` for faster exploratory runs |
| `BOWTIE2_K` | `30` | max alignments reported per read pair during IMGT competitive mapping |
| `ASSEMBLE_ALIGNER` | `parasail` | base-level scorer; `mappy` is faster but less exact |
| `ASSEMBLE_PREFILTER_TOP` | `200` | mappy prefilter size before parasail scoring; smaller is faster |
| `EM_REFINE_PER_GENE_CHI` | `0` | experimental; fixed pooled/global χ is the recommended default |
| `EM_REFINE_RECIPIENT_MINOR_RESCUE` | `1` | recover low-frequency recipient-only alleles when donor-major EM fitting collapses R/D to the donor-like pair |
| `QUARTET_OPTIMIZATION_PROFILE` | `normalized_joint_v1` | `normalized_joint_v1` applies the frozen optimizers; `shadow` audits without rewriting calls; `off` disables the stage |
| `REUSE_BINNING_ROOT` | empty | seed deduped FASTQs, DB BAM, per-gene FASTQs, and `header.sam` from a prior run |
| `REUSE_BINNING_CLEAN_DOWNSTREAM` | `0` | remove downstream outputs after seeding cache when intentionally recomputing calls |

The options above cover the recommended user-facing settings.

HLA-E/F/G/H and MICA/MICB are enabled by default and are included in the same
`final_calls.tsv` and `copy_calls.tsv` outputs. To run only the validated six
classical genes, explicitly disable the extended set:

```bash
EXTRA_TYPING_GENES="" \
bash polyphase_v2.sh
```

When the extended set is enabled, the driver automatically builds an augmented
reference, gene BED, and per-gene BWA references under
`EXTRA_TYPING_RESOURCE_ROOT` from the bundled IMGT-style FASTA. The global
pooled χ estimate remains restricted to the classical six contigs, so enabling
the extended set cannot change the established global estimator input. Extended
gene calls and local abundance QC should still be treated as exploratory until
validated with truth data for these loci.

Read-bin rescue is currently a validation-only diagnostic, not part of the
default production pipeline. Run `scripts/diagnostics/run_gendx_input_root_diagnostic.sh`
first to prove that strict gene binning is losing usable read evidence; only
after a rerun shows accuracy gain should rescue be promoted into the main flow.

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

The pipeline auto-estimates the lower mixture component from the data. The
pooled estimator is currently designed for mixture fractions in `[0.10, 0.50]`;
lower fractions require separate validation. For boundary cases:

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
    <SAMPLE>.copy_calls.tsv           ★ PRIMARY mixed-sample result (four sorted allele copies per gene)
    <SAMPLE>.copy_calls.compact.tsv   compact copy multiset/proportion/read-count result
    <SAMPLE>.final_calls.tsv          detailed R/D-slot aggregate result (one row per gene)
    <SAMPLE>.final_calls.compact.tsv  compact R/D-slot allele/proportion/read-count result
    <SAMPLE>.gene_abundance.tsv       15-gene mixture/copy-abundance QC report
    <SAMPLE>.cnv_loh.tsv              joint four-haplotype integer CNV/LOH states
    <SAMPLE>.exon_calls.tsv           exon-level G group diagnostic for high-mask genes
    <SAMPLE>.quartet_optimization.manifest.tsv  baseline/proposal/gate/application audit
    <gene_lc>/<HLA-X>/
        calls.tsv                     per-gene final 4-hap call (R/D-tagged)
        calls.baseline.tsv            baseline before EM refinement (if overridden)
        calls.quartet_optimization_input.tsv  pre-optimization input (if rewritten)
        hap{1..4}.fa                  per-haplotype masked FASTA

spechla_out/<SAMPLE>/                 intermediate alignments + variants
    <SAMPLE>.merge.bam, .freebayes.vcf.gz, .pooled_continuous.vcf.gz, ...
    <SAMPLE>.chimerism.txt            χ from AD-cluster estimator
    <SAMPLE>.chi_pooled.txt           χ from pooled-continuous (per gene)
    read_bin_rescue_manifest.tsv      rescue counts, retention gate, final status
    class2_joint_rescue_manifest.tsv  guarded class-II rescue audit trail
    em_refine/<gene>.{calls,summary,iterative}.tsv
    drb345/                          DRB3/4/5 linked add-on typing outputs
```

* `<SAMPLE>.final_calls.tsv` columns:
  `sample | gene | R1_full | R2_full | D1_full | D2_full | R1_2field | ... |
  R1_g_group | ... | R1_report | ... | R1_fraction | R2_fraction |
  D1_fraction | D2_fraction | R1_read_fraction | ... | R1_read_count | ... |
  R1_copy_fraction_fit | ... | copy_fit_error | copy_identifiability |
  copy_chi_r | allele_support_fraction_sum | source | mean_mask_fraction |
  report_level | warning`.
* `<SAMPLE>.final_calls.compact.tsv` columns:
  `sample | gene | R1_allele | R1_copy_fraction | R1_read_count |
  R2_allele | R2_copy_fraction | R2_read_count | D1_allele |
  D1_copy_fraction | D1_read_count | D2_allele | D2_copy_fraction |
  D2_read_count | copy_identifiability | copy_fit_error`.
* `<SAMPLE>.copy_calls.tsv` columns:
  `sample | gene | copy_id | copy_rank | allele | allele_2field |
  copy_fraction | allele_read_count | copy_read_count | proportion_source | legacy_slot |
  raw_copy_fraction | copy_identifiability | copy_fit_error`.
* `<SAMPLE>.copy_calls.compact.tsv` columns:
  `sample | gene | allele_multiset | allele_2field_multiset |
  copy_fractions | allele_read_counts | copy_read_counts | proportion_source | copy_identifiability |
  copy_fit_error`.
* `<SAMPLE>.gene_abundance.tsv` columns:
  `sample | gene | model | global_chi | global_chi_source |
  global_qc_status | global_ci95 | local_chi | n_af | median_dp | residual |
  delta_from_global | low_source_copies | high_source_copies |
  low_source_fraction | high_source_fraction | expected_gene_abundance |
  observed_low_fraction | observed_read_fraction | called_copies | status |
  reasons`.
  Fixed-diploid rows use `model=fixed_diploid`; HLA-DRB3/4/5 use
  `model=drb1_linked_conditional_copy`. `expected_gene_abundance` is normalized
  to the abundance of a diploid gene, so a standard two-copy fixed locus is
  `1.0`. DRB3/4/5 may be between `0.0` and `1.0` depending on source-specific
  presence and the mixture fraction.
* Abundance `status` values are `PASS`, `LOW_CONFIDENCE`, `MODEL_MISMATCH`,
  `NOT_ENABLED`, and `NOT_RUN`. A normal absent DRB3/4/5 locus is represented
  by `status=PASS`, `expected_gene_abundance=0`, and
  `reasons=structural_absence`. If pooled χ fails QC, the GT estimate can keep
  typing operational, but affected abundance rows remain `LOW_CONFIDENCE`.
* `<SAMPLE>.cnv_loh.tsv` includes normalized depth, four integer dosages,
  total copies, `event`, `confidence`, best/second/normal objectives, and
  `event_support`. `ASSIGNMENT_AMBIGUOUS` preserves a detected integer state
  while marking cross-source shared-allele attribution as non-identifiable.
* Per-gene `calls.tsv` columns:
  `global_hap | assignment(R/D) | allele | hap_fraction |
  allele_read_fraction | allele_read_count | em_weight` for EM-refined calls,
  or `... | hap_fraction | total_assembly_score` for baseline assembly calls.
  DRB345 add-on calls are stored under
  `asm_v2/<SAMPLE>/hla-drb345/HLA-DRB345/calls.tsv` with an extra
  `drb1_linked_locus` column.

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
