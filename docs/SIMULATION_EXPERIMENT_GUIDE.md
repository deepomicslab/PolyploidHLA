# PolyploidHLA Simulation Experiment Guide

## 1. Objective and scope

This benchmark evaluates three claims about PolyploidHLA:

1. **Accuracy:** the reported four-copy HLA genotype and copy proportions agree with predefined simulation truth.
2. **Breadth:** performance is characterized across loci, graft fractions, sequencing depths, and identifiable genotype structures relevant to solid-organ transplantation.
3. **Robustness:** common sequencing perturbations do not cause an unacceptable loss of performance.

The target application is mixed DNA from solid-organ transplant recipients, including lung, liver, and kidney transplantation. No organ-specific graft-fraction or sequencing-depth distributions are currently available. Therefore, this study uses one common solid-organ transplant model and must not claim organ-specific performance. Lung, liver, and kidney applicability will require independent real cohorts.

The recipient is the major source and the graft is the minor source. All production runs use:

```bash
RECIPIENT_MAJOR=1
```

Maternal-fetal mixtures are outside the scope of this benchmark.

## 2. Code and data locations

Keep source code and generated benchmark data separate.

```text
/data2/wangxuedong/polyploid_benchmark_v2/
  PolyploidHLA/                 # project source code; tracked by Git
  PolyploidHLA_simulation/      # generated simulation data; project sibling
```

Use the following fixed paths:

```bash
PROJECT_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA
BENCH_ROOT=/data2/wangxuedong/polyploid_benchmark_v2/PolyploidHLA_simulation
mkdir -p "$BENCH_ROOT"/{config,truth,allele_fasta,gene_reads,reads,runs,metrics,summaries,figures,logs}
```

All generated data must be written under `BENCH_ROOT`, never under the source repository:

```text
PolyploidHLA_simulation/
  config/              # frozen design manifests and software versions
  truth/               # allele and proportion truth; unavailable to the caller
  allele_fasta/        # selected truth sequences
  gene_reads/          # temporary per-gene/per-copy wgsim reads
  reads/               # merged individual-level R1/R2 FASTQs
  runs/                # PolyploidHLA outputs
  metrics/             # sample-by-locus results
  summaries/           # confidence intervals and summary tables
  figures/             # manuscript figures
  logs/                # simulation and pipeline logs
```

### 2.1 Software environment

Use the existing Conda environment defined by the project `environment.yml`:

```bash
conda activate polyploid-hla
export PATH="$CONDA_PREFIX/bin:$PATH"
```

Putting `$CONDA_PREFIX/bin` first is required on servers where user-installed `samtools` or `bcftools` directories precede Conda on `PATH`. Verify the actual executables, not only the Conda package records:

```bash
command -v python bowtie2 bwa samtools wgsim bcftools freebayes whatshap
python --version
bowtie2 --version | head -n 2
samtools --version | head -n 2
bcftools --version | head -n 2
freebayes --version
whatshap --version
```

The currently installed environment contains Python 3.10.14, Bowtie2 2.5.4, BWA 0.7.18, Samtools/wgsim 1.19.2, Bcftools 1.19, FreeBayes 1.3.6, and WhatsHap 2.8. Create the environment only if `conda env list` shows that `polyploid-hla` is absent:

```bash
cd "$PROJECT_ROOT"
conda env create -f environment.yml
```

Do not create a second environment for the benchmark when the existing `polyploid-hla` environment passes the executable checks above.

### 2.2 Executable smoke run

The smoke launcher automatically re-enters the existing `polyploid-hla` environment, puts its `bin` directory first, generates one six-locus dataset, and runs the production caller:

```bash
cd "$PROJECT_ROOT"
bash benchmark/slurm/run_simulation_smoke.sh
```

To generate and QC the merged FASTQs without running the caller:

```bash
RUN_CALLER=0 bash benchmark/slurm/run_simulation_smoke.sh
```

The underlying driver is `benchmark/scripts/simulate_wgsim_benchmark.py`. For example, a manifest-only preview of the main accuracy design is:

```bash
python benchmark/scripts/simulate_wgsim_benchmark.py \
  --experiment accuracy_main \
  --scenario distinct4 \
  --individuals 80 \
  --graft-fractions 0.10 0.20 0.30 0.40 0.50 \
  --coverages 300 \
  --dry-run
```

Remove `--dry-run` to generate FASTQs. Add `--run-caller` only when compute allocation and wall time are appropriate for all requested samples.

### 2.3 Formal one-day submission

Submit the prespecified Experiment 1 matrix as 20 concurrent Slurm array tasks:

```bash
cd "$PROJECT_ROOT"
bash benchmark/slurm/submit_formal_simulation.sh
```

Each task covers four consecutive individuals, runs all five graft fractions at 300x, and therefore performs 20 caller runs. The 20 tasks cover `SIM0001-SIM0080` without overlap and write to independent shards (`accuracy_main_v4_shard01` through `accuracy_main_v4_shard20`), preventing concurrent manifest writes. Together they produce 400 caller runs. The v4 benchmark generates reads directly from the production caller allele FASTA and uses the same file as the callable-allele database. It sets `CLASS2_DPB1_RARE_COLLAPSE=0` while retaining the other truth-free class-II rescue rules. Each task requests 4 CPUs and 16 GB for 36 hours; the expected elapsed time is 12-20 hours under the measured node load. Do not submit the script again unless the existing job has been cancelled or completed and a rerun is intentional.

Monitor the array tasks and inspect their scheduler logs with:

```bash
squeue -j <job_id>
tail -f "$BENCH_ROOT"/logs/slurm/formal_main_v4_<job_id>_*.out
```

## 3. Locus scope

Report loci in separate validation tiers:

| Tier | Loci | Role |
| --- | --- | --- |
| Primary | HLA-A, HLA-B, HLA-C, HLA-DRB1, HLA-DQB1, HLA-DPB1 | Prespecified primary endpoints |
| DRB345 add-on | HLA-DRB3, HLA-DRB4, HLA-DRB5 | Evaluate locus presence and allele accuracy using valid DRB1 linkage |
| Exploratory extensions | HLA-E, HLA-F, HLA-G, HLA-H, MICA, MICB | Report separately; do not pool with primary loci |

Do not report one pooled "15-locus accuracy." Promote an exploratory locus to validated status only after it passes the same prespecified QC and accuracy criteria as the primary loci.

## 4. Mixture model and minimum graft fraction

For graft fraction $f$, the four physical copy proportions are:

$$
\boldsymbol{p}=\left(\frac{1-f}{2},\frac{1-f}{2},\frac{f}{2},\frac{f}{2}\right),
$$

where the first two copies are recipient alleles and the last two are graft alleles.

Equal depth within an individual is a hard simulation constraint:

$$
p_{R1}=p_{R2}=\frac{1-f}{2},\qquad
p_{G1}=p_{G2}=\frac{f}{2}.
$$

The main study does not simulate within-person allele imbalance or copy-number abnormalities. If one source is homozygous, its two physical copies still have equal proportions, but their observable support is combined for proportion scoring.

### 4.1 Minimum fraction consistent with production code

The production defaults include:

```text
FB_MIN_AF=0.03
GT_DROP_FP_AF=0.05
```

`reassign_gt_chimeric.py` drops variants with observed AF below `GT_DROP_FP_AF`. The effective minimum tested single-copy AF is therefore 5%. Because each graft allele has expected fraction $f/2$:

$$
\frac{f_{\min}}{2}=0.05 \quad\Rightarrow\quad f_{\min}=0.10.
$$

The main benchmark consequently uses graft fractions at or above 10%. Keep `FB_MIN_AF=0.03` and `GT_DROP_FP_AF=0.05` fixed for every condition. Do not tune them using known simulation truth.

Testing graft fractions below 10% requires a separately frozen, truth-independent low-frequency mode and is not part of this benchmark.

### 4.2 Experimental class-I distinct-allele rescue

The production caller contains a default-off rescue for HLA-A, HLA-B, and HLA-C. It applies only when the EM quartet contains exactly three distinct 2-field alleles and replaces one duplicated copy with an unselected top-four EM candidate. The candidate must have EM fraction at least 0.005, EM weight at least 20, and a fourth-to-fifth candidate support ratio of at least 1.5. The replacement position is chosen by the lowest dosage residual, and the updated residual and rescue details are written to the per-gene summary.

Enable the frozen experimental profile with:

```bash
EM_REFINE_CLASS_I_DISTINCT_RESCUE=1
```

The rescue remains disabled by default until it is evaluated with independent seeds in `distinct4`, `shared1`, and `homozygous_graft` scenarios and on real validation samples. The existing v4 `distinct4` runs were used to derive the gate and are development data: replay improvements on those runs must not be reported as independent validation performance.

Validation must be paired on identical FASTQs. The rescue arm and baseline arm differ only in `EM_REFINE_CLASS_I_DISTINCT_RESCUE`; the baseline arm reuses the FASTQs generated for the rescue arm. Submit the three-scenario smoke validation with:

```bash
cd "$PROJECT_ROOT"
bash benchmark/slurm/submit_class_i_rescue_validation.sh smoke
```

The launcher submits a three-task rescue array and a dependent three-task baseline array. After both arrays complete, summarize final calls with:

```bash
python benchmark/scripts/summarize_class_i_rescue_validation.py \
  --bench-root "$BENCH_ROOT" \
  --experiment-glob 'class_i_rescue_smoke_v1_*' \
  --output "$BENCH_ROOT/metrics/class_i_rescue_smoke_v1/paired_locus_metrics.tsv"
```

Do not submit the formal matrix if either arm has missing outputs, or if `shared1` or `homozygous_graft` shows a rescue-associated regression or increased oversplitting in the smoke run. If the smoke run passes those plumbing and safety checks, submit the independent formal matrix with `bash benchmark/slurm/submit_class_i_rescue_validation.sh formal`. This creates 60 rescue tasks and a dependent 60-task paired baseline stage, covering 80 individuals per scenario at graft fractions 0.10, 0.20, and 0.30 with master seed `20260729`. Score it with experiment glob `class_i_rescue_validation_v1_*` and prespecify the final acceptance margins before inspecting formal truth-based results.

The paired smoke run completed on 2026-07-29 with all six caller jobs successful and no missing class-I outputs. The rescue improved `distinct4` HLA-C from 3/4 to 4/4, did not change `homozygous_graft`, but failed the prespecified `shared1` safety criterion. In `shared1` HLA-A it replaced a correct duplicated `A*24:29` copy with false `A*24:02`, reducing recall from 4/4 to 3/4 and increasing the predicted distinct-allele count from three to four. It also increased distinct-allele count in `shared1` HLA-C without improving copy recall. Across the nine class-I smoke sample-loci, baseline and rescue both recovered 27/36 copies: one locus improved and one regressed. Therefore, the formal 60+60 validation matrix was not submitted, the rescue remains disabled by default, and the current distinct-quartet gate must not be promoted to production.

### 4.3 Offline private-read replay

The diagnostic-only private-read replay is implemented in
`offline_class_i_private_rescue.py`. The frozen development profile enables the
four-distinct path, scans the top eight EM candidates, and retains the default
private-read thresholds of 30 incoming pairs, 10 weak-singleton pairs, and a
3-fold incoming-to-weak support ratio. Run the complete v4 replay with:

```bash
python offline_class_i_private_rescue.py \
  --bench-root ../PolyploidHLA_simulation \
  --experiment-glob 'accuracy_main_v4_shard*' \
  --enable-four-distinct \
  --four-distinct-top-n 8 \
  --support-cache ../PolyploidHLA_simulation/metrics/class_i_private_support_cache_v7 \
  --output ../PolyploidHLA_simulation/metrics/class_i_private_v7_optimized_full.tsv
```

The cache key includes the baseline quartet, incoming candidates, k-mer
parameters, and the paths, sizes, and modification times of both FASTQs and the
IMGT FASTA. Stale inputs therefore cannot silently reuse support. Cache files
are written atomically and a stopped run can resume from completed loci. On
2026-07-29, a full cache-hit replay took 1.35 seconds and was byte-identical to
the frozen v7 per-locus output: 4500/4800 class-I copies (93.75%), with 233
improved loci and no development-set regressions. This is a development replay,
not independent validation, and the rule remains outside the production caller.

### 4.4 Class-II post-rescue safety scope

The v4 per-locus audit showed that broad post-aggregation class-II rescue
reduced DPB1 from 1312/1600 to 1150/1600 correct copies and DQB1 from 1496/1600
to 1487/1600, while DRB1 remained unchanged at 1278/1600. The production
`truthfree_readsupport_class2` profile therefore defaults to
`CLASS2_RESCUE_GENES=HLA-DRB1`. This retains the non-regressing DRB1 branch and
leaves DQB1 and DPB1 at their pre-rescue calls.

With this scope, class-II recovery is 4086/4800 (85.125%), versus 3915/4800
(81.5625%) under the former broad profile. Across the original six-gene v4
calls, this changes recovery from 8182/9600 (85.229%) to 8353/9600 (87.010%).
These values use the original class-I calls; combining a separately evaluated
class-I replay with this class-II result must be reported as an offline composite,
not as one end-to-end caller run.

The historical broad eligibility can be reproduced explicitly with:

```bash
CLASS2_RESCUE_GENES="HLA-DRB1 HLA-DQB1 HLA-DPB1" bash polyphase_v2.sh ...
```

Do not re-enable DQB1 or DPB1 by default without independent, gene-specific
evidence and a zero-regression gate on held-out simulations.

### 4.5 DRB1 private-read candidate

The optimized offline private-read replay also accepts an explicit gene list.
The first DRB1 experiment reused the frozen class-I v7 thresholds without
DRB1-specific tuning: EM top eight, at least 30 candidate-private read pairs,
at most 10 pairs for the replaced singleton, and a candidate-to-replaced ratio
of at least 3. Reproduce it with:

```bash
python offline_class_i_private_rescue.py \
  --bench-root ../PolyploidHLA_simulation \
  --experiment-glob 'accuracy_main_v4_shard*' \
  --genes HLA-DRB1 \
  --enable-four-distinct \
  --four-distinct-top-n 8 \
  --support-cache ../PolyploidHLA_simulation/metrics/drb1_private_support_cache_v1 \
  --output ../PolyploidHLA_simulation/metrics/drb1_private_v1_full.tsv
```

On the complete v4 development set, DRB1 increased from 1278/1600 (79.875%)
to 1321/1600 (82.5625%): 43 loci improved, none regressed, and eight additional
accepted replacements were score-neutral. Net gains by graft fraction from
10% through 50% were +7, +9, +5, +5, and +17 copies. Combined with protected
DQB1 and DPB1 pre-rescue calls, this gives a development-only class-II result
of 4129/4800 (86.021%). Combining that result with the separate class-I v8a
replay gives 8640/9600 (90.000%), but this is an offline composite rather than
one end-to-end or independently validated result.

The DRB1 decision is truth-blind, but the reported gain was measured on the
same development simulation used to select this candidate. It is now part of
the integrated default profile, but publication claims still require a
prespecified independent simulation matrix with no meaningful regression in
shared-allele and homozygous-graft strata.

### 4.6 Main-pipeline integration and batch result

The main pipeline now runs `apply_private_read_rescue.py` after class-II rescue
and before final DRB3/4/5 typing. Its default frozen parameters are:

- Genes: HLA-A, HLA-B, HLA-C, and HLA-DRB1.
- EM candidate rank: top eight.
- Candidate-private support: at least 30 read pairs.
- Replaced singleton support: at most 10 read pairs.
- Candidate-to-replaced support ratio: at least 3.
- EM-gap override: enabled for HLA-A/B/C only; DRB1 uses the four-distinct path.
- Strict second pass: enabled for HLA-A/B/C only, after a successful first pass.
- Second-pass candidate support: at least 50 private read pairs.
- Second-pass replaced singleton: at most 5 private read pairs.
- Second-pass support ratio: at least 5; the first removed allele is blocked.

This is the `class_i_v8b + drb1_v1 + class2_safe` profile used by the optimized
development replay. The replay recovered 8663/9600 copies (90.240%) and 1850/2400
exact quartets (77.083%), versus 8182/9600 and 1452/2400 at baseline. These are
development replay results, not independent end-to-end validation results.

Set `PRIVATE_READ_RESCUE=0` for a complete ablation, or set
`PRIVATE_READ_RESCUE_SECOND_PASS_GENES=""` to disable only v8b's second pass.
Every run writes a truth-free
per-sample audit to `<ASM_ROOT>/<sample>/<sample>.private_read_rescue.tsv` and
backs up a changed gene call table as `calls.pre_private_rescue.tsv`.

After copy-level reporting, `update_batch_results.py` writes one row per sample
and gene to `BATCH_RESULTS_FILE`. File locking and atomic replacement allow
concurrent Slurm tasks to share the file; rerunning a sample replaces its
`experiment + condition + sample + gene` rows rather than duplicating them.
Simulation columns include scenario, graft fraction, coverage, read length,
insert mean and SD, error rate, and master seed.

The formal v4 launcher uses:

```text
PolyploidHLA_simulation/results/accuracy_main_v4.all_samples.tsv
```

The current backfill contains 2,800 rows from 400 completed sample-condition
runs, seven reported genes per run, and no duplicate keys. It records the
existing end-to-end calls and has SHA-256
`cfa4502d422baca99de7baadbb5d4eeb8e82451a9ab387cc997303047d038e65`.
It must not be described as a private-rescue rerun: future executions of the
integrated main pipeline will update the corresponding rows in this file.

## 5. Reproducible allele selection

### 5.1 Candidate filtering

For each gene:

1. Extract allele records from the frozen bundled IMGT-style FASTA.
2. Retain records with a parseable allele name, usable sequence over the typing region, and no more than 1% ambiguous bases.
3. Collapse full-resolution records into 2-field families.
4. Sample a 2-field family first, then select one full-resolution sequence from that family. This prevents families with many full-resolution records from being oversampled.
5. Save the candidate pool, exclusions, selected records, and `allele_seed` in the truth manifest.
6. Reuse the same individual-level allele combinations across graft fractions, depths, and paired perturbations.

The balanced sampling set tests allele diversity; it does not represent population frequency. Population-generalization claims require a separate frequency-weighted real-data study.

### 5.2 Genotype sets

Construct 80 complete simulated individuals for each scenario. Each individual contains a genotype at every target gene.

| Scenario | Per-gene construction | Purpose |
| --- | --- | --- |
| `distinct4` | Four different 2-field families | Primary accuracy and depth experiments |
| `shared1` | Three families; one is shared by recipient and graft | Shared-allele identifiability |
| `homozygous_graft` | One graft family duplicated; two different recipient families | Graft dosage and homozygosity |
| `near_neighbor` | One pair of different families in the lowest nonzero 5% sequence-distance stratum, plus two other families | Sequence ambiguity stress test |

For `near_neighbor`, compute distance over the jointly covered typing region, ignore unmatched terminal sequence, and require at least one observable difference. If a locus cannot provide 80 unique valid combinations, use all available combinations and report the actual sample size. Never duplicate combinations to reach a nominal count.

DRB3/4/5 must be selected using a frozen DRB1 linkage table rather than independent random sampling.

## 6. Read simulation with wgsim

Use the `wgsim` executable supplied by the pinned environment. Record its absolute path and version.

### 6.1 Coverage definition

Let $C$ be total locus coverage across all four copies, $p_i$ the copy proportion, $L_i$ the allele sequence length, and $r$ the read length. Generate:

$$
N_i=\operatorname{round}\left(\frac{C p_i L_i}{2r}\right)
$$

read pairs for copy $i$.

The two alleles from one person must have equal expected base coverage, not necessarily equal read-pair counts. Different sequence lengths require different pair counts. Verify before simulation:

$$
D_i=\frac{2rN_i}{L_i}\approx Cp_i.
$$

Within-source depth differences may only reflect integer rounding. Post-alignment depth imbalance is recorded as mapping or coverage bias and must not alter truth proportions.

### 6.2 Standard wgsim settings

The standard condition is PE150 with a 350-bp mean outer distance and 50-bp standard deviation. Use a fixed per-copy seed derived from the individual, gene, source, and haplotype identifiers.

```bash
wgsim \
  -N "${READ_PAIRS}" \
  -1 150 -2 150 \
  -d 350 -s 50 \
  -e 0.001 \
  -r 0 -R 0 -X 0 \
  -S "${COPY_SEED}" \
  truth_copy.fa copy.R1.fastq copy.R2.fastq
```

Setting `-r 0` prevents wgsim from introducing untracked germline variants into an allele sequence whose identity is already the truth. The sequencing-error rate is controlled by `-e`.

Rewrite read names after simulation so that production FASTQs do not expose gene, source, haplotype, or allele truth. Preserve the original-to-blinded name map under `truth/`; the caller must not access it.

### 6.3 Merge reads by simulated individual

Generate temporary reads independently for each gene and physical copy, then merge all target-gene reads for one individual and one condition into one paired FASTQ dataset:

```text
reads/<experiment>/<condition>/<sample>.R1.fastq.gz
reads/<experiment>/<condition>/<sample>.R2.fastq.gz
```

After read-name blinding, merge and shuffle as follows:

1. Read each R1/R2 input as synchronized four-line FASTQ records and fail if names or record counts do not match.
2. Concatenate the paired records from every simulated gene and physical copy for the individual.
3. Apply one seeded permutation to the pair indices; never shuffle R1 and R2 independently.
4. Write the permuted mates to the final R1/R2 paths with gzip compression.
5. Verify matching final record counts and mate names, then remove temporary gene-level FASTQs only after the final files pass QC.

The final FASTQs must preserve pairing but contain no truth labels. Each merged pair is processed once through the normal production entry point, preserving competitive gene binning and pooled chimerism estimation. Record the merge seed and the ordered list of input files in the sample manifest.

Example run:

```bash
FQ1="$BENCH_ROOT/reads/<experiment>/<condition>/<sample>.R1.fastq.gz" \
FQ2="$BENCH_ROOT/reads/<experiment>/<condition>/<sample>.R2.fastq.gz" \
SAMPLE=<sample> \
RECIPIENT_MAJOR=1 \
WORK_DIR="$BENCH_ROOT/runs/<experiment>/<condition>" \
bash "$PROJECT_ROOT/polyphase_v2.sh"
```

## 7. Depth groups

Use four total per-locus depths for every gene:

| Total depth | Interpretation |
| ---: | --- |
| $50\times$ | Low-coverage boundary |
| $100\times$ | Gives each graft copy about $5\times$ at $f=10\%$ |
| $300\times$ | Main condition; gives each graft copy about $15\times$ at $f=10\%$ |
| $1000\times$ | High-depth upper-bound condition |

Expected coverage per graft copy is $Cf/2$:

| Total depth | $f=10\%$ | $f=20\%$ | $f=30\%$ |
| ---: | ---: | ---: | ---: |
| $50\times$ | $2.5\times$ | $5\times$ | $7.5\times$ |
| $100\times$ | $5\times$ | $10\times$ | $15\times$ |
| $300\times$ | $15\times$ | $30\times$ | $45\times$ |
| $1000\times$ | $50\times$ | $100\times$ | $150\times$ |

Use remapped observed depth for analysis, while retaining planned depth in the manifest.

## 8. Core experiments

### 8.1 Experiment 1: absolute typing and fraction accuracy

- Scenario: `distinct4`
- Total depth: $300\times$
- Graft fractions: 0.10, 0.20, 0.30, 0.40, 0.50
- Replicates: 80 complete individuals reused across fractions
- Runs: $5\times80=400$

The primary typing endpoint is 2-field allele-set recall across the four-copy multiset for the six primary loci. Each sample-locus receives 0, 1, 2, 3, or 4 correct alleles, and the aggregate accuracy is the number of correctly recovered truth copies divided by the number of truth copies. Exact-quartet accuracy is a stringent secondary endpoint. Absolute graft-fraction error is the primary proportion endpoint; recovery of both graft alleles and no-call rate are secondary endpoints.

At $f=0.50$, source labels are exchangeable; score the unordered four-copy multiset and copy proportions, not recipient-versus-graft assignment.

### 8.2 Experiment 2: fraction-depth operating region

- Scenario: `distinct4`
- Graft fractions: 0.10, 0.20, 0.30
- Total depths: $50\times$, $100\times$, $300\times$, $1000\times$
- Replicates: 80 complete individuals reused across all 12 cells
- Runs: $3\times4\times80=960$

For each cell, report recovery of at least one and both graft alleles with two-sided Clopper-Pearson 95% confidence intervals. If all 80 individuals succeed, the two-sided exact lower confidence bound is approximately 95.5%. Do not fit or claim continuous LOD95 without denser fraction points and independent validation.

The non-300x cells were launched on 2026-07-30 with the integrated
`class_i_v8b + drb1_v1 + class2_safe` main-pipeline profile:

- Slurm array: `3151794` (`1-60%12`).
- Dependent scoring job: `3151795` (`afterok:3151794`).
- Depths: 50x, 100x, and 1000x; 20 shards per depth.
- Fractions: 0.10, 0.20, and 0.30; four individuals per shard.
- Total new caller runs: 720.
- Master seed: `20260728`; sample IDs and genotypes are paired with v4.
- Batch result: `PolyploidHLA_simulation/results/accuracy_depth_v1.all_samples.tsv`.
- Score output: `PolyploidHLA_simulation/metrics/accuracy_depth_v1/allele_set_summary.tsv`.

The depth launcher sets `--cleanup-caller-intermediates`. Cleanup occurs only
after a successful caller exit and removes regenerable input FASTQs, BAMs, and
augmented-reference indexes. It preserves truth/design manifests, caller logs,
ASM calls, final/copy result tables, EM evidence TSVs, rescue manifests, and the
shared batch result. `samples.tsv` records the number of cleaned files in
`caller_intermediates_cleaned`; all deleted inputs can be regenerated from the
recorded design and seeds.

### 8.3 Experiment 3: genotype identifiability

- Scenarios: `distinct4`, `shared1`, `homozygous_graft`, `near_neighbor`
- Graft fractions: 0.10, 0.20, 0.30
- Total depth: $300\times$
- Replicates: 80 complete individuals per scenario, reused across fractions

Reuse matching `distinct4` FASTQs from Experiment 1. Classify failures as allele-set error, copy-dose error, source-assignment error, or unidentifiable copy split. Merge proportions for indistinguishable repeated alleles before scoring.

### 8.4 Experiment 4: common solid-organ transplant matrix

Because no lung-, liver-, or kidney-specific graft-fraction and depth distributions are available, run one common technical matrix:

- Scenario: `distinct4`
- Graft fractions: 0.10, 0.20, 0.30, 0.40
- Total depths: $100\times$, $300\times$, $1000\times$
- Replicates: 80 complete individuals per cell

Do not duplicate these FASTQs under lung, liver, and kidney labels. Describe the result as a solid-organ transplant simulation. Organ-specific claims require real organ-specific data or defensible empirical distributions.

### 8.5 Experiment 5: technical robustness

Use 80 complete individuals at $f=0.10$ and $C=300\times$. Reuse the standard condition from Experiment 1 and generate four paired perturbations with the same allele truth and deterministic seeds:

| Condition | wgsim or post-processing setting |
| --- | --- |
| Standard | PE150, `-d 350 -s 50 -e 0.001` |
| Short reads | PE100, with pair counts recalculated to preserve total base coverage |
| High error | PE150, `-e 0.01` |
| PCR duplicates | Replace 20% of unique fragments with duplicated fragments while keeping total reads fixed |
| Background | Add a prespecified non-HLA background model while keeping HLA depth fixed |

Compare each individual with its standard-condition result. Prespecify non-inferiority margins before the full run; provisional margins are no more than a 5-percentage-point accuracy loss and no more than a 0.02 increase in graft-fraction absolute error.

## 9. Truth manifests and blinding

The individual manifest should contain:

```text
sample_id experiment condition scenario graft_fraction total_coverage read_length insert_mean insert_sd error_rate genotype_seed read_seed imgt_truth_version imgt_run_version recipient_major
```

The copy-level truth table should contain:

```text
sample_id gene source haplotype allele_full allele_2field g_group sequence_id sequence_length expected_fraction read_pairs copy_seed
```

Write truth only under `$BENCH_ROOT/truth`. PolyploidHLA commands, environment variables, and working directories must not contain truth paths. Score results only after calls have been frozen.

## 10. Evaluation

### 10.1 Typing

Use `<SAMPLE>.copy_calls.tsv` as the primary side-agnostic output. Normalize truth and predictions to 2-field resolution and compare them as multisets that preserve duplicate alleles.

For locus $g$:

$$
\mathrm{Recall}_g=\frac{|T_g\cap P_g|_{\mathrm{multiset}}}{4}.
$$

Retain the complete per-locus hit distribution rather than reducing partial calls to failures:

$$
H_g=|T_g\cap P_g|_{\mathrm{multiset}}\in\{0,1,2,3,4\}.
$$

The main aggregate typing accuracy is $\sum_g H_g/(4N)$, where $N$ is the number of evaluated sample-loci. Report the counts and proportions of 0/4 through 4/4 hits. Exact-quartet accuracy, $\Pr(H_g=4)$, remains a secondary all-or-none summary.

Report:

- allele-set recall and the 0/4 through 4/4 hit distribution;
- exact-quartet accuracy as a stringent secondary endpoint;
- recovery of both graft alleles;
- no-call rate;
- 2-field primary results and G-group/full-resolution supplementary results.

Program failures, missing samples, missing locus rows, no-calls, and partial calls remain in the denominator.

### 10.2 Proportions

Match predicted and truth alleles before comparing proportions. For shared or homozygous alleles, combine proportions that cannot be uniquely assigned. Report graft-fraction absolute error and calibration as primary proportion metrics; report copy-level mean absolute error as supplementary.

### 10.3 Uncertainty

Bootstrap complete simulated individuals, not loci, because all loci from one individual share the same merged FASTQ and pipeline run. Use exact binomial confidence intervals for condition-specific detection rates.

## 11. Quality control

Before running PolyploidHLA, require:

- `wgsim` command, version, and seed logged for every physical copy;
- read-pair counts matching the manifest;
- equal `expected_fraction` for the two alleles from one person;
- length-normalized expected depths equal within each person, up to pair-count rounding;
- merged R1 and R2 containing the same number of records with preserved pairing;
- no source, gene, haplotype, or allele labels in production read names;
- no truth paths in the caller environment;
- fixed `FB_MIN_AF=0.03`, `GT_DROP_FP_AF=0.05`, and `RECIPIENT_MAJOR=1`.

After running PolyploidHLA, require:

- a recorded exit status and resource usage;
- all six primary locus rows present or explicitly scored as no-call;
- nonnegative copy fractions summing to $1\pm10^{-4}$ when reported;
- estimated chimerism in $[0,1]$;
- per-locus assigned-read count, observed depth, mask fraction, and warning fields retained.

## 12. Figures and supported conclusions

The simulation study should produce five main figures:

1. Workflow, truth generation, merged-FASTQ design, and evaluation.
2. Truth-based typing and graft-fraction accuracy.
3. Fraction-depth operating-region heat maps with exact confidence intervals.
4. Locus and genotype-scenario performance, including primary, DRB345, and exploratory tiers.
5. Paired robustness changes under the four perturbations.

Supported claims must remain within the tested model:

- accurate four-copy typing and fraction estimation under closed-set simulation;
- an empirically defined operating region for graft fraction and depth;
- performance across tested loci and genotype structures;
- robustness within the tested wgsim perturbation range;
- technical relevance to solid-organ transplant mixtures, without organ-specific clinical claims.

Simulation alone cannot establish clinical validity. Independent lung, liver, and kidney transplant cohorts remain necessary for organ-specific claims.

## 13. Execution stages

1. Run five individuals per cell as a smoke/pilot study.
2. Verify manifests, wgsim generation, merged FASTQs, blinding, pipeline execution, and scoring.
3. Freeze code, thresholds, endpoints, exclusions, and seeds.
4. Expand each required cell to 80 complete individuals.
5. Generate all tables and figures from one version-controlled aggregation workflow.

The v2 repository now contains the production caller, manifest generator, `wgsim` driver, merged-FASTQ builder, Experiment 1 launcher, and frozen side-agnostic allele-set scorer. A frozen statistical summarizer is still required before manuscript tables and figures are generated. The remaining genotype-scenario and robustness generators must be implemented and validated before Experiments 3 and 5 are launched.