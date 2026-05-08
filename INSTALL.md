# Installation

This pipeline ships as a self-contained set of scripts plus bundled HLA
reference resources under `resources/spechla/`. Follow the steps below.

---

## 1. Conda environment

```bash
conda env create -f environment.yml
conda activate polyploid-hla
```

This installs python + the variant-calling toolchain (whatshap, pysam,
parasail, mappy, numpy, bowtie2, bwa, samtools, bcftools, freebayes, tabix)
plus `blastn` for optional exon-level G group fallback diagnostics. The
`samtools` package is also used for `wgsim`, which powers the benchmark read
simulation scripts. Version pins match what the pipeline was validated on.

If you cannot use conda, install the same set manually and ensure all
binaries are on `PATH`. The driver auto-discovers `python`, `whatshap`, and
`freebayes` from the active environment; override with `PYBIN=/path/to/python`,
`WHATSHAP=/path/to/whatshap`, or `FREEBAYES=/path/to/freebayes` if needed.

---

## 2. Bundled HLA resources

The files needed from SpecHLA are already copied into the project:

| Path | Contents |
| --- | --- |
| `resources/spechla/script/` | `uniq_read_name.py`, `assign_reads_to_genes.py` |
| `resources/spechla/db/ref/` | combined HLA ref, IMGT allele FASTA, bowtie2/BWA/faidx indexes |
| `resources/spechla/db/HLA/` | six per-gene refs, G group map, exon FASTAs |

No separate SpecHLA database download is required for the default workflow.
Set `SPECHLA=/path/to/custom/resources` only when intentionally testing a custom
or refreshed database laid out in the same `db/` and `script/` structure.

If the resource directory was copied without indexes, or if you replace it with
a custom database, rebuild indexes with:

```bash
conda activate polyploid-hla
bash build_resource_indexes.sh

# or for a custom resource directory:
bash build_resource_indexes.sh --resources /path/to/custom/resources
```

The script builds the required `samtools faidx`, BWA, and bowtie2 indexes. If
`makeblastdb` is available, it also prepares BLAST databases for exon-level
diagnostics.

---

## 3. Smoke test

Verify the environment is correctly resolved without running on real data:

```bash
conda activate polyploid-hla
bash -n polyphase_v2.sh                                        # syntax check
for t in python whatshap bowtie2 bowtie2-build bwa samtools wgsim bcftools freebayes tabix makeblastdb; do
    command -v "$t" >/dev/null && echo "OK   $t" || echo "MISS $t"
done
python -c "import pysam, parasail, mappy, numpy, tqdm; print('python deps OK')"
bash build_resource_indexes.sh                                 # index check / repair
ls "${SPECHLA:-./resources/spechla}/db/ref/hla.ref.extend.fa"  # resource check
```

All lines should print `OK` / file exists. A failed `MISS` or missing
reference means the environment / database step above did not complete.

A real end-to-end test requires paired short-read FASTQs and is documented
in [README.md](README.md) §3.

---

## 4. Troubleshooting

| Symptom | Likely cause / fix |
| ------- | ------------------ |
| `freebayes: --pooled-continuous: option requires an argument` (or unknown option) | freebayes < 1.3; reinstall via `conda install -c bioconda freebayes=1.3.*`. |
| `freebayes` aborts with a C++ vector assertion | use the pinned `freebayes=1.3.6` from `environment.yml`, or set `FREEBAYES=/path/to/freebayes-1.3.6`. |
| `whatshap: error: unrecognized arguments: --ploidy` | whatshap < 1.4; upgrade to ≥ 2.2 (in `environment.yml`). |
| Driver exits with `HLA_REF not found` | Bundled resources are missing, or `SPECHLA` points to an incomplete custom resource directory. |
| `bowtie2` / `bwa` reports missing index files | run `bash build_resource_indexes.sh --resources "${SPECHLA:-resources/spechla}"`. |
| `parasail`/`mappy` import error | conda env not activated, or installed via pip into wrong python. |
| Empty `calls.tsv` for a gene | usually no reads assigned at step 1; check `spechla_out/<SAMPLE>/<gene>.R1.fq.gz` size. |

---

## 5. Containerization (optional)

`environment.yml` works as the basis for a Docker / Singularity image; e.g.

```dockerfile
FROM continuumio/miniconda3
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy
SHELL ["conda", "run", "-n", "polyploid-hla", "/bin/bash", "-c"]
```

The bundled HLA resources under `resources/spechla/` are part of the
software image unless you deliberately replace them with a mounted custom
database through `SPECHLA=/path/to/resources`.
