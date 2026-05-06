# Installation

This pipeline ships as a self-contained set of scripts under `scripts/` plus
one external dependency (the **SpecHLA** reference database). Follow the
three steps below.

---

## 1. Conda environment

```bash
conda env create -f scripts/environment.yml
conda activate polyploid-hla
```

This installs python + the variant-calling toolchain (whatshap, pysam,
parasail, mappy, numpy, bowtie2, bwa, samtools, bcftools, freebayes, tabix)
with version pins matching what the pipeline was validated on.

If you cannot use conda, install the same set manually and ensure all
binaries are on `PATH`. The driver auto-discovers `python` and `whatshap`
from the active environment; override with `PYBIN=/path/to/python` or
`WHATSHAP=/path/to/whatshap` if needed.

---

## 2. SpecHLA reference database

The pipeline reuses the SpecHLA database (combined per-gene reference +
IMGT/HLA bowtie2 db + per-gene mini-refs). Get it from the upstream repo:

- https://github.com/deepomicslab/SpecHLA

```bash
git clone https://github.com/deepomicslab/SpecHLA.git
# Follow SpecHLA's own README to download db/ (HLA references + IMGT db).
```

Place (or symlink) the resulting `SpecHLA/` directory next to `scripts/`:

```
<repo-root>/
├── scripts/
└── SpecHLA/
    ├── db/
    │   ├── ref/hla.ref.extend.fa
    │   ├── ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta
    │   └── HLA/HLA_<gene>/HLA_<gene>.fa
    └── script/uniq_read_name.py, assign_reads_to_genes.py, ...
```

Or set `SPECHLA=/path/to/SpecHLA` before invoking the driver.

---

## 3. Smoke test

Verify the environment is correctly resolved without running on real data:

```bash
conda activate polyploid-hla
bash -n scripts/polyphase_v2.sh                                # syntax check
for t in python whatshap bowtie2 bwa samtools bcftools freebayes tabix; do
    command -v "$t" >/dev/null && echo "OK   $t" || echo "MISS $t"
done
python -c "import pysam, parasail, mappy, numpy, tqdm; print('python deps OK')"
ls "${SPECHLA:-./SpecHLA}/db/ref/hla.ref.extend.fa"            # database check
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
| `whatshap: error: unrecognized arguments: --ploidy` | whatshap < 1.4; upgrade to ≥ 2.2 (in `environment.yml`). |
| Driver exits with `HLA_REF not found` | SpecHLA db not at `<repo>/SpecHLA/` and `SPECHLA` env var not set. |
| `parasail`/`mappy` import error | conda env not activated, or installed via pip into wrong python. |
| Empty `calls.tsv` for a gene | usually no reads assigned at step 1; check `spechla_out/<SAMPLE>/<gene>.R1.fq.gz` size. |

---

## 5. Containerization (optional)

`environment.yml` works as the basis for a Docker / Singularity image; e.g.

```dockerfile
FROM continuumio/miniconda3
COPY scripts/environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy
SHELL ["conda", "run", "-n", "polyploid-hla", "/bin/bash", "-c"]
```

The SpecHLA database is large (~ a few GB) and is typically mounted as a
volume rather than baked into the image.
