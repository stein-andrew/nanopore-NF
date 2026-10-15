# ont-pipeline

ONT DNA (fiber-seq) and RNA (RNA004 direct RNA) processing, split so the GPU-heavy
half can run on whichever cluster has the GPUs.

| Stage | Script | Runs on | Does |
|---|---|---|---|
| 1 | `basecall.nf` | AICR (Parabricks) or Explorer | pod5 → dorado → align → restore mod tags → BAM per sample per reference |
| 2 | `analysis.nf` | Explorer | *not written yet*: DNA (Clair3, whatshap, fibertools, modkit) / RNA (modkit, gIVT, IsoQuant, StringTie2, isomatch, Clair3-RNA) |

Stage 1 writes `analysis_samplesheet.csv` with paths relative to its outdir, so
after an rsync it is Stage 2's `--input` unchanged.

## Quick start

```bash
# once per cluster: build the Stage 1 images into /scratch/$USER/ont-pipeline/containers
cd ~/ont-pipeline
sbatch -p short -c 8 --mem 32G -t 4:00:00 containers/build.sh dorado-2.1.2-dna dorado-1.3.0-rna

# Stage 1 on AICR (Parabricks) or Explorer (CPU minimap2 by default)
nextflow run basecall.nf -profile aicr,dna     --input runs.csv --outdir /scratch/$USER/ont-pipeline/results/run1
nextflow run basecall.nf -profile explorer,rna --input runs.csv --outdir /scratch/$USER/ont-pipeline/results/run1

# move to Explorer and verify
scripts/transfer.sh /scratch/$USER/ont-pipeline/results/run1 \
    stein.an@<explorer-login>:/scratch/stein.an/ont-pipeline/results/run1
```

Run Nextflow itself from a login node (Explorer) or a long, low-resource job; it
submits every task to SLURM. Requires Nextflow ≥ 25.10 (workflow outputs).

**On AICR the login node kills the JVM**, so the head job runs under SLURM:

```bash
scripts/aicr_launch.sh /scratch/$USER/ont-pipeline/runs/run1 \
    basecall.nf -profile aicr,dna --input runs.csv --outdir /scratch/$USER/ont-pipeline/results/run1
```

Dorado runs on `--gpu_queue` (default `b200-batch`; for the RTX PRO 6000 nodes use
`--gpu_queue rtx-batch --gpu_type rtx_pro_6000`). SLURM's start-time estimates on
AICR are unreliable, so don't choose a queue by them.
Parabricks runs on `--pbrun_queue` (default `b200-devel`): it is CPU-bound, so each
job gets 1 GPU and 32 CPUs, and a flowcell aligns in about 1-1.5 h (a second GPU was
no faster in benchmarks). Images and indexes live under
`/work/neu/p2026_0018_neu/ont-pipeline`, because AICR `/scratch` purges files older
than 30 days.

## Profiles

Pick one data type and one cluster.

| Profile | Sets |
|---|---|
| `dna` | dorado 2.1.2, `dna_r10.4.1_e8.2_400bps_hac@v6.0.0`, mods `5mCG_5hmCG@v1 6mA@v1`; moves (for Clair3 `_with_mv`); `lr:hq`; GRCh38 **and** HG002; keep supplementaries (`-F 260`) |
| `rna` | dorado 1.3.0, `rna004_130bps_sup@v5.2.0`, mods `inosine_m6A_2OmeA pseU_2OmeU m5C_2OmeC 2OmeG` (@v1); moves + poly-A; `splice -uf -k14`; GRCh38; primary only (`-F 2308`) |
| `explorer` | SLURM queues, reference paths, images on `/scratch/stein.an` |
| `aicr` | SLURM: GPU steps on `b200-batch` (switchable), CPU steps on `cpu`; Parabricks aligner; references under `/work/neu/p2026_0018_neu` |
| `stub` | local, no containers: `-stub-run` checks of the DAG |

Override anything in a `-params-file`, e.g. all-context 5mC instead of CpG:

```yaml
mod_models: ['5mC_5hmC@v1', '6mA@v1']   # one C model at most: dorado rejects two
```

## Samplesheet (Stage 1)

```csv
sample_id,run_id,pod5,condition,replicate
hg002_treated_rep1,PAY22766,/scratch/stein.an/HG002_data/ont/PAY22766,treated,1
```

* `sample_id` becomes the SM tag; rows sharing it are merged after alignment.
* `run_id` (one flowcell) becomes the RG ID; defaults to `sample_id`; must be unique.
* `pod5` is a file or a directory (searched recursively).
* Any other columns are carried through to `analysis_samplesheet.csv`.

## Why the tags are handled the way they are

Parabricks `--in-bam` drops MM/ML, and `--copy-comment` segfaults on large tags, so
after Parabricks `bin/copy_tags.py` grafts the carried tags (MM/ML/MN, `qs`, `mv`,
and `pt` for RNA) back from the uBAM, verbatim. The uBAM is name-sorted once per
run and shared by every reference that run was aligned to.

* **Supplementaries (DNA).** Parabricks has no `-Y`, so it hard-clips them, and
  positional tags (MM/ML/MN, `mv`) would point past the end of the clipped SEQ.
  `--unclip` puts the clipped bases back from the uBAM (H → S, reverse-complemented
  for reverse-strand records), which is what `minimap2 -Y` writes on the CPU path.
  Tested on HG002: the Parabricks output matched the CPU output record for record
  (1,176/1,176: position, CIGAR, SEQ, MM, ML, `mv`, MN), with nothing dropped.
  This replaced `modkit repair`, which aborts intermittently while finishing its
  output (modkit 0.5.0 and 0.6.4: `malloc(): unaligned fastbin chunk detected`;
  about half of repeated runs on one GM18865 input) and drops records it cannot place.
* **`mv` (DNA).** Clair3's dwell-time models
  (`r1041_e82_400bps_hac_v600_with_mv --enable_dwell_time`) read it on primary
  alignments; Clair3 skips supplementaries.
* **RNA.** Only primary alignments survive `-F 2308`, so nothing is clipped. `pt` is
  for IsoQuant. `ts` is never copied: dorado's `ts:i` would clash with minimap2's
  `ts:A`, which IsoQuant reads. (Clair3-RNA's `ont_dorado_drna004` model does not
  use `mv`.)

The CPU path (`--aligner minimap2`) carries tags through FASTQ comments
(`samtools fastq -T | minimap2 -y`) and uses `-Y` for DNA, so no transfer is needed.
Both paths give the same records.

Parabricks also ignores `--read-group-sm` and leaves records without `RG`; the
pipeline rebuilds the `@RG` (ID = run_id, SM = sample_id, plus dorado's PU/DT/DS
fields) and tags every record, as the CPU path does.

## Outputs (Stage 1)

```
<outdir>/
  analysis_samplesheet.csv   Stage 2 input; relative paths
  checksums.md5              md5sum -c from <outdir> after transfer
  bam/<ref>/<sample>.<ref>.bam(.bai)
  qc/<ref>/<sample>.<ref>.bam.flagstat
  ubam/<run_id>.ubam.bam     (params.publish_ubam)
  logs/                      pbrun logs (Parabricks path)
  pipeline_info/
```

Outputs are hard links into the work dir (`workflow.output.mode = 'link'`), so
publishing costs no space as long as outdir and workDir share a filesystem.
`.fai` and `.mmi` files are built once per reference and preset into
`params.index_dir` and reused.

## Containers

| Image | Contents |
|---|---|
| `dorado-2.1.2-dna.sif` | dorado 2.1.2 + HAC v6.0.0 and all four mod models, samtools 1.21, minimap2 2.28, modkit 0.6.4 (CPU) |
| `dorado-1.3.0-rna.sif` | dorado 1.3.0 + SUP v5.2.0 and five mod models, samtools 1.21, minimap2 2.28, modkit 0.5.0 |
| `clara-parabricks_4.7.0-1.sif` | NVIDIA NGC image (not built here) |

## Testing

```bash
nextflow lint .
nextflow run basecall.nf -profile explorer,dna,stub -stub-run --input sheet.csv --outdir stub-out
```
