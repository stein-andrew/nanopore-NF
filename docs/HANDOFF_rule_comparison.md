# Handoff: modification filtering rules, mappability, and the Dorado 2.2.0 comparison

Repo: `~/ont-pipeline-v2` on Explorer (branch `main`, pushed to `github.com/stein-andrew/nanopore-NF`). User: `stein.an`.
This file hands over work started in a chat session. Read it fully before running anything.

## Working rules (from the user)

- Never run pandas + matplotlib on the login node (it has run out of memory before). Light pandas-only summaries are fine there.
  Anything heavier goes to Slurm: `sbatch` (preferred for anything over a few minutes) or an `srun` allocation inside `tmux`.
- Interactive `srun` must be started on its own line; wait for the compute-node prompt before sending more commands.
- Exported variables do not survive a new allocation or shell: set them again every time (block below).
- `/scratch` is purged by mtime. Durable storage is `/projects/Genometechlab`. Flag anything purge-sensitive.
- Never delete or overwrite existing results. Ask before deleting anything, before using GPUs, and before anything destructive.
- Numbers in any summary must come from the output tables, never from memory. Say plainly when evidence is weak.
- Figures: lab style, no grids, no stats boxes (style sheets in `/projects/Genometechlab/matplotlib_stylesheets/`).
- The user prefers short answers that lead with the result, then commands that can be pasted as they are.

## Environment block (paste into every new shell)

```bash
mamba activate bioinf_tools
export REPO=~/ont-pipeline-v2 B=/scratch/stein.an/basecall_cmp \
       ALL=/scratch/stein.an/basecall_cmp/pileup064_nofilter Q40=/scratch/stein.an/basecall_cmp/pileup064_nofilter_q40_v2 \
       ORTHO=/projects/Genometechlab/andrew/orthogonal \
       FA=/projects/Genometechlab/referenceFastaFiles/human/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
       MASK=/scratch/stein.an/ont-pipeline/cohort-v1/psi_loci_check/clusters.bed \
       OUT=/projects/Genometechlab/andrew/rule_comparison
```

Since 24 Nov `bioinf_tools` has `PYTHONNOUSERSITE=1` (`conda env config vars`), so an activated shell ignores
`~/.local/lib/python3.12/site-packages` and uses the env's pandas 2.3.3, matplotlib 3.9.1, pyarrow 16.1.0, pysam 0.22.1 and
pod5 0.3.15 (`~/.local` has 3.0.2, 3.10.9, 22.0.0, 0.24.0 and 0.3.44). The 19 `~/.local` packages the env lacked were
installed into it at the same versions; `currentview` (genometechlab git package) was not, since it needs pod5 >= 0.3.23.
Jobs that call `.../envs/bioinf_tools/bin/python` without activating still see `~/.local`, so sbatch command lines now
start with `PYTHONNOUSERSITE=1` (from 24 Nov; earlier sweeps ran on the `~/.local` packages). On the env's own versions
the wilson sweep (mock), cmp_ortho truth-set loading and make_synthetic_results.py -> cohort_plots.py (all 18 figures)
give the same tables as on the `~/.local` versions. The env's state is recorded in `envs/bioinf_tools.conda.yml` (`conda
env export --no-builds`) and `envs/bioinf_tools.pip.txt` (`pip freeze`), both written with `PYTHONNOUSERSITE=1` on 24 Nov.
Not pure conda (decision: keep it so): pip-installed are openpyxl, et-xmlfile and pyxlsb (before 24 Nov, for the GLORI
.xlsx/.xlsb truth tables), and on 24 Nov the 19 packages only `~/.local` had (dash, dash-bootstrap-components, kaleido,
logomaker, umap-learn, pynndescent, matplotlib-venn, nbformat, numexpr, Bottleneck, pydantic, ...) with their 19
dependencies and typing_extensions 4.16.0 (over conda's 4.12.2), so that switching off the user site removed nothing the
`~/.local` setup provided.

## Data

- Same reads, three setups (identical pod5 subsets; aligned BAMs `$B/tower_ubam/out_<run>/bam/GRCh38/<sample>.GRCh38.bam`):
  `A_d212_v6` (Dorado 2.1.2, rna004_sup@v6.0.0 + inosine_m6A_2OmeA@v2, m5C_2OmeC@v2, pseU_2OmeU@v1, 2OmeG@v1; the two v2
  models came from a pre-release models folder), `C_d212_v52` (Dorado 2.1.2, rna004_130bps_sup@v5.2.0 + four v1 models),
  `B_d130_v52` (Dorado 1.3.0, same v5.2 models).
- Samples: `GM12878`, `08_01_25_RNA004_HEK293`, `09_09_25_RNA004_HEK293`, `10_06_25_RNA004_HEK293_CP` (truth maps are HEK293(T), so
  accuracy uses the three HEK293 runs).
- No-filter pileups (modkit 0.6.4, every primary alignment, fixed thresholds 0): `$ALL/<run>/<sample>.bed`.
- No-filter pileups with MAPQ >= 40 (same, `MIN_MAPQ=40`), HEK293 only, runs A and C: `$Q40/<run>/<sample>.bed` in
  `pileup064_nofilter_q40_v2` (jobs 12094113-12094119; all six beds end at chrUn). Each run also left `<sample>.mapq.bam`
  there: keep for now (user's decision, 24 Nov), delete only when the user says so.
- The first MAPQ-40 set (`pileup064_nofilter_q40`, jobs 12085729-12085734) was corrupted and deleted on 24 Nov: an
  interactive pileup started at 13:29 on 23 Nov wrote to the same paths as the Slurm jobs and finished last, so the C
  `08_01` and `09_09` beds stopped at chr8 (73.1M and 65.9M rows instead of 164.0M and 147.9M). Do not use any output read
  from it except `wilson_q40_A_d212_v6_m6A` (read before the overwrite; byte-identical to `wilson_q40v2_A_d212_v6_m6A`).
- Truth maps: `$ORTHO` (GLORI 1.0/2.0, m6A-SAC-seq, miCLIP lifted to hg38, DART-seq lifted, BID-seq, PRAISE, SLIC-seq,
  UBS-seq, NJU-seq). Loaded, reference-checked and combined by `cmp_ortho.py`. Main set per modification: m6A
  `GLORI1.0∩2.0`, inosine `SLIC-seq_2of3`, Psi `BID∩PRAISE`, m5C `UBS-seq`; none usable for 2'-OMe.
- Depth: the deepest position in the v5.2 08/01 pileup has 41,628 reads, below the 65,535 cap, so capping never applied.

## Scripts (in `$REPO/scripts`; md5 of the versions described here)

| Script | md5 | What it does |
|---|---|---|
| `cmp_ortho.py` | c6f76c9d103958d1a6c5962a3b876399 | truth-set loading, reference check, consensus sets, keys; per-run scoring |
| `cmp_rules.py` | 9d382e3f33edb0f5a5fa5c5015e3d8e9 | six site rules (pct20, or_k20, or_k20_f5, or_k20_f10, wilson10, bg5) on no-filter pileups: precision, recall, by depth, between runs, negatives. Since 24 Nov `read_pileup(..., fail_as_canonical=True)` reads N_fail (default unchanged). The version that produced every result before that is b12d19af5918e13a734c4eb90f721e24; the chat copy listed here first (cad681c773899abd82773e4291c8352d) is the same code with a final newline |
| `wilson_sweep.py` | d677514780bead36294d6d9808287f9d | Wilson lower-bound rules vs plain cutoffs: precision at matched recall, recall by truth level, agreement by depth. History: 300ebef5579a7e5c917bac90389d6b02 produced the all-reads and MAPQ-40 sweeps (the chat copy listed here first, d1bd89ad10d27378e62f8c92c83007e2, is the same code with a final newline); 4f269ae6569cf60801eb1e47092d6de6 added the opt-in wilsonbg rules (`--background empirical` and/or `--background givt07:<table>`) and produced `wilsonbg_q40v2_*`; d6775147 adds `--truth-set` (any cmp_ortho set) and `--fail-as-canonical` (depth N_valid + N_fail and occupancy N_mod / (N_valid + N_fail) for every rule, the 20-read floor and the depth bands; for fixed-threshold pileups). Each step leaves the default tables byte-identical on the mock data sets. Mock test: `tests/mock_wilsonbg.py` (27acf031848f390b8a90568646f72673) |
| `cmp_mappability.py` | f2a86999a385d231d33a7a94e8e37711 | calls by unique-read fraction (MAPQ-40 / all valid reads) and by static BEDs (GIAB); model gap within strata |
| `cmp_pileup.sh` | 217a84f99b3f9bbfd28177af9bf22fdb | modkit pileups: `THRESH_MODE=fixed THRESHOLDS='A:0 C:0 G:0 T:0'` for no filter, optional `MIN_MAPQ` |
| `make_short_figures.py` | 5d9bc2afde193ef29962ac54be1b4c6b | figures of the short v6 vs v5.2 report (needs `make_report_figures.py`, 475daa6f72da71b1f074d84830ecd4d7) |

Check: `cd $REPO/scripts && md5sum cmp_ortho.py cmp_rules.py wilson_sweep.py cmp_mappability.py`. `cmp_rules.py`,
`cmp_mappability.py`, `wilson_sweep.py` (4f269ae6) and `tests/mock_wilsonbg.py` were committed to main on 24 Nov.

## Results so far (all reads, no MAPQ filter, means of the three HEK293 runs)

Outputs: `$OUT/<model>_<mod>/` (cmp_rules), `$OUT/wilson_<model>_m6A/` (sweep, all reads),
`$OUT/wilson_q40v2_<model>_m6A/` (sweep, MAPQ >= 40; `wilson_q40_C_d212_v52_m6A` is empty, that job failed), logs in `$OUT/logs/`.

- m6A, v5.2: no OR-gate beats pct20. Extra calls over pct20 that are GLORI sites: or_k20 about 1%, or_k20_f5 about 3%,
  or_k20_f10 about 8% (background: 0.8% of covered A are GLORI sites). bg5 (binomial vs 5-mer background) failed on real
  data: precision 0.28, recall below pct20; background varies between sites within a context, which the binomial ignores.
- Agreement between runs is not a fair criterion for count-based rules (depth repeats between runs); compare it within depth
  strata, or only between occupancy rules.
- Background per 5-mer reaches 15.6% (m6A), 17% (m5C, 12 contexts above 10%), about 12% (2'-OMe-G and -U): fixed floors are
  unsafe there.
- Inosine, v5.2: no call at >= 100 reads in 08/01 is a SLIC-seq site, under any rule (about 300 calls).
- Wilson (m6A): wilson95_t10 needs 25% at 20 reads, 20% at 50, 16% at 100, 11.9% at 1,000. At matched recall it is
  +0.9 precision points over plain cutoffs for v5.2 and +0.8 for v6, positive in all three runs. Its extra recall
  (about 1.5 points) is all in GLORI's 10-20% band (v5.2: 31% -> 37% of those sites recovered). Its added deep calls repeat
  less between runs (>= 100 reads: 0.72 vs 0.77 for pct20, v5.2). Raising the depth floor (pct20 with >= 30 reads) is far
  worse: recall 0.46. GLORI lists no site below 10%.
- v6 vs v5.2 (m6A) is unchanged by the rule: v6 about 1 point more recall, 13-14 points less precision, lower agreement.

Pre-registered decision rule for any new site rule: highest recall whose precision is at least pct20's in all three runs,
with agreement compared within depth strata. A negative control (no-filter pileup of v5.2 gIVT RNA, `--negative-bed`) is
the decisive test; it does not exist yet.

### MAPQ >= 40 Wilson sweep (task 1, done 24 Nov; `wilson_q40v2_*` against `wilson_*`, m6A, means of three runs)

Decision: no Wilson rule passes the pre-registered rule, with or without the MAPQ filter. Keep pct20, with MAPQ >= 40.

- The filter changes little: pct20 precision v5.2 0.473 -> 0.481, v6 0.338 -> 0.346; recall unchanged (0.769, 0.779);
  the v6 vs v5.2 precision gap is unchanged (13.4 -> 13.5 points).
- pct20 precision at 300+ reads recovers only partly: v5.2 0.364 -> 0.394, v6 0.262 -> 0.273 (up in all three runs), still
  14 (v5.2) and 11 (v6) points below the 100-299 band. Calls at 300+ drop about 15% (v5.2 394 -> 334 per run).
- The Wilson gain at matched recall survives and grows a little: wilson95_t10 v5.2 +0.90 -> +1.25 points (worst run +0.81),
  v6 +0.83 -> +1.12 (worst run +0.81). Its extra recall is still in GLORI's 10-20% band (v5.2 30.5% -> 36.5%).
- But wilson95_t10 against pct20 at MAPQ 40 is +1.5 recall points and -1.2 precision points (both models, precision lower
  in all three runs); no Wilson rule in the sweep has higher recall with precision >= pct20 in all three runs.
- Wilson's deep calls still repeat less: >= 100 reads, v5.2 pct20 0.766 -> 0.768, wilson95_t10 0.719 -> 0.721; v6 pct20
  0.702 -> 0.700, wilson95_t10 0.669 -> 0.675. Multimappers do not explain the gap.

### Wilson over a background (wilsonbg, 24 Nov; MAPQ >= 40, m6A, `$OUT/wilsonbg_q40v2_<model>_m6A/`)

Rules: lower bound >= background + d (`wilsonbg<conf>_d<d>`), and >= max(10%, background + d) (`..._t10_d<d>`; identical to
`_d<d>` at d = 10 by construction). Backgrounds: empirical 5-mer from each run (both models), and the GM12878 gIVT 9-mer
table at the 0.70 threshold (v5.2 only; a lower bound on no-filter background, since the table counts only calls with
probability >= 0.70 as modified, over reference-matched 9-mers). Decision (user, 24 Nov): no hybrid rule; pct20 with
MAPQ >= 40 stays the default. givt07 is promising; follow-ups (v5.2 only): score it against maps not used to choose it,
and test it on 0.70-threshold pileups, where the background definition matches the table.

- Empirical 5-mer background: no rule passes the pre-registered rule in either model, and all 24 empirical wilsonbg rules
  lose precision at matched recall in every run (means -1.9 to -7.8 points for v6, -3.5 to -14.1 for v5.2). Median background, v5.2: DRACH 1.12-1.42% vs other 0.83-1.08% (true
  sites lift DRACH); v6: DRACH 1.41-1.56% vs other 1.92-2.30%.
- gIVT 0.70 background (v5.2): median 0.35% (DRACH 0.21%, other 0.37%), max 10.2%; the table covers 100% of sites (every
  m6A 9-mer has support >= 218). `wilsonbg_givt07_90_d10` passes the precision/recall part: precision +2.8 to +3.6 points
  and recall +3.0 to +3.6 points over pct20 in each run; +8.6 points at matched recall; GLORI 10-20% recall 0.41 vs 0.30.
  `wilsonbg_givt07_95_d10`: precision +7.1 to +8.2, recall +0.2 to +1.0. But the gain is all below 50 reads: precision at
  100-299 reads 0.44 vs 0.53 and at 300+ 0.24 vs 0.39 (pct20), and calls at >= 100 reads repeat less (0.729 vs 0.768,
  lower in every pair; below 100 reads 0.694 vs 0.667).

## Open tasks, in order

1. ~~Read the MAPQ-40 Wilson sweep~~: done, see the results above.
2. **Mappability layer.** Get GIAB GRCh38 v3.6 stratifications (low mappability `GRCh38_lowmappabilityall.bed.gz`,
   segmental duplications `GRCh38_segdups.bed.gz`; check whether the user already has them from the Fiber-seq work) into
   `/projects/Genometechlab/andrew/references/giab_strat_v3.6/`. Run `cmp_mappability.py` for m6A, Inosine, Psi, m5C with
   `--model v6=A_d212_v6 --model v5.2=C_d212_v52`, the three HEK293 samples, `--all-dir $ALL --q40-dir $Q40`. Read
   `gaps.tsv` (v6 minus v5.2 within each stratum) and the truth-density column: short-read maps are blind in
   low-mappability regions, so compare models within strata, not strata against each other.
3. **Rerun `cmp_rules.py` on the MAPQ-40 pileups** (all eight modifications, both models) so every rule result matches the
   pipeline's new default (`mods_unique_only = true`).
4. **Dorado 2.2.0 comparison** (released 23 Nov 2026). It adds `rna004_sup@v6.0.0_inosine_m6A_2OmeA@v2`, `_m5C_2OmeC@v2` and
   `_pseU_2OmeU@v2` (2OmeG stays at v1). Steps:
   - Download `https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.2.0-linux-x64.tar.gz` and the models
     (`dorado download --model <name> --models-directory <durable dir>`; no GPU needed): `rna004_sup@v6.0.0`, the three v2
     models, `rna004_sup@v6.0.0_2OmeG@v1`, and `rna004_130bps_sup@v5.2.0` with its four v1 models.
   - Model identity: `diff -rq` the released v2 models against the pre-release folder run A used (find it from run A's
     BAM `@PG` lines or the tower bundle config). Report the result: the short report's methods depend on it.
   - Plan two new basecalls of the same pod5 subsets: `D` = Dorado 2.2.0 + released v6 set (v2 x3, 2OmeG v1); `E` =
     Dorado 2.2.0 + v5.2 + v1 models. D vs E isolates the models on the new Dorado; E vs C checks Dorado 2.2.0 against
     2.1.2; D vs A shows what changed since the report (pseU v2, Dorado). Dorado 2.x needs CUDA >= 12.8. Earlier runs used
     the tower (4 A100), about 35-45 min per sample per setup. **Ask the user before any GPU job**, and confirm where.
   - Then no-filter MAPQ-40 pileups (`cmp_pileup.sh`) and the same scoring; the scorers use run labels such as
     `A_d212_v6`, so check how each script names and pairs runs before adding new labels.
5. **Pipeline (after the user decides):** report the Wilson lower bound as a column in the site tables, and optionally a
   `site_rule` parameter (`pct` default, `wilson` opt-in with target and confidence). Find where the 20% rule is applied
   (`conf/rna.config`, `bin/mod_stages.sh`, `bin/*.py`) and propose a patch before editing.
6. **Housekeeping:** copy the no-filter pileups and `clusters.bed` to `/projects` if the work continues more than a couple
   of weeks; commit `params/Ver005_wbRNA.yaml` and its samplesheet (Ver005 basecalling: Nextflow head job `nf-Ver005`,
   run dir `/scratch/stein.an/ont-pipeline/runs/Ver005_wbRNA`, outputs `/scratch/stein.an/ont-pipeline/results/11_18_26_RNA004_Ver005_wbRNA`;
   its 437 GB of pod5 under `/scratch/stein.an/11_18_26_RNA004_Ver005_wbRNA` need a durable copy).

## Related state

- Short report (4 pages, no filter, 2'-OMe naming) is on the user's Mac in `~/Downloads/v6_vs_v52_short`. Before it goes to
  ONT, the model-identity check in task 4 decides one methods sentence.
- modkit in the pipeline's `mods` image is 0.5.0 (no `--high-depth`); the comparison used 0.6.4. An upgrade to 0.6.4 is
  planned together with a switch to no-filter pileups; it needs `modkit_max_depth = 65535` plus `--high-depth` in
  `bin/validated_pileup.sh`, and a one-sample check of 0.5.0 vs 0.6.4 no-filter output first.
