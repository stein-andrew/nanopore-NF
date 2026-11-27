# Handoff: modification filtering rules, mappability, and the Dorado 2.2.0 comparison

Repo: `~/ont-pipeline-v2` on Explorer (branch `main`, pushed to `github.com/stein-andrew/nanopore-NF`). User: `stein.an`.
This file hands over work started in a chat session. Read it fully before running anything.

## Working rules (from the user)

- Never run pandas + matplotlib on the login node (it has run out of memory before). Light pandas-only summaries are fine there.
  Anything heavier goes to Slurm: `sbatch` (preferred for anything over a few minutes) or an `srun` allocation inside `tmux`.
- Interactive `srun` must be started on its own line; wait for the compute-node prompt before sending more commands.
- Exported variables do not survive a new allocation or shell: set them again every time (block below).
- `/scratch` is purged by mtime. Durable storage is `/projects/Genometechlab`. Flag anything purge-sensitive.
- Standing rule (24 Nov): submit CPU-only sbatch jobs and run analyses without asking, as long as they write only to new
  folders. Still ask before: any GPU job; deleting, moving or overwriting existing files; pushing to GitHub or merging
  into main; installing or changing packages; editing pipeline code on main. Never delete or overwrite existing results.
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
| `wilson_sweep.py` | bdfe601bfa2d9bff7df42ed0e4f67f4c | Wilson lower-bound rules vs plain cutoffs: precision at matched recall, recall by truth level, agreement by depth. History: 300ebef5579a7e5c917bac90389d6b02 produced the all-reads and MAPQ-40 sweeps (the chat copy listed here first, d1bd89ad10d27378e62f8c92c83007e2, is the same code with a final newline); 4f269ae6569cf60801eb1e47092d6de6 added the opt-in wilsonbg rules (`--background empirical` and/or `--background givt07:<table>`) and produced `wilsonbg_q40v2_*`; d6775147 adds `--truth-set` (any cmp_ortho set) and `--fail-as-canonical` (depth N_valid + N_fail and occupancy N_mod / (N_valid + N_fail) for every rule, the 20-read floor and the depth bands; for fixed-threshold pileups); 821efb0a adds `fpcorr_pct20` (with givt07: the pipeline's cov20_fpcorr_pct20, percent modified - 100 x the gIVT rate >= 20, uncorrectable sites never called, as bin/fp_correct.py + mod_stages.sh); bdfe601b adds `--fpcorr-threshold` (gIVT column for fpcorr_pct20, the threshold rounded down to 0.01, 0.70 at least) and d = 20 to the wilsonbg offsets. Each step leaves the default tables byte-identical on the mock data sets. Mock test: `tests/mock_wilsonbg.py` (a73ccc91ebebf4c1e50a1fbdbec06801) |
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

### gIVT follow-ups (24 Nov; v5.2, MAPQ >= 40, m6A; GLORI and the independent maps m6A-SAC-seq, miCLIP-CITS, miCLIP-CIMS)

Outputs: `$OUT/wilsonbg_q40v2_C_d212_v52_m6A_vs_<map>` and `..._fpc_vs_<map>` (no-filter pileups);
`$OUT/wilsonbg_f070q40_C_d212_v52_m6A_vs_<map>` and `..._fpc_vs_<map>` (fixed-0.70 pileups `$B/pileup064_fixed070_q40`,
`--fail-as-canonical`; jobs 12111072-74). `m6A_2+methods` counts GLORI as one of its methods, so it is not independent.

- Independent maps: the givt07 rules beat pct20 in every run on every map, by less in points because pct20's precision
  against these maps is 0.07-0.11. `givt07_90_d10` passes the pre-registered test on all three; `givt07_95_d10` on both
  miCLIP maps but not SAC-seq (recall -0.04, -0.23, +0.13). The deep-band weakness replicates on every map.
- Against the pipeline's correction: `fpcorr_pct20` already gains most of what gIVT adds over pct20 (GLORI, no filter:
  precision +6.9 to +7.9, recall -1.4; fixed 0.70: +3.3 to +4.2, about -1.4). `givt07_95_d10` matches or beats fpcorr's
  precision with 1.4-2.4 points more recall in every run on GLORI and miCLIP-CIMS (CITS: one run -0.05 precision);
  SAC-seq is a wash. The Wilson 95% gIVT curve at fpcorr's recall: GLORI +2.4 to +3.1 points (no filter), +0.9 to +1.9
  (fixed 0.70).
- The >= 100-read weakness stays when the background matches (fixed 0.70, fails canonical; GLORI, 100-299 / 300+ reads):
  pct20 0.840 / 0.809, fpcorr 0.873 / 0.848, givt07_95_d10 0.802 / 0.689 (262 calls at 300+ vs 156 for pct20); repeat at
  >= 100 reads 0.883, 0.892, 0.862. Same on every independent map.
- Occupancy definition: on the fixed-0.70 pileups pct20's GLORI precision is 0.754 / 0.780 / 0.717 (no filter: 0.468 /
  0.531 / 0.445) for 1.7-2.1 points less recall; calls in 08_01 drop from 19,735 to 11,968 with the same 12,049 GLORI
  sites covered.
- modkit omits nothing at 0.70: the >= 20-read universe equals the no-filter one, and N_valid + N_fail at 0.70 equals the
  no-filter call count at every truth site; 505-525 GLORI sites per run reach 20 reads only through N_fail.

### Occupancy definitions (24 Nov; v5.2, MAPQ >= 40, m6A; GLORI and the three independent maps)

Four definitions, each scored with pct_20, fpcorr_pct20 (gIVT column matching the threshold: 0.70 for 1-2, 0.81 for 3-4)
and wilsonbg_givt07_95_d20 (lower bound >= p0 + 0.20, pre-specified): (1) no filter; (2) fixed 0.70, fails canonical;
(3) pipeline thresholds A:0.8105 C:0.834 G:0.9512 T:0.9219, fails dropped (what the pipeline uses now); (4) pipeline
thresholds, fails canonical. Pileups `$B/pileup064_pipethr_q40` (jobs 12124039-41); sweeps `$OUT/wilsonbg_pthrq40_{dropfail,
failcanon}_C_d212_v52_m6A_vs_<map>`, `$OUT/wilsonbg_{q40v2,f070q40}_C_d212_v52_m6A_d20_vs_<map>`. Definition 3 here uses
modkit 0.6.4; on chr1 of 08_01 it matches the pipeline's own modkit 0.5.0 path (N_mod identical, N_valid/N_fail at 99.1%/98.2%
of sites, pct_20 calls differ at 4 of 171,606 sites).

- The definition matters more than the rule. GLORI, pct_20 (calls, true calls, precision range over runs): (1) 21,560,
  10,202, 0.445-0.531; (2) 13,366, 9,953, 0.717-0.780; (3) 12,403, 9,506, 0.743-0.795; (4) 12,156, 9,823, 0.786-0.828.
- (4) beats the current (3): pct_20 precision +3.3 to +4.7 points in each run and 317 more true GLORI calls (3% more);
  dropping fails leaves 6% fewer sites with >= 20 valid reads (12,423 vs 13,233 GLORI sites). Same order on every
  independent map (true calls (4) vs (3): SAC 1,835 vs 1,743, CITS 1,554 vs 1,500, CIMS 1,979 vs 1,938).
- fpcorr_pct20 adds little once calls are thresholded: gain at matched recall on GLORI (1) +5.2, (2) +2.6, (3) +1.1,
  (4) +0.8 points; on the independent maps under (3) and (4) -0.01 to +0.23.
- wilsonbg_givt07_95_d20 is strict: GLORI precision 0.746-0.899, recall 0.50-0.56, about 30% fewer true calls than
  pct_20; gain at matched recall GLORI +0.8 to +2.9, CITS +0.8 to +1.3, CIMS +1.0 to +1.8, SAC -0.24 to +0.10. At >= 100
  reads it is the most precise of the three under every definition (GLORI (4): 100-299 / 300+ 0.937 / 0.924 vs pct_20
  0.910 / 0.910); repeat at >= 100 reads (4): 0.939 vs 0.914.
- Every Wilson rule on definitions 2 and 4 (`$OUT/summaries/wilson_rules_def2_def4.tsv`: wilson90/95/99 t5-t20 and the 27
  gIVT rules; the empirical wilsonbg rules were not computed for these two definitions): none has precision >= definition-4
  pct_20 and more true calls in every run on every map. Definition 2 fails everywhere; on definition 4, givt07_99_d7.5,
  givt07_95_d10 and wilson95_t10 pass on GLORI and both miCLIP maps but not m6A-SAC-seq (ties within a few true calls,
  or slightly lower precision), and none keeps pct_20's precision at >= 100 reads. Largest gain at matched recall: the
  givt07 d = 5 rules (definition 4, GLORI +7.5), at 15-17 points less precision than pct_20.

**Conclusion (user, 24 Nov):** no Wilson rule beats definition-4 pct_20 on every map, so the m6A rule is definition 4
(pipeline thresholds, failed calls counted as canonical) + pct_20 (+ fpcorr_pct20, as now). Wilson rules are reporting
columns only.

### Definition 3 vs 4 beyond m6A (24 Nov; v5.2, MAPQ >= 40, pipeline thresholds)

Sweeps `$OUT/wilsonbg_pthrq40_{dropfail,failcanon}_C_d212_v52_{Inosine,Psi,m5C}` (fpcorr column 0.81 for inosine, 0.92
for Psi, 0.83 for m5C; table `$OUT/summaries/def3_def4_m5C.txt`, script `def3_def4_report.py`); calls and repeat for
the 2'-OMe codes in `$OUT/summaries/{calls,repeat}_def3_def4_other_mods.tsv`.

**Decision (user, 26 Nov): definition 4 for all modifications.** Validation status by code: m6A (four maps), inosine
(SLIC-seq) and Psi (BID∩PRAISE, weak: about 20 true calls) were scored; **m5C and the 2'-OMe codes (Am, Cm, Gm, Um) are
not validated by any map**: m5C's UBS-seq gives too few true calls (13-25 per run, precision about 0.5% under both
definitions) to validate either, and the 2'-OMe codes have no map. Definition 4 is applied to them by extension only.

- Inosine (SLIC-seq_2of3; 447 sites covered under 3, 569 under 4): pct_20 under 4 makes 789 calls per run, 121 true,
  precision 0.111-0.185, against 1,337, 110 and 0.053-0.108 under 3. No call at >= 100 reads is a SLIC-seq site under
  either definition (precision 0 at 100-299 and 300+), as before.
- Psi (BID∩PRAISE; 117 / 129 sites covered): under 4, 141 calls, 20 true, precision 0.138-0.150, repeat at >= 100
  reads 0.872; under 3, 355, 21, 0.047-0.073, 0.786. Weak evidence: about 20 true calls.
- m5C (UBS-seq; 686 / 757 of its 2,181 sites covered): pct_20 under 4 makes 2,921 calls per run, 16 true, precision
  0.005-0.006, recall 0.019-0.022; under 3, 5,296, 25, 0.004-0.006, 0.033-0.041. fpcorr_pct20: 1,930, 13, 0.006-0.008
  vs 4,164, 21, 0.005-0.006. Neither definition finds UBS-seq sites (precision about 0.5% under both); very weak
  evidence, 13-25 true calls.
- m5C, Am, Cm, Um, Gm (calls only): definition 4 cuts pct_20 calls per run to m5C 2,921 (from 5,296), Am 28 (88),
  Cm 204 (469), Um 264 (929), Gm 51 (171). Repeat below 100 reads is lower under 4 for every one of these codes (next
  point); at >= 100 reads it is mixed and rests on few calls under 4 (Am 2, Cm 9, Gm 5, Um 19 per run; m5C 167).
- Between-run agreement below 100 reads (pct_20 repeat, mean of the three run pairs, definition 4 vs 3): lower under 4
  for every modification except m6A. Inosine 0.444 vs 0.522, Psi 0.426 vs 0.472, m5C 0.406 vs 0.557, Am 0.224 vs 0.322,
  Cm 0.252 vs 0.395, Gm 0.357 vs 0.416, Um 0.254 vs 0.454; m6A goes up, 0.848 vs 0.832 (`between_runs.tsv` of the
  pthrq40 sweeps, the same for every map; 2'-OMe codes from `repeat_def3_def4_other_mods.tsv`). The 2'-OMe codes (Am, Cm,
  Gm, Um) have no map to check precision on; m5C's (UBS-seq, above) gives too few true calls to tell the definitions apart.

The switch, branch `mods-fail-canonical`, merged to main 26 Nov: `mods_fail_canonical` (default true) rewrites the
pileup counts before staging (`mod_stages.sh failcanon`: Nvalid_cov += Nfail, percent = Nmod / Nvalid_cov; raw.bed keeps
modkit's counts) and counts failed calls as canonical in MOD_FEATURES (`mod_feature_counts.sh --fail-canonical`);
`mods_wilson_columns` (default false) adds column 20 (95% Wilson lower bound) and 21 (givt_d20, using the threshold-matched
gIVT column) to the gIVT stages and `percent_wilson_lb95`, `givt_d20` to feature_mods.tsv.gz. Same branch: the boolean
params are string-safe (`--mods_unique_only false` from the command line used to count as true), and
`conf/explorer.config` restricts CPU tasks to AVX2 nodes (`params.slurm_cpu_constraint`; a run's broken nodes go in
`params.slurm_exclude`, not `process.clusterOptions`, which would replace the constraint).

Real runs on test-chr22 (3 samples, fixed pipeline thresholds; `/scratch/stein.an/ont-pipeline/test-chr22-failcanon/`,
`check_failcanon.txt`): main vs branch with both flags false: every output file identical except logs and the recorded
sample-probs estimate (those tasks reran, see below; the thresholds are fixed). Both flags true: the coverage stage equals
`cmp_rules.read_pileup(fail_as_canonical=True)` on raw.bed after the base check, 0 mismatches in Nvalid, Nmod and percent
for every code and sample; MOD_FEATURES pileup agreement 100% in all nine sample x run pairs; Wilson bound and givt_d20
match an independent computation. Definition 4 vs 3 on chr22: cov20 rows x1.06-1.07, cov20_fpcorr_pct20 m6A x0.92, Ino
x0.68, m5C x0.31, the 2'-OMe codes and Psi x0.10-0.20. Two bugs the real runs found, fixed before the merge:
feature_mods.py rejected the 21-column stage (867d103), and **modkit pileup puts a failed call in Nfail whatever the read's
base at the site** (a passing call on another base is Ndiff), so definition 4 (Nvalid + Nfail, here and in the rule
comparison) includes failed calls on reads with another base, about 0.5% of depth at chr22 candidate sites; MOD_FEATURES
now counts them the same way (8b32f32).

Resume caveat (Nextflow 26.04): a task's hash includes the bin scripts it calls, by path, size and date, so a run from a
worktree reruns every bin-using task, and in the main checkout any task whose bin script was edited reruns.

Cohort rerun, prepared 26 Nov, not launched: `cohort-v1/params_v3.yaml` (params_v2 + `mods_fail_canonical: true`, outdir
`results-v3`), from `cohort-v1` with `-resume 489ff548-3100-46cd-867a-566c928e2741`. Expected to rerun MOD_STAGES and
MOD_FEATURES (15 each; script and bin changes, including the mods-filter-mismatch base check) and COHORT_MATRICES; the rest
cached (all 200 work folders of the v2 run present on 26 Nov).

### MOD_STAGES reference-base check (24 Nov; cohort-v2 results, `cohort-v1/results-v2`)

Rows whose reference base on their strand is not the code's base: cov20 0.15-0.54% per code; cov20_pct20 0.2-0.7% (m6A)
to 7-25% (Ino, Am, Cm, Psi, Um, Gm); cov20_fpcorr_pct20 0 (gIVT cannot correct them). They reach cohort_plots.py only
through fig 3's table, fig 12's not-correctable share (they are 99.8-99.98% of uncorrectable rows) and fig 13's
cross-talk and background (0.24-0.56% of target rows). Merged to main 24 Nov (edf95ca, from branch `mods-filter-mismatch`):
fp_correct.py --filter_mismatch in MOD_STAGES (uncorrected stages rebuilt from the checked rows; raw.bed untouched) and
upper-cased base comparison; final stage unchanged (checked on ST1932 chr22); tests: compare_fpcorr.py --run, the polars
port harness (its one DIFFERENT, polya_features on real chr22, is in code this change does not touch), a stub run (200
processes, 0 failed). The cohort is not rerun for it; it goes in with the next cohort rerun.

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
   - Done 24 Nov: image `containers/dorado-2.2.0-rna.def` built (job 12124759; `%test` passed: dorado 2.2.0+645d4e4,
     samtools 1.21, minimap2 2.28-r1209, modkit 0.5.0, 12 model folders). Copies:
     `/scratch/stein.an/ont-pipeline/containers/dorado-2.2.0-rna.sif` and
     `/projects/Genometechlab/andrew/containers/dorado-2.2.0-rna.sif`, sha256
     84c42df869d7a5af7e77d33d835182f0f1bd5da3804c317ad44c5e97084c2351 (both; also in `dorado-2.2.0-rna.sif.sha256`).
   - Done 24 Nov, model identity (run A used `RNA004_SUP_v6.0.0_v2MODS/` copies: tower bundle
     `/data/stein_cmp/tower_bundle/models/...`, Explorer copies `/projects/Genometechlab/andrew/custom_models/` and
     `/scratch/stein.an/tower_bundle/models/`, identical to each other), the image's `/opt/models/SHA256SUMS` file by file:
     `inosine_m6A_2OmeA@v2` released = pre-release (24 of 24 files identical); `m5C_2OmeC@v2` released != pre-release
     (21 of 24 files differ: `config.toml` only in `creation_date`, 22 Sep 2026 pre-release vs 3 Nov 2026 released;
     20 of 23 tensors have other values at the same sizes, i.e. retrained weights, same architecture). `pseU_2OmeU`:
     run A used v1, the release has v2. For the short report: run A's m6A, inosine and 2'-OMe-A results are the
     released model's; its m5C and 2'-OMe-C results come from a pre-release m5C model that the release replaced.
   - 26 Nov, runs D, E, Cx (approved; all on one GPU type, A100 80GB, `$B/explorer_ubam/`): the tower's driver is
     550.90.07 (CUDA 12.4), too old for Dorado 2.2.0 (bundles cudart 12.8.57); Explorer d1028/d1029 have 570.86.15 (CUDA
     12.8), as do the H200s. The image's v5.2.0 set is byte-identical to run C's models and its v6.0.0 base and 2OmeG@v1 to
     run A's, so E vs C differs only in Dorado and GPU. Cx = Dorado 2.1.2 + run C's models on the same A100s, 09_09 only
     (GPU effect). One benchmark file for all three, `bench_explorer_a100.csv` (Dorado 2.2.0, both base models; job
     12178278 on d1029, NVIDIA A100-SXM4-80GB, driver 570.86.15). Basecalls: array 12178280 (D, E: multigpu, 4 GPUs, two
     at a time) and Cx via `submit_cx_when_free.sh` (the gpu QOS allows 1 GPU per job); each task records node, GPU model,
     driver and dorado's batch sizes in `<run>/<sample>.gpus.txt`. Then `chain_align_pg.sh`: `align_runs.sh` (basecall.nf
     ubam rows, as run C) and `compare_pg.sh` (minimap2/samtools @PG lines vs run C; stop if they differ). Pileups for the
     comparisons: no filter and definition 2 (fixed 0.70) at MAPQ 40 for D, E, Cx and (definition 2, missing) A; definition
     4 (pipeline thresholds) for E. The pod5 subsets stay on /scratch only (user's decision; /projects had 216 GB free for
     264 GB); their md5s: `$B/pod5_durable_check/scratch_pod5.md5`.
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
