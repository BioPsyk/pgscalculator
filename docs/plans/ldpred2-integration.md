% LDpred2 integration plan

**Status:** Draft / not started
**Target version:** v2.2.0 (additive; sBayesR-only configs remain valid)
**Owner:** TBD
**Last updated:** 2026-05-13

**Design decisions locked in (2026-05-13):**

- Container: **multi-stage `r_builder` pattern** (mirrors the existing `java_builder` / `rust_builder` stages); see [§8].
- Method execution: support both **parallel** (one submission, multiple methods) and
  **sequential / incremental** (run sBayesR today, run LDpred2 tomorrow, finalize
  picks up both). See [§6.4] and [§6.5].

---

## 0. TL;DR

Add **LDpred2** alongside sBayesR as a second weight-generation method. Reuse
everything from the existing pipeline that is method-agnostic (prep-genotypes,
format-sumstat, filter-variants up to LD-ref filtering, calc-score,
combine-scores, calc-benchmark, finalize-output). Add:

1. an R-based `calc-ldpred2` step (single **genome-wide** job, not chr-parallel
   — this is mandated by the LDpred2 method, see [§4]),
2. an LDpred2-aware LD reference prep (`prep-ldref-ldpred2`) that consumes the
   precomputed HM3+/HM3 reference (`LD_with_blocks_chr{N}.rds` +
   `map_hm3_plus.rds`),
3. a `methods:` config key (and a matching `--methods` CLI override) so users
   can choose `sbayesr`, `ldpred2`, or both — **in the same submission or
   across separate submissions** (see [§6.5] on incremental runs),
4. a third SLURM weights profile (`weights_ldpred2`) submitted in parallel with
   the existing `weights_sbayesr` and `weights_benchmark` arrays when both
   methods are requested in one submission,
5. minimal schema changes to outputs: `augmented_sumstat.gz` gains
   `postEffect_ldpred2`; one `scores_<method>.gz` per method (see [§9]).

LDpred2 has **fundamentally different parallelism** than sBayesR: it operates
on a single genome-wide sparse correlation matrix (SFBM). Running per chromosome
is documented by the LDpred2 authors as producing *less accurate* scores. We
therefore implement it as a single non-array SLURM job with higher memory and
multiple threads, not as a 22-task array.

---

## 1. Goals & non-goals

### Goals

- Allow users to choose `methods: [sbayesr]`, `methods: [ldpred2]`, or
  `methods: [sbayesr, ldpred2]` in `config.yaml`.
- Reuse the existing variant-map / sumstat / scoring / finalize machinery as
  much as possible — LDpred2 is a *plug-in* weight producer, not a parallel
  pipeline.
- Keep all current sBayesR runs binary-compatible: a config without `methods:`
  defaults to `methods: [sbayesr]` (current behavior).
- Produce per-method scores and a single unified `augmented_sumstat.gz` where
  posterior columns are clearly suffixed by method.
- Run LDpred2-auto by default (no tuning set needed → matches the current
  sBayesR philosophy of "fully automatic, no held-out validation").

### Non-goals (for this PR)

- LDpred2-grid (needs an external validation phenotype + individuals). Can be
  added later with a separate `validation` config block.
- lassosum2. It uses the same `corr` and `df_beta` as LDpred2 so it would be
  almost free, but we deliberately scope it out for the first PR; track as
  a follow-up (see [§13]).
- Custom (non-precomputed) LD reference for LDpred2 (i.e. computing
  `snp_cor()` from PLINK files). We rely on the precomputed HM3/HM3+ LD
  references from the bigsnpr authors.
- Cross-ancestry support. The precomputed LD ref is European-only; users with
  non-European GWAS get a logged warning, same as today.

---

## 2. Background: how the current pipeline is shaped

For context (and to make the "what's reusable" mapping unambiguous), here is
the current per-sumstat step flow as implemented today:

```text
prep (sumstat-agnostic, run once per project)
├── prep-genotypes        → prep/genotypes/                  (method-agnostic)
├── prep-ldref            → prep/ldref/, prep/ldref_augmented (sbayesR-specific)
└── prep-inclusion-list   → prep/variant_map.tsv            (sbayesR-specific row set)

per sumstat
├── format-sumstat        → work/formatted/chr{N}.tsv        (method-agnostic)
├── filter-variants       → work/filtered/chr{N}_filtered.tsv (sbayesR-specific row set)
│                            + chr{N}_matched.tsv
├── weights (group)
│   ├── calc-posteriors   → work/posteriors/chr{N}.snpRes    (sbayesR-specific)
│   ├── format-posteriors → work/posteriors_mapped/chr{N}.snpRes (sbayesR-specific)
│   └── calc-benchmark    → work/benchmark/chr{N}_bench.sscore (method-agnostic)
├── calc-score            → work/scores/chr{N}.sscore        (method-agnostic, takes posteriors_mapped/)
└── finalize
    ├── combine-scores    → scores.gz
    └── finalize-output   → augmented_sumstat.gz, variant_map.gz, bench_score.gz
```

SLURM driver currently submits, in order: `format-sumstat` directly →
`sumstat` array → `weights_sbayesr` array → `weights_benchmark` array →
`score` array → `finalize` job.

The plan below changes the highlighted "(sbayesR-specific)" rows to be
method-parameterised, and keeps everything labelled "(method-agnostic)"
unchanged.

---

## 3. What LDpred2 needs (a one-page distillation)

From the bigsnpr vignette
([privefl.github.io/bigsnpr/articles/LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2))
and the reference implementation by the bigsnpr author at
[privefl/paper-infer/code/example-with-provided-LD.R](https://github.com/privefl/paper-infer/blob/main/code/example-with-provided-LD.R):

### Sumstat schema

Required per-variant columns:

```
chr, pos, a0, a1, beta, beta_se, n_eff   (+ rsid if matching by RSID)
```

where `a1` is the **effect allele** and `a0` is the **reference/other allele**.
`n_eff` = `N` for quantitative traits, `4 / (1/N_case + 1/N_control)` for
binary traits.

### LD reference

For European GWAS, the bigsnpr authors provide **precomputed** LD references:

| Set     | Variants  | URL                                                | Files |
|---------|-----------|----------------------------------------------------|-------|
| HapMap3 | 1,054,330 | <https://doi.org/10.6084/m9.figshare.19213299>     | `LD_with_blocks_chr{1..22}.rds`, `map_hm3.rds` |
| HapMap3+ | 1,444,196 | <https://doi.org/10.6084/m9.figshare.21305061>    | `LD_with_blocks_chr{1..22}.rds`, `map_hm3_plus.rds` |

Each `LD_with_blocks_chr{N}.rds` is an R `dsCMatrix` (sparse, blockwise). The
`map_*.rds` file holds variant metadata: `chr, pos, a0, a1, rsid, af_UKBB, ld,
block_id, pos_hg18, pos_hg38`. (`pos` is GRCh37/hg19.)

**Recommendation:** default to HM3+. It has 37 % more variants and matches the
modern paper. Users can override via `ldpred2.ld_variant_set: hm3`.

### Compute model

LDpred2 builds a single genome-wide sparse correlation matrix (SFBM) on disk
by iterating chr 1..22 and concatenating per-chromosome blocks into one file.
Then `snp_ldpred2_auto` runs many parallel Gibbs chains on this single SFBM.

```r
for (chr in 1:22) {
  corr_chr <- readRDS(sprintf("LD_with_blocks_chr%d.rds", chr))[ind.chr3, ind.chr3]
  if (chr == 1) corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  else          corr$add_columns(corr_chr, nrow(corr))
}
ldsc <- snp_ldsc(df_beta$ld, ld_size = nrow(map_ldref),
                 chi2 = (df_beta$beta / df_beta$beta_se)^2,
                 sample_size = df_beta$n_eff, ncores = NCORES)
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES)
# Filter "stable" chains, average their beta_est
keep    <- which(sapply(multi_auto, function(a) diff(range(a$corr_est))) >
                 0.95 * quantile(sapply(multi_auto, ...), 0.95))
beta_auto <- rowMeans(sapply(multi_auto[keep], function(a) a$beta_est))
```

### Resource footprint (HM3+, genome-wide)

- **Disk:** the on-disk SFBM (`tmp.sbk`) is roughly the size of the loaded LD
  blocks restricted to the matched variants (~5–15 GB for HM3+ after sumstat
  intersection, much less for low-power GWAS).
- **Memory:** the bigsnpr docs say "60 GB should be enough for HM3 with one
  million variants"; in practice, with HM3+ and `ncores = 16`, **plan for
  64 GB**.
- **Time:** the vignette quotes <5 min on 15 cores for one Gibbs run; with
  30 chains × 1000 iters and HM3+, expect 30–90 min wall on 16 cores.
- **Threads:** scales near-linearly to ~16 cores per `ncores` argument of
  `snp_ldpred2_auto`. More than that gives diminishing returns.

### Critical parallelism note

> "*you should run LDpred2 genome-wide. Just build the SFBM (the sparse LD
> matrix on disk) so that it contains selected variants for all chromosomes
> at once.*" — bigsnpr vignette

We must therefore **not** spread LDpred2 across 22 array tasks. One job,
genome-wide, with `ncores` set to the slurm `weights_ldpred2.cpus` value.

---

## 4. Reuse map: existing components → LDpred2

| Existing artifact                                   | Reused as-is for LDpred2? | Notes |
|-----------------------------------------------------|---------------------------|-------|
| `pgscalculator-v2.sh` wrapper (mounts, sbatch dispatch) | **Yes** + 1 new bind for LDpred2 LD ref dir + new `weights_ldpred2` profile |
| `bin/pgscalculator` CLI dispatcher                  | **Yes** + new `calc-ldpred2`, `format-posteriors-ldpred2`, `prep-ldref-ldpred2` commands |
| `bin/lib/common.sh` parse_config, log helpers, get_chromosomes, etc. | **Yes** | Already handles arbitrary `prefix.subkey` keys → `CFG_LDPRED2_*` will Just Work |
| `prep-genotypes`                                    | **Yes**, unchanged | Genotype-only, already method-agnostic |
| `prep-ldref` (sBayesR HM3 LD .info parsing + liftover augmentation) | **No** for the LDpred2 LD ref itself, but stays as-is for sBayesR | We add a *new* parallel step `prep-ldref-ldpred2` rather than mutate this one |
| `prep-inclusion-list`                               | **Conditional** | See [§5.2] — we keep its current output (sBayesR variant map) and add a separate `prep-inclusion-list-ldpred2` that builds the LDpred2 variant map. Both can run; if `methods` contains the corresponding method we run the corresponding prep |
| `format-sumstat`                                    | **Yes**, unchanged | Per-chr split of `cleaned_GRCh38.gz`; method-agnostic |
| `filter-variants`                                   | **Method-parameterised** | Filters per LD ref; produces per-method `filtered/chr{N}_filtered.tsv`. See [§5.3] |
| `calc-posteriors` (sBayesR)                         | **Yes** (sBayesR path only). LDpred2 has its own step (`calc-ldpred2`) |
| `format-posteriors` (sBayesR)                       | **Yes** (sBayesR path only). LDpred2 has its own step (`format-posteriors-ldpred2`) which emits files in the same on-disk format so calc-score can read them |
| `calc-benchmark`                                    | **Yes**, unchanged | Uses observed GWAS effects + LD pruning + plink2 `--score`; method-agnostic |
| `calc-score` (plink2 `--score`)                     | **Yes**, unchanged | We point it at the per-method `posteriors_mapped/<method>/chr{N}.snpRes` directory |
| `combine-scores`                                    | **Method-parameterised** | Produces per-method `work/scores_<method>/...`; see [§9] |
| `finalize-output`                                   | **Extended** | Adds `postEffect_ldpred2` column to `augmented_sumstat.gz`; emits `scores_ldpred2.gz` (or merged `scores.gz` — see [§9]) |
| Docker image (`docker/Dockerfile`)                  | **Extended** | Add `bigsnpr`, `bigreadr`, `bigsparser`, `runonce`, `argparser`, `stringr`, `ggplot2` R packages. Bump container tag to `2.1.0` |

### Anti-patterns to avoid

- **Do not** chromosome-parallelise LDpred2. Even if it were correct, the
  per-chr SFBM construction overhead exceeds the gain because each chunk of
  the precomputed `LD_with_blocks_chr{N}.rds` is already block-diagonal.
- **Do not** reuse the sBayesR `posteriors_mapped/chr{N}.snpRes` directory
  for LDpred2 output. Separate them so the score step is unambiguous about
  which method it is scoring. Use `posteriors_mapped_<method>/`.
- **Do not** silently change the default `whichn` semantics. LDpred2 needs
  `n_eff`; existing code uses `whichn: totalN | effectiveN`. We will compute
  `n_eff` in the LDpred2 R script from the per-variant N column (or fall
  back to the metadata-derived effective N), so this stays compatible.

---

## 5. Concrete design

### 5.1 New CLI / step taxonomy

| Step group (`--steps`) | Concrete steps (driver order)                                                                 |
|------------------------|----------------------------------------------------------------------------------------------|
| `prep`                 | `prep-genotypes`, `prep-ldref` *(sBayesR LD)*, `prep-ldref-ldpred2`, `prep-inclusion-list`, `prep-inclusion-list-ldpred2` |
| `sumstat`              | `format-sumstat`, `filter-variants` *(per method)*                                            |
| `weights`              | sBayesR: `calc-posteriors`, `format-posteriors`<br>LDpred2: `calc-ldpred2`, `format-posteriors-ldpred2`<br>Benchmark: `calc-benchmark` |
| `score`                | `calc-score` *(once per method)*                                                              |
| `finalize`             | `combine-scores`, `finalize-output`                                                           |

Concrete steps gain `--method <sbayesr|ldpred2>` where applicable. The
driver passes this through.

> **Naming choice:** keep `calc-posteriors` as the sBayesR step name (it has
> shipped) and introduce `calc-ldpred2` for the new step rather than
> renaming `calc-posteriors → calc-posteriors-sbayesr`. This is least
> disruptive to existing scripts and docs. The general term "posteriors" is
> retained where it appears in user-facing output (`postEffect_sbayesr`,
> `postEffect_ldpred2`).

### 5.2 LDpred2 LD reference handling (`prep-ldref-ldpred2`)

**Input** (provided by user, downloaded once and stored on shared faststorage):

```
${CFG_LDPRED2_LD_DIR}/
├── LD_with_blocks_chr1.rds
├── LD_with_blocks_chr2.rds
├── ...
├── LD_with_blocks_chr22.rds
└── map_hm3_plus.rds            (or map_hm3.rds if hm3)
```

**Step responsibility:**

1. Validate that all 22 `LD_with_blocks_chr{N}.rds` files exist and that
   `map_hm3_plus.rds` is loadable (use a tiny R helper script: `Rscript -e
   'm <- readRDS("..."); cat(nrow(m))'`).
2. Convert `map_hm3_plus.rds` → a TSV (`prep/ldref_ldpred2/map.tsv`) with
   columns `chr, pos_b37, pos_b38, a0, a1, rsid, af_UKBB, ld, block_id`.
   This is the LDpred2 equivalent of `prep/ldref_augmented/` and is what we
   use for the LDpred2 variant map join.
3. Cache the conversion. If `prep/ldref_ldpred2/map.tsv` already exists and
   the `.rds` is older, skip.

We do **not** rewrite the `.rds` files (they stay binary and bound into the
container at runtime).

### 5.3 LDpred2 variant map (`prep-inclusion-list-ldpred2`)

Currently `prep-inclusion-list` joins `prep/ldref_augmented/chr{N}_ld_augmented.tsv`
with `prep/genotypes/chr{N}_pvar_fmt` on `chr:pos + alleles` to produce a
single `prep/variant_map.tsv`. The LDpred2 variant set is **different** (HM3+
has ~1.44M variants vs the sBayesR `band_ukb_10k_hm3` reference's ~1.13M),
so the matched/non-matched sets differ.

Solution: produce **two** variant maps in prep, keyed on the method:

```
prep/variant_map_sbayesr.tsv      (current; rename existing variant_map.tsv → variant_map_sbayesr.tsv)
prep/variant_map_ldpred2.tsv      (new; built from prep/ldref_ldpred2/map.tsv)
```

**No symlink for back-compat.** Symlinks don't survive cleanly through
Singularity/Docker bind mounts on every site (broken target resolution,
copy-on-export quirks). Instead, the rename is handled the same way the
sumstat-dir reshuffle already is: a `migrate_prep_variant_map_to_sbayesr()`
helper runs at the start of every `prep-inclusion-list` invocation and, if it
sees a legacy `prep/variant_map.tsv` (or `prep/variant_map/`) on disk with no
new sBayesR-named counterpart, renames it in place. Idempotent and safe to
re-run. All in-pipeline readers (`filter_variants.sh`, `format_posteriors.sh`,
`finalize_output.sh`) point at the new sBayesR-suffixed path directly — no
double-lookup with fallback.

The schema for `variant_map_ldpred2.tsv` mirrors `variant_map_sbayesr.tsv`:

```
chr  pos_b37  pos_b38  geno_snpid  geno_a1  geno_a2  ldref_snpid  ldref_a1  ldref_a2  ldref_a2freq  ld  block_id
```

Allele-convention reconciliation: the bigsnpr LDpred2 map uses `a0` (reference
allele) and `a1` (alternative / effect allele), and ships `af_UKBB` =
frequency of `a1`. The existing sBayesR variant map follows the GCTB
convention where `ldref_a1` is the effect allele and `ldref_a2` is the other
allele, with `ldref_a2freq` = frequency of the *other* allele. To keep
downstream readers method-agnostic we map:

| LDpred2 map column | variant_map_ldpred2.tsv column |
|--------------------|-------------------------------|
| `a1`               | `ldref_a1` (effect allele)    |
| `a0`               | `ldref_a2` (other allele)     |
| `af_UKBB`          | `1 - af_UKBB` → `ldref_a2freq`|

`ld` and `block_id` are carried through unchanged and are consumed by the
LDpred2 R script during LD-score regression and SFBM construction.

### 5.4 `filter-variants` becomes method-parameterised

Today the step writes to `work/filtered/chr{N}_filtered.tsv` (one place,
implicitly sBayesR-shaped). New shape:

```
work/filtered_sbayesr/chr{N}_filtered.tsv     (current location, renamed)
work/filtered_sbayesr/chr{N}_matched.tsv      (kept for finalize)
work/filtered_ldpred2/chr{N}_filtered.tsv     (new)
work/filtered_ldpred2/chr{N}_matched.tsv      (new)
```

The filter-variants step takes a `--method` flag and reads
`prep/variant_map_<method>.tsv`. The bulk of the logic (column derivation,
EAF/B/SE fill, bad-value drop) is unchanged: the awk parser already
auto-detects columns by name, so the 12-col LDpred2 variant map (`ld`,
`block_id` appended) is consumed transparently — the same code path picks
up the first 10 columns it knows about and ignores the trailing two.

The `--method ldpred2` invocation hard-fails if
`prep/variant_map_ldpred2.tsv` is missing, pointing the user at
`prep-inclusion-list-ldpred2` (which is itself opt-in via `ldpred2.ld_dir`,
mirroring the §5.3 prep steps).

> **Migration helper:** add a `migrate_filtered_dirs()` that, if
> `work/filtered/` exists but `work/filtered_sbayesr/` does not, renames the
> directory. Same pattern as the existing `migrate_sumstat_step_dir`. The
> helper also runs as part of `migrate_sumstat_all_step_dirs`, so any call
> site that already invokes that wrapper (e.g. `run_pipeline.sh`) gets the
> rename for free without per-step bookkeeping.

The driver runs `filter-variants` **once per requested method**. Both can
run inside the same `sumstat` array task (they're cheap; chr-parallel
already).

**Scope deferral (this commit, Phase 1 §5.4):**

- The per-sumstat-level `${sumstat_dir}/variant_map.tsv` filename stays
  unsuffixed for now. Splitting it into `variant_map_<method>.tsv` requires
  updating every downstream reader (`format_posteriors.sh`,
  `finalize_output.sh`, `calc_posteriors.sh`, etc.), which is part of the
  downstream-method work tracked in §5.6–§5.7 / §9. Since only one method
  runs per filter-variants invocation today, the file is owned by whichever
  method wrote it last; LDpred2 + sBayesR coexistence at the sumstat-mapfile
  level lands with the §9 finalize discovery work.
- The driver-side per-method dispatch (one `filter-variants` invocation per
  active method) belongs to the wrapper commit in §6.4 and is intentionally
  not threaded through `bin/lib/steps/run_pipeline.sh` here. For now the
  step-runner default is sbayesr; passing `CFG_METHOD=ldpred2` (or
  `--method ldpred2`) into a single `pgscalculator filter-variants` call
  exercises the LDpred2 path end-to-end.

### 5.5 The LDpred2 R script (`bin/lib/scripts/run_ldpred2.R`)

Adapted from the comorment/containers reference at
[containers/scripts/pgs/LDpred2/ldpred2.R](https://github.com/comorment/containers/blob/main/scripts/pgs/LDpred2/ldpred2.R)
plus the bigsnpr-author reference at
[paper-infer/code/example-with-provided-LD.R](https://github.com/privefl/paper-infer/blob/main/code/example-with-provided-LD.R).

We deliberately **do not** require a genotype `.rds`/`.bk` bigSNP file like
the comorment script. We score with PLINK2 downstream (same as sBayesR), so
the R script just needs to emit per-chr posterior effects keyed by
`ldref_snpid`. The downstream pipeline maps those to `geno_snpid` via the
LDpred2 variant map.

**Inputs (CLI args):**

```
--sumstat-dir       <work/filtered_ldpred2/>                  # has chr{N}_filtered.tsv
--ld-dir            <CFG_LDPRED2_LD_DIR>                       # has LD_with_blocks_chr{N}.rds
--ld-meta           <CFG_LDPRED2_LD_META_FILE>                 # map_hm3_plus.rds
--variant-map       <prep/variant_map_ldpred2.tsv>             # rsid<->geno_snpid + ld + block_id
--out-dir           <work/posteriors_ldpred2/>                  # writes chr{N}.snpRes
--mode              auto                                        # auto | inf
--shrink-corr       0.95
--allow-jump-sign   false
--hyper-p-max       0.2
--hyper-p-length    30
--burn-in           500
--num-iter          500
--use-mle           true
--seed              1
--ncores            16
--genotype-build    GRCh37                                      # so we know which pos column to match on
```

**Algorithm:**

1. `library(bigsnpr); library(bigreadr); library(bigsparser)`.
2. Read `map_<set>.rds` → `map_ldref`.
3. Read per-chr filtered sumstat files, `rbind` them → `sumstats` data.frame
   with columns `chr, pos, a0, a1, rsid, beta, beta_se, n_eff`.
   - We derive `n_eff` from the per-variant `N` column (or `Neff` if
     available). If both case/control counts are present in the cleansumstats
     metadata, we override `n_eff = 4 / (1/N_case + 1/N_control)`.
4. `df_beta <- snp_match(sumstats, map_ldref)` (match by chr/pos or rsid
   depending on `--merge-by-rsid`; default match by chr/pos using the
   build-appropriate column).
5. Standard QC (from the bigsnpr vignette):
   ```r
   sd_ldref <- sqrt(2 * df_beta$af_UKBB * (1 - df_beta$af_UKBB))
   sd_ss    <- 2 / sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
   is_bad   <- sd_ss < 0.5 * sd_ldref | sd_ss > sd_ldref + 0.1 |
               sd_ss < 0.05 | sd_ldref < 0.05
   df_beta  <- df_beta[!is_bad, ]
   ```
6. Build genome-wide SFBM by `for (chr in 1:22) corr$add_columns(...)`.
7. `snp_ldsc()` → starting `h2`.
8. `snp_ldpred2_auto()` (or `_inf`).
9. Filter "stable" chains, average their `beta_est` → vector `beta_auto`.
10. Write per-chr `posteriors_ldpred2/chr{N}.snpRes` with schema:
    ```
    SNP A1 A2 Freq Effect  (no SE/PIP — LDpred2-auto's posterior per-variant
    SD is not directly comparable to sBayesR's; we add a separate `PIP` column
    from `auto$postp_est` for symmetry but mark it as method-specific in docs)
    ```
    Actually keep the same schema as sBayesR's `posteriors/chr{N}.snpRes`
    (`SNP A1 A2 Freq Effect SE PIP`) — fill `SE = NA`, `PIP = postp_est`. The
    downstream `format-posteriors` step only consumes `SNP`, `A1`, `A2`,
    `Freq`, `Effect`, so this stays compatible.
11. Diagnostic plot of chain convergence → `logs/ldpred2_chains.png`.

**Failure modes** (all soft-fail, write placeholder, mirroring sBayesR
behavior):

| Marker file                                | Reason |
|-------------------------------------------|--------|
| `work/posteriors_ldpred2/FAILED_match`    | `snp_match` matched <100 variants |
| `work/posteriors_ldpred2/FAILED_sfbm`     | SFBM construction errored (corrupted .rds or out-of-disk) |
| `work/posteriors_ldpred2/FAILED_ldsc`     | LDSC returned NA / non-finite h2 |
| `work/posteriors_ldpred2/FAILED_ldpred2`  | All chains diverged (`keep` is empty) |

When any of these is written we still emit header-only
`chr{1..22}.snpRes` so the score step can produce header-only placeholder
scores (consistent with the existing sBayesR failure pattern documented in
`docs/plans/general-design.md` §"Failure handling").

### 5.6 `format-posteriors-ldpred2`

Almost identical to the existing `format-posteriors`. Reads
`work/posteriors_ldpred2/chr{N}.snpRes`, joins with
`prep/variant_map_ldpred2.tsv` to add `geno_snpid`, writes
`work/posteriors_mapped_ldpred2/chr{N}.snpRes`.

The cleanest implementation is to **refactor** the existing
`bin/lib/steps/format_posteriors.sh` so the variant_map path and the input
/ output directories are parameters of `run_format_posteriors`, then have
two thin step entry points (`run_format_posteriors_sbayesr`,
`run_format_posteriors_ldpred2`) that just supply the method-specific
paths. Less duplicated code, easier to maintain.

### 5.7 `calc-score` (method-aware, otherwise unchanged)

The score step today reads `work/posteriors_mapped/chr{N}.snpRes`. We
change it to read `work/posteriors_mapped_<method>/chr{N}.snpRes` and write
to `work/scores_<method>/chr{N}.sscore`. The driver invokes it once per
method.

`CFG_SCORE_COLUMNS` (`1 2 5` = ID, A1, Effect) stays the same because both
methods emit posteriors in the same column order.

---

## 6. Driver / SLURM topology

### 6.1 Current weights phase (one sumstat, --sbatch)

```
driver
  ├── (block) wait sumstat array (22 tasks)
  ├── (block) wait weights_sbayesr array (22 tasks)
  ├── (block) wait weights_benchmark array (22 tasks)
  └── (block) wait score array (22 tasks)
```

### 6.2 New weights phase

```
driver
  ├── (block) wait sumstat array (22 tasks, internally runs filter-variants for each method)
  ├── (block) wait weights_sbayesr array (22 tasks)        ┐
  ├── (block) wait weights_ldpred2 single job (no array)   ├── submitted concurrently
  ├── (block) wait weights_benchmark array (22 tasks)      ┘
  ├── (block) wait score_sbayesr array (22 tasks)          ┐
  ├── (block) wait score_ldpred2 array (22 tasks)          ┘ submitted after weights done
  └── (block) wait finalize single job
```

Submit-concurrently and wait-for-all is the same pattern the driver already
uses (just with one more job to track). The `submit_array_for_step()`
function in `pgscalculator-v2.sh:912` already supports per-`step_profile`
behavior — we add a new `weights_ldpred2` branch and a sibling
`submit_single_job_for_step()` for the non-array case (genome-wide).

### 6.3 SLURM resource defaults

Add to the template:

```yaml
slurm:
  # ... existing entries ...
  weights_ldpred2:    { mem: 64g, cpus: 16, time: '4:00:00' }    # single job, no array
  score_sbayesr:      { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }
  score_ldpred2:      { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }
```

For back-compat: if `score:` is present but `score_sbayesr` / `score_ldpred2`
are not, both inherit from `score:`. Keep the old key alive — same trick as
`genofile` ↔ `genotype_manifest`.

### 6.4 Method dispatch in the driver

The set of active methods for a given submission comes from (highest
precedence wins):

1. `--methods sbayesr,ldpred2` CLI flag on the wrapper (overrides config),
2. `methods: [sbayesr, ldpred2]` in the config file,
3. default: `[sbayesr]` (preserves current behavior for existing configs).

Parsing:

```bash
# CLI override (new)
if [[ -n "$methods_cli" ]]; then
  cfg_methods="$methods_cli"
else
  cfg_methods=$(parse_yaml_list "methods" "$config_file_host")
  [[ -z "$cfg_methods" ]] && cfg_methods="sbayesr"
fi
```

Then around `pgscalculator-v2.sh:1239`:

```bash
if has_method sbayesr "$cfg_methods"; then submit_array_for_step weights_sbayesr;     fi
if has_method ldpred2 "$cfg_methods"; then submit_single_job_for_step weights_ldpred2; fi
submit_array_for_step weights_benchmark   # always (orthogonal to method)

# Score phase: one array per active method
if has_method sbayesr "$cfg_methods"; then submit_array_for_step score_sbayesr; fi
if has_method ldpred2 "$cfg_methods"; then submit_array_for_step score_ldpred2; fi
```

No change for users on sBayesR-only configs.

### 6.5 Sequential / incremental method runs (small-HPC friendly)

Smaller HPC sites may not want to submit `weights_sbayesr` (22 tasks
× 20 GB) and `weights_ldpred2` (single job × 64 GB) at the same time. The
pipeline is designed so the user can run **methods one at a time across
separate submissions**, with the second submission only doing the
incremental work:

```bash
# Day 1: only sBayesR
./pgscalculator-v2.sh --config config.yaml \
  --steps sumstat,weights,score,finalize \
  --methods sbayesr \
  -i /path/to/sumstat_TRAIT --sbatch

# Day 2: add LDpred2 (sBayesR results untouched, finalize re-runs to merge)
./pgscalculator-v2.sh --config config.yaml \
  --steps sumstat,weights,score,finalize \
  --methods ldpred2 \
  -i /path/to/sumstat_TRAIT --sbatch
```

The second submission's behavior:

| Step                          | Behavior on Day-2 submission                                                                              |
|------------------------------|-----------------------------------------------------------------------------------------------------------|
| `format-sumstat`              | Skipped (existing `.completed` marker; no per-method state).                                              |
| `filter-variants --method sbayesr` | Not invoked (sBayesR not requested today).                                                          |
| `filter-variants --method ldpred2` | Runs (`work/filtered_ldpred2/` is missing or marker absent). Cheap (chr-parallel).                  |
| `calc-posteriors` (sBayesR)   | Not invoked (sBayesR not requested today). `work/posteriors_sbayesr/` is left untouched.                  |
| `format-posteriors` (sBayesR) | Not invoked.                                                                                              |
| `calc-ldpred2`                | Runs (single genome-wide job). Writes `work/posteriors_ldpred2/`.                                         |
| `format-posteriors-ldpred2`   | Runs. Writes `work/posteriors_mapped_ldpred2/`.                                                           |
| `calc-benchmark`              | Skipped (existing `.completed` marker; method-agnostic).                                                  |
| `calc-score` (sBayesR)        | Not invoked.                                                                                              |
| `calc-score` (LDpred2)        | Runs. Writes `work/scores_ldpred2/`.                                                                      |
| `combine-scores`              | Runs **per method**: regenerates `scores_ldpred2.gz`. Does **not** touch `scores_sbayesr.gz`.             |
| `finalize-output`             | Runs in **discovery mode**: scans `work/posteriors_mapped_*/`, regenerates `augmented_sumstat.gz` with **both** `postEffect_sbayesr` (from Day 1) **and** `postEffect_ldpred2` (from today). |

The mechanics that make this work:

1. **Per-method `.completed` markers.** Each method-specific step directory
   gets its own marker:
   - `work/filtered_<method>/.completed`
   - `work/posteriors_<method>/.completed`
   - `work/posteriors_mapped_<method>/.completed`
   - `work/scores_<method>/.completed`

   The existing `check_step_completed()` /`mark_step_completed()` helpers
   in `common.sh` already key off the step directory, so this works by
   construction once we route each method to its own directory (which we
   already planned to do in [§5]).

2. **Driver respects `--methods`/`methods:` as the *active set*, not as a
   re-run trigger.** If a method isn't in the active set, the driver
   simply skips dispatching its arrays/jobs. It doesn't delete or
   invalidate that method's existing outputs.

3. **`finalize-output` runs in discovery mode** (new in v2.2.0):
   - Reads `methods:` from the config as a *hint* for column ordering.
   - But for each method, **only includes the column if
     `work/posteriors_mapped_<method>/` contains non-empty outputs**.
   - For `scores_<method>.gz`, only emits files for methods that have
     non-empty `work/scores_<method>/`.

   This means:
   - Day 1 produces `scores_sbayesr.gz` + `augmented_sumstat.gz` with
     `postEffect_sbayesr` only.
   - Day 2 produces `scores_ldpred2.gz` + overwrites
     `augmented_sumstat.gz` to contain **both** `postEffect_sbayesr` and
     `postEffect_ldpred2` (using the Day-1 posteriors on disk).
   - `scores_sbayesr.gz` from Day 1 is left untouched.

4. **`benchEffect` is method-agnostic** (it's the LD-pruned GWAS effect),
   so it's filled identically in either run.

### 6.6 Re-running a single method

If the user wants to force a re-run of just LDpred2 (e.g. new LD ref):

```bash
./pgscalculator-v2.sh --config config.yaml \
  --steps weights,score,finalize \
  --methods ldpred2 \
  --force \
  -i /path/to/sumstat_TRAIT --sbatch
```

`--force` already exists and clears step markers. Combined with
`--methods ldpred2`, it only clears markers in method-specific
directories (`work/posteriors_ldpred2/`, `work/scores_ldpred2/`, etc.).
The sBayesR side is untouched.

### 6.7 Dispatch decision matrix

The driver's decision for any single submission boils down to two
inputs (`active_methods`, `requested_steps`) yielding the set of
SLURM jobs:

| `active_methods`        | `--steps`                       | Driver submits                                                                  |
|-------------------------|---------------------------------|--------------------------------------------------------------------------------|
| `[sbayesr]`             | `sumstat,weights,score,finalize` | sumstat array → weights_sbayesr + weights_benchmark → score_sbayesr → finalize  |
| `[ldpred2]`             | `sumstat,weights,score,finalize` | sumstat array → weights_ldpred2 + weights_benchmark → score_ldpred2 → finalize  |
| `[sbayesr,ldpred2]`     | `sumstat,weights,score,finalize` | sumstat array → weights_sbayesr ∥ weights_ldpred2 ∥ weights_benchmark → score_sbayesr ∥ score_ldpred2 → finalize |
| `[sbayesr]`             | `finalize`                       | finalize only (discovery picks up any prior outputs)                           |
| `[]` (or omitted method but explicit `--steps finalize`) | `finalize`                       | finalize only; discovery determines columns/files                              |

(The `∥` symbol denotes "submitted in parallel; driver waits for all".)

---

## 7. Config schema additions

Additive. No existing keys change.

```yaml
# Existing keys (unchanged)
input: /path/to/cleansumstats/output
outdir: /path/to/output
genodir: /path/to/genotypes
genofile: /path/to/genotype_manifest.tsv
lddir: /path/to/band_ukb_10k_hm3        # sBayesR LD reference (existing)
genotype_build: GRCh37
liftover_reference: /path/to/liftover.txt.gz
whichn: totalN
sbayesr: { ... }                         # existing
benchmark: { ... }                       # existing

# === NEW (additive) ===
# Which posterior methods to compute when this config is used.
# Omit or set to [sbayesr] for current behavior.
# Can be overridden per submission via --methods CLI flag (see §6.5 for
# the "run sBayesR today, LDpred2 tomorrow" workflow).
methods: [sbayesr, ldpred2]

ldpred2:
  mode: auto                             # auto | inf  (grid is future work)
  ld_variant_set: hm3_plus               # hm3 | hm3_plus
  ld_dir: /path/to/ldpred2_ref/ldref_hm3_plus    # contains LD_with_blocks_chr{N}.rds
  ld_meta_file: /path/to/ldpred2_ref/map_hm3_plus.rds
  shrink_corr: 0.95
  allow_jump_sign: false
  hyper_p_max: 0.2
  hyper_p_length: 30
  burn_in: 500
  num_iter: 500
  use_mle: true
  coef_shrink_keep: 0.95                 # for chain filtering
  merge_by_rsid: false                   # match sumstat<->map by chr:pos+alleles (recommended)
  seed: 1
  threads: 16

slurm:
  # ... existing entries ...
  weights_ldpred2:    { mem: 64g, cpus: 16, time: '4:00:00' }
  score_sbayesr:      { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }   # inherits from score: if missing
  score_ldpred2:      { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }
```

### Parsing in `common.sh`

`parse_config()` already turns nested keys into `CFG_<SECTION>_<KEY>`
variables, so we get `CFG_LDPRED2_MODE`, `CFG_LDPRED2_LD_DIR`,
`CFG_LDPRED2_SHRINK_CORR`, etc. for free.

The only addition needed in `common.sh` is `parse_yaml_list()` for the
`methods: [sbayesr, ldpred2]` line (inline array — the current parser only
handles scalar values). Reuse the existing `parse_yaml_list` from
`pgscalculator-v2.sh:255` (it handles `modules:` already, but block-style;
extend it to also handle inline `[a, b, c]` syntax).

---

## 8. Container changes (`docker/Dockerfile`)

### 8.1 Multi-stage `r_builder` (chosen design)

R-package installation is **isolated in a dedicated builder stage**, mirroring
the existing `java_builder` (Nextflow) and `rust_builder` (b3sum) stages. The
final runtime image just `COPY --from=r_builder`s the compiled R library tree.

**Why:** `bigsnpr` is a heavy install. From-source it pulls `Rcpp`,
`RcppArmadillo`, `bigstatsr`, `bigsparser`, `bigassertr`, `rmio`,
`bigparallelr`, `Matrix`, `data.table` and friends — easily 20–30 min on a
cold build. With the multi-stage pattern:

- changes to `apt-get install` lines in the final stage **do not bust** the
  R package cache,
- changes to `bin/`, `assets/`, plink/gctb/PRScs versions **do not bust** the
  R package cache,
- the R install only re-runs when *its own stage* changes (new R package,
  new system -dev lib needed for R, base image change).

Additional speed-up: use **`pak`** (parallel installer that pulls prebuilt
binaries from Posit Public Package Manager where available) instead of
`install.packages`. Empirically cuts the bigsnpr cold install from ~25 min
to ~3–5 min, and is much faster on re-installs too.

**Layer ordering inside the final stage** is also tightened: R-runtime libs
go in the same apt layer as everything else; the `COPY --from=r_builder`
goes just *after* that apt layer so it's stable, and the per-iteration
churn (`COPY bin/`, `COPY assets/`, `chmod`) stays at the very bottom.

### 8.2 Concrete Dockerfile sketch

```dockerfile
# ============================================================
# Existing builder stages (unchanged)
# ============================================================
FROM gradle:8.5-jdk11 AS java_builder
# ... (unchanged: Nextflow install) ...

FROM rust:1.84-slim-bookworm AS rust_builder
# ... (unchanged: b3sum install) ...

# ============================================================
# NEW: R packages builder
# ============================================================
FROM eclipse-temurin:11 AS r_builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base r-base-dev \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libbz2-dev liblzma-dev libgomp1 \
        build-essential pkg-config \
    && rm -rf /var/lib/apt/lists/*

# pak: parallel installer with binary repos (Posit PPM); much faster than
# install.packages() from source.
RUN R -e "install.packages('pak', repos = sprintf( \
      'https://r-lib.github.io/p/pak/stable/%s/%s/%s', \
      .Platform\$pkgType, R.Version()\$os, R.Version()\$arch))"

# All R packages needed for LDpred2.
# Pin bigsnpr >= 1.12.0 (validated LDpred2-auto) and bigsparser >= 0.6
# (compact SFBM format).
# pak version-pin syntax is 'pkg@>=X.Y.Z' (https://pak.r-lib.org/reference/pak_package_sources.html)
RUN R -e "pak::pkg_install(c( \
      'bigsnpr@>=1.12.0', \
      'bigsparser@>=0.6', \
      'bigreadr', 'runonce', 'argparser', 'stringr', \
      'ggplot2', 'cowplot', 'data.table', \
      'tibble', 'tidyr', 'dplyr'))"

# ============================================================
# Final runtime image
# ============================================================
FROM eclipse-temurin:11

COPY --from=java_builder /usr/local/bin/nextflow /usr/local/bin/nextflow
COPY --from=java_builder /root/.nextflow /root/.nextflow
COPY --from=rust_builder /usr/local/cargo/bin/b3sum /usr/bin/b3sum

WORKDIR /pgscalculator

# System libs. As-shipped in v2.2.0 we kept the existing -dev variants
# (libcurl4-openssl-dev, libssl-dev, libxml2-dev, libbz2-dev, liblzma-dev)
# from the pre-r_builder Dockerfile instead of switching to runtime-only
# (libcurl4 / libssl3 / libssl1.1 / libxml2 / libbz2-1.0 / liblzma5).
# Reason: `eclipse-temurin:11` is unpinned in this repo and may resolve to
# Focal (libssl1.1) or Jammy (libssl3) at build time; the -dev libs drag
# in the correct runtime variants regardless. Slimming to runtime-only is
# tracked as a follow-up (pin the base image + verify libssl name first).
# r-base-dev IS dropped from the final stage (compilers only needed in
# r_builder, since runtime R packages come pre-built via COPY).
RUN apt-get update && apt-get install -y \
        wget unzip curl gawk graphviz dos2unix pigz vim parallel git \
        r-base \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libbz2-dev liblzma-dev libgomp1 \
        --no-install-recommends \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Pull the prebuilt R library tree from the builder.
# Debian/Ubuntu R installs system packages to /usr/lib/R/site-library and
# user packages (what pak installs into without root config) to
# /usr/local/lib/R/site-library. Take both to be safe.
COPY --from=r_builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=r_builder /usr/lib/R/site-library       /usr/lib/R/site-library

# ... PRScs, gctb, plink, plink2 installs (unchanged) ...

# User setup (unchanged) ...

# pgscalculator v2 CLI (unchanged — stays at the bottom so dev iterations
# on bin/ never bust anything above)
COPY bin/   /pgscalculator/bin/
COPY assets/ /pgscalculator/assets/
COPY config.template.yaml /pgscalculator/
RUN chmod +x /pgscalculator/bin/pgscalculator
ENV PATH="/pgscalculator/bin:${PATH}"
```

### 8.3 Caveats

- **R version coupling.** Both `r_builder` and the final stage must end up
  with the same R version, otherwise compiled packages from the builder
  won't load (ABI mismatch). Since both `FROM eclipse-temurin:11` and both
  install Debian's `r-base`, this is automatic — but if we ever pin one
  base image, pin both.
- **`COPY` of `/usr/lib/R/site-library`** is only relevant if `pak` decided
  to install there for whatever reason. If pak only uses
  `/usr/local/lib/R/site-library`, the second COPY is a no-op. Harmless.
- **Don't `RUN R -e ...` for runtime smoke-test in the final stage.** Save
  that for the test phase (Phase 1's smoke test) — adding it as a RUN line
  in the Dockerfile would defeat caching.
- **Optional BuildKit cache mount.** If we feel pain even with `pak`, we
  can add `--mount=type=cache,target=/root/.cache/R` to the
  `pak::pkg_install` line so downloaded tarballs survive even a stage
  rebuild. Skip for v2.2.0 unless needed.

### 8.4 Tag bump & wrapper

Bump container tag to `2.1.0-amd64` (or `2.2.0-amd64` to align with
pipeline version — decide at release time). Update default in
`pgscalculator-v2.sh:1673` (`singularity_image_tag`).

### 8.5 Mounts

`pgscalculator-v2.sh` needs an additional bind for the LDpred2 LD ref dir:

```bash
ldpred2_lddir_host=$(realpath "${cfg_ldpred2_ld_dir}")
ldpred2_lddir_container="/pgscalculator/ldpred2_ref"
mount_opts="${mount_opts} ${mountflag} ${ldpred2_lddir_host}:${ldpred2_lddir_container}:ro"
```

Mount as **read-only** (`:ro`) — the LD ref is large and immutable.
`map_*.rds` is usually shipped in the same directory, so one bind covers
both.

---

## 9. Output schema changes

### `augmented_sumstat.gz`

Current schema (11 columns):

```
RSID  EffectAllele  OtherAllele  B  SE  Z  P  EAF  MAF  postEffect  benchEffect
```

New schema (13 columns when both methods run; 12 when only one runs):

```
RSID  EffectAllele  OtherAllele  B  SE  Z  P  EAF  MAF  postEffect_sbayesr  postEffect_ldpred2  benchEffect  postp_ldpred2
```

- `postEffect_sbayesr` ← `postEffect` of today (column renamed).
- `postEffect_ldpred2` ← averaged `beta_est` across kept LDpred2-auto chains
  (mapped to `ldref_snpid`).
- `postp_ldpred2` ← averaged `postp_est` (posterior inclusion probability;
  optional, gives users fine-mapping signal). Drop this if it's too noisy
  in practice.
- **Discovery-driven column inclusion** (see [§6.5]): `finalize-output`
  emits a `postEffect_<method>` column for every method that has non-empty
  outputs in `work/posteriors_mapped_<method>/`, regardless of what
  `methods:` says. This is what enables incremental runs (Day-1 sBayesR
  produces a 1-method file; Day-2 LDpred2 overwrites the file with both
  columns).
- For methods that *should* be there per config but have no outputs
  (because they haven't been run yet, or all chains diverged), the column
  is **omitted**. Header documents which methods are present.
  - **Alternative considered:** always emit all configured-methods
    columns, NA-fill the missing ones. Simpler for consumers but breaks
    the "Day-1 then Day-2 incrementally extends" property: the Day-1 file
    would carry an NA `postEffect_ldpred2` column that gets *replaced*
    rather than *added* on Day-2. **Recommendation:** omit-when-absent
    (matches discovery semantics); fix downstream consumers to read the
    header.

### Scores

**Option A** (recommended): one `scores.gz` per method.

```
scores_sbayesr.gz       (IID, SCORE_SUM, ALLELE_CT, N_VARIANTS)
scores_ldpred2.gz       (IID, SCORE_SUM, ALLELE_CT, N_VARIANTS)
bench_score.gz          (IID, ALLELE_CT, SCORE1_SUM)
```

For back-compat: when `methods == [sbayesr]`, emit both
`scores_sbayesr.gz` and a symlink `scores.gz → scores_sbayesr.gz`.

**Option B** (rejected): one combined `scores.gz` with
`SCORE_SBAYESR_SUM, SCORE_LDPRED2_SUM`. Reduces files but couples QC and
makes per-method N_VARIANTS columns awkward (different SNP sets per method).

### `variant_map.gz`

Unchanged. The user-facing variant_map.gz remains the sBayesR-shaped
mapfile. If the user wants the LDpred2 variant set, they can read
`prep/variant_map_ldpred2.tsv` directly. We do not need both as user
outputs.

### `details/steps.tsv`

The per-step variant-counts file gains rows for the LDpred2 pipeline
(`ldpred2-match`, `ldpred2-qc`, `ldpred2-sfbm`, `ldpred2-auto`,
`ldpred2-score`). Same TSV format.

---

## 10. Failure handling

Mirror the existing chromosome-level continue-on-failure policy
(`docs/plans/general-design.md` §"Failure handling"), adapted:

- **Per-method failure isolation.** If LDpred2 fails entirely (e.g. all
  chains diverge), the driver logs a warning, writes a global LDpred2
  `FAILED` marker, and *still produces sBayesR + benchmark outputs*.
  `scores_ldpred2.gz` is header-only and `augmented_sumstat.gz` simply
  omits the `postEffect_ldpred2` column (discovery sees empty outputs).
- **Per-chr score failures** are unchanged (same as today; the score array
  for LDpred2 follows the same pattern as for sBayesR).
- **prep-ldref-ldpred2 failure** is fatal at prep time (no point running
  per-sumstat if the LD ref is broken).
- **Incremental-run interaction with failures.** If Day-1 (sBayesR)
  produced full outputs and Day-2 (LDpred2) fails entirely, finalize on
  Day-2 still emits an `augmented_sumstat.gz` with the sBayesR column
  intact (discovery picks it up from Day-1's `posteriors_mapped_sbayesr/`).
  The user can retry just LDpred2 (`--methods ldpred2 --force`) without
  losing the sBayesR work.

---

## 11. Tests & smoke checks

### Unit tests (under `tests/unit/`)

| New script                                  | Tests |
|---------------------------------------------|-------|
| `test_ldpred2_match.sh`                     | `snp_match` against a tiny fake `map_hm3.rds` + tiny sumstat → expected row count |
| `test_format_for_ldpred2.sh`                | Per-chr filtered TSV → R input data.frame derivation (n_eff, chr coerced to int, etc.) |
| `test_format_posteriors_ldpred2.sh`         | LDpred2 `.snpRes` → mapped via `variant_map_ldpred2.tsv` |
| `test_prep_ldref_ldpred2.sh`                | `map_*.rds` → TSV conversion produces correct columns |
| `test_method_dispatch.sh`                   | Driver dispatch table: methods=[a], [b], [a,b] → expected set of submitted profiles |

### Smoke test (under `tests/smoke/v2.2-<date>/`)

Add a smoke test alongside the existing `v2.1-2026-01-11/`:

- `config.ldpred2.yaml`: methods=[ldpred2], small public sumstat, public
  HM3 LD ref.
- `config.both.yaml`: methods=[sbayesr, ldpred2].
- `run_local.sh`, `submit_slurm.sh`: same shape as v2.1.

Expected wall-clock for the smoke test on an interactive node with HM3 (not
HM3+) on a chr-22-only subset: ~2 min. For full genome HM3+ on a real
sumstat: ~30–60 min on 16 cores.

### Comparison test

Add `tests/smoke/v2.2-<date>/compare_methods.R` that loads
`scores_sbayesr.gz` and `scores_ldpred2.gz` and reports the
sample-level Pearson correlation. For well-powered sumstats this should
be 0.8–0.95; the test asserts >0.6 as a sanity check.

---

## 12. Implementation phases (recommended order)

Each phase is a self-contained PR/commit and the pipeline stays runnable
between phases.

| Phase | Scope                                                                                   | Approx. effort |
|-------|----------------------------------------------------------------------------------------|----------------|
| 1     | **Container (multi-stage `r_builder`)**: add `r_builder` stage with `pak` + R packages; final stage `COPY --from=r_builder`s the lib tree; runtime libs only in final stage; bump tag; smoke-test `Rscript -e 'library(bigsnpr)'`. See [§8] | 1 day          |
| 2     | Config schema: parse `methods:`, `ldpred2.*`, `slurm.weights_ldpred2`, **`--methods` CLI override**. No behavior change yet. Validation: fail fast if `ldpred2` is in active set but `ldpred2.ld_dir` / `ldpred2.ld_meta_file` missing | 0.5 day        |
| 3     | `prep-ldref-ldpred2`: `map_*.rds` → TSV converter. Unit test.                            | 0.5–1 day      |
| 4     | `prep-inclusion-list-ldpred2`: build `prep/variant_map_ldpred2.tsv`. Migration symlink for old `variant_map.tsv`. | 0.5 day        |
| 5     | `filter-variants --method`: parameterise output dirs; per-method `.completed` markers. Migration helper for legacy `work/filtered/` → `work/filtered_sbayesr/`. | 1 day          |
| 6     | `bin/lib/scripts/run_ldpred2.R` + `bin/lib/steps/calc_ldpred2.sh` shell wrapper (deps check, validation, log routing, FAILED markers). Single-machine test on a small sumstat | 2 days         |
| 7     | `format-posteriors-ldpred2` (refactor `format_posteriors.sh` to be method-parameterised). Per-method `posteriors_mapped_<method>/.completed`. | 1 day          |
| 8     | Driver: `submit_single_job_for_step()`, dispatch `methods` in weights & score phases, **`--methods` CLI flag plumbed through** | 1 day          |
| 9     | `combine-scores` + `finalize-output` extension: per-method `scores_<method>.gz`; **discovery-driven** `augmented_sumstat.gz` column inclusion (enables incremental runs from [§6.5]) | 1.5 days       |
| 10    | Unit tests, smoke test, comparison test, **incremental run integration test** (Day-1 sBayesR → Day-2 LDpred2 → assert augmented_sumstat has both columns and scores_sbayesr.gz is byte-identical to Day-1) | 1.5–2 days     |
| 11    | Docs: update `README-v2.md` (add "Choosing methods" + "Incremental runs" sections), `config.template.yaml`, `docs/plans/general-design.md` | 0.5–1 day      |
| 12    | Version bump → v2.2.0, CHANGELOG entry, container retag                                  | 0.5 day        |

**Total:** ~11–13 person-days.

**Critical-path note:** Phase 1 (container) is independent of everything
else and can be done first to de-risk the R install. Phases 2–9 have a
near-linear dependency chain. Phase 10's incremental-run test depends on
Phase 9 (discovery-mode finalize) — it's the acceptance test for [§6.5].

---

## 13. Open questions for review

1. **LD reference set:** default to **HM3+** (1.44M variants) or HM3 (1.05M)?
   Recommendation: **HM3+**. The bigsnpr paper authors recommend it for
   sufficiently powered GWAS, and most of our sumstats fit that profile.
2. **LDpred2 mode default:** **auto** (recommended; no tuning set), `inf`,
   or `grid`? Grid needs a validation phenotype which we currently don't
   plumb through. Recommendation: auto by default, expose `inf` as a flag,
   defer grid to a later PR.
3. **lassosum2:** lassosum2 is one line away once we have `corr` and
   `df_beta` in the R script. It would give a free third score column.
   Skip in v2.2.0 to keep the PR small, or include? Recommendation: skip,
   schedule for v2.3.0.
4. **Score columns vs files:** emit `scores_sbayesr.gz` + `scores_ldpred2.gz`
   ([§9] Option A) — yes/no? Recommendation: yes (Option A).
5. **Per-variant N source:** the LDpred2 R script needs `n_eff`. If the
   per-variant `N` column in the filtered sumstat is missing/zero for many
   variants, do we fall back to (a) the sumstat metadata `stats_TotalN`
   filled per row, (b) the bigsnpr trick `quantile(8 / beta_se^2, 0.999)`,
   or (c) error out? Recommendation: (a) for filled rows, (b) only for
   diagnostics in the log. Never (c) silently — log a clear warning.
6. **Removing the precomputed-LD-ref-only restriction:** if a user has a
   non-European cohort, do we support computing the LD ref from PLINK files
   (using `snp_cor` in R)? Cost: ~hours of one-time CPU + tens of GB.
   Recommendation: out of scope for v2.2.0; document that LDpred2 is
   European-only initially.
7. **Backwards-compatible variant_map:** rename today's
   `prep/variant_map.tsv` → `prep/variant_map_sbayesr.tsv` and symlink? Yes
   (safest). Or in-place keep the name + add `_ldpred2.tsv` next to it?
   Recommendation: rename + symlink. The symlink makes back-compat
   automatic; the rename keeps the naming consistent.

---

## 14. References

- LDpred2 vignette: <https://privefl.github.io/bigsnpr/articles/LDpred2>
- LDpred2 paper (Privé et al. 2020): <https://doi.org/10.1093/bioinformatics/btaa1029>
- LDpred2-auto inference paper (Privé et al. 2023): <https://doi.org/10.1016/j.ajhg.2023.10.010>
- HM3+ paper: <https://doi.org/10.1016/j.ajhg.2023.10.010>
- LDpred2-author reference script:
  <https://github.com/privefl/paper-infer/blob/main/code/example-with-provided-LD.R>
- Comorment standalone CLI:
  <https://github.com/comorment/containers/blob/main/scripts/pgs/LDpred2/ldpred2.R>
- Comorment LD ref repo: <https://github.com/comorment/ldpred2_ref>
- Precomputed LD ref (HM3): <https://doi.org/10.6084/m9.figshare.19213299>
- Precomputed LD ref (HM3+): <https://doi.org/10.6084/m9.figshare.21305061>
