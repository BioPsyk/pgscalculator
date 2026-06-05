# pgscalculator v2.2.0

Modular PGS calculation pipeline with step-by-step control for **sBayesR** and optional **LDpred2** polygenic scoring.

_Created by Jesper R. Gådin, Morten Dybdahl Krebs, and Andrew Schork (IBP)_

## What's New in v2.1

- **Config-first**: Reference paths in config.yaml, not CLI
- **Explicit steps**: `--steps` required — only runs what you specify
- **Prerequisite checks**: Helpful errors if prep/previous steps not done
- **Reusable prep**: Run prep once, reuse across multiple sumstats
- **Dual-position mapfile**: Prep builds a variant map with both GRCh37 and GRCh38 positions via `liftover_reference`
- **Finalize step**: Produces clean user-facing output files (`scores.gz`, `augmented_sumstat.gz`, `variant_map.gz`, `bench_score.gz`)
- **SLURM integration**: `--sbatch` flag auto-submits with config settings
- **SLURM driver jobs**: `--sbatch` submits one *driver job* per sumstat which runs `sumstat` and launches chromosome-parallel arrays for `weights` (sBayesR + benchmark, two arrays) and `score`, then a finalize job. Parallelism is controlled by `slurm.<step>.max_parallel` (set `max_parallel: 1` to disable parallelism).
- **Simplified CLI**: Just `--config`, `--steps`, and `-i`

## What's New in v2.2

- **LDpred2** (`calc-ldpred2`): genome-wide posterior estimation via `bigsnpr` (single job, not chr-parallel)
- **`methods:` config** and **`--methods` CLI**: choose `sbayesr`, `ldpred2`, or both per submission
- **Per-method outputs**: `scores_sbayesr.gz`, `scores_ldpred2.gz` (symlink `scores.gz` → `scores_sbayesr.gz` when sBayesR-only)
- **Per-method augmented sumstats**: each method gets its own self-contained `augmented_<method>.gz` (restricted to that method's LD-reference variant set, with its own `benchEffect` column inside); `augmented_sumstat.gz` is a back-compat symlink to the sBayesR file
- **Per-method benchmarks**: `bench_score_sbayesr.gz` / `bench_score_ldpred2.gz` (`bench_score.gz` → sBayesR, back-compat)
- **Discovery-driven, method-aware finalize**: finalize produces an augmented file for each method with mapped posteriors on disk and re-runs when a new method is added (incremental-safe)
- **LDpred2 diagnostics**: `details/ldpred2/summary.tsv` + `chains.png`
- **Incremental runs**: run sBayesR today and LDpred2 tomorrow in the same `outdir` without re-running sBayesR (see [Incremental runs](#incremental-runs))
- **New prep steps** (opt-in when `ldpred2.ld_dir` is set): `prep-ldref-ldpred2`, `prep-inclusion-list-ldpred2`
- **Smoke tests**: `tests/smoke/v2.2-2026-05-26/` (sBayesR-only, LDpred2-only, both, incremental Day-2 script + chr22 acceptance test)

## Quick Start

### Prerequisites

```bash
# Check required software
singularity --version
git --version

# Clone repository
git clone https://github.com/BioPsyk/pgscalculator.git
cd pgscalculator
```

### Pull Container Image

```bash
mkdir -p sif
singularity pull sif/ibp-pgscalculator-base_version-0.7.0.sif docker://biopsyk/ibp-pgscalculator:0.7.0-amd64
```

Image version (`docker/VERSION`, currently **0.7.0**) is **independent** of the pipeline version (`VERSION.v2`, **2.2.0**). Rebuild the image only when the Dockerfile changes (e.g. LDpred2 R packages added in Phase 1). sBayesR-only can use **0.6.0** until you enable LDpred2.

```bash
./scripts/docker-build.sh
singularity build sif/ibp-pgscalculator-base_version-0.7.0.sif docker-daemon://ibp-pgscalculator-base:0.7.0
./scripts/test-r-packages.sh --singularity sif/ibp-pgscalculator-base_version-0.7.0.sif
```

### Reference Data

The pipeline needs several reference data sets (liftover map, LD reference per
method, RSID map). Sources, download commands, and post-download sanity checks
for all of them live in [`docs/references.md`](docs/references.md). Install
what you need before running `prep`.

### Create Config File

Create `config.yaml` with your reference data paths (see `config.template.yaml` for all options):

```yaml
# config.yaml

# Optional: HPC modules to load before running (omit if singularity is already in PATH)
#modules:
#  - tools
#  - singularity/4.1.2

input: /path/to/cleansumstats/output
outdir: /path/to/output
genodir: /path/to/genotypes
genofile: /path/to/genotype_manifest.tsv
lddir: /path/to/ld_reference

# Genome build of the genotype files
genotype_build: GRCh37

# Liftover reference file (required for dual-position mapfile)
liftover_reference: /path/to/references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz

# Variant inclusion lists (optional)
filters:
  inclusion_list:
    gt: /path/to/genotype_snpids.txt
    ss: /path/to/sumstat_snpids.txt
    ld: /path/to/ldref_snpids.txt

whichn: totalN

methods: [sbayesr]

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  out_freq: 10
  p_value: 0.99
  rsq: 0.95
  threads: 6
  seed: 80851
  thin: 10
  exclude_mhc: true
  unscale_genotype: true
  no_mcmc_bin: false
  impute_n: false

# Optional LDpred2 (see docs/references.md §4)
#methods: [sbayesr, ldpred2]
#ldpred2:
#  mode: auto
#  ld_variant_set: hm3_plus
#  ld_dir: /path/to/references/ld-ldpred2/hm3_plus
#  threads: 16

# Scoring columns (variant_id allele effect)
score_columns: 1 2 5

# PLINK2 settings
plink:
  threads: 4

# Benchmark calculation (MAF and LD pruning for benchmark scores)
benchmark:
  maf_threshold: 0.05
  indep_pairwise: [250, 50, 0.25]   # window_kb, step, r2 for plink --indep-pairwise

# Optional: SLURM settings for --sbatch (weights = two array jobs, independently configured)
slurm:
  account: my_account
  partition: normal
  driver:             { mem: 1g, cpus: 1, time: '2:00:00' }
  prep:               { mem: 10g, cpus: 1, time: '1:00:00', max_parallel: 22 }
  sumstat:            { mem: 1g, cpus: 1, time: '0:30:00', max_parallel: 22 }
  weights_sbayesr:    { mem: 20g, cpus: 6, time: '2:00:00', max_parallel: 22 }
  weights_ldpred2:    { mem: 64g, cpus: 16, time: '4:00:00' }
  weights_benchmark:   { mem: 2g, cpus: 2, time: '0:30:00', max_parallel: 22 }
  score:              { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }
  finalize:           { mem: 16g, cpus: 1, time: '1:00:00' }
```

### Run Pipeline

```bash
# Step 1: Run prep (once per project)
./pgscalculator-v2.sh --config config.yaml --steps prep

# Step 2: Run per-sumstat steps (for each trait)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_TRAIT

# Or submit as SLURM jobs (uses slurm settings from config)
./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_TRAIT --sbatch

# Or submit as a SLURM driver job (recommended for running many sumstats in parallel)
# - prep must be run on its own
# - per sumstat: driver runs sumstat, weights (sBayesR + benchmark) and score arrays, then a separate finalize job (combine-scores + finalize-output)
# - parallelism is controlled via slurm.<step>.max_parallel (set 1 to serialize)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_TRAIT --sbatch

# You can also run subsets via the same driver mechanism:
./pgscalculator-v2.sh --config config.yaml --steps sumstat -i /path/to/sumstat_TRAIT --sbatch
./pgscalculator-v2.sh --config config.yaml --steps weights,score -i /path/to/sumstat_TRAIT --sbatch
```

### Batch Processing Multiple Sumstats

```bash
# Submit prep once
./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch

# Wait for prep to complete, then submit per-sumstat jobs
for sumstat in /path/to/sumstat_*; do
  ./pgscalculator-v2.sh --config config.yaml \
    --steps sumstat,weights,score,finalize \
    -i "$sumstat" \
    --sbatch
done

# Note: the per-sumstat output folder name matches the input sumstat folder basename
# (e.g. input: /path/to/sumstat_5759 -> outdir/sumstats/sumstat_5759/)
```

## Architecture

```
pgscalculator v2.2.0
├── pgscalculator-v2.sh      # Wrapper script (Singularity, SLURM, mounts)
├── bin/
│   ├── pgscalculator        # Main CLI entry point
│   └── lib/
│       ├── common.sh        # Shared functions
│       ├── scripts/run_ldpred2.R
│       └── steps/
│           ├── run_pipeline.sh       # Step orchestration
│           ├── prep_genotypes.sh
│           ├── prep_ldref.sh
│           ├── prep_ldref_ldpred2.sh
│           ├── prep_inclusion_list.sh
│           ├── prep_inclusion_list_ldpred2.sh
│           ├── format_sumstat.sh
│           ├── filter_variants.sh      # --method sbayesr|ldpred2
│           ├── calc_posteriors.sh
│           ├── calc_ldpred2.sh
│           ├── format_posteriors.sh    # method-parameterised
│           ├── calc_benchmark.sh
│           ├── calc_score.sh
│           ├── combine_scores.sh
│           ├── finalize_output.sh
│           └── status.sh
└── config.template.yaml     # Configuration template
```

## Pipeline Steps

### Preparation Steps (Run Once, Reuse Across Sumstats)

| Step | Command | Description |
|------|---------|-------------|
| 1 | `prep-genotypes` | Extract variant IDs from genotype .pvar files |
| 2 | `prep-ldref` | Augment sBayesR LD reference with dual positions (GRCh37/38) via liftover |
| 2b | `prep-ldref-ldpred2` | Validate LDpred2 LD RDS set and export `prep/ldref_ldpred2/map.tsv` (opt-in) |
| 3 | `prep-inclusion-list` | Create sBayesR variant map + inclusion list + compute MAF |
| 3b | `prep-inclusion-list-ldpred2` | Build `prep/variant_map_ldpred2.tsv` (opt-in) |

### Per-Sumstat Processing Steps

| Step | Command | Description |
|------|---------|-------------|
| 4 | `format-sumstat` | Add build coordinates, derive B/SE/EAF/N per chromosome |
| 5 | `filter-variants --method` | Filter sumstat per method (`work/filtered_sbayesr/`, `work/filtered_ldpred2/`) |
| 6 | `calc-posteriors` | Run sBayesR per chromosome (when `sbayesr` is active) |
| 6b | `calc-ldpred2` | Run LDpred2-auto genome-wide (when `ldpred2` is active) |
| 7 | `format-posteriors --method` | Map posteriors to genotype IDs per method |
| 8 | `calc-benchmark` | Benchmark scores (MAF filter + LD pruning; method-agnostic) |
| 9 | `calc-score --method` | Calculate PGS with plink2 per chromosome per method |
| 10 | `combine-scores` | Merge per-chr scores into `scores_<method>.gz` |
| 11 | `finalize-output` | Per-method `augmented_<method>.gz` + `bench_score_<method>.gz`, `variant_map.gz`, LDpred2 diagnostics (method-aware, incremental-safe) |

## Usage

### Option 1: Wrapper Script (Recommended for Batch Jobs)

The wrapper script (`pgscalculator-v2.sh`) uses a config-first approach:

```bash
# Run prep steps (once per project)
./pgscalculator-v2.sh --config config.yaml --steps prep

# Run per-sumstat steps
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_TRAIT

# Run specific steps only (prerequisite checking will warn if previous steps missing)
./pgscalculator-v2.sh --config config.yaml --steps weights,score,finalize -i /path/to/sumstat_TRAIT

# Submit as SLURM job (uses slurm settings from config)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_TRAIT --sbatch
```

#### Wrapper Script Options

| Option | Description |
|--------|-------------|
| `--config <file>` | **Required**: Path to config.yaml with reference paths |
| `--steps <list>` | **Required**: Steps to run: `prep`, `sumstat`, `weights`, `score`, `finalize` |
| `--methods <list>` | Active posterior methods for this submission: `sbayesr`, `ldpred2`, or `sbayesr,ldpred2` (overrides `methods:` in config) |
| `-i <dir>` | Path to sumstat folder (required for non-prep steps) |
| `-o <dir>` | Output directory (overrides config) |
| `--force` | Force re-run of steps even if already completed |
| `--sbatch` | Submit as SLURM job using slurm settings from config |
| `--cleanup` | Remove work/tmp folders after successful run |
| `-d` | Dev/verbose mode |
| `-v` | Show version |

#### Step Groups

| Group | Concrete steps | Description |
|-------|----------------|-------------|
| `prep` | `prep-genotypes`, `prep-ldref`, `prep-ldref-ldpred2`, `prep-inclusion-list`, `prep-inclusion-list-ldpred2` | Reference prep (LDpred2 prep steps run only when `ldpred2.ld_dir` is set) |
| `sumstat` | `format-sumstat`, `filter-variants` | Format sumstat; filter per active method |
| `weights` | sBayesR: `calc-posteriors`, `format-posteriors`; LDpred2: `calc-ldpred2`, `format-posteriors`; `calc-benchmark` | Posterior weights for active methods + benchmark |
| `score` | `calc-score` | PGS per chromosome per active method |
| `finalize` | `combine-scores`, `finalize-output` | Per-method score files, `augmented_<method>.gz`, `bench_score_<method>.gz` (method-aware) |

> **Note:** With `--sbatch`, sBayesR weights and benchmark run as chr-parallel arrays; LDpred2 weights run as a **single genome-wide** job (`weights_ldpred2`). Score arrays are submitted per active method (`score_sbayesr`, `score_ldpred2`).

#### Prerequisite Checking

The pipeline automatically checks that required outputs exist before running steps:

```
Error: Prep outputs not found.

Missing:
  - Genotype prep: /path/prep/genotypes/
  - Inclusion list: /path/prep/inclusion-list/inclusion_list.txt

Run prep first:
  ./pgscalculator-v2.sh --config config.yaml --steps prep

Then retry your command.
```

### Option 2: Direct CLI (Inside Container)

For interactive use or custom workflows:

```bash
# Start interactive container session
singularity shell --contain --cleanenv \
  -B /faststorage:/faststorage \
  sif/ibp-pgscalculator-base_version-0.7.0.sif

# Inside container - run individual steps
pgscalculator prep-genotypes --config /path/to/config.yaml
pgscalculator prep-ldref --config /path/to/config.yaml
pgscalculator prep-inclusion-list --config /path/to/config.yaml

pgscalculator format-sumstat --sumstat TRAIT --config /path/to/config.yaml
pgscalculator filter-variants --sumstat TRAIT --config /path/to/config.yaml
pgscalculator calc-posteriors --sumstat TRAIT --config /path/to/config.yaml
pgscalculator format-posteriors --sumstat TRAIT --config /path/to/config.yaml
pgscalculator calc-benchmark --sumstat TRAIT --config /path/to/config.yaml
pgscalculator calc-score --sumstat TRAIT --config /path/to/config.yaml
pgscalculator combine-scores --sumstat TRAIT --config /path/to/config.yaml
pgscalculator finalize-output --sumstat TRAIT --config /path/to/config.yaml

# Check status
pgscalculator status --config /path/to/config.yaml

# Run all steps at once
pgscalculator run --all --sumstat TRAIT --config /path/to/config.yaml
```

## Choosing methods

The `methods:` key (or `--methods` on the CLI) sets which posterior methods the **current submission** runs. Omitted `methods:` defaults to `[sbayesr]` (same as v2.1).

| Config / CLI | Behaviour |
|--------------|-----------|
| `methods: [sbayesr]` | sBayesR + benchmark only (default) |
| `methods: [ldpred2]` | LDpred2 + benchmark only; requires `ldpred2.ld_dir` |
| `methods: [sbayesr, ldpred2]` | Both methods in one submission (parallel SLURM weights jobs) |

LDpred2-specific settings live under `ldpred2:` in config (see `config.template.yaml`). Download and layout of the LD reference are documented in [`docs/references.md`](docs/references.md) (§4).

**Backwards compatibility:** configs without `methods:` behave as sBayesR-only. Existing `work/` layouts are migrated automatically (`work/filtered/` → `work/filtered_sbayesr/`, etc.).

## Incremental runs

On smaller HPC sites you can run **sBayesR and LDpred2 on different days** in the same `outdir` without re-running completed work.

**Day 1 — sBayesR only:**

```bash
./pgscalculator-v2.sh --config config.yaml \
  --steps sumstat,weights,score,finalize \
  --methods sbayesr \
  -i /path/to/sumstat_TRAIT
```

Produces `scores_sbayesr.gz` and `augmented_sbayesr.gz` (with `augmented_sumstat.gz` symlinked to it).

**Day 2 — add LDpred2** (same config `outdir`, same sumstat):

```bash
./pgscalculator-v2.sh --config config.yaml \
  --steps sumstat,weights,score,finalize \
  --methods ldpred2 \
  -i /path/to/sumstat_TRAIT
```

- **Include `sumstat`** on Day 2: `filter-variants` must run for LDpred2 (its LD-reference variant set differs from sBayesR's, so `work/filtered_ldpred2/` has to be built). `format-sumstat` is method-agnostic and is skipped as already complete.
- Runs the LDpred2 filter/weights/score path only; **does not** touch sBayesR `work/` outputs. `scores_sbayesr.gz` and `augmented_sbayesr.gz` stay **byte-identical** to Day 1.
- `finalize-output` is **discovery-driven and method-aware**: it writes `augmented_ldpred2.gz` (and `bench_score_ldpred2.gz`, `details/ldpred2/`) for the newly added method, and re-runs even though Day 1 already marked finalize complete. Each method keeps its own self-contained augmented file — there is no combined file with `postEffect_<method>` columns.

A scripted example (`run_incremental_day2.sh`) and an automated chr22 acceptance test with built-in byte-identity assertions (`integration_incremental_chr22.sh`) live under `tests/smoke/v2.2-2026-05-26/`.

To **re-run only LDpred2** (e.g. new LD ref), add `--force` with `--methods ldpred2`; sBayesR markers and outputs stay intact.

## Configuration

### config.yaml Format

See `config.template.yaml` for a complete annotated example. Key sections:

```yaml
# pgscalculator v2.2.0 Configuration

# Optional: HPC modules to load before running (omit if singularity is already in PATH)
#modules:
#  - tools
#  - singularity/4.1.2

input: /path/to/cleansumstats/output
outdir: /path/to/output
genodir: /path/to/genotypes
genofile: /path/to/genotype_manifest.tsv
lddir: /path/to/ld_reference

# Genome build of the genotype files (GRCh37 or GRCh38)
genotype_build: GRCh37

# Liftover reference file (required for dual-position mapfile)
# Pre-sorted on column 1 with LC_ALL=C
liftover_reference: /path/to/references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz

# Variant inclusion lists (optional)
filters:
  inclusion_list:
    gt: /path/to/genotype_snpids.txt
    ss: /path/to/sumstat_snpids.txt
    ld: /path/to/ldref_snpids.txt

whichn: totalN

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  out_freq: 10
  p_value: 0.99
  rsq: 0.95
  threads: 6
  seed: 80851
  thin: 10
  exclude_mhc: true
  unscale_genotype: true
  no_mcmc_bin: false
  impute_n: false

# Scoring columns (variant_id allele effect)
score_columns: 1 2 5

# PLINK2 settings
plink:
  threads: 4

# Benchmark calculation (MAF and LD pruning for benchmark scores)
benchmark:
  maf_threshold: 0.05
  indep_pairwise: [250, 50, 0.25]   # window_kb, step, r2 for plink --indep-pairwise

# Optional: SLURM settings for --sbatch
slurm:
  account: my_account
  partition: normal
  driver:             { mem: 1g, cpus: 1, time: '2:00:00' }
  prep:               { mem: 10g, cpus: 1, time: '1:00:00', max_parallel: 22 }
  sumstat:            { mem: 1g, cpus: 1, time: '0:30:00', max_parallel: 22 }
  weights_sbayesr:    { mem: 20g, cpus: 6, time: '2:00:00', max_parallel: 22 }
  weights_benchmark:   { mem: 2g, cpus: 2, time: '0:30:00', max_parallel: 22 }
  score:              { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }
  finalize:           { mem: 16g, cpus: 1, time: '1:00:00' }
```

## Testing

### Test Environment Setup

Request an interactive node:

```bash
srun --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --pty /bin/bash
```

### Run Smoke Test

**v2.1 (sBayesR only):** `tests/smoke/v2.1-2026-01-11/`

**v2.2 (LDpred2):** `tests/smoke/v2.2-2026-05-26/` — configs for sBayesR-only, LDpred2-only, both methods, and incremental Day-2.

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator
export LDPRED2_LD_DIR="${PWD}/references/ld-ldpred2/hm3"   # after download (see docs/references.md)

bash tests/smoke/v2.2-2026-05-26/run_local.sh config.sbayesr.yaml
bash tests/smoke/v2.2-2026-05-26/run_local.sh config.ldpred2.yaml
```

Unit tests (run all via `tests/run-unit-tests.sh`) cover methods config, driver dispatch + prep-step classifier, format-posteriors LDpred2, the LDpred2 variant-map join (`test_variant_map_for_ldpred2.sh`), phase-9 discovery, and the method-aware incremental finalize guard (`test_incremental_finalize.sh`). The chr22 incremental acceptance test is `tests/smoke/v2.2-2026-05-26/integration_incremental_chr22.sh`.

## Output Structure

```
outdir/
├── prep/                              # Preparation outputs (reusable across sumstats)
│   ├── genotypes/                     # Extracted genotype variant IDs
│   ├── ldref/                         # LD reference RSIDs per chromosome
│   ├── ldref_augmented/               # LD ref augmented with dual positions
│   ├── inclusion_list/                # Filtered variants + MAF
│   ├── variant_map/                   # Per-chromosome variant maps
│   ├── variant_map.tsv                # Combined variant map
│   ├── references/                    # Liftover reference copy
│   ├── details/                       # Prep step details
│   ├── logs/                          # Prep logs
│   └── tmp/                           # Temporary files (removed with --cleanup)
│
└── sumstats/
    └── sumstat_{name}/                # Per-sumstat outputs
        ├── scores_sbayesr.gz          # sBayesR PGS (IID, SCORE_SUM, ALLELE_CT, N_VARIANTS)
        ├── scores_ldpred2.gz          # LDpred2 PGS (when that method has been run)
        ├── scores.gz                  # Symlink → scores_sbayesr.gz (sBayesR-only / back-compat)
        ├── augmented_sbayesr.gz       # sBayesR set: RSID, alleles, B, SE, Z, P, EAF, MAF, postEffect, benchEffect
        ├── augmented_ldpred2.gz       # LDpred2 set: ... postEffect, postp_ldpred2, benchEffect (its own variant set)
        ├── augmented_sumstat.gz       # Symlink → augmented_sbayesr.gz (back-compat)
        ├── variant_map.gz             # Variant mapping (rsid ↔ genotype ID; sBayesR map)
        ├── bench_score_sbayesr.gz     # sBayesR benchmark PGS (IID, ALLELE_CT, SCORE1_SUM)
        ├── bench_score_ldpred2.gz     # LDpred2 benchmark PGS (its own variant set)
        ├── bench_score.gz             # Symlink → bench_score_sbayesr.gz (back-compat)
        ├── details/
        │   ├── steps.tsv              # Per-step variant counts
        │   ├── config.yaml            # Config used for this run
        │   ├── run_summary.txt        # Run summary
        │   └── ldpred2/               # LDpred2 diagnostics (when run)
        │       ├── summary.tsv        # h2, p, alpha, intercept, match/QC/chain counts
        │       └── chains.png         # LDpred2-auto p/h2 sampling paths (kept chains)
        ├── logs/                      # Per-sumstat logs
        └── work/                      # Intermediate files (removed with --cleanup)
            ├── formatted/             # Per-chr formatted sumstats
            ├── filtered_sbayesr/      # Per-chr filtered + matched (sBayesR)
            ├── filtered_ldpred2/      # Per-chr filtered + matched (LDpred2)
            ├── posteriors/            # Per-chr sBayesR posteriors
            ├── posteriors_ldpred2/    # LDpred2 .snpRes (genome-wide)
            ├── posteriors_mapped/     # sBayesR mapped posteriors
            ├── posteriors_mapped_ldpred2/
            ├── scores/                # Per-chr sBayesR plink2 scores
            ├── scores_ldpred2/        # Per-chr LDpred2 plink2 scores
            ├── scores_combined_sbayesr/
            ├── scores_combined_ldpred2/
            ├── benchmark_sbayesr/     # sBayesR benchmark intermediate files
            └── benchmark_ldpred2/     # LDpred2 benchmark intermediate files
```

## Troubleshooting

### "Required config key 'liftover_reference' missing"

Add the `liftover_reference` and `genotype_build` keys to your config.yaml. See `config.template.yaml`.

### "Missing augmented LD reference for chrN"

The liftover reference file could not be read or produced zero matches during `prep-ldref`. Check:
- The file path in `liftover_reference` is correct and readable
- The file is either plain text or gzipped (auto-detected)
- The file is pre-sorted on column 1 with `LC_ALL=C`

### "Config file not found"

The wrapper creates a container-specific config in the output directory. Check:
- Output directory is writable
- Container has write access to mounted output directory

### "No variants mapped"

The variant map couldn't match sbayesR RSIDs to genotype IDs. Check:
- Variant inclusion list was created successfully (`prep/inclusion_list/`)
- Genotype .pvar file uses expected ID format
- LD reference matches the expected HM3 variants

### "sbayesR failed on chromosome X"

Common causes:
- Too few variants after filtering
- Memory issues (increase `slurm.weights_sbayesr.mem` in config)
- LD matrix mismatch

Check logs in `outdir/sumstats/{name}/logs/`

### Scores are all zero or NA

- Posteriors may have failed (check `work/posteriors/chr*.snpRes`)
- Variant mapping failed (check `work/posteriors_mapped/`)
- Genotype IDs don't match (compare with `variant_map.gz`)

## Resource Requirements

Default SLURM resource settings (tunable via config):

| Step group | Memory | CPUs | Time |
|------------|--------|------|------|
| Prep | 10 GB | 1 | 1 h (once per project) |
| Sumstat | 1 GB | 1 | 30 min |
| Weights (sBayesR) | 20 GB | 6 | 2 h (per chr) |
| Weights (LDpred2) | 64 GB | 16 | 4 h (single genome-wide job) |
| Weights (benchmark) | 2 GB | 2 | 30 min (per chr) |
| Score | 10 GB | 4 | 30 min (per chr, per method) |
| Finalize | 16 GB | 1 | 1 h |

## More Documentation

- [General Design](docs/plans/general-design.md) — Pipeline architecture and design decisions
- [LDpred2 integration plan](docs/plans/ldpred2-integration.md) — Full design and implementation phases
- [Reference data](docs/references.md) — Download paths for genotypes, sBayesR LD, and LDpred2 LD

## Version History

- **v2.2.0** - LDpred2 method, `methods:` / `--methods`, per-method scores, discovery-mode finalize, incremental runs
- **v2.1.0** - Dual-position mapfile, finalize step, SLURM driver jobs
- **v2.0.0** - Modular CLI with step-by-step control
- **v1.x** - Nextflow-based monolithic pipeline (see [README.md](README.md))

## License

See LICENSE file for details.
