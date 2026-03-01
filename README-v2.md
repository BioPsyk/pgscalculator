# pgscalculator v2.1.0

Modular PGS calculation pipeline with step-by-step control for sbayesR-based polygenic scoring.

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
singularity pull sif/ibp-pgscalculator-base_version-2.0.0.sif docker://biopsyk/ibp-pgscalculator:2.0.0-amd64
```

### Create Config File

Create `config.yaml` with your reference data paths (see `config.template.yaml` for all options):

```yaml
# config.yaml
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

# Optional: SLURM settings for --sbatch (weights = two array jobs, independently configured)
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
pgscalculator v2.1.0
├── pgscalculator-v2.sh      # Wrapper script (Singularity, SLURM, mounts)
├── bin/
│   ├── pgscalculator        # Main CLI entry point
│   └── lib/
│       ├── common.sh        # Shared functions
│       └── steps/
│           ├── run_pipeline.sh       # Step orchestration
│           ├── prep_genotypes.sh
│           ├── prep_ldref.sh
│           ├── prep_inclusion_list.sh
│           ├── format_sumstat.sh
│           ├── filter_variants.sh
│           ├── calc_posteriors.sh
│           ├── format_posteriors.sh
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
| 2 | `prep-ldref` | Augment LD reference with dual positions (GRCh37/38) via liftover |
| 3 | `prep-inclusion-list` | Create variant inclusion list + variant map + compute MAF |

### Per-Sumstat Processing Steps

| Step | Command | Description |
|------|---------|-------------|
| 4 | `format-sumstat` | Add build coordinates, derive B/SE/EAF/N per chromosome |
| 5 | `filter-variants` | Filter sumstat to inclusion list variants per chromosome |
| 6 | `calc-posteriors` | Run sbayesR per chromosome |
| 7 | `format-posteriors` | Map posteriors to genotype variant IDs |
| 8 | `calc-benchmark` | Benchmark scores (MAF filter + LD pruning) |
| 9 | `calc-score` | Calculate PGS with plink2 per chromosome |
| 10 | `combine-scores` | Merge per-chromosome scores into `scores.gz` |
| 11 | `finalize-output` | Produce `augmented_sumstat.gz`, `variant_map.gz`, `bench_score.gz`, `steps.tsv` |

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
| `prep` | `prep-genotypes`, `prep-ldref`, `prep-inclusion-list` | Prepare genotypes, LD reference, and variant map (run once per project) |
| `sumstat` | `format-sumstat`, `filter-variants` | Format and filter sumstat per chromosome |
| `weights` | `calc-posteriors`, `format-posteriors`, `calc-benchmark` | sBayesR posteriors + benchmark weights |
| `score` | `calc-score` | Calculate PGS scores per chromosome |
| `finalize` | `combine-scores`, `finalize-output` | Merge scores and produce final output files |

> **Note:** `weights` runs both sBayesR (calc-posteriors, format-posteriors) and benchmark (calc-benchmark) as two separate SLURM array jobs when using `--sbatch`.

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
  sif/ibp-pgscalculator-base_version-2.0.0.sif

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

## Configuration

### config.yaml Format

See `config.template.yaml` for a complete annotated example. Key sections:

```yaml
# pgscalculator v2.1.0 Configuration
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

See `tests/smoke/v2.1-2026-01-11/` for a working smoke test setup:

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

# Local run (interactive node)
bash tests/smoke/v2.1-2026-01-11/run_local.sh

# SLURM run
bash tests/smoke/v2.1-2026-01-11/submit_slurm.sh
```

The smoke test config (`tests/smoke/v2.1-2026-01-11/config.yaml`) demonstrates all required config keys.

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
        ├── scores.gz                  # Final PGS scores (IID, SCORE_SUM, ALLELE_CT, N_VARIANTS)
        ├── augmented_sumstat.gz       # Variant-level results (RSID, B, SE, Z, P, EAF, MAF, postEffect, benchEffect, ...)
        ├── variant_map.gz             # Variant mapping (rsid ↔ genotype ID)
        ├── bench_score.gz             # Benchmark scores (IID, ALLELE_CT, SCORE1_SUM)
        ├── details/
        │   ├── steps.tsv              # Per-step variant counts
        │   ├── config.yaml            # Config used for this run
        │   └── run_summary.txt        # Run summary
        ├── logs/                      # Per-sumstat logs
        └── work/                      # Intermediate files (removed with --cleanup)
            ├── formatted/             # Per-chr formatted sumstats
            ├── filtered/              # Per-chr filtered sumstats
            ├── posteriors/            # Per-chr sbayesR posteriors
            ├── posteriors_mapped/     # Per-chr posteriors mapped to geno IDs
            ├── scores/               # Per-chr plink2 scores
            ├── scores_combined/      # Combined score file
            └── benchmark/            # Benchmark intermediate files
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
| Weights (benchmark) | 2 GB | 2 | 30 min (per chr) |
| Score | 10 GB | 4 | 30 min (per chr) |
| Finalize | 16 GB | 1 | 1 h |

## More Documentation

- [General Design](docs/plans/general-design.md) — Pipeline architecture and design decisions

## Version History

- **v2.0.0** - Modular CLI with step-by-step control
- **v1.x** - Nextflow-based monolithic pipeline (see [README.md](README.md))

## License

See LICENSE file for details.
