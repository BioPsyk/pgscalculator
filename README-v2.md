# pgscalculator v2.1.0

Modular PGS calculation pipeline with step-by-step control for sbayesR-based polygenic scoring.

_Created by Jesper R. Gådin, Morten Dybdahl Krebs, and Andrew Schork (IBP)_

## What's New in v2.1

- **Config-first**: Reference paths in config.yaml, not CLI
- **Explicit steps**: `--steps` required - only runs what you specify
- **Prerequisite checks**: Helpful errors if prep/previous steps not done
- **Reusable prep**: Run prep once, reuse across multiple sumstats
- **SLURM integration**: `--sbatch` flag auto-submits with config settings
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

Create `config.yaml` with your reference data paths:

```yaml
# config.yaml
outdir: /path/to/output
ld_reference: /path/to/ld-reference/band_ukb_10k_hm3
genotypes: /path/to/genotypes
genotype_manifest: /path/to/genotype_manifest.txt

info_threshold: 0.8
maf_threshold: 0.01

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 6
  seed: 80851
  exclude_mhc: true

# Optional: SLURM settings for --sbatch
slurm:
  account: my_account
  prep:       { mem: 10g, cpus: 6, time: '1:00:00' }
  posteriors: { mem: 20g, cpus: 8, time: '2:00:00' }
  score:      { mem: 10g, cpus: 4, time: '0:30:00' }
```

### Run Pipeline

```bash
# Step 1: Run prep (once per project)
./pgscalculator-v2.sh --config config.yaml --steps prep

# Step 2: Run per-sumstat steps (for each trait)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_TRAIT

# Or submit as SLURM jobs (uses slurm settings from config)
./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_TRAIT --sbatch
```

### Batch Processing Multiple Sumstats

```bash
# Submit prep once
./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch

# Wait for prep to complete, then submit per-sumstat jobs
for sumstat in /path/to/sumstat_*; do
  ./pgscalculator-v2.sh --config config.yaml \
    --steps sumstat,posteriors,score \
    -i "$sumstat" \
    --sbatch
done
```

## Architecture

```
pgscalculator v2.0.0
├── pgscalculator-v2.sh      # Wrapper script (v1 CLI compatible)
├── bin/
│   ├── pgscalculator        # Main CLI entry point
│   └── lib/
│       ├── common.sh        # Shared functions
│       └── steps/           # Individual step scripts
│           ├── prep_genotypes.sh
│           ├── prep_ldref.sh
│           ├── prep_inclusion_list.sh
│           ├── format_sumstat.sh
│           ├── filter_variants.sh
│           ├── calc_posteriors.sh
│           ├── format_posteriors.sh
│           ├── calc_score.sh
│           ├── combine_scores.sh
│           └── calc_benchmark.sh
└── config.template.yaml     # Configuration template
```

## Pipeline Steps

### Preparation Steps (Run Once, Reuse Across Sumstats)

| Step | Command | Description |
|------|---------|-------------|
| 1 | `prep-genotypes` | Extract variant IDs from genotype .pvar files |
| 2 | `prep-ldref` | Extract RSIDs from LD reference |
| 3 | `prep-inclusion-list` | Create variant inclusion list (INFO/MAF filtered, LD intersect) |

### Per-Sumstat Processing Steps

| Step | Command | Description |
|------|---------|-------------|
| 4 | `format-sumstat` | Add build coordinates, derive B/SE/EAF/N |
| 5 | `filter-variants` | Filter sumstat to inclusion list variants |
| 6 | `calc-posteriors` | Run sbayesR per chromosome |
| 7 | `format-posteriors` | Map posteriors to genotype variant IDs |
| 8 | `calc-score` | Calculate PGS with plink2 per chromosome |
| 9 | `combine-scores` | Merge per-chromosome scores |

### Optional

| Step | Command | Description |
|------|---------|-------------|
| - | `calc-benchmark` | Calculate benchmark scores (MAF filter + LD pruning) |

## Usage

### Option 1: Wrapper Script (Recommended for Batch Jobs)

The wrapper script (`pgscalculator-v2.sh`) uses a config-first approach:

```bash
# Run prep steps (once per project)
./pgscalculator-v2.sh --config config.yaml --steps prep

# Run per-sumstat steps
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_TRAIT

# Run specific steps only (prerequisite checking will warn if previous steps missing)
./pgscalculator-v2.sh --config config.yaml --steps posteriors,score -i /path/to/sumstat_TRAIT

# Limit to specific chromosomes (for testing)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_TRAIT --chr 21-22

# Submit as SLURM job (uses slurm settings from config)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_TRAIT --sbatch
```

#### Wrapper Script Options

| Option | Description |
|--------|-------------|
| `--config <file>` | **Required**: Path to config.yaml with reference paths |
| `--steps <list>` | **Required**: Steps to run: `prep`, `sumstat`, `posteriors`, `score` |
| `-i <dir>` | Path to sumstat folder (required for non-prep steps) |
| `-o <dir>` | Output directory (overrides config) |
| `--chr <range>` | Chromosome range (e.g., "21-22") |
| `--sbatch` | Submit as SLURM job using slurm settings from config |
| `-d` | Dev/verbose mode |
| `-v` | Show version |

#### Step Groups

| Step | Description |
|------|-------------|
| `prep` | Prepare genotypes and LD reference (run once per project) |
| `sumstat` | Format and filter sumstat |
| `posteriors` | Calculate posteriors with sbayesR |
| `score` | Calculate PGS scores |

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
pgscalculator calc-posteriors --sumstat TRAIT --config /path/to/config.yaml
pgscalculator calc-score --sumstat TRAIT --config /path/to/config.yaml

# Check status
pgscalculator status --config /path/to/config.yaml

# Run all steps at once
pgscalculator run --all --sumstat TRAIT --config /path/to/config.yaml
```

## Configuration

### config.yaml Format

```yaml
# pgscalculator v2.0.0 Configuration
input: /path/to/cleansumstats/output
outdir: /path/to/output
genodir: /path/to/genotypes
genofile: /path/to/genotype_manifest.tsv
lddir: /path/to/ld_reference

# Filtering thresholds
info_threshold: 0.8
maf_threshold: 0.01
whichn: totalN

# sbayesR parameters
sbayesr:
  gamma: 0.0,0.01,0.1,1
  pi: 0.95,0.02,0.02,0.01
  burn_in: 2000
  chain_length: 10000
  threads: 6
  seed: 80851
  exclude_mhc: true

# Scoring columns (variant_id allele effect)
score_columns: 2 5 9
```

## Testing

### Test Environment Setup

Request an interactive node:

```bash
srun --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --pty /bin/bash
```

### Test Paths (GDK Environment)

```bash
# Project location
PGSFOLD="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator"

# Sumstats (cleansumstats output)
SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.6.7"

# LD reference (sbayesR)
LD_REF="${PGSFOLD}/references/ld-sbayesr/ukb/band_ukb_10k_hm3"

# Test genotypes
GENO_DIR="${PGSFOLD}/references/genotypes_test/plink"
GENO_MANIFEST="${PGSFOLD}/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt"
```

### Run Test

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

# Test with example sumstat
./pgscalculator-v2.sh \
  -i ${SUMSTAT_DIR}/sumstat_TEST_ID \
  -l ${LD_REF} \
  -g ${GENO_DIR} \
  -f ${GENO_MANIFEST} \
  -c conf/sbayesr.config \
  -o ../test-zone/out_v2_test
```

### Verify Container Mounts

```bash
# Inside container
bash /pgscalculator/test-mounts.sh
```

### Batch Testing (SLURM)

See `test-zone-commands-v2.sh` for a template batch script:

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/test-zone
cp ${PGSFOLD}/test-zone-commands-v2.sh .

# Create a list of sumstat IDs to test
echo "TEST_ID_1" > short_list.txt
echo "TEST_ID_2" >> short_list.txt

# Submit batch jobs
bash test-zone-commands-v2.sh
```

## Output Structure

```
output_dir/
├── config.yaml              # Generated configuration
├── prep/                    # Preparation outputs (reusable)
│   ├── genotypes/
│   │   └── snplist_sorted   # All genotype variant IDs
│   ├── ldref/
│   │   └── chr*_ld_rsids    # LD reference RSIDs per chromosome
│   ├── inclusion_list/
│   │   └── variant_inclusion_list.tsv  # Filtered variants for analysis
│   └── variant_map.tsv      # Full rsid <-> genotype_id crosswalk
├── sumstat_{name}/          # Per-sumstat outputs
│   ├── formatted/
│   │   └── sumstat_formatted.tsv
│   ├── filtered/
│   │   └── sumstat_filtered.tsv
│   ├── posteriors/
│   │   └── chr*.snpRes      # sbayesR posteriors
│   ├── posteriors_mapped/
│   │   └── chr*.snpRes      # Mapped to genotype IDs
│   ├── scores/
│   │   └── chr*.sscore      # Per-chromosome scores
│   └── scores_combined/
│       └── merged.sscore    # Final combined scores
└── logs/                    # Execution logs
```

## Troubleshooting

### "pgscalculator: command not found"

Ensure the container was built with the v2 CLI:

```dockerfile
COPY bin/ /pgscalculator/bin/
ENV PATH="/pgscalculator/bin:${PATH}"
```

### "Config file not found"

The wrapper creates `config.yaml` in the output directory. Check:
- Output directory is writable
- Container has write access to mounted output directory

### "No variants mapped"

The variant map couldn't match sbayesR RSIDs to genotype IDs. Check:
- Variant inclusion list was created successfully
- Genotype .pvar file uses expected ID format
- LD reference matches the expected HM3 variants

### "sbayesR failed on chromosome X"

Common causes:
- Too few variants after filtering
- Memory issues (increase `--mem` in SLURM)
- LD matrix mismatch

Check logs in `{outdir}/sumstat_{name}/posteriors/chr{X}.log`

### Scores are all zero or NA

- Posteriors may have failed (check `*.snpRes` files)
- Variant mapping failed (check `posteriors_mapped/` files)
- Genotype IDs don't match (compare with `variant_map.tsv`)

## Resource Requirements

| Step | Memory | CPUs | Time (per sumstat) |
|------|--------|------|-------------------|
| Prep (all) | 10GB | 6 | 30 min (once) |
| Posteriors | 20GB | 22 | 30-60 min |
| Scoring | 10GB | 6 | 10-20 min |
| **Full pipeline** | **20GB** | **22** | **1-2 hours** |

## More Documentation

- [Testing v2 Setup](docs/testing-v2-setup.md) - Detailed testing guide
- [Technical Whitepaper](PGS_DST_whitepaper.md) - Full pipeline specification
- [Standalone Scoring](tmp/README_STANDALONE_SCORING.md) - Score external sbayesR outputs
- [SNP Inclusion List](docs/snp-inclusion-list.md) - Variant filtering details
- [FAQ](docs/FAQ.md) - Frequently asked questions

## Version History

- **v2.0.0** - Modular CLI with step-by-step control
- **v1.x** - Nextflow-based monolithic pipeline (see [README.md](README.md))

## License

See LICENSE file for details.
