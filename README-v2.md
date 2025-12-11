# pgscalculator v2.0.0

Modular PGS calculation pipeline with step-by-step control for sbayesR-based polygenic scoring.

_Created by Jesper R. Gådin, Morten Dybdahl Krebs, and Andrew Schork (IBP)_

## What's New in v2

- **Modular CLI**: Run individual steps or groups of steps
- **Reusable prep steps**: Run prep once, reuse across multiple sumstats
- **Step-by-step control**: Skip completed steps, resume failed runs
- **Better logging**: Track progress with status command
- **Backwards compatible**: Wrapper script supports v1 command-line interface

## Quick Start

### Prerequisites

```bash
# Check required software
singularity --version  # or docker --version
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

### Run Full Pipeline (v1-compatible wrapper)

```bash
./pgscalculator-v2.sh \
  -i /path/to/cleansumstats/output/sumstat_TRAIT \
  -l /path/to/ld-reference/band_ukb_10k_hm3 \
  -g /path/to/genotypes \
  -f /path/to/genotype_manifest.txt \
  -c conf/sbayesr.config \
  -o /path/to/output
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
│           ├── prep_whitelist.sh
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
| 3 | `prep-whitelist` | Create variant whitelist (INFO/MAF filtered, LD intersect) |

### Per-Sumstat Processing Steps

| Step | Command | Description |
|------|---------|-------------|
| 4 | `format-sumstat` | Add build coordinates, derive B/SE/EAF/N |
| 5 | `filter-variants` | Filter sumstat to whitelist variants |
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

The wrapper script (`pgscalculator-v2.sh`) provides backwards compatibility with v1:

```bash
# Run all steps
./pgscalculator-v2.sh \
  -i /path/to/sumstat_TRAIT \
  -l /path/to/ld-reference \
  -g /path/to/genotypes \
  -f /path/to/manifest.txt \
  -c conf/sbayesr.config \
  -o /path/to/output

# Run specific steps
./pgscalculator-v2.sh \
  -i /path/to/sumstat_TRAIT \
  -l /path/to/ld-reference \
  -c conf/sbayesr.config \
  -o /path/to/output \
  --steps prep,posteriors

# Skip prep if already completed
./pgscalculator-v2.sh \
  -i /path/to/sumstat_TRAIT \
  -l /path/to/ld-reference \
  -c conf/sbayesr.config \
  -o /path/to/output \
  --steps score --skip-prep
```

#### Wrapper Script Options

| Option | Description |
|--------|-------------|
| `-i <dir>` | Path to cleansumstats output folder |
| `-l <dir>` | LD reference directory |
| `-g <dir>` | Target genotypes directory |
| `-f <file>` | Genotype manifest file |
| `-c <file>` | Config file (sbayesr.config) |
| `-o <dir>` | Output directory |
| `-j <mode>` | Container mode: `docker`, `dockerhub_biopsyk`, or `singularity` (default) |
| `-d` | Dev mode (keep intermediates) |
| `-v` | Show version |
| `-1` | Disable posterior calculation |
| `-2` | Disable scoring |
| `--steps <list>` | Comma-separated steps: `prep`, `sumstat`, `posteriors`, `score` |
| `--skip-prep` | Skip prep steps if already completed |

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
pgscalculator prep-whitelist --config /path/to/config.yaml

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
│   └── whitelist/
│       ├── variant_whitelist.tsv
│       └── variant_map.tsv  # rsid <-> genotype_id mapping
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
- Variant whitelist was created successfully
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
