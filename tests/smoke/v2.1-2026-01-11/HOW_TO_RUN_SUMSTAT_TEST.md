# How to Run a Smoke Test for a Specific Sumstat

This document describes how to run the pgscalculator pipeline against a specific sumstat for testing/debugging purposes.

## Prerequisites

1. The pgscalculator environment is available
2. The sumstat has been processed by cleansumstats and exists in the sumstat clean library

## Quick Reference

### Sumstat Library Location

```
/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/
```

### Running a Test

From the pgscalculator root directory:

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

# Set the sumstat to test
export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_XXXX"

# Run interactively (recommended for debugging)
bash tests/smoke/v2.1-2026-01-11/run_local.sh

# Or submit to SLURM
bash tests/smoke/v2.1-2026-01-11/submit_slurm.sh
```

## Config Files Available

| Config File | Description |
|-------------|-------------|
| `config.yaml` | Quick test: chromosomes 21-22 only |
| `config.full.yaml` | Full genome: chromosomes 1-22 |
| `config.full.noimpute.yaml` | Full genome without imputation |
| `config.full.noimpute.v2parallel.yaml` | **Recommended**: Full genome, no imputation, v2.1 format |
| `config.sumstat_805.paramsearch_default.yaml` | Parameter search example |

**Recommended config**: `config.full.noimpute.v2parallel.yaml` - this is the most complete and up-to-date config, matching the README-v2.md specification with:
- `whichn: totalN`
- `score_columns: 1 2 5`
- `plink.threads: 4`
- `sbayesr.impute_n: false`

## Running with a Specific Config

To use a different config file, you can either:

1. **Modify run_local.sh temporarily** (not recommended for git-tracked changes)
2. **Run the wrapper directly:**

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_XXXX"
CONFIG="tests/smoke/v2.1-2026-01-11/config.yaml"

# Step 1: Run prep (only needed once per config, shared across sumstats)
./pgscalculator-v2.sh --config "$CONFIG" --steps prep

# Step 2: Run sumstat + posteriors + score for your sumstat
./pgscalculator-v2.sh --config "$CONFIG" --steps sumstat,posteriors,score -i "$SUMSTAT_DIR"
```

## Running Individual Steps

For debugging filtering issues, you may want to run only the `sumstat` step:

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_XXXX"
CONFIG="tests/smoke/v2.1-2026-01-11/config.yaml"

# Run only sumstat step (filter variants)
./pgscalculator-v2.sh --config "$CONFIG" --steps sumstat -i "$SUMSTAT_DIR"
```

## Output Location

All outputs go to the `outdir` specified in the config:

```
tests/smoke/v2.1-2026-01-11/out/
├── prep/                           # Shared prep outputs
└── sumstats/
    └── sumstat_XXXX/
        ├── logs/                   # Step logs
        │   └── slurm/              # SLURM logs (if using --sbatch)
        ├── work/                   # Intermediate files
        ├── tmp/                    # Temporary files
        └── results/                # Final results
```

## Checking Filtering Results

After running the `sumstat` step, check for filtering issues:

```bash
OUTDIR="tests/smoke/v2.1-2026-01-11/out/sumstats/sumstat_XXXX"

# Check per-chromosome variant counts
wc -l "$OUTDIR"/work/sumstat/*.filtered.gz 2>/dev/null || echo "No filtered files yet"

# Check logs for warnings/errors
grep -i "error\|warn\|0 variants" "$OUTDIR"/logs/*.log 2>/dev/null || true
```

## Example: Testing sumstat_5789

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5789"
CONFIG="tests/smoke/v2.1-2026-01-11/config.full.noimpute.v2parallel.yaml"

# Step 1: Run prep (only needed once)
./pgscalculator-v2.sh --config "$CONFIG" --steps prep

# Step 2: Run sumstat step to test filtering
./pgscalculator-v2.sh --config "$CONFIG" --steps sumstat -i "$SUMSTAT_DIR"

# Or run full pipeline (sumstat + posteriors + score)
./pgscalculator-v2.sh --config "$CONFIG" --steps sumstat,posteriors,score -i "$SUMSTAT_DIR"
```

## Cleanup

To remove test outputs:

```bash
rm -rf tests/smoke/v2.1-2026-01-11/out/sumstats/sumstat_XXXX
```

To remove all test outputs:

```bash
rm -rf tests/smoke/v2.1-2026-01-11/out/
```
