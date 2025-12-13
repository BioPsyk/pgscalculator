# SLURM Workflow Guide for pgscalculator v2.1

This guide shows how to run pgscalculator efficiently on a SLURM cluster.

## Key Concept: Prep Once, Run Many

The pipeline is designed so you can:
1. **Run prep steps once** - These process genotypes and LD reference (reusable)
2. **Run per-sumstat steps many times** - Each sumstat gets its own job

```
┌─────────────────────────────────────────────────────────────┐
│  PREP (run once per project)                                │
│  ├── prep-genotypes    → extracts variant IDs              │
│  ├── prep-ldref        → extracts LD RSIDs                 │
│  └── prep-inclusion    → creates variant filter list        │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  PER-SUMSTAT (run for each trait)                           │
│  ├── format-sumstat    → adds coordinates, derives stats   │
│  ├── filter-variants   → filters to inclusion list         │
│  ├── calc-posteriors   → runs sbayesR (slowest step)       │
│  ├── format-posteriors → maps to genotype IDs              │
│  ├── calc-score        → runs plink2 scoring               │
│  └── combine-scores    → merges chromosome scores          │
└─────────────────────────────────────────────────────────────┘
```

## CLI Overview

```bash
# Required: config file with paths
./pgscalculator-v2.sh --config config.yaml [options]

# Options:
#   -i <path>       Sumstat input path (absolute)
#   -o <path>       Output directory (overrides config)
#   --steps <list>  Steps to run (prep, sumstat, posteriors, score)
#   --skip-prep     Skip prep if already done
#   --chr <range>   Chromosomes to process
#   -d              Verbose/debug mode
```

---

## Quick Start

### Step 1: Create Project Config

Create a `config.yaml` file with all your paths:

```yaml
# config.yaml - PGS Project Configuration
# All paths are ABSOLUTE paths on the host system

# Project output directory
outdir: /faststorage/project/ibp_pipeline_pgscalculator/my_project

# Reference data paths (required)
ld_reference: /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/references/ld-sbayesr/ukb/band_ukb_10k_hm3
genotypes: /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/references/genotypes_test/plink2
genotype_manifest: /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/references/genotypes_test/mapfiles/plink2_genodir_genofiles.txt

# Filtering thresholds
info_threshold: 0.8
maf_threshold: 0.01

# sbayesR parameters
sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 6
  seed: 80851
  exclude_mhc: true

# SLURM settings for --sbatch flag
slurm:
  account: ibp_pipeline_cleansumstats
  partition: normal
  prep:       { mem: 10g, cpus: 6, time: '1:00:00' }
  sumstat:    { mem: 5g,  cpus: 2, time: '0:30:00' }
  posteriors: { mem: 20g, cpus: 8, time: '2:00:00' }
  score:      { mem: 10g, cpus: 4, time: '0:30:00' }
  default:    { mem: 20g, cpus: 8, time: '2:00:00' }

# Optional: limit chromosomes for testing
# chromosomes: "21-22"
```

```bash
# Create project directory and save config there
mkdir -p /faststorage/project/ibp_pipeline_pgscalculator/my_project
# Save the above YAML as config.yaml in that directory
```

### Step 2: Run Prep (Once)

**Option A: Using --sbatch (recommended)**

```bash
# Just add --sbatch flag - settings come from config.yaml
./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch
```

**Option B: Manual sbatch**

```bash
PGSFOLD="/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator"
CONFIG="/faststorage/project/ibp_pipeline_pgscalculator/my_project/config.yaml"

sbatch --mem=10g --cpus-per-task=6 --time=1:00:00 \
  --account=ibp_pipeline_cleansumstats \
  --job-name="pgs_prep" \
  --output="prep.out" \
  --error="prep.err" \
  --wrap="
${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} --steps prep
"
```

### Step 3: Run Per-Sumstat (For Each Trait)

**Option A: Using --sbatch (recommended)**

```bash
# Single sumstat - settings come from config.yaml
./pgscalculator-v2.sh --config config.yaml \
  -i /path/to/sumstat_814 \
  --skip-prep \
  --sbatch
```

**Option B: Manual sbatch**

```bash
PGSFOLD="/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator"
CONFIG="/faststorage/project/ibp_pipeline_pgscalculator/my_project/config.yaml"
SUMSTAT_LIB="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0"

# Single sumstat - use FULL PATH with -i
ID=814
sbatch --mem=20g --cpus-per-task=8 --time=2:00:00 \
  --account=ibp_pipeline_cleansumstats \
  --job-name="pgs_${ID}" \
  --output="pgs_${ID}.out" \
  --error="pgs_${ID}.err" \
  --wrap="
${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} -i ${SUMSTAT_LIB}/sumstat_${ID} --skip-prep
"
```

---

## Batch Processing Multiple Sumstats

### Option A: Loop with Delay

```bash
# Create list of FULL PATHS to sumstats
cat > sumstat_paths.txt << EOF
/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_814
/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_815
/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_816
/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5668
EOF

PGSFOLD="/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator"
CONFIG="/faststorage/project/ibp_pipeline_pgscalculator/my_project/config.yaml"

# Submit jobs with small delay
while read SUMSTAT_PATH; do
  ID=$(basename "$SUMSTAT_PATH" | sed 's/^sumstat_//')
  sbatch --mem=20g --cpus-per-task=8 --time=2:00:00 \
    --account=ibp_pipeline_cleansumstats \
    --job-name="pgs_${ID}" \
    --output="pgs_${ID}.out" \
    --error="pgs_${ID}.err" \
    --wrap="
${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} -i ${SUMSTAT_PATH} --skip-prep
"
  sleep 0.5  # Small delay between submissions
done < sumstat_paths.txt
```

### Option B: SLURM Job Array

```bash
# Create a batch script
cat > run_pgs_array.sh << 'SCRIPT'
#!/bin/bash
#SBATCH --array=1-4
#SBATCH --mem=20g
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --account=ibp_pipeline_cleansumstats
#SBATCH --job-name=pgs_array
#SBATCH --output=pgs_array_%a.out
#SBATCH --error=pgs_array_%a.err

PGSFOLD="/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator"
CONFIG="/faststorage/project/ibp_pipeline_pgscalculator/my_project/config.yaml"

# Get sumstat path from line number
SUMSTAT_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sumstat_paths.txt)

echo "Processing: ${SUMSTAT_PATH}"
${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} -i ${SUMSTAT_PATH} --skip-prep
echo "Completed: ${SUMSTAT_PATH}"
SCRIPT

# Submit array job
sbatch run_pgs_array.sh
```

---

## Resource Requirements

| Step | Memory | CPUs | Time (est.) |
|------|--------|------|-------------|
| `prep` | 10 GB | 6 | 30-60 min |
| `sumstat` | 5 GB | 2 | 5-10 min |
| `posteriors` | 20 GB | 8 | 30-90 min |
| `score` | 10 GB | 4 | 10-20 min |
| **Full pipeline** | **20 GB** | **8** | **1-2 hours** |

### Tips for Large Runs

1. **Use job dependencies** if prep not done:
   ```bash
   PREP_JOB=$(sbatch --parsable --wrap="${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} --steps prep")
   sbatch --dependency=afterok:${PREP_JOB} --wrap="${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} --sumstat 814 --skip-prep"
   ```

2. **Limit concurrent jobs** to avoid I/O bottlenecks:
   ```bash
   sbatch --array=1-100%10  # Max 10 concurrent
   ```

3. **Reduce chromosomes for testing** (in config.yaml):
   ```yaml
   chromosomes: "21-22"
   ```

---

## Checking Status

### Job Status
```bash
squeue -u $USER
```

### Pipeline Status (inside container)
```bash
singularity shell --contain --cleanenv \
  -B /faststorage:/faststorage \
  ${PGSFOLD}/sif/ibp-pgscalculator-base_version-2.0.0.sif

# Inside container:
pgscalculator status --config ${OUTDIR}/config.yaml
```

### Check Outputs
```bash
# List completed sumstats
ls -d ${OUTDIR}/sumstat_*/scores_combined/*.sscore

# Check a specific sumstat
ls -la ${OUTDIR}/sumstat_814/
```

---

## Troubleshooting

### Job Failed - How to Resume

```bash
# Check which step failed
tail -100 pgs_814.err

# Re-run just the failed step (e.g., posteriors)
sbatch --wrap="${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} --sumstat 814 --steps posteriors --skip-prep"

# Or force re-run all steps by removing output
rm -rf $(grep outdir ${CONFIG} | cut -d: -f2 | tr -d ' ')/sumstat_814
sbatch --wrap="${PGSFOLD}/pgscalculator-v2.sh --config ${CONFIG} --sumstat 814 --skip-prep"
```

### Out of Memory

Increase `--mem`:
```bash
sbatch --mem=40g ...
```

Or reduce sbayesR threads in config.yaml:
```yaml
sbayesr:
  threads: 4  # instead of 6
```

### Slow Jobs

Reduce chain length for testing in config.yaml:
```yaml
sbayesr:
  chain_length: 5000  # instead of 10000
```

Or run subset of chromosomes in config.yaml:
```yaml
chromosomes: "21-22"  # instead of all chromosomes
```

